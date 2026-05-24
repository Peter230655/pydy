"""The System class manages the simulation (integration) of a system whose
equations are given by :py:class:`KanesMethod`.

Many of the attributes are also properties, and can be directly modified.

Here is the procedure for using this class.

    1. specify your options either via the constructor or via the
       attributes.
    2. optionally, call :py:meth:`~pydy.system.System.generate_ode_function` if
       you want to customize how the ODE function is generated.
    3. call :py:meth:`~pydy.system.System.integrate` to simulate your system.

The simplest usage of this class is as follows. First, we need a
:external+sympy:py:class:`~sympy.physics.mechanics.kane.KanesMethod` object on
which we have already invoked
:external+sympy:py:meth:`~sympy.physics.mechanics.kane.KanesMethod.kanes_equations`::

    km = KanesMethod(...)
    km.kanes_equations(force_list, body_list)
    times = np.linspace(0, 5, 100)
    sys = System(km, times=times)
    sys.integrate()

In this case, we use defaults for the numerical values of the constants,
specified quantities, initial conditions, etc. You probably won't like
these defaults. You can also specify such values via constructor keyword
arguments or via the attributes::

    sys = System(km,
                 initial_conditions={dynamicsymbol('q1'): 0.5},
                 times=times)
    sys.constants = {symbol('m'): 5.0}
    sys.integrate()

To double-check the constants, specifieds, states and times in your problem,
look at these properties::

    sys.constants_symbols
    sys.specifieds_symbols
    sys.states
    sys.times

In this case, the System generates the numerical ode function for you
behind the scenes. If you want to customize how this function is generated,
you must call ``generate_ode_function`` on your own::

    sys = System(KM)
    sys.generate_ode_function(generator='cython')
    sys.integrate()

"""
import warnings
from itertools import repeat

import numpy as np
import sympy as sm
from sympy.physics.mechanics import dynamicsymbols, find_dynamicsymbols
from scipy.integrate import odeint
from scipy.optimize import root

from .codegen.ode_function_generators import generate_ode_function
from .utils import PyDyFutureWarning, PyDyUserWarning

SYMPY_VERSION = sm.__version__


warnings.simplefilter('once', PyDyFutureWarning)


class System(object):
    """Multibody dynamics system for simulation and numerical evaluation.

    See the class's attributes for a description of the arguments to this
    constructor.

    The parameters to this constructor are all attributes of the System.
    Actually, they are properties. With the exception of ``eom_method``,
    these attributes can be modified directly at any future point.

    Parameters
    ----------
    eom_method : sympy.physics.mechanics.KanesMethod
        You must have called ``KanesMethod.kanes_equations()`` *before*
        constructing this ``System``.
    constants : dict, optional (default: all 1.0)
        This dictionary maps SymPy Symbol objects to floats.
    specifieds : dict, optional (default: all 0.0)
        This dictionary maps SymPy Functions of time objects, or tuples of
        them, to floats, NumPy arrays, or functions of the state and time.
    ode_solver : function, optional (default: scipy.integrate.odeint)
        This function computes the derivatives of the states.
    initial_conditions : dict, optional (default: all zero)
        This dictionary maps SymPy Functions of time objects to floats.
    times : array_like, shape(n,), optional
        An array_like object, which contains time values over which
        equations are integrated. It has to be supplied before
        System.integrate() can be called.
    constraint_loads : iterable of Functions of time, optional
        If the ``eom_method`` includes noncontributig forces (Kane's method) or
        constraint forces from Lagrange multipliers, provide a list of variable
        names for these forces and they will be computed when evaluating the
        differential equations.
    outputs : dictionary
        Three possible entries for this dictionary:

            1 .{Function('y')(t): SymPy Expression}
            2. {(Function('y1')(t), Function('y2')(t): SymPy Column Matrix
               where y1 and y2 are not present}
            3. {(Function('f1')(t), Function('f2')(t): SymPy Column Matrix
               where both expressions are linear f1, f2 and possilby the time
               derivative of the speeds}

    """
    def __init__(self, eom_method, constants=None, specifieds=None,
                 ode_solver=None, initial_conditions=None, times=None,
                 outputs=None, constraint_loads=None):

        self._eom_method = eom_method
        self._extract_constraints()

        if outputs is None:
            outputs = dict()

        if self.constraints:
            self._con_syms = tuple(sm.Dummy('c' + str(i)) for i in
                                   range(self.num_constraints))
            outputs[self._con_syms] = self.constraints

        self.outputs = outputs  # calls _parse_outputs

        # NOTE: must be set before the state variables are intialized, so do
        # this first
        if constraint_loads is None:
            self._constraint_loads = []
        else:
            self.constraint_loads = constraint_loads  # calls parse_outputs again

        # TODO : What if user adds symbols after constructing a System?
        # TODO : For large equations of motion, these two methods can take
        # unncessary time. One option is to let the user optionally pass in
        # these lists of symbols or derive them from the constants and
        # specifieds dictionaries.
        self._constants_symbols = self._Kane_constant_symbols()
        self._specifieds_symbols = self._Kane_undefined_dynamicsymbols()

        if constants is None:
            self.constants = dict()
        else:
            self.constants = constants

        if specifieds is None:
            self.specifieds = dict()
        else:
            self.specifieds = specifieds

        if ode_solver is None:
            self.ode_solver = odeint
        else:
            self.ode_solver = ode_solver

        if initial_conditions is None:
            self.initial_conditions = dict()
        else:
            self.initial_conditions = initial_conditions

        if times is None:
            self.times = []
        else:
            self.times = times

        self._evaluate_ode_function = None

        self._needs_code_regeneration = True
        self._last_generated_ode_user_kwargs = {}

    @property
    def coordinates(self):
        """Returns a list of the symbolic functions of time representing the
        system's generalized coordinates."""
        return self.eom_method.q[:]

    @property
    def speeds(self):
        """Returns a list of the symbolic functions of time representing the
        system's generalized speeds."""
        return self.eom_method.u[:]

    @property
    def states(self):
        """Returns a list of the symbolic functions of time representing the
        system's states, i.e. generalized coordinates plus the generalized
        speeds. These are in the same order as used in integration (as
        passed into evaluate_ode_function) and match the order of the mass
        matrix and forcing vector.

        """
        # requires settable attributes : outputs

        if self._linear_outputs_symbols:
            speeds = self.speeds + self.dummy_states
            return self.coordinates + speeds
        else:
            return self.coordinates + self.speeds

    @property
    def eom_method(self):
        """This is a sympy.physics.mechanics.KanesMethod. The method used to
        generate the equations of motion. Read-only.

        """
        return self._eom_method

    @property
    def constants(self):
        """A dict that provides the numerical values for the constants in the
        problem (all non-dynamics symbols). Keys are the symbols for the
        constants, and values are floats. Constants that are not specified in
        this dict are given a default value of 1.0.

        """
        return self._constants

    @constants.setter
    def constants(self, constants):
        self._check_constants(constants)
        self._constants = constants

    @property
    def constants_symbols(self):
        """A set of the symbolic constants (not functions of time) in the
        system.
        """
        # requires settable attributes : constants
        return self._constants_symbols

    def _check_constants(self, constants):
        symbols = self.constants_symbols
        for k in constants.keys():
            if k not in symbols:
                raise ValueError("Symbol {} is not a constant.".format(k))

    def _constants_padded_with_defaults(self):
        d = dict(zip(self.constants_symbols,
                     repeat(1.0, len(self.constants_symbols))))
        d.update(self.constants)
        return d

    @property
    def _constants_array(self):
        p_dict = self._constants_padded_with_defaults()
        return np.array([p_dict[pi] for pi in self.constants_symbols])

    @property
    def specifieds(self):
        """A dict that provides numerical values for the specified quantities
        in the problem (all dynamicsymbols that are not defined by the
        equations of motion). There are two possible formats. (1) is more
        flexible, but (2) is more efficient (by a factor of 3).

        (1) Keys are the symbols for the specified quantities, or a tuple of
        symbols, and values are the floats, arrays of floats, or functions that
        generate the values. If a dictionary value is a function, it must have
        the same signature as ``f(x, t)``, the ode right-hand-side function
        (see the documentation for the ``ode_solver`` attribute). You needn't
        provide values for all specified symbols. Those for which you do not
        give a value will default to 0.0.

        (2) There are two keys: 'symbols' and 'values'. The value for 'symbols'
        is an iterable of *all* the specified quantities in the order that you
        have provided them in 'values'. Values is an ndarray, whose length is
        `len(sys.specifieds_symbols)`, or a function of x and t that returns an
        ndarray (also of length `len(sys.specifieds_symbols)`). NOTE: You must
        provide values for all specified symbols. In this case, we do *not*
        provide default values.

        NOTE: If you switch formats with the same instance of System, you
        *must* call `generate_ode_function()` before calling `integrate()`
        again.

        Examples
        --------
        Here are examples for (1). Keys can be individual symbols, or a tuple
        of symbols. Length of a value must match the length of the
        corresponding key. Values can be functions that return iterables::

            sys = System(km)
            sys.specifieds = {(a, b, c): np.ones(3), d: lambda x, t: -3 * x[0]}
            sys.specifieds = {(a, b, c): lambda x, t: np.ones(3)}

        Here are examples for (2)::

            sys.specifieds = {'symbols': (a, b, c, d),
                              'values': np.ones(4)}
            sys.specifieds = {'symbols': (a, b, c, d),
                              'values': lambda x, t: np.ones(4)}

        """
        return self._specifieds

    @specifieds.setter
    def specifieds(self, specifieds):
        self._check_specifieds(specifieds)
        self._specifieds = specifieds

    @property
    def specifieds_symbols(self):
        """A set of the dynamicsymbols you must specify."""
        # TODO : Eventually use a method in the KanesMethod class.
        return self._specifieds_symbols

    def _assert_is_specified_symbol(self, symbol, all_symbols):
        if symbol not in all_symbols:
            raise ValueError("Symbol {} is not a 'specified' symbol.".format(
                symbol))

    def _assert_symbol_appears_multiple_times(self, symbol, symbols_so_far):
        if symbol in symbols_so_far:
            raise ValueError("Symbol {} appears more than once.".format(
                symbol))

    def _specifieds_are_in_format_2(self, specifieds):
        keys = specifieds.keys()
        if ('symbols' in keys and 'values' in keys):
            return True
        else:
            return False

    def _check_specifieds(self, specifieds):
        symbols = self.specifieds_symbols

        symbols_so_far = list()

        if self._specifieds_are_in_format_2(specifieds):

            # The symbols must be specifieds.
            for sym in specifieds['symbols']:
                self._assert_is_specified_symbol(sym, symbols)

            # Each specified symbol can appear only once.
            for sym in specifieds['symbols']:
                self._assert_symbol_appears_multiple_times(sym, symbols_so_far)
                symbols_so_far.append(sym)

            # Must have provided all specifieds.
            for sym in self.specifieds_symbols:
                if sym not in specifieds['symbols']:
                    raise ValueError(
                        "Specified symbol {} is not provided.".format(sym))

        else:

            for k, v in specifieds.items():

                # The symbols must be specifieds.
                if isinstance(k, tuple):
                    for ki in k:
                        self._assert_is_specified_symbol(ki, symbols)
                else:
                    self._assert_is_specified_symbol(k, symbols)

                # Each specified symbol can appear only once.
                if isinstance(k, tuple):
                    for ki in k:
                        self._assert_symbol_appears_multiple_times(
                            ki, symbols_so_far)
                        symbols_so_far.append(ki)
                else:
                    self._assert_symbol_appears_multiple_times(
                        k, symbols_so_far)
                    symbols_so_far.append(k)

    def _symbol_is_in_specifieds_dict(self, symbol, specifieds_dict):
        for k in specifieds_dict.keys():
            if symbol == k or (isinstance(k, tuple) and symbol in k):
                return True
        return False

    def _specifieds_padded_with_defaults(self):
        d = dict(zip(self.specifieds_symbols,
                     repeat(0.0, len(self.specifieds_symbols))))
        d.update(self.specifieds)
        return d

    @property
    def times(self):
        """A 1D ndarray of monotonic time values over which the equations of
        motion are numerically integrated. Can be set with an array-like for a
        shape(n,) array."""
        return self._times

    @times.setter
    def times(self, new_times):
        times = np.asarray(new_times)
        self._check_times(times)
        self._times = np.asarray(new_times)

    def _check_times(self, times):

        # TODO : this check can probably be removed.
        if len(times.shape) == 0:
            raise TypeError("Times supplied should be in an array_like format.")

        if not np.all(times >= 0):
            raise ValueError("Times supplied must have positive values.")

        if not np.all(np.diff(times) >= 0):
            raise ValueError("Times supplied should be in an ascending order.")

        return True

    @property
    def ode_solver(self):
        """A function that performs forward integration. It must have the same
        signature as :py:func:`scipy.integrate.odeint`, which is::

            x_history = ode_solver(f, x0, t, args=f_args)

        where ``f`` is a function ``f(x, t, *f_args)``, ``x0`` are the initial
        conditions, ``x_history`` is the state time history, ``x`` is the
        state, ``t`` is the time, and ``args`` is a keyword argument takes
        arguments that are then passed to ``f``. The default solver is
        ``odeint``.

        Examples
        ========

        SciPy introduced a unified :py:func:`scipy.integrate.solve_ivp` API
        which can be used with PyDy. ``solve_ivp`` requires a function that has
        swapped first arguments and it returns a solution object where the
        trajectory is the transpose of what ``odeint`` outputs. You can make a
        custom ODE solver function to use ``solve_ivp`` like so:

        >>> from pydy.models import multi_mass_spring_damper
        >>> sys = multi_mass_spring_damper()
        >>> sys.initial_conditions[sys.coordinates[0]] = 1.0
        >>> sys.times = [1.0, 2.0, 3.0]
        >>> from scipy.integrate import solve_ivp
        >>> def custom_ode_solver(f, x0, ts, args=(), **kwargs):
        ...     return solve_ivp(lambda t, x: f(x, t, *args), ts[[0, -1]], x0,
        ...                      t_eval=ts, **kwargs).y.T
        >>> sys.ode_solver = custom_ode_solver

        This then allows one to easiliy change methods and settings following
        SciPy's API:

        >>> sys.integrate(method='LSODA', rtol=1e-10)
        array([[ 1.00000000e+00, -5.67952532e-17],
               [ 6.59700039e-01, -5.33506568e-01],
               [ 1.50574778e-01, -4.19279930e-01]])
        >>> sys.integrate(method='RK23', rtol=1e-12)
        array([[ 1.        ,  0.        ],
               [ 0.65970115, -0.53350624],
               [ 0.15057689, -0.41928088]])

        """
        return self._ode_solver

    @ode_solver.setter
    def ode_solver(self, ode_solver):
        if not hasattr(ode_solver, '__call__'):
            msg = "``ode_solver`` ({}) is not a function."
            raise ValueError(msg.format(ode_solver))
        self._ode_solver = ode_solver

    @property
    def initial_conditions(self):
        """Initial conditions for all states (coordinates and speeds). Keys
        are the symbols for the coordinates and speeds, and values are
        floats. Coordinates or speeds that are not specified in this dict
        are given a default value of 0.0.

        """
        return self._initial_conditions

    @initial_conditions.setter
    def initial_conditions(self, initial_conditions):
        self._check_initial_conditions(initial_conditions)
        self._initial_conditions = initial_conditions

    def _check_initial_conditions(self, initial_conditions):
        symbols = self.states
        for k in initial_conditions.keys():
            if k not in symbols:
                raise ValueError("Symbol {} is not a state.".format(k))

    def _initial_conditions_padded_with_defaults(self):
        d = dict(zip(self.states,
                     repeat(0.0, len(self.states))))
        d.update(self.initial_conditions)
        return d

    @property
    def _initial_conditions_array(self):
        x0_dict = self._initial_conditions_padded_with_defaults()
        return np.array([x0_dict[xi] for xi in self.states])

    @property
    def outputs(self):
        """Dictionary of functions of time or utple of functions of time mapped
        to SymPy expressions or iterables of expressions that represent extra
        functions of the state that should be evaluted alongside the ordinary
        differential equations. Acceptable key pairs for this dictionary take
        the following three forms:

        A single function of time mapped to a function of the state::

            outputs[p(t)] = k*x(t)

        A tuple of functions of time mapped to to functions of the state::

            outputs[(f1(t), f2(t))] = (k*x(t), c*v(t))

        A tuple of functions of time mapped to a system of linear equations in
        the functions and the time derivatives of the states::

            outputs[(m1(t), m2(t))] = (m1(t) - 4*m2(t) + k*x(t).diff(t) + 2,
                                       m1(t) + 3*m2(t) - theta.diff(t))

        If equations of the last form are provided, this linear system will be
        numerically solved alongside the ordinary differential equations.

        Notes
        =====

        If your system has configuration or motion constraints, these will
        automatically be added to the outputs dictionary. If your system has
        noncontributing forces exposed and you provide names for those forces,
        these will automatically be added to the outputs dictionary.

        """
        return self._outputs

    @outputs.setter
    def outputs(self, outputs):
        self._outputs = outputs
        self._parse_outputs()
        self._needs_code_regeneration = True

    def _parse_outputs(self):
        # Divide the equations into three types:
        #
        # 1. Equations that are simply a function of the state:
        # [y1] = [f1(t, x, r, p)]
        # [y2]   [f2(t, x, r, p)]
        # [y3]   [f3(t, x, r, p)]
        #
        # 2. Equations that are linear in the state derivatives and the new
        # output variables, for example noncontributing forces. The essential
        # equations of motion will be augmented with these equations and the
        # outputs y will be solved for along in the inversion of the new
        # augmented mass matrix.
        # [Md  0] [u'] = [Fd]
        # [Mu My] [y ]   [Fy]
        #
        # TODO : How could we also support Lagrange multipliers?:
        # [Md CT] [u'] = [Fd]
        # [C  Ml] [l ]   [Fl]
        #
        # TODO : support nonlinear functions of x' later.
        # 3. Equations that are nonlinear functions of the state and its state
        # derivatives.
        # [y1] = [f1(t, x', x, r, p)]
        # [y2]   [f2(t, x', x, r, p)]
        # [y3]   [f3(t, x', x, r, p)]
        #
        # The end goal is to have something like:
        #
        # def rhs(t, x, r, p):
        #     M, F, Y1 = eval_eqs(t, x, r, p)
        #     sol = solve(M, F)
        #     xdot = sol[:len(x)]
        #     Y2 = sol[len(x):]
        #     Y3 = eval_eqs2(t, x, xdot, r, p)
        #     Y = put_in_order((Y1, Y2, Y3))
        #     return xdot, Y
        #
        # We should retain the order of the outputs in the provided dictionary
        # for Y.
        #
        # The dictionary is iterated and for all ouputs:
        # Y = [y1, y2, y3, y4, y5, y6]
        # it is separated into simple outputs:
        # Y1 = [y1, y2, y3, y6]
        # and linear outputs:
        # Y2 = [y4, y5]
        # so to reconstruct Y in the order of the dictionary we need to store
        # the indices in Y so we can do:
        # Y[simple_idxs] = Y1
        # Y[linear_idxs] = Y2

        output_names_in_order = []

        funcs_of_x = []
        solved_eq_names = []

        funcs_of_xdot = []
        linear_eq_names = []

        for var, expr in self.outputs.items():
            if isinstance(var, tuple):
                for v, e in zip(var, expr):
                    output_names_in_order.append(v)
                    if e.has(sm.Derivative):
                        funcs_of_xdot.append(e)
                        linear_eq_names.append(v)
                    else:
                        funcs_of_x.append(e)
                        solved_eq_names.append(v)
            else:
                output_names_in_order.append(var)
                if expr.has(sm.Derivative):
                    funcs_of_xdot.append(expr)
                    linear_eq_names.append(var)
                else:
                    funcs_of_x.append(expr)
                    solved_eq_names.append(var)

        if len(set(output_names_in_order)) < len(output_names_in_order):
            raise ValueError('Duplicate names found in outputs keys.')

        self._num_simple_outputs = len(solved_eq_names)
        self._simple_outputs_symbols = solved_eq_names
        if funcs_of_x:
            self._simple_outputs_matrix = sm.Matrix(funcs_of_x)
        else:
            self._simple_outputs_matrix = funcs_of_x

        if funcs_of_xdot:
            funcs_of_xdot = sm.Matrix(funcs_of_xdot)
            xd = [ui.diff() for ui in self.speeds] + linear_eq_names
            # TOOD : linear_eq_to_matrix is ideal here but doesn't function in
            # oldest support sympy.
            mass_matrix_rows = funcs_of_xdot.jacobian(xd)
            forcing_rows = -funcs_of_xdot.xreplace({xdi: 0 for xdi in xd})
        else:
            mass_matrix_rows = sm.Matrix([])
            forcing_rows = sm.Matrix([])

        self.dummy_states = [sm.Symbol('∫ ' + s.name + ' dt') for s in
                             linear_eq_names]
        self._num_linear_outputs = len(linear_eq_names)
        self._linear_outputs_symbols = linear_eq_names
        self._linear_outputs_mass_matrix_rows = mass_matrix_rows
        self._linear_outputs_forcing_rows = forcing_rows

        if self.constraints:
            self._constraint_idxs = [self._simple_outputs_symbols.index(ci)
                                     for ci in self._con_syms]

        self.outputs_symbols = output_names_in_order
        self.num_outputs = len(output_names_in_order)

        self._simple_idxs = [output_names_in_order.index(si)
                             for si in solved_eq_names]
        self._linear_idxs = [output_names_in_order.index(si)
                             for si in linear_eq_names]

    def _augment_dynamical_diff_eqs(self):
        # [Md  0] [u'] = [Fd]  <- dynamical differential equations
        # [Mu My] [y ]   [Fy]  <- extra outputs y that are linear in u'
        Md = self.eom_method.mass_matrix
        Fd = self.eom_method.forcing
        MuMy = self._linear_outputs_mass_matrix_rows
        Fy = self._linear_outputs_forcing_rows
        Mz = sm.zeros(Md.shape[0], MuMy.shape[0])
        M = Md.row_join(Mz).col_join(MuMy)
        F = Fd.col_join(Fy)
        return M, F

    @property
    def constraint_loads(self):
        """List of symbolic functions of time representing the constraint loads
        (forces & torques) associated with noncontributing loads."""
        return self._constraint_loads

    @constraint_loads.setter
    def constraint_loads(self, constraint_loads):
        if not hasattr(self.eom_method, 'auxiliary_eqs'):
            msg = ('The KanesMethod object has no auxiliary equations and '
                   'thus constraint loads cannot be provided.')
            raise RuntimeError(msg)
        if len(constraint_loads) != len(self.eom_method._uaux):
            msg = ('You must provide symbols for {} constraint loads that '
                    'are present in the auxiliary equations.')
            raise ValueError(msg.format(len(self.eom_method._uaux)))
        # TODO : Check that the constraint loads are present in the auxiliary
        # equations and that they are not a coordinate, speed, or specified.
        # dj/dt = j' = constraint_load
        # create some dummy impulse states
        self._constraint_loads = list(constraint_loads)
        if tuple(constraint_loads) in self.outputs:
            raise ValueError('Constraint loads already present in outputs.')

        self.outputs[tuple(constraint_loads)] = self.eom_method.auxiliary_eqs
        self._parse_outputs()
        self._needs_code_regeneration = True

    @property
    def evaluate_ode_function(self):
        """A function generated by ``generate_ode_function`` that computes
        the state derivatives:

            x' = evaluate_ode_function(x, t, *args)

        This function is used by the ``ode_solver``.

        """
        # TODO : It would be more useful if the generated docstring was shown
        # when help(system.evaluate_ode_function) is called.
        return self._evaluate_ode_function

    def _args_for_gen_ode_func(self):
        """Returns a tuple of arguments in the form required by
        ``pydy.codegen.ode_function_generators.generate_ode_function``.

        """
        if self._linear_outputs_symbols:
            Fd = self.eom_method.forcing
            Fa = self._linear_outputs_forcing_rows
            # [Md  0] [u'] = [Fd]
            # [Mu Mj] [j']   [Fa]
            forcing = Fd.col_join(Fa)
            speeds = self.speeds + self.dummy_states
        else:
            forcing = self.eom_method.forcing
            speeds = self.speeds

        args = (forcing,
                self.coordinates,
                speeds,
                self.constants_symbols)

        return args

    def _kwargs_for_gen_ode_func(self):
        """Returns a dictionary of arguments in the form required by
        ``pydy.codegen.ode_function_generators.generage_ode_function``.

        """
        if self._specifieds_are_in_format_2(self.specifieds):
            specifieds = self.specifieds['symbols']
        else:
            specifieds = self.specifieds_symbols

        # generate_ode_func does not accept an empty tuple for the
        # specifieds, so set it to None
        if not specifieds:
            specifieds = None

        kin_diff_dict = self.eom_method.kindiffdict()
        kin_diff_rhs = sm.Matrix([kin_diff_dict[q.diff()] for q in
                                  self.coordinates])

        if self._linear_outputs_symbols:
            Md = self.eom_method.mass_matrix
            MuMj = self._linear_outputs_mass_matrix_rows
            Mz = sm.zeros(Md.shape[0], MuMj.shape[0])
            # [Md  0] [u'] = [Fd]
            # [Mu Mj] [j']   [Fa]
            mass_matrix = Md.row_join(Mz).col_join(MuMj)
        else:
            mass_matrix = self.eom_method.mass_matrix

        kwargs = {
            'mass_matrix': mass_matrix,
            'coordinate_derivatives': kin_diff_rhs,
            'specifieds': specifieds,
        }

        if self._simple_outputs_symbols:
            kwargs['outputs'] = self._simple_outputs_matrix

        return kwargs

    def _extract_constraints(self):
        """Extracts the configuration and motion constraints from the
        eom_method and stores them in attributes."""

        if self.eom_method._f_h:
            self.config_constraints = self.eom_method._f_h
        else:
            self.config_constraints = sm.Matrix([])

        if self.eom_method._k_nh:
            # rebuild the nonholonomic constraints from KanesMethod
            # TODO : KanesMethod and _Method should store the original
            # constraints passed by the user. Fix in sympy.physics.mechanics!
            self.motion_constraints = (
                self.eom_method._k_nh*self.eom_method.u +
                self.eom_method._f_nh)
        else:
            self.motion_constraints = sm.Matrix([])

        self.num_config_constraints = len(self.config_constraints)
        self.num_motion_constraints = len(self.motion_constraints)

    @property
    def constraints(self):
        """A column matrix of configuration and motion constraints expressions,
        ordered as stored in KanesMethod."""
        constraints = sm.Matrix([])

        if self.config_constraints or self.motion_constraints:
            if self.config_constraints:
                constraints = self.config_constraints

            if constraints and self.motion_constraints:
                constraints = constraints.col_join(self.motion_constraints)
            else:
                constraints = self.motion_constraints

        return constraints

    @property
    def num_constraints(self):
        """Total number of configuration and motion constaints."""
        return self.num_config_constraints + self.num_motion_constraints

    def set_dependent_initial_conditions(self, dep_vars=None, use_jac=False,
                                         **root_kwargs):
        """Sets the initial conditions of the dependent coordinates and
        dependent speeds using the holonomic and nonholonomic constraints,
        respectively.

        Parameters
        ==========
        dep_vars : iterable of Function()(t), optional
            Dependent coordinates and speeds to solve for. The number of
            coordinates should be equal to the number of holonomic constraints.
            The number of speeds should be equal to the number of nonholonic
            constraints. If None, the dependent coordinates and speeds are
            those used in KanesMethod instantiation.
        use_jac : boolean, optional
            If true the Jacobian of the constraint equations will be used to
            solve the constraint equations for the dependent states.
        root_kwargs
            Extra keyword arguments that are passed to ``scipy.optimize.root``.

        """
        # TODO : The nonholonomic constraints can be solved analytically, root
        # is only required for the coordinates.

        if self.num_constraints == 0:
            msg = ('This system does not have constraints, set all initial '
                   'conditions yourself.')
            raise ValueError(msg)
        else:
            num_holo = self.num_config_constraints
            num_nonh = self.num_motion_constraints

        # TODO : These variables should be publicly accessible on KanesMethod.
        if dep_vars is None:
            dep_vars = self.eom_method._qdep[:] + self.eom_method._udep[:]

        # TODO : Would be nice to check if the dependent variables are present
        # in the constraints and that the right number of coordinates and
        # speeds are each supplied.
        if len(dep_vars) != self.num_constraints:
            msg = (f'You must supply {num_holo} dependent coordinates and '
                   f'{num_nonh} dependent speeds.')
            raise ValueError(msg)

        x = self._initial_conditions_array
        p = self._constants_array
        dep_guess = [self.initial_conditions[xi] for xi in dep_vars]
        dep_idxs = [self.states.index(xi) for xi in dep_vars]

        if use_jac:
            jac = self.constraints.jacobian(dep_vars)
            eval_jac = sm.lambdify((self.states, self.constants_symbols), jac,
                                   cse=True)

            def eval_f(x_dep, p):
                x[dep_idxs] = x_dep
                return self.evaluate_constraints(x=x), eval_jac(x, p)

            fprime = True
        else:

            def eval_f(x_dep, p):
                x[dep_idxs] = x_dep
                return self.evaluate_constraints(x=x)

            fprime = False

        # not settable by the user, so overwrite:
        root_kwargs['args'] = (p, )
        root_kwargs['jac'] = fprime

        sol = root(eval_f, dep_guess, **root_kwargs)
        if not sol.success:
            msg = ('Failed to find a solution that meets tolerance. Maybe a '
                   'better guess will help or you may have to manually solve '
                   'for the dependent coordinates. SciPy root() failure '
                   'message: ' + sol.message)
            warnings.warn(msg, PyDyUserWarning, stacklevel=2)

        dep_vals = sol.x

        for si, vi in zip(dep_vars, dep_vals):
            self.initial_conditions[si] = vi

    def generate_ode_function(self, **kwargs):
        """Calls :py:func:`~generate_ode_function` with the appropriate
        arguments, and sets the ``evaluate_ode_function`` attribute to the
        resulting function.

        Parameters
        ----------
        kwargs
            All other kwargs are passed onto
            ``pydy.codegen.ode_function_generators.generate_ode_function()``.
            Don't specify the ``specifieds`` keyword argument though; the
            ``System`` class takes care of those.

        Returns
        -------
        evaluate_ode_function : function
            A function which evaluates the derivaties of the states.

        Notes
        -----

        If the Cython generator is selected and you have a custom
        ``ode_solver`` set, keyword argument ``force_c_contiguous`` will be
        automatically set to ``True``. You can disable this by setting it to
        ``False`` but you must ensure ensure that ode solver only passes C
        contiguous arrays to the generated ode function. Forcing C contiguous
        arrays introduces a small performance penalty due to the necessity of
        copying arrays.

        """
        self._parse_outputs()  # call before generating the args/kwargs below

        self._last_generated_ode_user_kwargs = kwargs.copy()

        if 'specified' in kwargs:
            kwargs.pop('specified')
            print("User supplied 'specified' kwarg was disregarded.")

        if 'specifieds' in kwargs:
            kwargs.pop('specifieds')
            print("User supplied 'specifieds' kwarg was disregarded.")

        kwargs.update(self._kwargs_for_gen_ode_func())

        if 'force_c_contiguous' not in kwargs:
            if 'generator' in kwargs:
                if (kwargs['generator'] == 'cython' and self.ode_solver is
                        odeint):
                    kwargs['force_c_contiguous'] = False
                # NOTE : This ensures that the arrays are forced to be C
                # contiguous if a user sets an integrator other than odeint,
                # but only if they are using the Cython generator. There is a
                # performance penalty but it avoids errors and they can
                # manually override this to False if they know their ode_solver
                # choice delivers C contiguous arrays.
                elif kwargs['generator'] == 'cython':
                    kwargs['force_c_contiguous'] = True

        self._evaluate_ode_function = generate_ode_function(
            *self._args_for_gen_ode_func(),
            **kwargs)

        return self.evaluate_ode_function

    def _prep_for_evaluate(self):
        # Users might have changed these properties by directly accessing the
        # dict, without using the setter. Before we integrate, make sure they
        # did not muck up these dicts.
        self._check_constants(self.constants)
        self._check_specifieds(self.specifieds)
        self._check_initial_conditions(self.initial_conditions)
        self._check_times(self.times)

        if self.evaluate_ode_function is None or self._needs_code_regeneration:
            self.generate_ode_function(**self._last_generated_ode_user_kwargs)
            self._needs_code_regeneration = False

        if self._specifieds_are_in_format_2(self.specifieds):
            specified_value = self.specifieds['values']
        else:
            specified_value = self._specifieds_padded_with_defaults()

        # If there are no specifieds then specified_value will be an empty
        # dict.
        if isinstance(specified_value, dict) and not specified_value:
            args = (self._constants_padded_with_defaults(),)
        else:
            args = (specified_value, self._constants_padded_with_defaults())

        return self._initial_conditions_array, args

    def _prep_x_t_overrides(self, x, t):
        x_default, args = self._prep_for_evaluate()

        if x is None:
            x = x_default
        x = np.asarray(x)

        # TODO : Maybe the default value of .times shouldn't be empty. Think
        # about this.
        if t is None:
            if len(x.shape) == 1:
                if self.times.size == 0:
                    t = 0.0
                else:
                    t = self.times[0]
            else:
                if self.times.size == 0:
                    t = np.zeros(x.shape[0])
                else:
                    t = self.times

        return x, t, args

    def evaluate_ode(self, x=None, t=None):
        """Returns the right hand side of the differential equations. The
        default is to evaluate at the set initial_conditions at the first time
        value. Pass in optional arguments to override using the initial state
        and time.

        Parameters
        ==========
        x : array_like, shape(n,) or shape(m, n), optional
            State values at time t.
        t : float or array_like, shape(m,), optional
            Time or m time values.

        Returns
        =======
        x_dot : ndarray, shape(n,) or shape(m, n)
           Time derivative of the states at time t.

        Notes
        =====

        This method is present for convenience, it is not designed to be used
        where performance matters, use :py:meth:`System.evaluate_ode_function`
        directly when performance is needed.

        To see the order of the state values use::

            >>> system = System(...)
            >>> system.states

        or::

            >>> rhs = system.generate_ode_function()
            >>> help(rhs)

        """
        x, t, args = self._prep_x_t_overrides(x, t)

        if len(x.shape) == 1 and not isinstance(t, float):
            raise ValueError('Time must be a float.')
        elif len(x.shape) == 1 and isinstance(t, float):
            res = self.evaluate_ode_function(x, t, *args)
            if self._simple_outputs_symbols:
                return res[0]
            else:
                return res
        # NOTE : I tried to make use of numpy.vectorize but it is not possible
        # due to args not being necessarily being comprised of arrays.
        elif len(x.shape) == 2:
            if x.shape[0] != len(t):
                raise ValueError('x trajectory must have same length as t.')
            xdot = np.zeros_like(x)
            for i, (ti, xi) in enumerate(zip(t, x)):
                res = self.evaluate_ode_function(xi, ti, *args)
                if self._simple_outputs_symbols:
                    xdot[i, :] = res[0]
                else:
                    xdot[i, :] = res
            return xdot

    def evaluate_outputs(self, x=None, t=None):
        """Returns an array of the evaluated outputs. The default is to
        evaluate at the initial conditions at the first time value. Pass in
        optional arguments to override the state or time.

        Parameters
        ==========
        x : array_like, shape(n,) or shape(m, n), optional
            State values at time t.
        t : float or array_like, shape(m,), optional
            Time or m time values.

        Returns
        =======
        y : ndarray, shape(o,) or shape(m, o)
           Output values at time t.

        Notes
        =====

        This method is present for convenience, it is not designed to be used
        where performance matters, use :py:meth:`System.evaluate_ode_function`
        directly when performance is needed.

        To see the order of the state values use::

            >>> system = System(...)
            >>> system.states

        or::

            >>> rhs = system.generate_ode_function()
            >>> help(rhs)

        """
        if not self.outputs:
            raise ValueError('This system has no outputs.')

        x, t, args = self._prep_x_t_overrides(x, t)

        if len(x.shape) == 1 and not isinstance(t, float):
            raise ValueError('Time must be a float.')
        elif len(x.shape) == 1 and isinstance(t, float):
            if self._linear_outputs_symbols and self._simple_outputs_symbols:
                y = np.zeros(self.num_outputs)
                xdot, y1 = self.evaluate_ode_function(x, t, *args)
                y[self._simple_idxs] = y1
                y[self._linear_idxs] = xdot[-len(self.dummy_states):]
                return y
            elif (self._linear_outputs_symbols and not
                  self._simple_outputs_symbols):
                xdot = self.evaluate_ode_function(x, t, *args)
                return xdot[-len(self.dummy_states):]
            else:
                return self.evaluate_ode_function(x, t, *args)[1]
        # NOTE : I tried to make use of numpy.vectorize but it is not possible
        # due to args not being necessarily being comprised of arrays.
        elif len(x.shape) == 2:
            if x.shape[0] != len(t):
                raise ValueError('x trajectory must have same length as t.')
            y = np.zeros((len(t), self.num_outputs))
            for i, (ti, xi) in enumerate(zip(t, x)):
                if (self._linear_outputs_symbols and
                    self._simple_outputs_symbols):
                    xdot, y1 = self.evaluate_ode_function(xi, ti, *args)
                    y[i, self._simple_idxs] = y1
                    y[i, self._linear_idxs] = xdot[-len(self.dummy_states):]
                elif (self._linear_outputs_symbols and not
                      self._simple_outputs_symbols):
                    xdot = self.evaluate_ode_function(xi, ti, *args)
                    y[i, :] = xdot[-len(self.dummy_states):]
                else:
                    y[i, :] = self.evaluate_ode_function(xi, ti, *args)[1]
            return y

    def evaluate_constraints(self, x=None, t=None):
        """Returns the values of the configuration and motion constraints at
        the initial condition or, alternatively, for the provided state vector.

        Parameters
        ==========
        x : array_like, shape(n,) or shape(m, n), optional
            State vector of n states or a series of m state vectors.
        t : float or array_like, shape(m,), optional
            Time or m time values.

        Returns
        =======
        ndarray, shape(o,) or shape(m, o)
            Constraint vector of o constraints or a series of m constraint
            vectors.

        Notes
        =====

        To see the order of the state values use::

            >>> system = System(...)
            >>> system.states

        or::

            >>> rhs = system.generate_ode_function()
            >>> help(rhs)

        """
        if self.num_constraints == 0:
            raise ValueError('This system has no constraints.')

        x, t, args = self._prep_x_t_overrides(x, t)

        if len(x.shape) == 1 and not isinstance(t, float):
            raise ValueError('Time must be a float.')
        elif len(x.shape) == 1 and isinstance(t, float):
            y = self.evaluate_ode_function(x, t, *args)[1]
            return y[self._constraint_idxs]
        elif len(x.shape) == 2:
            if x.shape[0] != len(t):
                raise ValueError('x trajectory must have same length as t.')
            con = np.zeros((x.shape[0], self.num_constraints))
            for i, (ti, xi) in enumerate(zip(t, x)):
                y = self.evaluate_ode_function(xi, ti, *args)[1]
                con[i, :] = y[self._constraint_idxs]
            return con

    def evaluate_config_constraints(self, x=None, t=None):
        """Returns the values of the configuration at the initial condition or,
        alternatively, for the provided state vector.

        Parameters
        ==========
        x : array_like, shape(n,) or shape(m, n), optional
            State vector of n states or a series of m state vectors.
        t : float or array_like, shape(m,), optional
            Time or m time values.

        Returns
        =======
        ndarray, shape(o,) or shape(m, o)
            Constraint vector of o constraints or a series of m constraint
            vectors.

        Notes
        =====

        To see the order of the state values use::

            >>> system = System(...)
            >>> system.states

        or::

            >>> rhs = system.generate_ode_function()
            >>> help(rhs)

        """
        if self.num_config_constraints == 0:
            raise ValueError('This system has no configuration constraints.')
        con = self.evaluate_constraints(x=x, t=t)
        if len(con.shape) == 1:
            return con[:self.num_config_constraints]
        else:
            return con[:, :self.num_config_constraints]

    def evaluate_motion_constraints(self, x=None, t=None):
        """Returns the values of the motion at the initial condition or,
        alternatively, for the provided state vector.

        Parameters
        ==========
        x : array_like, shape(n,) or shape(m, n), optional
            State vector of n states or a series of m state vectors.
        t : float or array_like, shape(m,), optional
            Time or m time values.

        Returns
        =======
        ndarray, shape(o,) or shape(m, o)
            Constraint vector of o constraints or a series of m constraint
            vectors.

        Notes
        =====

        To see the order of the state values use::

            >>> system = System(...)
            >>> system.states

        or::

            >>> rhs = system.generate_ode_function()
            >>> help(rhs)

        """
        if self.num_motion_constraints == 0:
            raise ValueError('This system has no motion constraints.')
        con = self.evaluate_constraints(x=x, t=t)
        if len(con.shape) == 1:
            return con[self.num_config_constraints:]
        else:
            return con[:, self.num_config_constraints:]

    def integrate(self, **solver_kwargs):
        """Integrates the equations ``evaluate_ode_function()`` using
        ``ode_solver``.

        It is necessary to have first generated an ode function. If you have
        not done so, we do so automatically by invoking
        ``generate_ode_function()``. However, if you want to customize how this
        function is generated (e.g., change the generator to cython), you can
        call ``generate_ode_function()`` on your own (before calling
        ``integrate()``).

        Parameters
        ----------
        **solver_kwargs
            Optional arguments that are passed on to the ``ode_solver``.

        Returns
        -------
        x_history : ndarray, shape(num_integrator_time_steps, 2)
            The trajectory of states (coordinates and speeds) through the
            requested time interval. num_integrator_time_steps is either
            len(times) if len(times) > 2, or is determined by the
            ``ode_solver``.

        """
        if len(self.times) < 2:
            raise ValueError('The times vector must be at least length 2.')

        x0, args = self._prep_for_evaluate()

        # NOTE : User cannot pass in args, System handles that.
        solver_kwargs.pop('args', None)


        # NOTE : skip the outputs if present
        def func(*args, **kwargs):
            return self.evaluate_ode_function(*args, **kwargs)[0]


        x_history = self.ode_solver(
            (func if self._simple_outputs_symbols else
             self.evaluate_ode_function),
            x0,
            self.times,
            args=args, **solver_kwargs)

        return x_history

    def _Kane_inlist_insyms(self):
        """TODO temporary."""

        uaux = self.eom_method._uaux[:]
        uauxdot = [sm.diff(i, dynamicsymbols._t) for i in uaux]

        # Checking for dynamic symbols outside the dynamic differential
        # equations; throws error if there is.

        # TODO : KanesMethod should provide public attributes for qdot,
        # udot, uaux, and uauxdot.
        # NOTE : There have been breaking changes in SymPy, e.g.
        # https://github.com/pydy/pydy/issues/395, so all of these need to be
        # converted to lists before summing.
        insyms = set(list(self.eom_method.q[:]) +
                     list(self.eom_method._qdot[:]) +
                     list(self.eom_method.u[:]) +
                     list(self.eom_method._udot[:]) +
                     list(uaux) + list(uauxdot))

        inlist = (self.eom_method.forcing_full[:] +
                  self.eom_method.mass_matrix_full[:])

        return inlist, insyms

    def _Kane_undefined_dynamicsymbols(self):
        """Similar to ``_find_dynamicsymbols()``, except that it checks all
        syms used in the system. Code is copied from ``linearize()``.

        TODO temporarily here until KanesMethod and Lagranges method have an
        interface for obtaining these quantities.

        """
        from_eoms, from_sym_lists = self._Kane_inlist_insyms()
        functions_of_time = set()
        for expr in from_eoms:
            functions_of_time = functions_of_time.union(
                find_dynamicsymbols(expr))
        return functions_of_time.difference(from_sym_lists)

    def _Kane_constant_symbols(self):
        """Similar to ``_find_othersymbols()``, except it checks all syms used in
        the system.

        Remove the time symbol.

        TODO temporary.

        """
        from_eoms, from_sym_lists = self._Kane_inlist_insyms()
        unique_symbols = set()
        for expr in from_eoms:
            unique_symbols = unique_symbols.union(expr.free_symbols)
        constants = unique_symbols
        constants.remove(dynamicsymbols._t)
        return constants
