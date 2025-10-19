import timeit

import numpy as np
import sympy as sm
import scipy as sp

from pydy import models
from pydy.codegen.ode_function_generators import (
    generate_ode_function, CythonODEFunctionGenerator)


sys = models.n_link_pendulum_on_cart(8, True, True)
constants = list(sm.ordered(sys.constants_symbols))
specifieds = list(sm.ordered(sys.specifieds_symbols))

p_array = np.random.random(len(constants))
r_array = np.random.random(len(specifieds))
x = np.random.random(len(sys.states))

itr = 10000

rhs_sym = generate_ode_function(
    sys.eom_method.forcing,
    sys.coordinates,
    sys.speeds,
    constants,
    specifieds=specifieds,
    mass_matrix=sys.eom_method.mass_matrix,
    coordinate_derivatives=sm.Matrix(sys.speeds),
    linear_sys_solver='sympy',
    constants_arg_type='array',
    specifieds_arg_type='array',
    generator='cython')

time = timeit.timeit(lambda: rhs_sym(x, 0.0, r_array, p_array), number=itr)
print('SymPy LU decomposition symbolic solve time: ', time)

g = CythonODEFunctionGenerator(
    sys.eom_method.forcing,
    sys.coordinates,
    sys.speeds,
    constants,
    specifieds=specifieds,
    mass_matrix=sys.eom_method.mass_matrix,
    coordinate_derivatives=sm.Matrix(sys.speeds),
    linear_sys_solver='numpy',
    constants_arg_type='array',
    specifieds_arg_type='array',
)
rhs_num = g.generate()

time = timeit.timeit(lambda: rhs_num(x, 0.0, r_array, p_array), number=itr)
print('numpy.linalg.solve time: ', time)


def numpy_umath_linalg_solve(A, b):
    return np.squeeze(np.linalg._umath_linalg.solve(A, np.atleast_2d(b).T))


g.linear_sys_solver = numpy_umath_linalg_solve
time = timeit.timeit(lambda: rhs_num(x, 0.0, r_array, p_array), number=itr)
print('numpy.linalg._umath_linalg.solve time: ', time)


def scipy_linalg_solve(A, b):
    return sp.linalg.solve(A, b, lower=True, assume_a='pos',
                           check_finite=False, overwrite_a=True,
                           overwrite_b=True)


g.linear_sys_solver = scipy_linalg_solve
time = timeit.timeit(lambda: rhs_num(x, 0.0, r_array, p_array), number=itr)
print('scipy.linalg.solve (Cholesky) time: ', time)


def scipy_linalg_lapack_dposv(A, b):
    sp.linalg.lapack.dposv(A, b, lower=True, overwrite_a=True,
                           overwrite_b=True)


g.linear_sys_solver = scipy_linalg_lapack_dposv
time = timeit.timeit(lambda: rhs_num(x, 0.0, r_array, p_array), number=itr)
print('scipy.linalg.lapack.dposv time: ', time)
