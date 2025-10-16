==================
3D N-Body Pendulum
==================

.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`3d-n-body-pendulum` or Jupyter notebook:
   :jupyter-download:notebook:`3d-n-body-pendulum`.


Objectives
----------

- Show how to use ``generate_ode_function`` as an alternative to ``lambdify``
- Show how to use ``PyDy Visualization`` to generate a 3D animation


Description
-----------

A pendulum consisting of n rods of length :math:`l` and mass
:math:`n_{\textrm{link}}`.
A ball of mass :math:`m` and radius :math:`r` is attached to each rod, such
that the rod goes through the center of each ball. The center of the ball is
fixed at the middle of the rod. The ball rotates around the rod.  If two balls
collide, they are ideally elastic, with spring constant :math:`k`.
The collision is modeled using the ``sm.Heaviside(..)`` function.
The balls are ideally slick, so collisions will not affect their rotation.
There may be speed dependent friction between
the rod and the ball, with coefficient of friction :math:`\textrm{reibung}`.
A particle with mass :math:`m_1` is attached to each ball.

Notes
-----

- ``generate_ode_function`` needs C - contiguous arrays as inputs, if
  generator = 'cython' is used. If in ``solve_ivp`` method :math:`\neq`'RK45',
  it will return non C contiguous ``y``. Then ``rhs_gen`` cannot be used
  directly, but a separate function must be defined, which must contain
  ``y = np.ascontiguousarray(y)``

.. jupyter-execute::

    import sympy as sm
    import sympy.physics.mechanics as me
    import time
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.integrate import solve_ivp, odeint
    from scipy.optimize import root
    from pydy.codegen.ode_function_generators import generate_ode_function
    from pydy.viz.shapes import Cylinder, Sphere
    from pydy.viz.scene import Scene
    from pydy.viz.visualization_frame import VisualizationFrame

Number of Bodies

Number of pendulum bodies, labelled 0, 1, .., n-1, must be two or more

.. jupyter-execute::

    n = 3

    if n <= 1 or not isinstance(n, int):
        raise ValueError('n must be an integer larger than 1')

Equations of Motion, Kane's Method
==================================

.. jupyter-execute::

    start = time.time()

    m, m1, m_link, g, r, l, reibung, k, t = sm.symbols(
        'm, m1, m_link, g, r, l, reibung, k, t')
    iXX, iYY, iZZ = sm.symbols('iXX, iYY, iZZ')
    t = me.dynamicsymbols._t

    q = []     # holds the generalized coordinates of each rod
    u = []     # generalized angular speeds

    A = []     # frames of each link

    Dmc = []       # geometric center of each body
    Dmc_link = []  # mass center of each link
    P = []         # points at end of each rod. Dmc_i is between P_i and P_i+1
    punkt = []     # marks a red dot on each ball, just used for animation
    rhs_subs = []  # substitutes for the rhs, to calculate the reaction forces


Virtual speeds and reaction forces needed to get the reaction forces at the
suspension point P0.

.. jupyter-execute::

    auxx, auxy, auxz, fx, fy, fz = me.dynamicsymbols(
        'auxx, auxy, auxz, fx, fy, fz')

Define the general coordinates, the gen. speeds, the frames and the points.

.. jupyter-execute::

    for i in range(n):
        for j in ('x', 'y', 'z'):
            q.append(me.dynamicsymbols('q' + j + str(i)))
            u.append(me.dynamicsymbols('u' + j + str(i)))
            rhs_subs.append(me.dynamicsymbols('rhs_subs' + j + str(i)))

        A.append(me.ReferenceFrame('A' + str(i)))
        Dmc.append(me.Point('Dmc' + str(i)))
        Dmc_link.append(me.Point('Dmc_link' + str(i)))
        P.append(me.Point('P' + str(i)))
        punkt.append(me.Point('punkt' + str(i)))

    N = me.ReferenceFrame('N')        # inertial frame
    P0 = me.Point('P0')
    P0.set_vel(N, 0)  # fixed in inertial frame

The lists rot, rot1 are needed for the kinematical equations, see below.

.. jupyter-execute::

    rot = []
    rot1 = []

It is very important, that the angular speeds be expressed in terms of
the 'child frame', otherwise the equations of motion become very large.

Note the virtual speeds at P[0], the suspension point.

.. jupyter-execute::

    A[0].orient_body_fixed(N, (q[0], q[1], q[2]), '123')
    rot.append(A[0].ang_vel_in(N))
    A[0].set_ang_vel(N, u[0]*A[0].x + u[1]*A[0].y + u[2]*A[0].z)
    rot1.append(A[0].ang_vel_in(N))

    for i in range(1, n):
        A[i].orient_body_fixed(A[i-1], (q[3*i], q[3*i+1], q[3*i+2]), '123')
        rot.append(A[i].ang_vel_in(N))   # needed for the kin. equations below
        A[i].set_ang_vel(N, u[3*i]*A[i].x + u[3*i+1]*A[i].y + u[3*i+2]*A[i].z)
        rot1.append(A[i].ang_vel_in(N))

Locate the various points, and define their speeds.

.. jupyter-execute::

    P[0].set_pos(P0, 0.)
    P[0].set_vel(N, auxx*N.x + auxy*N.y + auxz*N.z)  # Suspension point.
    Dmc[0].set_pos(P[0], l/2. * A[0].y)
    Dmc_link[0].set_pos(P[0], l/2. * A[0].y)
    Dmc_link[0].v2pt_theory(P[0], N, A[0])
    Dmc[0].v2pt_theory(P[0], N, A[0])
    punkt[0].set_pos(Dmc[0], r*A[0].z)  # only for the red dot in the animation
    punkt[0].v2pt_theory(Dmc[0], N, A[0])

    for i in range(1, n):
        P[i].set_pos(P[i-1], l * A[i-1].y)
        P[i].v2pt_theory(P[i-1], N, A[i-1])
        Dmc[i].set_pos(P[i], l/sm.S(2.) * A[i].y)
        Dmc[i].v2pt_theory(P[i], N, A[i])
        Dmc_link[i].set_pos(P[i], l/sm.S(2.) * A[i].y)
        Dmc_link[i].v2pt_theory(P[i], N, A[i])
        punkt[i].set_pos(Dmc[i], r*A[i].z)
        punkt[i].v2pt_theory(Dmc[i], N, A[i])

Make the list of the bodies.

.. jupyter-execute::

    balls = []
    points = []
    links = []
    for i in range(n):
        Inert = me.inertia(A[i], iXX, iYY, iZZ)
        balls.append(me.RigidBody('body' + str(i), Dmc[i], A[i], m,
                              (Inert, Dmc[i])))
        # the red dot may have a mass
        points.append(me.Particle('punct' + str(i), punkt[i], m1))
        inert_link = me.inertia(A[i], m_link*l**2/12., 0, m_link*l**2/12.)
        links.append(me.RigidBody('link' + str(i), Dmc_link[i], A[i], m_link,
                                  (inert_link, Dmc_link[i])))
    BODY = balls + points + links


Set up the forces.

There are:
- potential gravitational forces forces
- potential spring forces when the balls penetrate each other
- dissipative friction forces due to friction between rods and balls

Note how the reaction forces at P[0] are set up.

.. jupyter-execute::

    FG = ([(Dmc[i], -m*g*N.y) for i in range(n)] +
          [(punkt[i], -m1*g*N.y)for i in range(n)] +
          [(Dmc_link[i], -m_link*g*N.y) for i in range(n)])

    FB = [(P[0], fx*N.x + fy*N.y + fz*N.z)]
    for i in range(n):
        for j in range(i+1, n):
            aa = Dmc[j].pos_from(Dmc[i])
            bb = aa.magnitude()
            aa = aa.normalize()
            forceij = (Dmc[j],  k * (2 * r - bb) * aa *
                       sm.Heaviside(2 * r - bb))
            FB.append(forceij)
            forceji = (Dmc[i], -k * (2 * r - bb) * aa *
                       sm.Heaviside(2 * r - bb))
            FB.append(forceji)

        friction_i = (A[i], -reibung * u[3*i + 1] * A[i].y)  # around A[i].y
        FB.append(friction_i)

    FL = FG + FB  # list of forces

Kinematic equations.

Again it is very important that the frames A[i] be used below. Otherwise
the equations of motion become very large. Note how rot, rot1 from above
are used.

.. jupyter-execute::

    kd = []
    for i in range(n):
        for uv in A[i]:
            kd.append(me.dot(rot[i] - rot1[i], uv))

Kanes's Equations
-----------------

.. jupyter-execute::

    q1 = q
    u1 = u
    aux = [auxx, auxy, auxz]

    KM = me.KanesMethod(N, q_ind=q1, u_ind=u1, kd_eqs=kd, u_auxiliary=aux)
    fr, frstar = KM.kanes_equations(BODY, FL)

    react_forces = KM.auxiliary_eqs

The reaction forces (of course) contain accelerations.
They are replaced by ``rhs`` as place holders and will be calculated
numerically below. Symbolic calculation would be possible to too
time consuming.

.. jupyter-execute::

    react_forces = me.msubs(react_forces, {u[i].diff(t): rhs_subs[i]
                                           for i in range(len(u))})


``generate_ode_function`` needs the mass matrix and the forcing. The
kinematic differential equations to be supplied separately.

.. jupyter-execute::

    MM = KM.mass_matrix
    print('MM DS', me.find_dynamicsymbols(MM))
    print('MM free symbols', MM.free_symbols)
    print(f'MM contains {sm.count_ops(MM):,} operations', '\n')

    force = KM.forcing
    print('force DS', me.find_dynamicsymbols(force))
    print('force free symbols', force.free_symbols)
    print(f'force contains {sm.count_ops(force):,} operations', '\n')


Energy, Momentum
----------------

.. jupyter-execute::

    pot_energie = (sum(
        [m*g*me.dot(Dmc[i].pos_from(P[0]), N.y) for i in range(n)]) +
        sum([m1*g*me.dot(punkt[i].pos_from(P[0]), N.y)
        for i in range(n)]) +
        sum([m_link*g*me.dot(Dmc_link[i].pos_from(P[0]), N.y)
        for i in range(n)]))

    kin_energie = me.msubs(sum([BODY[i].kinetic_energy(N)
                                for i in range(3*n)]), {i: 0 for i in aux})
    spring_energie = sm.S(0.)
    for i in range(n):
        for j in range(i+1, n):
            aa = Dmc[j].pos_from(Dmc[i])
            bb = aa.magnitude()
            aa = aa.normalize()
            spring_energie += 0.5 * k * (2*r - bb)**2 * sm.Heaviside(2.*r - bb)

    aux_dict = {i: 0 for i in aux}
    ang_moment_x = sum([body.angular_momentum(P0, N).dot(N.x).subs(aux_dict)
                        for body in BODY])
    ang_moment_y = sum([body.angular_momentum(P0, N).dot(N.y).subs(aux_dict)
                        for body in BODY])
    ang_moment_z = sum([body.angular_momentum(P0, N).dot(N.z).subs(aux_dict)
                        for body in BODY])
    ang_momentum = [ang_moment_x, ang_moment_y, ang_moment_z]


Compilation
-----------

.. jupyter-execute::

    qL = q1 + u1
    pL = [m, m1, m_link, g, r, l, iXX, iYY, iZZ, reibung, k]

    specified = None
    constants = np.array(pL)

The kinematic equations must be solved for the derivatives of the coordinates.
The solution must be sorted so that it corresponds to the sequence in KM.q

.. jupyter-execute::

    kd = sm.Matrix(kd)
    loesung = sm.solve(kd, [q[i].diff(t) for i in range(3 * n)])

    schluessel = [i.diff(t) for i in KM.q]
    kin_eqs_solved = sm.Matrix([loesung[i] for i in schluessel])

    rhs_gen = generate_ode_function(
        force,
        KM.q,
        KM.u,
        constants=constants,
        mass_matrix=MM,
        specifieds=specified,
        coordinate_derivatives=kin_eqs_solved,  # rhs of kin. diff. equations
        generator='symjit',
        linear_sys_solver='numpy',
        constants_arg_type='array',
        specifieds_arg_type='array',
        time_first='True',
    )

Below lambdify is used as speed is of no concern.

Dmc_loc, punkt_loc are needed for the animation only.

.. jupyter-execute::

    punkt_loc = []
    Dmc_loc = []
    for i in range(n):
        punkt_loc.append([me.dot(punkt[i].pos_from(P[0]), uv) for uv in N])
        Dmc_loc.append([me.dot(Dmc[i].pos_from(P[0]), uv) for uv in N])


    pot_lam = sm.lambdify(qL + pL, pot_energie, cse=True)
    kin_lam = sm.lambdify(qL + pL, kin_energie, cse=True)
    spring_lam = sm.lambdify(qL + pL, spring_energie, cse=True)

    Dmc_loc_lam = sm.lambdify(qL + pL, Dmc_loc, cse=True)
    punkt_loc_lam = sm.lambdify(qL + pL, punkt_loc, cse=True)
    eingepraegt_lam = sm.lambdify([fx, fy, fz] + qL + pL + rhs_subs,
                                  react_forces, cse=True)
    ang_momentum_lam = sm.lambdify(qL + pL, ang_momentum, cse=True)



Numerical Integration
=====================

.. jupyter-execute::

    # Input values

    r1 = 1.5                       # radius of the ball
    m1 = 1.0                       # mass of the ball
    m11 = m1 / 5.                  # mass of the red dot
    m_link1 = 0.5 * m1             # mass of the link
    l1 = 6.0                       # length of the massless rod of the pendulum
    k1 = 1000.                     # 'spring constant' of the balls
    reibung1 = 0.0                 # friction in the joints
    q1x, q1y, q1z = 0.2, 0.2, 0.2  # initial deflection of the first rod

    omega1 = 7.5                  # initial rot. speed of ball_i around A[i].y
    u1x, u1y, u1z = 0.0, omega1, 0.0  # initial ang. velocity of the first rod
    intervall = 5.0

    schritte = 100 * int(intervall)
    times = np.linspace(0., intervall, schritte)
    iXX1 = 2./5. * m1 * r1**2
    iYY1 = iXX1
    iZZ1 = iXX1

    pL_vals = np.array([m1, m11, m_link1, 9.8, r1, l1, iXX1, iYY1, iZZ1,
                       reibung1, k1])

    y0 = [q1x, q1y, q1z] + [0.0, 0.0, 0.0] * (n-1) + [u1x, u1y,
        u1z] + [0.0, u1y, 0.0] * (n-1)

    t_span = (0.0, intervall)

    resultat1 = solve_ivp(rhs_gen, t_span, y0, t_eval=times, args=(pL_vals,),
                          method='Radau',
                          )

    resultat = resultat1.y.T

    print('Shape of resultat', resultat.shape)
    print(f"To numerically integrate an intervall of {intervall:.3f} sec "
          f"the routine cycled {resultat1.nfev:,} times")


Calculate the Reaction Forces at the Suspension Point
-----------------------------------------------------

The accelerations needed are calculated numerically and stored in ``RHS``


.. jupyter-execute::

    RHS = np.empty((resultat.shape))
    for i in range(resultat.shape[0]):
        res = np.array(resultat[i])
        RHS[i] = rhs_gen(0.0, res, pL_vals)

    react_x = np.empty(resultat.shape[0])
    react_y = np.empty(resultat.shape[0])
    react_z = np.empty(resultat.shape[0])


    def func_react(x0, args):
        return eingepraegt_lam(*x0, *args).squeeze()


    x0 = np.array([0., 0., 0.])
    for i in range(resultat.shape[0]):
        args = np.array([*resultat[i], *pL_vals, *RHS[i, 3*n:]])
        loesung = root(func_react, x0, args=args)
        react_x[i] = loesung.x[0]
        react_y[i] = loesung.x[1]
        react_z[i] = loesung.x[2]
        x0 = loesung.x


Plot Energy, Angular Speeds, Reaction Forces and Angular Momentum
-----------------------------------------------------------------

.. jupyter-execute::

    pot_np = np.empty(schritte)
    kin_np = np.empty(schritte)
    spring_np = np.empty(schritte)
    total_np = np.empty(schritte)

    for i in range(schritte):
        zeit = times[i]
        pot_np[i] = pot_lam(*[resultat[i, j] for j in range(resultat.shape[1])],
                            *pL_vals)
        kin_np[i] = kin_lam(*[resultat[i, j] for j in range(resultat.shape[1])],
                            *pL_vals)
        spring_np[i] = spring_lam(*[resultat[i, j]
                                    for j in range(resultat.shape[1])],
                                  *pL_vals)
        total_np[i] = pot_np[i] + kin_np[i] + spring_np[i]

    if reibung1 == 0.:
        total_max = np.max(total_np)
        total_min = np.min(total_np)
        print('deviation of total energy from constant is {:.5f} % of max. '
              'total energy'.format((total_max - total_min)/total_max*100))

    fig, ax = plt.subplots(4, 1, figsize=(8, 10), layout='constrained',
                           sharex=True)
    ax[0].plot(times, pot_np, label='gravitational potential energy')
    ax[0].plot(times, kin_np, label='kinetic energy')
    ax[0].plot(times, spring_np, label='spring potential energy')
    ax[0].plot(times, total_np, label='total energy')
    msg = r'$\mu$'
    ax[0].set_title(f"Energies of the system, {msg} = {reibung1}")
    _ = ax[0].legend()

    for i in range(n, 2*n):
        ax[1].plot(times, resultat[:, 3*i+1], label='rotational speed of '
                   f'body {i - n} in Y direction in its coordinate system')
    ax[1].set_title('Rotational speeds')
    ax[1].set_ylabel('Rotational speed')
    _ = ax[1].legend()

    ax[2].plot(times, react_x, label='reaction_X')
    ax[2].plot(times, react_y, label='reaction_Y')
    ax[2].plot(times, react_z, label='reaction_Z')
    ax[2].set_title('Reaction forces at the suspension point')
    ax[2].set_xlabel('Time')
    ax[2].set_ylabel('Reaction force [N]')
    _ = ax[2].legend()

    max_y = np.max(ang_momentum_lam(*[resultat[:, j]
                                      for j in range(resultat.shape[1])],
                                *pL_vals)[1])
    min_y = np.min(ang_momentum_lam(*[resultat[:, j]
                                      for j in range(resultat.shape[1])],
                                *pL_vals)[1])
    error = (max_y - min_y) / max_y * 100.
    print('deviation of Y - component of ang. momentum from being constant is '
          '{:.5f} % of max. angular momentum'.format(error))
    ax[3].plot(
        times, ang_momentum_lam(
            *[resultat[:, j] for j in range(resultat.shape[1])],
            *pL_vals)[0], label='Angular momentum X')
    ax[3].plot(
        times, ang_momentum_lam(
            *[resultat[:, j] for j in range(resultat.shape[1])],
            *pL_vals)[1], label='Angular momentum Y')
    ax[3].plot(
        times, ang_momentum_lam(
            *[resultat[:, j] for j in range(resultat.shape[1])],
            *pL_vals)[2], label='Angular momentum Z')
    ax[3].set_title('Angular momentum')
    ax[3].set_ylabel('Angular momentum [kg*m^2/s]')
    _ = ax[3].legend()

Animation with PyDy Visualization
=================================

``groesse`` is an empirical factor to get the right size of the animation,
found by trial and error.

.. jupyter-execute::

    groesse = 25.0

    farben = ['orange', 'blue', 'green', 'red', 'yellow']
    if n > len(farben):
        raise ValueError('More colors must be  given')
    viz_frames = []

    for i, (ball, point, link) in enumerate(zip(balls, points, links)):
        ball_shape = Sphere(name='sphere{}'.format(i),
                            radius=r,
                            color=farben[i])

        viz_frames.append(VisualizationFrame('ball_frame{}'.format(i),
                                             ball,
                                             ball_shape))

        point_shape = Sphere(name='point{}'.format(i), radius=0.5 / groesse,
                             color='black')

        viz_frames.append(VisualizationFrame('point_frame{}'.format(i),
                                             A[i],
                                             point,
                                             point_shape))

        link_shape = Cylinder(name='cylinder{}'.format(i),
                              radius=0.25 / groesse,
                              length=l,
                              color='red')

        viz_frames.append(VisualizationFrame('link_frame{}'.format(i),
                                             link,
                                             link_shape))

    scene = Scene(N, P0, *viz_frames)

    scene.times = times
    pL_vals_scene = [val / groesse for val in pL_vals]
    scene.constants = dict(zip(pL, pL_vals_scene))
    scene.states_symbols = q + u
    scene.states_trajectories = resultat

    scene.display_jupyter(axes_arrow_length=20)

Total running time

.. jupyter-execute::

    print(f"It took {time.time() - start:.2f} sec to run the simulation")
