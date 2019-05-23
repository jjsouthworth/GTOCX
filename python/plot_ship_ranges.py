import argparse
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

# Define constants
k8 = -1.94316e-12
k7 = 3.7516e-10
k6 = -2.70559e-8
k5 = 9.70521e-7
k4 = -1.88428e-5
k3 = 0.000198502
k2 = -0.0010625
k1 = 0.0023821
k0 = 0.00287729
kpc_per_myr = 30856775814671900. / 31557600e6


def circular_velocity(r):
    # Compute powers of the radius
    r2 = r*r
    r3 = r2*r
    r4 = r2*r2
    r5 = r3*r2
    r6 = r3*r3
    r7 = r4*r3
    r8 = r4*r4

    # Compute the circular velocity and mean motion.
    vc = 1.0 / (k8*r8 + k7*r7 + k6*r6 + k5*r5 + k4*r4 + k3*r3 + k2*r2 + k1*r + k0) / kpc_per_myr
    return vc


def equations_of_motion(_, x):
    r = np.linalg.norm(x[:3])
    vc = circular_velocity(r)
    fr = vc**2 / r

    acc = np.array(-x[:3]/r*fr)
    dx = np.concatenate((x[3:], acc))
    return dx


def elem_to_state(elem, t=0.0):
    # Extract the elements into variables.
    elem = np.array(elem)
    r = elem[:, 0]
    inc, omega, phi, theta = np.radians(elem[:, 1:]).transpose()

    # Compute the circular velocity and mean motion.
    vc = circular_velocity(r)
    n = vc/r

    # Compute repeated variables.
    nt_phi = n*t + phi
    c_nt_phi = np.cos(nt_phi)
    s_nt_phi = np.sin(nt_phi)
    c_omega = np.cos(omega)
    s_omega = np.sin(omega)
    c_inc = np.cos(inc)
    s_inc = np.sin(inc)

    # Compute the Cartesian state.
    x = r*(c_nt_phi*c_omega - s_nt_phi*c_inc*s_omega)
    y = r*(c_nt_phi*s_omega + s_nt_phi*c_inc*c_omega)
    z = r*(s_nt_phi*s_inc)

    vx = vc*(-s_nt_phi*c_omega - c_nt_phi*c_inc*s_omega)
    vy = vc*(-s_nt_phi*s_omega + c_nt_phi*c_inc*c_omega)
    vz = vc*(c_nt_phi*s_inc)

    return np.column_stack((x, y, z, vx, vy, vz))


def plot_stars(elems, t0=0.0, tf=90.0, outfile=None, draw=True):

    fig, ax = plt.subplots()

    # Create the time-vector
    tvec = np.arange(t0, tf, 1)
    if not tvec.size:
        tvec = np.array([t0])

    # Plot each star's path.
    for t in tvec:
        states = elem_to_state(elems, t)
        x, y, _, _, _, _ = states.transpose()
        ax.plot(x, y, 'k.', markersize=1)
        ax.plot(x[0], y[0], 'r*')

    ax.grid(True)
    ax.axis('equal')

    if draw:
        plt.show()

    if outfile is not None:
        fig.savefig(outfile)

    return ax


def rotz(x, theta):
    r = np.array([[np.cos(theta), -np.sin(theta), 0.0],
                  [np.sin(theta), np.cos(theta), 0.0],
                  [0.0, 0.0, 1.0]])

    return np.dot(r, x)


def unit(x):
    return x / np.linalg.norm(x)


def integrate_trajectories(initial_elems, dv, thetas):
    traj = list()

    # Create an integrator
    integrator = integrate.ode(equations_of_motion).set_integrator('dopri5')

    # Compute starting initial conditions and apply delta-v.
    x0 = elem_to_state([initial_elems], 0.0)[0]
    vhat = unit(x0[3:])

    for theta in thetas:
        # Compute and apply delta-v vector
        deltav = rotz(dv*vhat, theta)
        xi = x0 + np.concatenate(([0, 0, 0], deltav))
        integrator.set_initial_value(xi, t=0.0)

        # Integrate to every 1 Myr step.
        x = np.zeros((91, 6))
        times = range(91)
        for t in times:
            x[t, :] = integrator.integrate(t)

        traj.append(x)

    return traj


def plot_trajectories(ax, dv, traj, color='k'):
    for i, x in enumerate(traj):
        label = "{:.0f} km/s".format(dv) if i == 0 else None
        linespec = color + '-'
        ax.plot(x[:, 0], x[:, 1], linespec, label=label)

    return ax


def main(args):
    star_data = np.loadtxt(args.star_file, delimiter=',', skiprows=1)

    ids = star_data[:, 0].astype(int)
    star_elem = star_data[:, 1:]
    ax = plot_stars(star_elem, args.time, args.time, draw=False)

    # Generate trajectories with various delta-v's.
    dvs = np.array([100, 200, 300, 400, 750])  # km/s
    colors = ['b', 'r', 'g', 'm', 'y']
    for dv, c in zip(dvs, colors):
        thetas = np.linspace(0, 2*np.pi, 10)
        traj = integrate_trajectories(star_elem[0, :], dv / kpc_per_myr, thetas)

        # Plot trajectories.
        plot_trajectories(ax, dv, traj, c)
        ax.legend()
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Draw or animate the GTOCX galaxy.")
    parser.add_argument("star_file", help="Star data file.")
    parser.add_argument("-t", "--time", type=float, default=0.0, help="The stars will be drawn at the given time.")

    input_args = parser.parse_args()
    main(input_args)
