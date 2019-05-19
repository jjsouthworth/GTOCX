import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani


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

kpc = 30856775814671900  # km
yr = 31557600  # sec
Myrs = 1e6*yr  # sec


def elem_to_state(elem, t=0.0):
    # Extract the elements into variables.
    elem = np.array(elem)
    r = elem[:, 0]
    inc, omega, phi, theta = np.radians(elem[:, 1:]).transpose()

    # Compute powers of the radius
    r2 = r*r
    r3 = r2*r
    r4 = r2*r2
    r5 = r3*r2
    r6 = r3*r3
    r7 = r4*r3
    r8 = r4*r4

    # Compute the circular velocity and mean motion.
    vc = 1.0 / (k8*r8 + k7*r7 + k6*r6 + k5*r5 + k4*r4 + k3*r3 + k2*r2 + k1*r + k0)
    n = vc/(r*kpc) * Myrs

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


def animate_stars(elems, t0=0.0, tf=90.0, num_frames=500, outfile=None, draw=True):

    fig, ax = plt.subplots(1, 2, sharey=True, figsize=(15, 7.5))
    line1, = ax[0].plot([0, 0, -32, 32], [32, -32, 0, 0], 'k.', markersize=1)
    line2, = ax[1].plot([0, 0, -32, 32], [5, -5, 0, 0], 'k.', markersize=1)
    line3, = ax[0].plot([0, 0], [0, 0], 'r*', zorder=1)
    line4, = ax[1].plot([0, 0], [0, 0], 'r*', zorder=1)

    # Create the prescribed time-vector.
    tvec = np.linspace(t0, tf, num_frames)

    def update(i):
        t = tvec[i]
        states = elem_to_state(elems, t)
        line1.set_xdata(states[:, 0])
        line1.set_ydata(states[:, 1])
        line2.set_xdata(states[:, 2])
        line2.set_ydata(states[:, 0])
        line3.set_xdata(states[0, 0])
        line3.set_ydata(states[0, 1])
        line4.set_xdata(states[0, 2])
        line4.set_ydata(states[0, 0])

        return line1, line2, line3, line4

    update(0)
    for a in ax:
        a.axis('equal')
        a.grid(True)

    animation = ani.FuncAnimation(fig, update, len(tvec), blit=True, interval=1.0)
    if draw:
        plt.draw()
        plt.show()

    if outfile is not None:
        animation.save(outfile)


def main(args):
    draw = not args.suppress
    star_data = np.loadtxt(args.star_file, delimiter=',', skiprows=1)

    ids = star_data[:, 0].astype(int)
    star_elem = star_data[:, 1:]

    if args.animate:
        animate_stars(star_elem, args.time, outfile=args.outfile, draw=draw)
    else:
        plot_stars(star_elem, args.time, args.time, args.outfile, draw)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Draw or animate the GTOCX galaxy.")
    parser.add_argument("star_file", help="Star data file.")
    mut_exc_group = parser.add_mutually_exclusive_group()
    mut_exc_group.add_argument("-t", "--time", type=float, default=0.0, help="Time of the galaxy in Myrs. 0 <= t <= 90.")
    mut_exc_group.add_argument("-a", "--animate", action='store_true', help="Animate the trajectory.")
    parser.add_argument("-o", "--outfile", help="Save image or animation to given output file. Extension will "
                                                "determine the file type (.gif, .mp4, .png, .jpg, etc).")
    parser.add_argument("-s", "--suppress", action='store_true', help="Use to suppress drawing the plot/animation.")

    input_args = parser.parse_args()
    main(input_args)
