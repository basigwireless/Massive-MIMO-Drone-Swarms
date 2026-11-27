"""
Uniform distribution of drones inside a 3-D sphere: simulation vs. theory.

Python port of ThreeD_sphere_distribution_ThSim.m
Original MATLAB author: Prabhu Chandhar, Chandhar Research Labs, Chennai, India.

Drones are dropped uniformly inside a sphere of radius R centred on the
ground-station array. In spherical coordinates the three random variables are
independent with PDFs

    elevation  f_Theta(theta) = sin(theta) / 2,   0 <= theta <= pi
    azimuth    f_Phi(phi)     = 1 / (2*pi),       0 <= phi   <= 2*pi
    distance   f_D(r)         = 3 r^2 / R^3,      0 <= r     <= R

The script draws samples, compares their histograms with these PDFs and plots
a 3-D scatter of the drone positions.
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np


def sample_uniform_sphere(radius, num_samples, rng):
    """Draw points uniformly distributed inside a sphere.

    Returns elevation, azimuth, distance and Cartesian x, y, z arrays.
    """
    rvals = 2 * rng.random(num_samples) - 1
    el = np.arcsin(rvals) + np.pi / 2          # elevation angle in [0, pi]
    az = 2 * np.pi * rng.random(num_samples)   # azimuth angle in [0, 2*pi]
    dist = radius * rng.random(num_samples) ** (1 / 3)  # distance in [0, R]

    x = dist * np.cos(az) * np.sin(el)
    y = dist * np.sin(az) * np.sin(el)
    z = dist * np.cos(el)
    return el, az, dist, x, y, z


def empirical_pdf(samples, bins=1000):
    counts, edges = np.histogram(samples, bins=bins, density=True)
    centres = 0.5 * (edges[:-1] + edges[1:])
    return centres, counts


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--radius", type=float, default=1000.0, help="sphere radius R in m")
    parser.add_argument("--samples", type=int, default=1_000_000, help="number of drone positions")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--outdir", default=os.path.join(os.path.dirname(__file__), "..", "figures"))
    parser.add_argument("--show", action="store_true", help="display figures interactively")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    R = args.radius
    el, az, dist, x, y, z = sample_uniform_sphere(R, args.samples, rng)

    # Elevation angle PDF
    fig, ax = plt.subplots(figsize=(6, 4.2))
    b, a = empirical_pdf(el)
    th = np.linspace(0, np.pi, 400)
    ax.plot(b, a, color="tab:blue", lw=1, label="Simulation")
    ax.plot(th, np.sin(th) / 2, "r", lw=2, label="Theory: sin(θ)/2")
    ax.set(xlim=(0, np.pi), ylim=(0, 1), xlabel="Elevation angle θ (rad)",
           ylabel=r"PDF $f_\Theta(\theta)$", title="Elevation angle distribution")
    ax.grid(alpha=0.3); ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(args.outdir, "sphere_elevation_pdf.png"), dpi=150)

    # Azimuth angle PDF
    fig, ax = plt.subplots(figsize=(6, 4.2))
    b, a = empirical_pdf(az)
    ax.plot(b, a, color="tab:blue", lw=1, label="Simulation")
    ax.axhline(1 / (2 * np.pi), color="r", lw=2, label="Theory: 1/(2π)")
    ax.set(xlim=(0, 2 * np.pi), ylim=(0, 0.3), xlabel="Azimuth angle φ (rad)",
           ylabel=r"PDF $f_\Phi(\phi)$", title="Azimuth angle distribution")
    ax.grid(alpha=0.3); ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(args.outdir, "sphere_azimuth_pdf.png"), dpi=150)

    # Distance PDF
    fig, ax = plt.subplots(figsize=(6, 4.2))
    b, a = empirical_pdf(dist)
    r = np.linspace(0, R, 400)
    ax.plot(b, a, color="tab:blue", lw=1, label="Simulation")
    ax.plot(r, 3 * r ** 2 / R ** 3, "r", lw=2, label="Theory: 3r²/R³")
    ax.set(xlim=(0, R), xlabel="Distance r (m)", ylabel=r"PDF $f_D(r)$",
           title="Distance distribution")
    ax.grid(alpha=0.3); ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(args.outdir, "sphere_distance_pdf.png"), dpi=150)

    # 3-D scatter of the first 1000 drones
    fig = plt.figure(figsize=(6, 5.5))
    ax = fig.add_subplot(projection="3d")
    n = min(1000, args.samples)
    ax.scatter(x[:n], y[:n], z[:n], s=4, c=z[:n], cmap="viridis")
    ax.scatter([0], [0], [0], s=60, c="r", marker="^", label="Ground station array")
    ax.set(xlabel="x (m)", ylabel="y (m)", zlabel="z (m)",
           title=f"{n} drones uniformly distributed in a sphere (R = {R:g} m)")
    ax.legend(loc="upper left")
    fig.tight_layout(); fig.savefig(os.path.join(args.outdir, "sphere_drone_positions_3d.png"), dpi=150)

    print(f"Saved figures to {os.path.abspath(args.outdir)}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
