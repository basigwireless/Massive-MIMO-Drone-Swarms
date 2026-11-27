"""
Accuracy of near-field distance approximations for a large planar array.

Python port of ThreeD_Dist_Approximation_error.m
Original MATLAB author: Prabhu Chandhar, Chandhar Research Labs, Chennai, India.

A drone sits at distance R from the first element of an Mx x My uniform planar
array, at elevation theta and azimuth phi. The distance from element (p, q) to
the drone is compared for three models:

    exact              R * sqrt(1 + (a_pq^2)/R^2 - 2 sin(theta) b_pq / R)
    with higher order  R + a_pq^2 / (2R) - sin(theta) b_pq      (Fresnel / spherical wave)
    without h.o. term  R - sin(theta) b_pq                      (far-field / plane wave)

where a_pq^2 = ((p-1) dx)^2 + ((q-1) dy)^2 and
      b_pq   = (p-1) dx cos(phi) + (q-1) dy sin(phi).

Figures produced:
  * distance to each element for the last (theta, phi) on the grid
  * relative error (%) of the last element over all (theta, phi) for both
    approximations
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np


def element_distances(R, th, ph, Mx, My, dx, dy):
    """Return (exact, with_ho, without_ho) distances for every array element.

    th and ph may be arrays; the element dimensions are appended as the last axes.
    """
    p = np.arange(Mx)[:, None]  # (p-1)
    q = np.arange(My)[None, :]  # (q-1)
    th = np.asarray(th)[..., None, None]
    ph = np.asarray(ph)[..., None, None]
    a2 = (p * dx) ** 2 + (q * dy) ** 2
    b = p * dx * np.cos(ph) + q * dy * np.sin(ph)

    without_ho = R - np.sin(th) * b
    with_ho = R + a2 / (2 * R) - np.sin(th) * b
    exact = R * np.sqrt(1 + a2 / R ** 2 - 2 * np.sin(th) * b / R)
    return exact, with_ho, without_ho


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--freq", type=float, default=2.5e9, help="carrier frequency in Hz")
    parser.add_argument("--Mx", type=int, default=100, help="elements along x")
    parser.add_argument("--My", type=int, default=1, help="elements along y")
    parser.add_argument("--spacing", type=float, default=1.5, help="element spacing in wavelengths")
    parser.add_argument("--R", type=float, default=25.0, help="distance first element -> drone (m)")
    parser.add_argument("--outdir", default=os.path.join(os.path.dirname(__file__), "..", "figures"))
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    c = 3e8
    lam = c / args.freq
    dx = dy = args.spacing * lam
    R = args.R

    th_grid = np.arange(0, np.pi + 1e-12, np.pi / 20)       # elevation 0..180 deg
    ph_grid = np.arange(0, 2 * np.pi + 1e-12, np.pi / 20)   # azimuth   0..360 deg
    TH, PH = np.meshgrid(th_grid, ph_grid, indexing="ij")

    exact, with_ho, without_ho = element_distances(R, TH, PH, args.Mx, args.My, dx, dy)
    # Relative error of the farthest (last) element, as in the MATLAB script
    last = (..., -1, -1)
    err_with = np.abs((exact[last] - with_ho[last]) * 100 / exact[last])
    err_without = np.abs((exact[last] - without_ho[last]) * 100 / exact[last])

    # Figure 1: per-element distance for the last (theta, phi) on the grid
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    idx = np.arange(1, args.Mx * args.My + 1)
    sl = (-1, -1)
    ax.plot(idx, without_ho[sl].ravel(order="F"), lw=2, label="Without higher-order term")
    ax.plot(idx, with_ho[sl].ravel(order="F"), "g", lw=2, label="With higher-order term")
    ax.plot(idx, exact[sl].ravel(order="F"), "r--", lw=2, label="Exact")
    ax.set(xlabel="Antenna index", ylabel="Distance (m)",
           title=f"Element-to-drone distance (R = {R:g} m, θ = {np.degrees(th_grid[-1]):.0f}°, "
                 f"φ = {np.degrees(ph_grid[-1]):.0f}°)")
    ax.grid(alpha=0.3); ax.legend(loc="upper left")
    fig.tight_layout(); fig.savefig(os.path.join(args.outdir, "distance_vs_antenna_index.png"), dpi=150)

    # Figures 2 & 3: error surfaces
    for err, name, tag in [(err_with, "with higher-order term", "with_ho"),
                           (err_without, "without higher-order term", "without_ho")]:
        fig = plt.figure(figsize=(7, 5))
        ax = fig.add_subplot(projection="3d")
        surf = ax.plot_surface(np.degrees(PH), np.degrees(TH), err, cmap="viridis",
                               vmin=0, vmax=max(err_with.max(), err_without.max()),
                               edgecolor="k", linewidth=0.2)
        ax.set(xlabel="Azimuth angle (deg)", ylabel="Elevation angle (deg)",
               zlabel="Error (%)", xlim=(0, 360), ylim=(0, 180),
               zlim=(0, max(err_with.max(), err_without.max()) * 1.05),
               title=f"Relative distance error {name}")
        fig.colorbar(surf, ax=ax, shrink=0.6, pad=0.1, label="Error (%)")
        fig.tight_layout(); fig.savefig(os.path.join(args.outdir, f"distance_error_{tag}.png"), dpi=150)

    print(f"Max error with higher-order term:    {err_with.max():.2f} %")
    print(f"Max error without higher-order term: {err_without.max():.2f} %")
    print(f"Saved figures to {os.path.abspath(args.outdir)}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
