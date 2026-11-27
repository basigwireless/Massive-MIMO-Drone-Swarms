"""
Moments of the normalised interference for a drone swarm served by a
massive-MIMO uniform linear array (ULA) with spherical-wave channels.

Python port of ThreeD_sphere_Interference_Distribution_Moments_ThSim.m
Original MATLAB author: Prabhu Chandhar, Chandhar Research Labs, Chennai, India.

K drones are uniformly distributed inside a sphere of radius R. The channel
from the M-element ULA (spacing delta) to drone k with distance d_k and
direction cosine z_k = sin(theta_k) cos(phi_k) is modelled with the
second-order (Fresnel) phase approximation

    g_k[m] = exp(-j 2 pi / lambda * ( d_k + m^2 delta^2 / (2 d_k) - (m-1) delta z_k )),
    m = 1..M.

With matched-filter (conjugate) beamforming to drone 1, the normalised
interference from the other K-1 drones is

    U = (1 / M^2) * sum_{j=2..K} | g_j g_1^H |^2.

The script estimates E[U], E[U^2], E[U^3], E[U^4] by Monte-Carlo simulation
and compares them with the closed-form expressions from the MATLAB script.
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np

from sphere_distribution import sample_uniform_sphere


def channel(d, zc, M, delta, lam):
    """Spherical-wave ULA channel vectors, shape (N, M)."""
    m = np.arange(1, M + 1)[None, :]
    d = d[:, None]
    zc = zc[:, None]
    phase = d + (m ** 2 * delta ** 2) / (2 * d) - delta * (m - 1) * zc
    return np.exp(-1j * 2 * np.pi * phase / lam)


def simulate_moments(M_values, K, R, lam, delta, num_samples, rng, chunk=20000):
    """Monte-Carlo estimate of the first four moments of U for each M."""
    # Drone positions: K drones x num_samples independent drops
    d = np.empty((K, num_samples))
    zc = np.empty((K, num_samples))
    for k in range(K):
        el, az, dist, *_ = sample_uniform_sphere(R, num_samples, rng)
        d[k] = dist
        zc[k] = np.sin(el) * np.cos(az)

    moments = np.zeros((4, len(M_values)))
    for iM, M in enumerate(M_values):
        U = np.zeros(num_samples)
        for s in range(0, num_samples, chunk):
            e = min(s + chunk, num_samples)
            gk = channel(d[0, s:e], zc[0, s:e], M, delta, lam)
            for j in range(1, K):
                gj = channel(d[j, s:e], zc[j, s:e], M, delta, lam)
                U[s:e] += np.abs(np.sum(gj * np.conj(gk), axis=1)) ** 2
        U /= M ** 2
        moments[:, iM] = [np.mean(U ** n) for n in range(1, 5)]
        print(f"M = {M:4d}  E[U] = {moments[0, iM]:.4e}")
    return moments


def theoretical_moments(M_values, K):
    """Closed-form moments of U (from the original MATLAB script)."""
    out = np.zeros((4, len(M_values)))
    for i, M in enumerate(np.asarray(M_values, dtype=float)):
        a1 = (M ** 3 * 2 / 3) + (M / 3)
        a2 = M ** 2
        b1 = (M * (11 * M ** 4 + 5 * M ** 2 + 4)) / 20
        b2 = (M * M * (2 * M ** 2 + 1)) / 3
        b3 = M ** 3
        c1 = (M * (151 * M ** 6 + 70 * M ** 4 + 49 * M ** 2 + 45)) / 315
        c2 = (M * M * (11 * M ** 4 + 5 * M ** 2 + 4)) / 20
        c3 = ((M * (2 * M ** 2 + 1)) / 3) ** 2
        c4 = M * M * M * (2 * M ** 2 + 1) / 3
        c5 = M ** 4

        out[0, i] = (K - 1) / M
        out[1, i] = ((K - 1) * a1 + a2 * (K - 1) * (K - 2)) / M ** 4
        out[2, i] = ((K - 1) * b1 + 3 * (K - 1) * (K - 2) * b2
                     + (K - 1) * (K - 2) * (K - 3) * b3) / M ** 6
        out[3, i] = (c1 * (K - 1) + c2 * 4 * (K - 1) * (K - 2) + c3 * 3 * (K - 1) * (K - 2)
                     + c4 * 6 * (K - 1) * (K - 2) * (K - 3)
                     + c5 * (K - 1) * (K - 2) * (K - 3) * (K - 4)) / M ** 8
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--freq", type=float, default=2.4e9, help="carrier frequency in Hz")
    parser.add_argument("--R", type=float, default=1000.0, help="sphere radius in m")
    parser.add_argument("--K", type=int, default=10, help="number of drones")
    parser.add_argument("--spacing", type=float, default=0.5, help="element spacing in wavelengths")
    parser.add_argument("--samples", type=int, default=100_000, help="Monte-Carlo drops")
    parser.add_argument("--Mmin", type=int, default=10)
    parser.add_argument("--Mmax", type=int, default=100)
    parser.add_argument("--Mstep", type=int, default=10)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--outdir", default=os.path.join(os.path.dirname(__file__), "..", "figures"))
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    lam = 3e8 / args.freq
    delta = args.spacing * lam
    M_values = np.arange(args.Mmin, args.Mmax + 1, args.Mstep)
    rng = np.random.default_rng(args.seed)

    sim = simulate_moments(M_values, args.K, args.R, lam, delta, args.samples, rng)
    th = theoretical_moments(M_values, args.K)

    names = ["First", "Second", "Third", "Fourth"]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5))
    for n, ax in enumerate(axes.ravel()):
        ax.semilogy(M_values, sim[n], "bo-", ms=5, label="Simulation")
        ax.semilogy(M_values, th[n], "r--", lw=2, label="Theory")
        ax.set(xlabel="Number of antennas M", ylabel=f"E[U$^{n + 1}$]",
               title=f"{names[n]} moment of U")
        ax.grid(alpha=0.3, which="both"); ax.legend()
    fig.suptitle(f"Normalised interference moments (K = {args.K} drones, R = {args.R:g} m)")
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "interference_moments.png"), dpi=150)

    np.savetxt(os.path.join(args.outdir, "interference_moments.csv"),
               np.column_stack([M_values, sim.T, th.T]), delimiter=",",
               header="M,sim_m1,sim_m2,sim_m3,sim_m4,th_m1,th_m2,th_m3,th_m4", comments="")
    print(f"Saved figures to {os.path.abspath(args.outdir)}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
