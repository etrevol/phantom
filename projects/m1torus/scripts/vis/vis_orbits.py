#!/usr/bin/env python3
"""Orbital structure of the m = 1 mode, the paper's Figs. 2, 9, 13, B.1 and Sect. 7.2.

From the full dumps of a m1torus run:
  ecc_maps.png      face-on particles coloured by eccentricity at four times, with the
                    m = 1 arrow (paper Fig. 9)
  incl_peri.png     inclination against longitude of pericentre for e >= 0.5 at the
                    last dump (Fig. 13)
  ecc_hist.png      eccentricity histogram at the last dump with a Rayleigh fit (Fig. B.1)
  pattern.png       unwrapped phase Phi_1(t) and the pattern speed of the saturated mode
  corotating.png    surface density stacked in the frame of the mode (Fig. 2, lower
                    right) and its azimuthal profile with the Poisson-kernel fit (Sect. 7.2)

Elements are osculating, about the central mass. Prints the fitted numbers.

    vis_orbits.py [--run DIR] [--times 35,50,100,150] [--window 90,150]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.optimize import curve_fit  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "common"))
import phantom_io as pio
import labels  # noqa: E402

P_ORB = 2.0 * np.pi
plt.rcParams.update({"font.size": 10, "axes.labelsize": 11, "xtick.direction": "in",
                     "ytick.direction": "in", "xtick.top": True, "ytick.right": True})


def load(path):
    """Positions, velocities and elements about the central mass at one dump."""
    d = pio.read_dump(path)
    g, s = d["gas"], d["sinks"]
    keep = g["h"] > 0
    x = np.stack([g["x"], g["y"], g["z"]], 1)[keep]
    v = np.stack([g["vx"], g["vy"], g["vz"]], 1)[keep]
    xc = np.array([s["x"][0], s["y"][0], s["z"][0]])
    vc = np.array([s["vx"][0], s["vy"][0], s["vz"][0]])
    mu = s["m"][0]
    r, u = x - xc, v - vc
    rr = np.linalg.norm(r, axis=1)
    hv = np.cross(r, u)
    ev = np.cross(u, hv) / mu - r / rr[:, None]
    e = np.linalg.norm(ev, axis=1)
    inc = np.degrees(np.arccos(np.clip(hv[:, 2] / np.linalg.norm(hv, axis=1), -1, 1)))
    varpi = np.degrees(np.arctan2(ev[:, 1], ev[:, 0])) % 360.0
    return dict(t=d["header"]["time"] / P_ORB, x=x, xc=xc, e=e, inc=inc, varpi=varpi)


def nearest(dumps, times, t):
    return dumps[int(np.argmin(np.abs(times - t)))]


def poisson(psi, s0, s, k):
    return s0 + s * ((1 - k**2) / (1 - 2 * k * np.cos(psi) + k**2) - 1)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", default=".")
    ap.add_argument("--times", default="35,50,100,150")
    ap.add_argument("--window", default="90,150", help="saturated window, orbits")
    a = ap.parse_args(argv)
    run = Path(a.run).resolve()
    out = run / "figs"
    out.mkdir(exist_ok=True)
    dumps = sorted(run.glob("m1torus_[0-9][0-9][0-9][0-9][0-9]"))
    times = np.array([pio.read_dump(p)["header"]["time"] / P_ORB for p in dumps])
    m1 = pio.read_ev(run / "m1torus_m1.ev")
    tm = m1["time"] / P_ORB
    a1 = m1["A1"] / m1["sigma0"]
    w0, w1 = (float(x) for x in a.window.split(","))

    # Fig. 9: eccentricity maps
    ts = [float(x) for x in a.times.split(",")]
    fig, axs = plt.subplots(1, len(ts), figsize=(3.6 * len(ts) + 0.8, 4.0),
                            constrained_layout=True)
    for ax, t in zip(axs, ts):
        f = load(nearest(dumps, times, t))
        o = np.argsort(f["e"])
        sc = ax.scatter(f["x"][o, 0], f["x"][o, 1], c=f["e"][o], s=0.8, lw=0,
                        cmap="plasma", vmin=0, vmax=1, rasterized=True)
        k = int(np.argmin(np.abs(tm - f["t"])))
        L = 1.6 * a1[k] / 0.368
        ax.annotate("", xy=(L * np.cos(m1["Phi1"][k]), L * np.sin(m1["Phi1"][k])),
                    xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color="k", lw=1.6))
        ax.plot(*f["xc"][:2], "ko", ms=4)
        ax.set(xlim=(-2, 2), ylim=(-2, 2), aspect="equal", xlabel=r"$x/R_{\rm tor}$",
               title=rf"$t = {f['t']:.0f}\,P_{{\rm orb}}$")
        ax.set_xticks(np.arange(-1.5, 1.6, 0.5))
        ax.set_yticks(np.arange(-1.5, 1.6, 0.5))
    axs[0].set_ylabel(r"$y/R_{\rm tor}$")
    fig.colorbar(sc, ax=axs, shrink=0.8, pad=0.01, label="$e$")
    fig.suptitle(labels.wrapped(labels.model_label(run), 75), fontsize=10)
    fig.savefig(out / "ecc_maps.png", dpi=170)
    plt.close(fig)

    last = load(dumps[-1])

    # Fig. 13: inclination against longitude of pericentre, e >= 0.5
    sel = last["e"] >= 0.5
    k = int(np.argmin(np.abs(tm - last["t"])))
    phase = np.degrees(m1["Phi1"][k]) % 360.0
    fig, ax = plt.subplots(figsize=(6.2, 4.2), constrained_layout=True)
    ax.scatter(last["varpi"][~sel], last["inc"][~sel], s=0.6, lw=0, color="0.75",
               rasterized=True, label=r"$e<0.5$")
    ax.scatter(last["varpi"][sel], last["inc"][sel], s=1.2, lw=0, color="C1",
               rasterized=True, label=r"$e\geq0.5$")
    ax.axvline(phase, color="k", ls="--", lw=1, label=r"phase $\Phi_1$ of $m=1$")
    # orbits whose apocentre points along Phi_1 have their pericentre opposite
    ax.axvline((phase + 180.0) % 360.0, color="k", ls=":", lw=1,
               label=r"$\Phi_1 + 180^\circ$")
    ax.set(xlim=(0, 360), ylim=(0, 90), xlabel=r"longitude of pericentre $\varpi$ [deg]",
           ylabel=r"inclination $i$ [deg]", title=rf"$t = {last['t']:.0f}\,P_{{\rm orb}}$")
    ax.set_xticks(range(0, 361, 60))
    ax.legend(fontsize=8.5, loc="upper right", markerscale=6)
    fig.suptitle(labels.wrapped(labels.model_label(run), 75), fontsize=10)
    fig.savefig(out / "incl_peri.png", dpi=170)
    plt.close(fig)

    # Fig. B.1: eccentricity histogram with a Rayleigh fit
    hist, edges = np.histogram(last["e"], bins=50, range=(0, 1), density=True)
    cen = 0.5 * (edges[1:] + edges[:-1])
    (sig,), _ = curve_fit(lambda e, s: e / s**2 * np.exp(-e**2 / (2 * s**2)), cen, hist,
                          p0=[0.3])
    fig, ax = plt.subplots(figsize=(5.6, 3.9), constrained_layout=True)
    ax.stairs(hist, edges, color="0.3", lw=1.2, label="particles")
    ee = np.linspace(0, 1, 200)
    ax.plot(ee, ee / sig**2 * np.exp(-ee**2 / (2 * sig**2)), color="C3", lw=1.5,
            label=rf"Rayleigh, $\sigma = {sig:.3f}$")
    ax.plot(ee, ee / 0.27**2 * np.exp(-ee**2 / (2 * 0.27**2)), color="C3", ls=":", lw=1,
            label=labels.REF["bannikova26"] + r", $\sigma = 0.27$")
    ax.set(xlabel="$e$", ylabel="$dN/de$ (normalised)",
           title=rf"$t = {last['t']:.0f}\,P_{{\rm orb}}$")
    ax.legend(fontsize=9, frameon=False)
    fig.suptitle(labels.wrapped(labels.model_label(run), 75), fontsize=10)
    fig.savefig(out / "ecc_hist.png", dpi=170)
    plt.close(fig)

    # pattern speed of the saturated mode
    ph = np.unwrap(m1["Phi1"])
    w = (tm >= w0) & (tm <= w1)
    omega_p = np.polyfit(m1["time"][w], ph[w], 1)[0]       # rad per code time, Omega_0 = 1
    fig, ax = plt.subplots(figsize=(6.4, 3.9), constrained_layout=True)
    ax.plot(tm, np.degrees(ph), color="C0", lw=1)
    ax.plot(tm[w], np.degrees(np.polyval(np.polyfit(m1["time"][w], ph[w], 1),
                                         m1["time"][w])), color="C3", lw=1.5,
            label=rf"$\Omega_p = {omega_p:.4f}\,\Omega_0$, period ${2*np.pi/omega_p/P_ORB:.0f}\,P_{{\rm orb}}$")
    ax.set(xlabel=r"$t/P_{\rm orb}$", ylabel=r"$\Phi_1$ (unwrapped) [deg]")
    ax.legend(fontsize=9, frameon=False)
    fig.suptitle(labels.wrapped(labels.model_label(run), 75), fontsize=10)
    fig.savefig(out / "pattern.png", dpi=170)
    plt.close(fig)

    # stack in the frame of the mode: rotate each dump by -Phi_1
    nb, lim = 160, 2.0
    sigma = np.zeros((nb, nb))
    psi_all = []
    nst = 0
    for p, t in zip(dumps, times):
        if not (w0 <= t <= w1):
            continue
        f = load(p)
        k = int(np.argmin(np.abs(tm - t)))
        c, s = np.cos(-m1["Phi1"][k]), np.sin(-m1["Phi1"][k])
        xr = c * f["x"][:, 0] - s * f["x"][:, 1]
        yr = s * f["x"][:, 0] + c * f["x"][:, 1]
        hh, _, _ = np.histogram2d(xr, yr, bins=nb, range=[[-lim, lim], [-lim, lim]])
        sigma += hh
        psi_all.append(np.arctan2(yr, xr))
        nst += 1
    sigma /= max(nst, 1)
    psi = np.concatenate(psi_all)
    prof, pe = np.histogram(psi, bins=72, range=(-np.pi, np.pi))
    prof = prof / nst
    pc = 0.5 * (pe[1:] + pe[:-1])
    (s0, sfit, kfit), _ = curve_fit(poisson, pc, prof, p0=[prof.mean(), prof.mean() * 0.2,
                                                        0.2], bounds=([0, 0, 0], [np.inf, np.inf, 0.99]))
    fig = plt.figure(figsize=(11.5, 4.6), constrained_layout=True)
    ax1 = fig.add_subplot(1, 2, 1)
    im = ax1.imshow(np.log10(np.maximum(sigma.T, 0.5)), origin="lower",
                    extent=[-lim, lim, -lim, lim], cmap="turbo")
    ax1.set(xlabel=r"$x'/R_{\rm tor}$", ylabel=r"$y'/R_{\rm tor}$", aspect="equal",
            title=rf"stacked in the frame of $m=1$, $t = {w0:.0f}$-${w1:.0f}\,P_{{\rm orb}}$")
    fig.colorbar(im, ax=ax1, shrink=0.85, label=r"$\log_{10}$ particles per pixel")
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.plot(np.degrees(pc), prof / prof.mean(), color="0.3", lw=1.2, label="torus")
    ax2.plot(np.degrees(pc), poisson(pc, s0, sfit, kfit) / prof.mean(), color="C3", lw=1.5,
             label=rf"Poisson kernel, $k = {kfit:.3f}$")
    ax2.set(xlabel=r"$\psi = \phi - \Phi_1$ [deg]", ylabel=r"$\Sigma(\psi)/\Sigma_0$",
            xlim=(-180, 180))
    ax2.set_xticks(range(-180, 181, 60))
    ax2.legend(fontsize=9, frameon=False)
    fig.suptitle(labels.wrapped(labels.model_label(run), 75), fontsize=10)
    fig.savefig(out / "corotating.png", dpi=170)
    plt.close(fig)

    print(f"Rayleigh sigma at t = {last['t']:.0f}: {sig:.3f} (paper 0.27)")
    print(f"pattern speed {omega_p:.4f} Omega_0, period {2*np.pi/omega_p/P_ORB:.1f} orbits "
          f"(window {w0:.0f}-{w1:.0f})")
    print(f"Poisson-kernel k = {kfit:.3f} (paper 0.20-0.25); A2/A1 in the window "
          f"{np.mean(m1['A2'][w]) / np.mean(m1['A1'][w]):.3f}")
    dv = (last["varpi"][sel] - phase + 180.0) % 360.0 - 180.0
    mean_dv = np.degrees(np.angle(np.mean(np.exp(1j * np.radians(dv)))))
    print(f"e >= 0.5 at the end: {sel.sum()} of {sel.size}; their mean varpi - Phi_1 = "
          f"{mean_dv:.0f} deg")
    for name in ("ecc_maps", "incl_peri", "ecc_hist", "pattern", "corotating"):
        print(out / f"{name}.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
