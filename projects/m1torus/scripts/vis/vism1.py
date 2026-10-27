#!/usr/bin/env python3
"""The m = 1 mode of a m1torus run against Bannikova et al. (2026).

Reads <prefix>_m1.ev (live analysis) and phantom's .ev files; draws the panels the paper
uses - A_1..A_3 against log time (their Fig. 7), the central-mass and torus-barycentre
tracks (Fig. 3), the ratios k = A_2/A_1 and Q_3 (Fig. 18), the virial quantity (Fig. 4)
and the mean elements - and prints the saturated averages next to their Table 2.

    vism1.py [--run DIR] [--window 100,150]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "common"))
import phantom_io as pio  # noqa: E402

P_ORB = 2.0 * np.pi
N_PAPER = 131072

# Table 2 and Sect. 9 of the paper, run 0.5-60-0.3 (N = 128k, 1200 sectors)
PAPER = dict(A1=40.21, A2=8.46, A3=1.77, k=0.21, Q3=0.99, r_c=0.024, r_tb=0.24)


def find_run(start):
    p = Path(start).resolve()
    for q in (p, *p.parents[:3]):
        if list(q.glob("*_m1.ev")):
            return q
    raise FileNotFoundError(f"no *_m1.ev at or above {p}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", default=".")
    ap.add_argument("--window", default=None,
                    help="orbits t0,t1 for the saturated averages (default: last third)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)

    run = find_run(a.run)
    m1 = pio.read_ev(sorted(run.glob("*_m1.ev"))[0])
    par = pio.read_setup(sorted(run.glob("*.setup"))[0])
    ev = pio.read_ev_series(run, sorted(run.glob("*.setup"))[0].stem)
    t = m1["time"] / P_ORB
    s0 = m1["sigma0"]
    n = m1["n"][0]
    r_c = np.hypot(m1["xc"], m1["yc"])
    r_tb = np.hypot(m1["xtb"], m1["ytb"])
    a1, a2, a3 = (m1[f"A{m}"] for m in (1, 2, 3))

    t0, t1 = ((float(x) for x in a.window.split(",")) if a.window
              else (t[-1] * 2.0 / 3.0, t[-1]))
    w = (t >= t0) & (t <= t1)
    mean = {k: float(np.mean(v[w])) for k, v in
            dict(A1=a1, A2=a2, A3=a3, r_c=r_c, r_tb=r_tb, s0=s0).items()}
    k_ratio = mean["A2"] / mean["A1"]
    q3 = (mean["A3"] / mean["A1"]) / k_ratio**2

    fig, ax = plt.subplots(2, 3, figsize=(17, 10))
    x = ax[0, 0]
    tt = np.maximum(t, 1e-2)
    for amp, lab, c in ((a1, "A$_1$", "C0"), (a2, "A$_2$", "C3"), (a3, "A$_3$", "C2")):
        x.semilogx(tt, amp / s0, color=c, lw=1.2, label=lab)
    x.axhline(PAPER["A1"] / (N_PAPER - 1) * 1200, color="C0", ls=":", lw=1,
              label="paper A$_1$ (Table 2)")
    x.set_ylabel(r"$A_m/\Sigma_0$")
    x.set_title("Mode amplitudes (Fig. 7)", fontweight="bold")
    x.legend(fontsize=9)

    x = ax[0, 1]
    x.plot(m1["xc"], m1["yc"], color="C1", lw=0.6, label=r"$r_c$ central mass")
    x.plot(par.get("M_tor", 0.1) * m1["xtb"], par.get("M_tor", 0.1) * m1["ytb"],
           color="C0", lw=0.6, alpha=0.6, label=r"$M_{\rm tor} r_{tb}$")
    x.set_aspect("equal")
    x.set_title("Central mass and torus barycentre (Fig. 3)", fontweight="bold")
    x.legend(fontsize=9)

    x = ax[0, 2]
    x.plot(t, r_c, color="C1", lw=1, label=r"$|r_c|$")
    x.axhline(PAPER["r_c"], color="C1", ls=":", lw=1, label="paper")
    x.set_ylabel(r"$|r_c|$")
    x.set_title("Central-mass offset", fontweight="bold")
    x.legend(fontsize=9)

    x = ax[1, 0]
    with np.errstate(divide="ignore", invalid="ignore"):
        x.plot(t, a2 / a1, color="C3", lw=1, label=r"$k = A_2/A_1$")
        x.plot(t, (a3 / a1) / (a2 / a1) ** 2, color="C2", lw=1, alpha=0.7, label=r"$Q_3$")
    x.axhspan(0.20, 0.25, color="C3", alpha=0.12, label="paper k 0.20-0.25")
    x.axhline(1.0, color="C2", ls=":", lw=1)
    x.set_ylim(0, 2)
    x.set_title("Harmonic ratios (Fig. 18)", fontweight="bold")
    x.legend(fontsize=9)

    x = ax[1, 1]
    if ev:
        te = ev["time"] / P_ORB
        x.semilogx(np.maximum(te, 1e-2), 2 * ev["ekin"] + ev["epot"], color="C2", lw=1)
        x.set_ylabel(r"$2E_{\rm kin}+E_{\rm pot}$ (whole system)")
    x.set_title("Virial quantity (Fig. 4)", fontweight="bold")

    x = ax[1, 2]
    x.plot(t, m1["e_mean"], label=r"$\langle e\rangle$")
    x.plot(t, m1["e_sigma"], label=r"$\sigma_e$ (Rayleigh)")
    x.axhline(0.27, color="C1", ls=":", lw=1, label=r"paper $\sigma_e$")
    x2 = x.twinx()
    x2.plot(t, m1["i_mean"], color="C3", lw=0.8)
    x2.set_ylabel(r"$\langle i\rangle$ [deg]", color="C3")
    x.set_title("Mean elements", fontweight="bold")
    x.legend(fontsize=9)

    for x in ax.flat:
        x.grid(True, alpha=0.3)
    for x in (ax[0, 2], ax[1, 0], ax[1, 2]):
        x.set_xlabel(r"$t\;[P_{\rm orb}]$")
    ax[0, 0].set_xlabel(r"$t\;[P_{\rm orb}]$")
    ax[1, 1].set_xlabel(r"$t\;[P_{\rm orb}]$")
    info = (f"N = {int(n)}, M_tor = {par.get('M_tor')}, e_max = {par.get('e_max')}, "
            f"i_max = {par.get('i_max')}, R0 = {par.get('R0')}")
    fig.suptitle(f"m1torus: {info}", fontsize=13, fontweight="bold")
    fig.tight_layout()
    out = Path(a.out) if a.out else run / "figs"
    out.mkdir(parents=True, exist_ok=True)
    f = out / "m1.png"
    fig.savefig(f, dpi=170, bbox_inches="tight")
    plt.close(fig)

    scale = (N_PAPER - 1) / (n)          # particles per sector scale with N
    print(f"window {t0:.0f}-{t1:.0f} orbits, N = {int(n)} (paper amplitudes x {1/scale:.3f})")
    print(f"{'':10s} {'ours':>9s} {'paper':>9s}")
    print(f"{'A1/S0':10s} {mean['A1']/mean['s0']:9.3f} {PAPER['A1']/(N_PAPER-1)*1200:9.3f}")
    print(f"{'A2/S0':10s} {mean['A2']/mean['s0']:9.3f} {PAPER['A2']/(N_PAPER-1)*1200:9.3f}")
    print(f"{'A3/S0':10s} {mean['A3']/mean['s0']:9.3f} {PAPER['A3']/(N_PAPER-1)*1200:9.3f}")
    print(f"{'k':10s} {k_ratio:9.3f} {PAPER['k']:9.3f}")
    print(f"{'Q3':10s} {q3:9.3f} {PAPER['Q3']:9.3f}")
    print(f"{'|r_c|':10s} {mean['r_c']:9.4f} {PAPER['r_c']:9.4f}")
    print(f"{'|r_tb|':10s} {mean['r_tb']:9.4f} {PAPER['r_tb']:9.4f}")
    print(f)
    return 0


if __name__ == "__main__":
    sys.exit(main())
