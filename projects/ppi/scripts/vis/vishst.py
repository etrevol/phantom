#!/usr/bin/env python3
"""History of a PPI torus run: conservation, mass loss, shock capturing, centre of mass.

Reads phantom's <prefix>NN.ev (every two steps) and <prefix>_ppi.ev (every dtmax).

    vishst.py [--run DIR]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ppi_data as pdat  # noqa: E402


def rel(x):
    return x / x[0] - 1.0 if len(x) and x[0] != 0 else x


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", default=".")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)

    info = pdat.run_info(a.run)
    p_orb = info["model"].p_orb
    ev = pdat.read_ev_series(info["dir"], info["prefix"])
    md = pdat.read_ev(info["modes"]) if info["modes"].exists() else {}
    if not ev:
        print("no .ev files")
        return 1
    t = ev["time"] / p_orb

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    ax = axes[0, 0]
    ax.plot(t, rel(ev["etot"]), label=r"$E_{\rm tot}$")
    ax.plot(t, rel(ev["angtot"]), label=r"$L$")
    ax.set_ylabel("relative change")
    ax.set_title("Conservation", fontweight="bold")
    ax.legend()

    ax = axes[0, 1]
    ax.plot(t, ev["ekin"], label=r"$E_{\rm kin}$")
    ax.plot(t, ev["etherm"], label=r"$E_{\rm th}$")
    ax.set_title("Energies", fontweight="bold")
    ax.legend()

    ax = axes[0, 2]
    if md:
        tm = md["time"] / p_orb
        ax.plot(tm, rel(md["mass"]), color="k")
        ax.set_ylabel(r"$M/M_0 - 1$")
    ax.set_title("Torus mass (accretion through r_acc)", fontweight="bold")

    ax = axes[1, 0]
    ax.semilogy(t, ev["alpha max"], label=r"$\alpha_{\rm max}$")
    ax.semilogy(t, ev["rmsmach"], label="rms Mach")
    ax.set_title("Shock capturing", fontweight="bold")
    ax.legend()

    ax = axes[1, 1]
    ax.plot(t, ev["rho max"] / ev["rho max"][0], label=r"$\rho_{\max}/\rho_{\max}(0)$")
    ax.set_title("Peak density", fontweight="bold")
    ax.legend()

    ax = axes[1, 2]
    if md:
        r_cm = np.hypot(md["xcm"], md["ycm"])
        ax.semilogy(tm, np.maximum(r_cm, 1e-300), color="C3")
        ax.set_ylabel(r"$|\mathbf{r}_{\rm cm}|$")
    ax.set_title("Centre of mass offset (m = 1)", fontweight="bold")

    for ax in axes.flat:
        ax.set_xlabel(r"$t\;[P_{\rm orb}]$")
        ax.grid(True, alpha=0.3)
    fig.suptitle("History", fontsize=15, fontweight="bold")
    fig.text(0.012, 0.988, "\n".join(pdat.info_lines(info)), ha="left", va="top",
             fontsize=8.5, bbox=dict(boxstyle="round", fc="white", ec="#b8b8b4", alpha=0.9))
    fig.tight_layout()
    out_dir = Path(a.out) if a.out else info["dir"] / "figs"
    out_dir.mkdir(parents=True, exist_ok=True)
    f = out_dir / "history.png"
    fig.savefig(f, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f)
    return 0


if __name__ == "__main__":
    sys.exit(main())
