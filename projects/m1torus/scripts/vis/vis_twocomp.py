#!/usr/bin/env python3
"""Stars and gas of a two-component run (collisionless stars carrying the mode, colliding
gas clouds): the m = 1 amplitude of each with the gas noise level 2/sqrt(N), the phase of
the gas against the stars, the vertical dispersion of the gas, the radius holding 90% of
the bound gas, and the gas clouds within the torus region (r < rmax) over time.

    vis_twocomp.py --run DIR [--rmax 1.6] [--every 1]
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "obs"))
import labels as L  # noqa: E402
import twocomp as TC  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--rmax", type=float, default=None,
                    help="default 1.6 with inflow (leave out the clouds falling in), else all gas")
    ap.add_argument("--every", type=int, default=1)
    a = ap.parse_args()
    run = Path(a.run).resolve()
    if a.rmax is None:
        feed = float(str(L.pio.read_setup(run / "m1torus.in").get("feed_mdot", 0) or 0).split()[0])
        a.rmax = 1.6 if feed > 0 else np.inf
    rows = [TC.row(p, a.rmax) for p in sorted(run.glob("m1torus_[0-9][0-9][0-9][0-9][0-9]"))[::a.every]]
    t = np.array([r["t"] for r in rows])
    g = lambda c, k: np.array([r[c][k] for r in rows])
    dph = (g("gas", "phase") - g("stars", "phase") + 180) % 360 - 180

    fig, ax = plt.subplots(5, 1, figsize=(8.5, 10.5), sharex=True)
    ax[0].plot(t, g("stars", "A1"), color="0.2", label="stars")
    ax[0].plot(t, g("gas", "A1"), color="#d81bd8", label="gas clouds")
    ax[0].fill_between(t, 0, g("gas", "noise"), color="#d81bd8", alpha=0.15, lw=0, label=r"gas noise $2/\sqrt{N}$")
    ax[0].set(ylabel=r"$A_1/\Sigma_0$", ylim=(0, None))
    ax[0].legend(fontsize=8, ncol=3, loc="lower left")
    ax[1].plot(t, dph, color="#d81bd8")
    ax[1].axhline(0, color="0.6", lw=0.6)
    ax[1].set(ylabel=r"$\Phi_1^{\rm gas} - \Phi_1^{\rm stars}$ [deg]")
    ax[2].plot(t, g("gas", "sz"), color="#d81bd8")
    ax[2].plot(t, g("stars", "sz"), color="0.2")
    ax[2].set(ylabel=r"$\sigma_z / v_{\rm K}$ ($0.5 < R < 1.5$)", ylim=(0, None))
    ax[3].plot(t, g("gas", "R90b"), color="#d81bd8")
    ax[3].set(ylabel=r"$R_{90}$ of bound gas $[R_{\rm tor}]$")
    ax[4].plot(t, g("gas", "n"), color="#d81bd8")
    ax[4].set(ylabel=(rf"gas clouds, $r < {a.rmax:g}\,R_{{\rm tor}}$" if np.isfinite(a.rmax) else "gas clouds"),
              xlabel=r"$t / P_{\rm orb}$")
    name, pars, suffix = L.model_parts(run)
    fig.suptitle("Stars and Gas of the Torus\n" + L.wrapped(f"{name} ({', '.join(pars.values())}){suffix}", 95),
                 fontsize=9.5)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    out = run / "figs"
    out.mkdir(exist_ok=True)
    fig.savefig(out / "twocomp.png", dpi=120)
    print(out / "twocomp.png")


if __name__ == "__main__":
    main()
