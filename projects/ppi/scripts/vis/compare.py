#!/usr/bin/env python3
"""Several PPI torus runs on one figure: torus mass and |A_m|/A_0 against time.

    compare.py RUN_A RUN_B ... [--m 1,2] [--labels a,b] [--out FILE]

Labels default to the part of each directory name after the timestamp.
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
from vismodes import COLORS, LINEAR_CEILING  # noqa: E402

STYLES = ["-", "--", ":", "-."]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--m", default="1,2,4")
    ap.add_argument("--labels", default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)

    ms = [int(s) for s in a.m.split(",")]
    labels = a.labels.split(",") if a.labels else [
        Path(r).resolve().name.split("-", 2)[-1] for r in a.runs]
    fig, (axm, axa) = plt.subplots(1, 2, figsize=(16, 6.5))
    print(f"{'run':32s} {'orbits':>7s} {'dM/M':>9s} {'dM/M per orbit':>15s}")
    for k, (run, lab) in enumerate(zip(a.runs, labels)):
        info = pdat.run_info(run)
        ev = pdat.read_ev(info["modes"])
        t = ev["time"] / info["model"].p_orb
        dm = ev["mass"] / ev["mass"][0] - 1.0
        ls = STYLES[k % len(STYLES)]
        axm.plot(t, dm, color="k", ls=ls, lw=1.6, label=lab)
        for j, m in enumerate(ms):
            axa.semilogy(t, ev[f"A{m}"] / ev["A0"], color=COLORS[j % len(COLORS)], ls=ls,
                         lw=1.4, label=f"{lab}, m = {m}")
        print(f"{lab:32s} {t[-1]:7.2f} {dm[-1]:9.2e} {dm[-1] / max(t[-1], 1e-9):15.2e}")
    axm.set_ylabel(r"$M/M_0-1$")
    axm.set_title("Torus mass", fontweight="bold")
    axa.axhline(LINEAR_CEILING, color="k", lw=0.8, ls="--", alpha=0.6)
    axa.set_ylabel(r"$|A_m|/A_0$")
    axa.set_title("Azimuthal modes", fontweight="bold")
    for ax in (axm, axa):
        ax.set_xlabel(r"$t\;[P_{\rm orb}]$")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, ncol=2 if ax is axa else 1)
    fig.tight_layout()
    out = Path(a.out) if a.out else Path(a.runs[0]).resolve().parent / "compare.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
