#!/usr/bin/env python3
"""B35: how much of the broad emission at the NGC 613 nucleus is H13CN(4-3), which falls
at +396 km/s in the CO(3-2) cube. If the same gas emits both lines, the spectrum is
S(v) = C(v) + xi C(v - 396); xi follows from v > +300, where CO itself is weak, with C
taken from the spectrum. Apertures: the nucleus (1 and 2 beams about the kinematic
centre of B34) and the torus annulus.

    nucleus_spectrum.py [--out DIR]
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
sys.path.insert(0, str(HERE.parent / "vis"))
import ngc613_cube as NC  # noqa: E402

DV = 299792.458 * (345.79599 - 345.33976) / 345.33976    # H13CN(4-3) in the CO(3-2) frame
CENTRE = (0.0, 0.004)                                      # kinematic centre, B34 [arcsec]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=str(HERE.parents[3] / "results" / "m1torus" / "B" / "ngc613_cube"))
    a = ap.parse_args()
    cube, v, E, N, pix, beam = NC.load()
    cube = np.nan_to_num(cube.astype(float))
    EE, NN = np.meshgrid(E, N)
    r = np.hypot(EE - CENTRE[0], NN - CENTRE[1])
    print(f"H13CN(4-3) at {DV:+.1f} km/s")
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    for axi, (lab, sel) in zip(ax, [("nucleus, r < 1 beam", r < beam), ("nucleus, r < 2 beams", r < 2 * beam),
                                    ("torus, 0.03 < r < 0.1 arcsec", (r > 0.03) & (r < 0.1))]):
        S = cube[:, sel].sum(1)
        Ssh = np.interp(v - DV, v, S, left=0, right=0)
        k = (v > 300) & (v < 700)
        xi = float((S[k] * Ssh[k]).sum() / (Ssh[k] ** 2).sum())
        off = (v < -500) | (v > 720)
        noise = np.std(S[off])
        co = S - xi * Ssh
        core = (v > -300) & (v < 300)
        print(f"{lab:30s} xi {xi:.3f}; residual over +300..+700 {np.std(S[k] - xi * Ssh[k]) / noise:.2f} x noise; "
              f"CO at +230..+300 {co[(v > 230) & (v < 300)].sum() / co[core].sum() * 100:5.1f}% "
              f"(before removal {S[(v > 230) & (v < 300)].sum() / S[core].sum() * 100:5.1f}%), "
              f"CO below -200 {co[(v < -200) & (v > -450)].sum() / co[core].sum() * 100:5.1f}%")
        axi.step(v, S * 1e3, where="mid", lw=0.8, color="0.3", label="CO(3-2) cube")
        axi.step(v, xi * Ssh * 1e3, where="mid", lw=0.8, color="C1",
                 label=rf"H13CN(4-3): {xi:.2f} $\times$ CO shifted by {DV:+.0f} km/s")
        axi.step(v, co * 1e3, where="mid", lw=0.9, color="C0", label="CO after removal")
        axi.axhline(0, color="k", lw=0.4)
        axi.set(title=f"NGC 613, {lab}", xlabel=r"$v - 1481$ [km/s]", ylabel="flux [mJy]")
        axi.legend(fontsize=7)
    fig.suptitle("ALMA 2019.1.00273.S archive cube, CO(3-2) window", fontsize=9)
    fig.tight_layout()
    out = Path(a.out) / "ngc613_nucleus_spectrum.png"
    fig.savefig(out, dpi=110)
    print(out)


if __name__ == "__main__":
    main()
