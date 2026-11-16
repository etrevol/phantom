#!/usr/bin/env python3
"""Movie of the eccentric warped disc of a whole-cube fit (kinematic model, NGC 1672): the
clouds move along their fixed Keplerian ellipses for one orbit of the outer edge, seen
face-on in the inner disc plane and on the sky at the fitted view. The ellipses do not
precess (no self-gravity in this model), so the lopsided pattern stands while the clouds
stream through it: crowding at the aligned apocentres makes the bright side.

    anim_eccdisc.py --fit FIT_DIR [--frames 72] [--fps 12] [--out FILE]
"""
import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "obs"))
import fit_torus as T  # noqa: E402
import mockobs as mo  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--fit", required=True)
    ap.add_argument("--galaxy", default="ngc1672")
    ap.add_argument("--frames", type=int, default=72)
    ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--n", type=int, default=8000, help="clouds drawn")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    fit = json.load(open(Path(a.fit) / "fit.json"))
    b = dict(T.DEFAULT, **fit["best"])
    T.GAL = a.galaxy
    t = replace(T.target(), R_pc=b["R_pc"], inc=b["inc"], M_bh=10 ** b["logM"], pa=b["pa"])
    pcas = mo.pc_per_arcsec(t)
    v_kms = np.sqrt(4.30091e-3 * t.M_bh / b["R_pc"])                  # G in pc (km/s)^2 / Msun
    years = 2 * np.pi * b["R_pc"] * 3.0857e13 / v_kms / 3.156e7        # one orbit at the outer edge
    args = (b["a_in"], b["e0"], b["e1"], b["dw"], b["warp_i"], b["warp_n"], b["i_max"])

    fig, (af, asky) = plt.subplots(1, 2, figsize=(13, 6.6))
    fig.subplots_adjust(left=0.06, right=0.90, top=0.80, bottom=0.10, wspace=0.25)
    sf = af.scatter([], [], s=1.5, c=[], cmap="viridis", vmin=0, vmax=0.5, lw=0, rasterized=True)
    ss = asky.scatter([], [], s=1.5, c=[], cmap="viridis", vmin=0, vmax=0.5, lw=0, rasterized=True)
    af.set(xlim=(-1.6, 1.6), ylim=(-1.6, 1.6), aspect="equal", xlabel=r"$x/R$", ylabel=r"$y/R$",
           title="face-on, in the plane of the inner disc")
    lim = 1.25 * b["R_pc"] / pcas
    asky.set(xlim=(lim, -lim), ylim=(-lim, lim), aspect="equal", xlabel="RA offset [arcsec]",
             ylabel="Dec offset [arcsec]", title=f"on the sky as fitted: $i = {b['inc']:.0f}^\\circ$, PA ${b['pa'] % 360:.0f}^\\circ$")
    for ax in (af, asky):
        ax.plot(0, 0, "o", color="k", ms=5, zorder=5)
    cb = fig.colorbar(ss, cax=fig.add_axes([0.92, 0.18, 0.015, 0.55]))
    cb.set_label("eccentricity of the orbit")
    fig.text(0.5, 0.975, "NGC 1672: the eccentric warped disc of the whole-cube fit (kinematic model)",
             ha="center", va="top", fontsize=12)
    fig.text(0.5, 0.935, rf"$a = {b['a_in']:.2f}-1\,R$, $R = {b['R_pc']:.1f}$ pc, $e = {b['e0']:.2f} \to {b['e1']:.2f}$, "
             rf"twist ${b['dw']:.0f}^\circ$, warp ${b['warp_i']:.0f}^\circ$ at ${b['warp_n']:.0f}^\circ$, "
             rf"$i_{{\max}} = {b['i_max']:.0f}^\circ$, $M_{{\rm BH}} = {10 ** b['logM'] / 1e7:.2f}\times10^7\,M_\odot$; "
             "clouds move on fixed ellipses, black dot: the black hole", ha="center", va="top", fontsize=9)
    clock = fig.text(0.98, 0.975, "", ha="right", va="top", fontsize=11)
    from matplotlib.animation import FFMpegWriter
    out = Path(a.out) if a.out else Path(a.fit) / "anim_eccdisc.mp4"
    writer = FFMpegWriter(fps=a.fps, bitrate=5000)
    with writer.saving(fig, str(out), dpi=110):
        for k in range(a.frames):
            time = 2 * np.pi * k / a.frames
            p = T.eccdisc(*args, n=a.n, time=time)
            x, v = p["x"], p["v"]
            r = np.linalg.norm(x, axis=1)
            ecc = np.linalg.norm(np.cross(v, np.cross(x, v)) - x / r[:, None], axis=1)
            sf.set_offsets(x[:, :2])
            sf.set_array(ecc)
            east, north, _, _ = mo.to_sky(x, v, b["inc"], b["pa"], b["az"])
            ss.set_offsets(np.c_[east * b["R_pc"] / pcas, north * b["R_pc"] / pcas])
            ss.set_array(ecc)
            clock.set_text(rf"$t = {k / a.frames:.2f}\,P(R)$ = {k / a.frames * years / 1e3:.0f} kyr")
            writer.grab_frame()
    print(out)


if __name__ == "__main__":
    main()
