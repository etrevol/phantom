#!/usr/bin/env python3
"""Movie of a run as ALMA would see it, next to the NGC 613 archive maps.

Every dump is observed from one fixed direction (the view of a whole-cube fit, so the
mode precesses across the sky): CO(3-2) moments 0 and 1 with the archive beam, channels
and noise, next to the archive cube processed the same way and centred on the fitted
black hole; a strip follows A_1/Sigma_0 of the run with a cursor.

    anim_mockobs.py --run RUN_DIR --fit FIT_DIR [--every 1] [--fps 8] [--out FILE]
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
sys.path.insert(0, str(HERE.parents[2] / "common"))
import compare_obs as CO  # noqa: E402
import fit_torus as T  # noqa: E402
import labels as L  # noqa: E402
import mockobs as mo  # noqa: E402
import phantom_io as pio  # noqa: E402

LIM = 0.12      # half field [arcsec]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--fit", required=True, help="fit_torus output directory with fit.json (the view)")
    ap.add_argument("--every", type=int, default=1)
    ap.add_argument("--fps", type=int, default=8)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    run = Path(a.run).resolve()
    b = dict(T.DEFAULT, **json.load(open(Path(a.fit) / "fit.json"))["best"])
    T.GAL = "ngc613"
    T.EDGES = (0.01, 0.025, 0.04, 0.06, 0.08, LIM)
    t = replace(T.target(), R_pc=b["R_pc"], inc=b["inc"], M_bh=10 ** b["logM"], pa=b["pa"])
    dumps = sorted(run.glob("m1torus_[0-9][0-9][0-9][0-9][0-9]"))[::a.every]
    # one fixed observer: the fitted azimuth relative to the mode phase of the fitted dump
    p_fit = mo.load(json.load(open(Path(a.fit) / "fit.json"))["model"])
    phi = mo.m1_phase(p_fit["x"], p_fit["m"]) + b["az"]

    d = np.load(T.ROOT / "results/m1torus/B/ngc613_cube/ngc613_mask_moments.npz")
    de, dn = d["east"] - b["dx"], d["north"] - b["dy"]
    m1 = pio.read_ev(run / "m1torus_m1.ev")
    tm, a1 = m1["time"] / (2 * np.pi), m1["A1"] / m1["sigma0"]

    fig = plt.figure(figsize=(10, 10.6))
    axs = [[fig.add_axes([0.07 + 0.45 * j, 0.555 - 0.385 * i, 0.36, 0.33]) for j in range(2)] for i in range(2)]
    ax_s = fig.add_axes([0.07, 0.05, 0.81, 0.055])
    ims = []
    for i, k in enumerate(("m0", "m1")):
        cm, vmin, vmax = CO.PAL[k]
        data = d[k] - (b["dv"] if k == "m1" else 0.)
        ax = axs[i][1]
        ax.pcolormesh(de, dn, np.where(np.isfinite(data), data, np.nan), cmap=cm, vmin=vmin, vmax=vmax,
                      shading="auto", rasterized=True)
        ax.set_title(f"NGC 613, moment {i}: {CO.ALMA}", fontsize=8)
        im = axs[i][0].imshow(np.full((2, 2), np.nan), cmap=cm, vmin=vmin, vmax=vmax, origin="lower",
                              extent=(LIM, -LIM, -LIM, LIM), interpolation="nearest")
        ims.append(im)
        axs[i][0].set_title(f"Synthetic CO(3-2), moment {i}", fontsize=9)
        cax = fig.add_axes([0.89, 0.585 - 0.385 * i, 0.014, 0.27])
        fig.colorbar(im, cax=cax).set_label("moment 0 [Jy/beam km/s]" if i == 0 else "moment 1 [km/s]")
        for ax in axs[i]:
            ax.set(xlim=(LIM, -LIM), ylim=(-LIM, LIM), aspect="equal")
            ax.set_xlabel("RA offset from the black hole [arcsec]", fontsize=8)
            ax.set_ylabel("Dec offset [arcsec]", fontsize=8)
            ax.plot(0, 0, "+", color="k", ms=9, mew=1.5)
    ax_s.plot(tm, a1, color="C0", lw=1.0)
    ax_s.set(xlim=(tm[0], tm[-1]), ylim=(0, None), xlabel=r"$t\,[P_{\rm orb}]$", ylabel=r"$A_1/\Sigma_0$")
    cursor = ax_s.axvline(tm[0], color="k", lw=0.9)
    name, pars, suffix = L.model_parts(run)
    fig.text(0.5, 0.985, "NGC 613 in CO(3-2) at 0.014 arcsec: a synthetic observation over time",
             ha="center", va="top", fontsize=11)
    fig.text(0.5, 0.962, L.wrapped(f"{name} ({', '.join(pars.values())})", 110), ha="center", va="top", fontsize=8.5)
    fig.text(0.5, 0.925, rf"view of the whole-cube fit: $R_{{\rm tor}} = {b['R_pc']:.1f}$ pc, $i = {b['inc']:.0f}^\circ$, "
             rf"$M_{{\rm BH}} = {10 ** b['logM'] / 1e7:.2f}\times10^7\,M_\odot$, PA ${b['pa'] % 360:.0f}^\circ$; "
             "the observer is fixed, the mode turns", ha="center", va="top", fontsize=8.5)
    clock = fig.text(0.93, 0.985, "", ha="right", va="top", fontsize=11)

    from matplotlib.animation import FFMpegWriter
    out = Path(a.out) if a.out else run / "figs" / "anim_mockobs.mp4"
    out.parent.mkdir(parents=True, exist_ok=True)
    writer = FFMpegWriter(fps=a.fps, bitrate=6000)
    with writer.saving(fig, str(out), dpi=110):
        for k, f in enumerate(dumps):
            p = mo.load(f)
            cb = T.fast_cube(p, t, phi)
            m0, mm1, _ = mo.moments(cb, t, rms=cb["rms"])
            xs = cb["xs"]
            for im, m in zip(ims, (m0, mm1)):
                im.set_data(m)
                im.set_extent((xs[0], xs[-1], xs[0], xs[-1]))
            time = p["time"] / (2 * np.pi)
            cursor.set_xdata([time, time])
            clock.set_text(rf"$t = {time:5.1f}\,P_{{\rm orb}}$")
            writer.grab_frame()
            print(f"frame {k + 1}/{len(dumps)} t {time:.1f}", flush=True)
    print(out)


if __name__ == "__main__":
    main()
