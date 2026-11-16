#!/usr/bin/env python3
"""Channel-by-channel movie of a whole-cube fit: the archive cube, the best model with the
archive beam, channels and noise (times the fitted flux scale), and their difference,
centred on the fitted black hole, one velocity channel per frame.

    anim_channels.py --fit FIT_DIR [--galaxy ngc613] [--lim 0.12] [--fps 8] [--out FILE]
"""
import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.ndimage import map_coordinates  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "obs"))
import fit_torus as T  # noqa: E402
import labels as L  # noqa: E402
from galaxies import GALAXIES  # noqa: E402


def describe(fit, b):
    """The model in words, as on the fit figures."""
    src = fit["model"]
    if src == "eccdisc":
        s = (rf"Eccentric warped disc (kinematic): $a = {b['a_in']:.2f}-1\,R$, $e = {b['e0']:.2f} \to {b['e1']:.2f}$, "
             rf"warp ${b['warp_i']:.0f}^\circ$, $i_{{\max}} = {b['i_max']:.0f}^\circ$")
    elif src == "torus":
        s = (rf"Collisionless torus (phase-mixed, no mode): $R_0 = {b['R0']:.2f}$, $e_{{\max}} = {b['e_max']:.2f}$, "
             rf"$i_{{\max}} = {b['i_max']:.0f}^\circ$")
    else:
        s = L.model_label(Path(src).parent)
    if b.get("f_n", 0) > 0:
        s += rf"; inner component {100 * b['f_n']:.0f}% of the gas to ${b['R_n']:.2f}\,R_{{\rm tor}}$"
    return s + (rf"; $R_{{\rm tor}} = {b['R_pc']:.1f}$ pc, $i = {b['inc']:.0f}^\circ$, "
                rf"$M_{{\rm BH}} = {10 ** b['logM'] / 1e7:.2f}\times10^7\,M_\odot$")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--fit", required=True)
    ap.add_argument("--galaxy", default="ngc613", choices=sorted(GALAXIES))
    ap.add_argument("--lim", type=float, default=0.12, help="half field [arcsec]")
    ap.add_argument("--fps", type=int, default=8)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    fit = json.load(open(Path(a.fit) / "fit.json"))
    b = dict(T.DEFAULT, **fit["best"])
    T.GAL = a.galaxy
    T.EDGES = (0.01, a.lim)
    src = fit["model"] if fit["model"] in ("torus", "eccdisc") else str(Path(fit["model"]).parent)
    _, cb, t = T.model_metrics(b, src, "profiles+pv")      # noisy model cube, archive beam and channels
    model = 1e3 * fit.get("flux_scale", 1.0) * cb["cube"]      # mJy/beam
    xs, vs = cb["xs"], cb["vs"]
    T.observed("cube")
    dcb = T._data[T.GAL][1]
    dxs, dvs = dcb["xs"], dcb["vs"]
    # the archive cube on the model grid about the fitted black hole and systemic velocity
    E, N = np.meshgrid(xs, xs)
    rms = float(T.target().rms_mjy) * 1e-3
    name = GALAXIES[a.galaxy]["name"] if "name" in GALAXIES[a.galaxy] else a.galaxy.upper()
    vmax = 1e3 * float(np.nanpercentile(dcb["cube"], 99.9))

    fig, axs = plt.subplots(1, 3, figsize=(16.5, 5.6))
    fig.subplots_adjust(left=0.05, right=0.97, top=0.80, bottom=0.11, wspace=0.42)
    ims = []
    for ax, title, cmap, vr in zip(axs, (f"{name}: ALMA 2019.1.00273.S archive cube",
                                         r"best model $\times$ fitted flux scale, archive beam and noise",
                                         "difference (data - model)"),
                                   ("magma", "magma", "RdBu_r"), ((0, vmax), (0, vmax), (-vmax, vmax))):
        im = ax.imshow(np.zeros((len(xs), len(xs))), origin="lower", extent=(xs[0], xs[-1], xs[0], xs[-1]),
                       cmap=cmap, vmin=vr[0], vmax=vr[1], interpolation="nearest")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03).set_label("mJy/beam")
        ax.set(xlim=(a.lim, -a.lim), ylim=(-a.lim, a.lim), title=title,
               xlabel="RA offset from the black hole [arcsec]", ylabel="Dec offset [arcsec]")
        ax.title.set_fontsize(9)
        ax.plot(0, 0, "+", color="c", ms=9, mew=1.4)
        ims.append(im)
    fig.text(0.5, 0.975, f"{name} in CO(3-2), channel by channel: data, whole-cube fit, difference",
             ha="center", va="top", fontsize=12)
    fig.text(0.5, 0.935, L.wrapped(describe(fit, b), 150), ha="center", va="top", fontsize=8.5)
    clock = fig.text(0.97, 0.975, "", ha="right", va="top", fontsize=12)

    from matplotlib.animation import FFMpegWriter
    out = Path(a.out) if a.out else Path(a.fit) / "anim_channels.mp4"
    writer = FFMpegWriter(fps=a.fps, bitrate=6000)
    with writer.saving(fig, str(out), dpi=100):
        for k, v in enumerate(vs):
            iv = (v + b["dv"] - dvs[0]) / (dvs[1] - dvs[0])
            if iv < 0 or iv > len(dvs) - 1:
                continue
            idx = [np.full(E.shape, iv), (N + b["dy"] - dxs[0]) / (dxs[1] - dxs[0]),
                   (E + b["dx"] - dxs[0]) / (dxs[1] - dxs[0])]
            dat = 1e3 * map_coordinates(dcb["cube"], idx, order=1, mode="constant", cval=np.nan)
            for im, m in zip(ims, (dat, model[k], dat - model[k])):
                im.set_data(m)
            clock.set_text(rf"$v - v_{{\rm sys}} = {v:+5.0f}$ km/s")
            writer.grab_frame()
    print(out, f"rms {rms:.2e} Jy/beam")


if __name__ == "__main__":
    main()
