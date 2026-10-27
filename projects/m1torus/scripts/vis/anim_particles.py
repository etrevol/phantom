#!/usr/bin/env python3
"""Particle movie of a m1torus run, face-on and edge-on, in the style of Bannikova et al.

Each particle is coloured by its own density, rho = m (h_fact/h)^3 from its smoothing
length, on the logarithmic scale of the paper's density maps (1e-4 to 1e-1 M_c R_tor^-3);
dense particles are drawn last. Marked: the central mass, the torus barycentre and the
m = 1 mode (arrow along its phase, length proportional to A_1/Sigma_0). A strip below
follows A_1/Sigma_0 and |r_c| with a cursor at the frame time.

    anim_particles.py [--run DIR] [--every 1] [--fps 12] [--still T] [--corotate]

--corotate turns every frame by -Phi_1, so the mode stands still and its shape and the
orbits streaming through it are visible.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "common"))
import phantom_io as pio  # noqa: E402

P_ORB = 2.0 * np.pi
HFACT = 1.2
LIM = 2.0
NORM = LogNorm(1e-4, 1e-1)
CMAP = "turbo"
ARROW_FULL = 1.6        # arrow length at the paper's saturated A_1/Sigma_0 = 0.368
A1_PAPER = 0.368

plt.rcParams.update({"font.size": 10, "axes.labelsize": 11, "axes.linewidth": 0.8,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})


def frame_data(path):
    d = pio.read_dump(path)
    g, s, hdr = d["gas"], d["sinks"], d["header"]
    keep = g["h"] > 0
    return dict(t=hdr["time"] / P_ORB, x=g["x"][keep], y=g["y"][keep], z=g["z"][keep],
                h=g["h"][keep],
                xc=(s["x"][0], s["y"][0], s["z"][0]) if s else (0.0, 0.0, 0.0))


def build_figure(setup, n, tm, a1, rc, corotate=False):
    fig = plt.figure(figsize=(12.5, 7.6))
    ax_f = fig.add_axes([0.065, 0.31, 0.39, 0.62])
    ax_e = fig.add_axes([0.505, 0.31, 0.39, 0.62])
    cax = fig.add_axes([0.905, 0.40, 0.012, 0.44])
    ax_s = fig.add_axes([0.065, 0.075, 0.83, 0.14])

    for ax, ylab in ((ax_f, r"$y/R_{\rm tor}$"), (ax_e, r"$z/R_{\rm tor}$")):
        ax.set_xlim(-LIM, LIM)
        ax.set_ylim(-LIM, LIM)
        ax.set_aspect("equal")
        ticks = np.arange(-1.5, 1.6, 0.5)      # no labels at the corners, where they collide
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xlabel(r"$x/R_{\rm tor}$")
        ax.set_ylabel(ylab, labelpad=2)
    ax_f.text(0.03, 0.97, r"face-on, frame of $m=1$" if corotate else "face-on",
              transform=ax_f.transAxes, va="top", fontsize=10)
    ax_e.text(0.03, 0.97, "edge-on", transform=ax_e.transAxes, va="top", fontsize=10)

    # point size so the panels fill alike at any N (0.6 at the paper's 128k)
    size = min(3.0, 0.6 * (131072 / n) ** 0.5)
    sc_f = ax_f.scatter([], [], s=size, c=[], cmap=CMAP, norm=NORM, lw=0, rasterized=True)
    sc_e = ax_e.scatter([], [], s=size, c=[], cmap=CMAP, norm=NORM, lw=0, rasterized=True)
    cb = fig.colorbar(sc_f, cax=cax)
    cb.set_label(r"$\rho\;[M_c\,R_{\rm tor}^{-3}]$")
    cb.ax.tick_params(direction="out")

    bh_f, = ax_f.plot([], [], "o", color="k", ms=5, zorder=6)
    bh_e, = ax_e.plot([], [], "o", color="k", ms=5, zorder=6)
    tb_f, = ax_f.plot([], [], "x", color="k", ms=7, mew=1.6, zorder=6)
    handles = [Line2D([], [], ls="", marker="o", color="k", ms=5, label="central mass"),
               Line2D([], [], ls="", marker="x", color="k", ms=7, mew=1.6,
                      label="torus barycentre"),
               Line2D([], [], color="crimson", lw=2, marker=">", ms=6,
                      label=r"$m=1$: phase $\Phi_1$, length $\propto A_1/\Sigma_0$")]
    ax_f.legend(handles=handles, loc="lower left", fontsize=8.5, frameon=True,
                framealpha=0.85, edgecolor="0.7", handlelength=1.6, borderpad=0.5)

    ax_s.plot(tm, a1, color="C0", lw=1.0)
    ax_s.set_ylabel(r"$A_1/\Sigma_0$", color="C0")
    ax_s.tick_params(axis="y", colors="C0", right=False)
    ax_s.set_xlabel(r"$t/P_{\rm orb}$", labelpad=1)
    ax_s.set_xlim(tm[0], tm[-1])
    ax_s.set_ylim(0, max(0.4, 1.1 * np.nanmax(a1)))
    ax_r = ax_s.twinx()
    ax_r.plot(tm, rc, color="C1", lw=1.0)
    ax_r.set_ylabel(r"$|r_c|/R_{\rm tor}$", color="C1")
    ax_r.tick_params(axis="y", colors="C1", direction="in")
    ax_r.set_ylim(0, max(0.03, 1.1 * np.nanmax(rc)))
    cursor = ax_s.axvline(tm[0], color="k", lw=0.9)

    params = (rf"$N = {n:,}".replace(",", r"\,") + r",\;\;M_{\rm tor} = "
              + f"{setup['M_tor']:g}" + r"\,M_c,\;\;e_{\max} = " + f"{setup['e_max']:g}"
              + r",\;\;i_{\max} = " + f"{setup['i_max']:g}" + r"^\circ$")
    fig.text(0.065, 0.955, params, ha="left", va="center", fontsize=11)
    # the clock in three fixed pieces, digits of equal width, so nothing moves
    fig.text(0.80, 0.955, r"$t\,=$", ha="left", va="center", fontsize=11)
    clock = fig.text(0.868, 0.955, "", ha="right", va="center", fontsize=11)
    fig.text(0.873, 0.955, r"$P_{\rm orb}$", ha="left", va="center", fontsize=11)
    return dict(fig=fig, ax_f=ax_f, sc_f=sc_f, sc_e=sc_e, bh_f=bh_f, bh_e=bh_e,
                tb_f=tb_f, cursor=cursor, clock=clock, arrow=None)


def draw(art, path, pmass, tm, a1, m1, corotate=False):
    fr = frame_data(path)
    k = int(np.argmin(np.abs(tm - fr["t"])))
    if corotate:
        cr, sr = np.cos(-m1["Phi1"][k]), np.sin(-m1["Phi1"][k])
        fr["x"], fr["y"] = cr * fr["x"] - sr * fr["y"], sr * fr["x"] + cr * fr["y"]
        xc = fr["xc"]
        fr["xc"] = (cr * xc[0] - sr * xc[1], sr * xc[0] + cr * xc[1], xc[2])
    rho = pmass * (HFACT / fr["h"]) ** 3
    order = np.argsort(rho)
    x, y, z, c = fr["x"][order], fr["y"][order], fr["z"][order], rho[order]
    art["sc_f"].set_offsets(np.c_[x, y])
    art["sc_f"].set_array(c)
    art["sc_e"].set_offsets(np.c_[x, z])
    art["sc_e"].set_array(c)
    art["bh_f"].set_data([fr["xc"][0]], [fr["xc"][1]])
    art["bh_e"].set_data([fr["xc"][0]], [fr["xc"][2]])
    xtb, ytb = m1["xtb"][k], m1["ytb"][k]
    if corotate:
        xtb, ytb = cr * xtb - sr * ytb, sr * xtb + cr * ytb
    art["tb_f"].set_data([xtb], [ytb])
    if art["arrow"] is not None:
        art["arrow"].remove()
    length = ARROW_FULL * a1[k] / A1_PAPER
    phi = 0.0 if corotate else m1["Phi1"][k]
    art["arrow"] = art["ax_f"].annotate(
        "", xy=(length * np.cos(phi), length * np.sin(phi)), xytext=(0, 0), zorder=7,
        arrowprops=dict(arrowstyle="-|>", color="crimson", lw=2, mutation_scale=14))
    art["cursor"].set_xdata([fr["t"], fr["t"]])
    art["clock"].set_text(rf"$\mathtt{{{fr['t']:.1f}}}$")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", default=".")
    ap.add_argument("--every", type=int, default=1, help="use every n-th dump")
    ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--still", type=float, default=None, help="one PNG at this orbit")
    ap.add_argument("--corotate", action="store_true", help="frame turning with the mode")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)

    run = Path(a.run).resolve()
    dumps = sorted(run.glob("m1torus_[0-9][0-9][0-9][0-9][0-9]"))[::a.every]
    m1 = pio.read_ev(run / "m1torus_m1.ev")
    setup = pio.read_setup(run / "m1torus.setup")
    n = int(setup["np"])
    pmass = float(setup["M_tor"]) / n
    tm = m1["time"] / P_ORB
    a1 = m1["A1"] / m1["sigma0"]
    rc = np.hypot(m1["xc"], m1["yc"])
    art = build_figure(setup, n, tm, a1, rc, a.corotate)

    out = Path(a.out) if a.out else run / "figs"
    out.mkdir(parents=True, exist_ok=True)
    if a.still is not None:
        times = np.array([pio.read_dump(p)["header"]["time"] / P_ORB for p in dumps])
        draw(art, dumps[int(np.argmin(np.abs(times - a.still)))], pmass, tm, a1, m1,
             a.corotate)
        f = out / f"particles{'_corot' if a.corotate else ''}_t{a.still:06.1f}.png"
        art["fig"].savefig(f, dpi=170)
        print(f)
        return 0

    from matplotlib.animation import FFMpegWriter
    f = out / ("anim_particles_corot.mp4" if a.corotate else "anim_particles.mp4")
    writer = FFMpegWriter(fps=a.fps, bitrate=8000)
    with writer.saving(art["fig"], str(f), dpi=130):
        for p in dumps:
            draw(art, p, pmass, tm, a1, m1, a.corotate)
            writer.grab_frame()
    print(f"{len(dumps)} frames -> {f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
