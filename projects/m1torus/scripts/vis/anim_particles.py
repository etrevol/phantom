#!/usr/bin/env python3
"""Particle movie of a m1torus run, in the style of Bannikova et al.: face-on and edge-on.

Each particle is coloured by its own density, rho = m (h_fact/h)^3 from its smoothing
length, on the logarithmic scale of the paper's density maps (1e-4 to 1e-1 in units of
M_c / R_tor^3); dense particles are drawn last. The central mass is the black dot, the
torus barycentre a cross, and the arrow gives the phase and amplitude of the m = 1 mode.
A strip under the panels follows A_1/Sigma_0 and |r_c| with a marker at the frame time.

    anim_particles.py [--run DIR] [--every 1] [--fps 12] [--still T]
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

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "common"))
import phantom_io as pio  # noqa: E402

P_ORB = 2.0 * np.pi
HFACT = 1.2
LIM = 2.0
NORM = LogNorm(1e-4, 1e-1)
CMAP = "turbo"


def frame_data(path):
    d = pio.read_dump(path)
    g, s, hdr = d["gas"], d["sinks"], d["header"]
    keep = g["h"] > 0
    x, y, z = (g[k][keep] for k in ("x", "y", "z"))
    return dict(t=hdr["time"] / P_ORB, x=x, y=y, z=z, h=g["h"][keep],
                xc=(s["x"][0], s["y"][0], s["z"][0]) if s else (0.0, 0.0, 0.0))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", default=".")
    ap.add_argument("--every", type=int, default=1, help="use every n-th dump")
    ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--still", type=float, default=None, help="one PNG at this orbit")
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

    fig = plt.figure(figsize=(14, 8.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[4.3, 1.2], hspace=0.28, wspace=0.18)
    ax_f, ax_e = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    ax_s = fig.add_subplot(gs[1, :])
    for ax, lab in ((ax_f, "y"), (ax_e, "z")):
        ax.set_xlim(-LIM, LIM)
        ax.set_ylim(-LIM, LIM)
        ax.set_aspect("equal")
        ax.set_xlabel("x [NB]")
        ax.set_ylabel(f"{lab} [NB]")
        ax.grid(True, alpha=0.25, ls=":")
    ax_f.set_title("face-on", fontsize=11)
    ax_e.set_title("edge-on", fontsize=11)
    # point size so that the panels fill alike at any N (0.6 at the paper's 128k)
    size = min(3.0, 0.6 * (131072 / n) ** 0.5)
    sc_f = ax_f.scatter([], [], s=size, c=[], cmap=CMAP, norm=NORM, lw=0, rasterized=True)
    sc_e = ax_e.scatter([], [], s=size, c=[], cmap=CMAP, norm=NORM, lw=0, rasterized=True)
    bh_f, = ax_f.plot([], [], "ko", ms=5, zorder=5)
    bh_e, = ax_e.plot([], [], "ko", ms=5, zorder=5)
    tb_f, = ax_f.plot([], [], "kx", ms=7, mew=1.5, zorder=5)
    arrow = [None]
    cb = fig.colorbar(sc_f, ax=[ax_f, ax_e], orientation="horizontal", location="top",
                      fraction=0.04, pad=0.08, aspect=50)
    cb.set_label(r"density [$M_c/R_{\rm tor}^3$]")

    ax_s.plot(tm, a1, color="C0", lw=1, label=r"$A_1/\Sigma_0$")
    ax_s.set_ylabel(r"$A_1/\Sigma_0$", color="C0")
    ax_s.set_xlabel(r"$t\;[P_{\rm orb}]$")
    ax_s.set_xlim(tm[0], tm[-1])
    ax_s.grid(True, alpha=0.3)
    ax_r = ax_s.twinx()
    ax_r.plot(tm, rc, color="C1", lw=1)
    ax_r.set_ylabel(r"$|r_c|$", color="C1")
    mark = ax_s.axvline(0, color="k", lw=1)
    title = fig.suptitle("", fontsize=13, fontweight="bold")

    def draw(path):
        fr = frame_data(path)
        rho = pmass * (HFACT / fr["h"]) ** 3
        order = np.argsort(rho)
        x, y, z, c = fr["x"][order], fr["y"][order], fr["z"][order], rho[order]
        sc_f.set_offsets(np.c_[x, y])
        sc_f.set_array(c)
        sc_e.set_offsets(np.c_[x, z])
        sc_e.set_array(c)
        bh_f.set_data([fr["xc"][0]], [fr["xc"][1]])
        bh_e.set_data([fr["xc"][0]], [fr["xc"][2]])
        k = int(np.argmin(np.abs(tm - fr["t"])))
        tb_f.set_data([m1["xtb"][k]], [m1["ytb"][k]])
        if arrow[0] is not None:
            arrow[0].remove()
        phi, amp = m1["Phi1"][k], a1[k]
        arrow[0] = ax_f.annotate("", xy=(1.9 * amp * np.cos(phi) / 0.37, 1.9 * amp *
                                         np.sin(phi) / 0.37), xytext=(0, 0),
                                 arrowprops=dict(arrowstyle="-|>", color="crimson", lw=2))
        mark.set_xdata([fr["t"], fr["t"]])
        title.set_text(f"N = {n},  M_tor = {setup['M_tor']},  e_max = {setup['e_max']},  "
                       f"i_max = {setup['i_max']} deg     t = {fr['t']:6.1f} orbits")

    out = Path(a.out) if a.out else run / "figs"
    out.mkdir(parents=True, exist_ok=True)
    if a.still is not None:
        times = np.array([pio.read_dump(p)["header"]["time"] / P_ORB for p in dumps])
        draw(dumps[int(np.argmin(np.abs(times - a.still)))])
        f = out / f"particles_t{a.still:06.1f}.png"
        fig.savefig(f, dpi=150, bbox_inches="tight")
        print(f)
        return 0

    from matplotlib.animation import FFMpegWriter
    f = out / "anim_particles.mp4"
    writer = FFMpegWriter(fps=a.fps, bitrate=8000)
    with writer.saving(fig, str(f), dpi=120):
        for p in dumps:
            draw(p)
            writer.grab_frame()
    print(f"{len(dumps)} frames -> {f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
