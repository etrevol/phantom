#!/usr/bin/env python3
"""Synthetic NGC 613 maps of a dump next to the published ALMA panels (Combes et al. 2025).

Ours on the left, theirs on the right, on the same field, colour scale and clip:
  <stem>_moments.png   CO(3-2) moments 0, 1, 2 against their Fig. 6
  <stem>_pv.png        PV along PA -20 deg against their Fig. 12, Keplerian curve of 3.4e7 Msun
  <stem>_spectrum.png  integrated spectrum against their Fig. 13 (FWHM 144 km/s, Table 2)
The reference crops come from ref_figures.py (results/m1torus/B/reference).

    compare_obs.py DUMP [--phi DEG] [--label TEXT] [--out DIR]
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm, ListedColormap  # noqa: E402

import mockobs as mo  # noqa: E402

ROOT = Path(__file__).resolve().parents[4]
REF = ROOT / "results" / "m1torus" / "B" / "reference"
# GILDAS rainbow as in the published maps: black, blue, cyan, green, yellow, red, magenta, white
GILDAS = LinearSegmentedColormap.from_list(
    "gildas", ["#000000", "#0000ff", "#00ffff", "#00ff00", "#ffff00", "#ff0000", "#ff00ff", "#ffffff"])
GILDAS.set_bad("white")


def panel_ref(ax, name, title):
    f = REF / f"{name}.png"
    if f.exists():
        ax.imshow(plt.imread(f))
    else:
        ax.text(0.5, 0.5, f"run ref_figures.py\n({name})", ha="center", transform=ax.transAxes)
    ax.set_axis_off()
    ax.set_title(title, fontsize=9)


def beam(ax, t, xs):
    from matplotlib.patches import Circle
    r = 0.5 * t.beam_arcsec
    ax.add_patch(Circle((xs[-1] - 3 * r - 0.02, xs[0] + 3 * r + 0.02), r, color="k"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dump")
    ap.add_argument("--phi", type=float, default=None)
    ap.add_argument("--label", default="")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    t = mo.NGC613
    p = mo.load(a.dump)
    phase = mo.m1_phase(p["x"], p["m"])
    phi = phase + 90.0 if a.phi is None else a.phi
    cb = mo.cube(p, t, phi)
    m0, m1, m2 = mo.moments(cb, t)
    xs, vs = cb["xs"], cb["vs"]
    h = mo.hole(m0, xs, t)
    out = Path(a.out or Path(a.dump).parent / "mockobs")
    out.mkdir(parents=True, exist_ok=True)
    stem = Path(a.dump).name
    tag = (a.label + ", " if a.label else "") + f"t = {p['time'] / (2 * np.pi):.0f} orbits, phi = {phi:.0f} deg"
    ext = [xs[0], xs[-1], xs[0], xs[-1]]

    # moments: masked where moment 0 is below the lowest colour (white there, as theirs)
    show = m0 >= 0.05
    fig, axs = plt.subplots(3, 2, figsize=(9, 13))
    specs = [(np.where(m0 >= 0.04, m0, np.nan), 0.05, 0.28, "moment 0 [Jy/beam km/s]", "combes25_fig6_m0"),
             (np.where(show, m1, np.nan), -100, 100, "moment 1 [km/s]", "combes25_fig6_m1"),
             (np.where(show, m2, np.nan), 20, 190, "moment 2 [km/s]", "combes25_fig6_m2")]
    for row, (img, lo, hi, lab, ref) in enumerate(specs):
        ax = axs[row, 0]
        im = ax.imshow(img, origin="lower", extent=ext, cmap=GILDAS, vmin=lo, vmax=hi,
                       interpolation="nearest")
        ax.set_xlim(xs[-1], xs[0])
        ax.plot(0, 0, "w+", ms=8, mew=1.5)
        beam(ax, t, xs)
        ax.set_xlabel("RA offset [arcsec]")
        ax.set_ylabel("Dec offset [arcsec]")
        fig.colorbar(im, ax=ax, fraction=0.046, label=lab)
        ax.set_title(f"synthetic, {lab.split(' [')[0]}", fontsize=9)
        panel_ref(axs[row, 1], ref, f"Combes et al. (2025) Fig. 6, {lab.split(' [')[0]}")
    axs[0, 0].plot(*_xy(h, t), "wx", ms=9, mew=2)
    fig.suptitle(f"NGC 613, CO(3-2) at 0.015 arcsec: {tag}\n"
                 f"hole {h['offset_pc']:.1f} pc from the AGN at PA {h['pa']:.0f} deg "
                 f"(observed 3 pc, NW)", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / f"{stem}_moments.png", dpi=110)
    plt.close(fig)

    # PV along PA -20 deg (their axis: 0.3 .. 0 arcsec, AGN at 0.15)
    pvd = mo.pv(cb, t, pa=t.pa, halfwidth=t.beam_arcsec)
    pos = 0.15 - xs
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    levels = [0.0002, 0.0003, 0.0004, 0.0005, 0.0006, 0.0008, 0.001, 0.0011, 0.0013, 0.0015, 0.0017, 0.002]
    cmap = ListedColormap(["#808080", "#7f00ff", "#0000cd", "#6495ed", "#008000", "#32cd32",
                           "#ffff00", "#ffa500", "#ff8c00", "#ff4500", "#ff0000"])
    cmap.set_under("white")
    ax = axs[0]
    im = ax.pcolormesh(pos, vs, pvd, cmap=cmap, norm=BoundaryNorm(levels, cmap.N), shading="nearest")
    r = np.linspace(1e-3, 0.3, 300)
    pcas = mo.pc_per_arcsec(t)
    vk = np.sqrt(mo.G_PC * t.M_bh / (r * pcas)) * np.sin(np.radians(t.inc))
    ax.plot(0.15 + r, -vk, "r", 0.15 - r, vk, "r")
    ax.set_xlim(0.3, 0)
    ax.set_ylim(-300, 300)
    ax.set_xlabel("major axis position [arcsec]")
    ax.set_ylabel("velocity [km/s]")
    fig.colorbar(im, ax=ax, label="Jy/beam")
    ax.set_title(f"synthetic PV, PA {t.pa - 360:.0f} deg; red: Kepler, 3.4e7 Msun, i = 50", fontsize=9)
    panel_ref(axs[1], "combes25_fig12_pv", "Combes et al. (2025) Fig. 12")
    fig.suptitle(f"NGC 613: {tag}", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / f"{stem}_pv.png", dpi=110)
    plt.close(fig)

    # integrated spectrum within the 0.6 arcsec field [Jy]
    sb = 2 * np.pi * (t.beam_arcsec / mo.FWHM) ** 2 / t.pix_arcsec ** 2
    spec = cb["cube"].sum((1, 2)) / sb
    half = spec >= 0.5 * spec.max()
    fwhm = vs[half].max() - vs[half].min() + t.dv
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    axs[0].step(vs, spec, where="mid", color="k")
    axs[0].fill_between(vs, spec, step="mid", color="cyan")
    axs[0].set_xlim(-500, 500)
    axs[0].set_xlabel("velocity [km/s]")
    axs[0].set_ylabel("flux [Jy]")
    axs[0].set_title(f"synthetic, FWHM {fwhm:.0f} km/s (observed 144), peak {spec.max() * 1e3:.0f} mJy (observed 100)",
                     fontsize=9)
    panel_ref(axs[1], "combes25_fig13_spec", "Combes et al. (2025) Fig. 13, CO(3-2)")
    fig.suptitle(f"NGC 613: {tag}", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / f"{stem}_spectrum.png", dpi=110)
    plt.close(fig)
    print(out / f"{stem}_moments.png")
    print(f"hole {h}, FWHM {fwhm:.0f} km/s, peak {spec.max():.3f} Jy")


def _xy(h, t):
    r = h["offset_pc"] / mo.pc_per_arcsec(t)
    th = np.radians(h["pa"])
    return [r * np.sin(th)], [r * np.cos(th)]


if __name__ == "__main__":
    main()
