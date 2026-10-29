#!/usr/bin/env python3
"""Synthetic NGC 613 maps of a dump next to the published ALMA panels (Combes et al. 2025).

Ours on the left, theirs on the right, on the same field, colour scale and clip:
  <stem>_moments.png   CO(3-2) moments 0, 1, 2 against their Fig. 6
  <stem>_pv.png        PV along PA -20 deg against their Fig. 12, Keplerian curve of 3.4e7 Msun
  <stem>_spectrum.png  integrated spectrum against their Fig. 13 (FWHM 144 km/s, Table 2)
The reference crops come from ref_figures.py (results/m1torus/B/reference).

    compare_obs.py DUMP [--phi DEG] [--out DIR]

The model is named from its run directory (labels.py, FIGURES.md).
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm, ListedColormap  # noqa: E402

import mockobs as mo  # noqa: E402
import labels as L  # noqa: E402

ROOT = Path(__file__).resolve().parents[4]
REF = ROOT / "results" / "m1torus" / "B" / "reference"
# NGC 613 from the ALMA archive cube, processed like the synthetic maps (scripts/obs/ngc613_cube.py)
CUBE_DIR = ROOT / "results" / "m1torus" / "B" / "ngc613_cube"
ALMA = "ALMA 2019.1.00273.S archive cube, same processing"
# Palettes and value ranges read pixel by pixel from the colour bars of Combes et al.
# (2025) Figs. 6 and 12 (400 dpi render; ranges from the tick marks): moment 0
# 0.052-0.280 Jy/beam km/s, moment 1 -108..+109 km/s, moment 2 22-194 km/s, PV in ten
# discrete levels 0.0002-0.002 Jy/beam. Ours and theirs then share one colour scale.
_C = np.load(Path(__file__).with_name("combes25_colours.npz"))
PAL = {}
for k in ("m0", "m1", "m2"):
    cm = ListedColormap(_C[k], name=f"combes25_{k}")
    cm.set_bad("white")
    PAL[k] = (cm, float(_C["v" + k[1]][0]), float(_C["v" + k[1]][-1]))
PV_LEVELS = list(_C["pv_levels"])
PV_CMAP = ListedColormap(_C["pv_rgb"], name="combes25_pv")
PV_CMAP.set_under("white")
PV_CMAP.set_over(_C["pv_rgb"][-1])


def panel_ref(ax, name, title, extent=None):
    """A published panel; with extent [left, right, bottom, top] its frame is drawn on
    calibrated axes (Fig. 6: +-0.30 arcsec from the tick marks, see ref_figures.py)."""
    f = REF / f"{name}.png"
    if not f.exists():
        ax.text(0.5, 0.5, f"run ref_figures.py\n({name})", ha="center", transform=ax.transAxes)
        ax.set_axis_off()
    elif extent is None:
        ax.imshow(plt.imread(f))
        ax.set_axis_off()
    else:
        ax.imshow(plt.imread(f), extent=extent, aspect="equal")
        ax.set_xlabel("RA offset [arcsec]")
        ax.set_ylabel("Dec offset [arcsec]")
    ax.set_title(title, fontsize=9)


def beam(ax, t, xs):
    from matplotlib.patches import Circle
    r = 0.5 * t.beam_arcsec
    ax.add_patch(Circle((xs[-1] - 3 * r - 0.02, xs[0] + 3 * r + 0.02), r, color="k"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dump")
    ap.add_argument("--phi", type=float, default=None)
    ap.add_argument("--label", default="", help="ignored; the run directory names the model")
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
    model = L.wrapped(L.model_label(Path(a.dump).resolve().parent), 100)
    tag = (model + "\n" + rf"$t = {p['time'] / (2 * np.pi):.1f}\,P_{{\rm orb}}$, "
           + rf"viewing azimuth ${phi:.0f}^\circ$ from the torus $x$ axis")
    st = Path(a.dump).resolve().parent / "status"
    if st.exists() and "FAILED" in st.read_text():
        tag += "\n(the run stopped soon after this state; see the log)"
    c25 = L.REF["combes25"]
    ext = [xs[0], xs[-1], xs[0], xs[-1]]

    # moments wherever the mask kept emission; below the lowest colour it shows as the
    # lowest colour (black), as the faint specks of Combes et al. (2025), Fig. 6
    show = m0 > 0
    real = np.load(CUBE_DIR / "ngc613_mask_moments.npz") if (CUBE_DIR / "ngc613_mask_moments.npz").exists() else None
    fig, axs = plt.subplots(3, 3, figsize=(14.5, 13), layout="constrained")
    specs = [(np.where(show, m0, np.nan), "m0", "moment 0 [Jy/beam km/s]", "combes25_fig6_m0"),
             (np.where(show, m1, np.nan), "m1", "moment 1 [km/s]", "combes25_fig6_m1"),
             (np.where(show, m2, np.nan), "m2", "moment 2 [km/s]", "combes25_fig6_m2")]
    for row, (img, key, lab, ref) in enumerate(specs):
        ax = axs[row, 0]
        cmap, lo, hi = PAL[key]
        im = ax.imshow(img, origin="lower", extent=ext, cmap=cmap, vmin=lo, vmax=hi,
                       interpolation="nearest")
        ax.set_xlim(0.3, -0.3)
        ax.set_ylim(-0.3, 0.3)
        # above the published range the top colour is near white like "no data": outline it
        over = np.nan_to_num(img, nan=-np.inf) > hi
        if over.any():
            ax.contour(xs, xs, over.astype(float), levels=[0.5], colors="k", linewidths=0.7)
            ax.text(0.02, 0.97, f"above the scale inside the black line (max {np.nanmax(img):.3g})",
                    transform=ax.transAxes, va="top", fontsize=7.5)
        ax.plot(0, 0, "w+", ms=8, mew=1.5)
        beam(ax, t, xs)
        ax.set_xlabel("RA offset [arcsec]")
        ax.set_ylabel("Dec offset [arcsec]")
        ax.set_title(f"Synthetic CO(3-2), {lab.split(' [')[0]}", fontsize=9)
        # the archive cube with the same mask, beam and palette
        ax = axs[row, 1]
        if real is not None:
            rimg = real[key]
            rimg = np.where(real["m0"] > 0, rimg, np.nan)
            ax.imshow(rimg, origin="lower", extent=[real["east"][0], real["east"][-1],
                      real["north"][0], real["north"][-1]], cmap=cmap, vmin=lo, vmax=hi,
                      interpolation="nearest")
            ax.set_xlim(0.3, -0.3)
            ax.set_ylim(-0.3, 0.3)
            ax.plot(0, 0, "w+", ms=8, mew=1.5)
        ax.set_xlabel("RA offset [arcsec]")
        ax.set_ylabel("Dec offset [arcsec]")
        ax.set_title(f"NGC 613, {lab.split(' [')[0]}: {ALMA}", fontsize=8)
        panel_ref(axs[row, 2], ref, f"{c25}, Fig. 6, {lab.split(' [')[0]}",
                  extent=[0.3, -0.3, -0.3, 0.3])
        # one scale for the row: their palette and range, read from their colour bar
        fig.colorbar(im, ax=list(axs[row]), fraction=0.02, pad=0.01,
                     label=lab + f"\n(all panels; scale of {c25}, Fig. 6)")
    axs[0, 0].plot(*_xy(h, t), "wx", ms=9, mew=2)
    fig.suptitle(f"NGC 613, CO(3-2) at 0.015 arcsec\n{tag}\n"
                 rf"moment-0 hole {h['offset_pc']:.1f} pc from the AGN at PA ${h['pa']:.0f}^\circ$ "
                 f"(observed: 3 pc to the NW)", fontsize=10)
    fig.savefig(out / f"{stem}_moments.png", dpi=110)
    plt.close(fig)

    # PV along PA -20 deg (their axis: 0.3 .. 0 arcsec, AGN at 0.15)
    pvd = mo.pv(cb, t, pa=t.pa, halfwidth=t.beam_arcsec)
    pos = 0.15 - xs
    fig, axs = plt.subplots(1, 3, figsize=(18, 5.2), layout="constrained")
    ax = axs[0]
    im = ax.pcolormesh(pos, vs, pvd, cmap=PV_CMAP, norm=BoundaryNorm(PV_LEVELS, PV_CMAP.N),
                       shading="nearest")
    r = np.linspace(1e-3, 0.3, 300)
    pcas = mo.pc_per_arcsec(t)
    vk = np.sqrt(mo.G_PC * t.M_bh / (r * pcas)) * np.sin(np.radians(t.inc))
    ax.plot(0.15 + r, -vk, "r", 0.15 - r, vk, "r")
    ax.set_xlim(0.3, 0)
    ax.set_ylim(-300, 300)
    ax.set_xlabel("major axis position [arcsec]")
    ax.set_ylabel("velocity [km/s]")
    fig.colorbar(im, ax=list(axs), fraction=0.02, pad=0.01,
                 label=f"Jy/beam (all panels; levels of {c25}, Fig. 12)")
    cfile = CUBE_DIR / "ngc613_cube_win.npz"
    if cfile.exists():
        rc = np.load(cfile)
        e = rc["east"][::-1]                       # ascending east for the slit sampler
        cbr = dict(cube=rc["cube"][:, :, ::-1].astype(float), xs=e, vs=rc["v"])
        pvr = mo.pv(cbr, t, pa=t.pa, halfwidth=float(rc["beam"]))
        axs[1].pcolormesh(0.15 - e, rc["v"], pvr, cmap=PV_CMAP,
                          norm=BoundaryNorm(PV_LEVELS, PV_CMAP.N), shading="nearest")
        axs[1].plot(0.15 + r, -vk, "r", 0.15 - r, vk, "r")
        axs[1].set(xlim=(0.3, 0), ylim=(-300, 300), xlabel="major axis position [arcsec]",
                   ylabel="velocity [km/s]")
        axs[1].set_title(f"NGC 613 PV, same slit: {ALMA}", fontsize=8)
    ax.set_title(rf"Synthetic CO(3-2) PV along PA ${t.pa - 360:.0f}^\circ$; red: Keplerian, "
                 r"$3.4\times10^7\,M_\odot$, $i = 50^\circ$", fontsize=9)
    panel_ref(axs[2], "combes25_fig12_pv", f"{c25}, Fig. 12")
    fig.suptitle(f"NGC 613\n{tag}", fontsize=10)
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
    axs[0].set_title(f"Synthetic CO(3-2) spectrum: FWHM {fwhm:.0f} km/s (observed 144), "
                     f"peak {spec.max() * 1e3:.0f} mJy (observed 100)", fontsize=9)
    panel_ref(axs[1], "combes25_fig13_spec", f"{c25}, Fig. 13")
    fig.suptitle(f"NGC 613\n{tag}", fontsize=10)
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
