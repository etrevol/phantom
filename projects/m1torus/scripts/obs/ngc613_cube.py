#!/usr/bin/env python3
"""NGC 613 CO(3-2) from the ALMA archive cube: how hot is the torus, by three recipes.

The cube: project 2019.1.00273.S, pipeline spw25, SODA cutout of 0.35 arcsec radius
(materials/data/alma_ngc613/README.txt). Moments 0-2 by
  clip    the recipe Combes et al. (2025) state for Fig. 6: channels above 0.6 mJy/beam
          within +-300 km/s of the systemic velocity (1481 km/s)
  mask    smooth-and-mask (Dame 2011): the cube smoothed by one beam and 1.5 channels,
          above 3 times its rms, grown by one beam; the same as our synthetic maps
  gauss   one Gaussian per pixel with peak S/N > 5: sigma with the channel width
          (variance dv^2/12) removed
and the annulus numbers used for the models: |moment 1| (90th percentile) and moment 2
(median) in 0.02-0.05 and 0.05-0.1 arcsec. Figure in Combes' palettes.

    ngc613_cube.py [--out results/m1torus/B/ngc613_cube]
"""
import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.ndimage import gaussian_filter, binary_dilation  # noqa: E402
from scipy.optimize import curve_fit  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "vis"))
from fits_min import read_fits  # noqa: E402
import labels as L  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402

CUBE = ROOT / "materials" / "data" / "alma_ngc613" / "ngc613_co32_spw25_cut035.fits"
C_KMS = 299792.458
V_SYS = 1481.0          # Combes et al. (2025), Fig. 6 caption
FWHM = 2.0 * np.sqrt(2.0 * np.log(2.0))


def load():
    h, d = read_fits(CUBE)
    c = d[0]                                           # (nchan, ny, nx), Jy/beam
    n = h["NAXIS3"]
    f = h["CRVAL3"] + (np.arange(n) + 1 - h["CRPIX3"]) * h["CDELT3"]
    v = C_KMS * (h["RESTFRQ"] / f - 1.0) - V_SYS        # optical convention, about systemic
    order = np.argsort(v)
    c, v = c[order], v[order]
    pix = abs(h["CDELT2"]) * 3600.0
    nx = h["NAXIS1"]
    east = -(np.arange(nx) + 1 - h["CRPIX1"]) * pix     # RA offset, east positive
    north = (np.arange(nx) + 1 - h["CRPIX2"]) * pix
    beam = np.sqrt(h["BMAJ"] * h["BMIN"]) * 3600.0
    return c, v, east, north, pix, beam


def moments(w, v, dv):
    m0 = w.sum(0) * dv
    with np.errstate(invalid="ignore", divide="ignore"):
        m1 = (w * v[:, None, None]).sum(0) * dv / m0
        m2 = np.sqrt((w * (v[:, None, None] - m1) ** 2).sum(0) * dv / m0)
    m1[m0 <= 0] = np.nan
    m2[m0 <= 0] = np.nan
    return m0, m1, m2


def gauss_fit(c, v, rms, dv):
    g = lambda x, a, x0, s: a * np.exp(-0.5 * ((x - x0) / s) ** 2)  # noqa: E731
    ny, nx = c.shape[1:]
    amp, v0, sig = (np.full((ny, nx), np.nan) for _ in range(3))
    sel = np.abs(v) < 400
    for j in range(ny):
        for i in range(nx):
            y = c[sel, j, i]
            k = np.argmax(y)
            if y[k] < 5 * rms:
                continue
            try:
                p, _ = curve_fit(g, v[sel], y, p0=(y[k], v[sel][k], 30.0), maxfev=2000)
            except RuntimeError:
                continue
            s2 = p[2] ** 2 - dv ** 2 / 12.0
            if p[0] > 0 and 0 < s2 and abs(p[1]) < 400:
                amp[j, i], v0[j, i], sig[j, i] = p[0], p[1], np.sqrt(s2)
    return amp, v0, sig


def annuli(m1, m2, r, mask):
    out = []
    for lo, hi in ((0.02, 0.05), (0.05, 0.1)):
        s = (r > lo) & (r < hi) & mask & np.isfinite(m1) & np.isfinite(m2)
        out.append(dict(r=f"{lo}-{hi}", n=int(s.sum()),
                        m1_p90=float(np.percentile(np.abs(m1[s]), 90)) if s.any() else np.nan,
                        m2_med=float(np.median(m2[s])) if s.any() else np.nan))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=str(ROOT / "results" / "m1torus" / "B" / "ngc613_cube"))
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    c, v, east, north, pix, beam = load()
    dv = float(np.median(np.diff(v)))
    rms = float(np.std(c[np.abs(v) > 450]))
    E, N = np.meshgrid(east, north)
    r = np.hypot(E, N)
    win = np.abs(v) < 300

    w_clip = np.where((c > 6e-4) & win[:, None, None], c, 0.0)
    clip = moments(w_clip, v, dv)
    sb = beam / FWHM / pix
    sm = gaussian_filter(c, (1.5, sb, sb))
    rms_s = float(np.std(sm[np.abs(v) > 450]))
    msk = binary_dilation(sm > 3 * rms_s, iterations=max(1, int(round(sb)))) & win[:, None, None]
    w_mask = np.where(msk & (c > 6e-4), c, 0.0)
    mask = moments(w_mask, v, dv)
    amp, v0, sig = gauss_fit(c, v, rms, dv)

    res = dict(rms_mjy=rms * 1e3, dv=dv, beam_arcsec=beam,
               clip=annuli(clip[1], clip[2], r, clip[0] > 0),
               mask=annuli(mask[1], mask[2], r, mask[0] > 0),
               gauss=annuli(v0, sig, r, np.isfinite(sig)))
    (out / "ngc613_cube_numbers.json").write_text(json.dumps(res, indent=1))
    # the masked moments and the cube, for the comparison figures (compare_obs.py)
    np.savez_compressed(out / "ngc613_mask_moments.npz", m0=mask[0], m1=mask[1], m2=mask[2],
                        east=east, north=north, beam=beam)
    np.savez_compressed(out / "ngc613_cube_masked.npz", cube=np.where(msk, c, 0.0).astype(np.float32),
                        v=v, east=east, north=north, beam=beam, rms=rms)

    pal = np.load(HERE.parent / "vis" / "combes25_colours.npz")
    cm = {}
    for k in ("m0", "m1", "m2"):
        cmap = ListedColormap(pal[k])
        cmap.set_bad("white")
        cm[k] = (cmap, float(pal["v" + k[1]][0]), float(pal["v" + k[1]][-1]))
    ext = [east[0], east[-1], north[0], north[-1]]
    fig, axs = plt.subplots(3, 3, figsize=(13, 12.5), layout="constrained")
    rows = [("Stated recipe: clip at 0.6 mJy/beam within $\\pm300$ km/s", clip),
            ("Smooth-and-mask (Dame 2011)", mask),
            ("One Gaussian per pixel, peak S/N > 5", (np.where(np.isfinite(amp), amp * dv * np.sqrt(2 * np.pi) * sig / dv, np.nan), v0, sig))]
    names = ("moment 0 [Jy/beam km/s]", "moment 1 [km/s]", "moment 2 [km/s]")
    gnames = ("integrated flux [Jy/beam km/s]", "centroid [km/s]", r"$\sigma$ [km/s]")
    for i, (title, maps) in enumerate(rows):
        for k, key in enumerate(("m0", "m1", "m2")):
            ax = axs[i, k]
            cmap, lo, hi = cm[key]
            img = maps[k]
            im = ax.imshow(np.where(maps[0] > 0, img, np.nan) if key != "m0" else np.where(img > 0, img, np.nan),
                           origin="lower", extent=ext, cmap=cmap, vmin=lo, vmax=hi, interpolation="nearest")
            ax.set_xlim(0.3, -0.3)
            ax.set_ylim(-0.3, 0.3)
            ax.set_xlabel("RA offset [arcsec]")
            ax.set_ylabel("Dec offset [arcsec]")
            ax.set_title((gnames if i == 2 else names)[k], fontsize=9)
            if k == 2:
                fig.colorbar(im, ax=axs[i, :], fraction=0.02, pad=0.01,
                             label=f"scale of {L.REF['combes25']}, Fig. 6")
        axs[i, 0].text(-0.02, 1.13, title, transform=axs[i, 0].transAxes, fontsize=10.5,
                       fontweight="bold")
    fig.suptitle("NGC 613, CO(3-2), ALMA 2019.1.00273.S archive cube (pipeline, 6.8 km/s channels)\n"
                 f"rms {rms * 1e3:.2f} mJy/beam, beam {beam:.4f} arcsec", fontsize=11)
    fig.savefig(out / "ngc613_cube_moments.png", dpi=110)
    print(json.dumps(res, indent=1))
    print(out / "ngc613_cube_moments.png")


if __name__ == "__main__":
    main()
