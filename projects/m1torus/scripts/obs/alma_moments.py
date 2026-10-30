#!/usr/bin/env python3
"""Moments 0-2 of an archive cube by the recipe of the synthetic maps (mockobs.moments).

Writes <out>/<galaxy>_mask_moments.npz (m0, m1, m2, east, north, beam) and
<out>/<galaxy>_cube_win.npz (the cube within the velocity window, v, east, north, beam,
rms), the inputs of fit_torus.py, and a figure of the three maps.

    alma_moments.py ngc1672
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
from scipy.ndimage import gaussian_filter  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "vis"))
from fits_min import read_fits  # noqa: E402
from galaxies import GALAXIES, C_KMS  # noqa: E402
import mockobs as mo  # noqa: E402


def load(g):
    h, d = read_fits(g["cube"])
    c = d[0]
    n = h["NAXIS3"]
    f = h["CRVAL3"] + (np.arange(n) + 1 - h["CRPIX3"]) * h["CDELT3"]
    v = C_KMS * (h["RESTFRQ"] / f - 1.0) - g["v_sys"]    # optical convention, about systemic
    o = np.argsort(v)
    pix = abs(h["CDELT2"]) * 3600.0
    nx = h["NAXIS1"]
    east = -(np.arange(nx) + 1 - h["CRPIX1"]) * pix
    north = (np.arange(nx) + 1 - h["CRPIX2"]) * pix
    beam = np.sqrt(h["BMAJ"] * h["BMIN"]) * 3600.0
    return c[o], v[o], east, north, pix, beam


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("galaxy", choices=sorted(GALAXIES))
    a = ap.parse_args()
    g = GALAXIES[a.galaxy]
    out = Path(g["out"])
    out.mkdir(parents=True, exist_ok=True)
    c, v, east, north, pix, beam = load(g)
    dv = float(np.median(np.diff(v)))
    free = np.abs(v) > g["window"] + 50.
    rms = float(np.std(c[free]))
    sb = beam / mo.FWHM / pix
    rms_s = float(np.std(gaussian_filter(c, (1.5, sb, sb), mode="wrap")[free]))
    win = np.abs(v) < g["window"]
    tgt = replace(mo.NGC613, beam_arcsec=beam, dv=dv, pix_arcsec=pix, clip_mjy=g["clip_mjy"])
    m0, m1, m2 = mo.moments(dict(cube=c[win], vs=v[win], xs=north), tgt, rms=rms_s)
    np.savez_compressed(out / f"{a.galaxy}_mask_moments.npz", m0=m0, m1=m1, m2=m2, east=east, north=north, beam=beam)
    np.savez_compressed(out / f"{a.galaxy}_cube_win.npz", cube=c[win].astype(np.float32), v=v[win], east=east,
                        north=north, beam=beam, rms=rms)
    E, N = np.meshgrid(east, north)
    r = np.hypot(E, N)
    bpix = np.pi * beam ** 2 / (4 * np.log(2)) / pix ** 2
    info = dict(galaxy=g["name"], rms_mjy=rms * 1e3, dv=dv, beam=beam, pix=pix, window=g["window"],
                flux_jy_kms={f"r<{R}": float(np.nansum(np.where(r < R, m0, 0)) / bpix) for R in (0.1, 0.2, 0.3, 0.4)})
    (out / f"{a.galaxy}_moments.json").write_text(json.dumps(info, indent=1))
    from compare_obs import PAL
    fig, axs = plt.subplots(1, 3, figsize=(15, 5), layout="constrained")
    lim = g["field"]
    for ax, key, img in zip(axs, ("m0", "m1", "m2"), (m0, m1, m2)):
        cmap, lo, hi = PAL[key]
        if key == "m1":
            lo, hi = -g["window"] * 0.75, g["window"] * 0.75
        im = ax.imshow(np.where(m0 > 0, img, np.nan), origin="lower", cmap=cmap, vmin=lo, vmax=hi,
                       extent=[east[0], east[-1], north[0], north[-1]], interpolation="nearest")
        ax.set(xlim=(lim, -lim), ylim=(-lim, lim), xlabel="RA offset [arcsec]", ylabel="Dec offset [arcsec]")
        fig.colorbar(im, ax=ax, shrink=0.8, label=("moment 0 [Jy/beam km/s]", "moment 1 [km/s]",
                                                     "moment 2 [km/s]")[["m0", "m1", "m2"].index(key)])
    fig.suptitle(f"{g['name']}, CO(3-2), ALMA 2019.1.00273.S archive cube, smooth-and-mask moments "
                 f"(beam {beam:.4f} arcsec, rms {rms * 1e3:.2f} mJy/beam)", fontsize=10)
    fig.savefig(out / f"{a.galaxy}_moments.png", dpi=110)
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
