#!/usr/bin/env python3
"""B14: how much of the gap between the models and NGC 613 is viewing geometry.

Each snapshot is observed on a grid of inclination i and black-hole mass M_BH (the
velocity unit, at R_tor = 8 pc) from four azimuths about its m = 1 phase, with the beam,
channels and noise of the ALMA archive cube (2019.1.00273.S). Moments come from
mockobs.moments, the function the archive cube goes through (ngc613_cube.py), so the
recipe is the same on both sides. Metrics against the cube:
  |moment 1| (90th percentile) and moment 2 (median) in 0.02-0.05 and 0.05-0.1 arcsec,
  the axis ratio b/a of moment 0 (intensity-weighted second moments within 0.15 arcsec).
chi2 = sum of squared log ratios over sigma_log = 0.15 (the velocities move by ~0.1-0.25
between mask thresholds of 3 and 4.75 sigma, LOG 2026-10-29) and of (b/a) differences
over 0.1; the minimum over azimuths at each (i, M_BH).

    fit_geometry.py [--workers 4] [--out results/m1torus/B/fit_geometry]
"""
import argparse
import json
import sys
from dataclasses import replace
from multiprocessing import Pool
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE.parent / "vis"))
import labels as L  # noqa: E402
import mockobs as mo  # noqa: E402

A = ROOT / "results" / "m1torus" / "A" / "runs"
B = ROOT / "results" / "m1torus" / "B" / "runs"
SNAPS = {"A3": A / "20261027-141227-canon-16k-eps0" / "m1torus_00100",
         "B6b": B / "20261029-152300-clouds-from-a3t100-e05",
         "B8a": B / "20261029-163625-gas-from-a3t100-hor015",
         "B12": B / "20261029-181639-heated-s01-from-a3t100"}
INC = np.arange(15, 66, 5)
MBH = np.arange(1.5, 4.51, 0.25) * 1e7
AZ = (0, 90, 180, 270)               # azimuth about the m = 1 phase [deg]
ANN = ((0.02, 0.05), (0.05, 0.1))
SIG_LOG, SIG_Q = 0.15, 0.1
CUBE_DIR = ROOT / "results" / "m1torus" / "B" / "ngc613_cube"


def last_dump(p):
    p = Path(p)
    return p if p.is_file() else sorted(p.glob("m1torus_[0-9][0-9][0-9][0-9][0-9]"))[-1]


def target():
    """NGC 613 at the resolution of the archive cube."""
    c = np.load(CUBE_DIR / "ngc613_cube_win.npz")
    return replace(mo.NGC613, beam_arcsec=float(c["beam"]), dv=float(np.median(np.diff(c["v"]))),
                   pix_arcsec=float(abs(c["north"][1] - c["north"][0])), rms_mjy=float(c["rms"]) * 1e3)


def metrics(m0, m1, m2, east, north):
    E, N = np.meshgrid(east, north)
    r = np.hypot(E, N)
    out = []
    for lo, hi in ANN:
        s = (r > lo) & (r < hi) & (m0 > 0) & np.isfinite(m1) & np.isfinite(m2)
        out += [np.percentile(np.abs(m1[s]), 90) if s.sum() > 20 else np.nan,
                np.median(m2[s]) if s.sum() > 20 else np.nan]
    w = np.where((r < 0.15) & (m0 > 0), m0, 0.0)
    if w.sum() <= 0:
        return out + [np.nan]
    xe, yn = (w * E).sum() / w.sum(), (w * N).sum() / w.sum()
    cxx = (w * (E - xe) ** 2).sum()
    cyy = (w * (N - yn) ** 2).sum()
    cxy = (w * (E - xe) * (N - yn)).sum()
    lam = np.linalg.eigvalsh([[cxx, cxy], [cxy, cyy]])
    return out + [np.sqrt(lam[0] / lam[1])]


def observed():
    d = np.load(CUBE_DIR / "ngc613_mask_moments.npz")
    return np.array(metrics(d["m0"], d["m1"], d["m2"], d["east"], d["north"]))


def chi2(mod, obs):
    v = np.log(mod[:4] / obs[:4]) / SIG_LOG
    q = (mod[4] - obs[4]) / SIG_Q
    return float(np.nansum(v ** 2) + q ** 2) if np.all(np.isfinite(mod)) else np.inf


def _grid(inc, mbh):
    global INC, MBH
    INC, MBH = inc, mbh


def work(job):
    name, az = job
    t0 = target()
    p = mo.load(last_dump(SNAPS[name]))
    phi = mo.m1_phase(p["x"], p["m"]) + az
    res = np.full((len(INC), len(MBH), 5), np.nan)
    for a, inc in enumerate(INC):
        for b, mbh in enumerate(MBH):
            t = replace(t0, inc=float(inc), M_bh=float(mbh))
            cb = mo.cube(p, t, phi)
            m0, m1, m2 = mo.moments(cb, t)
            # cube columns run along east, rows along north, both on xs
            res[a, b] = metrics(m0, m1, m2, cb["xs"], cb["xs"])
    return name, az, res


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out", default=str(ROOT / "results" / "m1torus" / "B" / "fit_geometry"))
    ap.add_argument("--inc", default="15,65,5", help="lo,hi,step [deg]")
    ap.add_argument("--mbh", default="1.5,4.5,0.25", help="lo,hi,step [1e7 Msun]")
    ap.add_argument("--models", default=",".join(SNAPS))
    a = ap.parse_args()
    global INC, MBH
    lo, hi, st = map(float, a.inc.split(","))
    INC = np.arange(lo, hi + 1e-6, st).astype(int)
    lo, hi, st = map(float, a.mbh.split(","))
    MBH = np.arange(lo, hi + 1e-6, st) * 1e7
    for n in set(SNAPS) - set(a.models.split(",")):
        del SNAPS[n]
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    obs = observed()
    jobs = [(n, az) for n in SNAPS for az in AZ]
    with Pool(a.workers, initializer=_grid, initargs=(INC, MBH)) as pool:
        got = pool.map(work, jobs)
    grid = {n: np.full((len(AZ), len(INC), len(MBH), 5), np.nan) for n in SNAPS}
    for n, az, res in got:
        grid[n][AZ.index(az)] = res
    np.savez_compressed(out / "fit_geometry.npz", inc=INC, mbh=MBH, az=AZ, obs=obs,
                        **{n: g for n, g in grid.items()})

    summary = dict(observed=dict(zip(["m1_p90_in", "m2_med_in", "m1_p90_out", "m2_med_out", "b_a"],
                                     map(float, obs))))
    axs_n = len(SNAPS)
    fig, axs = plt.subplots(1, axs_n, figsize=(4.2 * axs_n + 1, 4.4), layout="constrained",
                            sharey=True, squeeze=False)
    axs = axs[0]
    for ax, (n, g) in zip(axs, grid.items()):
        c2 = np.array([[[chi2(g[k, i, j], obs) for j in range(len(MBH))] for i in range(len(INC))]
                       for k in range(len(AZ))])
        best_az = np.argmin(c2, axis=0)
        cmin = c2.min(axis=0)
        i, j = np.unravel_index(np.argmin(cmin), cmin.shape)
        k = best_az[i, j]
        i50, j34 = int(np.argmin(np.abs(INC - 50))), int(np.argmin(np.abs(MBH - 3.4e7)))
        summary[n] = dict(best_inc=int(INC[i]), best_mbh=float(MBH[j]), best_az=AZ[k],
                          chi2_best=float(cmin[i, j]), chi2_default=float(cmin[i50, j34]),
                          metrics_best=list(map(float, g[k, i, j])),
                          metrics_default=list(map(float, g[best_az[i50, j34], i50, j34])),
                          ratio_best=[float(g[k, i, j, 1] / g[k, i, j, 0]), float(g[k, i, j, 3] / g[k, i, j, 2])],
                          model=L.model_label(last_dump(SNAPS[n]).parent))
        im = ax.pcolormesh(MBH / 1e7, INC, np.log10(cmin), cmap="viridis_r", shading="nearest",
                           vmin=0, vmax=3)
        ax.contour(MBH / 1e7, INC, cmin - cmin[i, j], levels=[2.3, 6.2], colors="w",
                   linewidths=[1.2, 0.8])
        ax.plot(MBH[j] / 1e7, INC[i], "r*", ms=13)
        ax.plot(3.4, 50, "wo", mfc="none", ms=9, mew=1.5)
        ax.set_xlabel(r"$M_{\rm BH}$ [$10^7\,M_\odot$]")
        name = L.model_parts(last_dump(SNAPS[n]).parent)[0]
        ax.set_title(f"{name}\nbest $i = {INC[i]}^\\circ$, $M_{{\\rm BH}} = {MBH[j] / 1e7:.2f}$, "
                     f"$\\chi^2 = {cmin[i, j]:.1f}$ ($i = 50^\\circ$: {cmin[i50, j34]:.0f})", fontsize=8.5)
    axs[0].set_ylabel(r"inclination $i$ [deg]")
    fig.colorbar(im, ax=list(axs), label=r"$\log_{10}\chi^2$ (best of four viewing azimuths)")
    fig.suptitle("NGC 613: viewing geometry of the models against the ALMA 2019.1.00273.S archive cube\n"
                 r"$|$moment 1$|$ and moment 2 in two annuli, moment-0 axis ratio; "
                 r"star: best fit, circle: $i = 50^\circ$, $M_{\rm BH} = 3.4\times10^7\,M_\odot$ "
                 r"(Combes et al. 2025); white: $\Delta\chi^2 = 2.3$, 6.2", fontsize=10)
    fig.savefig(out / "fit_geometry.png", dpi=130)
    (out / "fit_geometry.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))
    print(out / "fit_geometry.png")


if __name__ == "__main__":
    main()
