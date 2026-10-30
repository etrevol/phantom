#!/usr/bin/env python3
"""B31: the lopsidedness of NGC 613 against tori with modes of different strength.

Metrics within 0.1 arcsec about the black hole (the continuum AGN in the data): the m = 1
harmonic A1 of moment 0 (0.02-0.1 arcsec), the moment-0 centroid offset and its PA.
Model spread: 5 snapshots near the end of each run x 3 azimuths (fitted +-10 deg,
relative to the mode phase) x own noise seeds; data error: 40 noise realisations at the
CO rms. chi2 of the three numbers with sigma^2 = model spread^2 + data error^2, added to
the kinematic chi2 of the snapshot fit (profiles + PV, fit_torus.py).

    b31_lopsided.py
"""
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "vis"))
import fit_torus as T  # noqa: E402
import mockobs as mo  # noqa: E402

EA, NA = 0.0100, 0.0025          # AGN: continuum peak (B18)
R_OUT, R_IN = 0.1, 0.02
FIT = ROOT / "results" / "m1torus" / "B" / "fit_torus"
RUNS = {"no mode (A6)": ("A/runs/20261028-144810-weak45-16k", range(392, 401, 2), "ngc613_B30t96_pv", 0.0),
        "weak (B25)": ("B/runs/20261030-103953-ngc613-best-16k", range(217, 226, 2), "ngc613_B25t225_pv", 0.067),
        "intermediate (B30)": ("B/runs/20261030-164437-ngc613-mid-16k", range(88, 97, 2), "ngc613_B30t96_pv", 0.17),
        "strong (A3)": ("A/runs/20261027-141227-canon-16k-eps0", range(92, 101, 2), "ngc613_A3t100_pv", 0.35)}


def lopsided(m0, E, N, e0, n0, pcas):
    m0 = np.where(np.isfinite(m0), m0, 0.0)
    r = np.hypot(E - e0, N - n0)
    w = np.where(r < R_OUT, m0, 0.0)
    de, dn = (w * E).sum() / w.sum() - e0, (w * N).sum() / w.sum() - n0
    s = (r > R_IN) & (r < R_OUT)
    a1 = abs((m0[s] * np.exp(1j * np.arctan2(N - n0, E - e0)[s])).sum()) / m0[s].sum()
    return np.array([a1, np.hypot(de, dn) * pcas, np.degrees(np.arctan2(de, dn)) % 360])


def data_metrics():
    c = np.load(ROOT / "results/m1torus/B/ngc613_cube/ngc613_cube_win.npz")
    d = np.load(ROOT / "results/m1torus/B/ngc613_cube/ngc613_mask_moments.npz")
    E, N = np.meshgrid(d["east"], d["north"])
    t = T.target()
    pcas = mo.pc_per_arcsec(t)
    base = lopsided(d["m0"], E, N, EA, NA, pcas)
    cube = c["cube"].astype(float)
    pix = t.pix_arcsec
    sb = t.beam_arcsec / mo.FWHM / pix
    rms = float(c["rms"])
    rms_s = float(np.std(gaussian_filter(cube, (1.5, sb, sb), mode="wrap")[:, :8, :]))
    rng = np.random.default_rng(11)
    reals = []
    for _ in range(40):
        nz = gaussian_filter(rng.standard_normal(cube.shape), (0, sb, sb), mode="wrap")
        nz *= rms / nz.std()
        m0 = mo.moments(dict(cube=cube + nz, vs=c["v"], xs=d["north"]), t, rms=rms_s)[0]
        reals.append(lopsided(m0, E, N, EA, NA, pcas))
    reals = np.array(reals)
    err = np.r_[reals[:, :2].std(0), np.degrees(np.std(np.unwrap(np.radians(reals[:, 2]))))]
    return base, err


def model_metrics(run, times, fitdir):
    b = json.load(open(FIT / fitdir / "fit.json"))["best"]
    t = replace(T.target(), R_pc=b["R_pc"], inc=b["inc"], M_bh=10 ** b["logM"])
    pcas = mo.pc_per_arcsec(t)
    out = []
    for k, ti in enumerate(times):
        p = mo.load(ROOT / "results/m1torus" / run / f"m1torus_{ti:05d}")
        ph = mo.m1_phase(p["x"], p["m"])
        for j, daz in enumerate((-10, 0, 10)):
            cb = mo.cube(p, t, ph + b["az"] + daz, seed=100 + 3 * k + j)
            E, N = np.meshgrid(cb["xs"], cb["xs"])
            out.append(lopsided(mo.moments(cb, t)[0], E, N, 0.0, 0.0, pcas))
    return np.array(out)


def main():
    T.GAL = "ngc613"
    T.EDGES = (0.01, 0.025, 0.04, 0.06, 0.08, 0.1)
    data, derr = data_metrics()
    print(f"NGC 613: A1 {data[0]:.3f} +- {derr[0]:.3f}, offset {data[1]:.2f} +- {derr[1]:.2f} pc, "
          f"PA {data[2]:.0f} +- {derr[2]:.0f} deg (noise only)")
    res = {}
    for name, (run, times, fitdir, a1s0) in RUNS.items():
        m = model_metrics(run, times, fitdir)
        mean = np.r_[m[:, :2].mean(0), np.degrees(np.angle(np.exp(1j * np.radians(m[:, 2])).mean())) % 360]
        spread = np.r_[m[:, :2].std(0), np.degrees(np.std(np.unwrap(np.radians(m[:, 2]))))]
        sig = np.hypot(spread, derr)
        dpa = (mean[2] - data[2] + 180) % 360 - 180
        c = ((mean[0] - data[0]) / sig[0]) ** 2 + ((mean[1] - data[1]) / sig[1]) ** 2 + (dpa / sig[2]) ** 2
        kin = json.load(open(FIT / fitdir / "fit.json"))["chi2"] if "A6" not in name else np.nan
        res[name] = dict(A1_S0_saturated=a1s0, mean=list(map(float, mean)), spread=list(map(float, spread)),
                         chi2_m1=float(c), chi2_kin=float(kin), chi2_total=float(kin + c))
        print(f"{name:20s} A1 {mean[0]:.3f} +- {spread[0]:.3f}, offset {mean[1]:.2f} +- {spread[1]:.2f} pc, "
              f"PA {mean[2]:.0f} +- {spread[2]:.0f}; chi2 m=1 {c:.1f}, kinematic {kin:.1f}, total {kin + c:.1f}")
    out = ROOT / "results/m1torus/B/b31_lopsided.json"
    out.write_text(json.dumps(dict(data=list(map(float, data)), data_err=list(map(float, derr)), models=res), indent=1))
    print(out)


if __name__ == "__main__":
    main()
