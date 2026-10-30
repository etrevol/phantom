#!/usr/bin/env python3
"""B19: fit a collisionless torus (or a run's snapshot) to NGC 613 in the ALMA archive cube.

Model sources
  torus      phase-mixed torus from the Keplerian elements of setup_m1torus.f90 (a uniform
             in [1 - R0, 1 + R0], e uniform in [0, e_max], i uniform in [-i_max, i_max],
             node random, periapsis on the node), mean anomaly uniform: the state a run
             reaches before the m = 1 mode, without self-gravity; same random numbers for
             every parameter set, so the objective is smooth
  torus+disc the torus plus a central disc of clouds: surface density ~ 1/R from
             0.02 R_tor to R_d, circular Keplerian orbits with a dispersion h_d v_K in
             every component, f_d of the mass, tilted by tilt_d about the sky major axis
             (its inclination is inc + tilt_d; the view is fixed at azimuth 0)
  snap:PATH  a dump (or a run directory: its last dump); only the view is fitted
Parameters  R0, e_max, i_max (torus); R_d, f_d, h_d, tilt_d (disc); R_pc (torus radius), inc, logM (log10 M_BH, the
            dynamical mass), az (view azimuth about the m = 1 phase, snapshots only)
Metrics     profiles (default): in six annuli 0.01-0.15 arcsec the share of moment 0,
            |moment 1| (90th percentile) and moment 2 (median), plus the moment-0 axis
            ratio: 19 numbers; basic: the 5 numbers of fit_geometry.py
Search      differential evolution (scipy) over the free parameters, all cores; every
            evaluation is kept, intervals are the ranges of evaluations within
            Delta chi2 = chi2_min / dof of the best (errors rescaled to a fit of
            chi2/dof = 1); a best value at a bound is flagged
Speed       generated tori have one particle size, so the cube is a histogram smoothed
            once by the beam and the line width, the same sum as the per-particle
            stamps of mockobs.cube; the noise is drawn once (same seed as mockobs.cube)

    fit_torus.py [--model torus|snap:PATH] [--free R0,e_max,i_max,R_pc,inc,logM]
                 [--fix name=value ...] [--bounds name=lo:hi ...] [--metrics profiles|basic]
                 [--maxiter 40] [--popsize 12] [--workers 8] [--out DIR] [--check]
"""
import argparse
import json
import sys
import time
from dataclasses import replace
from functools import partial
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.ndimage import gaussian_filter  # noqa: E402
from scipy.optimize import differential_evolution  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "vis"))
import fit_geometry as F  # noqa: E402
import labels as L  # noqa: E402
import mockobs as mo  # noqa: E402

N = 16384
SIZE = 0.05                     # particle size [R_tor] of generated tori (h/1.4 of the 16k runs)
EDGES = (0.01, 0.025, 0.04, 0.06, 0.08, 0.11, 0.15)    # annuli [arcsec]
SIG = dict(m0=0.2, m1=0.15, m2=0.15, ba=0.1)  # log errors (m0, m1, m2) and absolute (b/a)
MISS = 25.                      # chi2 of an annulus the model leaves empty
DEFAULT = dict(R0=0.3, e_max=0.5, i_max=60., R_pc=8., inc=50., logM=np.log10(3.4e7), az=90.,
               R_d=0.3, f_d=0.2, h_d=0.1, tilt_d=0.)
BOUNDS = dict(R0=(0.1, 0.8), e_max=(0.05, 0.8), i_max=(5., 89.), R_pc=(1.5, 10.),
              inc=(5., 80.), logM=(6.9, 8.0), az=(0., 360.),
              R_d=(0.05, 0.8), f_d=(0.0, 0.7), h_d=(0.02, 0.4), tilt_d=(-90., 90.))
CROP = 0.18                     # half-size of the field the metrics need [arcsec]
_noise = {}
_snap = {}


def torus(R0, e_max, i_max, n=N, seed=1):
    """Positions and velocities [code units] of a phase-mixed torus of Keplerian orbits."""
    rng = np.random.default_rng(seed)
    u = rng.random((6, n))
    a = (1. - R0) + 2. * R0 * u[0]
    e = e_max * u[1]
    inc = np.radians((2. * u[2] - 1.) * i_max)
    node = 2. * np.pi * u[3]
    M = 2. * np.pi * u[4]
    E = M.copy()
    for _ in range(40):                 # Kepler's equation by Newton steps
        E -= (E - e * np.sin(E) - M) / (1. - e * np.cos(E))
    nu = 2. * np.arctan2(np.sqrt(1. + e) * np.sin(E / 2.), np.sqrt(1. - e) * np.cos(E / 2.))
    p = a * (1. - e ** 2)
    r = p / (1. + e * np.cos(nu))
    xp, yp = r * np.cos(nu), r * np.sin(nu)
    vxp, vyp = -np.sin(nu) / np.sqrt(p), (e + np.cos(nu)) / np.sqrt(p)
    ci, si, cn, sn = np.cos(inc), np.sin(inc), np.cos(node), np.sin(node)

    def rot(x, y):                      # tilt by i about x, then turn by the node about z
        y, z = y * ci, y * si
        return np.c_[x * cn - y * sn, x * sn + y * cn, z]
    return dict(x=rot(xp, yp), v=rot(vxp, vyp), m=np.full(n, 1. / n), size=np.full(n, SIZE), m_tor=1.)


def disc(R_d, h_d, tilt_d, n, seed=2):
    """A central disc of clouds: Sigma ~ 1/R in [0.02, R_d], circular orbits about a unit
    mass with a dispersion h_d v_K, tilted by tilt_d [deg] about the y axis, which is
    the sky major axis for the view of fast_cube (azimuth 0): the disc's inclination
    becomes inc + tilt_d."""
    rng = np.random.default_rng(seed)
    u = rng.random((2, n))
    g = rng.standard_normal((4, n))
    R = 0.02 + (max(R_d, 0.021) - 0.02) * u[0]
    ph = 2. * np.pi * u[1]
    vk = 1. / np.sqrt(R)
    x = np.c_[R * np.cos(ph), R * np.sin(ph), h_d * R * g[0]]
    vR, vp, vz = h_d * vk * g[1], vk + h_d * vk * g[2], h_d * vk * g[3]
    v = np.c_[vR * np.cos(ph) - vp * np.sin(ph), vR * np.sin(ph) + vp * np.cos(ph), vz]
    c, s = np.cos(np.radians(tilt_d)), np.sin(np.radians(tilt_d))
    rot = np.array([[c, 0., -s], [0., 1., 0.], [s, 0., c]])
    return x @ rot.T, v @ rot.T


def noise(t, shape):
    """The noise of mockobs.cube (seed 1), drawn once per process."""
    if shape not in _noise:
        rng = np.random.default_rng(1)
        sb = t.beam_arcsec / mo.FWHM / t.pix_arcsec
        nz = gaussian_filter(rng.standard_normal(shape), (0, sb, sb), mode="wrap")
        nz = nz / nz.std() * t.rms_mjy * 1e-3
        # rms of the smoothed noise as mockobs.moments measures it (edge strip)
        sm = gaussian_filter(nz, (1.5, sb, sb), mode="wrap")
        _noise[shape] = (nz, float(np.std(sm[:, :8, :])))
    return _noise[shape]


def fast_cube(p, t, phi):
    """mockobs.cube for particles of one size: binned fluxes smoothed by beam and line width."""
    pcas = mo.pc_per_arcsec(t)
    east, north, _, vlos = mo.to_sky(p["x"], p["v"], t.inc, t.pa, phi)
    nx = int(round(2 * t.field_arcsec / t.pix_arcsec)) | 1
    nv = int(round(2 * t.vmax / t.dv)) | 1
    xs = (np.arange(nx) - nx // 2) * t.pix_arcsec
    vs = (np.arange(nv) - nv // 2) * t.dv
    flux = p["m"] / p["m_tor"] * t.m_gas * t.jy_kms_per_msun
    ex = np.r_[xs - t.pix_arcsec / 2, xs[-1] + t.pix_arcsec / 2]
    ev = np.r_[vs - t.dv / 2, vs[-1] + t.dv / 2]
    v = vlos * np.sqrt(mo.G_PC * t.M_bh / t.R_pc)
    h, _ = np.histogramdd(np.c_[v, north * t.R_pc / pcas, east * t.R_pc / pcas], bins=(ev, ex, ex),
                          weights=flux)
    size = float(np.median(p["size"]))
    sig_beam = t.beam_arcsec / mo.FWHM
    sig_p = np.hypot(size * t.R_pc / pcas, sig_beam)
    m_msun = t.m_gas / len(flux)
    sig_v = max(np.sqrt(mo.G_PC * m_msun / max(size * t.R_pc, 1e-6)), 0.5 * t.dv)
    c = gaussian_filter(h / (t.pix_arcsec ** 2 * t.dv),
                        (sig_v / t.dv, sig_p / t.pix_arcsec, sig_p / t.pix_arcsec), mode="constant")
    c *= 2 * np.pi * sig_beam ** 2                      # Jy/beam
    nz, rms = noise(t, c.shape)
    c += nz
    # the metrics need |x| < 0.15 arcsec only: crop, the noise level is known
    k = np.abs(xs) <= CROP
    return dict(cube=c[:, k][:, :, k], xs=xs[k], vs=vs, rms=rms)


def metrics(m0, m1, m2, east, north, kind="profiles"):
    if kind == "basic":
        return np.array(F.metrics(m0, m1, m2, east, north))
    E, N_ = np.meshgrid(east, north)
    r = np.hypot(E, N_)
    m0 = np.where(np.isfinite(m0), m0, 0.0)
    tot = m0[r < EDGES[-1]].sum()
    share, v1, v2 = [], [], []
    for lo, hi in zip(EDGES[:-1], EDGES[1:]):
        s = (r >= lo) & (r < hi)
        share.append(m0[s].sum() / tot if tot > 0 else np.nan)
        k = s & (m0 > 0) & np.isfinite(m1) & np.isfinite(m2)
        v1.append(np.percentile(np.abs(m1[k]), 90) if k.sum() > 20 else np.nan)
        v2.append(np.median(m2[k]) if k.sum() > 20 else np.nan)
    # axis ratio within the outer edge of the annuli (intensity-weighted second moments)
    w = np.where((r < EDGES[-1]) & (m0 > 0), m0, 0.0)
    if w.sum() > 0:
        xe, yn = (w * E).sum() / w.sum(), (w * N_).sum() / w.sum()
        cxy = (w * (E - xe) * (N_ - yn)).sum()
        lam = np.linalg.eigvalsh([[(w * (E - xe) ** 2).sum(), cxy], [cxy, (w * (N_ - yn) ** 2).sum()]])
        ba = float(np.sqrt(lam[0] / lam[1]))
    else:
        ba = np.nan
    return np.array(share + v1 + v2 + [ba])


def chi2(mod, obs, kind="profiles"):
    if kind == "basic":
        return F.chi2(mod, obs)
    n = len(EDGES) - 1
    sig = np.r_[np.full(n, SIG["m0"]), np.full(n, SIG["m1"]), np.full(n, SIG["m2"])]
    c = 0.0
    for k in range(3 * n):
        if not np.isfinite(obs[k]) or obs[k] <= 0:
            continue
        c += MISS if not np.isfinite(mod[k]) or mod[k] <= 0 else (np.log(mod[k] / obs[k]) / sig[k]) ** 2
    c += MISS if not np.isfinite(mod[-1]) else ((mod[-1] - obs[-1]) / SIG["ba"]) ** 2
    return float(c)


def observed(kind):
    d = np.load(F.CUBE_DIR / "ngc613_mask_moments.npz")
    return metrics(d["m0"], d["m1"], d["m2"], d["east"], d["north"], kind)


def model_metrics(par, src, kind):
    t = replace(F.target(), R_pc=float(par["R_pc"]), inc=float(par["inc"]), M_bh=10. ** float(par["logM"]))
    if src in ("torus", "torus+disc"):
        p = torus(par["R0"], par["e_max"], par["i_max"])
        if src == "torus+disc":
            # same particle mass: the disc takes f_d of the N clouds from the torus
            nd = int(round(par["f_d"] * N))
            if nd > 0:
                xd, vd = disc(par["R_d"], par["h_d"], par["tilt_d"], nd)
                p["x"][:nd], p["v"][:nd] = xd, vd
        cb = fast_cube(p, t, 0.0)
    else:
        if src not in _snap:
            _snap[src] = mo.load(F.last_dump(src))
        p = _snap[src]
        cb = mo.cube(p, t, mo.m1_phase(p["x"], p["m"]) + float(par["az"]))
    return metrics(*mo.moments(cb, t, rms=cb.get("rms")), cb["xs"], cb["xs"], kind), cb, t


def objective(x, names, fixed, src, kind, obs):
    par = dict(fixed, **dict(zip(names, x)))
    return chi2(model_metrics(par, src, kind)[0], obs, kind)


class Recorder:
    """Keeps every evaluation of the objective (the population of each generation)."""

    def __init__(self, names):
        self.names, self.x, self.f, self.gen = names, [], [], 0

    def __call__(self, intermediate_result):
        res = intermediate_result
        pop = getattr(res, "population", None)
        en = getattr(res, "population_energies", None)
        if pop is not None:
            self.x += list(pop)
            self.f += list(en)
        self.x.append(res.x)
        self.f.append(res.fun)
        self.gen += 1
        print(f"  generation {self.gen}: "
              f"best chi2 {res.fun:.2f} at " + ", ".join(f"{n} {v:.3g}" for n, v in zip(self.names, res.x)),
              flush=True)


def check():
    """fast_cube against mockobs.cube on the same torus and view."""
    p = torus(0.3, 0.5, 45.)
    for rpc, inc, m in ((8., 45., 2.5e7), (4., 20., 3.75e7)):
        t = replace(F.target(), R_pc=rpc, inc=inc, M_bh=m)
        t0 = time.time()
        a = metrics(*mo.moments(mo.cube(p, t, 0.), t), *(2 * [mo.cube(p, t, 0.)["xs"]]), "profiles")
        t1 = time.time()
        cb = fast_cube(p, t, 0.)
        b = metrics(*mo.moments(cb, t, rms=cb["rms"]), cb["xs"], cb["xs"], "profiles")
        t2 = time.time()
        print(f"R_pc {rpc}, i {inc}: slow {(t1 - t0) / 2:.2f} s, fast {t2 - t1:.2f} s")
        print("  slow", np.round(a, 3))
        print("  fast", np.round(b, 3))


def figures(out, best, src, kind, obs, X, fvals, names, bounds, label):
    m, cb, t = model_metrics(best, src, kind)
    # profiles: data against the best model
    if kind == "profiles":
        n = len(EDGES) - 1
        rc = 0.5 * (np.array(EDGES[:-1]) + np.array(EDGES[1:]))  # annulus centres
        fig, axs = plt.subplots(1, 3, figsize=(13, 4), layout="constrained")
        for ax, (k, lab) in zip(axs, ((0, "share of moment 0"), (1, r"$|$moment 1$|$, 90th percentile [km/s]"),
                                      (2, "moment 2, median [km/s]"))):
            ax.plot(rc, obs[k * n:(k + 1) * n], "ko-", label="NGC 613, ALMA 2019.1.00273.S archive cube")
            ax.plot(rc, m[k * n:(k + 1) * n], "s--", color="C3", label="best model, same processing")
            ax.set(xlabel="radius [arcsec]", ylabel=lab)
        axs[0].legend(fontsize=8, frameon=False)
        fig.suptitle(f"{label}\nmoment-0 axis ratio: NGC 613 {obs[-1]:.2f}, model {m[-1]:.2f}; "
                     rf"$\chi^2 = {chi2(m, obs, kind):.1f}$ for {len(obs)} numbers", fontsize=10)
        fig.savefig(out / "profiles.png", dpi=130)
        plt.close(fig)
    # pairs of free parameters: every evaluation coloured by chi2, the best on top
    k = len(names)
    if k >= 2:
        o = np.argsort(fvals)[::-1]
        fig, axs = plt.subplots(k - 1, k - 1, figsize=(2.4 * (k - 1) + 1.5, 2.4 * (k - 1)),
                                layout="constrained", squeeze=False)
        for i in range(1, k):
            for j in range(k - 1):
                ax = axs[i - 1, j]
                if j >= i:
                    ax.set_axis_off()
                    continue
                sc = ax.scatter(X[o, j], X[o, i], c=np.log10(np.array(fvals)[o]), s=4, cmap="viridis_r")
                ax.plot(best[names[j]], best[names[i]], "r*", ms=10)
                ax.set_xlim(bounds[j])
                ax.set_ylim(bounds[i])
                if i == k - 1:
                    ax.set_xlabel(names[j])
                if j == 0:
                    ax.set_ylabel(names[i])
        fig.colorbar(sc, ax=axs, label=r"$\log_{10}\chi^2$")
        fig.suptitle(f"{label}: every evaluation of the search", fontsize=10)
        fig.savefig(out / "pairs.png", dpi=120)
        plt.close(fig)
    # maps: the best model next to the archive cube, palettes of Combes et al. (2025)
    from compare_obs import PAL
    d = np.load(F.CUBE_DIR / "ngc613_mask_moments.npz")
    mm = mo.moments(cb, t, rms=cb.get("rms"))
    fig, axs = plt.subplots(2, 3, figsize=(13, 8.5), layout="constrained")
    for col, key in enumerate(("m0", "m1", "m2")):
        cmap, lo, hi = PAL[key]
        for row, (img, ex, east) in enumerate(((mm[col], cb["xs"], cb["xs"]), (d[key], d["east"], d["north"]))):
            m0 = mm[0] if row == 0 else d["m0"]
            ax = axs[row, col]
            im = ax.imshow(np.where(m0 > 0, img, np.nan), origin="lower", cmap=cmap, vmin=lo, vmax=hi,
                           extent=[ex[0], ex[-1], east[0], east[-1]], interpolation="nearest")
            ax.set(xlim=(0.15, -0.15), ylim=(-0.15, 0.15), xlabel="RA offset [arcsec]", ylabel="Dec offset [arcsec]")
            ax.set_title(("best model" if row == 0 else "NGC 613, archive cube") + f", moment {col}", fontsize=9)
        fig.colorbar(im, ax=list(axs[:, col]), orientation="horizontal", fraction=0.04,
                     label=("moment 0 [Jy/beam km/s]", "moment 1 [km/s]", "moment 2 [km/s]")[col])
    fig.suptitle(label, fontsize=10)
    fig.savefig(out / "maps.png", dpi=110)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--model", default="torus")
    ap.add_argument("--free", default=None, help="comma list; default all of the model source")
    ap.add_argument("--fix", nargs="*", default=[], help="name=value")
    ap.add_argument("--bounds", nargs="*", default=[], help="name=lo:hi")
    ap.add_argument("--metrics", default="profiles", choices=("profiles", "basic"))
    ap.add_argument("--edges", default=None, help="annulus edges [arcsec], e.g. 0.01,0.025,0.04,0.06,0.08,0.1")
    ap.add_argument("--maxiter", type=int, default=40)
    ap.add_argument("--popsize", type=int, default=12)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default=None)
    ap.add_argument("--check", action="store_true", help="fast cube against mockobs.cube, then stop")
    a = ap.parse_args()
    if a.check:
        check()
        return
    global EDGES
    if a.edges:
        EDGES = tuple(map(float, a.edges.split(",")))
    src = a.model if a.model in ("torus", "torus+disc") else str(Path(a.model.split(":", 1)[1]).resolve())
    allowed = {"torus": ["R0", "e_max", "i_max", "R_pc", "inc", "logM"],
               "torus+disc": ["R0", "e_max", "i_max", "R_d", "f_d", "h_d", "tilt_d", "R_pc", "inc", "logM"]
               }.get(src, ["R_pc", "inc", "logM", "az"])
    names = a.free.split(",") if a.free else allowed
    assert set(names) <= set(allowed), f"free parameters of this model: {allowed}"
    fixed = {k: v for k, v in DEFAULT.items()}
    for kv in a.fix:
        k, v = kv.split("=")
        fixed[k] = float(v)
    bnd = dict(BOUNDS)
    for kv in a.bounds:
        k, v = kv.split("=")
        bnd[k] = tuple(map(float, v.split(":")))
    bounds = [bnd[n] for n in names]
    tag = src.replace("+", "_") if src in ("torus", "torus+disc") else Path(src).parent.name + "_" + Path(src).name
    out = Path(a.out or ROOT / "results" / "m1torus" / "B" / "fit_torus" / f"{tag}_{a.metrics}")
    out.mkdir(parents=True, exist_ok=True)
    obs = observed(a.metrics)
    print(f"model {src}, free {names}, metrics {a.metrics} ({len(obs)} numbers)")
    rec = Recorder(names)
    t0 = time.time()
    res = differential_evolution(partial(objective, names=names, fixed=fixed, src=src, kind=a.metrics, obs=obs),
                                 bounds, maxiter=a.maxiter, popsize=a.popsize, seed=a.seed, tol=1e-3,
                                 workers=a.workers, updating="deferred", polish=False, callback=rec)
    X, fv = np.array(rec.x), np.array(rec.f)
    best = dict(fixed, **dict(zip(names, res.x)))
    dof = max(len(obs) - len(names), 1)
    delta = max(res.fun / dof, 1.0)          # Delta chi2 for errors rescaled to chi2/dof = 1
    near = fv <= res.fun + delta
    ranges = {n: [float(X[near, i].min()), float(X[near, i].max())] for i, n in enumerate(names)}
    edge = [n for i, n in enumerate(names)
            if min(res.x[i] - bounds[i][0], bounds[i][1] - res.x[i]) < 0.02 * (bounds[i][1] - bounds[i][0])]
    summary = dict(model=src, metrics=a.metrics, n_numbers=len(obs), free=names, chi2=float(res.fun),
                   chi2_per_dof=float(res.fun / dof), best={k: float(v) for k, v in best.items()},
                   M_bh=float(10 ** best["logM"]), ranges_delta_chi2=delta, ranges=ranges, at_bound=edge,
                   evaluations=int(res.nfev), minutes=round((time.time() - t0) / 60., 1),
                   observed=list(map(float, obs)))
    (out / "fit.json").write_text(json.dumps(summary, indent=1))
    np.savez_compressed(out / "evaluations.npz", names=names, x=X, chi2=fv)
    if src in ("torus", "torus+disc"):
        label = (rf"Collisionless torus (phase-mixed, no mode): $R_0 = {best['R0']:.2f}$, "
                 rf"$e_{{\max}} = {best['e_max']:.2f}$, $i_{{\max}} = {best['i_max']:.0f}^\circ$; ")
        if src == "torus+disc":
            label += (rf"central disc: $R_d = {best['R_d']:.2f}\,R_{{\rm tor}}$, {100 * best['f_d']:.0f}% of the mass, "
                      rf"$\sigma/v_{{\rm K}} = {best['h_d']:.2f}$, tilt ${best['tilt_d']:.0f}^\circ$;\n")
    else:
        label = L.wrapped(L.model_label(Path(src).parent), 110) + "\n"
    label += (rf"$R_{{\rm tor}} = {best['R_pc']:.1f}$ pc, $i = {best['inc']:.0f}^\circ$, "
              rf"$M_{{\rm BH}} = {10 ** best['logM'] / 1e7:.2f}\times10^7\,M_\odot$")
    figures(out, best, src, a.metrics, obs, X, fv, names, bounds, label)
    print(json.dumps(summary, indent=1))
    if edge:
        print(f"WARNING: best value at a bound for {edge}: widen --bounds")


if __name__ == "__main__":
    main()
