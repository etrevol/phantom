#!/usr/bin/env python3
"""B19, B26: fit a collisionless torus (or a run's snapshot) to a galaxy in its ALMA archive cube.

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
from galaxies import GALAXIES  # noqa: E402
import labels as L  # noqa: E402
import mockobs as mo  # noqa: E402

N = 16384
SIZE = 0.05                     # particle size [R_tor] of generated tori (h/1.4 of the 16k runs)
EDGES = (0.01, 0.025, 0.04, 0.06, 0.08, 0.11, 0.15)    # annuli [arcsec]
SIG = dict(m0=0.2, m1=0.15, m2=0.15, ba=0.1)  # log errors (m0, m1, m2) and absolute (b/a)
MISS = 25.                      # chi2 of an annulus the model leaves empty
SIG_PV = 20.                    # error of a PV velocity percentile [km/s] (~3 channels)
POS = (0.1, 0.2, 0.35, 0.5, 0.7, 0.9)   # PV positions along the major axis [x outer annulus edge]
DEFAULT = dict(R0=0.3, e_max=0.5, i_max=60., R_pc=8., inc=50., logM=np.log10(3.4e7), az=90.,
               R_d=0.3, f_d=0.2, h_d=0.1, tilt_d=0., f_ext=0., beta=1.5,
               dx=0., dy=0., dv=0., pa=np.nan,
               a_in=0.3, e0=0.1, e1=0.1, dw=0., warp_i=0., warp_n=0.)
# geometry of the data: black hole at (dx, dy) [arcsec east, north] of the map centre,
# systemic velocity dv [km/s] above the adopted one, sky PA of the model (nan: galaxies.py);
# the data metrics are taken about that centre and velocity, the model's about its own
GEO = ("dx", "dy", "dv", "pa")
BOUNDS = dict(R0=(0.1, 0.8), e_max=(0.05, 0.8), i_max=(5., 89.), R_pc=(1.5, 10.),
              inc=(5., 80.), logM=(6.9, 8.0), az=(0., 360.),
              R_d=(0.05, 0.8), f_d=(0.0, 0.7), h_d=(0.02, 0.4), tilt_d=(-90., 90.),
              f_ext=(0.0, 5.0), beta=(0.5, 3.0),
              dx=(-0.03, 0.03), dy=(-0.03, 0.03), dv=(-60., 60.), pa=(None, None),
              a_in=(0.05, 0.9), e0=(0., 0.7), e1=(0., 0.7), dw=(-90., 90.), warp_i=(-45., 45.),
              warp_n=(0., 180.))
GAL = "ngc613"                  # the galaxy fitted (galaxies.py), set by --galaxy
_noise = {}
_snap = {}
_data = {}


def extend(p, f_ext, beta):
    """Extended mass (stars, gas) with the black hole: M(<r) = M_BH [1 + f_ext (r/R_tor)^beta];
    beta = 3 - gamma for a density ~ r^-gamma. Velocities scale by sqrt(M(<r)/M_BH): exact
    for circular orbits, approximate for eccentric ones."""
    if f_ext > 0:
        r = np.linalg.norm(p["x"], axis=1)
        p["v"] = p["v"] * np.sqrt(1. + f_ext * r ** beta)[:, None]
    return p


def target():
    """The galaxy's Target: distances and numbers from galaxies.py, beam, channels and
    noise of its archive cube, total flux of its masked moment 0 within the field."""
    g = GALAXIES[GAL]
    out = Path(g["out"])
    c = np.load(out / f"{GAL}_cube_win.npz")
    pix = float(abs(c["north"][1] - c["north"][0]))
    if GAL == "ngc613":
        flux = 15.4 * g["m_gas"] / 9.7e6       # Combes et al. (2025), Table 2, as mockobs.NGC613
    else:
        d = np.load(out / f"{GAL}_mask_moments.npz")
        E, Nn = np.meshgrid(d["east"], d["north"])
        bpix = np.pi * float(d["beam"]) ** 2 / (4 * np.log(2)) / pix ** 2
        flux = float(np.nansum(np.where(np.hypot(E, Nn) < g["field"], d["m0"], 0.0)) / bpix)
    return replace(mo.NGC613, name=g["name"], R_pc=g["R_pc"], M_bh=g["M_bh"], D_Mpc=g["D_Mpc"], inc=g["inc"],
                   pa=g["pa"], beam_arcsec=float(c["beam"]), pix_arcsec=pix,
                   dv=float(np.median(np.diff(c["v"]))), vmax=g["window"], rms_mjy=float(c["rms"]) * 1e3,
                   clip_mjy=g["clip_mjy"], field_arcsec=max(g["field"], EDGES[-1] + 0.05),
                   jy_kms_per_msun=flux / g["m_gas"], m_gas=g["m_gas"])


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


def eccdisc(a_in, e0, e1, dw, warp_i, warp_n, i_max, n=N, seed=3):
    """A lopsided, warped thin disc of Keplerian orbits about a unit mass (kinematic, not
    self-consistent; as the eccentric-disc models of M31, Peiris & Tremaine 2003): semi-
    major axes uniform in [a_in, 1] (Sigma ~ 1/a, sharp edges, an inner hole), eccentricity
    e0 -> e1 linearly from the inner to the outer edge, pericentres along x turned by dw
    [deg] across the disc, uniform mean anomaly (time average: crowding at apocentre), a
    thickness of random tilts within +-i_max, and a warp: each orbit tilted by warp_i x
    [deg] (x = 0 inner edge, 1 outer edge) about the line of nodes at warp_n [deg]."""
    rng = np.random.default_rng(seed)
    u = rng.random((5, n))
    a = a_in + (1. - a_in) * u[0]
    x = (a - a_in) / max(1. - a_in, 1e-6)
    e = np.clip(e0 + (e1 - e0) * x, 0., 0.95)
    M = 2. * np.pi * u[1]
    E = M.copy()
    for _ in range(40):
        E -= (E - e * np.sin(E) - M) / (1. - e * np.cos(E))
    nu = 2. * np.arctan2(np.sqrt(1. + e) * np.sin(E / 2.), np.sqrt(1. - e) * np.cos(E / 2.))
    p = a * (1. - e ** 2)
    r = p / (1. + e * np.cos(nu))
    w = np.radians(dw) * x
    th = nu + w
    pos = np.c_[r * np.cos(th), r * np.sin(th), np.zeros(n)]
    vr, vt = e * np.sin(nu) / np.sqrt(p), (1. + e * np.cos(nu)) / np.sqrt(p)
    vel = np.c_[vr * np.cos(th) - vt * np.sin(th), vr * np.sin(th) + vt * np.cos(th), np.zeros(n)]

    def tilt(v, ang, node):             # rotate each vector by ang about the in-plane axis at node
        k = np.c_[np.cos(node), np.sin(node), np.zeros_like(node)]
        c, s = np.cos(ang)[:, None], np.sin(ang)[:, None]
        return v * c + np.cross(k, v) * s + k * (k * v).sum(1)[:, None] * (1. - c)
    thick, tnode = np.radians((2. * u[2] - 1.) * i_max), 2. * np.pi * u[3]
    pos, vel = tilt(pos, thick, tnode), tilt(vel, thick, tnode)
    wa, wn = np.radians(warp_i) * x, np.full(n, np.radians(warp_n))
    pos, vel = tilt(pos, wa, wn), tilt(vel, wa, wn)
    return dict(x=pos, v=vel, m=np.full(n, 1. / n), size=np.full(n, SIZE), m_tor=1.)


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
    # the metrics need |x| < the outer annulus only: crop, the noise level is known
    k = np.abs(xs) <= EDGES[-1] + 0.03
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


def pv_metrics(cb, t, centre=(0.0, 0.0), dv=0.0):
    """Velocity percentiles (10th, 90th) of the flux in the major-axis PV diagram (slit one
    beam wide) at +-POS x the outer annulus edge: the rotation curve, its asymmetry and the
    line wings, robust to noise (PV smoothed by one beam and 1.5 channels, kept above 4
    sigma, sigma from the median absolute deviation)."""
    pvd = mo.pv(cb, t, pa=t.pa, halfwidth=t.beam_arcsec / 2, centre=centre)
    xs, vs = cb["xs"], cb["vs"] - dv
    sm = gaussian_filter(pvd, (1.5, t.beam_arcsec / mo.FWHM / t.pix_arcsec))
    sig = 1.4826 * np.median(np.abs(sm - np.median(sm)))
    out = []
    for sgn in (1, -1):
        for f in POS:
            k = int(np.argmin(np.abs(xs - sgn * f * EDGES[-1])))
            w = np.where(sm[:, k] > 4 * sig, sm[:, k], 0.0)
            if w.sum() <= 0:
                out += [np.nan, np.nan]
                continue
            cdf = np.cumsum(w) / w.sum()
            out += [float(np.interp(0.1, cdf, vs)), float(np.interp(0.9, cdf, vs))]
    return np.array(out)


SIG_M1 = (0.03, 0.15, 0.15)      # errors of A1 and of the centroid offset east, north [pc]:
                                 # data noise (B31: 0.014, 0.07) and clump noise of the models


def lopsided(m0, east, north, t):
    """The m = 1 numbers of moment 0 about the black hole at the map origin: A1 of the
    annulus 0.02 arcsec - outer edge, and the centroid offset east and north [pc]."""
    E, N_ = np.meshgrid(east, north)
    m0 = np.where(np.isfinite(m0), m0, 0.0)
    r = np.hypot(E, N_)
    w = np.where(r < EDGES[-1], m0, 0.0)
    if w.sum() <= 0:
        return np.array([np.nan, np.nan, np.nan])
    s = (r > 0.02) & (r < EDGES[-1])
    a1 = abs((m0[s] * np.exp(1j * np.arctan2(N_, E)[s])).sum()) / max(m0[s].sum(), 1e-30)
    pcas = mo.pc_per_arcsec(t)
    return np.array([a1, (w * E).sum() / w.sum() * pcas, (w * N_).sum() / w.sum() * pcas])


def chi2(mod, obs, kind="profiles"):
    if kind == "profiles+pv+m1":
        c = chi2(mod[:-3], obs[:-3], "profiles+pv")
        for a, b, sg in zip(mod[-3:], obs[-3:], SIG_M1):
            c += MISS if not np.isfinite(a) else ((a - b) / sg) ** 2
        return float(c)
    if kind == "basic":
        return F.chi2(mod, obs)
    n = len(EDGES) - 1
    if kind == "profiles+pv":
        npr = 3 * n + 1
        c = chi2(mod[:npr], obs[:npr], "profiles")
        for a, b in zip(mod[npr:], obs[npr:]):
            if np.isfinite(b):
                c += MISS if not np.isfinite(a) else ((a - b) / SIG_PV) ** 2
        return float(c)
    sig = np.r_[np.full(n, SIG["m0"]), np.full(n, SIG["m1"]), np.full(n, SIG["m2"])]
    c = 0.0
    for k in range(3 * n):
        if not np.isfinite(obs[k]) or obs[k] <= 0:
            continue
        c += MISS if not np.isfinite(mod[k]) or mod[k] <= 0 else (np.log(mod[k] / obs[k]) / sig[k]) ** 2
    c += MISS if not np.isfinite(mod[-1]) else ((mod[-1] - obs[-1]) / SIG["ba"]) ** 2
    return float(c)


def observed(kind, geo=None):
    """The data metrics, about the black hole at (dx, dy) with systemic offset dv and the
    model's sky PA when geo (a parameter dict) is given, else about the map centre."""
    if GAL not in _data:
        d = dict(np.load(Path(GALAXIES[GAL]["out"]) / f"{GAL}_mask_moments.npz"))
        c = np.load(Path(GALAXIES[GAL]["out"]) / f"{GAL}_cube_win.npz")
        # ascending east for the slit sampler, as compare_obs.py
        cb = dict(cube=c["cube"][:, :, ::-1].astype(float), xs=c["east"][::-1], vs=c["v"])
        _data[GAL] = (d, cb)
    d, cb = _data[GAL]
    g = dict(dx=0., dy=0., dv=0., pa=np.nan)
    g.update({k: geo[k] for k in GEO if geo and k in geo})
    m = metrics(d["m0"], d["m1"] - g["dv"], d["m2"], d["east"] - g["dx"], d["north"] - g["dy"],
                "basic" if kind == "basic" else "profiles")
    if not kind.startswith("profiles+pv"):
        return m
    t = target()
    if np.isfinite(g["pa"]):
        t = replace(t, pa=float(g["pa"]))
    m = np.r_[m, pv_metrics(cb, t, centre=(g["dx"], g["dy"]), dv=g["dv"])]
    if kind == "profiles+pv+m1":
        m = np.r_[m, lopsided(d["m0"], d["east"] - g["dx"], d["north"] - g["dy"], t)]
    return m


def model_metrics(par, src, kind):
    t = replace(target(), R_pc=float(par["R_pc"]), inc=float(par["inc"]), M_bh=10. ** float(par["logM"]))
    if np.isfinite(par.get("pa", np.nan)):
        t = replace(t, pa=float(par["pa"]))
    if src == "eccdisc":
        p = extend(eccdisc(par["a_in"], par["e0"], par["e1"], par["dw"], par["warp_i"], par["warp_n"],
                           par["i_max"]), par["f_ext"], par["beta"])
        cb = fast_cube(p, t, float(par["az"]))
    elif src in ("torus", "torus+disc"):
        p = extend(torus(par["R0"], par["e_max"], par["i_max"]), par["f_ext"], par["beta"])
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
    mom = mo.moments(cb, t, rms=cb.get("rms"))
    m = metrics(*mom, cb["xs"], cb["xs"], "basic" if kind == "basic" else "profiles")
    if kind.startswith("profiles+pv"):
        m = np.r_[m, pv_metrics(cb, t)]
    if kind == "profiles+pv+m1":
        m = np.r_[m, lopsided(mom[0], cb["xs"], cb["xs"], t)]
    return m, cb, t


def objective(x, names, fixed, src, kind, obs, prior=None):
    par = dict(fixed, **dict(zip(names, x)))
    # Gaussian priors from independent measurements (e.g. the kinematic centre, B34)
    pen = sum(((par[k] - mu) / sig) ** 2 for k, (mu, sig) in (prior or {}).items())
    if any(k in names for k in GEO) or any(par[k] != 0 for k in ("dx", "dy", "dv")) or np.isfinite(par["pa"]):
        obs = observed(kind, par)
    return chi2(model_metrics(par, src, kind)[0], obs, kind) + pen


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
    if kind in ("profiles", "profiles+pv"):
        n = len(EDGES) - 1
        rc = 0.5 * (np.array(EDGES[:-1]) + np.array(EDGES[1:]))  # annulus centres
        npan = 4 if kind == "profiles+pv" else 3
        fig, axs = plt.subplots(1, npan, figsize=(4.4 * npan, 4), layout="constrained")
        if kind == "profiles+pv":
            npr = 3 * n + 1
            x = np.r_[np.array(POS), -np.array(POS)] * EDGES[-1]
            ax = axs[3]
            for arr, st, lab in ((obs, "ko", "data"), (m, "s", "model")):
                q = arr[npr:].reshape(-1, 2)
                ax.vlines(x, q[:, 0], q[:, 1], color="k" if lab == "data" else "C3", lw=3 if lab == "data" else 1.5,
                          alpha=0.6, label=lab)
            ax.axhline(0, color="0.7", lw=0.8)
            ax.set(xlabel="position along the major axis [arcsec]", ylabel="PV velocities, 10-90% of the flux [km/s]")
            ax.invert_xaxis()
            ax.legend(fontsize=8, frameon=False)
        for ax, (k, lab) in zip(axs[:3], ((0, "share of moment 0"), (1, r"$|$moment 1$|$, 90th percentile [km/s]"),
                                      (2, "moment 2, median [km/s]"))):
            ax.plot(rc, obs[k * n:(k + 1) * n], "ko-", label=f"{GALAXIES[GAL]['name']}, ALMA 2019.1.00273.S archive cube")
            ax.plot(rc, m[k * n:(k + 1) * n], "s--", color="C3", label="best model, same processing")
            ax.set(xlabel="radius [arcsec]", ylabel=lab)
        axs[0].legend(fontsize=8, frameon=False)
        iba = 3 * n
        fig.suptitle(f"{label}\nmoment-0 axis ratio: {GALAXIES[GAL]['name']} {obs[iba]:.2f}, model {m[iba]:.2f}; "
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
    d = np.load(Path(GALAXIES[GAL]["out"]) / f"{GAL}_mask_moments.npz")
    mm = mo.moments(cb, t, rms=cb.get("rms"))
    lim = EDGES[-1]
    fig, axs = plt.subplots(2, 3, figsize=(13, 8.5), layout="constrained")
    for col, key in enumerate(("m0", "m1", "m2")):
        cmap, lo, hi = PAL[key]
        for row, (img, ex, east) in enumerate(((mm[col], cb["xs"], cb["xs"]), (d[key], d["east"], d["north"]))):
            m0 = mm[0] if row == 0 else d["m0"]
            ax = axs[row, col]
            im = ax.imshow(np.where(m0 > 0, img, np.nan), origin="lower", cmap=cmap, vmin=lo, vmax=hi,
                           extent=[ex[0], ex[-1], east[0], east[-1]], interpolation="nearest")
            ax.set(xlim=(lim, -lim), ylim=(-lim, lim), xlabel="RA offset [arcsec]", ylabel="Dec offset [arcsec]")
            ax.set_title(("best model" if row == 0 else f"{GALAXIES[GAL]['name']}, archive cube") + f", moment {col}",
                         fontsize=9)
        fig.colorbar(im, ax=list(axs[:, col]), orientation="horizontal", fraction=0.04,
                     label=("moment 0 [Jy/beam km/s]", "moment 1 [km/s]", "moment 2 [km/s]")[col])
    fig.suptitle(label, fontsize=10)
    fig.savefig(out / "maps.png", dpi=110)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--galaxy", default="ngc613", choices=sorted(GALAXIES))
    ap.add_argument("--model", default="torus")
    ap.add_argument("--free", default=None, help="comma list; default all of the model source")
    ap.add_argument("--fix", nargs="*", default=[], help="name=value")
    ap.add_argument("--bounds", nargs="*", default=[], help="name=lo:hi")
    ap.add_argument("--prior", nargs="*", default=[], help="name=mean:sigma, a Gaussian prior added to chi2")
    ap.add_argument("--metrics", default="profiles", choices=("profiles", "profiles+pv", "profiles+pv+m1", "basic"))
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
    global EDGES, GAL
    GAL = a.galaxy
    if a.edges:
        EDGES = tuple(map(float, a.edges.split(",")))
    src = a.model if a.model in ("torus", "torus+disc", "eccdisc") else str(Path(a.model.split(":", 1)[1]).resolve())
    allowed = {"torus": ["R0", "e_max", "i_max", "R_pc", "inc", "logM", "f_ext", "beta"],
               "torus+disc": ["R0", "e_max", "i_max", "R_d", "f_d", "h_d", "tilt_d", "R_pc", "inc", "logM"],
               "eccdisc": ["a_in", "e0", "e1", "dw", "warp_i", "warp_n", "i_max", "R_pc", "inc", "logM", "az",
                           "f_ext", "beta"]
               }.get(src, ["R_pc", "inc", "logM", "az"]) + list(GEO)
    names = a.free.split(",") if a.free else [n for n in allowed if n not in GEO]
    assert set(names) <= set(allowed), f"free parameters of this model: {allowed}"
    fixed = {k: v for k, v in DEFAULT.items()}
    for kv in a.fix:
        k, v = kv.split("=")
        fixed[k] = float(v)
    bnd = dict(BOUNDS)
    pa0 = GALAXIES[GAL]["pa"]
    bnd["pa"] = (pa0 - 45., pa0 + 45.)
    for kv in a.bounds:
        k, v = kv.split("=")
        bnd[k] = tuple(map(float, v.split(":")))
    bounds = [bnd[n] for n in names]
    prior = {kv.split("=")[0]: tuple(map(float, kv.split("=")[1].split(":"))) for kv in a.prior}
    tag = src.replace("+", "_") if src in ("torus", "torus+disc", "eccdisc") else F.last_dump(src).parent.name + "_" + F.last_dump(src).name
    out = Path(a.out or ROOT / "results" / "m1torus" / "B" / "fit_torus" / f"{GAL}_{tag}_{a.metrics}")
    out.mkdir(parents=True, exist_ok=True)
    obs = observed(a.metrics)
    print(f"model {src}, free {names}, metrics {a.metrics} ({len(obs)} numbers)")
    rec = Recorder(names)
    t0 = time.time()
    res = differential_evolution(partial(objective, names=names, fixed=fixed, src=src, kind=a.metrics, obs=obs, prior=prior),
                                 bounds, maxiter=a.maxiter, popsize=a.popsize, seed=a.seed, tol=1e-3,
                                 workers=a.workers, updating="deferred", polish=False, callback=rec)
    X, fv = np.array(rec.x), np.array(rec.f)
    best = dict(fixed, **dict(zip(names, res.x)))
    dof = max(len(obs) - len(names), 1)
    delta = max(res.fun / dof, 1.0)          # Delta chi2 for errors rescaled to chi2/dof = 1
    near = fv <= res.fun + delta
    ranges = {n: [float(X[near, i].min()), float(X[near, i].max())] for i, n in enumerate(names)}
    obs = observed(a.metrics, best)           # the data about the fitted centre
    edge = [n for i, n in enumerate(names)
            if min(res.x[i] - bounds[i][0], bounds[i][1] - res.x[i]) < 0.02 * (bounds[i][1] - bounds[i][0])]
    pen = sum(((best[k] - mu) / sig) ** 2 for k, (mu, sig) in prior.items())
    summary = dict(model=src, metrics=a.metrics, n_numbers=len(obs), free=names, chi2=float(res.fun),
                   prior={k: list(v) for k, v in prior.items()}, chi2_prior=float(pen), chi2_data=float(res.fun - pen),
                   chi2_per_dof=float(res.fun / dof), best={k: float(v) for k, v in best.items()},
                   M_bh=float(10 ** best["logM"]), ranges_delta_chi2=delta, ranges=ranges, at_bound=edge,
                   evaluations=int(res.nfev), minutes=round((time.time() - t0) / 60., 1),
                   observed=list(map(float, obs)))
    (out / "fit.json").write_text(json.dumps(summary, indent=1))
    np.savez_compressed(out / "evaluations.npz", names=names, x=X, chi2=fv)
    if src == "eccdisc":
        label = (rf"Eccentric warped disc (kinematic): $a = {best['a_in']:.2f}-1$, $e = {best['e0']:.2f} \to {best['e1']:.2f}$, "
                 rf"twist ${best['dw']:.0f}^\circ$, warp ${best['warp_i']:.0f}^\circ$ at ${best['warp_n']:.0f}^\circ$, "
                 rf"$i_{{\max}} = {best['i_max']:.0f}^\circ$; ")
    elif src in ("torus", "torus+disc"):
        label = (rf"Collisionless torus (phase-mixed, no mode): $R_0 = {best['R0']:.2f}$, "
                 rf"$e_{{\max}} = {best['e_max']:.2f}$, $i_{{\max}} = {best['i_max']:.0f}^\circ$; ")
        if src == "torus+disc":
            label += (rf"central disc: $R_d = {best['R_d']:.2f}\,R_{{\rm tor}}$, {100 * best['f_d']:.0f}% of the mass, "
                      rf"$\sigma/v_{{\rm K}} = {best['h_d']:.2f}$, tilt ${best['tilt_d']:.0f}^\circ$;\n")
    else:
        label = L.wrapped(L.model_label(F.last_dump(src).parent), 110) + "\n"
    if best.get("f_ext", 0) > 0:
        label += rf"extended mass $f_{{\rm ext}} = {best['f_ext']:.2f}$, $\beta = {best['beta']:.1f}$; "
    label += (rf"$R_{{\rm tor}} = {best['R_pc']:.1f}$ pc, $i = {best['inc']:.0f}^\circ$, "
              rf"$M_{{\rm BH}} = {10 ** best['logM'] / 1e7:.2f}\times10^7\,M_\odot$")
    if any(k in names for k in GEO):
        label += (rf"; centre $({1e3 * best['dx']:+.0f}, {1e3 * best['dy']:+.0f})$ mas, "
                  rf"$\Delta v_{{\rm sys}} = {best['dv']:+.0f}$ km/s, PA ${best['pa'] if np.isfinite(best['pa']) else pa0:.0f}^\circ$")
    m1k = a.metrics == "profiles+pv+m1"
    if m1k:
        summary["m1_data"] = list(map(float, obs[-3:]))
        summary["m1_model"] = list(map(float, model_metrics(best, src, a.metrics)[0][-3:]))
        (out / "fit.json").write_text(json.dumps(summary, indent=1))
    figures(out, best, src, "profiles+pv" if m1k else a.metrics, obs[:-3] if m1k else obs, X, fv, names, bounds, label)
    print(json.dumps(summary, indent=1))
    if edge:
        print(f"WARNING: best value at a bound for {edge}: widen --bounds")


if __name__ == "__main__":
    main()
