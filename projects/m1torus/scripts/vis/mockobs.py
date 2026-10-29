#!/usr/bin/env python3
"""Synthetic ALMA observation of a m1torus dump, in the units of a target galaxy.

Each particle emits in proportion to its mass (optically thin, uniform excitation),
spread over its own size (cloud radius, or h/2 for SPH gas) and an internal line width,
then seen through a Gaussian beam with the target's inclination, position angle,
distance, pixel, channel width and noise. From the cube: moments 0, 1, 2 with the
observers' clip, a PV diagram along the major axis, the integrated spectrum, the hole
(deepest minimum of moment 0 inside the torus) and A1(R) of the deprojected moment 0.

Units: code G = M_BH = 1, R = 1 is the target's R_unit; the AGN is the sink. The
observer azimuth phi (the direction of the line of sight projected on the torus plane,
measured from the code x axis) is free; --phi mode puts it relative to the m = 1 phase.

    mockobs.py DUMP [--target ngc613|circinus] [--phi DEG] [--out DIR]
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "common"))
import phantom_io as pio  # noqa: E402

G_PC = 4.30091e-3          # G in pc (km/s)^2 / Msun
ISTAR, IGAS = 4, 1
FWHM = 2.0 * np.sqrt(2.0 * np.log(2.0))


@dataclass
class Target:
    name: str
    R_pc: float            # code length unit
    M_bh: float            # code mass unit [Msun]
    D_Mpc: float
    inc: float             # deg, 0 = face-on
    pa: float              # deg E of N, receding side of the kinematic major axis
    beam_arcsec: float     # FWHM
    pix_arcsec: float
    dv: float              # channel width [km/s]
    vmax: float            # cube half-range [km/s]
    rms_mjy: float         # noise per channel [mJy/beam]
    clip_mjy: float        # moment clip [mJy/beam]
    field_arcsec: float    # half-size of the map
    jy_kms_per_msun: float  # flux calibration
    m_gas: float           # gas mass the torus stands for [Msun]
    pv_halfwidth_arcsec: float


# NGC 613: Combes et al. (2025): D 17.2 Mpc, i 50, kinematic PA -20 (their moment 1
# recedes to the north),
# beam 0.015", 5.9 km/s channels, rms 0.4 mJy/beam, moments clipped at 0.6 mJy/beam,
# 15.4 Jy km/s for 9.7e6 Msun (Table 2), M_BH 3.4e7, torus radius 8 pc.
NGC613 = Target("NGC 613", 8.0, 3.4e7, 17.2, 50.0, 340.0, 0.015, 0.0029, 5.9, 300.0,
                0.4, 0.6, 0.3, 15.4 / 9.7e6, 1.0e7, 0.0075)
# Circinus: Izumi et al. (2023), Tristram et al. (2022), Lopez-Rodriguez et al. (2026):
# D 4.2 Mpc, edge-on maser disc, PA ~ 50, M_BH 1.7e6, dust disc 5 x 3 pc; beam ~0.1".
CIRCINUS = Target("Circinus", 2.5, 1.7e6, 4.2, 88.0, 50.0, 0.10, 0.01, 5.0, 200.0,
                  0.4, 1.2, 0.35, 1.0e-5, 3.0e5, 0.05)
TARGETS = {"ngc613": NGC613, "circinus": CIRCINUS}


def pc_per_arcsec(t: Target) -> float:
    return t.D_Mpc * 1.0e6 * np.pi / (180.0 * 3600.0)


def load(dump):
    """Particles relative to the sink, in code units, with per-particle size."""
    d = pio.read_dump(dump)
    g, s = d["gas"], d["sinks"]
    mtype = d["header_lists"]["massoftype"]
    itype = g["itype"] if "itype" in g else np.full(len(g["x"]), ISTAR)
    alive = g["h"] > 0
    m = np.array([mtype[it - 1] for it in itype])
    x = np.c_[g["x"], g["y"], g["z"]] - [s["x"][0], s["y"][0], s["z"][0]]
    v = np.c_[g["vx"], g["vy"], g["vz"]] - [s["vx"][0], s["vy"][0], s["vz"][0]]
    hdr = d["header"]
    # clouds: their radius (softening h = 1.4 eps = 1.4 r_cloud); gas: h/2
    size = np.where(itype == IGAS, 0.5 * g["h"], np.abs(g["h"]) / 1.4)
    return dict(x=x[alive], v=v[alive], m=m[alive], size=size[alive], time=hdr["time"],
                m_tor=m[alive].sum())


def m1_phase(x, m):
    """Phase of the m = 1 surface-density harmonic in the torus plane."""
    ph = np.arctan2(x[:, 1], x[:, 0])
    return np.degrees(np.arctan2((m * np.sin(ph)).sum(), (m * np.cos(ph)).sum()))


def to_sky(x, v, inc, pa, phi):
    """Code frame -> sky: E (left positive), N, line of sight (away), and v_los.

    The observer sits at azimuth phi in the torus plane, inclined by inc from the pole;
    the projected angular-momentum axis then points to PA - 90 so that the receding
    half of the major axis lies at PA.
    """
    c, s = np.cos(np.radians(phi)), np.sin(np.radians(phi))
    # frame with x' toward the observer's azimuth
    xr = np.c_[c * x[:, 0] + s * x[:, 1], -s * x[:, 0] + c * x[:, 1], x[:, 2]]
    vr = np.c_[c * v[:, 0] + s * v[:, 1], -s * v[:, 0] + c * v[:, 1], v[:, 2]]
    ci, si = np.cos(np.radians(inc)), np.sin(np.radians(inc))
    # line of sight n = (si, 0, ci) toward the observer; sky axes a = y', b = (ci, 0, -si)
    los = -(xr[:, 0] * si + xr[:, 2] * ci)            # distance away from the observer
    vlos = -(vr[:, 0] * si + vr[:, 2] * ci)           # positive = receding
    a = xr[:, 1]
    b = xr[:, 0] * ci - xr[:, 2] * si
    # rotate (a, b) so the receding side of the major axis (a > 0 for prograde) is at PA
    t = np.radians(pa)
    east = a * np.sin(t) + b * np.cos(t)
    north = a * np.cos(t) - b * np.sin(t)
    return east, north, los, vlos


def cube(p, t: Target, phi, seed=1):
    """Spectral cube [Jy/beam] on a grid of (v, north, east), with beam and noise."""
    pcas = pc_per_arcsec(t)
    vunit = np.sqrt(G_PC * t.M_bh / t.R_pc)
    east, north, _, vlos = to_sky(p["x"], p["v"], t.inc, t.pa, phi)
    east_as = east * t.R_pc / pcas
    north_as = north * t.R_pc / pcas
    v = vlos * vunit
    nx = int(round(2 * t.field_arcsec / t.pix_arcsec)) | 1
    nv = int(round(2 * t.vmax / t.dv)) | 1
    xs = (np.arange(nx) - nx // 2) * t.pix_arcsec
    vs = (np.arange(nv) - nv // 2) * t.dv
    # flux of a particle: its share of the gas mass the torus stands for
    flux = p["m"] / p["m_tor"] * t.m_gas * t.jy_kms_per_msun        # Jy km/s
    sig_beam = t.beam_arcsec / FWHM
    sig_p = np.sqrt((p["size"] * t.R_pc / pcas) ** 2 + sig_beam ** 2)
    # internal line width: a self-gravitating cloud, sqrt(G m / r), at least half a channel
    m_msun = p["m"] / p["m_tor"] * t.m_gas
    sig_v = np.maximum(np.sqrt(G_PC * m_msun / np.maximum(p["size"] * t.R_pc, 1e-6)), 0.5 * t.dv)
    c = np.zeros((nv, nx, nx))
    ix0 = np.searchsorted(xs, east_as)
    iy0 = np.searchsorted(xs, north_as)
    iv0 = np.searchsorted(vs, v)
    for k in range(len(flux)):
        rx = int(np.ceil(3 * sig_p[k] / t.pix_arcsec)) + 1
        rv = int(np.ceil(3 * sig_v[k] / t.dv)) + 1
        x0, x1 = max(ix0[k] - rx, 0), min(ix0[k] + rx, nx)
        y0, y1 = max(iy0[k] - rx, 0), min(iy0[k] + rx, nx)
        v0, v1 = max(iv0[k] - rv, 0), min(iv0[k] + rv, nv)
        if x0 >= x1 or y0 >= y1 or v0 >= v1:
            continue
        gx = np.exp(-0.5 * ((xs[x0:x1] - east_as[k]) / sig_p[k]) ** 2)
        gy = np.exp(-0.5 * ((xs[y0:y1] - north_as[k]) / sig_p[k]) ** 2)
        gv = np.exp(-0.5 * ((vs[v0:v1] - v[k]) / sig_v[k]) ** 2)
        # normalised so that the stamp integrates to the particle's flux per beam
        norm = flux[k] / (2 * np.pi * sig_p[k] ** 2) / (np.sqrt(2 * np.pi) * sig_v[k])
        c[v0:v1, y0:y1, x0:x1] += norm * gv[:, None, None] * gy[None, :, None] * gx[None, None, :]
    beam_area = 2 * np.pi * sig_beam ** 2                      # arcsec^2
    c *= beam_area                                             # Jy/beam
    # noise correlated over the beam, rms as observed
    rng = np.random.default_rng(seed)
    n = gaussian_filter(rng.standard_normal(c.shape), (0, sig_beam / t.pix_arcsec, sig_beam / t.pix_arcsec))
    c += n / n.std() * t.rms_mjy * 1e-3
    return dict(cube=c, xs=xs, vs=vs, east=east_as, north=north_as, v=v, vunit=vunit)


def moments(cb, t: Target):
    c, vs = cb["cube"], cb["vs"]
    w = np.where(c > t.clip_mjy * 1e-3, c, 0.0)
    m0 = w.sum(0) * t.dv
    with np.errstate(invalid="ignore", divide="ignore"):
        m1 = (w * vs[:, None, None]).sum(0) * t.dv / m0
        m2 = np.sqrt((w * (vs[:, None, None] - m1) ** 2).sum(0) * t.dv / m0)
    m1[m0 <= 0] = np.nan
    m2[m0 <= 0] = np.nan
    return m0, m1, m2


def pv(cb, t: Target, pa=None, halfwidth=None):
    """Position-velocity slice along PA through the AGN, summed over the slit width."""
    pa = t.pa if pa is None else pa
    hw = t.pv_halfwidth_arcsec if halfwidth is None else halfwidth
    xs = cb["xs"]
    th = np.radians(pa)
    ds = xs[1] - xs[0]
    offs = np.arange(-hw, hw + 1e-9, ds)
    out = np.zeros((len(cb["vs"]), len(xs)))
    for o in offs:
        e = xs * np.sin(th) + o * np.cos(th)
        nn = xs * np.cos(th) - o * np.sin(th)
        ie = (e - xs[0]) / ds
        jn = (nn - xs[0]) / ds
        for iv in range(len(cb["vs"])):
            out[iv] += map_coordinates(cb["cube"][iv], [jn, ie], order=1, mode="constant")
    return out / len(offs)


def hole(m0, xs, t: Target, rmax_pc=None):
    """Deepest minimum of the smoothed moment 0 inside the torus: offset [pc] and PA."""
    pcas = pc_per_arcsec(t)
    rmax = (rmax_pc or 0.6 * t.R_pc) / pcas
    sm = gaussian_filter(m0, 0.5 * t.beam_arcsec / FWHM / (xs[1] - xs[0]))
    X, Y = np.meshgrid(xs, xs)
    r = np.hypot(X, Y)
    inside = (r < rmax) & (sm > 0)
    # a hole is a minimum surrounded by emission: compare to the median of its ring
    k = np.argmin(np.where(inside, sm, np.inf))
    jy, jx = np.unravel_index(k, sm.shape)
    off = np.hypot(X[jy, jx], Y[jy, jx]) * pcas
    pa = np.degrees(np.arctan2(X[jy, jx], Y[jy, jx])) % 360
    return dict(offset_pc=float(off), pa=float(pa), depth=float(sm[jy, jx] / np.nanmax(sm)))


def deproject_a1(m0, xs, t: Target, nr=8):
    """A1/A0 of moment 0 in annuli of the deprojected plane, about the AGN."""
    pcas = pc_per_arcsec(t)
    X, Y = np.meshgrid(xs, xs)
    th = np.radians(t.pa)
    maj = X * np.sin(th) + Y * np.cos(th)
    mnr = (X * np.cos(th) - Y * np.sin(th)) / np.cos(np.radians(min(t.inc, 85.0)))
    R = np.hypot(maj, mnr) * pcas
    ph = np.arctan2(mnr, maj)
    edges = np.linspace(0, 1.2 * t.R_pc, nr + 1)
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (R >= a) & (R < b)
        w = m0[sel]
        a0 = w.sum()
        a1 = np.hypot((w * np.cos(ph[sel])).sum(), (w * np.sin(ph[sel])).sum())
        out.append((0.5 * (a + b), a1 / a0 if a0 > 0 else np.nan))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dump")
    ap.add_argument("--target", default="ngc613", choices=TARGETS)
    ap.add_argument("--phi", type=float, default=None,
                    help="observer azimuth [deg]; default: 90 deg from the m = 1 phase")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    t = TARGETS[a.target]
    p = load(a.dump)
    phase = m1_phase(p["x"], p["m"])
    phi = phase + 90.0 if a.phi is None else a.phi
    cb = cube(p, t, phi)
    m0, m1, m2 = moments(cb, t)
    res = dict(target=t.name, time=p["time"], phase=phase, phi=phi, hole=hole(m0, cb["xs"], t),
               a1=deproject_a1(m0, cb["xs"], t).tolist(), vunit=cb["vunit"])
    out = Path(a.out or Path(a.dump).parent / "mockobs")
    out.mkdir(parents=True, exist_ok=True)
    stem = f"{Path(a.dump).name}_{a.target}"
    np.savez_compressed(out / f"{stem}.npz", m0=m0, m1=m1, m2=m2, xs=cb["xs"], vs=cb["vs"],
                        pv=pv(cb, t), spec=cb["cube"].sum((1, 2)), target=json.dumps(asdict(t)))
    (out / f"{stem}.json").write_text(json.dumps(res, indent=1))
    print(json.dumps({k: res[k] for k in ("target", "time", "phase", "phi", "hole")}, indent=1))


if __name__ == "__main__":
    main()
