#!/usr/bin/env python3
"""Linear Papaloizou-Pringle modes of the constant-l torus, in 3D and in the 2D cylinder.

PP84 eq. (4.1). With perturbations ~ exp[i(m phi + sigma t)] and W = p'/[rho (sigma + m Omega)],

    (1/R) d/dR (rho R dW/dR) + d/dz (rho dW/dz) - (m^2/R^2) rho W
        = - (sigma + m Omega)^2 rho^2 W / (gamma p)

a quadratic eigenvalue problem sigma^2 A W + sigma B W + C W = 0 with A = rho^2/(gamma p),
B = 2 m Omega A, C = L + m^2 Omega^2 A, solved here as the 2N generalised problem
sigma N x = M x, x = [sigma W; W]. A mode grows when Im(sigma) < 0, at rate -Im(sigma);
its pattern speed is -Re(sigma)/m. dim = 2 drops d/dz and uses the mid-plane profile:
that is the z-independent cylinder of the athena-collab runs, and the same discretisation
as vvlib/ppi.py on branch athena-vv.

Cells are cell-centred on a uniform (R, z) box around the upper half of the torus; a cell
is active where rho > 0 at its centre. Fluxes through faces into inactive cells are
dropped: rho -> 0 at the surface, so that is the regularity condition PP84 (4.25) with
an O(rho_face) error. The mid-plane face carries dW/dz = 0 (even modes, the fastest in
PP84) or W = 0 (odd modes).

Units: G M = 1, R_0 = 1, so rates are in Omega_0 and e-folds per orbit = 2 pi rate.
"""

from __future__ import annotations

import argparse
import math
import sys
import time
from pathlib import Path

import numpy as np
import scipy.linalg
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, str(Path(__file__).resolve().parent))
import torus_model as tm  # noqa: E402


def _operators(model, m, nr, nz, dim, parity):
    r_in, r_out = model.r_in, model.r_out
    dr = (r_out - r_in) / nr
    R = r_in + (np.arange(nr) + 0.5) * dr
    Rf = r_in + np.arange(nr + 1) * dr
    if dim == 2:
        nz, dz = 1, 1.0
        z = np.zeros(1)
    else:
        dz = model.z_top / nz
        z = (np.arange(nz) + 0.5) * dz
    zf = np.arange(nz + 1) * dz

    rr, zz = np.meshgrid(R, z, indexing="ij")
    rho = model.rho(rr, zz)
    act = rho > 0.0
    idx = -np.ones((nr, nz), dtype=int)
    idx[act] = np.arange(act.sum())
    n = int(act.sum())

    pres = rho * model.p_over_rho(rr, zz)
    omega = np.sqrt(model.r_center) / rr**2           # l = sqrt(G M R_0)
    a_diag = np.where(act, rho**2 / (model.gamma * np.where(act, pres, 1.0)), 0.0)

    rows, cols, vals = [], [], []

    def put(r, c, v):
        rows.append(np.asarray(r)), cols.append(np.asarray(c)), vals.append(np.asarray(v))

    ii, jj = np.nonzero(act)
    k = idx[ii, jj]
    put(k, k, -(m**2 / R[ii] ** 2) * rho[ii, jj])

    # radial faces
    for s in (-1, 1):
        fi = ii + (s > 0)                                   # face index in Rf
        rho_f = model.rho(Rf[fi], z[jj]) if dim == 3 else model.rho(Rf[fi], 0.0)
        coef = rho_f * Rf[fi] / (dr * dr * R[ii])
        nb = ii + s
        ok = (nb >= 0) & (nb < nr)
        ok[ok] &= act[nb[ok], jj[ok]]
        put(k, k, -np.where(ok, coef, 0.0))
        put(k[ok], idx[nb[ok], jj[ok]], coef[ok])

    if dim == 3:
        for s in (-1, 1):
            fj = jj + (s > 0)
            coef = model.rho(R[ii], zf[fj]) / (dz * dz)
            nb = jj + s
            ok = (nb >= 0) & (nb < nz)
            ok[ok] &= act[ii[ok], nb[ok]]
            put(k, k, -np.where(ok, coef, 0.0))
            put(k[ok], idx[ii[ok], nb[ok]], coef[ok])
            if s < 0 and parity == "odd":
                # W_{-1} = -W_0 across the mid-plane
                bottom = jj == 0
                put(k[bottom], k[bottom], -2.0 * coef[bottom])

    lop = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                        shape=(n, n))

    a = a_diag[act]
    om = omega[act]
    return dict(lop=lop, a=a, om=om, n=n, R=R, z=z, act=act, idx=idx)


def _pencil(op, m):
    """Sparse M, N of sigma N x = M x."""
    n, a, om = op["n"], op["a"], op["om"]
    c_mat = op["lop"] + sp.diags(m * m * om * om * a)
    ident = sp.identity(n, format="csr")
    big_m = sp.bmat([[sp.diags(-2.0 * m * om * a), -c_mat], [ident, None]], format="csc")
    big_n = sp.bmat([[sp.diags(a), None], [None, ident]], format="csc")
    return big_m, big_n


def spectrum(model, m, nr=60, nz=30, dim=3, parity="even", vectors=False):
    """All complex eigenfrequencies sigma (and eigenvectors W if asked), dense."""
    op = _operators(model, m, nr, nz, dim, parity)
    big_m, big_n = (x.toarray() for x in _pencil(op, m))
    n = op["n"]
    if vectors:
        sig, vec = scipy.linalg.eig(big_m, big_n)
        ok = np.isfinite(sig)
        return sig[ok], vec[n:, ok], op
    sig = scipy.linalg.eig(big_m, big_n, right=False)
    return sig[np.isfinite(sig)]


def near(model, m, sigma0, nr, nz, dim=3, parity="even", k=8):
    """The k eigenfrequencies nearest sigma0, by sparse shift-invert."""
    op = _operators(model, m, nr, nz, dim, parity)
    big_m, big_n = _pencil(op, m)
    # (M - s N)^-1 N x = x / (sigma - s): the largest of these are the nearest sigma
    lu = spla.splu((big_m - sigma0 * big_n).tocsc())
    opinv = spla.LinearOperator(big_m.shape, matvec=lambda x: lu.solve(big_n @ x),
                                dtype=complex)
    mu = spla.eigs(opinv, k=k, which="LM", return_eigenvectors=False)
    return sigma0 + 1.0 / mu


DENSE_CELLS = 1000      # active cells up to which the full spectrum is affordable


def _pick(sig, m, min_rate=1e-6):
    grow = -sig.imag
    if not (grow > min_rate).any():
        return None
    return sig[np.argmax(grow)]


def ladder(model, m, nr=60, nz=30, dim=3, parity="even"):
    """Fastest mode on a sequence of grids ending at (nr, nz), finest last.

    The full spectrum is taken on the coarsest grid; each finer grid only follows the
    fastest mode, by shift-invert around the previous estimate. Returns a list of
    (nr, nz, sigma) with sigma None where no mode grows.
    """
    grids = [(nr, nz)]
    while True:
        gr, gz = grids[-1]
        if gr // 2 < 8:
            break
        n_act = _operators(model, m, gr, gz, dim, parity)["n"]
        if n_act <= DENSE_CELLS:
            break
        grids.append((gr // 2, max(gz // 2, 1)))
    grids.reverse()
    out, sig = [], None
    for k, (gr, gz) in enumerate(grids):
        if k == 0:
            sig = _pick(spectrum(model, m, gr, gz, dim, parity), m)
        elif sig is not None:
            sig = _pick(near(model, m, sig, gr, gz, dim, parity), m)
        out.append((gr, gz, sig))
    return out


def fastest(model, m, nr=60, nz=30, dim=3, parity="even"):
    """(rate, pattern speed) of the fastest growing mode, or (0, nan)."""
    sig = ladder(model, m, nr, nz, dim, parity)[-1][2]
    if sig is None:
        return 0.0, float("nan")
    return float(-sig.imag), float(-sig.real / m)


def converged(model, m, nr, nz, dim=3, parity="even"):
    """Fastest mode at (nr, nz), with the relative change from the grid half as fine."""
    steps = ladder(model, m, nr, nz, dim, parity)
    sig = steps[-1][2]
    if len(steps) == 1:
        coarse = _pick(spectrum(model, m, nr // 2, max(nz // 2, 1), dim, parity), m)
    else:
        coarse = steps[-2][2]
    rate = -sig.imag if sig is not None else 0.0
    rate_c = -coarse.imag if coarse is not None else 0.0
    rel = abs(rate - rate_c) / rate if rate > 0 else float("nan")
    patt = -sig.real / m if sig is not None else float("nan")
    return dict(m=m, rate=rate, pattern=patt, rate_coarse=rate_c, rel_change=rel,
                efolds_per_orbit=2.0 * math.pi * rate / model.omega_0)


def check_pp84(nr=60, nz=30):
    """Fastest even modes against PP84 Table 1 (growth time 1/rate in 1/Omega_0)."""
    rows = []
    for (preset, m), (ta, tb) in tm.PP84_TABLE1.items():
        model = tm.model_from_preset(preset)
        t0 = time.time()
        res = converged(model, m, nr, nz)
        tau = 1.0 / res["rate"] if res["rate"] > 0 else float("inf")
        rows.append((preset, m, tau, ta, tb, res["rel_change"], time.time() - t0))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--preset", default="athena", choices=sorted(tm.PRESETS))
    ap.add_argument("--m", default="1,2,3,4", help="comma-separated modes")
    ap.add_argument("--dim", type=int, default=3, choices=(2, 3))
    ap.add_argument("--nr", type=int, default=60)
    ap.add_argument("--nz", type=int, default=30)
    ap.add_argument("--parity", default="even", choices=("even", "odd"))
    ap.add_argument("--pp84", action="store_true", help="compare with PP84 Table 1")
    a = ap.parse_args(argv)

    if a.pp84:
        print(f"{'model':8s} {'m':>2s} {'tau':>7s} {'PP84 A':>7s} {'PP84 B':>7s} "
              f"{'d(rate)':>8s} {'sec':>5s}")
        for p, m, tau, ta, tb, rel, sec in check_pp84(a.nr, a.nz):
            print(f"{p:8s} {m:2d} {tau:7.2f} {ta:7.1f} {tb:7.1f} {rel:8.3f} {sec:5.0f}")
        return 0

    model = tm.model_from_preset(a.preset)
    print(f"{a.preset}: C' = {model.C_prime}, gamma = {model.gamma:.4g}, dim = {a.dim}, "
          f"grid {a.nr} x {a.nz if a.dim == 3 else 1}, {a.parity}")
    print(f"{'m':>2s} {'rate/Om0':>9s} {'e-fold/orb':>10s} {'pattern':>8s} {'d(rate)':>8s}")
    for m in [int(s) for s in a.m.split(",")]:
        r = converged(model, m, a.nr, a.nz, a.dim, a.parity)
        print(f"{m:2d} {r['rate']:9.4f} {r['efolds_per_orbit']:10.3f} {r['pattern']:8.4f} "
              f"{r['rel_change']:8.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
