"""Readers for the PPI torus runs: phantom dumps, .ev files, and the live-analysis output.

Not run directly. Imported by the vis scripts, which find it beside themselves.

    read_dump(path)      -> header dict and particle arrays, double precision as written
    read_ev(path)        -> dict of columns from a phantom .ev or a <prefix>_ppi.ev file
    read_maps(path)      -> list of frames from <prefix>_ppimaps.bin
    run_info(run_dir)    -> setup parameters, model and file names of one run
"""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
for _p in (_HERE, _HERE.parent / "theory"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

# --------------------------------------------------------------------------- dumps

# phantom's ndatatypes order; "real" is the default kind, 8 bytes with DOUBLEPRECISION
_DTYPES = ["<i4", "<i1", "<i2", "<i4", "<i8", "<f8", "<f4", "<f8"]


def _records(buf):
    """Fortran sequential records, 4-byte markers."""
    pos, n = 0, len(buf)
    while pos < n:
        ln = int.from_bytes(buf[pos:pos + 4], "little", signed=True)
        yield buf[pos + 4:pos + 4 + ln]
        pos += ln + 8


def read_dump(path):
    """Header and arrays of a full phantom dump (single MPI block layout)."""
    buf = Path(path).read_bytes()
    rec = _records(buf)
    first = next(rec)
    if np.frombuffer(first[:4], "<i4")[0] not in (60769, 690706):
        raise ValueError(f"{path}: not a phantom dump")
    fileid = next(rec).decode("ascii", "replace").strip()
    header = {}
    for dt in _DTYPES:
        n = int(np.frombuffer(next(rec), "<i4")[0])
        if n <= 0:
            continue
        tags = next(rec)
        vals = np.frombuffer(next(rec), dt)
        for k in range(n):
            tag = tags[16 * k:16 * (k + 1)].decode().strip()
            header.setdefault(tag, vals[k].item())
    nblocks = int(np.frombuffer(next(rec), "<i4")[0])
    blocks = []
    for _ in range(nblocks):
        r = next(rec)
        number = int(np.frombuffer(r[:8], "<i8")[0])
        nums = np.frombuffer(r[8:8 + 32], "<i4")
        blocks.append((number, nums))
    arrays = []
    for number, nums in blocks:
        out = {}
        for t, dt in enumerate(_DTYPES):
            for _ in range(int(nums[t])):
                tag = next(rec).decode().strip()
                out[tag] = np.frombuffer(next(rec), dt)
        arrays.append(out)
    return dict(fileid=fileid, header=header, gas=arrays[0],
                sinks=arrays[1] if len(arrays) > 1 else {})


# --------------------------------------------------------------------------- ev files

def read_ev(path):
    """Columns of a phantom .ev ('[01 time]' header) or of <prefix>_ppi.ev ('# time ...')."""
    path = Path(path)
    with path.open() as fh:
        head = fh.readline()
    names = re.findall(r"\[\s*\d+\s+([^\]]+?)\s*\]", head)
    if not names:
        names = head.lstrip("#").split()
    data = np.loadtxt(path, ndmin=2)
    if data.size == 0:
        return {n: np.array([]) for n in names}
    t = data[:, 0]
    # a restart appends from the dump time: keep the last occurrence of each time
    _, last = np.unique(t[::-1], return_index=True)
    keep = np.sort(len(t) - 1 - last)
    data = data[keep]
    return {n: data[:, i] for i, n in enumerate(names[:data.shape[1]])}


def read_ev_series(run_dir, prefix="ppitorus"):
    """All <prefix>NN.ev of a run, joined in order."""
    files = sorted(Path(run_dir).glob(f"{prefix}[0-9][0-9].ev"))
    parts = [read_ev(f) for f in files]
    if not parts:
        return {}
    keys = parts[0].keys()
    return {k: np.concatenate([p[k] for p in parts if k in p]) for k in keys}


# --------------------------------------------------------------------------- maps

_MAPS_MAGIC = 1984208
_HDR = np.dtype([("magic", "<i4"), ("version", "<i4"), ("time", "<f8"),
                 ("nr", "<i4"), ("nphi", "<i4"), ("nx", "<i4"), ("nz", "<i4"),
                 ("r_lo", "<f8"), ("r_hi", "<f8"), ("x_half", "<f8"), ("z_half", "<f8")])


def read_maps(path, every=1, tmin=None, tmax=None):
    """Frames of <prefix>_ppimaps.bin as dicts; arrays are views into one buffer.

    polar fields have shape (nphi, nr): rho, p, v_R, v_phi on the mid-plane;
    'slice' is rho(z, x) in y = 0; 'avg' is the azimuthal mean of rho(z, R).
    """
    buf = Path(path).read_bytes()
    frames, pos, k = [], 0, 0
    while pos + _HDR.itemsize <= len(buf):
        h = np.frombuffer(buf, _HDR, count=1, offset=pos)[0]
        if h["magic"] != _MAPS_MAGIC:
            raise ValueError(f"{path}: bad frame at byte {pos}")
        nr, nphi, nx, nz = int(h["nr"]), int(h["nphi"]), int(h["nx"]), int(h["nz"])
        pos += _HDR.itemsize
        npol, nsl, nav = 4 * nr * nphi, nx * nz, nr * nz
        if pos + 4 * (npol + nsl + nav) > len(buf):
            break                                      # frame still being written
        t = float(h["time"])
        if (k % every == 0 and (tmin is None or t >= tmin)
                and (tmax is None or t <= tmax)):
            pol = np.frombuffer(buf, "<f4", npol, pos).reshape(4, nphi, nr)
            sl = np.frombuffer(buf, "<f4", nsl, pos + 4 * npol).reshape(nz, nx)
            av = np.frombuffer(buf, "<f4", nav, pos + 4 * (npol + nsl)).reshape(nz, nr)
            r_lo, r_hi = float(h["r_lo"]), float(h["r_hi"])
            x_half, z_half = float(h["x_half"]), float(h["z_half"])
            frames.append(dict(
                time=t, rho=pol[0], p=pol[1], v_r=pol[2], v_phi=pol[3], slice=sl, avg=av,
                r=r_lo + (np.arange(nr) + 0.5) * (r_hi - r_lo) / nr,
                phi=(np.arange(nphi) + 0.5) * 2.0 * math.pi / nphi,
                x=-x_half + (np.arange(nx) + 0.5) * 2.0 * x_half / nx,
                z=-z_half + (np.arange(nz) + 0.5) * 2.0 * z_half / nz,
                r_edges=np.linspace(r_lo, r_hi, nr + 1),
                phi_edges=np.linspace(0.0, 2.0 * math.pi, nphi + 1),
                x_edges=np.linspace(-x_half, x_half, nx + 1),
                z_edges=np.linspace(-z_half, z_half, nz + 1)))
        pos += 4 * (npol + nsl + nav)
        k += 1
    # restarts append frames from the dump time: keep the last occurrence of each time
    seen, out = {}, []
    for f in frames:
        seen[f["time"]] = f
    for t in sorted(seen):
        out.append(seen[t])
    return out


# --------------------------------------------------------------------------- runs

def read_setup(path):
    """key = value pairs of a phantom .setup or .in file."""
    out = {}
    for line in Path(path).read_text().splitlines():
        line = line.split("!")[0].strip()
        if "=" not in line or line.startswith("#"):
            continue
        k, v = (s.strip() for s in line.split("=", 1))
        try:
            out[k] = int(v) if re.fullmatch(r"[+-]?\d+", v) else float(v)
        except ValueError:
            out[k] = v
    return out


def find_run_dir(start="."):
    """The run directory: one holding a .setup, at start or above it."""
    p = Path(start).resolve()
    for q in (p, *p.parents[:3]):
        if list(q.glob("*.setup")):
            return q
    raise FileNotFoundError(f"no .setup at or above {p}")


def run_info(run_dir="."):
    """Parameters, model and files of one run."""
    import torus_model as tm

    run_dir = find_run_dir(run_dir)
    setup = sorted(run_dir.glob("*.setup"))[0]
    prefix = setup.stem
    par = read_setup(setup)
    model = tm.TorusModel(r_center=par.get("r_center", 1.0),
                          C_prime=par.get("C_prime", 0.2),
                          gamma=par.get("gamma", 1.3), M_bh=par.get("M_bh", 4.5e7),
                          chi=par.get("chi", 500.0), rho_0=par.get("rho_0", 1e-13))
    return dict(dir=run_dir, prefix=prefix, par=par, model=model,
                modes=run_dir / f"{prefix}_ppi.ev",
                maps=run_dir / f"{prefix}_ppimaps.bin",
                dumps=sorted(run_dir.glob(f"{prefix}_[0-9][0-9][0-9][0-9][0-9]")))


def info_lines(info):
    """Model parameters for figure annotations, one item per entry."""
    par, m = info["par"], info["model"]
    seed = par.get("pert_kind", "none")
    if seed != "none":
        seed = f"{seed}, A = {par.get('pert_amp', 0):g}" + (
            f", m = {par.get('pert_m')}" if seed in ("single", "multi") else "")
    return [f"C' = {m.C_prime:g},  gamma = {m.gamma:.4g}",
            f"N = {par.get('np', 0):.3g},  {par.get('lattice', '?')}, "
            f"nsym = {par.get('nsym', 1)}",
            f"seed: {seed}",
            f"P_orb = {m.p_orb * m.t_unit / 3.15576e7:.3g} yr"]
