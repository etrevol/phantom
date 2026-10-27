"""Readers shared by every project: phantom dumps, .ev files, .setup and .in files.

    read_dump(path)            -> header dict and particle arrays, double precision as written
    read_ev(path)              -> dict of columns from a phantom .ev or any '# name ...' table
    read_ev_series(dir, pre)   -> all <pre>NN.ev of a run joined in order
    read_setup(path)           -> key = value pairs of a .setup or .in file
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

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
