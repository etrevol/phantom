"""A minimal reader of a primary-HDU FITS image (no astropy here): header and data.

    hdr, data = read_fits(path)   # data in numpy order (NAXIS4, NAXIS3, NAXIS2, NAXIS1)
"""
from __future__ import annotations

import numpy as np

_BITPIX = {8: ">u1", 16: ">i2", 32: ">i4", 64: ">i8", -32: ">f4", -64: ">f8"}


def read_fits(path):
    hdr = {}
    with open(path, "rb") as fh:
        nblock = 0
        done = False
        while not done:
            block = fh.read(2880)
            nblock += 1
            for i in range(36):
                card = block[80 * i:80 * (i + 1)].decode("ascii", "replace")
                key = card[:8].strip()
                if key == "END":
                    done = True
                    break
                if card[8:10] != "= ":
                    continue
                val = card[10:].split("/")[0].strip() if not card[10:].strip().startswith("'") \
                    else card[10:].split("'")[1].strip()
                try:
                    val = int(val)
                except ValueError:
                    try:
                        val = float(val)
                    except ValueError:
                        val = {"T": True, "F": False}.get(val, val)
                hdr[key] = val
        shape = [hdr[f"NAXIS{k}"] for k in range(hdr["NAXIS"], 0, -1)]
        data = np.frombuffer(fh.read(int(np.prod(shape)) * abs(hdr["BITPIX"]) // 8),
                             dtype=_BITPIX[hdr["BITPIX"]]).reshape(shape)
    data = data.astype(np.float64)
    if "BSCALE" in hdr or "BZERO" in hdr:
        data = data * hdr.get("BSCALE", 1.0) + hdr.get("BZERO", 0.0)
    return hdr, data
