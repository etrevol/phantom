#!/usr/bin/env python3
"""Read the r_c(t) curves of the paper's Fig. A.1 off the PDF, by colour.

Renders page 13 at 400 dpi, crops the figure and, for each visible run, takes the median
row of its colour in every pixel column. Prints the times at which r_c first crosses
0.005, 0.01, 0.015, 0.02 and the plateau after t = 200. The axis calibration (frame
columns, tick rows) was read off this crop by eye and holds for this PDF only.

    digitise_figA1.py PAPER.pdf
"""

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image

X0, X1, Y0, Y25 = 188, 1238, 801, 102        # t = 0.01 and 1000 columns, r_c = 0, 0.025 rows
RUNS = {"64k": (141, 181, 61), "128k, eps 1e-4": (139, 128, 190),
        "256k": (232, 166, 48), "sym": (199, 110, 36)}


def main(pdf):
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(["pdftoppm", "-f", "13", "-l", "13", "-r", "400", "-png", pdf,
                        f"{tmp}/p"], check=True)
        im = Image.open(next(Path(tmp).glob("p*.png"))).convert("RGB")
    w, h = im.size
    a = np.asarray(im.crop((int(0.08 * w), int(0.40 * h), int(0.49 * w), int(0.60 * h))),
                   dtype=float)
    for name, rgb in RUNS.items():
        m = np.sqrt(((a - np.array(rgb)) ** 2).sum(2)) < 45
        m[:, :X0 + 5] = m[:, X1 - 3:] = False
        m[:90] = m[Y0 + 8:] = False
        m[150:640, 280:620] = False                 # legend
        t, r = [], []
        for col in range(X0 + 5, X1 - 3):
            rows = np.nonzero(m[:, col])[0]
            if rows.size:
                t.append(10 ** (-2 + 5 * (col - X0) / (X1 - X0)))
                r.append(0.025 * (Y0 - np.median(rows)) / (Y0 - Y25))
        t, r = np.array(t), np.array(r)
        cross = [t[(r >= lev) & (t > 5)][0] for lev in (0.005, 0.01, 0.015, 0.02)]
        print(f"{name:15s} r_c 0.005/0.01/0.015/0.02 at t = "
              + " ".join(f"{c:6.1f}" for c in cross) + f"   plateau {r[t > 200].mean():.4f}")


if __name__ == "__main__":
    main(sys.argv[1])
