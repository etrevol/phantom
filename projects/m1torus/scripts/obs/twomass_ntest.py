#!/usr/bin/env python3
"""A17: the rates of mass segregation in two-mass tori of different particle number.

From twomass.npz of each run (vis_twomass.py): Delta a = <a>_light - <a>_heavy and
Delta i (inclination), their slopes after the mode saturates (A1 of all particles first
above 80% of its mean over the last 30 orbits), the spreading of the light particles
(slope of <a>_light) and the saturated A1; then the power p of rate ~ N^-p between runs
(p = 1 two-body, 0.5 resonant, 0 collective).

    twomass_ntest.py RUN [RUN ...]
"""
import sys
from pathlib import Path

import numpy as np


def rates(run):
    d = np.load(Path(run) / "twomass.npz")
    t, li, he, al = d["t"], d["light"], d["heavy"], d["all"]
    plateau = al[t >= t[-1] - 30, 0].mean()
    tsat = t[np.argmax(al[:, 0] > 0.8 * plateau)]
    k = t >= tsat
    out = dict(tsat=tsat, A1=plateau)
    for name, col in (("da", 2), ("di", 4)):
        y = li[:, col] - he[:, col]
        p, cov = np.polyfit(t[k], y[k], 1, cov=True)
        out[name] = (p[0] * 100, np.sqrt(cov[0, 0]) * 100, y[-1] - y[0])
    p, cov = np.polyfit(t[k], li[k, 2], 1, cov=True)
    out["spread"] = (p[0] * 100, np.sqrt(cov[0, 0]) * 100)
    return out


def main():
    rows = []
    for run in sys.argv[1:]:
        setup = open(Path(run) / "m1torus.setup").read()
        n = int(setup.split("np =")[1].split()[0])
        r = rates(run)
        r["N"] = n
        rows.append(r)
        print(f"N {n:6d}: mode saturated from t {r['tsat']:.0f} at A1/S0 {r['A1']:.3f} | per 100 orbits: "
              f"Delta a {r['da'][0]:+.4f} +- {r['da'][1]:.4f} (total {r['da'][2]:+.4f}), "
              f"Delta i {r['di'][0]:+.2f} +- {r['di'][1]:.2f} deg (total {r['di'][2]:+.2f}), "
              f"light spreading {r['spread'][0]:+.4f} +- {r['spread'][1]:.4f}")
    if len(rows) > 1:
        N = np.log([r["N"] for r in rows])
        for key in ("da", "di", "spread"):
            y = np.array([r[key][0] for r in rows])
            if np.all(y > 0):
                p = -np.polyfit(N, np.log(y), 1)[0]
                print(f"{key}: rate ~ N^-p with p = {p:.2f}")


if __name__ == "__main__":
    main()
