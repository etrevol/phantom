#!/usr/bin/env python3
"""B32: the m = 1 mode of the stars (dark-matter type) and of the gas clouds (star type)
in a two-component torus, per dump: amplitude A1/S0 and phase about the black hole,
the gas dispersion, unbound share and spread (stop criteria of LOG 2026-10-30).

    twocomp.py RUN_DIR [--every 2]
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "common"))
import phantom_io as pio  # noqa: E402


def mode(x, w):
    ph = np.arctan2(x[:, 1], x[:, 0])
    z = (w * np.exp(1j * ph)).sum() / w.sum()
    return 2 * abs(z), np.degrees(np.angle(z)) % 360


def row(path):
    d = pio.read_dump(path)
    g, s = d["gas"], d["sinks"]
    it = g["itype"]
    live = g["h"] > 0
    xb = np.r_[s["x"][0], s["y"][0], s["z"][0]]
    vb = np.r_[s["vx"][0], s["vy"][0], s["vz"][0]]
    M = s["m"][0]
    x = np.c_[g["x"], g["y"], g["z"]] - xb
    v = np.c_[g["vx"], g["vy"], g["vz"]] - vb
    out = dict(t=d["header"]["time"] / (2 * np.pi))
    for name, k in (("stars", 5), ("gas", 4)):
        sel = live & (it == k)
        xs, vs = x[sel], v[sel]
        a1, ph = mode(xs, np.ones(len(xs)))
        r = np.linalg.norm(xs, axis=1)
        R = np.hypot(xs[:, 0], xs[:, 1])
        E = 0.5 * (vs ** 2).sum(1) - M / r
        kk = (R > 0.5) & (R < 1.5)
        vk = np.sqrt(M / R[kk])
        out[name] = dict(n=int(sel.sum()), A1=a1, phase=ph, unbound=float(np.mean(E > 0)),
                         R90=float(np.percentile(R, 90)), R90b=float(np.percentile(R[E < 0], 90)),
                         sz=float(np.std(vs[kk, 2] / vk)),
                         noise=2 / np.sqrt(max(sel.sum(), 1)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run")
    ap.add_argument("--every", type=int, default=2)
    a = ap.parse_args()
    dumps = sorted(Path(a.run).glob("m1torus_[0-9][0-9][0-9][0-9][0-9]"))[::a.every]
    for p in dumps:
        r = row(p)
        st, ga = r["stars"], r["gas"]
        dph = (ga["phase"] - st["phase"] + 180) % 360 - 180
        print(f"t {r['t']:6.1f} | stars A1 {st['A1']:.3f} phase {st['phase']:5.0f} sz {st['sz']:.2f} | "
              f"gas n {ga['n']:5d} A1 {ga['A1']:.3f} (noise {ga['noise']:.3f}) phase-stars {dph:+5.0f} "
              f"sz {ga['sz']:.2f} unbound {ga['unbound']:.3f} R90 {ga['R90']:.2f} R90b {ga['R90b']:.2f}")


if __name__ == "__main__":
    main()
