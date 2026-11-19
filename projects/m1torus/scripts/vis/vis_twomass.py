#!/usr/bin/env python3
"""A two-mass torus (setup f_heavy, mu_heavy): light (star type) and heavy (dark-matter
type) particles over time - the m = 1 amplitude and phase of each species and of the
mass-weighted total, with the noise 2/sqrt(N) of each, and the mass segregation of
Foote et al. (2020): mean semi-major axis, eccentricity, inclination and |z| of each.

    vis_twomass.py --run DIR [--ref RUN] [--every 1]

--ref overplots the mass-weighted A1 of a one-mass run (A14: same draw and softening).
Writes figs/twomass.png and twomass.npz in the run directory.
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[2] / "common"))
import labels as L  # noqa: E402
import phantom_io as pio  # noqa: E402

LIGHT, HEAVY = 4, 5
P_ORB = 2.0 * np.pi


def mode(x, w):
    z = (w * np.exp(1j * np.arctan2(x[:, 1], x[:, 0]))).sum() / w.sum()
    return 2 * abs(z), np.degrees(np.angle(z)) % 360


def row(path):
    d = pio.read_dump(path)
    g, s = d["gas"], d["sinks"]
    it = g["itype"]
    mt = d["header_lists"]["massoftype"]
    m = np.array([mt[k - 1] for k in it])
    live = g["h"] > 0
    # the mode about the system barycentre (the analysis convention), elements about the sink;
    # the barycentre is measured, since a two-mass run need not keep it at the origin
    xb = np.r_[s["x"][0], s["y"][0], s["z"][0]]
    xs = ((m[live, None] * np.c_[g["x"], g["y"], g["z"]][live]).sum(0) + s["m"][0] * xb) / (m[live].sum() + s["m"][0])
    x = np.c_[g["x"], g["y"], g["z"]] - xs
    xb = xb - xs
    vb = np.r_[s["vx"][0], s["vy"][0], s["vz"][0]]
    r = x - xb
    v = np.c_[g["vx"], g["vy"], g["vz"]] - vb
    M = s["m"][0]
    rr = np.linalg.norm(r, axis=1)
    v2 = (v ** 2).sum(1)
    a = 1.0 / (2.0 / rr - v2 / M)
    ev = ((v2 - M / rr)[:, None] * r - (r * v).sum(1)[:, None] * v) / M
    h = np.cross(r, v)
    inc = np.degrees(np.arccos(np.clip(h[:, 2] / np.linalg.norm(h, axis=1), -1, 1)))
    out = dict(t=d["header"]["time"] / P_ORB, rbh=float(np.hypot(*xb[:2])), drift=float(np.linalg.norm(xs)))
    out["all"] = mode(x[live], m[live])
    for name, k in (("light", LIGHT), ("heavy", HEAVY)):
        sel = live & (it == k)
        b = sel & (a > 0)
        out[name] = (*mode(x[sel], np.ones(sel.sum())), a[b].mean(), np.linalg.norm(ev[b], axis=1).mean(),
                     inc[b].mean(), np.abs(r[sel, 2]).mean(), 2 / np.sqrt(sel.sum()))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--ref", default=None)
    ap.add_argument("--every", type=int, default=1)
    a = ap.parse_args()
    run = Path(a.run).resolve()
    rows = [row(p) for p in sorted(run.glob("m1torus_[0-9][0-9][0-9][0-9][0-9]"))[::a.every]]
    t = np.array([q["t"] for q in rows])
    li = np.array([q["light"] for q in rows])
    he = np.array([q["heavy"] for q in rows])
    al = np.array([q["all"] for q in rows])
    rbh = np.array([q["rbh"] for q in rows])
    drift = np.array([q["drift"] for q in rows])
    np.savez(run / "twomass.npz", t=t, light=li, heavy=he, all=al, rbh=rbh, drift=drift,
             cols="A1 phase a e inc |z| noise")

    fig, ax = plt.subplots(3, 2, figsize=(11, 9.5), sharex=True)
    ax = ax.ravel()
    cl, ch = "C0", "C3"
    ax[0].plot(t, li[:, 0], color=cl, label="light")
    ax[0].plot(t, he[:, 0], color=ch, label=r"heavy ($\times 10$ mass)")
    ax[0].plot(t, al[:, 0], color="k", lw=1.2, label="all, mass-weighted")
    if a.ref:
        m1 = pio.read_ev(Path(a.ref) / "m1torus_m1.ev")
        ax[0].plot(m1["time"] / P_ORB, m1["A1"] / m1["sigma0"], color="0.6", lw=0.9, ls="--",
                   label="one mass, same draw")
    ax[0].fill_between(t, 0, he[:, 6], color=ch, alpha=0.12, lw=0)
    ax[0].fill_between(t, 0, li[:, 6], color=cl, alpha=0.2, lw=0)
    ax[0].set(ylabel=r"$A_1/\Sigma_0$", ylim=(0, None))
    ax[0].legend(fontsize=8, loc="upper left")
    dph = (he[:, 1] - li[:, 1] + 180) % 360 - 180
    ax[1].plot(t, dph, color="k")
    ax[1].axhline(0, color="0.6", lw=0.6)
    ax[1].set(ylabel=r"$\Phi_1^{\rm heavy} - \Phi_1^{\rm light}$ [deg]", ylim=(-180, 180))
    for k, (j, lab) in enumerate(((2, r"mean $a$ $[R_{\rm tor}]$"), (3, r"mean $e$"),
                                  (4, r"mean inclination [deg]"), (5, r"mean $|z|$ $[R_{\rm tor}]$"))):
        axi = ax[2 + k]
        axi.plot(t, li[:, j], color=cl)
        axi.plot(t, he[:, j], color=ch)
        axi.set(ylabel=lab)
    for axi in ax[4:]:
        axi.set_xlabel(r"$t / P_{\rm orb}$")
    name, pars, suffix = L.model_parts(run)
    fig.suptitle("Light and Heavy Particles of the Torus\n" + L.wrapped(f"{name} ({', '.join(pars.values())}){suffix}", 110),
                 fontsize=9.5)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    (run / "figs").mkdir(exist_ok=True)
    fig.savefig(run / "figs" / "twomass.png", dpi=120)
    print(run / "figs" / "twomass.png")


if __name__ == "__main__":
    main()
