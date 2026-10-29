#!/usr/bin/env python3
"""Time series of a B run (colliding clouds or gas): how the torus dissipates.

Panels: the m = 1 mode (A1/Sigma0, |r_c|) from the live analysis; radiated energy
E_diss/|E0| and the mass fraction in the sink; thickness (median |z|/R) and
sigma_z/v_phi; the 10/50/90% radii of the mass; and, from synthetic NGC 613 maps of
every dump (mockobs), the offset of the moment-0 hole from the AGN against the
observed 3 pc. Also writes the series as JSON beside the figure.

    vis_b.py --run DIR [--every 1]
"""
import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "common"))
import phantom_io as pio  # noqa: E402
import mockobs as mo  # noqa: E402

P_ORB = 2 * np.pi


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--every", type=int, default=1)
    a = ap.parse_args()
    run = Path(a.run)
    m1 = pio.read_ev(run / "m1torus_m1.ev")
    ev = pio.read_ev_series(run, "m1torus")
    e0 = abs(ev["etot"][0])
    dumps = sorted(run.glob("m1torus_[0-9][0-9][0-9][0-9][0-9]"))[::a.every]
    rows = []
    for d in dumps:
        dd = pio.read_dump(d)
        p = mo.load(d)
        x, v, m = p["x"], p["v"], p["m"]
        R = np.hypot(x[:, 0], x[:, 1])
        vphi = (x[:, 0] * v[:, 1] - x[:, 1] * v[:, 0]) / np.maximum(R, 1e-9)
        q = np.quantile(R, [0.1, 0.5, 0.9])
        msink = dd["sinks"]["m"][0]
        cb = mo.cube(p, mo.NGC613, mo.m1_phase(x, m) + 90.0)
        m0, _, _ = mo.moments(cb, mo.NGC613)
        h = mo.hole(m0, cb["xs"], mo.NGC613)
        rows.append(dict(t=p["time"] / P_ORB, ediss=dd["header"].get("ediss", 0.0) / e0,
                         macc=(msink - 1.0) / (msink - 1.0 + m.sum()) if m.sum() > 0 else np.nan,
                         zr=float(np.median(np.abs(x[:, 2]) / np.maximum(R, 1e-9))),
                         sz=float(np.std(v[:, 2]) / np.mean(np.abs(vphi))),
                         r10=q[0], r50=q[1], r90=q[2], hole_pc=h["offset_pc"]))
    s = {k: np.array([r[k] for r in rows]) for k in rows[0]}
    (run / "figs").mkdir(exist_ok=True)
    (run / "figs" / "b_series.json").write_text(json.dumps({k: v.tolist() for k, v in s.items()}))

    fig, axs = plt.subplots(3, 2, figsize=(12, 11), sharex=True)
    t1 = m1["time"] / P_ORB
    ax = axs[0, 0]
    ax.plot(t1, m1["A1"] / m1["sigma0"], "C0")
    ax.axhline(np.sqrt(np.pi / m1["n"][0]), color="C0", ls="--", lw=0.8)
    ax.set_yscale("log")
    ax.set_ylabel("A1/Sigma0 (dashed: noise)")
    ax2 = ax.twinx()
    ax2.plot(t1, np.hypot(m1["xc"], m1["yc"]), "C3", lw=0.8)
    ax2.set_ylabel("|r_c|", color="C3")
    ax.set_title("m = 1 mode and black-hole offset", fontsize=10)
    ax = axs[0, 1]
    ax.plot(s["t"], s["ediss"], "k.-", label="E_diss/|E0|")
    ax.plot(s["t"], s["macc"], "C1.-", label="mass fraction in the sink")
    ax.legend(fontsize=8)
    ax.set_title("dissipation and inflow", fontsize=10)
    ax = axs[1, 0]
    ax.plot(s["t"], s["zr"], "k.-", label="median |z|/R")
    ax.plot(s["t"], s["sz"], "C2.-", label="sigma_z / <v_phi>")
    ax.legend(fontsize=8)
    ax.set_title("thickness", fontsize=10)
    ax = axs[1, 1]
    ax.fill_between(s["t"], s["r10"], s["r90"], color="0.85", label="10-90%")
    ax.plot(s["t"], s["r50"], "k.-", label="median R")
    ax.legend(fontsize=8)
    ax.set_ylabel("R [R_tor]")
    ax.set_title("radial distribution", fontsize=10)
    ax = axs[2, 0]
    ax.plot(s["t"], s["hole_pc"], "k.-")
    ax.axhline(3.0, color="r", ls="--", label="NGC 613: 3 pc")
    ax.set_ylabel("hole offset [pc]")
    ax.legend(fontsize=8)
    ax.set_title("synthetic NGC 613: moment-0 hole from the AGN", fontsize=10)
    axs[2, 1].axis("off")
    for ax in axs[-1, :1]:
        ax.set_xlabel("t [orbits]")
    axs[1, 1].set_xlabel("t [orbits]")
    fig.suptitle(run.name, fontsize=11)
    fig.tight_layout()
    fig.savefig(run / "figs" / "b_series.png", dpi=110)
    print(run / "figs" / "b_series.png")


if __name__ == "__main__":
    main()
