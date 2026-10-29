#!/usr/bin/env python3
"""Several m1torus runs on one figure, against Bannikova et al. (2026).

Panels: A_1/Sigma_0 on a log scale with each run's shot-noise floor sqrt(pi/N), the
central-mass offset |r_c|, k = A_2/A_1 and the energy error. Prints, per run, the onset
(A_1/Sigma_0 first above 0.1), the e-folding time of the growth, and the averages over a
saturated window, with A_3 and Q_3 also corrected for the noise bias <A^2> = A^2 + 2 s^2.

    compare_runs.py --runs DIR [DIR ...] [--labels 16k,64k] [--window 110,150] [--out F]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "common"))
import phantom_io as pio
import labels as L  # noqa: E402

P_ORB = 2.0 * np.pi
N_SEC = 1200
N_PAPER = 131072
PAPER = dict(A1=40.21, A2=8.46, A3=1.77, k=0.21, Q3=0.99, r_c=0.024, r_tb=0.24)
plt.rcParams.update({"font.size": 10, "axes.labelsize": 11, "xtick.direction": "in",
                     "ytick.direction": "in", "xtick.top": True, "ytick.right": True})


def running(y, n):
    # mean over n samples, the window shortened at the ends instead of padded with zeros
    k = np.ones(n)
    return np.convolve(y, k, mode="same") / np.convolve(np.ones_like(y), k, mode="same")


def load(run):
    run = Path(run)
    m1 = pio.read_ev(sorted(run.glob("*_m1.ev"))[0])
    ev = pio.read_ev_series(run, sorted(run.glob("*.setup"))[0].stem)
    t = m1["time"] / P_ORB
    return dict(t=t, s0=m1["sigma0"], n=int(m1["n"][0]),
                smooth=2 * max(1, int(round(1.0 / np.median(np.diff(t))))),  # 2 orbits
                a=[m1[f"A{m}"] for m in (1, 2, 3)],
                r_c=np.hypot(m1["xc"], m1["yc"]), r_tb=np.hypot(m1["xtb"], m1["ytb"]),
                te=ev["time"] / P_ORB, etot=ev["etot"])


def summary(d, t0, t1):
    t, s0 = d["t"], d["s0"]
    a1, a2, a3 = (x / s0 for x in d["a"])
    smooth = running(a1, d["smooth"])
    above = np.nonzero(smooth > 0.1)[0]
    onset = t[above[0]] if above.size else np.nan
    grow = (smooth > 0.05) & (smooth < 0.2) & (t < (onset if above.size else t[-1]) + 30)
    tau = 1.0 / np.polyfit(t[grow], np.log(smooth[grow]), 1)[0] if grow.sum() > 5 else np.nan
    w = (t >= t0) & (t <= t1)
    noise2 = 4.0 / d["n"]                    # 2 sigma^2 / S0^2 with sigma/S0 = sqrt(2/N)
    m = {k: float(np.mean(v[w])) for k, v in dict(A1=a1, A2=a2, A3=a3).items()}
    a3c = np.sqrt(max(np.mean(a3[w] ** 2) - noise2, 0.0))
    a2c = np.sqrt(max(np.mean(a2[w] ** 2) - noise2, 0.0))
    k = m["A2"] / m["A1"]
    return dict(onset=onset, tau=tau, **m, k=k, Q3=(m["A3"] / m["A1"]) / k**2,
                A3c=a3c, Q3c=(a3c / m["A1"]) / (a2c / m["A1"]) ** 2,
                r_c=float(np.mean(d["r_c"][w])), r_tb=float(np.mean(d["r_tb"][w])))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", nargs="+", required=True)
    ap.add_argument("--labels", default=None,
                    help="comma separated, one per run; default: the parameters that differ")
    ap.add_argument("--window", default="110,150", help="saturated window, orbits")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)

    runs = [load(r) for r in a.runs]
    title, auto = L.legend_labels([Path(r) for r in a.runs])
    labels = a.labels.split(",") if a.labels else auto
    ref = L.REF["bannikova26"]
    t0, t1 = (float(x) for x in a.window.split(","))
    a1_paper = PAPER["A1"] / (N_PAPER - 1) * N_SEC

    fig, ax = plt.subplots(2, 2, figsize=(12, 8.4), sharex=True)
    for i, (d, lab) in enumerate(zip(runs, labels)):
        c = f"C{i}"
        a1 = d["a"][0] / d["s0"]
        ax[0, 0].semilogy(d["t"], a1, color=c, lw=0.6, alpha=0.35)
        ax[0, 0].semilogy(d["t"], running(a1, d["smooth"]), color=c, lw=1.4, label=lab)
        ax[0, 0].axhline(np.sqrt(np.pi / d["n"]), color=c, ls="--", lw=0.9)
        ax[0, 1].plot(d["t"], d["r_c"], color=c, lw=1.0, label=lab)
        with np.errstate(divide="ignore", invalid="ignore"):
            ax[1, 0].plot(d["t"], running(d["a"][1] / d["a"][0], d["smooth"]), color=c,
                          lw=1.0, label=lab)
        ax[1, 1].plot(d["te"], (d["etot"] - d["etot"][0]) / abs(d["etot"][0]), color=c,
                      lw=1.0, label=lab)
    ax[0, 0].axhline(a1_paper, color="k", ls=":", lw=1, label=f"{ref}, Table 2")
    ax[0, 0].set(ylim=(3e-3, 0.8), ylabel=r"$A_1/\Sigma_0$",
                 title=r"Growth of the $m=1$ mode (dashed: noise $\sqrt{\pi/N}$)")
    ax[0, 1].axhline(PAPER["r_c"], color="k", ls=":", lw=1, label=ref)
    ax[0, 1].set(ylabel=r"$|r_c|$", title="Central-mass offset")
    ax[1, 0].axhspan(0.20, 0.25, color="0.5", alpha=0.15, label=f"{ref}, 0.20-0.25")
    ax[1, 0].set(ylim=(0, 0.6), ylabel=r"$k = A_2/A_1$", title="Harmonic ratio (2-orbit mean)")
    ax[1, 1].set(ylabel=r"$\Delta E/|E_0|$", title="Energy conservation")
    for x in ax.flat:
        x.grid(True, alpha=0.3)
        x.legend(fontsize=9)
        x.axvspan(t0, t1, color="0.85", alpha=0.35, lw=0)
    for x in ax[1]:
        x.set_xlabel(r"$t\;[P_{\rm orb}]$")
    fig.suptitle(L.wrapped(title, 110) + f"\nagainst {ref}; shaded: window "
                 + rf"${t0:.0f}$-${t1:.0f}\,P_{{\rm orb}}$", fontsize=12)
    fig.tight_layout()
    out = Path(a.out) if a.out else Path(a.runs[0]).resolve().parent / "compare_runs.png"
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)

    rows = [(lab, summary(d, t0, t1)) for d, lab in zip(runs, labels)]
    print(f"window {t0:.0f}-{t1:.0f} orbits")
    print(f"{'':12s}" + "".join(f"{'run ' + str(i + 1):>9s}" for i in range(len(rows)))
          + f"{'B+26':>9s}")
    for i, (lab, _) in enumerate(rows):
        print(f"  run {i + 1}: {lab}")
    paper = dict(onset=np.nan, tau=np.nan, A1=a1_paper, A2=PAPER["A2"] / (N_PAPER - 1) * N_SEC,
                 A3=PAPER["A3"] / (N_PAPER - 1) * N_SEC, k=PAPER["k"], Q3=PAPER["Q3"],
                 A3c=np.nan, Q3c=np.nan, r_c=PAPER["r_c"], r_tb=PAPER["r_tb"])
    names = dict(onset="onset t", tau="e-fold [orb]", A1="A1/S0", A2="A2/S0", A3="A3/S0",
                 k="k", Q3="Q3", A3c="A3/S0 corr", Q3c="Q3 corr", r_c="|r_c|", r_tb="|r_tb|")
    for key, name in names.items():
        vals = [s[key] for _, s in rows] + [paper[key]]
        fmt = "9.1f" if key in ("onset", "tau") else "9.4f"
        print(f"{name:12s}" + "".join(f"{v:{fmt}}" for v in vals))
    print(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
