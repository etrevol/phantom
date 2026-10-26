#!/usr/bin/env python3
"""Azimuthal mode growth of a PPI torus run against the linear theory.

Reads <prefix>_ppi.ev (written every dtmax by the live analysis), plots |A_m|/A_0 per
orbit on a log axis, fits the straightest exponential window below the linearity
ceiling, and compares the rate and pattern speed with ppi3d.py. Linear theory is cached
in <run>/ppi_theory.json, since a converged 3D spectrum costs about a minute per mode.

    vismodes.py                      # m = 1..4, figure and table
    vismodes.py --m 1,2,3,4,8,16 --no-theory
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ppi_data as pdat  # noqa: E402

#: |A_m| / A_0 above which the perturbation is not small against the axisymmetric
#: background and linear theory stops describing it (vvlib/ppi.py, athena-vv)
LINEAR_CEILING = 0.02
COLORS = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#8c564b", "#e377c2",
          "#7f7f7f", "#bcbd22", "#17becf"]


def fit_growth(t, a, lo, hi):
    """Least-squares d ln A / dt over [lo, hi], with R^2."""
    sel = (t >= lo) & (t <= hi) & (a > 0)
    if sel.sum() < 4:
        return dict(rate=np.nan, r2=np.nan, n=int(sel.sum()), window=(lo, hi))
    x, y = t[sel], np.log(a[sel])
    slope, icept = np.polyfit(x, y, 1)
    resid = y - (slope * x + icept)
    ss = ((y - y.mean()) ** 2).sum()
    r2 = 1.0 - (resid**2).sum() / ss if ss > 0 else np.nan
    return dict(rate=float(slope), r2=float(r2), n=int(sel.sum()),
                window=(float(x[0]), float(x[-1])), icept=float(icept))


def best_window(t, a, background=None, min_span=2.0, min_points=6, ceiling=LINEAR_CEILING,
                floor=None, stride=1):
    """Straightest window of at least min_span, below the ceiling and above the floor.

    As in vvlib/ppi.py: scanning finds the exponential phase from the data instead of
    assuming one, and the ceiling keeps the scan out of the weakly non-linear top where
    a saturating curve is straightest but the rate already depressed. The floor keeps
    it out of the noise a seed has to rise through first.
    """
    good = a > 0
    if background is not None:
        rel = a / np.where(background > 0, background, np.inf)
        good &= rel < ceiling
        if floor is not None:
            good &= rel > floor
    tt, aa = t[good], a[good]
    best = dict(rate=np.nan, r2=-np.inf, n=0, window=(np.nan, np.nan))
    for i in range(0, len(tt), stride):
        for j in range(i + min_points - 1, len(tt), stride):
            if tt[j] - tt[i] < min_span:
                continue
            fit = fit_growth(tt, aa, tt[i], tt[j])
            if fit["r2"] > best["r2"] + 1e-6 or (
                    abs(fit["r2"] - best["r2"]) <= 1e-6 and fit["n"] > best["n"]):
                best = fit
    return best


def pattern_speed(t, phase, m, window):
    """-d(phase)/dt / m over the window; phase is arg sum exp(-i m phi) at r_center."""
    sel = (t >= window[0]) & (t <= window[1])
    if sel.sum() < 4:
        return np.nan
    ph = np.unwrap(phase[sel])
    return float(-np.polyfit(t[sel], ph, 1)[0] / m)


def theory(info, ms, nr=80, nz=40):
    """Linear rates per orbit and pattern speeds from ppi3d, cached in the run."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "theory"))
    import ppi3d

    cache = info["dir"] / "ppi_theory.json"
    data = json.loads(cache.read_text()) if cache.exists() else {}
    model = info["model"]
    key = f"C{model.C_prime:g}_g{model.gamma:.6g}_{nr}x{nz}"
    tab = data.setdefault(key, {})
    for m in ms:
        if str(m) in tab:
            continue
        r = ppi3d.converged(model, m, nr, nz)
        tab[str(m)] = dict(rate_orbit=r["efolds_per_orbit"], pattern=r["pattern"],
                           rel_change=r["rel_change"], grid=[nr, nz])
        cache.write_text(json.dumps(data, indent=1))
    return {int(k): v for k, v in tab.items()}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", default=".")
    ap.add_argument("--m", default="1,2,3,4")
    ap.add_argument("--no-theory", action="store_true")
    ap.add_argument("--floor", type=float, default=None,
                    help="ignore |A_m|/A_0 below this in the fit (noise of the seed)")
    ap.add_argument("--min-span", type=float, default=2.0, help="orbits")
    ap.add_argument("--flux", action="store_true", help="also plot the v_R mass-flux modes")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)

    info = pdat.run_info(a.run)
    ev = pdat.read_ev(info["modes"])
    p_orb = info["model"].p_orb
    t = ev["time"] / p_orb
    ms = [int(s) for s in a.m.split(",")]
    th = {} if a.no_theory else theory(info, [m for m in ms if m <= 6])

    out_dir = Path(a.out) if a.out else info["dir"] / "figs"
    out_dir.mkdir(parents=True, exist_ok=True)
    nrow = 2 if a.flux else 1
    fig, axes = plt.subplots(nrow, 1, figsize=(11, 6.5 * nrow), squeeze=False)
    ax = axes[0, 0]
    a0 = ev["A0"]
    rows = []
    for k, m in enumerate(ms):
        c = COLORS[k % len(COLORS)]
        am = ev[f"A{m}"]
        rel = am / a0
        ax.semilogy(t, rel, color=c, lw=1.6, label=f"m = {m}")
        fit = best_window(t, am, a0, min_span=a.min_span, floor=a.floor,
                          stride=max(1, len(t) // 300))
        tr = th.get(m)
        # a rate that moves by more than 5% on halving the grid is not a prediction
        if tr and not (tr["rate_orbit"] > 0 and tr["rel_change"] < 0.05):
            tr = None
        patt = (pattern_speed(ev["time"], ev[f"phi{m}"], m,
                              tuple(w * p_orb for w in fit["window"]))
                if m <= 4 and np.isfinite(fit["rate"]) else np.nan)
        rows.append((m, fit, tr, patt))
        if np.isfinite(fit["rate"]):
            w = np.array(fit["window"])
            ax.semilogy(w, np.exp(fit["icept"] + fit["rate"] * w) / np.interp(w, t, a0),
                        color=c, lw=5, alpha=0.25)
        if tr and np.isfinite(fit["rate"]):
            # theory slope through the middle of the fitted window
            tm_ = 0.5 * (fit["window"][0] + fit["window"][1])
            amid = np.exp(fit["icept"] + fit["rate"] * tm_) / np.interp(tm_, t, a0)
            tt = np.linspace(t[0], t[-1], 50)
            ax.semilogy(tt, amid * np.exp(tr["rate_orbit"] * (tt - tm_)), color=c, ls=":",
                        lw=1.2)
    ax.axhline(LINEAR_CEILING, color="k", ls="--", lw=0.8, alpha=0.6)
    ax.set_xlabel(r"$t\;[P_{\rm orb}]$")
    ax.set_ylabel(r"$|A_m|/A_0$")
    ax.grid(True, alpha=0.3)
    ax.legend(ncol=min(len(ms), 4), fontsize=9, loc="lower right")
    ax.set_title("Azimuthal modes of the surface density"
                 + ("  (dotted: linear theory, PP84 eq. 4.1)" if th else ""),
                 fontsize=12, fontweight="bold")
    if a.flux:
        axf = axes[1, 0]
        for k, m in enumerate(ms):
            axf.semilogy(t, ev[f"F{m}"] / np.maximum(ev["A0"], 1e-300),
                         color=COLORS[k % len(COLORS)], lw=1.4, label=f"m = {m}")
        axf.set_xlabel(r"$t\;[P_{\rm orb}]$")
        axf.set_ylabel(r"$|F_m|/A_0$  ($\Sigma v_R$ modes)")
        axf.grid(True, alpha=0.3)
        axf.legend(ncol=4, fontsize=9)
    fig.text(0.012, 0.988, "\n".join(pdat.info_lines(info)), ha="left", va="top",
             fontsize=8.5, bbox=dict(boxstyle="round", fc="white", ec="#b8b8b4", alpha=0.9))
    fig.tight_layout()
    f = out_dir / "modes.png"
    fig.savefig(f, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"{'m':>3} {'fit/orb':>8} {'theory':>8} {'ratio':>6} {'R2':>7} {'window [orb]':>15}"
          f" {'patt':>7} {'patt th':>7}")
    for m, fit, tr, patt in rows:
        r_th = tr["rate_orbit"] if tr else np.nan
        p_th = tr["pattern"] if tr else np.nan
        ratio = fit["rate"] / r_th if tr and r_th > 0 else np.nan
        w = fit["window"]
        print(f"{m:3d} {fit['rate']:8.3f} {r_th:8.3f} {ratio:6.2f} {fit['r2']:7.4f} "
              f"{w[0]:6.2f}-{w[1]:6.2f}  {patt / info['model'].omega_0:7.3f} {p_th:7.3f}")
    print(f)
    return 0


if __name__ == "__main__":
    sys.exit(main())
