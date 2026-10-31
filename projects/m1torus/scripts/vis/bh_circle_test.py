#!/usr/bin/env python3
"""Do the oscillations of the black-hole offset come from a drifting centre of its circle?

In the saturated m = 1 state the black hole moves on a circle of radius A about the
centre of mass of the bound system, turning with the pattern speed Omega_p. If that
centre C(t) drifts away from the reference point (the origin), the distance from the
origin is
    |r_bh|^2 = A^2 + |C|^2 + 2 A |C| cos(Omega_p t - phi_C),
so |r_bh| oscillates with the pattern period between A - |C| and A + |C|. The test needs
only the trajectory (t, x_bh, y_bh): in windows of one pattern period it fits a circle
(algebraic least squares) for its centre C and radius A, then checks
  1. A(t) stays constant (the mode does not change),
  2. half the swing of |r_bh| in each window equals |C| (the oscillation is the drift),
  3. |r_bh - C(t)| no longer oscillates.
A real oscillation of the mode itself would show up as a changing A or as a swing of
|r_bh - C| left over after the drift is removed.

    bh_circle_test.py FILE [--cols t,x,y] [--period P] [--tmin T] [--out PNG]
    bh_circle_test.py --demo            # synthetic drift, as a check of the method
FILE: whitespace columns, time in orbits (or --time-unit 2pi for N-body time).
"""
import argparse
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def fit_circle(x, y):
    """Centre and radius from x^2 + y^2 = 2 a x + 2 b y + c (Kasa)."""
    M = np.c_[2 * x, 2 * y, np.ones_like(x)]
    (a, b, c), *_ = np.linalg.lstsq(M, x ** 2 + y ** 2, rcond=None)
    return a, b, np.sqrt(c + a ** 2 + b ** 2)


def period_from_phase(t, x, y, a, b):
    """Pattern period from the unwrapped angle of the black hole about the centre."""
    ph = np.unwrap(np.arctan2(y - b, x - a))
    return abs(2 * np.pi / np.polyfit(t, ph, 1)[0])


def analyse(t, x, y, period=None, tmin=None, step=0.25):
    k = t >= (tmin if tmin is not None else t[0])
    t, x, y = t[k], x[k], y[k]
    if period is None:          # one global fit over the last third for a first guess
        s = t > t[0] + 2 * (t[-1] - t[0]) / 3
        a, b, _ = fit_circle(x[s], y[s])
        period = period_from_phase(t[s], x[s], y[s], a, b)
    rows = []
    t0 = t[0]
    while t0 + period <= t[-1]:
        s = (t >= t0) & (t < t0 + period)
        if s.sum() > 20:
            a, b, A = fit_circle(x[s], y[s])
            r = np.hypot(x[s], y[s])
            rc = np.hypot(x[s] - a, y[s] - b)
            rows.append((t0 + period / 2, a, b, A, 0.5 * np.ptp(r), np.hypot(a, b), 0.5 * np.ptp(rc), r.mean()))
        t0 += step * period
    return period, np.array(rows), (t, x, y)


def demo():
    """Circle of radius 0.0225, period 54 orbits, centre drifting from t = 250 to
    (-0.003, -0.02) at t = 1000, plus noise: like the eps = 1e-4 run of the authors."""
    rng = np.random.default_rng(1)
    t = np.arange(100, 1000, 0.1)
    om = 2 * np.pi / 54.
    f = np.clip((t - 250) / 750, 0, 1)
    cx, cy = -0.003 * f, -0.020 * f
    x = cx + 0.0225 * np.cos(om * t) + 4e-4 * rng.standard_normal(t.size)
    y = cy + 0.0225 * np.sin(om * t) + 4e-4 * rng.standard_normal(t.size)
    return t, x, y


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("file", nargs="?")
    ap.add_argument("--cols", default="0,1,2", help="columns of t, x, y (0-based)")
    ap.add_argument("--time-unit", default="orbit", choices=("orbit", "2pi"))
    ap.add_argument("--period", type=float, default=None, help="pattern period [orbits]")
    ap.add_argument("--tmin", type=float, default=None)
    ap.add_argument("--demo", action="store_true")
    ap.add_argument("--out", default="bh_circle_test.png")
    a = ap.parse_args()
    if a.demo:
        t, x, y = demo()
    else:
        c = [int(i) for i in a.cols.split(",")]
        d = np.loadtxt(a.file, usecols=c)
        t, x, y = d[:, 0], d[:, 1], d[:, 2]
        if a.time_unit == "2pi":
            t = t / (2 * np.pi)
    P, w, (t, x, y) = analyse(t, x, y, a.period, a.tmin)
    if len(w) < 4:
        sys.exit(f"only {len(w)} window(s) of one pattern period ({P:.1f} orbits): too short a series")
    tm, cx, cy, A, half_r, C, half_rc, rmean = w.T
    late = tm > tm[0] + 0.5 * (tm[-1] - tm[0])
    print(f"pattern period {P:.1f} orbits; windows {len(tm)}")
    print(f"circle radius A: {A.mean():.4f} +- {A.std():.4f} (whole run), late half {A[late].mean():.4f}")
    print(f"half swing of |r_bh| from the origin, late half: {half_r[late].mean():.4f}")
    print(f"drift of the centre |C|, late half:              {C[late].mean():.4f}")
    print(f"half swing of |r_bh - C| (what is left), late half: {half_rc[late].mean():.4f}")
    print(f"correlation of the swing of |r_bh| with |C| over the windows: {np.corrcoef(half_r, C)[0, 1]:.3f}")

    Ci = np.c_[np.interp(t, tm, cx), np.interp(t, tm, cy)]
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8), layout="constrained")
    ax[0].plot(x, y, lw=0.4, color="C2")
    ax[0].plot(cx, cy, "k.-", ms=3, lw=1, label="fitted centre C(t)")
    ax[0].set(aspect="equal", xlabel=r"$x_{\rm bh}$", ylabel=r"$y_{\rm bh}$",
              title="Black-hole trajectory and the centre of its circle")
    ax[0].legend(fontsize=8, frameon=False)
    ax[1].plot(t, np.hypot(x, y), lw=0.5, color="0.3", label=r"$|r_{\rm bh}|$ from the origin")
    ax[1].plot(tm, A + C, "C3", tm, A - C, "C3", lw=1.2)
    ax[1].plot([], [], "C3", label=r"$A \pm |C|$ predicted from the drift")
    ax[1].set(xlabel="t [orbits]", ylabel="offset", title="Oscillation and its prediction from the drift")
    ax[1].legend(fontsize=8, frameon=False)
    ax[2].plot(t, np.hypot(x - Ci[:, 0], y - Ci[:, 1]), lw=0.5, color="C0",
               label=r"$|r_{\rm bh} - C(t)|$")
    ax[2].plot(t, np.hypot(x, y), lw=0.3, color="0.75", label="from the origin (for comparison)")
    ax[2].set(xlabel="t [orbits]", ylabel="offset", title="Offset from the drifting centre")
    ax[2].legend(fontsize=8, frameon=False)
    fig.savefig(a.out, dpi=130)
    print(a.out)


if __name__ == "__main__":
    sys.exit(main())
