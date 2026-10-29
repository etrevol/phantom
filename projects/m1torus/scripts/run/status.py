#!/usr/bin/env python3
"""Health and progress of m1torus runs, one line each.

    status.py [RUN_DIR ...]          # default: every run under results/m1torus/ not finished
    status.py --all                  # finished runs too

Per run: status, orbits done of the total, wall rate (since start and since the previous
call), ETA, energy and angular-momentum drift, particles left, A_1/Sigma_0 and |r_c|,
and the problems below. Exit code 1 if any run has a problem.

Problems flagged
  STALL      no new analysis row for 10 min while RUNNING (a hang, or a full disk)
  SLOWDOWN   recent rate below a third of the average (timestep collapse)
  ENERGY     |dE/E| > 1e-2 (integration has gone wrong)
  ACCRETED   particles lost to the central mass (> 0.1%)
  ANGMOM     runs with colliding clouds, which radiate energy and feed the black hole by
             design, in place of ENERGY and ACCRETED: |dL/L| > 1e-3
  LOG        NaN, ERROR or FATAL in run.log
  DISK       less than 20 GB free on the host drive
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "common"))
import phantom_io as pio  # noqa: E402

P_ORB = 2.0 * np.pi
ROOT = Path(__file__).resolve().parents[4]


def fmt_h(seconds):
    if not np.isfinite(seconds) or seconds < 0:
        return "?"
    h, m = divmod(int(seconds) // 60, 60)
    return f"{h}h{m:02d}m" if h else f"{m}m"


def check(run: Path):
    status = (run / "status").read_text().strip() if (run / "status").exists() else "?"
    # tmax of the .in, not n_orbits of the .setup: a restart ends where its .in says
    total = float(pio.read_setup(run / "m1torus.in").get("tmax", np.nan)) / P_ORB
    m1f = run / "m1torus_m1.ev"
    problems = []
    if not m1f.exists():
        return status, f"{run.name}: {status}, no analysis output yet", problems
    m1 = pio.read_ev(m1f)
    t = m1["time"][-1] / P_ORB
    age = time.time() - m1f.stat().st_mtime
    start = datetime.fromisoformat(json.loads((run / "provenance.json").read_text())
                                   ["created"]).timestamp()
    t_first = m1["time"][0] / P_ORB
    # after continue.sh the rate counts from the restart, not across the pause before it
    clog = run / "continued.log"
    if clog.exists() and clog.read_text().strip():
        stamp, rest = clog.read_text().strip().splitlines()[-1].split(" ", 1)
        start = datetime.fromisoformat(stamp).timestamp()
        t_first = pio.read_dump(run / rest.split()[2])["header"]["time"] / P_ORB
    wall = (m1f.stat().st_mtime if status != "RUNNING" else time.time()) - start
    rate = (t - t_first) / wall * 3600.0 if wall > 0 else np.nan

    # recent rate from the previous call, kept beside the run
    state = run / ".status_prev"
    recent = np.nan
    if state.exists():
        w0, t0 = map(float, state.read_text().split())
        if time.time() - w0 > 60 and t > t0:
            recent = (t - t0) / (time.time() - w0) * 3600.0
    state.write_text(f"{time.time()} {t}")
    eta = (total - t) / (recent if np.isfinite(recent) else rate) * 3600.0

    ev = pio.read_ev_series(run, "m1torus")
    de = abs(ev["etot"][-1] / ev["etot"][0] - 1.0) if ev else np.nan
    dl = abs(ev["angtot"][-1] / ev["angtot"][0] - 1.0) if ev else np.nan
    n0, n = m1["n"][0], m1["n"][-1]
    # colliding clouds radiate energy and fall into the black hole by design: watch L,
    # which collisions and accretion keep, and report the radiated energy instead
    clouds = float(pio.read_setup(run / "m1torus.in").get("r_cloud", 0) or 0) > 0
    ediss = 0.0
    if clouds:
        dumps = sorted(run.glob("m1torus_[0-9][0-9][0-9][0-9][0-9]"))
        if dumps:
            ediss = pio.read_dump(dumps[-1])["header"].get("ediss", 0.0) / abs(ev["etot"][0])
    a1 = m1["A1"][-1] / m1["sigma0"][-1]
    rc = float(np.hypot(m1["xc"][-1], m1["yc"][-1]))

    if status == "RUNNING" and age > 600:
        problems.append(f"STALL({fmt_h(age)})")
    if np.isfinite(recent) and recent < rate / 3.0:
        problems.append("SLOWDOWN")
    if clouds:
        if dl > 1e-3:
            problems.append("ANGMOM")
    else:
        if de > 1e-2:
            problems.append("ENERGY")
        if n < 0.999 * n0:
            problems.append("ACCRETED")
    log = run / "run.log"
    # builds before ebce1648e logged an ERROR for the absent optional key; harmless
    text = log.read_text(errors="ignore").replace("ERROR: hsoft_nbody not found", "") if log.exists() else ""
    if re.search(r"\bNaN\b|ERROR|FATAL", text):
        problems.append("LOG")
    free = shutil.disk_usage("/mnt/d" if Path("/mnt/d").exists() else "/").free / 1e9
    if free < 20:
        problems.append(f"DISK({free:.0f}GB)")

    line = (f"{run.name}: {status} {t:6.1f}/{total:.0f} orb, "
            f"{rate:5.1f} orb/h (recent {recent:5.1f}), ETA {fmt_h(eta) if status == 'RUNNING' else '-'}, "
            + (f"E_diss/|E0| {ediss:.2f}, " if clouds else f"dE/E {de:.1e}, ")
            + f"dL/L {dl:.1e}, N {int(n)}/{int(n0)}, "
            f"A1/S0 {a1:.3f}, |r_c| {rc:.4f}"
            + (f"  !! {' '.join(problems)}" if problems else "  ok"))
    return status, line, problems


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="*")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args(argv)
    runs = [Path(r) for r in a.runs] or sorted(
        p.parent for p in (ROOT / "results" / "m1torus").glob("*/*/*/status"))
    bad = False
    for run in runs:
        st = (run / "status").read_text().strip() if (run / "status").exists() else "?"
        if not a.all and not a.runs and (st in ("COMPLETED", "FAILED") or st.startswith("STOPPED")):
            continue
        status, line, problems = check(run)
        print(line)
        bad |= bool(problems) or status == "FAILED"
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
