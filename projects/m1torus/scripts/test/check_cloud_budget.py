#!/usr/bin/env python3
"""Energy and angular-momentum budget of a run with colliding clouds.

Reads the energy from phantom's .ev files and the radiated energy (header 'ediss')
from the full dumps; prints, at each dump, E, E_diss, (E + E_diss - E0)/|E0| and
(L - L0)/L0. Collisions remove energy and keep L exactly, so E + E_diss and L
should drift only by the integrator error of a collisionless run.

    check_cloud_budget.py RUN_DIR
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "common"))
import phantom_io as pio  # noqa: E402

run = Path(sys.argv[1])
ev = pio.read_ev_series(run, prefix="m1torus")
t_ev, e_ev, l_ev = ev["time"], ev["etot"], ev["angtot"]
e0, l0 = e_ev[0], l_ev[0]
print(f"{'t':>8} {'E':>14} {'E_diss':>12} {'dE+Ed/|E0|':>12} {'dL/L0':>12}")
for d in sorted(run.glob("m1torus_[0-9][0-9][0-9][0-9][0-9]")):
    hdr = pio.read_dump(d)["header"]
    t, ed = hdr["time"], hdr.get("ediss", 0.0)
    k = np.argmin(abs(t_ev - t))
    print(f"{t:8.3f} {e_ev[k]:14.7e} {ed:12.5e} {(e_ev[k] + ed - e0) / abs(e0):12.3e}"
          f" {(l_ev[k] - l0) / l0:12.3e}")
