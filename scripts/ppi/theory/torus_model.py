#!/usr/bin/env python3
"""Papaloizou-Pringle torus: the model, its units, and the .setup file for phantom.

Running this file writes the .setup that `phantomsetup` reads:

    python3 scripts/ppi/theory/torus_model.py                       # RUN_SETUP below
    python3 scripts/ppi/theory/torus_model.py --preset pp84-2 --np 3e5 -o run.setup

The physics lives in the TorusModel field defaults and in PRESETS, nowhere else.
RUN_SETUP holds run mechanics only. Anything that needs the model imports this module.

Code units are those of the athena-collab torus: length L_0 = chi r_g, mass M_bh,
G = 1, so R_0 = r_center = 1 and the orbit at the density maximum lasts 2 pi.
"""

from __future__ import annotations

import argparse
import dataclasses
import math
import sys
from pathlib import Path

import numpy as np

# CGS, as in acc_disk_visc.cpp
C_LIGHT = 2.99792458e10
G_GRAV = 6.67430e-8
M_SUN = 1.98841e33
YR = 3.15576e7


@dataclasses.dataclass
class TorusModel:
    r_center: float = 1.0      # R_0, density maximum [L_0]
    C_prime: float = 0.2       # surface constant, 0 < C' < 1/2
    gamma: float = 1.3
    M_bh: float = 4.5e7        # [M_sun]
    chi: float = 500.0         # L_0 / r_g
    rho_0: float = 1.0e-13     # density at R_0 [g/cm^3]; sets the mass scale only

    # ---------------------------------------------------------------- shape

    @property
    def n(self):
        return 1.0 / (self.gamma - 1.0)

    @property
    def f_c(self):
        return 0.5 - self.C_prime

    @property
    def r_in(self):
        return self.r_center / (1.0 + math.sqrt(1.0 - 2.0 * self.C_prime))

    @property
    def r_out(self):
        return self.r_center / (1.0 - math.sqrt(1.0 - 2.0 * self.C_prime))

    @property
    def omega_0(self):
        return self.r_center**-1.5

    @property
    def p_orb(self):
        return 2.0 * math.pi / self.omega_0

    def f(self, R, z=0.0):
        R = np.asarray(R, float)
        return (self.r_center / np.sqrt(R * R + z * z)
                - 0.5 * (self.r_center / R) ** 2 - self.C_prime)

    def rho(self, R, z=0.0):
        """rho / rho(R_0, 0); zero outside the surface."""
        fv = np.clip(self.f(R, z), 0.0, None)
        return (fv / self.f_c) ** self.n

    def p_over_rho(self, R, z=0.0):
        return np.clip(self.f(R, z), 0.0, None) / ((self.n + 1.0) * self.r_center)

    def cs(self, R, z=0.0):
        return np.sqrt(self.gamma * self.p_over_rho(R, z))

    def z_surface(self, R):
        """Height of the surface above the mid-plane at cylindrical radius R."""
        R = np.asarray(R, float)
        s = 1.0 / (0.5 * self.r_center / R**2 + self.C_prime / self.r_center)
        return np.sqrt(np.clip((self.r_center * s) ** 2 - R * R, 0.0, None))

    @property
    def z_top(self):
        R = np.linspace(self.r_in, self.r_out, 4001)
        return float(self.z_surface(R).max())

    @property
    def h_over_r(self):
        """c_s / v_phi at (R_0, 0), = sqrt((gamma-1)(1/2-C')) as in athena-collab."""
        return math.sqrt((self.gamma - 1.0) * self.f_c)

    def volume_integral(self, nr=800, nz=400):
        """int (rho/rho_c) dV over the torus, by the midpoint rule in (R, z)."""
        R = np.linspace(self.r_in, self.r_out, nr + 1)
        Rc, dR = 0.5 * (R[1:] + R[:-1]), np.diff(R)
        z = np.linspace(0.0, self.z_top, nz + 1)
        zc, dz = 0.5 * (z[1:] + z[:-1]), np.diff(z)
        rr, zz = np.meshgrid(Rc, zc, indexing="ij")
        return float(4.0 * math.pi * (self.rho(rr, zz) * rr * dR[:, None]
                                       * dz[None, :]).sum())

    # ---------------------------------------------------------------- units

    @property
    def r_g(self):
        return 2.0 * G_GRAV * self.M_bh * M_SUN / C_LIGHT**2

    @property
    def L_0(self):
        return self.chi * self.r_g

    @property
    def t_unit(self):
        return math.sqrt(self.L_0**3 / (G_GRAV * self.M_bh * M_SUN))

    @property
    def rho_unit(self):
        return self.M_bh * M_SUN / self.L_0**3

    @property
    def rho_c_code(self):
        return self.rho_0 / self.rho_unit

    @property
    def m_torus(self):
        """Torus mass in code units (M_bh)."""
        return self.rho_c_code * self.volume_integral()

    def check(self):
        if not 0.0 < self.C_prime < 0.5:
            raise ValueError(f"C' = {self.C_prime} outside (0, 1/2): no bounded torus")
        if not 1.0 < self.gamma < 2.0:
            raise ValueError(f"gamma = {self.gamma} outside (1, 2)")


# PP84 section 6 uses n = 3; d = 1/(2C') is their distortion parameter.
PRESETS = {
    "athena": {},
    "pp84-1": dict(C_prime=0.495, gamma=4.0 / 3.0),
    "pp84-2": dict(C_prime=0.444, gamma=4.0 / 3.0),
    "pp84-3": dict(C_prime=0.25, gamma=4.0 / 3.0),
}

# PP84 Table 1, growth time in units of 1/Omega_0, methods A and B.
PP84_TABLE1 = {
    ("pp84-1", 2): (8.0, 8.5), ("pp84-1", 4): (3.9, 4.7), ("pp84-1", 6): (2.9, 4.0),
    ("pp84-2", 2): (4.0, 4.6), ("pp84-3", 2): (8.5, 10.7),
}

RUN_SETUP = dict(
    preset="athena",
    np=100_000,           # target particle number
    lattice="rings",      # rings | random
    nsym=32,              # ring particle counts are multiples of nsym; 1 = no symmetry
    pert_kind="none",     # none | single | multi | noise
    pert_amp=0.0,         # delta v_R in units of the local c_s
    pert_m=2,             # mode (single) or highest mode (multi)
    pert_seed=1,
    relax=True,           # stage 1 at rest in the effective potential, then spin up
    n_relax=2.0,          # relaxation time [orbits]
    damp_relax=0.3,       # velocity damping rate [code units]
    r_acc=0.35,           # accretion radius of the central mass [L_0]
    n_orbits=40.0,
    n_analysis_orbit=20,  # live analysis calls per orbit (sets dtmax)
    n_dumps_orbit=2,      # full dumps per orbit; must divide n_analysis_orbit
)


def model_from_preset(name):
    if name not in PRESETS:
        raise ValueError(f"unknown preset '{name}', one of {sorted(PRESETS)}")
    m = TorusModel(**PRESETS[name])
    m.check()
    return m


def _fmt(v):
    if isinstance(v, str):
        return v
    if isinstance(v, (int, np.integer)):
        return f"{v:d}"
    return f"{v:.10g}"


def setup_lines(model, run):
    """The .setup file, keyword for keyword what setup_ppitorus.f90 reads."""
    if run["n_analysis_orbit"] % run["n_dumps_orbit"]:
        raise ValueError("n_dumps_orbit must divide n_analysis_orbit")
    if run["relax"] and model.r_center != 1.0:
        raise ValueError("relax needs r_center = 1 (iexternalforce = 5 is written for it)")
    if run["r_acc"] >= model.r_in:
        raise ValueError(f"r_acc = {run['r_acc']} must lie inside r_in = {model.r_in:.4f}")
    rows = [
        ("# units", None, None),
        ("M_bh", model.M_bh, "central mass [M_sun]"),
        ("chi", model.chi, "length unit L_0 = chi r_g"),
        ("rho_0", model.rho_0, "density at the maximum [g/cm^3]"),
        ("# torus", None, None),
        ("r_center", model.r_center, "radius of the density maximum [L_0]"),
        ("C_prime", model.C_prime, "surface constant (0 < C_prime < 0.5)"),
        ("gamma", model.gamma, "adiabatic index"),
        ("r_acc", run["r_acc"], "accretion radius of the central mass [L_0]"),
        ("# particles", None, None),
        ("np", int(run["np"]), "target number of particles"),
        ("lattice", run["lattice"], "rings | random"),
        ("nsym", int(run["nsym"]), "ring counts are multiples of nsym (1 = no symmetry)"),
        ("# seed perturbation", None, None),
        ("pert_kind", run["pert_kind"], "none | single | multi | noise"),
        ("pert_amp", run["pert_amp"], "delta v_R in units of the local c_s"),
        ("pert_m", int(run["pert_m"]), "azimuthal mode (single) or highest mode (multi)"),
        ("pert_seed", int(run["pert_seed"]), "random seed (multi, noise, random)"),
        ("# relaxation (stage 1, at rest in the effective potential)", None, None),
        ("relax", "T" if run["relax"] else "F", "relax before spinning up (needs r_center = 1)"),
        ("n_relax", run["n_relax"], "relaxation time [orbits at r_center]"),
        ("damp_relax", run["damp_relax"], "velocity damping rate [code units]"),
        ("# run", None, None),
        ("n_orbits", run["n_orbits"], "run length [orbits at r_center]"),
        ("n_analysis_orbit", int(run["n_analysis_orbit"]), "analysis calls per orbit"),
        ("n_dumps_orbit", int(run["n_dumps_orbit"]), "full dumps per orbit"),
    ]
    out = ["# input file for the Papaloizou-Pringle torus setup (ppitorus)"]
    for key, val, com in rows:
        if val is None:
            out += ["", key]
        else:
            out.append(f"{key:>18} = {_fmt(val):<14}   ! {com}")
    return out


def summary(model, run=None):
    rows = [
        f"preset parameters : C' = {model.C_prime}, gamma = {model.gamma:.6g}, n = {model.n:.4g}",
        f"surface           : r_in = {model.r_in:.6f}, r_out = {model.r_out:.6f}, "
        f"z_top = {model.z_top:.6f}",
        f"thickness         : c_s/v_phi at R_0 = {model.h_over_r:.4f}, "
        f"d = 1/(2C') = {0.5 / model.C_prime:.4f}",
        f"orbit at R_0      : P = {model.p_orb:.6f} code = "
        f"{model.p_orb * model.t_unit / YR:.4g} yr",
        f"units             : L_0 = {model.L_0:.6e} cm, t = {model.t_unit:.6e} s, "
        f"rho = {model.rho_unit:.6e} g/cm^3",
        f"torus mass        : {model.m_torus:.6e} M_bh = "
        f"{model.m_torus * model.M_bh:.6g} M_sun",
    ]
    if run is not None:
        rows.append(f"particle mass     : {model.m_torus / run['np']:.6e} (np = {run['np']:.3g})")
    return "\n".join(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--output", default="ppitorus.setup")
    ap.add_argument("--preset", default=RUN_SETUP["preset"], choices=sorted(PRESETS))
    for key, val in RUN_SETUP.items():
        if key == "preset":
            continue
        if isinstance(val, bool):
            ap.add_argument("--" + key.replace("_", "-"), dest=key, default=val,
                            action=argparse.BooleanOptionalAction)
            continue
        typ = float if isinstance(val, float) or key == "np" else type(val)
        ap.add_argument("--" + key.replace("_", "-"), dest=key, type=typ, default=val)
    a = ap.parse_args(argv)
    run = {k: getattr(a, k) for k in RUN_SETUP}
    run["np"] = int(run["np"])
    model = model_from_preset(a.preset)
    lines = setup_lines(model, run)
    Path(a.output).write_text("\n".join(lines) + "\n")
    print(summary(model, run))
    print(f"wrote {a.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
