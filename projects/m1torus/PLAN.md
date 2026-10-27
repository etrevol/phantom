# Plan

Ground truth here is the paper and our own measurements. The side projects on this
laptop (athena-torus, athena-vv, phantom-sgtorus) are reference material only: their
code and numbers are unverified and nothing is taken from them as correct. Of
athena-collab only the branch `artem-develop` counts as checked.

## A. Reproduction of the collisionless mode

### What the paper does (Sect. 2, Table 1, App. A)

- N = 128k equal-mass particles, torus mass M_tor = 0.1, central mass M_c = 1, units
  G = M_c = R_tor = 1, orbit at R_tor = 2 pi; time quoted in orbits.
- Keplerian elements about M_c: a uniform in [1 - R0, 1 + R0], e uniform in [0, e_max],
  i uniform in [-i_max, i_max], node and true anomaly uniform in [0, 2 pi], argument of
  periapsis 0. Canonical run `0.5-60-0.3`: e_max = 0.5, i_max = 60 deg, R0 = 0.3.
- Softening eps = 1e-2 (1e-3 and 1e-4 give the same result, App. A); Hermite 4th order,
  block steps, eta = 0.007.
- Diagnostics: counts in 1200 azimuthal sectors over the whole torus, least-squares
  harmonics up to m = 5; amplitude A_m in particles per sector, Sigma_0 = (N - 1)/1200;
  barycentre of the torus r_tb and of the central mass r_c = -M_tor r_tb; osculating
  elements; k = A_2/A_1, Q_3 = (A_3/A_1)/(A_2/A_1)^2.

### Checks (falsifiable; canonical run unless stated)

| quantity | paper | where |
|---|---|---|
| onset of m = 1 | after virial equilibrium; axisymmetric at t = 35, lopsided at t = 51 | Fig. 4, 9 |
| A_1 (A_1/Sigma_0) | 40.2 (0.377) | Table 2 |
| A_2, A_3 | 8.46, 1.77 | Table 2 |
| k, Q_3 | 0.21, 0.99 | Table 2 |
| central-mass offset r_c | ~0.024 | Sect. 9.3, Fig. 3 |
| torus barycentre r_tb | 0.24 | Sect. 9.2 |
| eccentricity distribution | Rayleigh, sigma = 0.27 | App. B.1 |
| pattern | prograde, period >> orbital | Sect. 3.1 |
| thin control `0.5-30-0.3` | no mode above noise | Sect. 5, Table 2 |
| mass scaling | r_c linear in M_tor (0.10, 0.06, 0.01) | Fig. 6 |

### Differences from the paper that must be measured, not assumed

- Gravity by tree (approximate force) instead of direct summation.
- Softening of collisionless particles is adaptive (from h) in phantom, fixed in the paper.
- Integrator: leapfrog with individual steps, not 4th-order Hermite.

### Steps

1. `setup_m1torus.f90`: elements to positions and velocities, collisionless particles,
   the central mass as a sink; `analysis_m1torus.f90`: the paper's diagnostics, live.
2. Numerics at 16k: energy and angular momentum drift, cost per orbit, effect of the
   tree opening angle and of softening.
3. Cost at 32k and 64k, extrapolation to 128k x 1000 orbits. **Report before any run
   longer than 30 min.**
4. Canonical run, then the thin control and one torus mass for the scaling.

## B. Dissipation (after A)

A fraction f_gas of the torus as SPH gas with short cooling (near isothermal), the fluid
limit of colliding clouds; collisionless rest unchanged. Does the mode survive, at what
f_gas does it die, where does the gas go? Viscosity options and their meaning: local
notes, `.notes/projects/m1torus.md`.

## C. Nested gas disc (after B)

A gas disc at 0.1-0.3 R_tor inside the lopsided torus, coplanar and tilted by 10-30 deg.
Measure forced eccentricity, inflow, precession, warp, tearing. Literature anchors:
Hopkins & Quataert (2010) on m = 1 driven inflow; Winkel et al. (2025, arXiv:2506.19902)
for NGC 4593; nested misaligned rings in NGC 1068 (Vermot 2022, arXiv:2204.03621); the
Laplace surface (Tremaine, Touma & Kazandjian 2009) - to be checked when C starts.
