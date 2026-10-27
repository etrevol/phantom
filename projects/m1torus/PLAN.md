# Plan

The single tracked list of tasks of this project. Status: `[x]` done, `[~]` in progress,
`[ ]` open. Each task has an acceptance check and a cost on this laptop (8 cores).
Results and decisions, dated, go to [`LOG.md`](LOG.md); this file says what is left.

Ground truth is the paper and our own measurements. The side projects on this laptop
(athena-torus, athena-vv, phantom-sgtorus) are reference only; of athena-collab only
`artem-develop` counts as checked.

Cost reference: 16k ~45 s/orbit (150 orbits ~2 h), 64k ~4.5 min/orbit (110 orbits ~9 h,
150 ~12-14 h). Ask the user before anything over 30 min.

## A. Reproduction of the collisionless mode

| id | status | task | acceptance | cost |
|---|---|---|---|---|
| A1 | [x] | setup, live analysis, run tooling, figures | runs end to end | - |
| A2 | [x] | numerics at 16k: cost, conservation, energy-drift diagnosis | dE/E < 2e-3 over 150 orbits | - |
| A3 | [x] | canonical `0.5-60-0.3`, 16k, 150 orbits | mode appears and holds | 2 h |
| A4 | [ ] | canonical at 64k, 110-150 orbits (the user starts it overnight) | saturated A1/S0, k, r_c within 10% of the paper; 16k -> 64k converged | 9-14 h |
| A5 | [ ] | thin control `0.5-30-0.3`, 16k, 150 orbits | no m = 1 above noise (paper Sect. 5) | 2 h |
| A6 | [ ] | weak case `0.5-45-0.3`, 16k, 150 orbits | A1/S0 ~0.09, k ~0.08 (Table 2) | 2 h |
| A7 | [ ] | torus mass M_tor = 0.06 and 0.01, 16k | r_c linear in M_tor with r_tb fixed (Fig. 6, Sect. 3.3); pattern speed linear in M_tor (Sect. 10); onset may come later at low mass, check Fig. 6 before setting the length | 2 x 2 h |
| A8 | [ ] | `0.2-60-0.3`, 16k (optional) | stronger mode, A1/S0 ~0.46 | 2 h |
| A9 | [x] | pericentres of e >= 0.5 sit at Phi_1 + 180 deg here; the paper says "near the phase": settle the convention | explanation in LOG | - |
| A10 | [ ] | short write-up of A: methods, differences from the paper, results | LaTeX, like `projects/ppi/formulary` | - |
| A11 | [ ] | M_tor = 0.3 (the NGC 613 mass ratio), 16k | does r_tb stay ~0.24 (Sect. 3.3) or grow (Sect. 9.2); cavity-AGN offset (1 + M_tor) r_tb against the observed 3/8 = 0.375 R_tor | 2 h |

### What the paper does (Sect. 2, Table 1, App. A)

- N = 128k equal-mass particles, M_tor = 0.1, M_c = 1, units G = M_c = R_tor = 1, the
  orbit at R_tor lasts 2 pi; time quoted in orbits.
- Keplerian elements about M_c: a uniform in [1 - R0, 1 + R0], e uniform in [0, e_max],
  i uniform in [-i_max, i_max], node and true anomaly uniform in [0, 2 pi], argument of
  periapsis 0. Canonical `0.5-60-0.3`: e_max = 0.5, i_max = 60 deg, R0 = 0.3.
- Softening eps = 1e-2 (1e-3 and 1e-4 give the same result, App. A); 4th-order Hermite,
  block steps, eta = 0.007.
- Diagnostics: counts in 1200 azimuthal sectors, least-squares harmonics to m = 5; A_m in
  particles per sector, Sigma_0 = (N - 1)/1200; r_tb and r_c = -M_tor r_tb; osculating
  elements; k = A_2/A_1, Q_3 = (A_3/A_1)/(A_2/A_1)^2.

### Checks (canonical run)

| quantity | paper | 16k, eps_c = 0 (A3) | where |
|---|---|---|---|
| onset of m = 1 | after virial equilibrium; axisymmetric t = 35, lopsided t = 51 | t ~ 35-40 | Fig. 4, 9 |
| A_1/Sigma_0 | 0.368 (A_1 = 40.2) | 0.349 | Table 2 |
| A_2/Sigma_0, A_3/Sigma_0 | 0.077, 0.016 | 0.073, 0.020 | Table 2 |
| k, Q_3 | 0.21, 0.99 | 0.209 (Poisson fit 0.200), 1.34 (A_3 at the 16k noise floor) | Table 2, Sect. 7.2 |
| r_c | ~0.024 | 0.0217, flat to t = 150 | Sect. 9.3, Fig. 3 |
| r_tb | 0.24 | 0.222 | Sect. 9.2 |
| eccentricities | Rayleigh, sigma = 0.27 | sigma = 0.29 | App. B.1 |
| pattern | prograde, period >> orbital | 0.0186 Omega_0, period 54 orbits | Sect. 3.1 |

### Differences from the paper (measured, not assumed)

- Central mass unsoftened (eps_c = 0): phantom softens the particle-sink force over
  max(eps_c, h_i) with the adaptive h_i, which does not conserve energy (LOG 2026-10-27).
  The paper's eps = 0.01 is negligible against pericentres >= 0.05.
- Gravity by tree, adaptive softening between torus particles (h ~ 0.08), leapfrog with
  individual steps; the paper: direct summation, eps = 0.01, Hermite.

## B. Dissipation (after A)

A fraction f_gas of the torus particles becomes SPH gas, the rest stays collisionless. At
what f_gas does the mode die, and where does the gas go? Design (B1, LOG 2026-10-27): gas
locally isothermal with H/R = 0.1, the effective turbulent sound speed of Hopkins &
Quataert (2010); no explicit viscosity, angular momentum moves by the resolved
gravitational torques and shocks (shock viscosity with its switch); gas self-gravity on;
the sink accretes gas only, r_acc = 0.05. SPH is the frequent-collision limit of a cloud
ensemble; slower collisions would need sticky particles, which phantom lacks.

| id | status | task | acceptance | cost |
|---|---|---|---|---|
| B1 | [x] | design with the user | decision in LOG | - |
| B2 | [ ] | moddump: f_gas of the particles to gas, each split into lighter gas particles (N_gas ~ 32k); control run of the gas without the torus | dE/E as in A2; the control gives the numerical inflow; cost per orbit | 1 h |
| B3 | [ ] | scan f_gas = 0.05, 0.1, 0.3 from the A3 dump at t = 100, ~50 orbits each | stellar A1/S0 and r_c against f_gas | 3 x 2-3 h |
| B4 | [ ] | at the threshold f_gas: c_s = 0.3 v_K; one beta-cooling run (Gammie 2001); one run from t = 0 | does the threshold depend on the gas model; does the mode form at all | 3 x 2-4 h |
| B5 | [ ] | where the gas goes: disc radius, eccentricity, alignment with the mode, inflow to the sink | figures, LOG; compare with Hopkins & Quataert (2010, 2011) | - |

Viscosity and dissipation options in phantom: local note `projects/m1torus/CLAUDE.md`.
Literature (in `materials/papers/`): Hopkins & Quataert 2010, 2011; Hopkins et al. 2012;
Krolik & Begelman 1988; Gammie 2001; Meru & Bate 2011; Lodato & Price 2010; Hönig 2019;
Combes et al. 2019.

## C. Nested gas disc (after B)

A gas disc at 0.1-0.3 R_tor inside the lopsided torus, coplanar and tilted by 10-30 deg:
forced eccentricity, inflow, precession, warp, tearing.

| id | status | task | acceptance | cost |
|---|---|---|---|---|
| C1 | [ ] | literature: Hopkins & Quataert 2010, Laplace surface (Tremaine et al. 2009), disc tearing | short summary in LOG | - |
| C2 | [ ] | setup: disc inside the evolved A torus (moddump onto an A dump) | stable disc without the torus | 1 h |
| C3 | [ ] | coplanar and tilted runs, diagnostics | figures, LOG | TBD |

Anchors: Winkel et al. (2025, arXiv:2506.19902) for NGC 4593; nested misaligned rings in
NGC 1068 (Vermot 2022, arXiv:2204.03621).
