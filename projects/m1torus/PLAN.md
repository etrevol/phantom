# Plan

The single tracked list of tasks of this project. Status: `[x]` done, `[~]` in progress,
`[ ]` open. Each task has an acceptance check and a cost on this laptop (8 cores).
Results and decisions, dated, go to [`LOG.md`](LOG.md); this file says what is left.

## Goal and principles (agreed with the user 2026-10-29)

Goal: a physical model of the molecular torus of NGC 613 that matches the ALMA data
well (A, the reproduction of Bannikova et al. 2026, is done and is the starting point),
and through it a model plus theory that predicts observations: predict a measurable
signature for a case and check it against what the observations of that case show.
Deliverables: a journal paper, a technical report (LaTeX), and the tools (the cloud
collision module in phantom, the synthetic-observation pipeline, the fit to ALMA data).
NGC 613 is approached in stages: cheap checks on existing data first, then a fit of
the best physically grounded model if the checks justify it.

Principles:
- Physics over fit: a physically grounded model with an honest mismatch beats a tuned
  one with a good match; parameters set by hand are named as such.
- Novelty through the literature: before every conclusion, check what is known, from
  downloaded papers; keep "ours" and "known" apart.
- Our own measurements of the data: compare with the raw archive data (ALMA) processed
  exactly like the models, rather than with published figures.
- Time economy: a stop criterion before every run, cheap checks on existing data first,
  no runs for their own sake.

Ground truth is the paper and our own measurements. The side projects on this laptop
(athena-torus, athena-vv, phantom-sgtorus) are reference only; of athena-collab only
`artem-develop` counts as checked.

Cost reference: 16k ~45 s/orbit (150 orbits ~2 h), 64k 4.7 min/orbit (150 orbits 12 h). Ask the user before anything over 30 min.

## Roadmap (agreed with the user 2026-10-30)

Where B stands: collisional cloud tori cannot stay hot for long whatever heats them
(known since KB88 for damping; our three heating laws quantify it; a speed-dependent
restitution gives a stable balance with sigma ~ v_A); the collisionless torus fits the
archive cube best once its size is free (R_tor ~6.6 pc, i ~60, M_dyn ~1.4e7) but lacks
the central component (40% of the CO within 4 pc, round core) and drives e_max to 1.
Order (the user's choice 2026-10-30, after B19-B26: consolidate before predicting;
no 64k run for now):
1. B27: the black-hole mass, ours 2-2.5 times below Combes et al. (2025) in both
   galaxies (1.5 vs 3.4e7, 1.9 vs 4.5e7): their PV method on our models, PV metrics in
   the fit.
2. B28: extended mass in the potential (stars, gas) where the rotation rises outward
   (NGC 1672, 107 -> 127 km/s at 1.6-6 pc); it also sets the precession of the mode.
3. B29: the m = 1 tension in NGC 613: m0 A1 0.17 against A2, A3 ~0.10 within 0.1
   arcsec and a torus centroid ~1 pc from the continuum AGN, closer to the strong mode
   (A3: r_tb ~0.24 R) than to the weak one of the fitted torus (B25: 0.17 pc); torus
   mass, orbit parameters, or another cause (feeding, clumps, outflow).
4. B20 only then, blind: predictions written before looking, for quantities not yet
   seen (HCO+(4-3), other galaxies with high-resolution archive data).
5. Write-up after the model: A10 report with the B chapters; papers decided later.

## A. Reproduction of the collisionless mode

| id | status | task | acceptance | cost |
|---|---|---|---|---|
| A1 | [x] | setup, live analysis, run tooling, figures | runs end to end | - |
| A2 | [x] | numerics at 16k: cost, conservation, energy-drift diagnosis | dE/E < 2e-3 over 150 orbits | - |
| A3 | [x] | canonical `0.5-60-0.3`, 16k, 150 orbits | mode appears and holds | 2 h |
| A4 | [x] | canonical at 64k, 150 orbits | saturated A1/S0, k, r_c within 10% of the paper; 16k -> 64k converged | 12 h |
| A5 | [x] | thin control `0.5-30-0.3`, 16k, 150 orbits | no m = 1 above noise (paper Sect. 5) | 2 h |
| A6 | [x] | weak case `0.5-45-0.3`, 16k, to 400 orbits | rises clearly above noise by t ~ 300-400 as in Fig. 14 (A1/S0 ~0.07-0.09); 1000 only if 400 leaves doubt | done: no mode, A1 ~ A2 ~ A3 at the noise level (LOG 2026-10-29) |
| A7 | [x] | torus mass M_tor = 0.06 (needed) and 0.01 (optional: 0.06, 0.1, 0.3 already span a factor 5), 16k | r_c linear in M_tor with r_tb fixed (Fig. 6, Sect. 3.3); pattern speed linear in M_tor (Sect. 10); growth slows with mass, check Fig. 6 before setting the length | 0.06: ~250 orbits, 3 h; 0.01: >= 1000 orbits, 13-17 h |
| A8 | [ ] | `0.2-60-0.3`, 16k (optional, no new test) | stronger mode, A1/S0 ~0.46 | 2 h |
| A9 | [x] | pericentres of e >= 0.5 sit at Phi_1 + 180 deg here; the paper says "near the phase": settle the convention | explanation in LOG | - |
| A10 | [ ] | short write-up of A: methods, differences from the paper, results | LaTeX, like `projects/ppi/formulary` | - |
| A11 | [x] | M_tor = 0.3 (the NGC 613 mass ratio), 16k, 120 orbits | does r_tb stay ~0.24 (Sect. 3.3) or grow (Sect. 9.2); cavity-AGN offset (1 + M_tor) r_tb against the observed 3/8 = 0.375 R_tor | 2-2.5 h |
| A12 | [x] | digitise the paper's Fig. A.1: onset and growth rate of r_c for 64k, 128k, 256k | their growth rate and onset against ours (64k starts ~20-30 orbits later here) | 30 min |
| A13 | [x] | can phantom soften star particles with a fixed length instead of the adaptive h | answer from the code | 30 min |
| A14 | [x] | fixed softening at 16k (setup `eps_star`, merged 2026-10-28): A3 again with the same seed and eps = 0.01 as the paper, after 16k tests of cost and energy | A1/S0, r_c, onset and growth against A3 (adaptive h ~ 0.08) and the paper | tests 30 min, run 2-15 h |
| A15 | [ ] | 16k canonical run continued to 1000 orbits (A3 from t = 150) | the mode lives to t = 1000 (paper Sect. 3.1) | ~14 h, overnight |
| A16 | [ ] | (lowest) oscillation of the central-mass radius against softening: an unpublished figure of the authors shows |r_c| at 64k oscillating with a ~55-orbit period and growing amplitude after t ~ 200 for eps = 1e-3, 1e-4, flat for 1e-2 | from A15 (our sink unsoftened, torus h ~ 0.08) and, if needed, fixed eps_star 1e-3 and 1e-2 runs (C = 0.04, ~24 orbits/h at 16k): which softening, sink-particle or particle-particle, drives it | after A15 |

Order (2026-10-29, agreed with the user): by day A11, A7 at 0.06, A6 on to 400; then A10
(LaTeX write-up of A, edited later if needed); A15 last of the needed tasks, overnight.
Optional, only if results call for them: A7 at 0.01, A8, A6 to 1000; A16 last of all.
All at 16k: 64k moved the saturated values by ~3% and the signals sought are 3-25 times
the 16k noise floor.

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

Saturated averages over t = 110-150 for both runs (the 16k numbers in LOG for t = 90-110 are
within 1%).

| quantity | paper | 16k, eps_c = 0 (A3) | 64k (A4) | where |
|---|---|---|---|---|
| onset of m = 1 | after virial equilibrium; axisymmetric t = 35, lopsided t = 51 | A1/S0 > 0.1 at t = 52 | A1/S0 > 0.1 at t = 69; same growth, e-fold 16 orbits | Fig. 4, 9 |
| A_1/Sigma_0 | 0.368 (A_1 = 40.2) | 0.344 | 0.356 | Table 2 |
| A_2/Sigma_0, A_3/Sigma_0 | 0.077, 0.016 | 0.072, 0.020 | 0.074, 0.016 | Table 2 |
| k, Q_3 | 0.21, 0.99 | 0.209, 1.33 (noise-corrected 1.07) | 0.207 (Poisson fit 0.204), 1.05 (0.99) | Table 2, Sect. 7.2 |
| r_c | ~0.024 | 0.0222 | 0.0229 | Sect. 9.3, Fig. 3 |
| r_tb | 0.24 | 0.218 | 0.226 | Sect. 9.2 |
| eccentricities | Rayleigh, sigma = 0.27 | sigma = 0.28 | sigma = 0.28 | App. B.1 |
| pattern | prograde, period >> orbital | 0.0186 Omega_0, period 54 orbits | 0.0190 Omega_0, 53 orbits | Sect. 3.1 |
| energy drift, 150 orbits | not quoted | 1.1e-3 | 5.1e-4 | - |

### Differences from the paper (measured, not assumed)

- Central mass unsoftened (eps_c = 0): phantom softens the particle-sink force over
  max(eps_c, h_i) with the adaptive h_i, which does not conserve energy (LOG 2026-10-27).
  The paper's eps = 0.01 is negligible against pericentres >= 0.05.
- Gravity by tree, adaptive softening between torus particles (h ~ 0.08), leapfrog with
  individual steps; the paper: direct summation, eps = 0.01, Hermite.

## B. Dissipation: the molecular torus as colliding clouds (after A)

Question: the paper applies its collisionless torus to the molecular torus of NGC 613
(Sect. 9.2) and blames its 5.8 pc against the observed 3 pc on dissipation between clouds,
left for future work. Does a torus of colliding clouds form and keep an m = 1 mode, does
the black hole still move, and do synthetic maps look like the ALMA maps? Targets and
their numbers: LOG 2026-10-29. Main target NGC 613 (Combes et al. 2025), second Circinus
(Izumi et al. 2023; Lopez-Rodriguez et al. 2026; Tristram et al. 2022); M31 optional.

Models (decision 2026-10-29, LOG):
- **Main: sticky clouds.** Each particle is one cloud: the paper's self-gravitating
  collisionless particles, softened by the cloud radius r_cl (fixed `hsoft_nbody`), plus
  partially inelastic collisions when two clouds overlap and approach. The impulse acts
  along the line of centres with restitution e, so momentum and angular momentum are exact
  and the lost energy is booked as radiated. Constant cloud masses (no coalescence or
  fragmentation, unlike Vollmer et al. 2022). r_cl follows from the covering number N0 of
  an equatorial ray, N0 = n pi r_cl^2 L: N0 ~ 4 (Circinus clumpy fit) with 16k clouds in a
  disc of H/R 0.2 gives r_cl ~ 0.01 R (0.08 pc at NGC 613), under the tidal radius
  R (m_cl / 3 M_BH)^(1/3) ~ 0.02 R; clouds of ~600 Msun and that size have n(H2) ~ 5e6
  cm^-3, the density of the HCO+(4-3) gas. New code: phantom has no collisions.
- **Second: fluid limit.** SPH gas, the limit of frequent collisions: effective turbulent
  sound speed with Q_min ~ 1.5, locally isothermal about the sink (`ieos = 6`), as
  Hopkins & Quataert (2010). No explicit viscosity, shock viscosity with its switch, gas
  self-gravity on, the sink accretes gas only, r_acc = 0.05.
- Recorded, later if the two disagree or call for it: thin isothermal disc with sink
  creation (Q < 1), beta cooling, ISM cooling with star formation and HII feedback.
  Agreement of the two limits makes a result robust; disagreement is itself the result.

Units and observables (code: G = M_BH = R = 1):

| | NGC 613 | Circinus |
|---|---|---|
| M_BH [Msun] | 3.4e7 (range 2.5-3.7e7) | 1.7e6 maser, 2.0e6 HCN |
| R = 1 is | 8 pc, the torus radius | 2.5 pc, half the 5 x 3 pc dust disc |
| v unit, P(R = 1) | 135 km/s, 0.36 Myr | 54 km/s, 0.28 Myr |
| M_gas / M_BH | 0.27-0.56 at Galactic alpha_CO, ~0.1 if alpha_CO is 3-5 times lower, as often in centres | < 0.5 (dynamical mass 3e6 within 2 pc) |
| thickness | unknown (i = 50 deg from the oval, assuming a circular torus) | dense gas h/r < 0.1-0.3, sigma ~ 15 km/s |
| stability | Q ~ (H/R)(M_BH/M_gas): 0.3 at H/R 0.1, M 0.3 | Q < 1 at r > 1 pc (Izumi) |
| asymmetry | hole r = 0.5 pc at 3 pc NW of the AGN (0.375 R); AGN in emission; two spirals | Keplerian side stronger on the blue side |
| beam | 0.015 arcsec = 1.2 pc = 0.15 R | ~2 pc = 0.8 R (ALMA 700 um); AMI 1.9 pc |
| other | PV on PA -20 deg; flux 15.4 Jy km/s for 9.7e6 Msun; stars within 10 pc 2.7e6 | inflow 0.2-0.34 Msun/yr at 0.27 pc = 0.005-0.009 in code units |

| id | status | task | acceptance | cost |
|---|---|---|---|---|
| B1 | [x] | first design with the user (2026-10-27) | LOG | - |
| B2 | [x] | targets and their numbers (2026-10-29) | LOG; papers in `materials/papers/` | - |
| B3 | [x] | models (2026-10-29): sticky clouds the main, SPH fluid limit the second, three more recorded | LOG | - |
| B4 | [x] | collision module: pairs from the tree within r_i + r_j, approaching only, normal restitution e, dissipated energy summed into the energy budget; hook after the step, global steps first with dt < C_coll r_cl / v_rel; `.in` keys for r_cl and e, off by default | two-cloud collision against the analytic result; p and L to round-off, E + E_diss to the integrator error; collision rate in a uniform cloud box against n sigma v_rel; e = 1 conserves E; r_cl = 0 reproduces A3 bit for bit | 2-4 days |
| B5 | [~] | sticky tests: cost per orbit at 16k with global steps; dt convergence; N0 = 4 and 10 | cost known; result unchanged when C_coll halves | 2-4 h |
| B6 | [~] | sticky runs from the paper's torus: M_tor 0.1 (pair of A3) and 0.3 (pair of A11) at e = 0.5; one at e = 0; one from A3 at t = 100 (collisions switched on with the mode in place) | mode growth and saturation, r_c, thickness and dispersion against A; dissipated power | 4-5 x TBD |
| B7 | [~] | fluid limit, code: gas option of `setup_m1torus` through `set_disc` (enclosed-mass rotation, `ieos = 6`, c_s from Q_min); moddump turning cloud particles into split gas particles for the same initial torus | a disc at Q >= 2 stays axisymmetric 5 orbits; gas-only control gives the numerical inflow | 1 day + 2 h tests |
| B8 | [ ] | fluid runs: the paper's torus as gas (M_tor 0.1, 0.3) and the NGC 613 disc (r 0.1-1.2, M_gas 0.1, 0.3) | same diagnostics as B6 | 4 x TBD |
| B9 | [x] | synthetic observations: `mockobs.py` (units, projection at i and PA, cube with the beam, moments 0-2, PV, integrated spectrum, hole offset, A1 of moment 0; clouds emit with their mass spread over r_cl); `ref_figures.py` (crops of the published panels with calibrated axes, kept local); digitised Combes Fig. 6 moment 0 | a mock of an axisymmetric Keplerian disc gives the Keplerian PV and no hole offset | 1-2 days |
| B10 | [ ] | comparison figures, ours left and the published panel right, same field, beam, threshold and colour scale, profile strip below where digitised; list below | figures, LOG, handbook | 1 day |
| B11 | [ ] | optional, M31: a few per cent of clouds in the A3 torus at t = 100 (gas of stellar mass loss, Chang et al. 2007) | does the stellar mode notice them | 2-3 h |
| B12 | [x] | heated colliding-cloud model (thermostat s = sigma/v_K, clump sinks): s = 0.1 (done, old thermostat holding sigma_R,phi,z) and 0.3 (thermostat on sigma_R, z after the old one unbound the torus; drains and evaporates in ~5 orbits), from A3 at t = 100, 20 orbits each; compare with B6b (no heating) and NGC 613 moment 2 | does the mode survive at the observed dispersion; moment 2 and thickness against Combes 2025 and Izumi 2023 | 2 x ~2.5 h |
| B13 | [x] | is NGC 613 dynamically hot? (a) the published moment 2 against our synthetic maps made with their stated recipe (plain clip 0.6 mJy/beam over +-300 km/s) and with smooth-and-mask; (b) the ALMA archive cube of project 2019.1.00273.S (CO(3-2), 0.015 arcsec): our own moments with a known recipe and Gaussian fits per pixel | moment 2 / moment 1 of NGC 613 from data we measured ourselves; which recipe made Fig. 6 | (a) 1-2 h; (b) download size decides |
| B14 | [x] | fit of the viewing geometry: inclination 20-60 deg and M_BH 1.5-4e7 Msun against |moment 1|, moment 2, moment-0 axis ratio and hole offset, on the existing snapshots (no new runs) | what part of the gap is i and M_BH | 2-3 h; also on A5, A6, A7, A11: collisionless tori fit best (M_tor 0.3: chi2 12) |
| B15 | [x] | if NGC 613 is hot: heated clouds at s 0.3 and 0.5 (mode survival at sigma/v ~ 0.5); clouds with rarer collisions, r_cloud for N0 ~ 0.5 and 1 (the paper's collisionless limit approached from the dissipative side) | survival of the mode and moment 2 / moment 1 against B13 | done 2026-10-29: s 0.3 drains (N0 2.6) or evaporates (N0 0.5); a hot torus that lasts needs N0 <~ 0.04 (LOG) |
| B16 | [ ] | only if B15 fails: non-circular motions from feeding (inflow along the nuclear spiral) and outflow, a torus fed from outside | moment 2 and PV against Combes 2025 | days (overlaps C) |
| B17 | [x] | literature on the dispersion of cloud tori: collisions damp it in orbits (Krolik & Begelman 1988), clumpy-torus covering numbers (Nenkova et al. 2008), tori as disc winds (Elitzur & Shlosman 2006), heated cloud models (Vollmer et al.), the molecular outflow of NGC 613 (Audibert et al. 2019; Combes et al. 2019, 2025) | what of B12-B15 is known and what is new; does an outflow take part of moment 2 | 1-2 h |
| B18 | [x] | the ALMA data in depth: continuum subtraction of the pipeline cube against Combes Fig. 6 and 12 (the stripe at the AGN, their hotter central moment 2), the continuum image, total flux against their Table 2, HCO+(4-3) | why the archive cube and the published figures differ; which to trust for kinematics | 1-2 h |
| B19 | [x] | fit of the collisionless torus to the archive cube on initial states (phantomsetup, seconds each): i_max, e_max, M_tor, R_0 with i and M_BH; one full run of the best (ask first) for the mode and the evolved maps | best physically grounded model and its mismatch | done 2026-10-30: torus region, R_tor 5.9 pc, i 48, M_dyn 1.5e7, chi2 1.8/dof |
| B20 | [ ] | prediction and check: a measurable signature of the model (e.g. black-hole offset from the torus centre, lopsided moment 0, PV asymmetry) for a case with observations, predicted before looking | the prediction against the data | open |
| B21 | [x] | heating by partially elastic collisions in the shear (Krolik & Begelman 1988, Sect. IVb): clouds with e_rest 0.9, no thermostat, from A3 at t = 100; extreme first, stop at t = 103 if the torus flattens as B6b; then 0.7-0.8 for the threshold and 20 orbits with the mode and synthetic maps | does a cloud torus heat itself to sigma ~ 0.3-0.5 v_K, and the e threshold against KB88 eq. 24 | done: constant e has no balance (0.5 collapses, 0.9 runs away); next a speed-dependent e (KB88 eq. 25), ask the user |
| B22 | [x] | the molecular outflow of NGC 613 (Audibert et al. 2019) added kinematically to the synthetic cube: how much of the nuclear moment 2 and the PV wing it explains | moment 2 at the nucleus and PV against the archive cube | done 2026-10-30: the outflow does not set the annulus moment 2; feeding needs the drained gas to leave |
| B23 | [x] | magnetised clouds: restitution from the collision speed (KB88 eq. 25, `v_alfven`), no thermostat, from A3 at t = 100; v_A 0.2 first, stop criterion at t = 103 | a stable hot torus (sigma ~0.25-0.4 v_K predicted), then the mode and the synthetic maps against the archive cube | done 2026-10-30: v_A 0.2 bound but cold (sigma_z 0.11), v_A 0.6 at sigma_z 0.28 drains and evaporates |
| B24 | [x] | torus plus central disc fitted together to the archive cube (`fit_torus.py` with a disc component: radius, mass share, thickness) | the central moment-0 share, central moment 2 and the round core without near-radial orbits | done 2026-10-30: the disc does not improve the fit per degree of freedom; dropped |
| B25 | [x] | phantom run of the best physically grounded model, stop criterion, three-column maps | does the mode grow in an NGC 613-like torus; evolved maps against the archive cube | done 2026-10-30: weak slow mode (A1/S0 0.05-0.06, locked at t ~ 78, still rising); evolved torus chi2 25.6 |
| B26 | [x] | other tori of Combes et al. (2025), archive cubes of 2019.1.00273.S, same processing and fit | the conclusion beyond one galaxy | per galaxy ~1-2 h |
| B27 | [ ] | black-hole mass: our fits 2-2.5 x below Combes et al. (2025) in both galaxies; their PV method (peak projected velocity at ~1 pc) applied to our models, PV metrics added to the fit | the source of the factor and a consistent M_BH | ~2 h |
| B28 | [ ] | extended mass (stellar cluster, gas) in the generator potential and in phantom (external potential) where the rotation rises outward (NGC 1672) | rotation profile shape; its effect on the precession and the mode | ~2-3 h |
| B29 | [ ] | NGC 613 m = 1 tension: m0 A1 0.17 > A2, A3 ~0.10, torus centroid ~1 pc from the AGN; which torus (M_tor, e_max, i_max) or which other cause gives it | a model of the observed lopsidedness or a reason it is not the slow mode | ~2-4 h |

Order (2026-10-29, the user's choice): B13b first (the only independent check of the
number the conclusion rests on), B13a and B14 alongside, s 0.3 after A6.


Order: B4 and B9 code by day while A runs; B5 tests when cores are free; B6 before B8.
Runs over 30 min are asked for first.

B10 comparisons, for both models. Viewing angles are the observed ones; the observer azimuth relative to
the mode phase is the one free parameter, fitted to the hole position and stated. Each
quantity is shown over time, not only at the best-matching snapshot.

| # | published panel | ours | measured |
|---|---|---|---|
| 1 | Combes 2025 Fig. 6 left: CO(3-2) moment 0, beam 0.015 arcsec, field +-0.3 arcsec | moment 0 at i = 50 deg, PA -20 deg | emission radius (8 pc), oval axis ratio (i against thickness), peak and hole positions from the AGN (hole 0.5 pc at 3 pc), A1(R) of the deprojected map against the digitised one |
| 2 | Fig. 6 middle: moment 1, +-100 km/s | moment 1 | kinematic PA, amplitude, twist of the isovelocity lines (eccentric orbits) |
| 3 | Fig. 6 right: moment 2, 50-300 km/s | moment 2 smeared by the same beam | dispersion in the outer disc, where smearing is small: tests H/R |
| 4 | Fig. 12: PV along PA -20 deg with the Keplerian curve of 3.4e7 Msun | PV along the same slit | envelope, asymmetry between the two sides |
| 5 | Fig. 13 / Table 2: integrated CO(3-2) spectrum, FWHM 144 km/s | integrated spectrum | FWHM, horn asymmetry (lopsidedness) |
| 6 | Lopez-Rodriguez 2026 Fig. 2 (700 um dust, HCN disc) and Fig. 3e (dust of the clumpy-torus fit) | column map edge-on at PA 50 deg, their beam | extent and aspect ratio; coarse, 1-2 beams across |
| 7 | Izumi 2023 Fig. 2C: HCN and [CI] PV along the major axis | PV | Keplerian envelope, blue/red asymmetry |
| 8 | Izumi 2023 Fig. 2E (sigma/V_rot, a proxy of h/r) and Fig. 4 (Toomre Q) | the same profiles | h/r(r), Q(r) |
| 9 | Izumi 2023: inflow 0.2-0.34 Msun/yr at 0.27 pc | sink accretion rate | in physical units |

Not in B: the AMI images at 3.8-4.8 um (Lopez-Rodriguez Fig. 1) trace AGN-heated dust and
need radiative transfer (RADMC-3D or HyperCAT); AGN feedback, which may make the hole
itself (Garcia-Burillo et al. 2024); feeding from the nuclear spiral. Synthetic maps take
line emission proportional to column (optically thin, uniform excitation); CO(3-2) is
partly thick, so moment 0 is compared in shape, not in peak value.

Literature (in `materials/papers/`): Combes et al. 2019, 2025; Izumi et al. 2023;
Tristram et al. 2022; Lopez-Rodriguez et al. 2026; Gamez Rosas et al. 2025;
Garcia-Burillo et al. 2024; Pereira-Santaella et al. 2026; Vollmer et al. 2008, 2022;
Adams, Ruden & Shu 1989; Chang et al. 2007; Hopkins & Quataert 2010, 2011; Hopkins et
al. 2012; Wada 2002, 2012; Krolik & Begelman 1988; Gammie 2001; Kratter & Lodato 2016;
Lodato & Price 2010; Hönig 2019. Viscosity options in phantom: local note.

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
