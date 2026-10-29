# Log

Dated results and decisions, newest last. Task ids refer to [`PLAN.md`](PLAN.md). Runs are
in `results/m1torus/<campaign>/` (gitignored); figures in each run's `figs/`.

## 2026-10-27

- **A1.** Setup (Keplerian elements, live sink, barycentre at rest), live analysis in the
  paper's units, `run.sh`, `status.py`/`watch.sh`, `vism1.py`, `anim_particles.py`.
- **A2, cost.** 16k 39 s/orbit, 32k 104 s/orbit (~N^1.4); 64k 4.5 min/orbit measured.
- **A3, first 16k run** (`20261027-095659-canon-16k`, eps_c = 0.01). Mode appears at
  t ~ 35-40, saturates at A1/S0 0.34; after t ~ 100 energy rises (2.7% by t = 150),
  ~1% of particles leave, |r_c| falls from 0.025 to 0.008.
- **A2, diagnosis.** Restarts from t = 100, 5 orbits, one change each (dE/E per orbit):
  base -1.4e-4; tree_accuracy 0.3 -2.3e-4; C_force 0.1 -3.0e-4; C_cour 0.1 -2.8e-4;
  tolh 1e-6 -1.5e-4; global steps -7.3e-5; eps_c 0.05 -1.8e-4; fixed point mass
  -1.2e-5; eps_c = 0 -1.7e-5. Direct-summation energy matched phantom's to 1e-5: the
  drift was real, carried by the torus-central mass term. Cause: phantom softens the
  particle-sink force over max(eps_c, h_i) (`ptmass.F90`); with adaptive h_i (~0.2 near
  the centre) that potential is not conserved. **Decision: eps_c = 0** (default now).
  Two-body relaxation ruled out: ~4e4 orbits at 16k.
- **A3, 16k with eps_c = 0** (`20261027-141227-canon-16k-eps0`): mode holds to t = 150;
  t = 90-110: A1/S0 0.349, A2/S0 0.073, k 0.209, |r_c| 0.0217, |r_tb| 0.222; energy
  drift 1.1e-3 over 150 orbits; Rayleigh sigma 0.29 at t = 150; pattern speed 0.0186
  Omega_0 (period 54 orbits); profile in the frame of the mode fits a Poisson kernel
  with k = 0.200. Pericentres of the e >= 0.5 orbits sit at Phi_1 + 180 deg (apocentres
  build the overdensity) - open point A9. Q3 1.34, A3 near the 16k noise.
  Figures: `m1.png`, `ecc_maps.png`, `incl_peri.png`, `ecc_hist.png`, `pattern.png`,
  `corotating.png`, `anim_particles.mp4`, `anim_particles_corot.mp4`;
  comparison of both 16k runs: `results/m1torus/A/runs/compare_eps_16k.png`.
- **A9, settled on our side.** With varpi = Omega + omega as in the paper, the e >= 0.5
  pericentres sit at Phi_1 + 180 deg (within 5 deg) at t = 60-150, mean resultant length
  0.50-0.67; the projected azimuth agrees within 2 deg; omega alone is uniform
  (R <= 0.03), so the axis of Fig. 13 is varpi. This is the apocentre effect the paper
  cites in Sect. 4.2 (an e = 0.5 orbit spends 80% of its time in the apocentre half).
  "Near the phase of the pattern" fits either the 180 deg ambiguity of
  Phi_m = arctan(b_m/a_m) (a plain arctan flips our Phi_1 in 60% of the rows after
  t = 50) or a phase taken from the eccentricity vector. Question for the authors; A10
  states both angles.
- **Q3 at 16k is noise-limited.** Harmonic noise sigma/S0 = sqrt(2/N) = 0.011 and
  <A^2> = A^2 + 2 sigma^2, so a true A3/S0 = 0.016 reads 0.022; we measure 0.020.
- **Paper, Sect. 9.2 (NGC 613).** r_tb is scaled with M_tor (0.24 -> 0.72 R_tor), but
  Sect. 3.3 has r_c linear in M_tor with the same participating fraction, i.e. r_tb
  independent of M_tor. With r_tb = 0.24 the cavity-AGN distance (1 + M_tor) r_tb is
  2.5 pc against 3 pc observed, without dissipation. Self-similarity is tested only to
  M_tor = 0.1: task A11. **Sect. 9.1 (M31):** eq. 12 with dmu = -0.3 and S1/S2 = 7-10
  gives J = 9-13, which matches the quoted A~ = 0.80-0.86; the quoted J = 4-6 would give
  0.6-0.7.
- **Phantom.** Sinks accrete gas and dust only (`is_accretable` in `part.F90`); star
  particles pass through h_acc.
- **B1, decision.** Literature: thermal pressure alone cannot keep a torus thick
  (Krolik & Begelman 1988; Hopkins et al. 2012); simulations on these scales carry the
  unresolved support as an effective turbulent sound speed, 10-50 km/s observed
  (Hopkins & Quataert 2010), 0.07-0.35 v_K for NGC 613; beyond ~0.01-0.1 pc gravitational
  torques, not local viscosity, move angular momentum, and near the black hole the
  eccentric stellar pattern torques the gas through shocks (Hopkins & Quataert 2010,
  2011). Decided: gas from the A3 dump at t = 100 (f_gas 0.05, 0.1, 0.3) plus one run from
  t = 0 at the threshold; locally isothermal H/R = 0.1, then c_s = 0.3 v_K and one
  beta-cooling run at the threshold; no explicit viscosity, gas-only control for the
  numerical inflow. Defaults: gas self-gravity on, r_acc = 0.05, 16k stars plus ~32k
  lighter gas particles.
- **Not cited by the paper:** Hopkins & Quataert (2010a, b) and Hopkins et al. (2012), who
  found m = 1 eccentric gas+stellar discs near black holes and proposed the M31 disc as
  a fossil of such an inflow.

## 2026-10-28

- **A4, 64k** (`20261027-215121-canon-64k`, 150 orbits, 21:51-09:40, 4.7 min/orbit):
  saturated t = 110-150, against 16k over the same window and the paper: A1/S0 0.356
  (16k 0.344, paper 0.368), A2/S0 0.074 (0.072, 0.077), A3/S0 0.016 (0.020, 0.016), k
  0.207 (0.209, 0.21), Q3 1.05, noise-corrected 0.99 (16k 1.33 and 1.07, paper 0.99),
  |r_c| 0.0229 (0.0222, 0.024), |r_tb| 0.226 (0.218, 0.24), Rayleigh sigma 0.28,
  pattern 0.0190 Omega_0 (53 orbits), pericentres of e >= 0.5 at Phi_1 + 180 deg. Within
  6% of the paper and 3% of 16k: converged. Q3 = 1 at 64k confirms the 16k Q3 was noise.
  The mode grows at the same rate (e-fold 16 orbits at both N) but starts later at 64k
  (A1/S0 > 0.1 at t = 69 against 52): the seed is the shot noise, sqrt(pi/N). Energy
  drift 5.1e-4 over 150 orbits (16k 1.1e-3), dL/L ~1.4e-3, no particle lost.
  Figures in the run's `figs/`; comparison `results/m1torus/A/runs/compare_16k_64k.png`
  (`scripts/vis/compare_runs.py`).
- **Incident.** The run was extended from 110 to 150 orbits with `continue.sh`. phantom
  names the .ev after `logfile` and rewrote `m1torus01.ev` (and the sink .ev) from
  t = 110, losing 0-110. The energy was rebuilt from the 220 full dumps exactly as
  `energies.f90` sums it (kinetic + sink potential + poten), agreeing with phantom's
  etot to 7e-8 where both exist, into `m1torus00.ev`; the mode diagnostics
  (`m1torus_m1.ev`) were never affected. `continue.sh` now moves `logfile` to the next
  free number.
- **A4, reading corrected.** "Converged" was too strong: from 16k to 64k A1/S0 and |r_c|
  rose by ~3% towards the paper, so they still drift with N. The paper's own Fig. A.1
  shows its onset moving later with N too (256k and the symmetrised run start last),
  with the same saturated r_c (plateau ~0.023 on the figure). Our 64k starts ~20-30
  orbits after their 64k, i.e. grows slower. Suspect: phantom's adaptive softening,
  h ~ 0.08 against eps = 0.01, weakens small-scale self-gravity; their eps tests
  (1e-2-1e-4) are all far below our h. New tasks A12-A15, run before the controls.
- **A12, Fig. A.1 digitised** (`scripts/vis/digitise_figA1.py`). Smoothed r_c reaches
  0.005 / 0.015 / plateau at: paper 64k 39 / 61 / 0.0233, 128k 38 / 58 / 0.0241, 256k
  60 / 79 / 0.0239, symmetrised 104 / 125 / 0.0230; ours 16k 49 / 69 / 0.0222, 64k
  66 / 85 / 0.0229. The rise 0.005 -> 0.015 takes 19-22 orbits in every run, theirs and
  ours: the growth rate agrees, the earlier "grows slower" was wrong. The onset scatters
  by 65 orbits among their own runs (seed realisation) and ours lies inside it. Their
  plateau scatters 0.0230-0.0241; our 64k 0.0229 is at its lower edge, ~3% under their
  mean, the same size as their run-to-run spread.
- **A13.** phantom has no fixed softening for star or dark-matter particles: their h comes
  from the neighbour sum (`dens.F90`); only sinks take a fixed hsoft. A fixed eps would
  need a code change; hfact rescales h.
- **A5, thin control** (`20261028-121511-thin30-16k`, i_max = 30 deg, 16k, 150 orbits;
  interrupted at t = 148 and finished with `continue.sh`, whose new logfile kept the .ev).
  No mode: over t = 110-150 A1/S0 0.019, A2/S0 0.015, A3/S0 0.015, the flat spectrum of
  shot noise (floor sqrt(pi/N) = 0.014); |r_c| 0.0010, |r_tb| 0.009, the noise level of
  a barycentre of 16k particles. The paper's 0.5-30-0.3 gives A1 = 1.11, A1/S0 = 0.010,
  twice its 128k floor: the same verdict. `results/m1torus/A/runs/compare_thin30.png`.
- **A14, tests** (16k, eps_star = 0.01, h = 0.014 fixed, 3 orbits, branch `m1torus/fixedsoft`).
  At phantom's step factors (C_cour 0.3, C_force 0.25) and at 0.1 the small fixed
  softening does not resolve close passes: dE/E slope -4.2e-4 per orbit at 0.1. At 0.04
  the slope is 4.7e-6 per orbit (~7e-4 over 150 orbits, as A3's 1.1e-3), at ~18 orbits/h
  on 4 threads. The run uses C_cour = C_force = 0.04.
- **A14, fixed softening** (`20261028-150605-canon-16k-fixeps001`: A3's initial state,
  eps = 0.01 Plummer-equivalent fixed, h = 0.014, C_cour = C_force = 0.04, 150 orbits,
  6 h at 24 orbits/h; dE/E 1.8e-3). Over t = 110-150 against A3 (adaptive h ~ 0.08):
  A1/S0 0.337 (0.344), A2/S0 0.070 (0.072), k 0.206 (0.209), Q3 noise-corrected 1.03
  (1.07), |r_c| 0.0212 (0.0222), |r_tb| 0.212 (0.218), pattern 0.0183 Omega_0 (0.0186),
  Rayleigh sigma 0.30 (0.29); growth e-fold 17.8 (16.8) orbits, onset 3 orbits later.
  The paper's softening does not move the run towards the paper (0.368, 0.024): the
  adaptive h is not the cause of the 3-6% gap. What remains is the N trend (16k -> 64k
  +3%, extrapolating to ~0.36 at 128k) and the run-to-run spread the paper shows itself
  (r_c plateau 0.0230-0.0241). Comparison: `results/m1torus/A/runs/compare_softening.png`.
- **A6 and the paper's Fig. 14.** The weak case 0.5-45-0.3 grows slowly in the paper: A1
  rises from t ~ 50 and saturates (A1/S0 ~ 0.09) only at t ~ 500-1000, about 0.04-0.05
  at t = 150. A6 is at the noise floor at t = 81; it goes on to 1000 orbits after 150.

## 2026-10-29

- **A11 start** (`20261029-100846-mtor03-16k`, M_tor = 0.3; the paper ran 0.1, 0.06, 0.01
  only). The Keplerian start is far from equilibrium at this mass: m = 2 grows to A2/S0
  0.37 by t = 3.8 and decays by t = 8, the torus relaxes to <a> 1.14, <e> 0.36.
- **B, what the targets are** (papers in `materials/papers/`). Observed molecular tori:
  R 6-27 pc, M(H2) 0.7-3.9e7 Msun (Combes et al. 2019), thin where seen edge-on (NGC 1672
  axis ratio 6.5-10), clumpy, warped, slightly lopsided (Combes et al. 2025). NGC 613: ring
  R = 8 pc, M(H2) 1.0-1.4e7, M_BH 3.4-4.5e7, so M_gas/M_BH 0.2-0.4 before the X_CO
  uncertainty; bulge mass within 10 pc 2.7e6: the torus is gas, the stars a small
  spheroid. Hole 0.5 pc offset 3 pc from the AGN (0.375 R), two nascent spirals; Combes
  estimates a torus life of a few Myr (P = 0.29 Myr at 8 pc). Circinus: dust disk 5 x 3 pc
  fit by clumpy models with ~4 clouds along an equatorial ray (Lopez-Rodriguez et al.
  2026), so clouds collide every ~1/(2 pi N0 H/R) ~ 0.4 orbits at H/R = 0.1: collisional.
  The paper itself attributes its 5.8 pc against 3 pc for NGC 613 to cloud dissipation and
  defers it (Sect. 9.2). Gas analogue: m = 1 instabilities of discs with M_D/M_* ~ 0.1-1
  grow on the outer dynamical time and move the star off the centre of mass (Adams, Ruden
  & Shu 1989). Stellar nuclei (M31): gas comes from stellar mass loss of the eccentric
  disc (Chang et al. 2007). Alternative for central holes: AGN feedback (Garcia-Burillo
  et al. 2024). Design consequences not yet agreed.
- **B, checks of the numbers.** Combes et al. (2025) quote 250 km/s at 1 pc and i = 50 deg
  for NGC 613, which gives 2.5e7 Msun, not their 3.4e7 (M-sigma 3.7e7): M_BH 2.5-3.7e7.
  Their gas mass is 1.0e7 (Table 3), 1.4e7 (text), 0.97e7 within 0.2 arcsec (Table 2),
  all at the Galactic alpha_CO = 4.36; in galaxy centres alpha_CO is often several times
  lower, so M_gas/M_BH spans ~0.1-0.5. i = 50 deg comes from the oval assuming a circular
  torus; an eccentric torus would bias it. The AGN sits in emission (the HCO+ peak is at
  the core) and the hole of 0.5 pc lies 3 pc away, so the hole is a void on one side of
  a lopsided disc, not a central cavity around the black hole: the paper's reading of the
  offset as the torus barycentre is doubtful, and B compares maps, not barycentres.
  Moment 2 reaches 300 km/s on the major axis by beam smearing, so the true dispersion is
  unknown; a mock map smeared the same way is the fair comparison. Stability: a uniform
  Keplerian disc has Q ~ (H/R)(M_BH/M_gas); NGC 1672 (M_gas/M_BH 0.04, H/R 0.10-0.15)
  gives Q 2.3-3.4, NGC 613 at H/R 0.1 gives Q ~ 0.3, a fragmenting disc. Circinus (Izumi
  et al. 2023): M_BH 2.0e6 from HCN, dense gas h/r ~ 0.1, sigma ~ 15 km/s, Q < 1 beyond
  1 pc, Keplerian emission stronger on the blue side, inflow 0.2-0.34 Msun/yr at 0.27 pc,
  ten times the viscous estimate; Tristram et al. (2022): thin cold disc, dynamical mass
  ~3e6 within 2 pc. phantom offers `set_disc` with the enclosed-mass correction of the
  rotation, `ieos = 6` (locally isothermal about the sink; `ieos = 3` is about the origin,
  wrong once the black hole moves) and sink creation.
- **B3, decision** (with the user). Three ways to treat a massive gas torus whose thin
  state has Q < 1: (1) an effective turbulent sound speed set so that Q_min ~ 1.5, locally
  isothermal about the sink (Hopkins & Quataert 2010 carry unresolved support this way);
  H/R comes out ~0.15 at M_gas 0.1 and ~0.45 at 0.3; (2) a thin disc, H/R = 0.1, with sink
  creation where clumps collapse (closest to the clumpy Q < 1 disc of Circinus, but the
  clumps then move as collisionless points); (3) beta cooling (Gammie 2001), the disc
  regulates itself to Q ~ 1 and its thickness follows, at the price of a free beta,
  while real molecular gas cools fast enough to fragment. Chosen: (1) as the model,
  (2) as the check, (3) recorded only. Limits of (1): it cannot make clumps, it predicts
  a thick torus at M_gas 0.3 (sigma ~ 60 km/s at 8 pc), which the oval of NGC 613 can
  test, and it disagrees with the thin Q < 1 dense gas of Circinus, which (2) covers. Its
  pressure also makes a fluid m = 1 mode precess, so the outcome is not given in advance.
