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
- **B3 revised: sticky clouds** (with the user, same day). Hopkins & Quataert (2010) already
  ran GADGET with an effective sound speed of 20-100 km/s and star formation around a
  3e7 Msun black hole and found lopsided eccentric gas+star discs whenever more than
  ~0.1 M_BH of gas reaches 1-10 pc; the fluid model alone repeats their method. Vollmer et
  al. (2022) ran partially inelastic cloud collisions for the Galactic Centre and NGC 1068,
  but in a fixed potential, without self-gravity or a live black hole, so with no m = 1
  mode and no black-hole offset. Chosen: self-gravitating sticky clouds (the paper's
  particles plus the dissipation it invokes) as the main model, the SPH fluid limit as the
  second; thin disc with sinks, beta cooling and ISM cooling with star formation are
  recorded for later. phantom has ISM cooling (Koyama & Inutsuka), sink star formation and
  HII feedback, but supernova injection only as a test problem. Estimate for the clouds:
  N0 ~ 4 with 16k clouds in a disc of H/R 0.2 gives r_cl ~ 0.01 R, below the tidal radius
  ~0.02 R; ~600 Msun in 0.08 pc is n(H2) ~ 5e6 cm^-3, as the HCO+(4-3) gas.
- **A11, M_tor = 0.3** (`20261029-100846-mtor03-16k`, 16k, 120 orbits, 2 h 10 min). After
  the m = 2 transient the torus relaxes to <a> 1.14, <e> 0.36; m = 1 rises from t ~ 40 to
  A1/S0 0.11 by t = 60 and creeps to 0.12-0.14 by t = 120, not clearly saturated. Over
  t = 80-120 against A3 (M_tor 0.1) in the same window: A1/S0 0.132 (0.340), k 0.15
  (0.21), |r_c| 0.0210 (0.0213), |r_tb| 0.070 (0.216), M_tor r_tb = r_c as momentum
  requires; pattern 0.050 Omega_0 (period 20 orbits, A3 0.019); Rayleigh sigma 0.33;
  e >= 0.5 pericentres at Phi_1 - 176 deg as in A3. The black-hole offset stays at
  ~0.021, the barycentre offset drops threefold: the paper's linear scaling r_tb ~ M_tor
  (0.72 R at 0.3, Sect. 9.2) fails here, and (1 + M_tor) r_tb = 0.09 R = 0.7 pc at NGC 613
  against the observed 3 pc. Caveats: one seed, the start is far from equilibrium at this
  mass, the mode may still grow. dE/E 3.0e-3 and dL/L 6.1e-3 grow linearly, three and
  fifteen times A3; cause not yet checked. `results/m1torus/A/runs/compare_mtor.png`.
- **A7 start** (`20261029-122207-mtor006-16k`, M_tor = 0.06, 16k, 250 orbits; the paper's
  Fig. 6 runs to t = 1400, ~220 orbits). Tests the self-similar claim r_c ~ M_tor with
  r_tb fixed from below, after A11 broke it from above.
- **B4, collision module** (branch `m1torus/sticky`, `src/main/cloud_collisions.f90`). Pairs
  from a sorted grid of 2 r_cloud cells, impulse along the line of centres with
  restitution e_rest, lost energy in `ediss_coll` and the dump header; global timesteps,
  dt < C_coll r_cloud / v_max; sinks accrete clouds when collisions are on. Unit tests
  (`scripts/test/run_test_cloudcollisions.sh`): head-on and oblique pairs exact to 1e-14 in
  v, p, L and E + E_diss; collision rate of 20k elastic clouds in a box with walls
  1.00, 1.03, 1.00 of n sigma <v_rel> (three seeds; free streaming without walls gave
  0.89, a thinning box). r_cloud = 0 reproduces master bit for bit. In phantom (4k clouds,
  r_cloud 0.02, sink softened 0.05): e = 1 keeps E to 4e-5 over 1.5 orbits and L to 2e-6;
  e = 0.5 closes E + E_diss to 1.2e-4 until clouds reach the sink, whose mergers remove
  energy phantom does not tally; L stays to 1e-5 with 352 clouds accreted. Two traps:
  without accretion a plunging cloud forced 1e-9 substeps (now sinks swallow clouds), and
  a central clump deepens phantom's tree beyond 2 maxp nodes at small N: build with
  NCELLSMAX=1048576.
- **B4 merged, B5 first numbers.** Branch merged into master after a second regression
  (2 orbits at 2k: `.ev`, `m1.ev` and every physical array of the last dump identical;
  only the unused `alpha` of star particles differs, by 5e-44). Covering number of the
  paper's torus at 16k: r_cloud 0.01 gives N0 = 2.6 along an equatorial ray at t = 0
  (range 0-7), 1.5 at t = 100; r_cloud 0.005 gives 0.6. Cost: 16k clouds, r_cloud 0.01,
  r_acc 0.05, global steps, 9.2 orbits/h on 6 threads beside A7. First physics: the
  paper's torus (M_tor 0.1, e = 0.5) radiates 0.044 = 78% of |E0| in half an orbit and
  2% of the clouds fall into the black hole; a thick torus of colliding clouds collapses
  within about one orbit. `mockobs.py` maps A3 at t = 100 onto NGC 613 as a ring of
  radius ~0.1 arcsec with a moment-0 peak 0.23 Jy/beam km/s (observed scale 0.05-0.25),
  but its PV diagram lacks the +-250 km/s near the centre: the paper's torus has an
  empty middle, NGC 613 has gas down to ~1 pc.
- **B9, synthetic NGC 613 maps of the paper's torus** (A3 at t = 100, `compare_obs.py`,
  `results/m1torus/B/compare_A3/`, against Combes et al. 2025 Figs. 6, 12, 13 cropped by
  `ref_figures.py`). The ring of radius ~0.1 arcsec and its moment-0 level match; the
  rest does not: the observed PV keeps emission within 0.05 arcsec of the AGN and rises
  towards +-250 km/s there, the paper's torus puts it at +-0.1 arcsec on the Keplerian
  curve with nothing inside; moment 2 of the outer torus is 80-120 km/s against the
  observed 40-80 (a hotter, thicker torus than NGC 613's, beam smearing alike); the
  integrated line is 195 km/s wide against 144. The observed moment 1 recedes to the
  north (PA -20 deg), which the mock now follows.
- **B6 start** (`results/m1torus/B/runs/20261029-133711-clouds-m01-e05-16k`): the paper's
  canonical torus as 16k colliding clouds, M_tor 0.1, r_cloud 0.01 (N0 ~ 2.6), e 0.5,
  r_acc 0.05, global steps, 30 orbits on 6 threads beside A7.
- **B7, fluid limit code.** `moddump_m1torus` mode = gas turns the torus particles into SPH
  gas, ieos = 6 about the sink with c_s = hor v_K (q = 1/2), r_acc set; polyk, qfacdisc
  and isink travel in the dump header. Test: the paper's torus at 4k as gas with H/R
  0.15 runs 2 orbits in 106 s on 2 threads; 12% of the gas falls into the sink in that
  time, the shocked torus settles. Queue from here: one A slot (A7, then A6 to 400) and
  one B slot (B6 clouds, then B8 fluid runs), two runs at a time.
- **B7, gas disc** (`setup_m1torus` geometry = 1, `set_disc`). At M_gas 0.1 and H/R 0.15
  set_disc reports Q_min = 2.13 (the uniform-disc estimate gave 1.5); Q_min 1.5 needs
  H/R ~ 0.105 at M 0.1 and ~0.32 at M 0.3. Test, 4k, r 0.1-1.2, Sigma ~ 1/R, 3 orbits:
  axisymmetric (A1, A2 <= 0.03), |z|/R 0.10, but 19% of the gas reaches the sink at
  r_acc 0.05: inflow from the inner edge at low resolution; the 16k gas-only control
  (B7 acceptance) will measure it before the B8 disc runs.
- **B6 after one orbit**: the cloud torus has collapsed to a thin disc, median |z|/R
  0.27 -> 0.056, sigma_z 0.44 -> 0.15, median R 0.84 -> 0.60, E_diss 0.075 (1.3 |E0|),
  4% of the clouds in the sink. Its synthetic PV on NGC 613 now follows the Keplerian
  curve of 3.4e7 Msun from 0.03 to 0.08 arcsec at 100-200 km/s, brighter and faster
  than the observed PV, whose emission at +-0.05 arcsec stays within ~100 km/s (lower
  than Kepler for 3.4e7 at i = 50; with the 2.5e7 that Combes' own numbers give, or a
  smaller inclination, the gap shrinks).
- **B9 done**: `mockobs.py`, `ref_figures.py` (Combes 2025 Figs. 6, 12, 13; Lopez-Rodriguez
  2026 Figs. 2f, 2l, 3d, 3e; Izumi 2023 Figs. 2C, 2E), `compare_obs.py` (NGC 613),
  `compare_circinus.py` (dust column and axis ratio, PV, sigma_z/v_phi, accretion in
  Msun/yr). On B6 at t = 3 orbits scaled to Circinus: axis ratio 0.50 of the dust map
  (their disc 0.6), accretion 0.017 Msun/yr against the measured 0.2-0.34; the torus
  mass 0.1 M_BH = 1.7e5 Msun is an assumption, not a Circinus measurement.
- **Hole measure reworked** (`mockobs.hole`): the lowest ratio of the smoothed moment 0 to
  its local maximum over 4 beams, inside the emission. The paper's torus (A3, t = 100)
  on NGC 613: hole 2.5 pc from the AGN (observed 3 pc). B6 after its collapse: the hole
  sits on the sink (0.3-0.5 pc, one pixel), the cavity of the accretion radius, with no
  offset. A7 at t = 148: A1/S0 0.345, |r_c| 0.013 = (0.06/0.1) x 0.022: the paper's
  self-similar scaling holds from below, while A11 breaks it from above.
- **moddump mode = clouds**: switches collisions on in an evolved torus (r_cloud, e_rest,
  r_acc, softening 1.4 r_cloud). Checked on A3 at t = 100 (0.1 orbit, `B/tests/`): the
  new .in carries r_cloud, e_rest, C_coll and hsoft_nbody; the mode is in place
  (A1/S0 0.35). Next in the B slot after B6: this restart for 20 orbits (does the mode
  survive collisions?), then the fluid runs B8.
- **Snapshot 14:46.** A7 (M_tor 0.06, t = 176): onset t = 83 against 52 at 0.1, e-fold 28.5
  orbits against 16.8 (ratio 1.70, the mass ratio is 1.67: growth rate ~ M_tor),
  saturating at A1/S0 0.35 like 0.1, |r_c| 0.0134 = 0.6 x 0.022 (`compare_mass.png`).
  B6 (t = 15): |z|/R 0.015 (the cloud radius), sigma_z/v_phi 0.05, 24% of the mass in
  the sink, E_diss 4.4 |E0|, no m = 1 (A1/S0 at the noise floor), median R back to 0.86
  as the disc spreads. On NGC 613 it is a filled disc with a centred hole of the
  accretion radius; moment 2 (40-100 km/s, high at the centre) now resembles the
  observed map, moment 1 reaches +-100 km/s against the observed ~+-50.
- **B6 ends at t = 28.3 by fragmentation** (tree overflow even at NCELLSMAX = 2^20). After
  the collapse the cold disc of clouds (sigma_R ~ 0.04, Q_stellar 1.1-1.5 at R 0.5-1.2)
  grows multi-armed spirals and from t ~ 24 fragments: at t = 28, 357 clouds within one
  cloud radius, several clumps at R 1.3-1.7; equal-radius clouds that never merge pile
  up without limit, which is the model failing, not the physics. The clumps raise
  A1/S0 to 0.19 and |r_c| to 0.008 in three orbits: a lopsided state from fragments,
  not the paper's mode, which never appeared (A1/S0 at the noise floor to t = 24).
  Final: E_diss 5.3 |E0|, 31% of the clouds in the sink, L to 1.3e-4.
  `figs/particles_t0028.0.png`, `figs/b_series.png`. The particle movie colours by
  h, which is fixed here, so its colour scale carries no information for clouds.
- **B6b start** (`20261029-152300-clouds-from-a3t100-e05`): collisions switched on in A3
  at t = 100 with the mode in place, 15 orbits (fragmentation took ~24 in B6). A first
  launch wrote dumps every 5 orbits (wrong dtmax) and was stopped after 7 s, kept in
  `B/tests/`.
- **B6b, first 3.5 orbits**: with collisions switched on in the saturated A3 torus, the
  mode survives the collapse so far: A1/S0 0.33 -> 0.23 (half an orbit) -> 0.28-0.35,
  |r_c| 0.021 -> 0.014 -> 0.021, while E_diss reaches 2.25 |E0| and 12% of the clouds
  fall into the sink. B6 from t = 0 collapsed before any mode grew; with the mode in
  place the dissipating torus keeps its lopsidedness (preliminary, 3.5 orbits).
- **A7 done** (M_tor 0.06, 250 orbits, dE/E 9.7e-4): over t = 200-250 A1/S0 0.357, k 0.21,
  |r_c| 0.0143, |r_tb| 0.233 (0.1: 0.344, 0.022, 0.218), pattern 0.0108 Omega_0 (0.0186 at
  0.1: ratio 0.58 for a mass ratio 0.6), e-fold 28.5 orbits (16.8). Below 0.1 the paper's
  self-similar picture holds: fixed amplitude and r_tb, r_c, growth rate and pattern
  speed linear in M_tor. Above (A11, 0.3) it fails. `compare_mass.png`, run `figs/`.
- **B6b done** (collisions on at t = 100, 15 orbits, e 0.5, r_cloud 0.01): the torus
  collapses to a thin disc within an orbit (|z|/R 0.23 -> 0.015, sigma_z/v_phi 0.41 ->
  0.045), radiates 4.3 |E0|, puts 24% of the clouds into the sink, and keeps the mode:
  A1/S0 0.33 -> 0.25-0.30, |r_c| 0.016-0.022, a lopsided disc with spiral arms, no
  fragmentation in 15 orbits. Order matters: dissipation from t = 0 (B6) left no mode,
  dissipation after the mode (B6b) keeps it in the collapsed disc.
- **Axis check of the published panels**: Fig. 6 ticks every 34.3 px = 0.05 arcsec, half
  width 0.299 arcsec: the +-0.3 arcsec of the synthetic maps holds to 1%; the panels are
  now drawn on the same calibrated axes.
- **B8a start** (`20261029-163625-gas-from-a3t100-hor015`): A3 at t = 100 as isothermal gas,
  H/R 0.15, r_acc 0.05, 15 orbits, the fluid pair of B6b. A first launch inherited
  r_cloud from the moddump default and stopped (fixed, kept in `B/tests/`).
- **B8a done** (A3 at t = 100 as isothermal gas, H/R 0.15, 15 orbits): the fluid limit kills
  the mode: A1/S0 0.33 -> 0.072 at t = 109 -> 0.037 at 115 (e-fold ~6 orbits); the gas
  settles to |z|/R 0.10 (its c_s), sigma_z/v_phi 0.02, spreads viscously (R90 1.37 ->
  2.0) and feeds 31% to the sink. L changes by 1.4e-3 while the sink spin takes 8% of
  L, so phantom books the accreted L; the ANGMOM alarm at 1e-3 is tight for gas.
  Against B6b (same start, clouds, mode kept at 0.25-0.30): the two limits of the cloud
  ensemble disagree. Plausible cause: pressure (c_s = 0.15 v_K) makes eccentric orbits
  precess differentially and phase-mixes the mode; the cold cloud disc (sigma ~ 0.045)
  barely precesses. If so, mode survival is set by the effective dispersion, which the
  calibrated heating (s = 0.1, 0.3) will test.
- **B8b, gas at H/R 0.08** (same start as B8a): fails after half an orbit, the tree
  overflows as the thin isothermal gas (|z|/R 0.03) fragments, Q < 1; the dispersion
  test needs sink creation for the clumps, which (b) of the B3 extension supplies.
- **Comparison figures fixed** (user's check): the palette was a hand-made GILDAS
  rainbow and the ranges rounded. Now the palettes are read pixel by pixel from the
  colour bars of Combes Fig. 6 and Fig. 12 at 400 dpi (`vis/combes25_colours.npz`), the
  ranges from their tick marks (moment 0 0.052-0.280 Jy/beam km/s, moment 1 -108..+109
  km/s, moment 2 22-194 km/s, PV levels 0.0002-0.002 Jy/beam), and each pair of panels
  carries one shared colour bar. Synthetic noise is now periodic (no bright rim at the
  map edge). Particle movies of cloud runs colour by the density of the 32 nearest
  neighbours, since a fixed softening carries no density.
- **Figure rules** (with the user): `FIGURES.md` and `scripts/vis/labels.py`. Every figure
  and movie names the model it shows (Collisionless torus model, Colliding-cloud model,
  Isothermal gas model, with their parameters in mathtext and the restart time), never
  "paper", "ours" or a run id; published data are cited in full. `regen_all.sh`
  regenerates every figure and movie of the finished runs.
- **Moment 2 of the synthetic maps was wrong** (user's check: our noise specks were pink,
  theirs black). A plain clip at 1.5 rms passes noise in ~7% of the channels over +-300
  km/s, so a noise pixel had moment 2 ~ 150 km/s, and the noise also inflated the disc's
  own moment 2. Now a smooth-and-mask selection (Dame 2011): noise specks have sigma ~ 0
  as theirs, and the cold cloud disc of B6b has moment 2 ~ 10-20 km/s against the
  observed 40-80 km/s. The earlier statement that the cloud model's moment 2 "resembles
  the observed map" was the artefact: NGC 613's gas is hotter than the collapsed cloud
  disc, which argues for the heating of the next step.
- **Heated colliding-cloud model** (B3 extension, with the user: (c) heating + (b) clump
  sinks). Code: `heat_clouds` (thermostat to sigma = heat_s v_K per annulus, relaxed over
  heat_tau local orbits, kicks keep momentum and L_z, energy in `eheat`) and
  `form_clump_sinks` (a cloud overlapped by sf_n clouds in a bound group becomes a sink
  at the centre of mass; mass, momentum, L with spin kept; internal energy in `esf`).
  Unit tests pass. Phantom test (paper torus from t = 0, s 0.1, tau 1, 3 orbits): no
  clump sinks yet (fragmentation took ~24 orbits in B6), but tau = 1 is too slow:
  collisions radiated 2.7 |E0| against 0.22 added, sigma_z/v_phi 0.076 instead of 0.1,
  |z|/R 0.04. The thermostat must act faster than collisions (~0.2 orbits): tau = 0.1
  in the runs; s stays the calibrated parameter, tau a numerical one.
- **B12 start** (`20261029-181639-heated-s01-from-a3t100`): heated clouds, s 0.1, tau 0.1,
  sf_n 16, from A3 at t = 100 with the mode in place, 20 orbits. s 0.3 follows when A6
  ends.
- **How far the models are from NGC 613, in numbers.** Combes et al. (2025) Fig. 6
  digitised with the palettes of their colour bars: in the annulus 0.02-0.05 arcsec
  (1.7-4 pc) |moment 1| (90th percentile) 86 km/s and moment 2 (median) 83 km/s, at
  0.05-0.1 arcsec 55 and 51 km/s: moment 2 / |moment 1| ~ 0.95. The synthetic maps:
  collisionless torus (A3, t 100) 107/43 and 96/34 (ratio 0.35-0.40); clouds (B6b)
  144/24 and 117/11 (0.10-0.16); gas H/R 0.15 (B8a) 127/23 and 106/13 (0.12-0.18);
  heated clouds s 0.1 (B12, t 110) 143/31 and 109/17 (0.15-0.21). Two mismatches: all
  models rotate 1.5-2 times faster in projection (M_BH sin^2 i too large for the
  assumed 3.4e7 and i = 50, both uncertain: Combes' own numbers give 2.5e7, and i comes
  from an oval assumed circular), and all are 2.5-9 times too cold relative to their
  rotation; the ratio rises as 1/sin i, so a lower inclination closes part of the gap.
  NGC 613's torus is dynamically hot (sigma/v ~ 0.5-1 for i = 25-50 deg) or its moment 2
  carries non-circular motions. A fit over i and M_BH plus heated runs at s 0.3 (and
  higher if needed) will separate the two.
- **Plan for the hotness gap** (with the user): PLAN B13-B16. First B13b, the ALMA archive
  cube of NGC 613, because Combes' stated recipe (plain clip at 1.5 rms over +-300 km/s)
  is the one that inflated our own moment 2 by noise, while their black noise specks
  suggest a mask: from the paper alone the recipe of Fig. 6 cannot be settled.
- **B13b, NGC 613 from the ALMA archive** (project 2019.1.00273.S, MOUS X3435, pipeline
  cube spw25 = CO(3-2), 6.8 km/s channels, beam 0.0137 arcsec, rms 0.42 mJy/beam; SODA
  cutout of 0.35 arcsec, `materials/data/alma_ngc613/`; `scripts/obs/ngc613_cube.py`,
  `results/m1torus/B/ngc613_cube/`). Combes' stated recipe (clip at 0.6 mJy/beam over
  +-300 km/s) gives a noise-dominated pink moment 2 everywhere, unlike their Fig. 6;
  smooth-and-mask reproduces Fig. 6 (black specks, blue-cyan torus, green centre), so
  Fig. 6 was masked and our synthetic maps are processed like the data. Numbers, annuli
  0.02-0.05 and 0.05-0.1 arcsec: mask |m1| (p90) 77 and 59, m2 (median) 47 and 32 km/s,
  ratio 0.61 and 0.54; one Gaussian per pixel (S/N > 5, few pixels) sigma 41 and 37,
  centroid p90 63 and 60. The digitised Fig. 6 (ratio ~0.95) overstated it. NGC 613 is
  dynamically warm, sigma/|v_proj| ~ 0.55-0.6; against our models by the same mask:
  collisionless torus 0.35-0.40, heated clouds s 0.1 0.15-0.21, cold clouds and gas
  0.10-0.18. The ratio scales as 1/sin i: the collisionless torus would match at
  i ~ 30 deg, the cloud discs only at i ~ 10-15 deg; B14 fits i and M_BH with the axis
  ratio of moment 0 as the second constraint.
- **A6 to 400 orbits** (`20261028-144810-weak45-16k`, weak case `0.5-45-0.3`, 16k): no
  m = 1 mode. Over 350-400 A1/S0 0.016-0.018 with A2 and A3 at the same level and the
  noise floor of 16k at 0.014; |r_c| ~ 0.002. Bannikova et al. (2026), Fig. 14 has
  A1/S0 ~0.04-0.05 already at 150 orbits and a clear rise by 300-400. Not confirmed at
  16k. Plausible causes, none tested: the weak case is closer to the threshold and the
  16k noise (and its relaxation) damps it; a different spread in i_max sampling. A6 to
  1000 or at 64k would decide; not scheduled (optional in PLAN).
- **B12 result** (heated clouds s 0.1, tau 0.1, sf_n 16, t 100-120): the mode survives
  20 orbits, A1/S0 0.33 -> 0.20, r_c grows 0.022 -> 0.037; 41% of the clouds end in clump
  sinks (N 9596/16384); E_diss 8.4 |E0| radiated and resupplied by the thermostat;
  dL/L 9e-7. sigma_z/v_K ~0.09, |z|/R 0.06. On NGC 613 by the mask recipe moment 2 /
  |moment 1| 0.15-0.21 against 0.55-0.6 observed: still too cold, and the moment-0 map
  is faint because the sinks took the dense gas.
- **B12 at s 0.3** (`20261029-194616-heated-s03-from-a3t100`) started 19:46, same
  settings, ~8 orbits/h.
- **Comparison figures in three columns**: synthetic map | NGC 613 from the archive cube
  with the same mask, beam and palette | Combes et al. (2025) Fig. 6 (PV: same slit
  through the archive cube, then Fig. 12). The archive column reproduces Fig. 6, which
  checks the axes, palettes and recipe of the whole pipeline at once. Movie titles
  shrink to end before the clock; movies labelled like figures (FIGURES.md).
- **First two columns made identical** (the user saw them differ). The archive cube had
  its own mask (3-D dilation, rms from far channels), the synthetic one a Gaussian-grown
  mask, and both panels hid m0 < 0.04-0.05, which Fig. 6 shows black. Now one function,
  `mockobs.moments`, for both; all masked pixels shown; same limits; PV unmasked in both.
  The threshold 4.75 sigma reproduces the speck coverage of Fig. 6 in 0.15-0.3 arcsec
  (15%; 3 sigma gives 41%). NGC 613 by this mask: |m1| p90 73 and 56, m2 median 39 and
  25 km/s, ratio 0.54 and 0.44 (3 sigma: 0.63 and 0.55; Gaussian sigma 0.62-0.65): the
  recipe moves the ratio by ~0.1, so models are measured with the same function only.
  Unexplained: m2 inside the torus in Fig. 6 (cyan-green, 60-100 km/s) is about twice the
  archive cube's by any of our masks, and Fig. 12 has a vertical stripe at the AGN over
  +-300 km/s that the pipeline cube lacks; their own reduction (continuum, weighting,
  channels) differs from the pipeline, not tested.
- **s 0.3 run invalid; thermostat corrected** (`20261029-194616-heated-s03-from-a3t100`,
  status INVALID). The torus flew apart: median R of the live clouds 0.95 -> 11 in 6.5
  orbits, 30% unbound by t = 102, 38% into the black hole, which took their momentum
  and drifted (|r_c| 0.40; total momentum kept). Cause: the thermostat held each of
  sigma_R, sigma_phi, sigma_z at s v_K, but Keplerian epicycles fix sigma_phi/sigma_R =
  kappa/2 Omega = 0.5 (Binney & Tremaine 2008, Sect. 4.4.3); the start (A3, t 100)
  already had sigma_R 0.30-0.35, sigma_phi 0.16-0.17 (ratio 0.48-0.53), and pumping
  sigma_phi to 0.27 drove sigma_R past 0.5. B12 (s 0.1) is consistent (ratio 0.53-0.59)
  but its effective state is sigma_R ~0.18, sigma_phi 0.10, sigma_z 0.10 v_K: labelled
  now as the floor sigma_{R,phi,z} >= 0.1 v_K. With the user: the thermostat holds
  sigma_R and sigma_z only (`heat_phi = F`, default; T keeps the old mode); unit test
  covers both modes. Restarted s 0.3 as `20261029-201331-heated-rz-s03-from-a3t100`.
- **B14, viewing geometry against the archive cube** (`scripts/obs/fit_geometry.py`,
  `results/m1torus/B/fit_geometry*/`). i 15-65 deg, M_BH 1.5-4.5e7 (A3 and B12 also
  i 15-45, M_BH 3-8e7), four azimuths about the m = 1 phase, beam, channels and noise
  of the archive cube, moments by the same function. Metrics: |m1| p90 and m2 median in
  0.02-0.05 and 0.05-0.1 arcsec, moment-0 b/a within 0.15 arcsec; observed 73/56, 39/25,
  b/a 0.90 (0.82-0.90 for any aperture and threshold). No geometry fits any model: best
  chi2 (5 metrics) 26 for the collisionless torus (i 55, M_BH 2.0e7; i = 50, 3.4e7: 37),
  104 clouds (B6b), 91 gas (B8a), 84 heated clouds (B12, at the grid edge i 65,
  M_BH 1.5e7; default 194). Geometry is not the main cause of the gap. What stays at
  every geometry: masked m2 of the inner annulus too low (A3 24, clouds and gas <= 19
  against 39), moment 0 lopsided (b/a 0.46-0.83, the m = 1 clump) against a round
  observed core. Caveat: the mask cuts the broad, multi-peaked model lines harder than
  the observed one at the same peak S/N (5.2): noise-free m2 of A3 at its best geometry
  is 49/35, as hot as NGC 613's Gaussian sigma 41/37; the heated clouds stay cold (26/17).
  So the collisionless torus is hot enough but too lopsided and has the wrong line
  profile; the dissipative models are too cold; B15 (hotter clouds) is the lever, not
  geometry.
- **s 0.3 with the corrected thermostat** (`20261029-201331-heated-rz-s03-from-a3t100`,
  sigma_R, sigma_z >= 0.3 v_K, t 100-120): the torus does not survive. 32.5% of the
  clouds fall into the black hole within ~5 orbits, then the inflow stops because the
  inner edge recedes (r10 0.39 -> 1.24); the thermostat puts 3.7 |E0| in during the
  first 5 orbits and unbinds the outer part: r50 0.79 -> 6.7, r90 1.5 -> 70 R. A1/S0
  0.33 -> 0.11, |r_c| 0.21; the synthetic NGC 613 maps are empty (the emission spreads
  far beyond 0.3 arcsec). dL/L 9e-5.
- **B12 re-read in the same terms**: no clump sinks formed in either run; B12 drains like
  a viscous disc, 41% into the black hole in 20 orbits, r10 0.36 -> 0.61, r50 0.70 ->
  1.54, r90 1.4 -> 3.4; it radiated 8.4 |E0| against 1.4 |E0| from the thermostat, so
  most of its dissipation is paid by the inflow. Early inflow rates 4.3 (B12) and 7 %
  per orbit (s 0.3): ratio 1.7, against ~3.8 from sigma_R^2 if collisional viscosity
  alone set them; at s 0.3 evaporation, not viscosity, dominates. Conclusion
  (plausible): a torus of colliding clouds with r_cl = 0.01 cannot be held at
  sigma ~ 0.3 v_K by random kicks; it drains (sigma_R ~ 0.2) or evaporates (0.35).
  NGC 613's sigma/|v| ~ 0.5 then needs rarer collisions (B15, r_cloud), resupply, or
  non-circular motions (B16). Caveat: the kick thermostat also heats clouds already
  leaving the torus (annuli out to 8 R), which may speed up the evaporation.
- **B15 start: rarer collisions at s 0.3.** N0 (covering number of an equatorial ray)
  scales as r_cl^2: 2.6 at 0.01, ~1 at 0.0062, ~0.5 at 0.0044. Only the extreme run
  (`20261029-212729-heated-rz-s03-rcl0044-from-a3t100`) goes: if even N0 ~ 0.5 drains
  or evaporates, the intermediate one is not needed (`...-rcl0062-...` stopped at its
  start). Stop criterion fixed in advance, at t = 103 against r_cl 0.01 (in the black
  hole 19.5%, R90 6.5): stop if R90 > 4 or > 15% accreted; on to 110 if R90 < 3 and
  < 10%; between, one more orbit.
- **B15 result (stopped at t = 103 by the criterion, R90 4.7 > 4).** N0 ~ 0.5 at s 0.3:
  collisions radiate 1.06 |E0| in 3 orbits against 5.75 at N0 2.6, and the inflow
  drops 19.5% -> 7.3%, but clouds still escape (6.4% unbound against 10.2%). Where
  they start: 6-8% of every initial radius 0.3-2 R, not the colder edges. Controls
  from existing runs at t = 103: no thermostat (B6b) 0.4% unbound, R90 1.5; s 0.1
  (B12) 0.5%, R90 1.7. So the escape comes from re-randomising the velocities at
  sigma 0.3 v_K (a Jeans-like tail refilled every step), nearly independent of the
  collision rate; the inflow comes from the collisions. Cooling without heating (B6b,
  N0 2.6): sigma_z 0.33 -> 0.09 in 2 orbits, |z|/R 0.23 -> 0.025; with t_cool ~ 1/N0,
  a torus that stays hot for ~100 orbits (36 Myr at 8 pc) needs N0 <~ 0.04, clouds
  that practically never collide: the collisionless limit, which B14 found hot enough
  (noise-free m2) but too lopsided. Plausible reading: NGC 613's dispersion is either
  that of a nearly collisionless cloud ensemble, a short-lived or fed state, or
  non-circular motions (B16). The intermediate run r_cl 0.0062 is not needed.
