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
