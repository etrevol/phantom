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
