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
