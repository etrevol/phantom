# Slow m = 1 mode in a self-gravitating torus (SETUP=m1torus)

Bannikova et al. (2026, A&A, arXiv:2604.11528) find, with direct N-body (phi-GPU, N = 128k),
that a thick collisionless torus of mass 0.1 M_c around a central mass grows a global,
slowly precessing m = 1 mode on its own, keeps it for 1000 orbits, and displaces the
central mass from the barycentre. This project asks three questions with phantom:

| | question | needs |
|---|---|---|
| **A** | does an independent code reproduce the mode: onset, amplitude, harmonic spectrum, central-mass offset? | collisionless particles, self-gravity, a live sink |
| **B** | how much dissipation (a gas fraction that cools) does the mode survive, and what does the gas do? | A, plus SPH gas with cooling |
| **C** | how does a circumnuclear gas disc at 0.1-0.3 R_tor respond to the lopsided torus: forced eccentricity and inflow (coplanar), precession, warp or tearing (tilted)? | A, B's gas machinery |

Order: A, then B, then C. Tasks, status and budgets: [`PLAN.md`](PLAN.md); dated results
and decisions: [`LOG.md`](LOG.md).

Paper: `materials/papers/2026_bannikova_2604.11528_*.pdf`.

Status: see PLAN.md.
