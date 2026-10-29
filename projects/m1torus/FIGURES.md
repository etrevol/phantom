# Figures and movies of m1torus: naming and style

Agreed with the user on 2026-10-29. Every figure and movie follows these rules; a new
script follows them from its first version. The names come from
`scripts/vis/labels.py`, never typed by hand.

## Names of our models

A figure names the physical model a run computes, with the parameters that set it, as
in a journal. It never says "paper", "ours", "phantom run" or a run id: our
implementation differs from Bannikova et al. (2026) and is not "the paper".

| model | name on figures | parameters shown |
|---|---|---|
| collisionless N-body torus | Collisionless torus model | M_tor, N; e_max, i_max, R_0 when not 0.5, 60 deg, 0.3; fixed epsilon; epsilon_c |
| colliding (sticky) clouds | Colliding-cloud model | M_tor, N, e, r_cl |
| SPH gas converted from the torus | Isothermal gas model | M_tor, N, H/R |
| SPH gas disc set up as a disc | Isothermal gas-disc model | M_tor, N, R range, H/R |
| heated clouds with star formation (next) | Heated colliding-cloud model | as clouds, plus the heating and sink parameters |

A run restarted from another model says so: "..., collisions from t = 100 P_orb of the
collisionless model". Symbols are mathtext: $M_{\rm tor}$, $r_{\rm cl}$, $i_{\max}$,
$P_{\rm orb}$. Several runs on one figure: the title carries what they share, the
legend only what differs (`legend_labels`).

## Published data

Cited in full, with the table or figure: "Bannikova et al. (2026), Table 2",
"Combes et al. (2025), Fig. 6", "López-Rodríguez et al. (2026), Fig. 2f",
"Izumi et al. (2023), Fig. 2C". The digitised curves of Bannikova et al. Fig. A.1 keep
their own run names (their N and epsilon).

## Layout

- English text, capitalised titles as in a journal, units in square brackets, time in
  $P_{\rm orb}$.
- No run ids or directory names in figures or movies; file names may carry them.
- Title: the subject on the first line, the model name on the second, then the time.
- Synthetic observations: our map on the left ("Synthetic CO(3-2), moment 0"), the
  published panel on the right on the same calibrated axes ("Combes et al. (2025),
  Fig. 6, moment 0"), one colour bar for the pair with the palette and range read from
  the published colour bar. The noise, beam, channels and clip are the observers'; the
  moments use a smooth-and-mask selection (Dame 2011).
- Movies: the model name across the top, the clock at the right; cloud runs colour by
  the density of the 32 nearest neighbours.
- After a change of these rules, regenerate every figure and movie of the finished runs.
