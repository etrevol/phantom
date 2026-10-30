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
| colliding (sticky) clouds | Colliding-cloud model | M_tor, N, e (or "magnetised, $v_{\rm A} = ...\,v_{\rm K}(R_{\rm tor})$" when e follows the collision speed), r_cl |
| SPH gas converted from the torus | Isothermal gas model | M_tor, N, H/R |
| SPH gas disc set up as a disc | Isothermal gas-disc model | M_tor, N, R range, H/R |
| collisionless stars plus colliding gas clouds | Two-component model | as clouds, plus "gas clouds f%, collisionless stars"; with inflow "inflow $\dot M = ...\,M_c\,P_{\rm orb}^{-1}$ at ...$\,R_{\rm tor}$" |
| clouds with the thermostat (and clump sinks) | Heated colliding-cloud model | as clouds, plus the thermostat floor $\sigma_{R,z} \geq s\,v_{\rm K}$ ($\sigma_{R,\phi,z}$ for runs before 2026-10-29 21:00), $\tau_{\rm heat}$ if not 1, the clump size that makes a sink |

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
- Synthetic observations of NGC 613: three columns on the same calibrated axes - our
  map ("Synthetic CO(3-2), moment 0"), the ALMA archive cube processed the same way
  ("NGC 613, moment 0: ALMA 2019.1.00273.S archive cube, same processing"), the
  published panel ("Combes et al. (2025), Fig. 6, moment 0"); the PV diagram likewise
  (same slit through the archive cube, then Fig. 12). One colour bar per row with the
  palette and range read from the published colour bar. The noise, beam, channels and
  clip are the observers'; the moments use a smooth-and-mask selection (Dame 2011),
  one function (`mockobs.moments`, 4.75 sigma, calibrated on the speck coverage of
  Fig. 6) for the synthetic and the archive cube; every pixel the mask keeps is shown,
  below the scale in its lowest colour, as in Fig. 6. PV diagrams are unmasked in all
  columns, as Fig. 12. Other targets without a cube we hold (Circinus): two columns, ours and theirs.
- Movies are labelled like figures: the model name with its parameters across the top
  (from `labels.py`; a restart suffix on the second line), the clock in $P_{\rm orb}$ at
  the right, axes with units, the frame (inertial or "corotating with the mode")
  stated. No run ids in the frame. Every finished run gets its movie
  (`anim_particles.py`; collisionless runs also `--corotate`), made when the run ends
  and again by `regen_all.sh`. Cloud runs colour by the density of the 32 nearest
  neighbours. Two-component runs colour the stars by density and draw the gas clouds in
  one colour (magenta) over them; their stars-and-gas time series is `vis_twocomp.py`.
- After a change of these rules, regenerate every figure and movie of the finished runs.
