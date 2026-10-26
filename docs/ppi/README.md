# Papaloizou-Pringle torus (SETUP=ppitorus)

The constant specific angular momentum torus of Papaloizou & Pringle (1984, MNRAS 208,
721) around a fixed point mass, in 3D SPH, with and without a controlled seed for the
non-axisymmetric instability. Parameters and units mirror the 2D Athena++ torus of
athena-collab (C' = 0.2, gamma = 1.3, L_0 = 500 r_g, M = 4.5e7 M_sun), with presets for
PP84's Models 1-3. Derivations: [`formulary/formulary-en.pdf`](formulary/formulary-en.pdf).

## Quick start

```bash
scripts/ppi/run/run.sh --name unseeded -- --np 1e5 --n-orbits 12
scripts/ppi/run/run.sh --name m2 -- --np 1e5 --pert-kind single --pert-m 2 --pert-amp 1e-2
cd results/runs/<date>-m2
python3 ../../../scripts/ppi/vis/vismodes.py      # A_m/A_0 against linear theory
python3 ../../../scripts/ppi/vis/vis2d.py         # polar maps, deviation, meridional
python3 ../../../scripts/ppi/vis/vis2d.py --animate --anim deviation
python3 ../../../scripts/ppi/vis/vishst.py
```

`scripts/ppi/theory/torus_model.py --help` lists every setup option; `--preset pp84-2`
switches to PP84's Model 2. `scripts/ppi/theory/ppi3d.py` gives the linear growth rates
(`--pp84` reproduces their Table 1, `--dim 2` the cylinder of the Athena++ runs).

## Variants

| what | options |
|---|---|
| unseeded, nsym-fold symmetric lattice | `--lattice rings --nsym 32 --pert-kind none` |
| stock-like rings, no symmetry | `--nsym 1` |
| Monte Carlo positions | `--lattice random` |
| one mode | `--pert-kind single --pert-m M --pert-amp A` |
| modes 1..M, random phases | `--pert-kind multi --pert-m M --pert-amp A` |
| white noise in v_R | `--pert-kind noise --pert-amp A` |

The seed is `delta v_R = A c_s (f/f_c) S(phi)`, as in athena-vv.

## Output

The live analysis (built into phantom) writes every `dtmax` (20 per orbit), whether or
not a dump is written:

- `ppitorus_ppi.ev`: torus mass, L_z, centre of mass, mode amplitudes A_m and F_m
  (m = 0..16) of the surface density and of the radial mass flux, phases of m = 1..4;
- `ppitorus_ppimaps.bin`: mid-plane (R, phi) maps of rho, p, v_R, v_phi; rho in y = 0;
  the azimuthal mean rho(R, z).

Full dumps are written `n_dumps_orbit` times per orbit and read with splash or
`scripts/ppi/vis/ppi_data.py`.
