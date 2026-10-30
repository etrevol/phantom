"""The galaxies whose archive cubes the models are fitted to, and their adopted numbers.

Combes et al. (2025), Tables 1 and 3 (distance, systemic redshift, torus radius, gas
mass, inclination, position angle, black-hole mass); the cubes are the ALMA pipeline
products of project 2019.1.00273.S (materials/data/alma_<name>/README.txt). The
velocity window covers the line: NGC 1672 is edge-on with M_BH 4.5e7, ~440 km/s at 1 pc.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
C_KMS = 299792.458

GALAXIES = {
    "ngc613": dict(name="NGC 613", cube=ROOT / "materials/data/alma_ngc613/ngc613_co32_spw25_cut035.fits",
                   v_sys=1481.0, D_Mpc=17.2, R_pc=8.0, m_gas=1.0e7, inc=50.0, pa=340.0, M_bh=3.4e7,
                   window=300.0, clip_mjy=0.6, field=0.3,
                   out=ROOT / "results/m1torus/B/ngc613_cube"),
    "ngc1672": dict(name="NGC 1672", cube=ROOT / "materials/data/alma_ngc1672/ngc1672_spw25_cube_cut05.fits",
                    v_sys=0.00444 * C_KMS, D_Mpc=11.4, R_pc=14.0, m_gas=0.2e7, inc=90.0, pa=70.0, M_bh=4.5e7,
                    window=450.0, clip_mjy=0.55, field=0.45,
                    out=ROOT / "results/m1torus/B/ngc1672_cube"),
}
