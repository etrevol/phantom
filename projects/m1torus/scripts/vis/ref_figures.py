#!/usr/bin/env python3
"""Crops of the published panels that B compares with, kept local (results/ is ignored).

Each entry: PDF in materials/papers, page, box in fractions of the page (left, top,
right, bottom), and the field the panel shows, so a synthetic map can be drawn on the
same axes next to it. Renders at 200 dpi with pdftoppm.

    ref_figures.py [--out results/m1torus/B/reference]
"""
import argparse
import subprocess
import tempfile
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parents[4]
PAPERS = ROOT / "materials" / "papers"
DPI = 200

# boxes measured on the 200 dpi renders (page 1654 x 2339 px for A&A)
PANELS = {
    # Combes et al. (2025), NGC 613, Fig. 6: CO(3-2) moments 0, 1, 2, field +-0.3 arcsec
    "combes25_fig6_m0": ("2026_combes_2509.08635v2", 6, (185, 243, 596, 655), "arcsec +-0.3"),
    "combes25_fig6_m1": ("2026_combes_2509.08635v2", 6, (597, 243, 1007, 655), "arcsec +-0.3"),
    "combes25_fig6_m2": ("2026_combes_2509.08635v2", 6, (1008, 243, 1418, 655), "arcsec +-0.3"),
    # the whole figure with axes and colour bars, for reading the scales
    "combes25_fig6": ("2026_combes_2509.08635v2", 6, (100, 180, 1500, 760), "full"),
    # Fig. 12: PV along PA -20 deg, 0.3 arcsec strip; Fig. 13: integrated CO and HCO+ spectra
    # with their velocity axis. Fig. 6 field: ticks every 34.3 px = 0.05 arcsec on the
    # frame, half-width 205 px = 0.299 arcsec (measured on the render)
    "combes25_fig12_pv": ("2026_combes_2509.08635v2", 8, (855, 150, 1545, 710), "full"),
    "combes25_fig13_spec": ("2026_combes_2509.08635v2", 8, (845, 960, 1185, 1555), "full"),
    # Circinus, Lopez-Rodriguez et al. (2026), Fig. 2f ALMA 700 um dust, Fig. 2l HCN(3-2),
    # Fig. 3d 700 um with the clumpy-torus contours, Fig. 3e its dust distribution
    # (AASTeX pages, 1700 x 2200 px); 2 pc bars in the panels
    "lr26_fig2f_700um": ("2026_lopez-rodriguez_2506.08077v1", 4, (575, 490, 865, 775), "2 pc bar"),
    "lr26_fig2l_hcn": ("2026_lopez-rodriguez_2506.08077v1", 4, (1278, 812, 1562, 1090), "2 pc bar"),
    "lr26_fig3d_700um": ("2026_lopez-rodriguez_2506.08077v1", 6, (1040, 425, 1305, 705), "2 pc bar"),
    "lr26_fig3e_dust": ("2026_lopez-rodriguez_2506.08077v1", 6, (1305, 425, 1565, 705), "2 pc bar"),
    # Izumi et al. (2023), Fig. 2C major-axis PV of [CI] and HCN(3-2), Fig. 2E sigma/V_rot
    "izumi23_fig2c_pv": ("2023_izumi_2305.03993v2", 11, (1092, 140, 1604, 667), "full"),
    "izumi23_fig2e_hr": ("2023_izumi_2305.03993v2", 11, (827, 680, 1604, 1100), "full"),
}


def render(pdf, page, tmp):
    stem = Path(tmp) / f"p{page}"
    subprocess.run(["pdftoppm", "-f", str(page), "-l", str(page), "-r", str(DPI), "-png",
                    str(pdf), str(stem)], check=True)
    return next(Path(tmp).glob(f"p{page}*.png"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=str(ROOT / "results" / "m1torus" / "B" / "reference"))
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    cache = {}
    with tempfile.TemporaryDirectory() as tmp:
        for name, (stem, page, box, _) in PANELS.items():
            pdf = next(PAPERS.glob(stem + "*.pdf"))
            key = (pdf, page)
            if key not in cache:
                sub = Path(tmp) / f"{len(cache)}"
                sub.mkdir()
                cache[key] = Image.open(render(pdf, page, sub)).copy()
            cache[key].crop(box).save(out / f"{name}.png")
            print(out / f"{name}.png")


if __name__ == "__main__":
    main()
