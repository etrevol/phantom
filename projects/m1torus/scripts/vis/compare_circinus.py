#!/usr/bin/env python3
"""Synthetic Circinus maps of a dump next to the published panels.

Ours on the left, theirs on the right:
  <stem>_circ_dust.png  gas column seen edge-on at PA 50 deg through a 2 pc beam (dust
                        mass, optically thin), against Lopez-Rodriguez et al. (2026)
                        Figs. 2f (ALMA 700 um) and 3e (dust of the clumpy-torus fit);
                        axis ratio of the map against their 5 x 3 pc disc
  <stem>_circ_pv.png    PV along the major axis against Izumi et al. (2023) Fig. 2C
  <stem>_circ_hr.png    sigma_z / v_phi of the particles against r (their Fig. 2E proxy
                        of h/r), with the sink's accretion rate in Msun/yr against
                        their 0.2-0.34 at 0.27 pc
Units: R = 2.5 pc, M_BH = 1.7e6 Msun (mockobs.CIRCINUS).

    compare_circinus.py DUMP [--out DIR]
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.ndimage import gaussian_filter  # noqa: E402

import mockobs as mo  # noqa: E402
import labels as L  # noqa: E402
from compare_obs import panel_ref  # noqa: E402


def column_map(p, t, phi, beam_pc=2.0, field_pc=7.0, pix_pc=0.1):
    east, north, _, _ = mo.to_sky(p["x"], p["v"], t.inc, t.pa, phi)
    e, n = east * t.R_pc, north * t.R_pc
    edges = np.arange(-field_pc, field_pc + 1e-9, pix_pc)
    img, _, _ = np.histogram2d(n, e, bins=[edges, edges], weights=p["m"])
    img = gaussian_filter(img, beam_pc / mo.FWHM / pix_pc)
    return img, edges


def axis_ratio(img, edges):
    """Minor/major axis from the second moments of the brighter half of the map."""
    c = 0.5 * (edges[1:] + edges[:-1])
    Y, X = np.meshgrid(c, c, indexing="ij")
    w = np.where(img > 0.5 * img.max(), img, 0.0)
    w = w / w.sum()
    mx, my = (w * X).sum(), (w * Y).sum()
    cov = np.array([[(w * (X - mx) ** 2).sum(), (w * (X - mx) * (Y - my)).sum()],
                    [(w * (X - mx) * (Y - my)).sum(), (w * (Y - my) ** 2).sum()]])
    ev = np.sort(np.linalg.eigvalsh(cov))
    return float(np.sqrt(ev[0] / ev[1]))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dump")
    ap.add_argument("--phi", type=float, default=None)
    ap.add_argument("--label", default="", help="ignored; the run directory names the model")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    t = mo.CIRCINUS
    p = mo.load(a.dump)
    phase = mo.m1_phase(p["x"], p["m"])
    phi = phase + 90.0 if a.phi is None else a.phi
    out = Path(a.out or Path(a.dump).parent / "mockobs")
    out.mkdir(parents=True, exist_ok=True)
    stem = Path(a.dump).name
    tag = (L.wrapped(L.model_label(Path(a.dump).resolve().parent), 110) + "\n"
           + rf"$t = {p['time'] / (2 * np.pi):.0f}\,P_{{\rm orb}}$")
    lr, iz = L.REF["lopezrodriguez26"], L.REF["izumi23"]

    # dust column, edge-on
    img, edges = column_map(p, t, phi)
    q = axis_ratio(img, edges)
    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    axs[0].imshow(img, origin="lower", extent=[edges[0], edges[-1]] * 2, cmap="magma")
    axs[0].set_xlim(edges[-1], edges[0])
    axs[0].set_xlabel("E offset [pc]")
    axs[0].set_ylabel("N offset [pc]")
    axs[0].set_title(f"Synthetic dust column, 2 pc beam; axis ratio {q:.2f} (their disc 3/5 = 0.6)",
                     fontsize=9)
    panel_ref(axs[1], "lr26_fig2f_700um", rf"{lr}, Fig. 2f, ALMA 700 $\mu$m")
    panel_ref(axs[2], "lr26_fig3e_dust", f"{lr}, Fig. 3e, dust of the clumpy-torus fit")
    fig.suptitle(f"Circinus\n{tag}", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / f"{stem}_circ_dust.png", dpi=110)
    plt.close(fig)

    # PV along the major axis, HCN-like beam 0.03 arcsec, velocities about systemic 446 km/s
    tp = mo.Target(**{**t.__dict__, "beam_arcsec": 0.03, "pix_arcsec": 0.012, "field_arcsec": 0.9,
                      "rms_mjy": 0.0, "clip_mjy": 0.0})
    cb = mo.cube(p, tp, phi)
    pvd = mo.pv(cb, tp, pa=tp.pa, halfwidth=tp.beam_arcsec)
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    axs[0].pcolormesh(cb["xs"], cb["vs"] + 446.0, pvd, cmap="Greys", shading="nearest")
    r = np.linspace(1e-3, 0.3, 300)
    vk = np.sqrt(mo.G_PC * t.M_bh / (r * mo.pc_per_arcsec(t))) * np.sin(np.radians(t.inc))
    axs[0].plot(r, 446 + vk, "r", -r, 446 - vk, "r", lw=0.8)
    axs[0].set_ylim(300, 600)
    axs[0].set_xlim(-0.9, 0.9)
    axs[0].set_xlabel("offset [arcsec]")
    axs[0].set_ylabel("velocity [km/s]")
    axs[0].set_title(r"Synthetic PV along the major axis; red: Keplerian, $1.7\times10^6\,M_\odot$",
                     fontsize=9)
    panel_ref(axs[1], "izumi23_fig2c_pv", f"{iz}, Fig. 2C (offsets to $\\pm0.9$ arcsec)")
    fig.suptitle(f"Circinus\n{tag}", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / f"{stem}_circ_pv.png", dpi=110)
    plt.close(fig)

    # sigma_z / v_phi against radius, from the particles; accretion rate of the sink
    x, v = p["x"], p["v"]
    R = np.hypot(x[:, 0], x[:, 1]) * t.R_pc
    vphi = (x[:, 0] * v[:, 1] - x[:, 1] * v[:, 0]) / np.hypot(x[:, 0], x[:, 1])
    bins = np.arange(0.25, 3.0 * t.R_pc, 0.5)
    prof = []
    for a0, b0 in zip(bins[:-1], bins[1:]):
        s = (R >= a0) & (R < b0)
        if s.sum() > 20:
            prof.append((0.5 * (a0 + b0), np.std(v[s, 2]) / np.mean(np.abs(vphi[s]))))
    prof = np.array(prof)
    mdot = sink_mdot(Path(a.dump).parent, t)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
    axs[0].plot(prof[:, 0], prof[:, 1], "ko-")
    axs[0].set_xlabel("radius [pc]")
    axs[0].set_ylabel("sigma_z / v_phi")
    axs[0].set_ylim(0, 0.7)
    axs[0].set_xlim(0, 16)
    axs[0].set_title(r"$\sigma_z/v_\phi$ of all mass; accretion "
                     + (rf"{mdot:.2f} $M_\odot\,$yr$^{{-1}}$" if np.isfinite(mdot) else "n/a")
                     + f" ({iz}: 0.2-0.34 at 0.27 pc)", fontsize=9)
    panel_ref(axs[1], "izumi23_fig2e_hr", f"{iz}, Fig. 2E")
    fig.suptitle(f"Circinus\n{tag}", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / f"{stem}_circ_hr.png", dpi=110)
    plt.close(fig)
    print(out / f"{stem}_circ_dust.png")
    print(f"axis ratio {q:.2f}, accretion {mdot:.3f} Msun/yr")


def sink_mdot(run, t, window=1.0):
    """Mean accretion rate of the sink over the last window orbits [Msun/yr]."""
    import phantom_io as pio
    files = sorted(run.glob("m1torusSink0001N[0-9][0-9].ev"))
    if not files:
        return np.nan
    ev = {}
    for f in files:
        d = pio.read_ev(f)
        for k, v in d.items():
            ev.setdefault(k, []).append(v)
    ev = {k: np.concatenate(v) for k, v in ev.items()}
    key = next((k for k in ev if "macc" in k.lower()), None)
    if key is None or len(ev["time"]) < 2:
        return np.nan
    tt, m = ev["time"], ev[key]
    sel = tt >= tt[-1] - window * 2 * np.pi
    if sel.sum() < 2:
        return np.nan
    rate = (m[sel][-1] - m[sel][0]) / (tt[sel][-1] - tt[sel][0])   # code mass / code time
    tunit_yr = t.R_pc / np.sqrt(mo.G_PC * t.M_bh / t.R_pc) * 0.9778e6
    return float(rate * t.M_bh / tunit_yr)


if __name__ == "__main__":
    main()
