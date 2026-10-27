#!/usr/bin/env python3
"""2D maps of a PPI torus run: mid-plane polar views, meridional cuts, animations.

Reads <prefix>_ppimaps.bin, which the live analysis writes every dtmax, so a movie
costs no extra dumps. Run from inside a run directory or point --run at one.

    vis2d.py                                   # polar, deviation, meridional at the last frame
    vis2d.py --mode polar --time 5             # nearest frame to 5 orbits
    vis2d.py --mode polar --normalize azimuthal
    vis2d.py --animate --anim deviation --every 2
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LogNorm, Normalize  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ppi_data as pdat  # noqa: E402

CONFIG = dict(dpi=200, fps=10, grid_alpha=0.3)

FIELDS = {
    "rho": dict(label=r"$\rho$ (Density)", cmap="inferno", log=True),
    "p": dict(label=r"$P$ (Pressure)", cmap="plasma", log=True),
    "v_r": dict(label=r"$v_R$ (Radial Velocity)", cmap="RdBu_r", log=False),
    "v_phi": dict(label=r"$v_\phi$ (Azimuthal Velocity)", cmap="coolwarm", log=False),
}
# Fields whose azimuthal mean is positive; the deviation is shown as a fraction of it.
# <v_R> is near zero, so dividing by it would turn a small deviation into noise.
RELATIVE = ("rho", "p")


def deviation(name, a):
    """f/<f>_phi - 1 for rho and p, f - <f>_phi otherwise; a has shape (nphi, nr)."""
    mean = a.mean(axis=0, keepdims=True)
    if name in RELATIVE:
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(mean > 0, a / mean - 1.0, 0.0)
    return a - mean


def _norm(name, a, normalize, pct=(2.0, 98.0)):
    info = FIELDS[name]
    if normalize:
        v = np.percentile(np.abs(a), 99.5) or 1.0
        return Normalize(-v, v), "RdBu_r"
    if info["log"]:
        pos = a[a > 0]
        if pos.size == 0:
            return Normalize(0, 1), info["cmap"]
        vmax = pos.max()
        return LogNorm(vmin=max(pos.min(), vmax * 1e-6), vmax=vmax), info["cmap"]
    if name == "v_phi":
        lo, hi = np.percentile(a[a != 0], pct) if np.any(a != 0) else (0, 1)
        return Normalize(lo, hi), info["cmap"]
    v = np.abs(a).max() or 1.0
    return Normalize(-v, v), info["cmap"]


def _info_box(fig, lines):
    if lines:
        fig.text(0.012, 0.988, "\n".join(lines), ha="left", va="top", fontsize=8.5,
                 color="#1a1a19", linespacing=1.5,
                 bbox=dict(boxstyle="round", fc="white", ec="#b8b8b4", alpha=0.9))


def _title(frame, p_orb, what):
    return f"{what} - t = {frame['time']:.3f} ({frame['time'] / p_orb:.2f} orbits)"


def _polar_axes(fig):
    return [fig.add_subplot(2, 2, k + 1, projection="polar") for k in range(4)]


def _draw_polar(axes, frame, normalize, norms=None):
    """pcolormesh on four polar axes; returns the mesh artists."""
    out = []
    for k, (ax, name) in enumerate(zip(axes, FIELDS)):
        a = frame[name]
        a = deviation(name, a) if normalize else a
        if norms is None:
            norm, cmap = _norm(name, a, normalize)
        else:
            norm, cmap = norms[k]
        mesh = ax.pcolormesh(frame["phi_edges"], frame["r_edges"], a.T, cmap=cmap,
                             norm=norm, shading="flat")
        ax.set_ylim(0, frame["r_edges"][-1])
        ax.grid(True, alpha=CONFIG["grid_alpha"])
        label = FIELDS[name]["label"]
        if normalize:
            label = (r"$\delta$" + label + r"$/\langle\cdot\rangle$" if name in RELATIVE
                     else r"$\delta$" + label)
        ax.set_title(label, fontsize=12, fontweight="bold", pad=16)
        out.append(mesh)
    return out


def plot_polar(frame, info, out_dir, normalize=False):
    p_orb = info["model"].p_orb
    fig = plt.figure(figsize=(14, 10))
    axes = _polar_axes(fig)
    meshes = _draw_polar(axes, frame, normalize)
    for ax, mesh in zip(axes, meshes):
        fig.colorbar(mesh, ax=ax, pad=0.1)
    what = "Mid-plane deviation from the azimuthal mean" if normalize else "Mid-plane"
    fig.suptitle(_title(frame, p_orb, what), fontsize=16, fontweight="bold")
    _info_box(fig, pdat.info_lines(info))
    fig.tight_layout()
    tag = "deviation" if normalize else "polar"
    f = out_dir / f"{tag}_t{frame['time'] / p_orb:06.2f}.png"
    fig.savefig(f, dpi=CONFIG["dpi"], bbox_inches="tight")
    plt.close(fig)
    return f


def plot_meridional(frame, info, out_dir):
    p_orb = info["model"].p_orb
    fig, axes = plt.subplots(2, 1, figsize=(12, 9))
    for ax, key, xe, title in (
            (axes[0], "slice", frame["x_edges"], r"$\rho$ in the plane $y = 0$"),
            (axes[1], "avg", frame["r_edges"], r"$\langle\rho\rangle_\phi(R, z)$")):
        a = frame[key]
        pos = a[a > 0]
        norm = LogNorm(vmin=pos.max() * 1e-6, vmax=pos.max()) if pos.size else None
        mesh = ax.pcolormesh(xe, frame["z_edges"], a, cmap="inferno", norm=norm,
                             shading="flat")
        ax.set_aspect("equal")
        ax.set_xlabel("x" if key == "slice" else "R")
        ax.set_ylabel("z")
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.grid(True, alpha=CONFIG["grid_alpha"])
        fig.colorbar(mesh, ax=ax, pad=0.02)
    fig.suptitle(_title(frame, p_orb, "Meridional cuts"), fontsize=15, fontweight="bold")
    _info_box(fig, pdat.info_lines(info))
    fig.tight_layout()
    f = out_dir / f"meridional_t{frame['time'] / p_orb:06.2f}.png"
    fig.savefig(f, dpi=CONFIG["dpi"], bbox_inches="tight")
    plt.close(fig)
    return f


def animate(frames, info, out_dir, kind="deviation", fps=CONFIG["fps"]):
    """mp4 of the four mid-plane fields (polar), their deviation, or both densities."""
    from matplotlib.animation import FFMpegWriter, FuncAnimation

    p_orb = info["model"].p_orb
    last = frames[-1]
    if kind in ("polar", "deviation"):
        normalize = kind == "deviation"
        fig = plt.figure(figsize=(14, 10))
        axes = _polar_axes(fig)
        # colour scales fixed from the last frame, where the mode is strongest
        norms = []
        for name in FIELDS:
            a = deviation(name, last[name]) if normalize else last[name]
            norms.append(_norm(name, a, normalize))
        meshes = _draw_polar(axes, frames[0], normalize, norms)
        for ax, mesh in zip(axes, meshes):
            fig.colorbar(mesh, ax=ax, pad=0.1)

        def arrays(fr):
            return [(deviation(n, fr[n]) if normalize else fr[n]).T.ravel() for n in FIELDS]
    else:
        fig = plt.figure(figsize=(16, 8))
        axes = [fig.add_subplot(1, 2, k + 1, projection="polar") for k in range(2)]
        r0 = last["rho"]
        dv = deviation("rho", r0)
        norms = [_norm("rho", r0, False), _norm("rho", dv, True)]
        meshes = []
        for ax, a, (norm, cmap), t in zip(axes, (frames[0]["rho"],
                                                deviation("rho", frames[0]["rho"])),
                                          norms, (FIELDS["rho"]["label"],
                                                  r"$\delta\rho/\langle\rho\rangle_\phi$")):
            m = ax.pcolormesh(last["phi_edges"], last["r_edges"], a.T, cmap=cmap, norm=norm,
                              shading="flat")
            ax.set_ylim(0, last["r_edges"][-1])
            ax.grid(True, alpha=CONFIG["grid_alpha"])
            ax.set_title(t, fontsize=13, fontweight="bold", pad=16)
            fig.colorbar(m, ax=ax, pad=0.1)
            meshes.append(m)

        def arrays(fr):
            return [fr["rho"].T.ravel(), deviation("rho", fr["rho"]).T.ravel()]

    title = fig.suptitle("", fontsize=15, fontweight="bold")
    _info_box(fig, pdat.info_lines(info))
    fig.tight_layout()

    def update(k):
        fr = frames[k]
        for m, a in zip(meshes, arrays(fr)):
            m.set_array(a)
        title.set_text(f"Frame {k + 1:4d}/{len(frames)} | t = {fr['time']:.3f} "
                       f"({fr['time'] / p_orb:.2f} orbits)")
        return meshes + [title]

    anim = FuncAnimation(fig, update, frames=len(frames), blit=False)
    f = out_dir / f"anim_{kind}.mp4"
    anim.save(f, writer=FFMpegWriter(fps=fps, bitrate=6000), dpi=120)
    plt.close(fig)
    return f


def pick(frames, t_orbits, p_orb):
    if t_orbits is None:
        return frames[-1]
    times = np.array([f["time"] for f in frames]) / p_orb
    return frames[int(np.argmin(np.abs(times - t_orbits)))]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", default=".", help="run directory (default: here)")
    ap.add_argument("--mode", default="all", choices=["all", "polar", "meridional"])
    ap.add_argument("--normalize", default="none", choices=["none", "azimuthal"])
    ap.add_argument("--time", type=float, default=None, help="orbits; default last frame")
    ap.add_argument("--animate", action="store_true")
    ap.add_argument("--anim", default="deviation", choices=["polar", "deviation", "density"])
    ap.add_argument("--every", type=int, default=1, help="use every n-th frame")
    ap.add_argument("--fps", type=int, default=CONFIG["fps"])
    ap.add_argument("--out", default=None, help="output directory (default <run>/figs)")
    a = ap.parse_args(argv)

    info = pdat.run_info(a.run)
    out_dir = Path(a.out) if a.out else info["dir"] / "figs"
    out_dir.mkdir(parents=True, exist_ok=True)
    every = a.every if a.animate else 1
    frames = pdat.read_maps(info["maps"], every=every)
    if not frames:
        print(f"no frames in {info['maps']}")
        return 1
    p_orb = info["model"].p_orb

    if a.animate:
        print(f"{len(frames)} frames -> ", end="", flush=True)
        print(animate(frames, info, out_dir, a.anim, a.fps))
        return 0
    fr = pick(frames, a.time, p_orb)
    done = []
    if a.mode in ("all", "polar"):
        if a.mode == "all" or a.normalize == "none":
            done.append(plot_polar(fr, info, out_dir, normalize=False))
        if a.mode == "all" or a.normalize == "azimuthal":
            done.append(plot_polar(fr, info, out_dir, normalize=True))
    if a.mode in ("all", "meridional"):
        done.append(plot_meridional(fr, info, out_dir))
    for f in done:
        print(f)
    return 0


if __name__ == "__main__":
    sys.exit(main())
