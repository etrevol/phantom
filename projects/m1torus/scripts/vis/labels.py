"""Names on every m1torus figure and movie; the rules are in projects/m1torus/FIGURES.md.

Our runs are named by the model they compute, never "paper", "ours" or a run id: the
collisionless torus model, the colliding-cloud model, the isothermal gas model, each
with the parameters that set it, in mathtext as in a journal. Published data are cited
in full ("Bannikova et al. (2026), Table 2").

    model_label(run_dir)        "Colliding-cloud model ($M_{\\rm tor} = 0.1$, $N = 16\\,384$,
                                $e = 0.5$, $r_{\\rm cl} = 0.01$), collisions from
                                $t = 100\\,P_{\\rm orb}$ of the collisionless model"
    model_parts(run_dir)        (name, {key: text}, suffix)
    legend_labels([run_dirs])   a common title and, per run, the parameters that differ
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "common"))
import phantom_io as pio  # noqa: E402

REF = {
    "bannikova26": "Bannikova et al. (2026)",
    "combes25": "Combes et al. (2025)",
    "lopezrodriguez26": "López-Rodríguez et al. (2026)",
    "izumi23": "Izumi et al. (2023)",
}
# the canonical torus 0.5-60-0.3 of Bannikova et al. (2026); other geometries are spelled out
CANON = {"e_max": 0.5, "i_max": 60.0, "R0": 0.3}
TEX = {"e_max": r"e_{\max}", "i_max": r"i_{\max}", "R0": r"R_0"}
P_ORB = 6.283185307179586


def _num(x):
    return float(str(x).split()[0])


def _fmt(x):
    return f"{x:g}"


def _math(s):
    return "$" + s + "$"


def model_parts(run):
    """Model name, its parameters as mathtext, and how it started."""
    run = Path(run)
    setup = pio.read_setup(run / "m1torus.setup")
    inp = pio.read_setup(run / "m1torus.in")
    prov_f = run / "provenance.json"
    prov = json.loads(prov_f.read_text()) if prov_f.exists() else {}
    r_cloud = _num(inp.get("r_cloud", 0) or 0)
    ieos = str(inp.get("ieos", "")).strip()
    geometry = int(_num(setup.get("geometry", 0) or 0))
    n = f"{int(_num(setup['np'])):,}".replace(",", r"\,")
    pars = {"M_tor": _math(r"M_{\rm tor} = " + _fmt(_num(setup["M_tor"]))),
            "N": _math("N = " + n)}
    if geometry == 1:
        name = "Isothermal gas-disc model"
        pars["disc"] = (_math("R = " + _fmt(_num(setup["r_in"])) + "-" + _fmt(_num(setup["r_out"])))
                        + ", " + _math("H/R = " + _fmt(_num(setup["hor"]))))
    elif r_cloud > 0:
        name = "Colliding-cloud model"
        pars["e"] = _math("e = " + _fmt(_num(inp.get("e_rest", 0.5))))
        pars["r_cl"] = _math(r"r_{\rm cl} = " + _fmt(r_cloud))
    elif ieos in ("3", "6"):
        name = "Isothermal gas model"
        opts = run / "moddump_m1torus.opts"
        hor = pio.read_setup(opts).get("hor") if opts.exists() else None
        pars["H/R"] = _math("H/R = " + _fmt(_num(hor))) if hor is not None else "isothermal"
    else:
        name = "Collisionless torus model"
        eps = _num(setup.get("eps_star", 0) or 0)
        if eps > 0:
            pars["eps"] = "fixed " + _math(r"\varepsilon = " + _fmt(eps))
        eps_c = _num(setup.get("eps_c", 0) or 0)
        if eps_c > 0:
            pars["eps_c"] = _math(r"\varepsilon_c = " + _fmt(eps_c))
    if geometry != 1:
        for k, v0 in CANON.items():
            v = _num(setup[k])
            if abs(v - v0) > 1e-9:
                pars[k] = _math(TEX[k] + " = " + _fmt(v) + (r"^\circ" if k == "i_max" else ""))
    suffix = ""
    if prov.get("restart"):
        num = Path(prov["restart"]).name.rsplit("_", 1)[-1]
        dump = run / f"m1torus_{num}"
        if dump.exists():
            t0 = pio.read_dump(dump)["header"]["time"] / P_ORB
            when = _math(f"t = {t0:.0f}" + r"\,P_{\rm orb}")
            what = {"Colliding-cloud model": "collisions", "Isothermal gas model": "gas"}.get(name)
            suffix = (f", {what} from {when} of the collisionless model" if what
                      else f", restarted at {when}")
    return name, pars, suffix


def model_label(run):
    name, pars, suffix = model_parts(run)
    return f"{name} ({', '.join(pars.values())}){suffix}"


def legend_labels(runs):
    """A title from what the runs share and, per run, only what differs."""
    parts = [model_parts(r) for r in runs]
    names = {p[0] for p in parts}
    keys = []
    for _, pars, _ in parts:
        keys += [k for k in pars if k not in keys]
    common = [k for k in keys if len({p[1].get(k) for p in parts}) == 1]
    title = parts[0][0] if len(names) == 1 else "Models"
    shared = [parts[0][1][k] for k in common if parts[0][1].get(k)]
    if shared:
        title += " (" + ", ".join(shared) + ")"
    labels = []
    for name, pars, suffix in parts:
        lab = ", ".join(pars[k] for k in keys if k not in common and k in pars)
        if len(names) > 1:
            lab = name + (f" ({lab})" if lab else "")
        labels.append((lab or name) + suffix)
    return title, labels


def wrapped(text, width=70):
    """Break a label into lines of about width characters at ', ' or ' (' outside $...$."""
    pieces, cur, depth = [], "", 0
    i = 0
    while i < len(text):
        ch = text[i]
        if ch == "$":
            depth ^= 1
        brk = depth == 0 and (text.startswith(", ", i) or text.startswith(" (", i))
        if brk:
            pieces.append(cur + (", " if text.startswith(", ", i) else ""))
            cur = "" if text.startswith(", ", i) else "("
            i += 2
            continue
        cur += ch
        i += 1
    pieces.append(cur)
    lines, line = [], ""
    for p in pieces:
        if line and len(line) + len(p) > width:
            lines.append(line.rstrip())
            line = ""
        line += p if not (line and p.startswith("(")) else " " + p
    lines.append(line)
    return "\n".join(lines)
