"""The shapes to send to CFD: which ones, their coordinate files, and a note on how to run them.

    from morphgeom.cfd import cfd_cases, write_cfd_set
    cases, skipped = cfd_cases(cfg)       # the baseline first, then every group in cfg.cfd_set
    write_cfd_set(cases, skipped, cfg, outdir)

Only the trailing edge morphs, with the rear spar at the hinge (as in flap_study.py).
Angles are equivalent flap angles (deg): + = droop, - = reflex.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from . import export, plots
from .aero import flight_condition, lift_per_span
from .config import copy_config
from .deflection import _baseline, _short_failure, delta_for, try_angle
from .pipeline import make_shape
from .surfaces import loop


def case_name(af_name: str, case: int, angle_deg=None, x_h=None, n=None, v2=False) -> str:
    """File name stem: sc20714_cfd01_baseline, sc20714_cfd07_flap+5deg_xh0.650_n3_constlen."""
    if angle_deg is None:
        return f"{af_name}_cfd{case:02d}_baseline"
    return f"{af_name}_cfd{case:02d}_flap{angle_deg:+g}deg_xh{x_h:.3f}_n{n:g}" + ("_constlen" if v2 else "")


def cfd_cases(cfg, af=None) -> tuple:
    """Every shape in cfg.cfd_set, in run order. Returns (cases, skipped).

    cases    one dict per shape that passes every check: case (1, 2, ...), group,
             name, x_h, n, angle_deg, delta (tip drop, fraction of chord) and res
             (the Result). Case 1 is the baseline (x_h, n and angle_deg None).
    skipped  shapes that fail a check: the same keys without case and name, plus
             failed (the checks that fail).
    A shape listed twice (same hinge, n and angle) is kept once, in its first
    group. Angle 0 is the baseline, so it is not repeated.
    """
    c = copy_config(cfg)
    c.morph.method = cfg.cfd_set.method
    c.morph.x_rear_spar = None
    c.morph.delta_TE = c.morph.delta_LE = 0.0
    af = af if af is not None else _baseline(c)
    v2 = c.morph.method == "constant_length"
    cases = [{"case": 1, "group": "baseline", "name": case_name(af.name, 1), "x_h": None, "n": None,
              "angle_deg": None, "delta": 0.0, "res": make_shape(c, af)}]
    skipped, seen = [], set()
    for g in cfg.cfd_set.groups:
        for x_h in g["x_h"]:
            for n in g["n"]:
                for a in g["angles"]:
                    x_h, n, a = float(x_h), float(n), float(a)
                    key = (round(x_h, 6), round(n, 6), round(a, 6))
                    if a == 0 or key in seen:
                        continue
                    seen.add(key)
                    s = copy_config(c)
                    s.morph.x_h, s.morph.n = x_h, n
                    res, failed = try_angle(s, af, "TE", a)
                    entry = {"group": g["name"], "x_h": x_h, "n": n, "angle_deg": a,
                             "delta": delta_for(a, s.morph, "TE"), "res": res}
                    if failed:
                        entry["failed"] = [_short_failure(f) for f in failed]
                        skipped.append(entry)
                    else:
                        entry["case"] = len(cases) + 1
                        entry["name"] = case_name(af.name, entry["case"], a, x_h, n, v2)
                        cases.append(entry)
    return cases, skipped


def case_rows(cases, cfg) -> list:
    """One row per case for cfd_cases.csv (and the printed table)."""
    fc = flight_condition(cfg.flight)
    chord = cfg.structure.chord_m
    rows = []
    for c in cases:
        r = c["res"].row
        base = c["angle_deg"] is None
        rows.append({
            "case": c["case"], "group": c["group"], "file": f"{c['name']}.dat",
            "hinge_x_c": "-" if base else c["x_h"], "n": "-" if base else c["n"],
            "flap_angle_deg": 0.0 if base else c["angle_deg"],
            "tip_drop_pct_chord": 100 * c["delta"],
            "strain_max_pct": r["strain_max_pct"],
            "lift_change_cl_estimate": r["dcl_estimate"],
            "lift_change_n_per_m_estimate": lift_per_span(r["dcl_estimate"], fc["q_pa"], chord),
        })
    return rows


def _source(cfg) -> str:
    b = cfg.baseline
    if b.source == "dat":
        return f"coordinates from {Path(b.dat_file).name}"
    return f"NACA 6-series thickness ({Path(b.thickness_file).name}) + a = 1.0 mean line, design cl {b.cl_design:g}"


def _list(values, fmt) -> str:
    return ", ".join(fmt(float(v)) for v in values)


def note(cases, skipped, cfg) -> str:
    """The text of README_cfd.txt: how the shapes are set up, the flight
    condition, what is still open and what to send back."""
    fl, st, pts = cfg.flight, cfg.structure, cfg.points
    fc = flight_condition(fl)
    base = cases[0]["res"].base
    gap = float(np.hypot(base.xu[-1] - base.xl[-1], base.yu[-1] - base.yl[-1]))
    v2 = cfg.cfd_set.method == "constant_length"
    order = "TE -> upper surface -> LE -> lower surface -> TE (Selig)" if pts.order == "selig" \
        else "TE -> lower surface -> LE -> upper surface -> TE"
    gap_rho = fc["density_from_p_t"] / fl.density_kg_m3 - 1
    L = [
        f"Shapes for CFD: {cfg.baseline.name} baseline + {len(cases) - 1} morphed trailing edges",
        f"Made by scripts/cfd_shapes.py from the cfd_set settings (config/default.yaml). Baseline: {_source(cfg)}.",
        "",
        "FILES",
        "  <name>.dat       one shape each. Line 1 is the name, then one x y pair per line.",
        f"                   Chord = 1, 6 decimals, {pts.n_per_surface} points per surface"
        f" ({pts.spacing} spacing{', bunched at the LE and TE' if pts.spacing == 'cosine' else ''}).",
        f"                   Order: {order}.",
        "  cfd_cases.csv    every case: group, hinge, shape exponent n, flap angle, tip drop, rough lift estimate",
        "  cfd_shapes.png   every shape drawn against the baseline",
        "",
        "HOW THE SHAPES ARE SET UP",
        "  - Only the trailing edge changes, behind the hinge (\"xh\" in the file name, x/c). Everything ahead",
        "    of the hinge is the same as the baseline. The rear spar is at the hinge.",
        "  - The flap angle (\"flap+5deg\") is the equivalent flap angle: a plain flap hinged at the same point",
        "    with the same trailing-edge drop would sit at this angle. + = trailing edge down (more lift), - = up.",
        "  - n is the shape exponent: how the bend is spread between the hinge and the trailing edge (the drop",
        "    grows as distance^n). Higher n keeps the skin straighter near the hinge and bends it more near the tip.",
        ("  - The skins bend but keep their length (the lower skin slides a few mm into the trailing edge)."
         if v2 else "  - The camber line is shifted and the thickness kept (the skins stretch slightly)."),
        "  - The shapes are NOT re-normalised. The baseline chord line stays on the x-axis and a drooped",
        "    trailing edge simply sits lower. Run every shape at the same angle of attack, measured from the",
        "    x-axis, so they are all compared at the same wing attitude.",
    ]
    if gap > 1e-6:
        L += [f"  - Blunt trailing edge: {100 * gap:.2f}% chord thick ({gap * st.chord_m * 1000:.1f} mm at the placeholder "
              f"{st.chord_m:g} m chord).",
              "    The first and last points are the top and bottom of the base; they are not joined."]
    else:
        L += ["  - Sharp trailing edge: the first and last points are the same point."]
    L += [
        "",
        "FLIGHT CONDITION (config/default.yaml, flight:)",
        f"  Altitude {fl.altitude_ft:,.0f} ft. Velocity {fl.velocity_m_s:g} m/s, speed of sound {fl.speed_of_sound_m_s:g} m/s:"
        f" Mach {fc['mach']:.3f}.",
        f"  Temperature {fc['temperature_k']:.1f} K (from the speed of sound). Pressure {fl.pressure_pa:,.0f} Pa."
        f" Density {fl.density_kg_m3:g} kg/m3 (as given).",
        f"  Reynolds number {fc['reynolds_per_m'] / 1e6:.2f} million per metre of chord"
        f" ({fc['reynolds_per_m'] * st.chord_m / 1e6:.2f} million at the placeholder {st.chord_m:g} m chord).",
        f"  Dynamic pressure {fc['q_pa']:,.0f} Pa.",
        "",
        "STILL OPEN (question numbers from the README)",
        f"  - Q18 wing sweep: the section sees about Mach {fc['mach']:.2f} x cos(sweep). Run at that Mach if the wing is swept.",
    ]
    if abs(gap_rho) > 0.02:
        L += [f"  - Q17 density: the pressure and speed of sound give {fc['density_from_p_t']:.3f} kg/m3, not "
              f"{fl.density_kg_m3:g}. That would change the Reynolds number by {100 * gap_rho:+.0f}%."]
    L += [
        f"  - Q4 chord: {st.chord_m:g} m is a placeholder, and the Reynolds number scales with it.",
        "  - Q5 point format and how to close or mesh the trailing edge.",
        "  - Angle of attack: to agree with aero (the cruise angle of attack, or a few angles around it).",
        "",
        "PLEASE SEND BACK, for every shape and angle of attack",
        "  - lift, drag and pitching moment coefficients (cl, cd, cm about the quarter chord)",
        "  - the pressure distribution, where the shock sits, and whether the flow separates",
        "",
        "THE GROUPS",
    ]
    for g in cfg.cfd_set.groups:
        ids = [c["case"] for c in cases if c["group"] == g["name"]]
        L += [f"  {g['name']} (cases {_ranges(ids)}): hinge {_list(g['x_h'], lambda v: f'{v:.3f}')}, "
              f"n {_list(g['n'], lambda v: f'{v:g}')}, flap angles {_list(g['angles'], lambda v: f'{v:+g}')} deg"]
        if g.get("why"):
            L += [f"    to show {g['why']}"]
    if skipped:
        L += ["", "LEFT OUT (fail a buildability check with the current placeholder skin, chord and strain limit)"]
        L += [f"  {s['group']}: {s['angle_deg']:+g} deg, hinge {s['x_h']:.3f}, n {s['n']:g} ({', '.join(s['failed'])})"
              for s in skipped]
    L += ["", "The lift numbers in cfd_cases.csv are a rough low-speed estimate (thin-airfoil theory). They are not",
          "reliable at this Mach number: replacing them is the point of these runs.", ""]
    return "\n".join(L)


def _ranges(ids) -> str:
    """[2, 3, 4, 7] -> '2-4, 7'."""
    if not ids:
        return "none"
    out, start = [], ids[0]
    for a, b in zip(ids, ids[1:] + [None]):
        if b != a + 1:
            out.append(f"{start}" if start == a else f"{start}-{a}")
            start = b
    return ", ".join(out)


def write_cfd_set(cases, skipped, cfg, outdir) -> dict:
    """Write the CFD folder. Old <airfoil>_cfd*.dat files in it are removed
    first, so it only ever holds the current set."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    af_name = cases[0]["name"].rsplit("_cfd", 1)[0]
    for old in outdir.glob(f"{af_name}_cfd*.dat"):
        old.unlink()
    for c in cases:
        export.write_selig(outdir / f"{c['name']}.dat", loop(c["res"].shape, cfg.points.order), c["name"])
    files = {"csv": outdir / "cfd_cases.csv", "png": outdir / "cfd_shapes.png", "note": outdir / "README_cfd.txt"}
    export.write_summary_csv(files["csv"], case_rows(cases, cfg))
    plots.cfd_shapes(cases, f"Shapes for CFD - {af_name}", files["png"])
    files["note"].write_text(note(cases, skipped, cfg))
    return files

