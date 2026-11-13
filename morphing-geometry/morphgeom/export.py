"""Output files and reports: Selig .dat, CAD XYZ (mm), summary rows / CSV, printed report."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from .morph import rear_spar


def shape_name(base_name: str, mp) -> str:
    """File name stem that spells out the morph settings.

        naca641412_TE+0.020_xh0.650_n3
        naca641412_LE+0.010_xle0.150_TE+0.020_xh0.650_n3   (both edges)
        naca641412_baseline                                (no morph)
        naca641412_TE+0.020_xh0.650_n3_constlen            (v2 method)
        naca641412_TE+0.020_xh0.700_rs0.600_n3             (hinge behind the rear spar)
    A part is left out when its delta is 0. rs (rear spar) only appears when
    the spar is set apart from the hinge.
    """
    parts = [base_name]
    if mp.delta_LE != 0:
        parts.append(f"LE{mp.delta_LE:+.3f}_xle{mp.x_LE:.3f}")
    if mp.delta_TE != 0:
        parts.append(f"TE{mp.delta_TE:+.3f}_xh{mp.x_h:.3f}")
        if rear_spar(mp) != mp.x_h:
            parts.append(f"rs{rear_spar(mp):.3f}")
    if len(parts) == 1:
        parts.append("baseline")
    else:
        parts.append(f"n{mp.n:g}")
        if getattr(mp, "method", "camber") == "constant_length":
            parts.append("constlen")
    return "_".join(parts)


def _clean(a: np.ndarray) -> np.ndarray:
    a = np.round(a, 6)
    a[a == 0] = 0.0          # no "-0.000000"
    return a


def write_selig(path, xy: np.ndarray, name: str) -> None:
    """Write a coordinate loop as a Selig-style .dat file.

    Line 1: name. Then one "x y" pair per line, chord = 1, 6 decimals, in the
    order given (surfaces.loop makes TE -> upper -> LE -> lower -> TE).
    XFOIL and most meshers read this directly.
    """
    xy = _clean(np.array(xy, float))
    with open(path, "w") as f:
        f.write(f"{name}\n")
        for x, y in xy:
            f.write(f"{x:10.6f} {y:10.6f}\n")


def write_cad_xyz(path, xy: np.ndarray, chord_m: float) -> None:
    """Write a coordinate loop as X Y Z in millimetres for CAD.

    One ordered loop (same order as the .dat), scaled by the chord, Z = 0,
    tab separated, no header. In SolidWorks: Insert > Curve > Curve Through
    XYZ Points > Browse (set the part units to mm first). For a closed
    trailing edge the first and last points are the same, so the curve
    closes. CATIA and most CAD tools read the same three-column text.
    """
    P = _clean(np.array(xy, float) * chord_m * 1000)
    with open(path, "w") as f:
        for x, y in P:
            f.write(f"{x:.6f}\t{y:.6f}\t{0.0:.6f}\n")


def summary_row(res) -> dict:
    """One flat row of settings, check values and pass/fail for a shape.

    Units are picked so three decimals are meaningful: mm, percent, degrees,
    1/m. Pass/fail fields are bool (written as PASS/FAIL). Counts are int.
    """
    mp, R = res.morph, res.checks
    row = {
        "name": res.name, "method": res.mr.method,
        "delta_TE": mp.delta_TE, "delta_LE": mp.delta_LE, "x_h": mp.x_h, "x_rear_spar": rear_spar(mp),
        "x_LE": mp.x_LE, "n": float(mp.n),
        "te_angle_deg": res.angles["te_deg"], "le_angle_deg": res.angles["le_deg"],
        "wingbox_max_mm": R["wingbox_max_mm"], "wingbox": R["wingbox_pass"],
    }
    for k in ("te_upper", "te_lower", "le_upper", "le_lower", "le_wrap"):
        row[f"arc_{k}_pct"] = R[f"arc_{k}_pct"]
    for k in ("te_upper", "te_lower", "le_upper", "le_lower", "le_wrap"):
        row[f"arc_{k}_mm"] = R[f"arc_{k}_mm"]
    row.update({
        "arc_upper_total_pct": R["arc_upper_total_pct"], "arc_lower_total_pct": R["arc_lower_total_pct"],
        "arc_max_pct": R["arc_max_pct"], "skin_length": R["arc_pass"],
        "te_lower_slide_mm": R["te_lower_slide_mm"],
        "strain_te_upper_pct": R["strain_te_upper_pct"], "strain_te_lower_pct": R["strain_te_lower_pct"],
        "strain_le_upper_pct": R["strain_le_upper_pct"], "strain_le_lower_pct": R["strain_le_lower_pct"],
        "strain_max_pct": R["strain_max_pct"], "strain": R["strain_pass"],
        "slope_jump_xh_deg": R["slope_jump_xh_deg"], "slope_jump_xle_deg": R["slope_jump_xle_deg"],
        "curv_jump_xh_per_m": R["curv_jump_xh_per_m"], "curv_jump_xle_per_m": R["curv_jump_xle_per_m"],
        "slope_join": R["slope_join_pass"], "curvature_join": R["curv_join_pass"],
        "skin_slope_jump_deg": R["skin_slope_jump_deg"], "skin_curv_jump_per_m": R["skin_curv_jump_per_m"],
        "min_thickness_pct": R["min_thickness_pct"], "thickness": R["thickness_pass"],
        "thickness_change_max_mm": R["thickness_change_max_mm"],
        "self_intersections": int(R["self_intersections"]), "no_self_crossing": R["intersect_pass"],
        "te_gap_change_mm": R["te_gap_change_mm"], "te_closure": R["te_pass"],
        "le_new_inflections": int(max(0, R["le_inflections_morph"] - R["le_inflections_base"])),
        "nose_radius_base_mm": R["nose_radius_base_mm"], "nose_radius_morph_mm": R["nose_radius_morph_mm"],
        "le_shape": R["le_pass"],
        "te_move_x_mm": R["te_move_x_mm"], "te_move_y_mm": R["te_move_y_mm"],
        "nose_drop_mm": R["nose_drop_mm"], "off_v1_mm": R["off_v1_mm"], "dcl_estimate": R["dcl_estimate"],
        "all_checks": R["all_pass"],
    })
    return row


def _cell(v) -> str:
    if isinstance(v, (bool, np.bool_)):
        return "PASS" if v else "FAIL"
    if isinstance(v, str):
        return v
    if isinstance(v, (int, np.integer)):
        return f"{v:d}"
    if math.isnan(v):
        return "NaN"
    s = f"{v:.3f}"
    return "0.000" if s == "-0.000" else s


def write_summary_csv(path, rows: list[dict]) -> None:
    """Write summary rows (same keys) to a CSV file.
    Numbers: 3 decimals. Pass/fail: PASS or FAIL. Counts: whole numbers."""
    names = list(rows[0].keys())
    with open(path, "w") as f:
        f.write(",".join(names) + "\n")
        for r in rows:
            f.write(",".join(_cell(r[k]) for k in names) + "\n")


def print_report(res, cfg) -> None:
    """Print one shape's settings and check results, 3 decimals."""
    r, ck = res.row, cfg.checks
    v2 = r["method"] == "constant_length"

    def pf(ok):
        return "PASS" if ok else "FAIL"

    def show(label, value, limit=math.nan, result=""):
        lim = "" if math.isnan(limit) else f"{limit:.3f}"
        v = f"{value:.3f}"
        v = "0.000" if v == "-0.000" else v
        print(f"  {label:<36s} {v:>12s} {lim:>12s}   {result}")

    print(f"\n{r['name']}")
    if v2:
        print("  method: constant_length (v2) - skins bend but keep their length; lower TE skin slides")
    else:
        print("  method: camber (v1) - camber line shifted, thickness kept; skins may stretch")
    spar = "" if r["x_rear_spar"] == r["x_h"] else f", rear spar {r['x_rear_spar']:.3f}"
    print(f"  TE: delta = {r['delta_TE']:+.3f} c ({r['te_angle_deg']:.3f} deg), hinge x_h = {r['x_h']:.3f}{spar}"
          f"   |   LE: delta = {r['delta_LE']:+.3f} c ({r['le_angle_deg']:.3f} deg), front spar x_LE = {r['x_LE']:.3f}"
          f"   |   n = {r['n']:g}")
    print(f"  chord = {cfg.structure.chord_m:.3f} m, skin = {cfg.structure.t_skin_mm:.3f} mm (placeholders unless set)")
    print(f"  {'check':<36s} {'value':>12s} {'limit':>12s}   result")
    show("wing box: max skin movement (mm)", r["wingbox_max_mm"], ck.tol_wingbox_mm, pf(r["wingbox"]))
    show("skin length change TE upper (%)", r["arc_te_upper_pct"], ck.tol_arc_pct, pf(r["skin_length"]))
    if v2:
        show("                   TE lower (%)", r["arc_te_lower_pct"], result="slides")
    else:
        show("                   TE lower (%)", r["arc_te_lower_pct"], ck.tol_arc_pct)
    show("                   LE upper (%)", r["arc_le_upper_pct"], ck.tol_arc_pct)
    show("                   LE lower (%)", r["arc_le_lower_pct"], ck.tol_arc_pct)
    show("                   LE wrap  (%)", r["arc_le_wrap_pct"], result="info")
    show("bending strain TE upper (%)", r["strain_te_upper_pct"], ck.tol_strain_pct, pf(r["strain"]))
    show("               TE lower (%)", r["strain_te_lower_pct"], ck.tol_strain_pct)
    show("               LE upper (%)", r["strain_le_upper_pct"], ck.tol_strain_pct)
    show("               LE lower (%)", r["strain_le_lower_pct"], ck.tol_strain_pct)
    what = "skin" if v2 else "camber"      # v2 is defined on the upper TE skin and the LE skin
    show(f"{what} slope jump, hinge (deg)", r["slope_jump_xh_deg"], ck.tol_slope_deg, pf(r["slope_join"]))
    show("                  front spar (deg)", r["slope_jump_xle_deg"], ck.tol_slope_deg)
    show(f"{what} curv. jump, hinge (1/m)", r["curv_jump_xh_per_m"], ck.tol_curv_per_m, pf(r["curvature_join"]))
    show("                  front spar (1/m)", r["curv_jump_xle_per_m"], ck.tol_curv_per_m)
    show("any skin slope jump, max (deg)", r["skin_slope_jump_deg"], result="info")
    show("any skin curvature jump, max (1/m)", r["skin_curv_jump_per_m"], result="info")
    show("min half-thickness (% chord)", r["min_thickness_pct"], result=pf(r["thickness"]))
    show("self-crossings (count)", r["self_intersections"], result=pf(r["no_self_crossing"]))
    show("TE gap change (mm)", r["te_gap_change_mm"], ck.tol_te_gap_mm, pf(r["te_closure"]))
    show("new LE inflection points (count)", r["le_new_inflections"], result=pf(r["le_shape"]))
    print("  info")
    print(f"  {'nose radius':<36s} {r['nose_radius_base_mm']:12.3f} -> {r['nose_radius_morph_mm']:.3f} mm")
    show("nose drop at 1% chord (mm)", r["nose_drop_mm"])
    show("TE moves: x (mm, - = forward)", r["te_move_x_mm"])
    show("          y (mm, - = down)", r["te_move_y_mm"])
    show("max thickness change (mm)", r["thickness_change_max_mm"])
    show("lift change, rough estimate (cl)", r["dcl_estimate"])
    if v2:
        show("lower TE skin slides into TE (mm)", r["te_lower_slide_mm"])
        show("max distance from v1 shape (mm)", r["off_v1_mm"])
    print(f"  ALL CHECKS: {pf(r['all_checks'])}")
