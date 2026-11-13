"""Deflection range: the largest droop and reflex that still pass every check.

    from morphgeom.deflection import find_range
    rng = find_range(cfg, af, "TE")        # or "LE"
    rng["droop"].angle_deg, rng["droop"].stopped_by, rng["droop"].res

Deflections are given as the equivalent flap angle (morph.equivalent_angles):
the angle of a straight flap, hinged at x_h, with the same tip drop.
+ = droop, - = reflex. The other edge is held unmorphed.

How it searches: it tries bigger and bigger angles in steps of
cfg.range.step_deg until a check fails (or cfg.range.max_deg is reached),
then halves the gap between the last pass and the first fail until it is
smaller than cfg.range.tol_deg. The answer is the last angle that PASSED,
so the true limit is at most tol_deg further out.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .config import copy_config
from .morph import _validate
from .pipeline import Result, load_baseline, make_shape
from .surfaces import build, loop, stations

# pass/fail checks, by the column names used in summary.csv
CHECKS = {"wingbox": "wingbox_pass", "skin_length": "arc_pass", "strain": "strain_pass",
          "slope_join": "slope_join_pass", "curvature_join": "curv_join_pass",
          "thickness": "thickness_pass", "no_self_crossing": "intersect_pass",
          "te_closure": "te_pass", "le_shape": "le_pass"}


@dataclass
class Limit:
    direction: str                 # "droop" or "reflex"
    angle_deg: float               # largest angle that passes (signed); 0 if even tol_deg fails
    delta: float                   # the same as a tip drop, fraction of chord (signed)
    stopped_by: list               # checks that fail just past it ([] = none failed up to max_deg)
    res: Optional[Result]          # the shape at the limit (None if angle_deg is 0)
    tried: list = field(default_factory=list)   # (angle_deg, summary row or None, failed checks), every shape tried


def delta_for(angle_deg: float, mp, edge: str) -> float:
    """Tip drop (fraction of chord) that gives this equivalent flap angle."""
    arm = 1 - mp.x_h if edge == "TE" else mp.x_LE
    return arm * math.tan(math.radians(angle_deg))


def failed_checks(R: dict, ignore=()) -> list:
    return [name for name, key in CHECKS.items() if name not in ignore and not R[key]]


def try_angle(cfg, af, edge: str, angle_deg: float, ignore=()):
    """Make the shape at one angle: (Result or None, list of failed checks)."""
    c = copy_config(cfg)
    c.morph.delta_TE = c.morph.delta_LE = 0.0
    setattr(c.morph, f"delta_{edge}", delta_for(angle_deg, c.morph, edge))
    try:
        res = make_shape(c, af)
    except (RuntimeError, ValueError, IndexError, ArithmeticError) as e:
        # e.g. the v2 solver cannot reach a very large deflection
        return None, [f"morph did not solve ({e})"]
    return res, failed_checks(res.checks, ignore)


def find_limit(cfg, af, edge: str, sign: int) -> Limit:
    """Largest deflection one way (sign +1 = droop, -1 = reflex) that passes."""
    rc = cfg.range
    ignore = list(rc.ignore or [])
    unknown = [k for k in ignore if k not in CHECKS]
    if unknown:
        raise KeyError(f"Unknown check name(s) in range.ignore: {unknown}. Use: {list(CHECKS)}.")
    if edge not in ("TE", "LE"):
        raise ValueError(f'edge must be "TE" or "LE" (got "{edge}").')
    c = copy_config(cfg)
    c.morph.delta_TE = c.morph.delta_LE = 0.0
    _validate(c.morph)                 # bad settings: stop now, not inside the search

    tried = []

    def attempt(a):
        res, failed = try_angle(cfg, af, edge, sign * a, ignore)
        tried.append((sign * a, None if res is None else res.row, failed))
        return res, failed

    lo, lo_res = 0.0, None             # last angle that passed
    hi, hi_failed = None, []           # first angle that failed
    k = 1
    while k * rc.step_deg <= rc.max_deg + 1e-9:
        a = k * rc.step_deg
        res, failed = attempt(a)
        if failed:
            hi, hi_failed = a, failed
            break
        lo, lo_res = a, res
        k += 1
    if hi is not None:
        while hi - lo > rc.tol_deg:
            a = (lo + hi) / 2
            res, failed = attempt(a)
            if failed:
                hi, hi_failed = a, failed
            else:
                lo, lo_res = a, res
    tried.sort(key=lambda t: t[0])
    return Limit(direction="droop" if sign > 0 else "reflex", angle_deg=sign * lo,
                 delta=delta_for(sign * lo, cfg.morph, edge), stopped_by=hi_failed, res=lo_res, tried=tried)


def find_range(cfg, af, edge: str = "TE") -> dict:
    """Largest droop and largest reflex: {"droop": Limit, "reflex": Limit}."""
    return {"droop": find_limit(cfg, af, edge, +1), "reflex": find_limit(cfg, af, edge, -1)}


def wing_box(af, x_LE: float, x_rs: float, n: int = 2001) -> dict:
    """Size of the wing box between the spars (fractions of chord):
    length, front and rear spar heights (vertical, skin to skin), and the
    box's share of the whole section area."""
    S = build(af.cam, af.thk, stations(n, "cosine", [x_LE, x_rs]))
    up, lo = S.xu > 0.02, S.xl > 0.02           # away from the nose, x runs one way along each skin

    def height(x):
        return np.interp(x, S.xu[up], S.yu[up]) - np.interp(x, S.xl[lo], S.yl[lo])

    xb = np.linspace(x_LE, x_rs, 2001)
    hb = height(xb)
    box = float(np.sum((hb[1:] + hb[:-1]) / 2 * np.diff(xb)))
    P = loop(S)
    section = 0.5 * abs(float(np.dot(P[:, 0], np.roll(P[:, 1], -1)) - np.dot(P[:, 1], np.roll(P[:, 0], -1))))
    return {"length": x_rs - x_LE, "front_height": float(height(x_LE)), "rear_height": float(height(x_rs)),
            "area_share": box / section}


def lift_per_degree(lim: Limit) -> float:
    """Rough lift change per degree of deflection, from the smallest
    deflection the search tried (the first step). Thin-airfoil theory is
    good at small deflections, so this is the trustworthy part of the
    estimate; at the limit itself the flow would separate."""
    small = [(a, r) for a, r, f in lim.tried if r is not None and a != 0]
    if not small:
        return float("nan")
    a, r = min(small, key=lambda t: abs(t[0]))
    return r["dcl_estimate"] / a


def describe(lim: Limit, max_deg: float) -> str:
    """What stopped the search, in words."""
    if not lim.stopped_by:
        return f"nothing failed up to {max_deg:g} deg"
    return ", ".join(lim.stopped_by)


# ---- one setting at a time (scripts/variable_study.py) ----

def get_setting(cfg, path: str):
    """cfg value by its dotted name, e.g. "structure.t_skin_mm"."""
    obj = cfg
    for part in path.split("."):
        if not hasattr(obj, part):
            raise KeyError(f'Unknown setting "{path}" (check the spelling against config/default.yaml).')
        obj = getattr(obj, part)
    return obj


def set_setting(cfg, path: str, value) -> None:
    """Set a cfg value by its dotted name."""
    head, _, last = path.rpartition(".")
    get_setting(cfg, path)                       # unknown name: stop with a clear message
    setattr(get_setting(cfg, head) if head else cfg, last, value)


def bend_radius_limit(cfg) -> float:
    """Tightest extra bend the skin may be given, as a fraction of chord:
    skin thickness / (2 x strain limit), over the chord. The strain check
    sees chord, skin thickness and strain limit only through this number."""
    return cfg.structure.t_skin_mm / (2 * cfg.checks.tol_strain_pct / 100) / (1000 * cfg.structure.chord_m)


def value_text(v) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if v is None:
        return "null"
    if isinstance(v, str):
        return v
    return f"{v:g}"


def _key(v):
    """Comparable form of a setting value: 3 and 3.0 match, True and 1 do not."""
    if isinstance(v, (int, float)) and not isinstance(v, bool):
        return float(v)
    return repr(v)


def vary_one(cfg, af, edge: str, entry: dict, cache: Optional[dict] = None) -> list:
    """Find the range at each value of one setting, all else as in cfg.

    entry: {"setting": "structure.t_skin_mm", "values": [1, 2, 4],
            "label": "skin thickness (mm)",       # optional, for tables and plots
            "hold": {"morph.x_rear_spar": 0.65}}  # optional, other settings fixed meanwhile
    Returns one row per value. cache (a dict) skips repeat runs of the same settings.
    """
    path, label = entry["setting"], entry.get("label", entry["setting"])
    hold = entry.get("hold") or {}
    rows = []
    for v in entry["values"]:
        c = copy_config(cfg)
        for p, hv in hold.items():
            set_setting(c, p, hv)
        set_setting(c, path, v)
        changed = tuple(sorted((p, _key(get_setting(c, p))) for p in list(hold) + [path]
                               if _key(get_setting(c, p)) != _key(get_setting(cfg, p))))
        key = (edge, changed)                    # () = the same as cfg
        if cache is not None and key in cache:
            rng = cache[key]
        else:
            rng = find_range(c, af, edge)
            if cache is not None:
                cache[key] = rng
        d, r = rng["droop"], rng["reflex"]
        rows.append({"setting": path, "label": label, "value": value_text(v),
                     "is_default": "yes" if not changed else "no", "droop_max_deg": d.angle_deg,
                     "reflex_max_deg": r.angle_deg, "range_deg": d.angle_deg - r.angle_deg,
                     "droop_tip_mm": d.delta * 1000 * c.structure.chord_m,
                     "reflex_tip_mm": r.delta * 1000 * c.structure.chord_m,
                     "lift_per_deg_estimate": lift_per_degree(d),
                     "bend_radius_limit_c": bend_radius_limit(c),
                     "droop_stopped_by": describe(d, c.range.max_deg).replace(",", ";"),
                     "reflex_stopped_by": describe(r, c.range.max_deg).replace(",", ";"),
                     "hit_search_cap": "yes" if not (d.stopped_by and r.stopped_by) else "no"})
    return rows


def unsolved(row: dict) -> bool:
    """True if the morph could not be made at any angle with this value
    (a solver limit, e.g. too few LE shape terms), so it says nothing about
    the setting's effect on the range."""
    return row["range_deg"] == 0 and row["droop_stopped_by"].startswith("morph did not solve")


def spread(rows: list) -> float:
    """How much one setting moves the total range (droop - reflex), degrees,
    over the values that could be solved."""
    v = [r["range_deg"] for r in rows if not unsolved(r)]
    return max(v) - min(v) if v else 0.0


# ---- every flap angle x hinge x shape exponent (scripts/flap_study.py) ----

_BASELINES = {}


def _baseline(cfg):
    """The baseline, loaded once per process. (Worker processes load their
    own: the loaded airfoil cannot be sent between processes.)"""
    b = cfg.baseline
    key = (b.source, b.dat_file, b.thickness_file, b.cl_design)
    if key not in _BASELINES:
        _BASELINES[key] = load_baseline(cfg)
    return _BASELINES[key]


def _short_failure(name: str) -> str:
    return "morph did not solve" if name.startswith("morph did not solve") else name


def flap_point(cfg, x_h: float, n: float, angles, af=None) -> tuple:
    """One hinge and shape exponent (rear spar at the hinge): every TE angle in
    `angles` (equivalent flap angles, deg), and the exact largest droop and
    reflex (find_range). Returns (one row per angle, one row of limits)."""
    c = copy_config(cfg)
    c.morph.x_h, c.morph.x_rear_spar, c.morph.n = x_h, None, n
    af = af if af is not None else _baseline(c)
    c_mm = 1000 * c.structure.chord_m
    rows = []
    for a in angles:
        res, failed = try_angle(c, af, "TE", a)
        r = res.row if res is not None else {}
        rows.append({"x_h": x_h, "n": float(n), "angle_deg": a, "tip_mm": delta_for(a, c.morph, "TE") * c_mm,
                     "passes": not failed, "failed": "; ".join(_short_failure(f) for f in failed),
                     "strain_max_pct": r.get("strain_max_pct", float("nan")),
                     "dcl_estimate": r.get("dcl_estimate", float("nan"))})
    rng = find_range(c, af, "TE")
    d, r = rng["droop"], rng["reflex"]
    limits = {"x_h": x_h, "n": float(n), "droop_max_deg": d.angle_deg, "reflex_max_deg": r.angle_deg,
              "range_deg": d.angle_deg - r.angle_deg, "droop_tip_mm": d.delta * c_mm, "reflex_tip_mm": r.delta * c_mm,
              "lift_per_deg_estimate": lift_per_degree(d),
              "droop_stopped_by": "; ".join(_short_failure(f) for f in d.stopped_by) or "nothing (search cap)",
              "reflex_stopped_by": "; ".join(_short_failure(f) for f in r.stopped_by) or "nothing (search cap)"}
    return rows, limits
