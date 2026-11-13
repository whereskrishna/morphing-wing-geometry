"""Automatic checks: wing box, skin length, bending strain, smooth joins, sanity."""
from __future__ import annotations

import numpy as np

from .aero import lift_change
from .morph import Morph, deflect_camber, rear_spar
from .surfaces import Airfoil, build, loop, point_to_polyline, stations


def curvature(xy: np.ndarray) -> np.ndarray:
    """Signed curvature at each point of a polyline, from its neighbours.

    Uses the circle through each point and its two neighbours. Positive =
    the path turns left (counter-clockwise). An airfoil loop in Selig order
    runs counter-clockwise, so on it convex skin = positive curvature.
    The first and last points have only one neighbour and get nan.
    """
    a, b, c = xy[:-2], xy[1:-1], xy[2:]
    ab, bc, ac = b - a, c - b, c - a
    cross2 = ab[:, 0] * bc[:, 1] - ab[:, 1] * bc[:, 0]
    k = 2 * cross2 / (np.linalg.norm(ab, axis=1) * np.linalg.norm(bc, axis=1) * np.linalg.norm(ac, axis=1))
    return np.concatenate([[np.nan], k, [np.nan]])


def bending_strain(t_skin_mm, dkappa_per_m):
    """Skin bending strain from a change in curvature: eps = (t/2) * dkappa.

    The outer fibre of a skin of thickness t sits t/2 from its middle, so
    bending it by dkappa stretches that fibre by (t/2) * dkappa.

    Worked example: a 2.000 mm skin bent from flat to a 0.500 m radius:
    dkappa = 1/0.500 = 2.000 1/m, eps = 0.001 m * 2.000 1/m = 0.002 (0.200%).
    """
    return (t_skin_mm / 1000 / 2) * np.asarray(dkappa_per_m)


def joint_jumps(shape, x0: float, h: float = 1e-4) -> dict:
    """Jump in slope and in curvature across station x0 (a spar).

    shape  function x -> Surfaces (e.g. lambda x: build(cam, thk, x), or the
           eval of a morph from morph.apply)

    Measured numerically on three curves: the camber line (x, yc), the
    upper surface and the lower surface. On each side of x0 the slope angle
    and the curvature are estimated with one-sided, second-order differences
    that use only points on that side. The jump is right minus left.

    slope_deg  |jump in slope angle|, degrees  [camber, upper, lower]
    curv       |jump in curvature|, 1/chord    [camber, upper, lower]
    """
    k = np.arange(4)
    SL = shape(x0 - h * k)        # left: x0, x0-h, ...
    SR = shape(x0 + h * k)        # right: x0, x0+h, ...
    curves = [("x", "yc"), ("xu", "yu"), ("xl", "yl")]
    slope, curv = np.zeros(3), np.zeros(3)
    for i, (cx, cy) in enumerate(curves):
        aL, kL = _one_sided(getattr(SL, cx), getattr(SL, cy), -h)
        aR, kR = _one_sided(getattr(SR, cx), getattr(SR, cy), h)
        slope[i] = abs(np.degrees(aR - aL))
        curv[i] = abs(kR - kL)
    return {"slope_deg": slope, "curv": curv}


def _one_sided(X, Y, h):
    """Derivatives at the first point from points at steps 0, h, 2h, 3h
    (h may be negative). Second-order accurate."""
    def d1(f):
        return (-3 * f[0] + 4 * f[1] - f[2]) / (2 * h)

    def d2(f):
        return (2 * f[0] - 5 * f[1] + 4 * f[2] - f[3]) / h ** 2

    x1, y1, x2, y2 = d1(X), d1(Y), d2(X), d2(Y)
    return np.arctan2(y1, x1), (x1 * y2 - y1 * x2) / (x1 ** 2 + y1 ** 2) ** 1.5


def self_intersections(xy: np.ndarray) -> int:
    """Number of places where a point loop crosses itself.

    Tests every pair of segments that do not share a point. Segments that
    only touch are not counted. If the first and last points are the same
    (closed trailing edge) the first and last segments count as neighbours.
    """
    P, Q = xy[:-1], xy[1:]
    m = P.shape[0]
    closed = bool(np.array_equal(xy[0], xy[-1]))
    lo, hi = np.minimum(P, Q), np.maximum(P, Q)

    def orient(a, b, c):
        # which side of line a->b the point(s) c are on (sign of the cross product)
        return (b[..., 0] - a[..., 0]) * (c[..., 1] - a[..., 1]) - (b[..., 1] - a[..., 1]) * (c[..., 0] - a[..., 0])

    n = 0
    for i in range(m - 2):
        j = np.arange(i + 2, m)
        if closed and i == 0:
            j = j[j != m - 1]
        # quick reject: bounding boxes must overlap
        j = j[(hi[j, 0] >= lo[i, 0]) & (lo[j, 0] <= hi[i, 0]) & (hi[j, 1] >= lo[i, 1]) & (lo[j, 1] <= hi[i, 1])]
        if j.size == 0:
            continue
        d1 = orient(P[i], Q[i], P[j])
        d2 = orient(P[i], Q[i], Q[j])
        d3 = orient(P[j], Q[j], P[i])
        d4 = orient(P[j], Q[j], Q[i])
        n += int(np.count_nonzero((d1 * d2 < 0) & (d3 * d4 < 0)))
    return n


def run_all(af: Airfoil, mr: Morph, cfg) -> dict:
    """Every automatic check for one morphed shape: a value and pass/fail for each.

    1. Wing box unchanged   largest skin movement between the spars, mm
    2. Skin length          % change of each skin piece
    3. Bending strain       (t/2) x curvature change, matched by arc-length fraction
    4. Smooth joins         slope and curvature jumps at the front spar and the hinge
    5. Sanity               thickness >= 0, no self-crossing, TE gap kept,
                            no new wiggles (inflection points) in the LE skin
    Info (no pass/fail): thickness change, TE and nose movement, nose radius,
    a rough lift change (thin-airfoil estimate, aero.py), and for v2 how far
    the shape is from the v1 shape and how far the lower TE skin slides.

    Skin pieces (one per surface per morph region): TE upper and TE lower
    (hinge -> TE), LE upper and LE lower (front spar -> nose tip). The
    whole LE skin round the nose ("LE wrap", spar to spar) is also reported,
    because constant-arc-length LE designs hold that total fixed. In v2 the
    lower TE skin is designed to slide into the TE, so its length change is
    reported as the slide, not judged.

    The checks use their own dense points (cfg.checks.n_points per surface,
    cosine spacing, with a station exactly at each spar), so the results do
    not depend on how many points go into the output files.
    R['dense'] keeps those points and curvatures for the plots.
    """
    mp = mr.morph
    xrs = rear_spar(mp)
    v2 = mr.method == "constant_length"
    ck = cfg.checks
    c_m = cfg.structure.chord_m
    c_mm = 1000 * c_m
    R = {}

    x = stations(ck.n_points, "cosine", [mp.x_LE, mp.x_h])
    N = x.size
    ile = int(np.nonzero(x == mp.x_LE)[0][0])
    ih = int(np.nonzero(x == mp.x_h)[0][0])
    B = build(af.cam, af.thk, x)
    M = mr.eval(x)
    loopB, loopM = loop(B), loop(M)
    kB, kM = curvature(loopB), curvature(loopM)      # 1/chord, convex = positive

    # Loop order is TE -> upper -> LE -> lower -> TE, so station i sits at
    # loop index N-1-i on the upper surface and N-1+i on the lower surface.
    station = np.concatenate([np.arange(N - 1, -1, -1), np.arange(1, N)])
    pieces = {"te_upper": np.arange(0, N - ih),
              "te_lower": np.arange(N - 1 + ih, 2 * N - 1),
              "le_upper": np.arange(N - 1 - ile, N),
              "le_lower": np.arange(N - 1, N + ile)}
    wrap = np.arange(N - 1 - ile, N + ile)             # LE skin, upper spar -> nose -> lower spar

    # 1. Wing box unchanged (front spar to rear spar)
    dU = np.hypot(M.xu - B.xu, M.yu - B.yu)
    dL = np.hypot(M.xl - B.xl, M.yl - B.yl)
    box = (x >= mp.x_LE) & (x <= xrs)
    R["wingbox_station_max"] = float(max(dU[box].max(), dL[box].max()))   # fraction of chord; should be exactly 0
    # A few skin points just outside the spar stations sit physically inside
    # the box (the thickness is not laid vertically). Count them too.
    inU = (B.xu >= mp.x_LE) & (B.xu <= xrs)
    inL = (B.xl >= mp.x_LE) & (B.xl <= xrs)
    R["wingbox_max_mm"] = float(max(dU[inU].max(), dL[inL].max())) * c_mm
    R["wingbox_pass"] = R["wingbox_station_max"] <= 1e-12 and R["wingbox_max_mm"] <= ck.tol_wingbox_mm

    # 2. Skin length
    def length(P, idx):
        return float(np.sum(np.linalg.norm(np.diff(P[idx], axis=0), axis=1)))

    def pct(idx):
        return 100 * (length(loopM, idx) - length(loopB, idx)) / length(loopB, idx)

    def dmm(idx):
        return (length(loopM, idx) - length(loopB, idx)) * c_mm

    R["arc_upper_total_pct"] = pct(np.arange(0, N))
    R["arc_lower_total_pct"] = pct(np.arange(N - 1, 2 * N - 1))
    R["arc_le_wrap_pct"] = pct(wrap)
    R["arc_le_wrap_mm"] = dmm(wrap)
    worst = 0.0
    for name, idx in pieces.items():
        R[f"arc_{name}_pct"] = pct(idx)
        R[f"arc_{name}_mm"] = dmm(idx)
        if not (v2 and name == "te_lower"):          # v2: this skin slides
            worst = max(worst, abs(pct(idx)))
    R["arc_max_pct"] = worst
    R["arc_pass"] = worst <= ck.tol_arc_pct
    R["te_lower_slide_mm"] = -R["arc_te_lower_mm"] if v2 else float("nan")   # + = slides into the TE

    # 3. Bending strain
    # Each skin point keeps its fraction of the piece's length (it can bend but
    # barely stretch). The curvature change is taken between matching points.
    # The last cfg.checks.te_exclude of chord (the solid TE strip) is skipped:
    # the a = 1.0 mean line's curvature grows without limit in the last ~1% of
    # chord. That is a property of the formula, not of a buildable skin.
    tip = x[station] > 1 - ck.te_exclude
    emax = 0.0
    for name, idx in pieces.items():
        fb, fm = _cumlen(loopB[idx]), _cumlen(loopM[idx])
        dk = np.interp(fb, fm, kM[idx]) - kB[idx]        # 1/chord
        use = np.isfinite(dk) & ~tip[idx]
        eps_pct = 100 * np.abs(bending_strain(cfg.structure.t_skin_mm, dk[use] / c_m))
        j = int(np.argmax(eps_pct))
        R[f"strain_{name}_pct"] = float(eps_pct[j])
        R[f"strain_{name}_at_x"] = float(x[station[idx[use][j]]])
        emax = max(emax, float(eps_pct[j]))
    R["strain_max_pct"] = emax
    R["strain_pass"] = emax <= ck.tol_strain_pct

    # 4. Smooth joins at the front spar and the hinge (where the TE starts to bend)
    # Pass/fail is on the curves the morph is defined on, where the n rule
    # holds exactly (n >= 2: slope matches, n >= 3: curvature matches):
    #   v1: the camber line;  v2: the upper TE skin and the LE skin.
    # The other skins follow those curves at a distance (thickness, ribs) and
    # need one order more; their jumps are reported as info.
    Jh = joint_jumps(mr.eval, mp.x_h)
    Jle = joint_jumps(mr.eval, mp.x_LE)
    drive_h, drive_le = ([1], [1, 2]) if v2 else ([0], [0])
    R["slope_jump_xh_deg"] = float(Jh["slope_deg"][drive_h].max())
    R["slope_jump_xle_deg"] = float(Jle["slope_deg"][drive_le].max())
    R["curv_jump_xh_per_m"] = float(Jh["curv"][drive_h].max()) / c_m
    R["curv_jump_xle_per_m"] = float(Jle["curv"][drive_le].max()) / c_m
    R["skin_slope_jump_deg"] = float(np.max(np.concatenate([Jh["slope_deg"][1:], Jle["slope_deg"][1:]])))
    R["skin_curv_jump_per_m"] = float(np.max(np.concatenate([Jh["curv"][1:], Jle["curv"][1:]]))) / c_m
    R["slope_join_pass"] = max(R["slope_jump_xh_deg"], R["slope_jump_xle_deg"]) <= ck.tol_slope_deg
    R["curv_join_pass"] = max(R["curv_jump_xh_per_m"], R["curv_jump_xle_per_m"]) <= ck.tol_curv_per_m

    # 5. Sanity
    tB = np.hypot(B.xu - B.xl, B.yu - B.yl)    # local thickness: upper to lower point at each station
    tM = np.hypot(M.xu - M.xl, M.yu - M.yl)
    R["min_thickness_pct"] = 100 * float(tM.min()) / 2   # as a half-thickness, like the NACA tables
    R["thickness_pass"] = bool(tM.min() >= 0)
    R["thickness_change_max_mm"] = float(np.max(np.abs(tM - tB))) * c_mm

    R["self_intersections"] = self_intersections(loopM)
    R["intersect_pass"] = R["self_intersections"] == 0

    gapB = np.hypot(B.xu[-1] - B.xl[-1], B.yu[-1] - B.yl[-1])
    gapM = np.hypot(M.xu[-1] - M.xl[-1], M.yu[-1] - M.yl[-1])
    R["te_gap_change_mm"] = float(abs(gapM - gapB)) * c_mm
    R["te_pass"] = R["te_gap_change_mm"] <= ck.tol_te_gap_mm

    R["le_inflections_base"] = _sign_changes(kB[wrap])
    R["le_inflections_morph"] = _sign_changes(kM[wrap])
    R["nose_radius_base_mm"] = c_mm / float(np.nanmax(kB[wrap]))
    R["nose_radius_morph_mm"] = c_mm / float(np.nanmax(kM[wrap]))
    R["le_pass"] = R["le_inflections_morph"] <= R["le_inflections_base"]

    # Info: where things moved
    R["te_move_x_mm"] = float(M.xu[-1] - B.xu[-1]) * c_mm             # - = forward
    R["te_move_y_mm"] = float(M.yu[-1] - B.yu[-1]) * c_mm             # - = down
    R["nose_drop_mm"] = (_nose_mid(B, 0.01) - _nose_mid(M, 0.01)) * c_mm   # + = down
    R["dcl_estimate"] = lift_change(B, M)
    R["off_v1_mm"] = 0.0
    if v2:
        V1 = build(deflect_camber(af.cam, mp), af.thk, x)
        moved = ~((x >= mp.x_LE) & (x <= mp.x_h))[station]
        R["off_v1_mm"] = float(point_to_polyline(loopM[moved], loop(V1)).max()) * c_mm

    R["all_pass"] = all(bool(R[k]) for k in (
        "wingbox_pass", "arc_pass", "strain_pass", "slope_join_pass", "curv_join_pass",
        "thickness_pass", "intersect_pass", "te_pass", "le_pass"))

    R["dense"] = {"x": x, "base": B, "morph": M, "loop_base": loopB, "loop_morph": loopM,
                  "k_base": kB, "k_morph": kM, "station": station, "N": N}
    return R


def _cumlen(P):
    """Arc-length fraction (0 to 1) along a polyline."""
    s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))])
    return s / s[-1]


def _sign_changes(k):
    """Inflection points: sign changes of curvature, ignoring near-flat spots
    (|k| < 1e-3 per chord is a radius of 1000 chords)."""
    s = np.sign(k[np.isfinite(k) & (np.abs(k) > 1e-3)])
    return int(np.count_nonzero(np.diff(s) != 0))


def _nose_mid(S, x0):
    """Height of the point halfway between the upper and lower surfaces where
    they cross the vertical line x = x0 (first crossing from the nose)."""
    return (_cross_y(S.xu, S.yu, x0) + _cross_y(S.xl, S.yl, x0)) / 2


def _cross_y(X, Y, x0):
    k = np.nonzero(((X[:-1] - x0) * (X[1:] - x0) <= 0) & (X[1:] != X[:-1]))[0][0]
    w = (x0 - X[k]) / (X[k + 1] - X[k])
    return float(Y[k] + w * (Y[k + 1] - Y[k]))
