"""The morphs: v1 'camber' (shift the camber line, keep the thickness) and
v2 'constant_length' (skins bend but keep their length).

Morph settings (cfg.morph) are read as attributes:
    x_LE, x_h, x_rear_spar, delta_LE, delta_TE, n, method, n_modes, le_smoothing, le_full_droop

x_h is the hinge: where the TE starts to bend. x_rear_spar is the back of
the wing box. The hinge can be anywhere from the rear spar back; the skin
between them does not move. x_rear_spar = None puts the spar at the hinge.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np
from scipy.interpolate import CubicSpline

from .surfaces import Airfoil, Camber, Surfaces, build, stations


# --------------------------------------------------------------------- v1

def camber_delta(x, mp):
    """Camber change from the LE and TE morph, and its slope: (dy, ddy).

    Trailing edge, x > x_h (hinge, at or behind the rear spar):
        s  = (x - x_h) / (c - x_h)
        dy = -delta_TE * c * s^n
    Leading edge, x < x_LE (front spar):
        s  = (x_LE - x) / x_LE
        dy = -delta_LE * c * s^n
    Between the front spar and the hinge dy and its slope are exactly 0
    (fixed wing box, plus any fixed skin between the rear spar and the hinge).

    delta > 0 = droop (edge down, more camber); delta < 0 = reflex (edge up).
    n >= 2: slope matches at the spar. n >= 3: curvature matches too.
    """
    x = np.asarray(x, float)
    dy = np.zeros_like(x)
    ddy = np.zeros_like(x)

    te = x > mp.x_h
    s = (x[te] - mp.x_h) / (1 - mp.x_h)
    dy[te] = -mp.delta_TE * s ** mp.n
    ddy[te] = -mp.delta_TE * mp.n * s ** (mp.n - 1) / (1 - mp.x_h)

    le = x < mp.x_LE
    s = (mp.x_LE - x[le]) / mp.x_LE
    dy[le] = -mp.delta_LE * s ** mp.n
    ddy[le] = mp.delta_LE * mp.n * s ** (mp.n - 1) / mp.x_LE
    return dy, ddy


def rear_spar(mp) -> float:
    """Rear spar station: x_rear_spar, or the hinge x_h if that is not set."""
    x = getattr(mp, "x_rear_spar", None)
    return mp.x_h if x is None else x


def _validate(mp):
    for f in ("x_LE", "x_h", "delta_LE", "delta_TE", "n"):
        v = getattr(mp, f, None)
        if v is None or not np.isscalar(v) or not math.isfinite(v):
            raise ValueError(f'Morph setting "{f}" is missing or not a finite number.')
    xrs = rear_spar(mp)
    if not np.isscalar(xrs) or not math.isfinite(xrs):
        raise ValueError('Morph setting "x_rear_spar" must be a number or null.')
    if not (0 < mp.x_LE < xrs <= mp.x_h < 1):
        raise ValueError(f"Need 0 < x_LE < x_rear_spar <= x_h < 1: the hinge sits at or behind the rear spar "
                         f"(got x_LE = {mp.x_LE:.3f}, x_rear_spar = {xrs:.3f}, x_h = {mp.x_h:.3f}).")
    if mp.n < 1:
        raise ValueError(f"Need n >= 1 (got {mp.n:.3f}).")


def deflect_camber(base: Camber, mp) -> Camber:
    """Morphed camber line = baseline camber + LE/TE deflection.

    The thickness is not touched. Build the morphed surfaces with
    surfaces.build(cam, af.thk, x), which lays the unchanged thickness
    perpendicular to the new camber line. Shapes stay in the baseline frame
    (no re-normalising), so the wing box does not move.
    """
    _validate(mp)

    def yc(x):
        return base.yc(x) + camber_delta(x, mp)[0]

    def dyc(x):
        return base.dyc(x) + camber_delta(x, mp)[1]

    return Camber(name=f"{base.name}, morphed", yc=yc, dyc=dyc,
                  le_blend_x=base.le_blend_x, base=base, morph=mp)


def equivalent_angles(mp) -> dict:
    """Morph size as a plain flap-style deflection angle, degrees (+ = droop).
        te_deg = atan(delta_TE * c / (c - x_h)),  le_deg = atan(delta_LE * c / x_LE)
    i.e. the angle of a straight flap hinged at x_h with the same tip drop.
    """
    return {"te_deg": math.degrees(math.atan(mp.delta_TE / (1 - mp.x_h))),
            "le_deg": math.degrees(math.atan(mp.delta_LE / mp.x_LE))}


# --------------------------------------------------------------------- v2

@dataclass
class ConstantLengthSolution:
    x: np.ndarray
    eval: Callable[[np.ndarray], Surfaces]
    te: Optional[dict] = None
    le: Optional[dict] = None


def constant_length(af: Airfoil, mp, cfg) -> ConstantLengthSolution:
    """v2 morph: skins bend but do not stretch.

    sol.eval(x) returns the morphed surfaces at any chord stations x, in the
    same form as surfaces.build. Each station is a skin point: sol.eval says
    where the baseline skin point at station x ends up.

    Each skin is cut into many short straight segments. The morph only TURNS
    segments, never stretches them, so the skin keeps its length exactly.
    The turn angle dphi grows smoothly from the spar; sigma = distance from
    the spar / skin length.

    Trailing edge (sliding-TE layout)
      A closed, sharp TE whose skins cannot stretch, held apart by ribs that
      keep the thickness, cannot droop: the outer skin needs more length on
      a bend than the inner one. Real morphing TEs let something give.
      Here the UPPER skin runs unbroken to the tip and keeps its length; the
      LOWER skin slides a little into the TE (the checks report how far).
      - Upper skin: dphi = a * sigma^(n-1), the same curvature law as v1
        (n = 2: slope matches at the spar, n = 3: curvature too). a is set
        so the TE drops by delta_TE. The TE moves forward a little.
      - Lower skin: each rib (baseline line from an upper skin point to the
        lower skin point at the same station) keeps its length and its angle
        to the upper skin. So the thickness is unchanged, the skins cannot
        cross, and the TE stays closed.

    Leading edge (constant-arc-length nose, as in the Swansea morphing-LE work)
      One skin wraps round the nose, clamped at both front-spar points
      (position and slope; curvature too for n >= 3), and keeps its length.
      dphi = sigma^(n-1) (1-sigma)^(n-1) * (c_0 + c_1 sigma + ...), n_modes
      terms (few terms keep the skin's curvature change smooth, and so the
      strain low).
      - Full droop (le_full_droop = True, default): the nose must drop
        exactly as far as the v1 nose, measured as the height of the point
        midway between the upper and lower surfaces at 1% chord. v1 itself
        is NOT a constant-length shape (it needs about 1% more skin for a 1%
        droop), so this gives a new droop shape with the same skin length.
        The skin bends mostly just ahead of the front spar and the nose
        swings down nearly as a whole.
      - Within that, the shape stays as close to the v1 droop as it can (sum
        of squared distances from the v1 surface) while bending smoothly
        (curvature change squared); le_smoothing sets the balance.
      - le_full_droop = False drops the droop condition: closest shape to
        v1, which reaches only about half the droop.
      Why not simply reshape v1 to the right length (e.g. by pulling the
      nose back a little)? That keeps the total length but makes the skin
      slide ~5 mm round the nose, so tip material ends up where the skin is
      much flatter: about 7% strain for a 1% droop. Bending the skin through
      smooth modes avoids that (about 0.3%).
    """
    n_modes = getattr(mp, "n_modes", 4)
    smoothing = getattr(mp, "le_smoothing", 0.01)
    full_droop = getattr(mp, "le_full_droop", True)
    if mp.n < 2:
        raise ValueError("The constant-length morph needs n >= 2 (n = 1 puts a kink in the skin at the spar).")

    # Work on a grid twice as fine as the checks use: the segment chain is
    # second-order accurate, so this keeps its error well below what the
    # checks can see.
    x = stations(2 * cfg.checks.n_points - 1, "cosine", [mp.x_LE, mp.x_h])
    N = x.size
    ile = int(np.nonzero(x == mp.x_LE)[0][0])
    ih = int(np.nonzero(x == mp.x_h)[0][0])
    B = build(af.cam, af.thk, x)
    U = np.column_stack([B.xu, B.yu])        # skin points by station, LE -> TE
    L = np.column_stack([B.xl, B.yl])
    DU = np.zeros((N, 2))                    # how far each skin point moves
    DL = np.zeros((N, 2))

    sol_te = sol_le = None
    if mp.delta_TE != 0:
        Qu, Ql, sol_te = _sliding_te(U[ih:], L[ih:], mp.n, mp.delta_TE)
        DU[ih:] = Qu - U[ih:]
        DL[ih:] = Ql - L[ih:]
    if mp.delta_LE != 0:
        T = build(deflect_camber(af.cam, mp), af.thk, x)              # v1 shape
        wrap = np.vstack([U[ile::-1], L[1:ile + 1]])                   # upper spar -> nose -> lower spar
        twrap = np.vstack([np.column_stack([T.xu, T.yu])[ile::-1],
                           np.column_stack([T.xl, T.yl])[1:ile + 1]])
        mid = None
        if full_droop:
            # reach the same nose height as v1 at 1% chord
            Z = np.zeros((twrap.shape[0], 1))
            mid = {"x": 0.01}
            mid["y"] = (_crossing(twrap, Z, Z, mid["x"], "first")[0]
                        + _crossing(twrap, Z, Z, mid["x"], "last")[0]) / 2
        Qw, sol_le = _nose_wrap(wrap, twrap, ile, mp.n, n_modes, smoothing, mp.delta_LE, mid)
        DU[:ile + 1] = Qw[ile::-1] - U[:ile + 1]
        DL[:ile + 1] = Qw[ile:] - L[:ile + 1]

    # Smooth fits of the movement in each region, in u = sqrt(x) (smooth at the
    # nose), so the shape can be evaluated at any stations. Between the spars
    # the movement is exactly zero.
    fits = {}
    if mp.delta_TE != 0:
        fits["teu"] = CubicSpline(np.sqrt(x[ih:]), DU[ih:], axis=0)
        fits["tel"] = CubicSpline(np.sqrt(x[ih:]), DL[ih:], axis=0)
    if mp.delta_LE != 0:
        fits["leu"] = CubicSpline(np.sqrt(x[:ile + 1]), DU[:ile + 1], axis=0)
        fits["lel"] = CubicSpline(np.sqrt(x[:ile + 1]), DL[:ile + 1], axis=0)
    closed_te = bool(np.array_equal(U[-1], L[-1]))
    return ConstantLengthSolution(x=x, eval=lambda xq: _evaluate(xq, af, mp, fits, closed_te),
                                  te=sol_te, le=sol_le)


def _sliding_te(Pu, Pl, n, delta):
    """TE: upper skin keeps its length, ribs carry the lower skin."""
    V, s_mid, Lu = _segments(Pu)
    s_pt = np.concatenate([[0.0], np.cumsum(np.linalg.norm(V, axis=1))])
    b_seg = (s_mid / Lu) ** (n - 1)          # turn-angle shape at segment middles
    b_pt = (s_pt / Lu) ** (n - 1)            # ... and at the skin points (for the ribs)
    y_target = Pu[-1, 1] - delta

    a = 0.0                                   # Newton on the one unknown
    for it in range(1, 51):
        Qu, R = _bend(Pu[0], V, a * b_seg)
        f = Qu[-1, 1] - y_target
        if abs(f) < 1e-14:
            break
        a -= f / np.sum(b_seg * R[:, 0])      # d(end y)/da
    if abs(f) >= 1e-12:
        raise RuntimeError(f"Constant-length TE morph did not converge (delta_TE = {delta:.3f} may be too large).")
    rib = Pl - Pu                             # baseline ribs
    ang = a * b_pt
    Ql = Qu + np.column_stack([rib[:, 0] * np.cos(ang) - rib[:, 1] * np.sin(ang),
                               rib[:, 0] * np.sin(ang) + rib[:, 1] * np.cos(ang)])
    info = {"a": a, "iterations": it, "L_upper": Lu, "tip_move": Qu[-1] - Pu[-1],
            "lower_slide": np.sum(np.linalg.norm(np.diff(Pl, axis=0), axis=1))
                           - np.sum(np.linalg.norm(np.diff(Ql, axis=0), axis=1))}
    return Qu, Ql, info


def _nose_wrap(P, T, itip, n, K, smoothing, delta, mid):
    """LE: one skin round the nose, clamped at both front-spar points.

    P: baseline wrap points (upper spar -> nose -> lower spar), T: the v1 wrap.
    mid: None or {'x', 'y'}: the nose must reach mid-height y at x (droop).
    """
    V, s_mid, Lw = _segments(P)
    Bm = _basis(s_mid / Lw, n, K)
    # Powers of sigma look alike, which makes the fit swing wildly. Swap them
    # for independent (orthonormal) combinations: same shapes, stable fit.
    l = np.linalg.norm(V, axis=1)
    _, Rb = np.linalg.qr(np.sqrt(l / Lw)[:, None] * Bm, mode="reduced")
    Rinv = np.linalg.inv(Rb)
    Bm = Bm @ Rinv
    H = Rinv.T @ _energy(n, K) @ Rinv / Lw    # bending energy = c' H c
    w = (np.concatenate([l, [0.0]]) + np.concatenate([[0.0], l])) / 2 / Lw   # every bit of skin counts the same
    ncon = 2 + (mid is not None)

    def conditions(Q, R):
        return _nose_conditions(Q, R, Bm, P[-1], mid)

    # Gauss-Newton steps on the distance from the v1 shape, with the
    # conditions held to first order. A step is only taken if it lowers
    # "misfit + rho * condition error", otherwise it is halved; no segment
    # turns more than 0.1 rad further in one step. The bending term (weight =
    # smoothing, relative to the misfit) keeps the skin smooth.
    c = np.zeros(K)
    match = np.arange(P.shape[0])
    rho = 0.0
    Q, R = _bend(P[0], V, Bm @ c)
    e, Je, match, _ = _misfit(Q, R, Bm, T, w, match)
    g, Jg = conditions(Q, R)
    prev = math.inf
    for it in range(1, 151):
        A = Je.T @ Je
        mh = max(smoothing, 1e-6) * np.trace(A) / np.trace(H)
        kkt = np.block([[A + mh * H, Jg.T], [Jg, np.zeros((ncon, ncon))]])
        sol = np.linalg.solve(kkt, np.concatenate([-(Je.T @ e) - mh * (H @ c), -g]))
        step = sol[:K]
        turn = np.max(np.abs(Bm @ step))
        step = step * min(1.0, 0.1 / turn)
        rho = max(rho, 2 * np.max(np.abs(sol[K:])))

        def score(e_, c_, g_):
            return 0.5 * (e_ @ e_) + 0.5 * mh * (c_ @ H @ c_) + rho * np.sum(np.abs(g_))

        s0 = score(e, c, g)
        a = 1.0
        for _ls in range(30):
            ct = c + a * step
            Qt, Rt = _bend(P[0], V, Bm @ ct)
            et, Jet, mt, _ = _misfit(Qt, Rt, Bm, T, w, match)
            gt, Jgt = conditions(Qt, Rt)
            if score(et, ct, gt) <= s0:
                break
            a /= 2
        c, Q, R, e, Je, match, g, Jg = ct, Qt, Rt, et, Jet, mt, gt, Jgt
        s1 = score(e, c, g)
        if a * min(turn, 0.1) < 1e-10 or abs(prev - s1) <= 1e-13 * max(s1, np.finfo(float).eps):
            break
        prev = s1
    for _k in range(20):                      # land exactly on the conditions
        if np.max(np.abs(g)) < 1e-14:
            break
        c = c - Jg.T @ np.linalg.solve(Jg @ Jg.T, g)
        Q, R = _bend(P[0], V, Bm @ c)
        g, Jg = conditions(Q, R)
    if np.max(np.abs(g)) >= 1e-12:
        raise RuntimeError(f"Constant-length LE morph did not converge (delta_LE = {delta:.3f} may be too large).")
    Q[-1] = P[-1]                             # exactly back on the lower spar point
    _, _, _, d = _misfit(Q, R, Bm, T, w, match)
    info = {"c": c, "iterations": it, "L_wrap": Lw, "tip_move": Q[itip] - P[itip],
            "max_off_v1": float(np.max(np.abs(d)))}
    return Q, info


def _nose_conditions(Q, R, Bm, spar_end, mid):
    """The wrap ends back on the lower spar point; optionally the nose reaches
    the droop target (mid-height between the surfaces at x = mid['x'])."""
    Jx, Jy = _point_jacobians(R, Bm)
    g = Q[-1] - spar_end
    Jg = np.vstack([Jx[-1], Jy[-1]])
    if mid is not None:
        yu, Ju = _crossing(Q, Jx, Jy, mid["x"], "first")   # upper side
        yl, Jl = _crossing(Q, Jx, Jy, mid["x"], "last")    # lower side
        g = np.concatenate([g, [(yu + yl) / 2 - mid["y"]]])
        Jg = np.vstack([Jg, (Ju + Jl) / 2])
    return g, Jg


def _crossing(Q, Jx, Jy, x0, which):
    """Height where the wrap crosses the vertical line x = x0, and how it moves
    with the mode amplitudes. 'first' = from the upper spar, 'last' = from the lower."""
    f = Q[:, 0] - x0
    idx = np.nonzero((f[:-1] * f[1:] <= 0) & (f[:-1] != f[1:]))[0]
    if idx.size == 0 or not np.min(Q[:, 0]) < x0:
        raise RuntimeError(f"The nose moved behind {100 * x0:.0f}% chord: the LE morph is too large.")
    k = idx[0] if which == "first" else idx[-1]
    w = -f[k] / (f[k + 1] - f[k])
    y = Q[k, 1] + w * (Q[k + 1, 1] - Q[k, 1])
    slope = (Q[k + 1, 1] - Q[k, 1]) / (Q[k + 1, 0] - Q[k, 0])
    J = (1 - w) * Jy[k] + w * Jy[k + 1] - slope * ((1 - w) * Jx[k] + w * Jx[k + 1])
    return y, J


def _point_jacobians(R, Bm):
    """How each skin point moves per unit change of each mode amplitude."""
    K = Bm.shape[1]
    Jx = np.vstack([np.zeros((1, K)), np.cumsum(-R[:, 1:2] * Bm, axis=0)])
    Jy = np.vstack([np.zeros((1, K)), np.cumsum(R[:, 0:1] * Bm, axis=0)])
    return Jx, Jy


def _misfit(Q, R, Bm, T, w, match):
    """Weighted distance of each skin point from the v1 curve T, measured
    square to the curve (sliding along it is free), and its derivative."""
    nrm, match, d = _project(Q, T, match)
    Jx, Jy = _point_jacobians(R, Bm)
    sw = np.sqrt(w)
    e = sw * d
    Je = sw[:, None] * (nrm[:, 0:1] * Jx + nrm[:, 1:2] * Jy)
    return e, Je, match, d


def _project(Q, T, seg0):
    """Closest point on polyline T for each point of Q, searching near the
    segment it matched last time. Returns the unit normal there, the segment,
    and the signed distance along the normal."""
    W = 300
    ns = T.shape[0] - 1
    cand = np.clip(seg0[:, None] + np.arange(-W, W + 1)[None, :], 0, ns - 1)
    A = T[cand]
    Bv = T[cand + 1] - A
    qx = Q[:, 0:1] - A[..., 0]
    qy = Q[:, 1:2] - A[..., 1]
    l2 = np.maximum(Bv[..., 0] ** 2 + Bv[..., 1] ** 2, np.finfo(float).eps)
    tt = np.clip((qx * Bv[..., 0] + qy * Bv[..., 1]) / l2, 0, 1)
    k = np.argmin((qx - tt * Bv[..., 0]) ** 2 + (qy - tt * Bv[..., 1]) ** 2, axis=1)
    seg = cand[np.arange(cand.shape[0]), k]
    bx = T[seg + 1, 0] - T[seg, 0]
    by = T[seg + 1, 1] - T[seg, 1]
    bl = np.maximum(np.hypot(bx, by), np.finfo(float).eps)
    nrm = np.column_stack([-by / bl, bx / bl])
    d = np.sum((Q - T[seg]) * nrm, axis=1)
    return nrm, seg, d


def _segments(P):
    V = np.diff(P, axis=0)
    l = np.linalg.norm(V, axis=1)
    return V, np.cumsum(l) - l / 2, float(np.sum(l))


def _bend(P0, V, dphi):
    """Turn every segment by its own angle and chain them from P0. Lengths stay."""
    R = np.column_stack([V[:, 0] * np.cos(dphi) - V[:, 1] * np.sin(dphi),
                         V[:, 0] * np.sin(dphi) + V[:, 1] * np.cos(dphi)])
    return np.vstack([P0, P0 + np.cumsum(R, axis=0)]), R


def _basis(sigma, n, K):
    """Nose turn-angle modes: sigma^(n-1+j) (1-sigma)^(n-1), zero at both spars."""
    sigma = np.asarray(sigma, float)
    return np.column_stack([sigma ** (n - 1 + j) * (1 - sigma) ** (n - 1) for j in range(K)])


def _energy(n, K):
    """G[j, k] = integral over sigma of (d mode_j/d sigma)(d mode_k/d sigma).
    Divided by the skin length this is integral (curvature change)^2 ds."""
    M = 20000
    s = (np.arange(1, M + 1) - 0.5) / M      # midpoint rule
    b = n - 1
    D = np.column_stack([(n - 1 + j) * s ** (n - 2 + j) * (1 - s) ** b
                         - b * s ** (n - 1 + j) * (1 - s) ** (b - 1) for j in range(K)])
    return D.T @ D / M


def _evaluate(xq, af, mp, fits, closed_te):
    """The morphed surfaces at any stations."""
    S = build(af.cam, af.thk, xq)
    x = S.x
    te = x > mp.x_h
    le = x < mp.x_LE
    if "teu" in fits:
        _move(S, te, fits["teu"](np.sqrt(x[te])), fits["tel"](np.sqrt(x[te])))
    if "leu" in fits:
        _move(S, le, fits["leu"](np.sqrt(x[le])), fits["lel"](np.sqrt(x[le])))
    # shared points stay shared: nose tip (x = 0) and a closed TE (x = 1)
    f = x == 0
    S.xl[f], S.yl[f] = S.xu[f], S.yu[f]
    if closed_te:
        f = x == 1
        S.xl[f], S.yl[f] = S.xu[f], S.yu[f]
    # there is no camber line any more: keep the mid-line and half-gap instead
    S.yc = (S.yu + S.yl) / 2
    S.yt = np.hypot(S.xu - S.xl, S.yu - S.yl) / 2
    S.dyc = np.full_like(x, np.nan)
    S.theta = np.full_like(x, np.nan)
    return S


def _move(S, f, du, dl):
    S.xu[f] += du[:, 0]
    S.yu[f] += du[:, 1]
    S.xl[f] += dl[:, 0]
    S.yl[f] += dl[:, 1]


# --------------------------------------------------------------- choose

@dataclass
class Morph:
    method: str
    morph: object
    eval: Callable[[np.ndarray], Surfaces]
    cam: Optional[Camber] = None
    sol: Optional[ConstantLengthSolution] = None


def apply(af: Airfoil, mp, cfg) -> Morph:
    """Morph the baseline with the chosen method.

    mp.method  'camber'          v1: shift the camber line, keep the thickness
                                 (skins stretch a little; see the checks)
               'constant_length' v2: skins bend but keep their length
    """
    method = getattr(mp, "method", "camber")
    if method == "camber":
        cam = deflect_camber(af.cam, mp)
        return Morph(method, mp, lambda x: build(cam, af.thk, x), cam=cam)
    if method == "constant_length":
        _validate(mp)                          # same setting checks as v1
        sol = constant_length(af, mp, cfg)
        return Morph(method, mp, sol.eval, sol=sol)
    raise ValueError(f'Unknown morph method "{method}". Use camber or constant_length.')
