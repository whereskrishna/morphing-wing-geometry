"""Camber + thickness -> upper/lower surfaces, point loops, and the .dat reader/splitter.

The shared data types live here too:
    Camber      a camber line: yc(x) and its slope dyc(x)
    Thickness   a thickness form: half-thickness yt(x)
    Airfoil     a named camber + thickness pair
    Surfaces    surfaces built at a set of chord stations
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

import numpy as np
from scipy.interpolate import CubicSpline

Fn = Callable[[np.ndarray], np.ndarray]


@dataclass
class Camber:
    """A camber line. yc(x) and dyc(x) take chord stations x in [0, 1].

    le_blend_x > 0 tells build() how to handle a camber line whose slope runs
    to infinity at the LE (the a = 1.0 mean line); see build().
    """
    name: str
    yc: Fn
    dyc: Fn
    le_blend_x: float = 0.0
    base: Optional["Camber"] = None
    morph: object = None


@dataclass
class Thickness:
    """A thickness form: yt(x) is the half-thickness.

    construction = 'perpendicular' (NACA standard) or 'vertical' (airfoils
    loaded from a .dat file, see split_dat).
    """
    name: str
    yt: Fn
    table: Optional[np.ndarray] = None
    le_radius: float = math.nan
    le_radius_implied: float = math.nan
    te_half_thickness: float = 0.0
    construction: str = "perpendicular"
    source: str = ""


@dataclass
class Airfoil:
    name: str
    cam: Camber
    thk: Thickness
    description: str = ""
    source: str = ""
    split_error: float = math.nan     # loaded .dat only: largest distance file -> rebuilt shape
    split_rms: float = math.nan       # ... and the typical (rms) distance
    frame: str = ""                   # ... whether the file's own frame was kept
    decimals: int = -1                # ... how many decimals the file was published to


@dataclass
class Surfaces:
    """Surfaces at chord stations x, all arrays ordered LE -> TE."""
    x: np.ndarray
    yc: np.ndarray
    dyc: np.ndarray
    yt: np.ndarray
    theta: np.ndarray
    xu: np.ndarray
    yu: np.ndarray
    xl: np.ndarray
    yl: np.ndarray

    def copy(self) -> "Surfaces":
        return Surfaces(**{k: np.array(v, copy=True) for k, v in self.__dict__.items()})


def stations(n: int, spacing: str = "cosine", snap=()) -> np.ndarray:
    """Chord stations from 0 to 1, n points including both ends.

    spacing  'cosine'  points bunch up at the LE and TE (default for output)
             'uniform' equal steps
    snap     stations that must appear exactly (e.g. the spars). The nearest
             interior station is moved onto each one, and its neighbours shift
             with it by less and less (over 10 stations each side), so the
             spacing stays smooth and the point count stays the same.
    """
    assert n >= 3, "Need at least 3 stations."
    i = np.arange(n)
    if spacing == "cosine":
        x = (1 - np.cos(np.pi * i / (n - 1))) / 2
    elif spacing == "uniform":
        x = i / (n - 1)
    else:
        raise ValueError(f'Unknown spacing "{spacing}". Use cosine or uniform.')
    x[0] = 0.0
    x[-1] = 1.0

    m = 10
    for xs in np.atleast_1d(snap):
        k = int(np.argmin(np.abs(x[1:-1] - xs))) + 1
        j = np.arange(max(k - m + 1, 1), min(k + m - 1, n - 2) + 1)   # never move the end points
        shift = xs - x[k]
        x[j] = x[j] + shift * (1 - np.abs(j - k) / m)
        x[k] = xs
    assert np.all(np.diff(x) > 0), "Snapped stations are too close together. Use more points."
    return x


def thickness_fit(x, yt, le_radius: float = math.nan) -> Thickness:
    """Smooth half-thickness function yt(x) through tabulated points.

    x, yt      table as fractions of chord, x from 0 to 1, yt(0) = 0
    le_radius  leading-edge radius (fraction of chord), or nan if unknown

    Method: yt = sqrt(x) * g(x), with g a cubic spline in x.
    - A round nose behaves like yt = sqrt(2 r x), so g is smooth and nearly
      flat at the LE. Splining g (not yt) avoids wiggles near the nose and
      gives yt(0) = 0 exactly.
    - Slope and curvature are continuous everywhere (no jumps at the table
      points), which the curvature and strain checks rely on.
    - g(0) = sqrt(2 r_LE). A given le_radius is used as an extra data point.
      Without one, the spline extrapolates g(0) from the table. Either way
      the radius the fit ends up with is in le_radius_implied.
      (Check: on the PDAS 64(1)-012 table, extrapolating gives 1.002% chord
      against the published 1.007%.)
    """
    x = np.asarray(x, float).ravel()
    yt = np.asarray(yt, float).ravel()
    assert x[0] == 0 and x[-1] == 1 and yt[0] == 0, \
        "Thickness table must run from x = 0 (with yt = 0) to x = 1."
    assert np.all(np.diff(x) > 0), "Thickness table x values must increase."

    xg = x[1:]
    g = yt[1:] / np.sqrt(xg)
    if not math.isnan(le_radius):
        xg = np.concatenate([[0.0], xg])
        g = np.concatenate([[math.sqrt(2 * le_radius)], g])
    spl = CubicSpline(xg, g, bc_type="not-a-knot")
    yt_te = yt[-1]

    def evaluate(xq):
        xq = np.asarray(xq, float)
        y = np.sqrt(xq) * spl(xq)
        return np.where(xq == 1, yt_te, y)   # exact TE value

    return Thickness(name="tabulated thickness", yt=evaluate, table=np.column_stack([x, yt]),
                     le_radius=le_radius, le_radius_implied=float(spl(0.0)) ** 2 / 2,
                     te_half_thickness=float(yt_te))


def build(cam: Camber, thk: Thickness, x) -> Surfaces:
    """Upper and lower surfaces from a camber line and a thickness form.

    Standard NACA construction: the half-thickness is laid off perpendicular
    to the camber line.
        theta = atan(dyc/dx)
        xu = x - yt sin(theta),  yu = yc + yt cos(theta)
        xl = x + yt sin(theta),  yl = yc - yt cos(theta)
    Airfoils loaded from a .dat file use thk.construction = 'vertical'
    (see split_dat). Then theta is taken as 0: yu = yc + yt, yl = yc - yt.

    If the camber line sets le_blend_x (the a = 1.0 mean line does, at
    0.005), theta ahead of that station is a parabola matching theta, its
    slope and its curvature there, instead of following the mean-line slope
    to infinity. That keeps the nose smooth (no kink, no curvature spikes).
    Stations at or behind le_blend_x are not affected.
    """
    x = np.asarray(x, float).ravel()
    yc = cam.yc(x)
    dyc = cam.dyc(x)
    yt = thk.yt(x)
    if thk.construction == "vertical":
        theta = np.zeros_like(x)
    else:
        theta = np.arctan(dyc)
        if cam.le_blend_x > 0:
            # ahead of xb, continue theta as a parabola that matches its value,
            # slope and curvature at xb: smooth, and it keeps falling steadily
            xb = cam.le_blend_x
            h = 1e-3 * xb
            t3 = np.arctan(cam.dyc(np.array([xb - h, xb, xb + h])))
            d1 = (t3[2] - t3[0]) / (2 * h)
            d2 = (t3[2] - 2 * t3[1] + t3[0]) / h ** 2
            inside = x < xb
            theta[inside] = t3[1] + d1 * (x[inside] - xb) + d2 / 2 * (x[inside] - xb) ** 2
    s = np.sin(theta)
    c = np.cos(theta)
    return Surfaces(x=x, yc=yc, dyc=dyc, yt=yt, theta=theta,
                    xu=x - yt * s, yu=yc + yt * c, xl=x + yt * s, yl=yc - yt * c)


def loop(S: Surfaces, order: str = "selig") -> np.ndarray:
    """Join the two surfaces into one ordered point loop, an (m, 2) array.

    order  'selig'   TE -> upper -> LE -> lower -> TE
           'reverse' TE -> lower -> LE -> upper -> TE
    The LE point (station x = 0) appears once. For a closed trailing edge the
    first and last points are the same point.
    """
    assert S.x[0] == 0, "The first station must be the LE (x = 0)."
    up = np.column_stack([S.xu, S.yu])
    lo = np.column_stack([S.xl, S.yl])
    if order == "selig":
        return np.vstack([up[::-1], lo[1:]])
    if order == "reverse":
        return np.vstack([lo[::-1], up[1:]])
    raise ValueError(f'Unknown point order "{order}". Use selig or reverse.')


@dataclass
class DatFile:
    name: str
    xy: np.ndarray
    format: str


def read_dat(path) -> DatFile:
    """Read an airfoil coordinate file in Selig or Lednicer format.

    Returns name (first line), xy (one loop in Selig order:
    TE -> upper -> LE -> lower -> TE) and format ('selig' or 'lednicer').

    Selig:    name line, then x y pairs around the loop.
    Lednicer: name line, a line with the two point counts (e.g. "26. 26."),
              then the upper surface LE -> TE, then the lower surface LE -> TE.
              (UIUC stores most NACA sections this way.)
    Windows line endings are fine.
    """
    lines = Path(path).read_text().replace("\r", "").split("\n")
    name = lines[0].strip()
    rows = []
    for line in lines[1:]:
        parts = line.split()
        if len(parts) == 2:
            try:
                rows.append([float(parts[0]), float(parts[1])])
            except ValueError:
                pass
    rows = np.array(rows, float).reshape(-1, 2)
    assert rows.shape[0] >= 5, f"Found fewer than 5 points in {path}."

    if np.all(rows[0] > 1):
        # Lednicer: first numeric line holds the point counts
        nu, nl = int(round(rows[0, 0])), int(round(rows[0, 1]))
        pts = rows[1:]
        assert pts.shape[0] == nu + nl, \
            f"Point counts in {path} ({nu} + {nl}) do not match the data ({pts.shape[0]})."
        up, lo = pts[:nu], pts[nu:]
        if np.array_equal(up[0], lo[0]):
            lo = lo[1:]          # shared LE point appears once
        return DatFile(name, np.vstack([up[::-1], lo]), "lednicer")
    return DatFile(name, rows, "selig")


def split_dat(path, name: str, smooth: bool = True, n_stations: int = 401) -> Airfoil:
    """Load any Selig/Lednicer .dat airfoil and split it into camber + thickness.

    Returns the same Airfoil as naca.naca6_section, so a loaded airfoil
    (SC(2)-0714, RAE 2822, ...) can be morphed the same way.

    THIS SPLIT IS APPROXIMATE. A coordinate file has no separate camber and
    thickness tables, so they are estimated with a VERTICAL split:
      1. Frame: if the file already has its LE at (0, 0) and its TE at x = 1
         (NASA and UIUC files do), its own frame is kept, so angles of attack
         stay comparable with published test data. Otherwise: LE = the point
         farthest from the TE midpoint; shift, rotate and scale so the LE is
         at (0, 0) and the TE midpoint at (1, 0).
      2. At each station x: camber = midpoint of the two surfaces,
         half-thickness = half the vertical gap.
      3. Smooth curves through those (see below).
    Because of this, loaded airfoils lay their thickness VERTICALLY when the
    surfaces are rebuilt (construction = 'vertical'), not perpendicular to
    the camber line like the NACA sections. A perpendicular split was tried
    and rejected: near the nose a tiny camber change rotates the cut lines
    enough that the fit runs away.

    Smoothing (smooth=True, default). Published coordinates are rounded
    (SC(2) files: 4 decimals). A curve passed exactly through rounded points
    has curvature that swings up and down at every point, which would swamp
    the strain and inflection checks and show as small pressure wiggles in
    CFD. So camber and thickness are fitted as smoothly as possible while
    the misfit stays at the size of that rounding (the number of decimals is
    read from the file). The nose point and the TE (height and gap) are kept
    exactly. smooth=False passes exactly through every point instead.

    split_error is the largest distance (fraction of chord) between the
    file's points and the rebuilt baseline; split_rms the typical one.
    Check them before using the shape.
    """
    d = read_dat(path)
    decimals = _decimals(path)

    # 1. frame
    xy = d.xy
    te = (xy[0] + xy[-1]) / 2
    ile = int(np.argmax(np.linalg.norm(xy - te, axis=1)))
    le = xy[ile]
    if np.allclose(le, 0, atol=1e-9) and abs(xy[0, 0] - 1) < 1e-9 and abs(xy[-1, 0] - 1) < 1e-9:
        P, frame = xy.copy(), "file frame kept"
    else:
        chord = float(np.linalg.norm(te - le))
        ang = math.atan2(te[1] - le[1], te[0] - le[0])
        rot = np.array([[math.cos(ang), -math.sin(ang)], [math.sin(ang), math.cos(ang)]])
        P = (xy - le) @ rot / chord          # rotate by -ang, scale to unit chord
        frame = f"normalised: turned {math.degrees(-ang):+.3f} deg, scaled x{1 / chord:.6f}"
    up = P[: ile + 1][::-1]                  # both surfaces LE -> TE
    lo = P[ile:]
    if up[:, 1].mean() < lo[:, 1].mean():
        up, lo = lo, up                      # file was in reverse order

    if smooth:
        cam, thk = _smooth_split(up, lo, decimals)
    else:
        cam, thk = _exact_split(up, lo, n_stations)
    cam.name = f"camber split from {d.name}"
    thk.name = f"thickness split from {d.name}"
    thk.construction = "vertical"
    thk.source = str(path)

    af = Airfoil(name=name, cam=cam, thk=thk, source=str(path),
                 description=f"vertical split of {path} ({'smoothed' if smooth else 'exact'}, {frame})")
    # how well does camber + thickness rebuild the file?
    S = build(cam, thk, stations(4001, "cosine"))
    dist = point_to_polyline(P, loop(S))
    af.split_error = float(dist.max())
    af.split_rms = float(np.sqrt(np.mean(dist ** 2)))
    af.frame = frame
    af.decimals = decimals
    return af


def _decimals(path) -> int:
    """How many decimals the coordinates were really given to (trailing zeros
    don't count: 0.048400 is 4 decimals)."""
    vals = []
    for line in Path(path).read_text().replace("\r", "").split("\n")[1:]:
        parts = line.split()
        if len(parts) == 2:
            try:
                vals += [float(parts[0]), float(parts[1])]
            except ValueError:
                pass
    vals = np.array(vals)
    vals = vals[np.abs(vals) <= 1]           # skip a Lednicer count line
    for dec in range(0, 9):
        if np.all(np.abs(vals * 10 ** dec - np.round(vals * 10 ** dec)) < 1e-6):
            return dec
    return 9


def _station_split(up, lo):
    """Upper and lower heights at every station either surface was given at.
    Where a surface has a published point at a station, that value is used;
    elsewhere it is read off a smooth curve through that surface."""
    x = np.unique(np.concatenate([up[:, 0], lo[:, 0]]))
    x = x[(x >= 0) & (x <= min(up[-1, 0], lo[-1, 0]))]
    return x, _heights(up, x), _heights(lo, x)


def _heights(P, x):
    y = _vertical_cut(_densify(P), x)
    own = np.isin(x, P[:, 0])
    first = {xi: yi for xi, yi in zip(P[::-1, 0], P[::-1, 1])}   # first published point at each x
    y[own] = [first[xi] for xi in x[own]]
    return y


def _smooth_split(up, lo, decimals):
    """Camber and thickness fitted as smoothly as the file's rounding allows."""
    from scipy.interpolate import UnivariateSpline

    x, yu, yl = _station_split(up, lo)
    yt = (yu - yl) / 2
    yc = (yu + yl) / 2
    # How far the points scatter about a smooth curve. Measured from the data
    # (published decimals can overstate the real accuracy: SC(2)-0714 is
    # printed to 5 decimals but scatters like 4). Never less than the
    # rounding itself: +-half a unit spread evenly, sd 0.5e-d / sqrt(3),
    # and camber / half-thickness average two values.
    floor = 0.5 * 10.0 ** (-decimals) / math.sqrt(3) / math.sqrt(2)
    sig_c = max(_scatter(x, yc), floor)
    sig_t = max(_scatter(x, yt), floor)
    x0 = x[1:]
    g = yt[1:] / np.sqrt(x0)                  # as in thickness_fit: yt = sqrt(x) g(x)
    gs = UnivariateSpline(x0, g, w=np.sqrt(x0) / sig_t, k=3, s=x0.size)
    cs = UnivariateSpline(x, yc, w=np.full_like(x, 1 / sig_c), k=3, s=x.size)
    dcs = cs.derivative()
    # pin the nose (camber 0) and the TE (height and gap) to the file exactly;
    # the straight-line corrections are at the size of the rounding
    yc_le, yc_te = float(cs(0.0)), float(cs(1.0))
    t_te = float(yt[-1])
    g_te = float(gs(1.0))

    def cam_y(q):
        q = np.asarray(q, float)
        return cs(q) - (1 - q) * yc_le + q * (yc[-1] - yc_te)

    def cam_dy(q):
        q = np.asarray(q, float)
        return dcs(q) + yc_le + (yc[-1] - yc_te)

    def thk_y(q):
        q = np.asarray(q, float)
        return np.sqrt(q) * (gs(q) + q * (t_te - g_te))

    cam = Camber(name="", yc=cam_y, dyc=cam_dy)
    thk = Thickness(name="", yt=thk_y, table=np.column_stack([x, yt]),
                    le_radius_implied=float(gs(0.0)) ** 2 / 2, te_half_thickness=t_te)
    return cam, thk


def _scatter(x, y, m: int = 7) -> float:
    """Typical scatter of the points about a smooth curve: fit a cubic through
    each run of m neighbouring points (x from 5% to 95% chord, away from the
    nose and TE) and see how far the middle point is off it. Only runs no
    wider than 10% of chord count: over a wider run a cubic can't follow the
    real shape, and that misfit would be taken for scatter. 0 if the file
    has too few close points to tell."""
    sel = (x >= 0.05) & (x <= 0.95)
    x, y = x[sel], y[sel]
    r2 = []
    c = m // 2
    for i in range(x.size - m + 1):
        xs, ys = x[i:i + m], y[i:i + m]
        if xs.max() - xs.min() > 0.1:
            continue
        A = np.vander((xs - xs.mean()) / (xs.max() - xs.min()), 4)
        Hm = A @ np.linalg.pinv(A)            # the local cubic fit, as a matrix
        r2.append((ys[c] - Hm[c] @ ys) ** 2 / (1 - Hm[c, c]))
    if len(r2) < 5:
        return 0.0
    return float(np.sqrt(np.mean(r2)))


def _exact_split(up, lo, n_stations):
    """Camber and thickness passing exactly through dense vertical cuts."""
    x = stations(n_stations, "cosine")
    yu = _vertical_cut(_densify(up), x)
    yl = _vertical_cut(_densify(lo), x)
    yc = (yu + yl) / 2
    yt = np.maximum((yu - yl) / 2, 0)
    yc[0] = 0.0
    yt[0] = 0.0                              # LE is the origin
    spl = CubicSpline(x, yc, bc_type="not-a-knot")
    dspl = spl.derivative()
    cam = Camber(name="", yc=lambda xq: spl(np.asarray(xq, float)), dyc=lambda xq: dspl(np.asarray(xq, float)))
    return cam, thickness_fit(x, yt)


def _densify(P: np.ndarray) -> np.ndarray:
    """Smooth curve through the points, parameterised by distance along them."""
    s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))])
    keep = np.concatenate([[True], np.diff(s) > 0])
    P, s = P[keep], s[keep]
    q = s[-1] * (1 - np.cos(np.pi * np.arange(4001) / 4000)) / 2
    return np.column_stack([CubicSpline(s, P[:, 0])(q), CubicSpline(s, P[:, 1])(q)])


def _vertical_cut(C: np.ndarray, x: np.ndarray) -> np.ndarray:
    """y where curve C crosses each vertical line x = const (first crossing from
    the curve's start; nearest point if it never crosses)."""
    y = np.zeros_like(x)
    for i, xi in enumerate(x):
        f = C[:, 0] - xi
        k = np.nonzero(f[:-1] * f[1:] <= 0)[0]
        if k.size == 0:
            y[i] = C[np.argmin(np.abs(f)), 1]
            continue
        k = k[0]
        if f[k] == f[k + 1]:
            y[i] = C[k, 1]
        else:
            w = f[k] / (f[k] - f[k + 1])
            y[i] = C[k, 1] + w * (C[k + 1, 1] - C[k, 1])
    return y


def point_to_polyline(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Distance from each point of P to the polyline through the points Q."""
    A = Q[:-1]
    AB = Q[1:] - A
    len2 = np.maximum(np.sum(AB ** 2, axis=1), np.finfo(float).eps)
    dist = np.zeros(P.shape[0])
    for i, p in enumerate(P):
        t = np.clip(((p[0] - A[:, 0]) * AB[:, 0] + (p[1] - A[:, 1]) * AB[:, 1]) / len2, 0, 1)
        dist[i] = np.min(np.hypot(A[:, 0] + t * AB[:, 0] - p[0], A[:, 1] + t * AB[:, 1] - p[1]))
    return dist
