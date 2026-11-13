"""NACA building blocks: the a = 1.0 mean line, 6-series assembly, thickness
tables, and NACA 4-digit sections (used by the tests)."""
from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np

from .surfaces import Airfoil, Camber, Thickness, thickness_fit


def camber_a10(cl_i: float) -> Camber:
    """NACA a = 1.0 mean line (uniform loading).

        yc/c   = -(cl_i/(4*pi)) * [ (1-x) ln(1-x) + x ln(x) ]
        dyc/dx =  (cl_i/(4*pi)) * ln( (1-x)/x )

    Endpoints: yc is exactly 0 at x = 0 and x = 1 (x ln x -> 0). The slope
    really is infinite there (+inf at the LE, -inf at the TE). arctan() turns
    that into +/-90 deg, and the half-thickness is 0 at both ends, so the
    surface points at x = 0 and x = 1 still come out exact.

    Nose: laying thickness along that run-away slope puts a microscopic kink
    in the nose (below 0.001% chord), with fake curvature spikes. NACA's
    tables have no points between x = 0 and 0.005 (only the LE radius), so
    le_blend_x = 0.005 tells surfaces.build to continue the thickness
    direction smoothly forward from x = 0.005 instead. Every tabulated
    station (0.5% and up) is untouched.
    """
    k = cl_i / (4 * math.pi)

    def xlogx(x):
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(x == 0, 0.0, x * np.log(x))

    def yc(x):
        x = np.asarray(x, float)
        return -k * (xlogx(1 - x) + xlogx(x))

    def dyc(x):
        x = np.asarray(x, float)
        if k == 0:
            return np.zeros_like(x)     # avoid 0*inf = nan at the ends
        with np.errstate(divide="ignore"):
            return k * (np.log1p(-x) - np.log(x))

    return Camber(name=f"a=1.0 mean line, cl_i = {cl_i:.3f}", yc=yc, dyc=dyc, le_blend_x=0.005)


def camber_4digit(m: float, p: float) -> Camber:
    """NACA 4-digit mean line: m = max camber, p = its location (fractions of
    chord). NACA 4412 -> m = 0.04, p = 0.4."""

    def yc(x):
        x = np.asarray(x, float)
        if m == 0:
            return np.zeros_like(x)
        return np.where(x < p, m / p ** 2 * (2 * p * x - x ** 2),
                        m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x - x ** 2))

    def dyc(x):
        x = np.asarray(x, float)
        if m == 0:
            return np.zeros_like(x)
        return np.where(x < p, 2 * m / p ** 2 * (p - x), 2 * m / (1 - p) ** 2 * (p - x))

    return Camber(name=f"NACA 4-digit mean line, m = {m:.3f}, p = {p:.3f}", yc=yc, dyc=dyc)


def thickness_4digit(t: float) -> Thickness:
    """NACA 4-digit thickness form, t = max thickness (NACA 0012 -> 0.12).
    Standard coefficients, so the TE is slightly open (yt(1) = 0.00126 for t = 0.12)."""

    def yt(x):
        x = np.asarray(x, float)
        return 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2
                        + 0.2843 * x ** 3 - 0.1015 * x ** 4)

    return Thickness(name=f"NACA 4-digit thickness, t = {t:.3f}", yt=yt)


def thickness_table(path) -> Thickness:
    """Read a thickness-form table file and make it a smooth yt(x).

    File format: two columns, x/c and yt/c (half-thickness), from x = 0 to
    x = 1. Lines starting with '#' are comments. An optional comment line
    "# le_radius = <value>" gives the leading-edge radius (fraction of chord).
    See data/naca64_1-012_thickness_abbott.dat.
    """
    path = Path(path)
    text = path.read_text()
    table = np.loadtxt(path, comments="#")
    tok = re.search(r"^#\s*le_radius\s*=\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)\s*$", text, re.M)
    le_radius = float(tok.group(1)) if tok else math.nan
    thk = thickness_fit(table[:, 0], table[:, 1], le_radius)
    thk.name = path.name
    thk.source = str(path)
    return thk


def naca6_section(thickness_file, cl_i: float, name: str) -> Airfoil:
    """6-series section from a thickness table + the a = 1.0 mean line.

    Example, NACA 64(1)-412:
        af = naca6_section('data/naca64_1-012_thickness_abbott.dat', 0.4, 'naca641412')

    Camber and thickness stay separate so the morph can change camber only.
    Surfaces come from surfaces.build.
    """
    cam = camber_a10(cl_i)
    thk = thickness_table(thickness_file)
    return Airfoil(name=name, cam=cam, thk=thk, description=f"{thk.name} + {cam.name}")
