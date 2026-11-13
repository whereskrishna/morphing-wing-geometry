"""Rough lift change from a shape change, by thin-airfoil theory, and the flight condition.

    dcl = lift_change(base, morph)     # surfaces from surfaces.build or a morph

Thin-airfoil theory (thin, inviscid, low speed) gives the lift coefficient
of a mean line z(x), chord 0 to 1, at angle of attack alpha (radians,
measured from the x-axis):

    cl = 2 pi alpha + cl_0,    cl_0 = -4 * integral from 0 to 1 of sqrt(x / (1 - x)) dz

So a morph changes the lift at a fixed angle of attack by the change in
cl_0. The weight sqrt(x / (1 - x)) grows towards the trailing edge, which is
why trailing-edge flaps are so effective. For a straight flap hinged at x_h
and turned down by a small angle d this gives the textbook result
dcl = 2 d (pi - t + sin t), with cos t = 1 - 2 x_h.

This is an ESTIMATE for comparing shapes with each other (more or less lift,
how much a hinge position costs). It ignores thickness, viscosity and
compressibility, and it reads too high at large deflections, where the flow
separates. XFOIL or CFD give the real numbers.
"""
from __future__ import annotations

import numpy as np

from .surfaces import Surfaces


def _weight_integral(u):
    """Integral of sqrt(u / (1 - u)) from 0 to u."""
    u = np.clip(u, 0.0, 1.0)
    return np.arcsin(np.sqrt(u)) - np.sqrt(u * (1 - u))


def zero_lift_cl(S: Surfaces) -> float:
    """cl_0 of a shape: its lift coefficient at zero angle of attack (x-axis).

    The mean line is the midpoint of the upper and lower skin points at each
    station (exactly the camber line for v1 and the baseline). Its x range,
    nose to TE, is mapped onto 0 to 1 (v2 moves the TE forward a little).
    z is taken as straight between points, which makes the integral exact
    for that piecewise-straight line.
    """
    x = (S.xu + S.xl) / 2
    z = (S.yu + S.yl) / 2
    u = (x - x[0]) / (x[-1] - x[0])
    du = np.diff(u)
    dF = np.diff(_weight_integral(u))
    mid = np.clip((u[:-1] + u[1:]) / 2, 0.0, 1 - 1e-15)
    # mean weight on each piece; the weight itself where a piece has no width
    w = np.where(np.abs(du) > 1e-14, dF / np.where(du == 0, 1, du), np.sqrt(mid / (1 - mid)))
    return float(-4 * np.sum(w * np.diff(z)))


def lift_change(base: Surfaces, morph: Surfaces) -> float:
    """Change in lift coefficient at the same angle of attack (+ = more lift)."""
    return zero_lift_cl(morph) - zero_lift_cl(base)


# ---- flight condition (config flight:) ----

GAMMA = 1.4          # air, ratio of specific heats
R_AIR = 287.05287    # air, gas constant J/(kg K)


def flight_condition(fl) -> dict:
    """Numbers that follow from the flight settings.

    mach                 velocity / speed of sound
    q_pa                 dynamic pressure, 1/2 rho V^2 (uses the density as given)
    temperature_k        from the speed of sound, a = sqrt(1.4 R T)
    density_from_p_t     pressure / (R T): what the density should be if the
                         pressure and speed of sound are right (a cross-check)
    reynolds_per_m       rho V / mu, for one metre of chord (Sutherland viscosity)
    prandtl_glauert      1 / sqrt(1 - M^2): the simple compressibility factor
                         on lift. It is only good well below the speed where
                         shocks form, so at Mach 0.8 it is a rough guide only.
    """
    a, V, rho, p = fl.speed_of_sound_m_s, fl.velocity_m_s, fl.density_kg_m3, fl.pressure_pa
    M = V / a
    T = a ** 2 / (GAMMA * R_AIR)
    mu = 1.458e-6 * T ** 1.5 / (T + 110.4)
    return {"mach": M, "q_pa": 0.5 * rho * V ** 2, "temperature_k": T,
            "density_from_p_t": p / (R_AIR * T), "reynolds_per_m": rho * V / mu,
            "prandtl_glauert": 1 / np.sqrt(1 - M ** 2) if M < 1 else float("nan")}


def lift_per_span(dcl: float, q_pa: float, chord_m: float) -> float:
    """Lift change per metre of span (N/m) from a lift coefficient change."""
    return dcl * q_pa * chord_m
