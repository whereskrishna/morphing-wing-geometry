"""PNG plots: baseline vs morphed overlay, camber (mean) lines, skin curvature,
and the baseline against published coordinates.

Colours: baseline = blue, morphed = orange. Both pass colour-blind separation
and 3:1 contrast on the light surface colour (dataviz palette validator).
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .morph import rear_spar  # noqa: E402
from .surfaces import Surfaces, loop  # noqa: E402

SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
TEXT2 = "#52514e"
GRID = "#e4e3df"
SPAR = "#a3a29d"
BASE = "#2a78d6"
MORPH = "#eb6834"
LINE = 2.0
FONT = 10


def _figure(w_in, h_in, rows=1):
    fig, axes = plt.subplots(rows, 1, figsize=(w_in, h_in), facecolor=SURFACE, squeeze=False)
    return fig, axes[:, 0]


def _style(ax, legend=True):
    """Recessive axes: hairline solid grid, muted ticks, no box."""
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, linewidth=0.6, linestyle="-")
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(TEXT2)
        ax.spines[side].set_linewidth(0.6)
    ax.tick_params(colors=TEXT2, labelsize=FONT, direction="out")
    ax.xaxis.label.set_color(TEXT2)
    ax.yaxis.label.set_color(TEXT2)
    ax.title.set_color(TEXT)
    if legend:
        ax.legend(frameon=False, labelcolor=TEXT, fontsize=FONT)


def _spars(ax, mp):
    """Thin grey vertical lines at the front and rear spar (and the hinge, if
    it is behind the rear spar), labelled."""
    xrs = rear_spar(mp)
    lines = [(mp.x_LE, "front spar"), (xrs, "rear spar" if xrs != mp.x_h else "rear spar / hinge")]
    if xrs != mp.x_h:
        lines.append((mp.x_h, "hinge"))
    for xs, label in lines:
        ax.axvline(xs, color=SPAR, linewidth=0.8, zorder=0)
        ax.text(xs + 0.004, 0.02, label, rotation=90, color=SPAR, fontsize=FONT - 1,
                transform=ax.get_xaxis_transform(), va="bottom")


def _save(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def overlay(res, path) -> None:
    """Baseline and morphed airfoil on top of each other, true scale."""
    fig, (ax,) = _figure(11, 3.8)
    B, M = loop(res.base), loop(res.shape)
    ax.plot(B[:, 0], B[:, 1], color=BASE, linewidth=LINE, label="baseline")
    ax.plot(M[:, 0], M[:, 1], color=MORPH, linewidth=LINE, label="morphed")
    _spars(ax, res.morph)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlim(-0.03, 1.03)
    ax.set_xlabel("x / c")
    ax.set_ylabel("y / c")
    ax.set_title(res.name, fontsize=FONT + 1)
    _style(ax)
    _save(fig, path)


def camber_lines(res, path) -> None:
    """Baseline and morphed mean lines. Vertical scale is stretched.

    The mean line is the midpoint between the upper and lower skin points at
    each station. For v1 (and the baseline) that is exactly the camber line.
    v2 has no camber line of its own, so its mean line is shown instead.
    """
    D = res.checks["dense"]

    def mid(S: Surfaces):
        return (S.xu + S.xl) / 2, (S.yu + S.yl) / 2

    bx, by = mid(D["base"])
    mx, my = mid(D["morph"])
    fig, (ax,) = _figure(11, 4.2)
    ax.plot(bx, by, color=BASE, linewidth=LINE, label="baseline")
    ax.plot(mx, my, color=MORPH, linewidth=LINE, label="morphed")
    _spars(ax, res.morph)
    ax.set_xlim(min(bx.min(), mx.min()), 1)
    ax.set_xlabel("x / c")
    ax.set_ylabel("y / c (vertical scale stretched)")
    ax.set_title(f"Camber (mean) lines - {res.name}", fontsize=FONT + 1)
    _style(ax)
    _save(fig, path)


def curvature(res, cfg, path) -> None:
    """Skin curvature vs x for the upper and lower surfaces.

    Curvature in 1/m (uses the chord in cfg.structure.chord_m), convex = +.
    The nose is far more curved than the rest of the skin, so the y range is
    set from x >= 0.02 and the nose runs off the top of the plot.
    """
    D = res.checks["dense"]
    N, x = D["N"], D["x"]
    c_m = cfg.structure.chord_m
    x_from, x_to = 0.02, 1 - cfg.checks.te_exclude
    surf = [("Upper surface", np.arange(N - 1, -1, -1), D["base"].xu, D["morph"].xu),
            ("Lower surface", np.arange(N - 1, 2 * N - 1), D["base"].xl, D["morph"].xl)]
    fig, axes = _figure(11, 7, rows=2)
    for k, (title, idx, xb, xm) in enumerate(surf):
        ax = axes[k]
        kb = D["k_base"][idx] / c_m
        km = D["k_morph"][idx] / c_m
        ax.plot(xb, kb, color=BASE, linewidth=LINE, label="baseline")
        ax.plot(xm, km, color=MORPH, linewidth=LINE, label="morphed")
        _spars(ax, res.morph)
        win = (x >= x_from) & (x <= x_to)
        vals = np.concatenate([kb[win], km[win]])
        lo, hi = np.nanmin(vals), np.nanmax(vals)
        pad = 0.1 * max(hi - lo, np.finfo(float).eps)
        ax.set_ylim(lo - pad, hi + pad)
        ax.set_xlim(0, 1)
        ax.set_ylabel("curvature (1/m), convex +")
        ax.set_title(f"{title}   (nose, x < {x_from:.2f}, off scale)", fontsize=FONT + 1)
        _style(ax, legend=(k == 0))
        if k == 1:
            ax.set_xlabel("x / c")
    fig.suptitle(f"Skin curvature - {res.name}", color=TEXT, fontsize=FONT + 2)
    _save(fig, path)


def baseline_check(S: Surfaces, ref_xy: np.ndarray, title: str, path, label: str = "published (PDAS)") -> None:
    """Built baseline (line) against published coordinates (dots)."""
    fig, (ax,) = _figure(11, 3.8)
    L = loop(S)
    ax.plot(L[:, 0], L[:, 1], color=BASE, linewidth=LINE, label="built here")
    ax.plot(ref_xy[:, 0], ref_xy[:, 1], "o", markersize=4, markerfacecolor=MORPH,
            markeredgecolor=SURFACE, markeredgewidth=1, label=label)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlim(-0.03, 1.03)
    ax.set_xlabel("x / c")
    ax.set_ylabel("y / c")
    ax.set_title(title, fontsize=FONT + 1)
    _style(ax)
    _save(fig, path)


REFLEX = "#1baf7a"        # third categorical slot (aqua); droop uses MORPH
SHADE = "#efeeea"
# what each stopping check measures: (summary column, limit setting, axis label)
_MEASURE = {"strain": ("strain_max_pct", "tol_strain_pct", "peak bending strain (%)"),
            "skin_length": ("arc_max_pct", "tol_arc_pct", "skin length change, worst piece (%)")}


def deflection_range(rng, base: Surfaces, cfg, edge: str, title: str, path) -> None:
    """One page for a deflection range: the two limit shapes over the baseline,
    the check that limits it against the deflection angle, and the rough lift
    change. Droop = orange, reflex = aqua; the allowed range is shaded."""
    lim_d, lim_r = rng["droop"], rng["reflex"]
    stop = next((s for s in lim_d.stopped_by + lim_r.stopped_by if s in _MEASURE), "strain")
    col, tol_key, ylabel = _MEASURE[stop]
    tol = getattr(cfg.checks, tol_key)
    ang = f"{edge.lower()}_angle_deg"

    fig = plt.figure(figsize=(11, 10), facecolor=SURFACE)
    gs = fig.add_gridspec(3, 1, height_ratios=[1.2, 1, 1])

    ax = fig.add_subplot(gs[0])
    B = loop(base)
    ax.plot(B[:, 0], B[:, 1], color=BASE, linewidth=LINE, label="baseline")
    for lim, color in ((lim_d, MORPH), (lim_r, REFLEX)):
        if lim.res is not None:
            M = loop(lim.res.shape)
            ax.plot(M[:, 0], M[:, 1], color=color, linewidth=LINE,
                    label=f"largest {lim.direction}: {lim.angle_deg:+.3f} deg")
    _spars(ax, (lim_d.res or lim_r.res).morph if (lim_d.res or lim_r.res) else cfg.morph)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlim(-0.03, 1.03)
    ax.set_xlabel("x / c")
    ax.set_ylabel("y / c")
    ax.set_title(title, fontsize=FONT + 1)
    _style(ax, legend=False)
    ax.legend(frameon=False, labelcolor=TEXT, fontsize=FONT, loc="upper center",
              bbox_to_anchor=(0.5, -0.22), ncol=3)       # below the axis, clear of the airfoil and spar labels

    for k, (ycol, label) in enumerate(((col, ylabel), ("dcl_estimate", "lift change, rough estimate (cl)"))):
        ax = fig.add_subplot(gs[k + 1])
        ax.axvspan(lim_r.angle_deg, lim_d.angle_deg, color=SHADE, zorder=0, linewidth=0)
        for lim, color in ((lim_d, MORPH), (lim_r, REFLEX)):
            pts = [(0.0, 0.0)] + [(r[ang], r[ycol]) for a, r, f in lim.tried if r is not None]
            pts.sort()
            a, v = np.array(pts).T
            ax.plot(a, v, color=color, linewidth=LINE, marker="o", markersize=4,
                    markeredgecolor=SURFACE, markeredgewidth=1, label=lim.direction)
        if k == 0:
            ax.axhline(tol, color=TEXT2, linewidth=1, linestyle="--")
            ax.text(ax.get_xlim()[0], tol, f"  limit {tol:.3f}", color=TEXT2, fontsize=FONT - 1, va="bottom")
            ax.set_title("Allowed range shaded. Points: every shape the search tried", fontsize=FONT)
        else:
            ax.axhline(0, color=SPAR, linewidth=0.8)
            ax.set_title("Rough thin-airfoil estimate: good near 0, too high at large deflections "
                         "(the real flow separates)", fontsize=FONT)
            ax.set_xlabel(f"{edge} deflection, equivalent flap angle (deg), + = droop")
        ax.set_ylabel(label)
        _style(ax, legend=(k == 0))
    _save(fig, path)


def spar_study(rows, cfg, answer, title: str, path) -> None:
    """Deflection range, rough lift change at the limits, and rear spar height,
    against the rear spar station (hinge at the spar)."""
    st = cfg.spar_study
    x = np.array([r["x_rear_spar"] for r in rows])
    fig, axes = _figure(11, 10, rows=3)

    ax = axes[0]
    ax.plot(x, [r["droop_max_deg"] for r in rows], color=MORPH, linewidth=LINE, marker="o", markersize=6,
            markeredgecolor=SURFACE, markeredgewidth=1, label="largest droop")
    ax.plot(x, [r["reflex_max_deg"] for r in rows], color=REFLEX, linewidth=LINE, marker="o", markersize=6,
            markeredgecolor=SURFACE, markeredgewidth=1, label="largest reflex")
    for need, label in ((st.need_droop_deg, "droop needed"), (st.need_reflex_deg, "reflex needed")):
        ax.axhline(need, color=TEXT2, linewidth=1, linestyle="--")
        ax.text(x.min(), need, f" {label} {need:+.1f} deg (placeholder)", color=TEXT2, fontsize=FONT - 1,
                va="bottom" if need > 0 else "top")
    ax.axhline(0, color=SPAR, linewidth=0.8)
    if answer is not None:
        ax.axvline(answer, color=TEXT, linewidth=1)
        ax.text(answer, 0.97, f" furthest aft\n x = {answer:.3f}", color=TEXT, fontsize=FONT - 1,
                transform=ax.get_xaxis_transform(), va="top")
    ax.set_ylabel("equivalent flap angle (deg)")
    ax.set_title(title, fontsize=FONT + 1)
    _style(ax)
    ax.legend(frameon=False, labelcolor=TEXT, fontsize=FONT, loc="lower right")

    ax = axes[1]
    ax.plot(x, [r["lift_per_deg_estimate"] for r in rows], color=BASE, linewidth=LINE, marker="o", markersize=6,
            markeredgecolor=SURFACE, markeredgewidth=1)
    ax.set_ylim(bottom=0)
    ax.set_ylabel("lift change per degree (cl/deg)")
    ax.set_title("Rough thin-airfoil estimate at small deflections: a shorter morphing TE gives less lift per degree",
                 fontsize=FONT)
    _style(ax, legend=False)

    ax = axes[2]
    ax.plot(x, [r["rear_spar_height_mm"] for r in rows], color=BASE, linewidth=LINE, marker="o", markersize=6,
            markeredgecolor=SURFACE, markeredgewidth=1)
    ax.set_ylabel("rear spar height (mm)")
    ax.set_xlabel("rear spar station x / c   (the TE bends from the spar back)")
    ax.set_title("Wing box: a spar further back is shorter but the box is bigger", fontsize=FONT)
    _style(ax, legend=False)
    for a in axes:
        a.set_xlim(x.min() - 0.02, x.max() + 0.02)
    _save(fig, path)


def _numeric(values):
    try:
        return [float(v) for v in values], None
    except ValueError:
        return list(range(len(values))), list(values)


def _real_stop(text: str) -> bool:
    """True if a check stopped the search (not the solver or the search cap)."""
    return not (text.startswith("morph did not solve") or text.startswith("nothing failed"))


def _short_stop(text: str) -> str:
    if text.startswith("morph did not solve"):
        return "solver limit"
    if text.startswith("nothing failed"):
        return "search cap"
    return text.replace("_", " ")


def variable_study(groups, cfg, edge: str, default: dict, title: str, path, need=None) -> None:
    """One page for a one-setting-at-a-time study.

    Top: how far each setting moves the largest droop (orange) and reflex
    (aqua) over the values tried, biggest effect first. Below: one small
    panel per setting, all on the same angle scale; a point stopped by a
    different check from the rest of its panel is labelled, and a hollow
    point is a solver limit or the search cap, not a check. For the TE, a
    last panel puts the chord, skin thickness and strain limit runs against
    the one number the strain check sees (bend_radius_limit_c).
    groups: [(label, rows)], rows from deflection.vary_one.
    default: {"droop": deg, "reflex": deg} at the settings in cfg.
    need: (droop deg, reflex deg) to mark, or None.
    """
    from .deflection import spread, unsolved

    order = sorted(groups, key=lambda g: spread(g[1]))   # biggest effect last = top row
    strain_like = [g for g in groups if g[1][0]["setting"] in
                   ("structure.t_skin_mm", "structure.chord_m", "checks.tol_strain_pct")]
    panels = len(groups) + (1 if edge == "TE" and len(strain_like) >= 2 else 0)
    ncol = 3
    nrow = -(-panels // ncol)
    top_h = 0.42 * len(groups) + 2.2
    fig = plt.figure(figsize=(12, top_h + 3.1 * nrow), facecolor=SURFACE, layout="constrained")
    sub_top, sub_grid = fig.subfigures(2, 1, height_ratios=[top_h, 3.1 * nrow])   # separate layouts:
    sub_top.set_facecolor(SURFACE)                       # the long labels on top don't squeeze the panels
    sub_grid.set_facecolor(SURFACE)
    kw = dict(color=TEXT, fontsize=FONT - 1)

    # --- top: spread of each setting's effect (a bar per setting and direction) ---
    ax = sub_top.subplots()
    for i, (label, rows) in enumerate(order):
        for key, color in (("droop_max_deg", MORPH), ("reflex_max_deg", REFLEX)):
            v = [r[key] for r in rows if not unsolved(r)] or [0.0]
            ax.hlines(i, min(v), max(v), color=color, linewidth=7, capstyle="round")
            if max(v) - min(v) < 0.2:                     # no effect: show a dot, not nothing
                ax.plot([min(v)], [i], marker="o", markersize=8, color=color, markeredgecolor=SURFACE)
    ax.set_yticks(range(len(order)), [g[0] for g in order])
    top = len(order) - 0.4
    ax.set_ylim(-0.7, top + 1.2)                           # two headroom rows for the labels
    frac = (top + 0.7) / (top + 1.9)                       # vertical lines stop below the labels
    for key in ("droop", "reflex"):
        ax.axvline(default[key], color=TEXT, linewidth=1, ymax=frac)
        ax.text(default[key], top + 0.9, f" default {default[key]:+.1f} ", va="center",
                ha="left" if key == "droop" else "right", **kw)
    if need is not None:
        for v in need:
            ax.axvline(v, color=TEXT2, linewidth=1, linestyle="--", ymax=frac)
            ax.text(v, top + 0.35, f" needed {v:+.0f} (placeholder) ", va="center", color=TEXT2,
                    fontsize=FONT - 1, ha="left" if v > 0 else "right")
    ax.axvline(0, color=SPAR, linewidth=0.8)
    ax.set_xlabel(f"largest {edge} deflection that passes every check, equivalent flap angle (deg)")
    ax.set_title(title + "\nEach bar: the spread over the values tried, biggest effect at the top. "
                 "Orange = droop, aqua = reflex", fontsize=FONT + 1, loc="left", color=TEXT)
    _style(ax, legend=False)
    ax.grid(False, axis="y")

    # --- one panel per setting, same angle scale ---
    grid = sub_grid.subplots(nrow, ncol, squeeze=False).ravel()
    axes = list(grid[:panels])
    for spare in grid[panels:]:
        spare.axis("off")
    lo = min(min(r["reflex_max_deg"] for r in g[1]) for g in groups)
    hi = max(max(r["droop_max_deg"] for r in g[1]) for g in groups)
    pad = 0.12 * (hi - lo)
    for k, (ax, (label, rows)) in enumerate(zip(axes, groups)):
        x, cats = _numeric([r["value"] for r in rows])
        for key, stop, color in (("droop_max_deg", "droop_stopped_by", MORPH),
                                 ("reflex_max_deg", "reflex_stopped_by", REFLEX)):
            y = [r[key] for r in rows]
            ax.plot(x, y, color=color, linewidth=LINE, zorder=2)
            for xi, yi, r in zip(x, y, rows):
                real = _real_stop(r[stop])
                ax.plot([xi], [yi], marker="o", markersize=6 if real else 7, zorder=3,
                        color=color if real else SURFACE, markeredgecolor=SURFACE if real else color,
                        markeredgewidth=1 if real else 1.6)
            # label what stopped the search where it differs from the rest of the panel,
            # once for each run of neighbouring points with the same reason
            usual = max({r[stop] for r in rows}, key=[r[stop] for r in rows].count)
            runs = []
            for xi, yi, r in zip(x, y, rows):
                if r[stop] == usual:
                    continue
                if runs and runs[-1][0] == r[stop] and runs[-1][-1] == rows.index(r) - 1:
                    runs[-1][1].append(xi), runs[-1][2].append(yi)
                    runs[-1][-1] = rows.index(r)
                else:
                    runs.append([r[stop], [xi], [yi], rows.index(r)])
            up = key[0] == "d"
            for text, xs, ys, _ in runs:
                ax.annotate(_short_stop(text), (float(np.mean(xs)), max(ys) if up else min(ys)),
                            xytext=(0, 9 if up else -9), textcoords="offset points", ha="center",
                            va="bottom" if up else "top", color=TEXT2, fontsize=FONT - 2)
        for xi, r in zip(x, rows):
            if r["is_default"] == "yes":
                ax.axvline(xi, color=SPAR, linewidth=0.8, zorder=0)
                ax.text(xi, 0.98, " default", color=TEXT2, fontsize=FONT - 2,
                        transform=ax.get_xaxis_transform(), va="top")
        if cats is not None:
            ax.set_xticks(x, cats)
            ax.set_xlim(-0.5, len(x) - 0.5)
        if need is not None:
            for v in need:
                ax.axhline(v, color=TEXT2, linewidth=1, linestyle="--")
        ax.axhline(0, color=SPAR, linewidth=0.8)
        ax.set_ylim(lo - pad, hi + pad)
        ax.set_title(label, fontsize=FONT, loc="left", color=TEXT)
        _style(ax, legend=False)
        if k % ncol == 0:
            ax.set_ylabel("largest angle (deg)")

    if panels > len(groups):
        ax = axes[-1]
        marks = {"structure.t_skin_mm": ("o", "changing skin thickness"), "structure.chord_m": ("s", "changing chord"),
                 "checks.tol_strain_pct": ("^", "changing strain limit")}
        for label, rows in strain_like:
            m, name = marks[rows[0]["setting"]]
            xb = [r["bend_radius_limit_c"] for r in rows]
            for key, stop, color in (("droop_max_deg", "droop_stopped_by", MORPH),
                                     ("reflex_max_deg", "reflex_stopped_by", REFLEX)):
                for xi, r in zip(xb, rows):
                    real = _real_stop(r[stop])
                    ax.plot([xi], [r[key]], linestyle="none", marker=m, markersize=7,
                            color=color if real else SURFACE, markeredgecolor=TEXT2 if real else color,
                            markeredgewidth=0.6 if real else 1.6)
            ax.plot([], [], linestyle="none", marker=m, markersize=7, color=SURFACE, markeredgecolor=TEXT2, label=name)
        ax.axhline(0, color=SPAR, linewidth=0.8)
        ax.set_ylim(lo - pad, hi + pad)
        ax.set_xlabel("skin thickness / (2 x strain limit) / chord")
        ax.set_title("The same runs against one number: one curve", fontsize=FONT, loc="left", color=TEXT)
        _style(ax, legend=False)
        where = grid[panels] if panels < len(grid) else ax          # the spare cell, if there is one
        h, lab = ax.get_legend_handles_labels()
        where.legend(h, lab, frameon=False, labelcolor=TEXT, fontsize=FONT - 1, loc="center left",
                     title="marker = which setting was changed\nhollow (any panel) = solver limit or\n"
                           "search cap, not a check", title_fontsize=FONT - 2, alignment="left")
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


# ---- every flap angle x hinge x shape exponent (scripts/flap_study.py) ----

def _ramp(hex_color):
    """Sequential colour map: one hue, from near the surface colour to full strength."""
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("ramp", ["#f3f2ee", hex_color], N=256)


def _heat(ax, xs, ns, vals, cmap, fmt, title, outline=None):
    """Hinge (columns) x n (rows) map with the value written in each cell."""
    V = np.array(vals, dtype=float)
    lo, hi = np.nanmin(V), np.nanmax(V)
    hi = hi if hi > lo else lo + 1
    ax.imshow(V, cmap=cmap, vmin=lo, vmax=hi, aspect="auto", origin="lower", interpolation="nearest")
    for i in range(len(ns)):
        for j in range(len(xs)):
            v = V[i, j]
            if np.isnan(v):
                ax.text(j, i, "fails", ha="center", va="center", fontsize=FONT - 2, color=TEXT2)
                continue
            r, g, b, _ = cmap((v - lo) / (hi - lo))
            dark = 0.2126 * r + 0.7152 * g + 0.0722 * b < 0.45          # light text on dark cells only
            ax.text(j, i, fmt(v), ha="center", va="center", fontsize=FONT - 2,
                    color=SURFACE if dark else TEXT)
            if outline is not None and outline[i][j]:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor=TEXT, linewidth=1.6))
    ax.set_xticks(range(len(xs)), [f"{x:.3f}".rstrip("0") for x in xs], fontsize=FONT - 2)
    ax.set_yticks(range(len(ns)), [f"{n:g}" for n in ns], fontsize=FONT - 1)
    ax.set_xlabel("hinge x / c   (rear spar at the hinge)", color=TEXT2, fontsize=FONT - 1)
    ax.set_ylabel("shape exponent n", color=TEXT2, fontsize=FONT - 1)
    ax.set_title(title, fontsize=FONT, color=TEXT, loc="left")
    for side in ax.spines.values():
        side.set_visible(False)
    ax.tick_params(colors=TEXT2, length=0)


def flap_study_limits(limits, cfg, fc, title: str, path) -> None:
    """Largest droop, largest reflex, total range and lift per degree at
    cruise, as maps over hinge x n."""
    xs = sorted({r["x_h"] for r in limits})
    ns = sorted({r["n"] for r in limits})
    by = {(r["x_h"], r["n"]): r for r in limits}

    def grid(key):
        return [[by[(x, n)][key] for x in xs] for n in ns]

    need = cfg.spar_study
    reach = [[by[(x, n)]["reaches_need"] for x in xs] for n in ns]
    fig, axes = plt.subplots(2, 2, figsize=(15, 9.5), facecolor=SURFACE, layout="constrained")
    _heat(axes[0, 0], xs, ns, grid("droop_max_deg"), _ramp(MORPH), lambda v: f"{v:.1f}",
          "Largest droop (deg): passes every check up to here")
    _heat(axes[0, 1], xs, ns, np.abs(grid("reflex_max_deg")), _ramp(REFLEX),
          lambda v: f"-{v:.1f}" if v >= 0.05 else "0.0",
          "Largest reflex (deg)")
    _heat(axes[1, 0], xs, ns, grid("range_deg"), _ramp(BASE), lambda v: f"{v:.1f}",
          f"Total range, droop to reflex (deg). Outlined: reaches the needed {need.need_droop_deg:+g} / "
          f"{need.need_reflex_deg:+g} deg (placeholder)", outline=reach)
    step = cfg.flap_study.angle_step_deg
    lift = [[by[(x, n)]["lift_per_deg_n_per_m"] if by[(x, n)]["droop_max_deg"] >= step else np.nan for x in xs]
            for n in ns]                                  # blank where even the first step fails a check
    _heat(axes[1, 1], xs, ns, lift, _ramp(BASE), lambda v: f"{v:.0f}",
          f"Lift change at cruise, N per metre of span per degree (rough: low-speed theory,\n"
          f"first few degrees only; q = {fc['q_pa']:.0f} Pa, chord {cfg.structure.chord_m:g} m)")
    fig.suptitle(f"{title}\nchord {cfg.structure.chord_m:g} m, skin {cfg.structure.t_skin_mm:g} mm, strain limit "
                 f"{cfg.checks.tol_strain_pct:g} % (placeholders); cruise Mach {fc['mach']:.2f}",
                 color=TEXT, fontsize=FONT + 1)
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


# what stopped a shape, for the map: (label, colour); the first that applies is shown
_STOPS = [("passes every check", "#b9d3f0"), ("morph did not solve", SPAR), ("curvature_join", "#7f56d9"),
          ("strain", MORPH), ("other check", TEXT2)]


def _stop_kind(row) -> int:
    if row["passes"]:
        return 0
    failed = row["failed"]
    for i, (key, _) in enumerate(_STOPS[1:-1], 1):
        if key in failed:
            return i
    return len(_STOPS) - 1


def flap_study_map(grid, limits, cfg, title: str, path) -> None:
    """Every shape tried, one panel per n: angle against hinge, coloured by
    whether it passes or what stopped it, with the exact limits on top."""
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    xs = sorted({r["x_h"] for r in grid})
    ns = sorted({r["n"] for r in grid})
    angles = sorted({r["angle_deg"] for r in grid})
    step = min(np.diff(angles)) if len(angles) > 1 else 1.0
    a_all = np.arange(angles[0], angles[-1] + step / 2, step)
    ia = {a: i for i, a in enumerate(np.round(a_all, 9))}
    cols = 4
    rows = int(np.ceil(len(ns) / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(15, 4.6 * rows + 0.8), facecolor=SURFACE, layout="constrained",
                             sharex=True, sharey=True, squeeze=False)
    cmap = ListedColormap([c for _, c in _STOPS])
    need = cfg.spar_study
    for k, n in enumerate(ns):
        ax = axes[k // cols, k % cols]
        Z = np.full((len(a_all), len(xs)), np.nan)
        for r in grid:
            if r["n"] == n:
                Z[ia[round(r["angle_deg"], 9)], xs.index(r["x_h"])] = _stop_kind(r)
        dx = np.diff(xs).min() if len(xs) > 1 else 0.05
        ax.pcolormesh(np.r_[np.array(xs) - dx / 2, xs[-1] + dx / 2], np.r_[a_all - step / 2, a_all[-1] + step / 2],
                      np.ma.masked_invalid(Z), cmap=cmap, vmin=-0.5, vmax=len(_STOPS) - 0.5,
                      edgecolor=SURFACE, linewidth=0.2)
        lim = sorted((r for r in limits if r["n"] == n), key=lambda r: r["x_h"])
        for key in ("droop_max_deg", "reflex_max_deg"):
            ax.plot([r["x_h"] for r in lim], [r[key] for r in lim], color=TEXT, linewidth=1.4, marker="o",
                    markersize=3)
        for v in (need.need_droop_deg, need.need_reflex_deg):
            ax.axhline(v, color=TEXT, linewidth=0.9, linestyle="--")
        ax.set_title(f"n = {n:g}", fontsize=FONT, color=TEXT, loc="left")
        ax.set_facecolor(SURFACE)
        ax.tick_params(colors=TEXT2, labelsize=FONT - 1)
        for side in ax.spines.values():
            side.set_visible(False)
        if k % cols == 0:
            ax.set_ylabel("equivalent flap angle (deg)", color=TEXT2, fontsize=FONT - 1)
        if k // cols == rows - 1:
            ax.set_xlabel("hinge x / c", color=TEXT2, fontsize=FONT - 1)
    for k in range(len(ns), rows * cols):
        axes[k // cols, k % cols].set_visible(False)
    seen = {_stop_kind(r) for r in grid}
    handles = [Patch(facecolor=c, label=lab.replace("_", " ")) for i, (lab, c) in enumerate(_STOPS) if i in seen]
    handles += [plt.Line2D([], [], color=TEXT, linewidth=1.4, marker="o", markersize=3,
                           label="exact largest droop / reflex"),
                plt.Line2D([], [], color=TEXT, linewidth=0.9, linestyle="--",
                           label=f"needed {need.need_droop_deg:+g} / {need.need_reflex_deg:+g} deg (placeholder)")]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, frameon=False, labelcolor=TEXT, fontsize=FONT - 1)
    fig.suptitle(f"{title}\nEach cell is one shape. Colour: passes, or the check that stops it (solver first, "
                 f"then curvature join, then strain)", color=TEXT, fontsize=FONT + 1)
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
