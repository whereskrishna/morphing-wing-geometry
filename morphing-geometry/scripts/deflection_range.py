"""Deflection range: how far the morphing edge can droop and reflex before a check fails.

Change the settings below, then:
    python scripts/deflection_range.py                         # NASA SC(2)-0714 (the default)
    python scripts/deflection_range.py config/naca641412.yaml  # any case file, e.g. the NACA 64(1)-412

For the shape settings (method, hinge, spars, n) it tries bigger and bigger
deflections, first droop, then reflex, and finds the largest of each that
still passes every check. Deflections are equivalent flap angles: a straight
flap hinged at x_h with the same tip drop. The search settings are under
range: in config/default.yaml.

Writes <output dir>/range/:
    <shape>.dat, <shape>_cad_mm.txt   the two limit shapes (largest droop, largest reflex)
    <name>.csv                        every shape the search tried, with all check values
    <name>.png                        limit shapes, the limiting check and the rough lift change
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from morphgeom import load_baseline, load_config, write_outputs  # noqa: E402
from morphgeom import plots  # noqa: E402
from morphgeom.deflection import describe, find_range, lift_per_degree  # noqa: E402
from morphgeom.export import write_summary_csv  # noqa: E402
from morphgeom.morph import rear_spar  # noqa: E402
from morphgeom.surfaces import build, stations  # noqa: E402

cfg = load_config(sys.argv[1] if len(sys.argv) > 1 else None)   # optional case file
cfg.morph.method = "constant_length"   # "constant_length" (v2, skins keep their length) or "camber" (v1)
edge = "TE"                            # "TE" or "LE"
# cfg.morph.x_h = 0.65                 # hinge: where the TE starts to bend
# cfg.morph.x_rear_spar = 0.60         # rear spar, at or ahead of the hinge (None = at the hinge)
# cfg.morph.n = 3
# cfg.range.ignore = ["skin_length"]   # e.g. to see where v1 would stop if its skins could stretch

af = load_baseline(cfg)
mp = cfg.morph
v2 = mp.method == "constant_length"
where = f"hinge x_h {mp.x_h:.3f}" if edge == "TE" else f"front spar x_LE {mp.x_LE:.3f}"
if edge == "TE" and rear_spar(mp) != mp.x_h:
    where += f", rear spar {rear_spar(mp):.3f}"
name = (f"{af.name}_{edge}range_" + (f"xh{mp.x_h:.3f}" if edge == "TE" else f"xle{mp.x_LE:.3f}")
        + (f"_rs{rear_spar(mp):.3f}" if edge == "TE" and rear_spar(mp) != mp.x_h else "")
        + f"_n{mp.n:g}" + ("_constlen" if v2 else ""))

rng = find_range(cfg, af, edge)
outdir = Path(cfg.output.dir) / "range"
c_mm = 1000 * cfg.structure.chord_m

print(f"\nDeflection range: {af.name}, {edge}, {'constant_length (v2)' if v2 else 'camber (v1)'}, {where}, n = {mp.n:g}")
print(f"  chord {cfg.structure.chord_m:.3f} m, skin {cfg.structure.t_skin_mm:.3f} mm, "
      f"strain limit {cfg.checks.tol_strain_pct:.3f} %, skin length limit {cfg.checks.tol_arc_pct:.3f} %"
      "  (placeholders unless set)")
if cfg.range.ignore:
    print(f"  not stopping at: {', '.join(cfg.range.ignore)}")
print(f"  {'':<16s} {'angle (deg)':>12s} {'tip move (mm)':>14s} {'peak strain (%)':>16s}   stopped by")
for lim in rng.values():
    if lim.res is not None:
        write_outputs(lim.res, cfg, outdir, make_plots=False)
        vals = f"{lim.angle_deg:12.3f} {lim.delta * c_mm:14.3f} {lim.res.row['strain_max_pct']:16.3f}"
    else:
        vals = f"{0.0:12.3f} {'-':>14s} {'-':>16s}"
    print(f"  largest {lim.direction:<8s} {vals}   {describe(lim, cfg.range.max_deg)}")
print(f"  (each angle passes every check; the true limit is less than {cfg.range.tol_deg:g} deg further out;")
print("   tip move: + = down)")
print(f"  lift change per degree, small deflections: {lift_per_degree(rng['droop']):.3f} cl/deg "
      "(rough thin-airfoil estimate, README explains)")

outdir.mkdir(parents=True, exist_ok=True)
rows = sorted((r for lim in rng.values() for a, r, f in lim.tried if r is not None),
              key=lambda r: r[f"{edge.lower()}_angle_deg"])
write_summary_csv(outdir / f"{name}.csv", rows)
base = build(af.cam, af.thk, stations(cfg.points.n_per_surface, cfg.points.spacing))
plots.deflection_range(rng, base, cfg, edge, f"Deflection range - {name}", outdir / f"{name}.png")
print(f"\nFiles written to {outdir}")
