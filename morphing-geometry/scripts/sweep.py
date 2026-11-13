"""Run every combination of the sweep settings and write summary.csv.

Set the lists in config/default.yaml (sweep:), or below, then:
    python scripts/sweep.py                         # NASA SC(2)-0714 (the default)
    python scripts/sweep.py config/naca641412.yaml  # any case file, e.g. the NACA 64(1)-412

Writes <output dir>/sweep/: one .dat (+ CAD file) per shape and summary.csv with
every check value and PASS/FAIL.
"""
import itertools
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from morphgeom import copy_config, load_baseline, load_config, make_shape, write_outputs  # noqa: E402
from morphgeom.export import write_summary_csv  # noqa: E402

cfg = load_config(sys.argv[1] if len(sys.argv) > 1 else None)   # optional case file
# cfg.sweep.delta_TE = [-0.02, -0.01, 0.0, 0.01, 0.02, 0.03]   # example of changing a list here

af = load_baseline(cfg)          # load once, reuse for every shape
outdir = Path(cfg.output.dir) / "sweep"
sw = cfg.sweep
# same order as before: delta_TE changes fastest, method slowest
grid = [combo[::-1] for combo in itertools.product(sw.method, sw.n, sw.x_rear_spar, sw.x_h, sw.x_LE,
                                                    sw.delta_LE, sw.delta_TE)]
print(f"\nSweep: {len(grid)} shapes")
print(f"  {'name':<57s} {'TE (deg)':>9s} {'LE (deg)':>9s} {'max skin len (%)':>16s} "
      f"{'max strain (%)':>16s} {'off v1 (mm)':>12s}   all checks")

rows = []
for d_te, d_le, x_le, x_h, x_rs, n, method in grid:
    c = copy_config(cfg)
    c.morph.delta_TE, c.morph.delta_LE = d_te, d_le
    c.morph.x_LE, c.morph.x_h, c.morph.x_rear_spar, c.morph.n, c.morph.method = x_le, x_h, x_rs, n, method
    res = make_shape(c, af)
    write_outputs(res, c, outdir, sw.plots)
    rows.append(res.row)
    r = res.row
    print(f"  {r['name']:<57s} {r['te_angle_deg']:9.3f} {r['le_angle_deg']:9.3f} {r['arc_max_pct']:16.3f} "
          f"{r['strain_max_pct']:16.3f} {r['off_v1_mm']:12.3f}   {'PASS' if r['all_checks'] else 'FAIL'}")

csv = outdir / "summary.csv"
write_summary_csv(csv, rows)
print(f"\nWrote {len(rows)} shapes and {csv}")
