"""Make one morphed shape: coordinate file, CAD file, plots, and the check report.

Change the settings below, then:
    python scripts/make_shape.py                         # NASA SC(2)-0714 (the default)
    python scripts/make_shape.py config/naca641412.yaml  # any case file, e.g. the NACA 64(1)-412
Everything not set here comes from config/default.yaml (or the case file).

Writes <output dir>/shapes/<name>.dat, <name>_cad_mm.txt and three PNG plots.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from morphgeom import load_config, make_shape, write_outputs  # noqa: E402
from morphgeom.export import print_report  # noqa: E402

cfg = load_config(sys.argv[1] if len(sys.argv) > 1 else None)   # optional case file
cfg.morph.method = "camber"      # "camber" (v1) or "constant_length" (v2: skins keep their length)
cfg.morph.delta_TE = 0.02        # TE tip drop, fraction of chord (+ = droop)
cfg.morph.delta_LE = 0.00        # LE tip drop, fraction of chord (+ = droop)
cfg.morph.n = 3                  # 3 = curvature matches at the spars
# cfg.morph.x_h = 0.65           # hinge: where the TE starts to bend (PLACEHOLDER in default.yaml)
# cfg.morph.x_rear_spar = 0.60   # rear spar, at or ahead of the hinge (None = at the hinge)
# cfg.morph.x_LE = 0.15          # front spar (PLACEHOLDER in default.yaml)

res = make_shape(cfg)
outdir = Path(cfg.output.dir) / "shapes"
write_outputs(res, cfg, outdir)
print_report(res, cfg)
print(f"\nFiles written to {outdir}")
