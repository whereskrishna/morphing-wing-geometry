"""Build the unmorphed baseline, write its files, check it against published data.

    python scripts/make_baseline.py                         # NASA SC(2)-0714 (the default)
    python scripts/make_baseline.py config/naca641412.yaml  # any case file, e.g. the NACA 64(1)-412

Writes <output dir>/baseline/:
    <name>_baseline.dat          coordinates for XFOIL / CFD (phase 2 baseline)
    <name>_baseline_cad_mm.txt   SolidWorks XYZ file
    <name>_baseline_check.png    built shape vs the published points
Prints how far the built baseline is from the published coordinates:
    NACA 6-series: against the PDAS 64(1)-412 table in data/
    loaded .dat:   against the file's own points (after smoothing)
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from morphgeom import load_config, make_shape, write_outputs  # noqa: E402
from morphgeom import plots  # noqa: E402
from morphgeom.naca import naca6_section  # noqa: E402
from morphgeom.surfaces import build, loop, read_dat  # noqa: E402

cfg = load_config(sys.argv[1] if len(sys.argv) > 1 else None)
cfg.morph.delta_LE = 0.0
cfg.morph.delta_TE = 0.0
res = make_shape(cfg)
af = res.af
outdir = Path(cfg.output.dir) / "baseline"
write_outputs(res, cfg, outdir, make_plots=False)
mm = 1000 * cfg.structure.chord_m

if cfg.baseline.source == "naca6":
    # Compare with the published table, point by point, at its 26 stations.
    ref = read_dat(ROOT / "data" / "naca64_1-412_pdas_reference.dat")
    xref = np.concatenate([[0, 0.005, 0.0075, 0.0125, 0.025, 0.05, 0.075, 0.1], np.linspace(0.15, 1.0, 18)])
    print("\nNACA 64(1)-412 baseline vs PDAS published coordinates (26 stations per surface)")
    print(f"  {'thickness table used':<46s} {'max diff (% chord)':>18s} {'max diff (mm)':>14s}")
    for table in (cfg.baseline.thickness_file, ROOT / "data" / "naca64_1-012_thickness_pdas.dat"):
        a = naca6_section(table, cfg.baseline.cl_design, cfg.baseline.name)
        S = build(a.cam, a.thk, xref)
        d = np.max(np.linalg.norm(loop(S) - ref.xy, axis=1))
        print(f"  {a.thk.name:<46s} {100 * d:18.3f} {d * mm:14.3f}")
    print("  (the first line is the baseline used everywhere; the second checks the mean line")
    print("   and the construction, using the same thickness table PDAS used)")
    title, label = "NACA 64(1)-412 baseline: built here vs published (PDAS)", "published (PDAS)"
else:
    ref = read_dat(cfg.baseline.dat_file)
    print(f"\n{ref.name}: built baseline vs the file's own points ({ref.xy.shape[0]} points)")
    print(f"  {af.frame}; coordinates given to {af.decimals} decimals")
    print(f"  largest distance from a file point: {100 * af.split_error:.4f} % chord ({af.split_error * mm:.3f} mm)")
    print(f"  typical (rms) distance:             {100 * af.split_rms:.4f} % chord ({af.split_rms * mm:.3f} mm)")
    print("  (smoothed to the scatter in the published points; see README, 'Loaded .dat files')")
    print(f"  trailing-edge gap: {200 * af.thk.te_half_thickness:.3f} % chord ({2 * af.thk.te_half_thickness * mm:.3f} mm)")
    title, label = f"{ref.name}: built here vs file points", "published (file)"
print(f"  leading-edge radius in the fit: {100 * af.thk.le_radius_implied:.3f} % chord")

plots.baseline_check(res.base, ref.xy, title, outdir / f"{res.name}_check.png", label=label)
print(f"\nFiles written to {outdir}")
