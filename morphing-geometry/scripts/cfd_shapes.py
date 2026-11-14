"""The shapes to send to CFD: coordinate files, a list of cases, a picture and a note for the CFD team.

    python scripts/cfd_shapes.py                         # NASA SC(2)-0714 (the default)
    python scripts/cfd_shapes.py config/naca641412.yaml  # any case file, e.g. the NACA 64(1)-412

Which shapes: cfd_set in config/default.yaml. Case 1 is the baseline; then each
group, every hinge x shape exponent n x flap angle (TE only, rear spar at the
hinge). Shapes that fail a check are left out and listed.

Writes <output dir>/cfd_shapes/ (old .dat files there are removed first, so the
folder only ever holds the current set):
    <name>.dat        one per shape: chord 1, Selig order, 6 decimals. Send these.
    cfd_cases.csv     every case: group, hinge, n, flap angle, tip drop, rough lift estimate
    cfd_shapes.png    every shape drawn against the baseline
    README_cfd.txt    for the CFD team: how the shapes are set up, the flight condition,
                      what is still open and what to send back
"""
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from morphgeom import load_config  # noqa: E402
from morphgeom.cfd import case_rows, cfd_cases, write_cfd_set  # noqa: E402

cfg = load_config(sys.argv[1] if len(sys.argv) > 1 else None)   # optional case file
t0 = time.time()
cases, skipped = cfd_cases(cfg)
outdir = Path(cfg.output.dir) / "cfd_shapes"
files = write_cfd_set(cases, skipped, cfg, outdir)

print(f"\nShapes for CFD: {cfg.baseline.name}, {cfg.cfd_set.method}, TE only, rear spar at the hinge")
print(f"  {'case':>4s}  {'group':<15s} {'flap (deg)':>10s} {'hinge':>6s} {'n':>4s} {'tip drop (% c)':>15s} "
      f"{'strain (%)':>11s} {'lift est. (cl)':>15s}   file")
for r in case_rows(cases, cfg):
    hinge = r["hinge_x_c"] if isinstance(r["hinge_x_c"], str) else f"{r['hinge_x_c']:.3f}"
    n = r["n"] if isinstance(r["n"], str) else f"{r['n']:g}"
    angle = f"{r['flap_angle_deg']:+g}" if r["case"] > 1 else "0"
    print(f"  {r['case']:>4d}  {r['group']:<15s} {angle:>10s} {hinge:>6s} {n:>4s} "
          f"{r['tip_drop_pct_chord']:>15.3f} {r['strain_max_pct']:>11.3f} {r['lift_change_cl_estimate']:>15.3f}   {r['file']}")
if skipped:
    print("\n  Left out (fail a check with the current placeholders):")
    for s in skipped:
        print(f"    {s['group']}: {s['angle_deg']:+g} deg, hinge {s['x_h']:.3f}, n {s['n']:g} ({', '.join(s['failed'])})")
print(f"\n{len(cases)} shapes written to {outdir} ({time.time() - t0:.0f} s)")
print("Send the folder to the CFD team: README_cfd.txt explains it.")
