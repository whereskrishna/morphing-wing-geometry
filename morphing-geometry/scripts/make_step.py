"""Turn CAD point files into STEP files that SolidWorks opens directly (File > Open).

    python scripts/make_step.py outputs/sc20714_tunnel/baseline/sc20714_baseline_cad_mm.txt --span 203.2
    python scripts/make_step.py outputs/sc20714/shapes/*_cad_mm.txt          # 2D faces only

Takes any <name>_cad_mm.txt written by make_baseline.py, make_shape.py or
sweep.py (X Y Z in mm) and writes next to it:
    <name>_2d.step                 the closed profile as a flat face at z = 0 (2D CFD, sketches)
    <name>_3d_span<S>mm.step       a solid: that face extruded S mm, centred on z = 0 (only with --span)
Then reads each file back and prints its size, area and volume, so you can
see it came out right before opening it in CAD.

The span is not a design number this project sets: pass it with --span (mm).
For the NASA SC(2)-0714 tunnel model it is 203.2 mm (README, "Tunnel scale").

Needs pythonocc-core (in environment.yml).
FIND DEFLECTION RANGE
hinge point
after last spar free to morph wherever
find a way to test how far back the last spar should go -> spar positioning


"""
import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from morphgeom.step import check_step, write_step  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
ap.add_argument("cad_files", nargs="+", help="one or more <name>_cad_mm.txt files")
ap.add_argument("--span", type=float, default=0.0, help="span of the 3D solid in mm (leave out for 2D only)")
args = ap.parse_args()

for path in map(Path, args.cad_files):
    xy = np.loadtxt(path)[:, :2]
    stem = path.stem.removesuffix("_cad_mm")
    x, y = xy[:, 0], xy[:, 1]
    area = 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))   # area inside the points
    jobs = [(path.with_name(f"{stem}_2d.step"), 0.0)]
    if args.span > 0:
        jobs.append((path.with_name(f"{stem}_3d_span{args.span:g}mm.step"), args.span))

    print(f"\n{path.name}: {xy.shape[0]} points, chord {x.max() - x.min():.3f} mm, "
          f"area inside the points {area:.3f} mm2")
    for out, span in jobs:
        write_step(out, xy, span_mm=span, name=stem)
        c = check_step(out)
        x0, y0, z0, x1, y1, z1 = (round(v, 3) + 0.0 for v in c["bbox_mm"])   # + 0.0: no "-0.000"
        print(f"  {out.name}")
        body = f"1 solid with {c['faces']} faces" if c["solids"] == 1 else f"{c['faces']} flat face, no solid"
        print(f"    geometry valid: {'yes' if c['valid'] else 'NO'};  {body}")
        print(f"    size: x {x0:.3f} to {x1:.3f}, y {y0:.3f} to {y1:.3f}, z {z0:.3f} to {z1:.3f} mm")
        if span:
            print(f"    volume {c['volume_mm3']:.3f} mm3 = {c['volume_mm3'] / span:.3f} mm2 x span "
                  f"({100 * (c['volume_mm3'] / span / area - 1):+.3f} % from the points)")
        else:
            print(f"    face area {c['face_areas_mm2'][0]:.3f} mm2 "
                  f"({100 * (c['face_areas_mm2'][0] / area - 1):+.3f} % from the points)")
print("\nOpen in SolidWorks with File > Open, file type STEP (*.step; *.stp).")
