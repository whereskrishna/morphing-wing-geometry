"""How far back can the rear spar go? Spar positioning study.

    python scripts/spar_study.py                         # NASA SC(2)-0714 (the default)
    python scripts/spar_study.py config/naca641412.yaml  # any case file, e.g. the NACA 64(1)-412

A rear spar further back makes a bigger wing box, but leaves less skin
behind it to bend, so the morphing TE reaches less deflection before the
skin is over-strained. For each rear spar station in spar_study.x_rear_spar
(config/default.yaml) it puts the hinge at the spar, so everything behind
the spar is free to morph, and finds the largest droop and reflex that pass
every check (the same search as deflection_range.py). Then it reports the
furthest-aft station that still reaches the needed range
(spar_study.need_droop_deg / need_reflex_deg: PLACEHOLDERS).

Writes <output dir>/spar_study/:
    <name>.csv   one row per station: wing box size, range, lift per degree, what stopped it
    <name>.png   range, rough lift per degree and rear spar height against the station
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from morphgeom import copy_config, load_baseline, load_config  # noqa: E402
from morphgeom import plots  # noqa: E402
from morphgeom.deflection import describe, find_range, lift_per_degree, wing_box  # noqa: E402
from morphgeom.export import write_summary_csv  # noqa: E402

cfg = load_config(sys.argv[1] if len(sys.argv) > 1 else None)   # optional case file
st = cfg.spar_study
# st.x_rear_spar = [0.60, 0.625, 0.65, 0.675, 0.70]   # e.g. a finer look around one region
cfg.morph.method = st.method

af = load_baseline(cfg)
c_mm = 1000 * cfg.structure.chord_m
v2 = st.method == "constant_length"
print(f"\nRear spar study: {af.name}, {'constant_length (v2)' if v2 else 'camber (v1)'}, n = {cfg.morph.n:g}, "
      f"front spar x_LE = {cfg.morph.x_LE:.3f}, hinge at the rear spar")
print(f"  chord {cfg.structure.chord_m:.3f} m, skin {cfg.structure.t_skin_mm:.3f} mm, "
      f"strain limit {cfg.checks.tol_strain_pct:.3f} %  (placeholders unless set)")
print(f"  needed: droop {st.need_droop_deg:+.3f} deg, reflex {st.need_reflex_deg:+.3f} deg  (PLACEHOLDERS)")
print(f"  {'rear spar':>9s} {'box length':>11s} {'spar height':>12s} {'box share':>10s} "
      f"{'largest droop':>14s} {'largest reflex':>15s} {'lift per deg':>13s}   reaches need   stopped by")
print(f"  {'x / c':>9s} {'(mm)':>11s} {'(mm)':>12s} {'(% area)':>10s} {'(deg)':>14s} {'(deg)':>15s} "
      f"{'(cl, est.)':>13s}")

rows = []
for xs in sorted(st.x_rear_spar):
    c = copy_config(cfg)
    c.morph.x_h, c.morph.x_rear_spar = xs, None
    rng = find_range(c, af, "TE")
    d, r = rng["droop"], rng["reflex"]
    box = wing_box(af, cfg.morph.x_LE, xs)
    margin = min(d.angle_deg - st.need_droop_deg, st.need_reflex_deg - r.angle_deg)
    row = {"x_rear_spar": xs, "box_length_mm": box["length"] * c_mm,
           "rear_spar_height_mm": box["rear_height"] * c_mm, "box_area_share_pct": 100 * box["area_share"],
           "droop_max_deg": d.angle_deg, "reflex_max_deg": r.angle_deg,
           "droop_tip_mm": d.delta * c_mm, "reflex_tip_mm": r.delta * c_mm,
           "lift_per_deg_estimate": lift_per_degree(d),
           "droop_stopped_by": describe(d, cfg.range.max_deg).replace(",", ";"),
           "reflex_stopped_by": describe(r, cfg.range.max_deg).replace(",", ";"),
           "margin_deg": margin, "reaches_need": margin >= 0}
    rows.append(row)
    stops = row["droop_stopped_by"] if row["droop_stopped_by"] == row["reflex_stopped_by"] \
        else f"{row['droop_stopped_by']} / {row['reflex_stopped_by']}"
    print(f"  {xs:9.3f} {row['box_length_mm']:11.3f} {row['rear_spar_height_mm']:12.3f} "
          f"{row['box_area_share_pct']:10.3f} {d.angle_deg:14.3f} {r.angle_deg:15.3f} "
          f"{row['lift_per_deg_estimate']:13.3f}   "
          f"{'yes' if margin >= 0 else 'NO ':<12s}   {stops}")

# Furthest-aft station that still reaches the need. Between the last station
# that does and the first behind it that does not, the margin (degrees to
# spare, the smaller of droop and reflex) is taken as straight.
answer = None
ok = [i for i, row in enumerate(rows) if row["reaches_need"]]
if not ok:
    print("\n  No tested station reaches the needed range: move the spar forward, or relax the limits.")
elif ok[-1] == len(rows) - 1:
    answer = rows[-1]["x_rear_spar"]
    print(f"\n  Every station up to x = {answer:.3f} reaches the needed range. Test further aft to find the limit.")
else:
    a, b = rows[ok[-1]], rows[ok[-1] + 1]
    answer = a["x_rear_spar"] + (b["x_rear_spar"] - a["x_rear_spar"]) * a["margin_deg"] / (a["margin_deg"] - b["margin_deg"])
    print(f"\n  Furthest-aft rear spar that still reaches the needed range: x = {a['x_rear_spar']:.3f} (tested);")
    print(f"  about x = {answer:.3f} between the tested stations. Behind that, the morphing TE runs out of skin.")
print("  These depend on the placeholder chord, skin thickness, strain limit and needed range (README).")
print("  Lift per degree: rough thin-airfoil estimate at small deflections. A spar further back")
print("  leaves a shorter morphing TE, so each degree of it gives less lift.")

name = f"{af.name}_spar_study_n{cfg.morph.n:g}" + ("_constlen" if v2 else "")
outdir = Path(cfg.output.dir) / "spar_study"
outdir.mkdir(parents=True, exist_ok=True)
write_summary_csv(outdir / f"{name}.csv", rows)
plots.spar_study(rows, cfg, answer, f"Rear spar study - {name}", outdir / f"{name}.png")
print(f"\nFiles written to {outdir}")
