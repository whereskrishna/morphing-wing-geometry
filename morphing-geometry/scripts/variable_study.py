"""Which settings matter? Change ONE setting at a time and find the range each time.

    python scripts/variable_study.py                         # NASA SC(2)-0714 (the default)
    python scripts/variable_study.py config/naca641412.yaml  # any case file, e.g. the NACA 64(1)-412

For each setting listed under variable_study (config/default.yaml), it tries
each value with everything else left at its default, and finds the largest
droop and reflex that pass every check (the same search as
deflection_range.py). Then it ranks the settings by how much they move the
range. The values are test values, chosen wide enough to show which
settings matter. They are not design choices.

Takes a few minutes (about 5 s per value). Writes <output dir>/variable_study/:
    <name>.csv   one row per value tried: largest droop and reflex, tip drop, what stopped it
    <name>.png   the ranking at the top, then one small plot per setting
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from morphgeom import load_baseline, load_config  # noqa: E402
from morphgeom import plots  # noqa: E402
from morphgeom.deflection import bend_radius_limit, describe, find_range, spread, unsolved, vary_one  # noqa: E402
from morphgeom.export import write_summary_csv  # noqa: E402

cfg = load_config(sys.argv[1] if len(sys.argv) > 1 else None)   # optional case file
vs = cfg.variable_study
# vs.edges = ["TE"]                                              # e.g. the TE only
# vs.TE = [{"setting": "morph.x_h", "values": [0.60, 0.62, 0.64, 0.66], "label": "hinge x/c"}]
cfg.morph.method = vs.method

af = load_baseline(cfg)
v2 = vs.method == "constant_length"
outdir = Path(cfg.output.dir) / "variable_study"
outdir.mkdir(parents=True, exist_ok=True)
mp, st = cfg.morph, cfg.structure
print(f"\nVariable study: {af.name}, {'constant_length (v2)' if v2 else 'camber (v1)'}")
print(f"  defaults: front spar {mp.x_LE:.3f}, hinge {mp.x_h:.3f}, n = {mp.n:g}, chord {st.chord_m:.3f} m, "
      f"skin {st.t_skin_mm:.3f} mm, strain limit {cfg.checks.tol_strain_pct:.3f} %  (placeholders)")

cache = {}
for edge in vs.edges:
    entries = getattr(vs, edge)
    n_runs = sum(len(e["values"]) for e in entries)
    print(f"\n=== {edge}: {len(entries)} settings, {n_runs} values ===")
    base = find_range(cfg, af, edge)
    cache[(edge, ())] = base
    default = {"droop": base["droop"].angle_deg, "reflex": base["reflex"].angle_deg}
    need = (cfg.spar_study.need_droop_deg, cfg.spar_study.need_reflex_deg) if edge == "TE" else None
    print(f"  at the defaults: droop {default['droop']:+.3f} deg, reflex {default['reflex']:+.3f} deg "
          f"(stopped by {describe(base['droop'], cfg.range.max_deg)} / {describe(base['reflex'], cfg.range.max_deg)})")
    if need:
        print(f"  needed: droop {need[0]:+.3f} deg, reflex {need[1]:+.3f} deg  (PLACEHOLDERS)")

    groups, all_rows = [], []
    for e in entries:
        rows = vary_one(cfg, af, edge, e, cache)
        label = rows[0]["label"]
        groups.append((label, rows))
        all_rows += rows
        print(f"\n  {label}")
        print(f"    {'value':>8s} {'largest droop':>14s} {'largest reflex':>15s} {'total range':>12s} "
              f"{'lift per deg':>13s}" + ("   reaches need" if need else "") + "   stopped by")
        for r in rows:
            reach = r["droop_max_deg"] >= need[0] and r["reflex_max_deg"] <= need[1] if need else None
            stops = r["droop_stopped_by"] if r["droop_stopped_by"] == r["reflex_stopped_by"] \
                else f"{r['droop_stopped_by']} / {r['reflex_stopped_by']}"
            print(f"    {r['value']:>8s} {r['droop_max_deg']:14.3f} {r['reflex_max_deg']:15.3f} "
                  f"{r['range_deg']:12.3f} {r['lift_per_deg_estimate']:13.3f}"
                  + (f"   {'yes' if reach else 'NO ':<12s}" if need else "")
                  + f"   {stops}" + ("   (default)" if r["is_default"] == "yes" else "")
                  + ("   (does not solve at any angle: left out of the ranking)" if unsolved(r) else ""))

    # ranking: how much each setting moves the total range (droop - reflex)
    print(f"\n  Which settings move the {edge} range most (over the values tried):")
    for label, rows in sorted(groups, key=lambda g: -spread(g[1])):
        ok = [r for r in rows if not unsolved(r)]
        lo, hi = min(ok, key=lambda r: r["range_deg"]), max(ok, key=lambda r: r["range_deg"])
        note = "no effect" if spread(rows) < 0.05 else \
            f"total range {lo['range_deg']:.1f} deg at {lo['value']} to {hi['range_deg']:.1f} deg at {hi['value']}"
        print(f"    {label:<38s} spread {spread(rows):7.3f} deg   {note}")

    # chord, skin thickness and strain limit: the strain check only sees one number
    trio = [r for r in all_rows if r["setting"] in ("structure.t_skin_mm", "structure.chord_m", "checks.tol_strain_pct")]
    if edge == "TE" and trio:
        same = {}
        for r in trio:
            same.setdefault(round(r["bend_radius_limit_c"], 6), []).append(r)
        worst = max((max(abs(a["droop_max_deg"] - b["droop_max_deg"]) for a in g for b in g) for g in same.values()
                     if len(g) > 1), default=None)
        print(f"\n  Chord, skin thickness and strain limit act through ONE number: skin thickness / (2 x strain limit)")
        print(f"  / chord = {bend_radius_limit(cfg):.3f} at the defaults (the tightest extra bend allowed, in chords).")
        if worst is not None:
            print(f"  Runs with the same number give the same range (largest difference {worst:.3f} deg).")
        print("  So the team only has to agree that one number to fix the strain limit on the range.")

    name = f"{af.name}_{edge}_variable_study" + ("_constlen" if v2 else "")
    write_summary_csv(outdir / f"{name}.csv", all_rows)
    plots.variable_study(groups, cfg, edge, default, f"{edge} variable study - {name}", outdir / f"{name}.png",
                         need=need)

print("\n  These depend on the placeholder settings (README open questions). Each setting was changed")
print("  on its own, so effects of changing two at once are not shown; rerun with a case file for that.")
print(f"\nFiles written to {outdir}")
