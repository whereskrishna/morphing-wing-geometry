"""Every flap angle x hinge x shape exponent: which TE shapes pass, and how much lift.

    python scripts/flap_study.py                         # NASA SC(2)-0714 (the default)
    python scripts/flap_study.py config/naca641412.yaml  # any case file, e.g. the NACA 64(1)-412

For every hinge station and shape exponent n listed under flap_study
(config/default.yaml), with the rear spar at the hinge, it makes the morphed
TE at every equivalent flap angle from angle_min_deg to angle_max_deg and
runs every check. It also finds the exact largest droop and reflex (the same
search as deflection_range.py). Then it turns the rough lift estimate into a
force at the cruise condition under flight: (config/default.yaml).

The shapes are made on all processor cores at once: about 1.5 minutes on 15 cores.
Writes <output dir>/flap_study/:
    <name>_grid.csv     one row per shape tried: passes or what failed, peak strain, lift
    <name>_limits.csv   one row per hinge and n: largest droop and reflex, lift per degree
    <name>_limits.png   largest droop, largest reflex, total range and lift, as hinge x n maps
    <name>_map.png      every shape tried, one panel per n: passes, or what stopped it
"""
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from morphgeom import load_config  # noqa: E402
from morphgeom import plots  # noqa: E402
from morphgeom.aero import flight_condition, lift_per_span  # noqa: E402
from morphgeom.deflection import flap_point  # noqa: E402
from morphgeom.export import write_summary_csv  # noqa: E402


def table(title, limits, xs, ns, value, fmt):
    """Print one value as a table: a row per n, a column per hinge."""
    by = {(r["x_h"], r["n"]): r for r in limits}
    print(f"\n  {title}")
    print("    hinge x/c " + "".join(f"{x:7.3f}" for x in xs))
    for n in ns:
        print(f"    n = {n:<6g}" + "".join(fmt(value(by[(x, n)])) for x in xs))


if __name__ == "__main__":       # needed for the worker processes
    cfg = load_config(sys.argv[1] if len(sys.argv) > 1 else None)   # optional case file
    fs = cfg.flap_study
    # fs.x_h = [0.60, 0.65, 0.70]; fs.n = [3, 4]                     # e.g. a smaller study
    cfg.morph.method = fs.method
    v2 = fs.method == "constant_length"
    st, need = cfg.structure, cfg.spar_study

    xs = sorted(float(x) for x in fs.x_h)
    ns = sorted(float(n) for n in fs.n)
    k = np.arange(np.ceil(fs.angle_min_deg / fs.angle_step_deg), np.floor(fs.angle_max_deg / fs.angle_step_deg) + 1)
    angles = [float(a) for a in k * fs.angle_step_deg if abs(a) > 1e-9]

    fc = flight_condition(cfg.flight)
    fl = cfg.flight
    print(f"\nFlap study: {cfg.baseline.name}, {'constant_length (v2)' if v2 else 'camber (v1)'}, TE, "
          f"rear spar at the hinge, front spar {cfg.morph.x_LE:.3f}")
    print(f"  {len(xs)} hinges x {len(ns)} values of n x {len(angles)} angles "
          f"({angles[0]:+g} to {angles[-1]:+g} deg) = {len(xs) * len(ns) * len(angles)} shapes, plus the exact limits")
    print(f"  chord {st.chord_m:.3f} m, skin {st.t_skin_mm:.3f} mm, strain limit {cfg.checks.tol_strain_pct:.3f} %"
          f"  (placeholders unless set)")
    print(f"\n  Cruise: {fl.altitude_ft:g} ft, {fl.velocity_m_s:g} m/s, Mach {fc['mach']:.3f}, "
          f"density {fl.density_kg_m3:g} kg/m3, pressure {fl.pressure_pa:g} Pa")
    print(f"    dynamic pressure q = 1/2 rho V^2 = {fc['q_pa']:.0f} Pa;  Reynolds number "
          f"{fc['reynolds_per_m'] / 1e6:.2f} million per metre of chord")
    gap = fl.density_kg_m3 / fc["density_from_p_t"] - 1
    if abs(gap) > 0.02:
        print(f"    CHECK: the pressure and speed of sound give a density of {fc['density_from_p_t']:.4f} kg/m3 "
              f"(temperature {fc['temperature_k']:.1f} K),")
        print(f"    {100 * abs(gap):.0f}% {'lower' if gap > 0 else 'higher'} than the "
              f"{fl.density_kg_m3:g} given. Lift forces below use the given density; with "
              f"{fc['density_from_p_t']:.4f}")
        print(f"    they would change by {100 * (fc['density_from_p_t'] / fl.density_kg_m3 - 1):+.0f}%.")

    workers = fs.workers or os.cpu_count()
    t0 = time.time()
    grid, limits = [], []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        jobs = [pool.submit(flap_point, cfg, x, n, angles) for n in ns for x in xs]
        for i, job in enumerate(as_completed(jobs), 1):
            rows, lim = job.result()
            grid += rows
            limits.append(lim)
            print(f"\r  done {i} of {len(jobs)} hinge / n pairs ({time.time() - t0:.0f} s)", end="", flush=True)
    print()
    grid.sort(key=lambda r: (r["n"], r["x_h"], r["angle_deg"]))
    limits.sort(key=lambda r: (r["n"], r["x_h"]))

    # lift as a force at cruise: per metre of span, with the chord above
    for r in grid:
        r["lift_change_n_per_m"] = lift_per_span(r["dcl_estimate"], fc["q_pa"], st.chord_m)
    for r in limits:
        r["lift_per_deg_n_per_m"] = lift_per_span(r["lift_per_deg_estimate"], fc["q_pa"], st.chord_m)
        r["lift_per_deg_n_per_m_compressible"] = r["lift_per_deg_n_per_m"] * fc["prandtl_glauert"]
        r["reaches_need"] = r["droop_max_deg"] >= need.need_droop_deg and r["reflex_max_deg"] <= need.need_reflex_deg

    table("Largest droop (deg)", limits, xs, ns, lambda r: r["droop_max_deg"], lambda v: f"{v:7.1f}")
    table("Largest reflex (deg)", limits, xs, ns, lambda r: r["reflex_max_deg"], lambda v: f"{v:7.1f}")
    table(f"Reaches the needed {need.need_droop_deg:+g} / {need.need_reflex_deg:+g} deg (PLACEHOLDERS)?",
          limits, xs, ns, lambda r: r["reaches_need"], lambda v: f"{'yes' if v else '-':>7s}")
    table("Lift change per degree, low-speed thin-airfoil estimate (cl per deg)", limits, xs, ns,
          lambda r: r["lift_per_deg_estimate"], lambda v: f"{v:7.3f}")
    table(f"The same at cruise: newtons of lift per metre of span per degree (q = {fc['q_pa']:.0f} Pa, "
          f"chord {st.chord_m:g} m)", limits, xs, ns, lambda r: r["lift_per_deg_n_per_m"], lambda v: f"{v:7.0f}")

    ok = [r for r in limits if r["reaches_need"]]
    print(f"\n  {len(ok)} of {len(limits)} hinge / n pairs reach the needed range (placeholder need).")
    if ok:
        far = max(ok, key=lambda r: (r["x_h"], r["range_deg"]))
        best = max(ok, key=lambda r: r["range_deg"])
        print(f"    furthest-aft hinge that does: x/c = {far['x_h']:.3f} with n = {far['n']:g} "
              f"({far['droop_max_deg']:+.1f} / {far['reflex_max_deg']:+.1f} deg)")
        print(f"    widest range: x/c = {best['x_h']:.3f}, n = {best['n']:g} "
              f"({best['droop_max_deg']:+.1f} / {best['reflex_max_deg']:+.1f} deg)")
    stops = {}
    for r in limits:
        for s in (r["droop_stopped_by"], r["reflex_stopped_by"]):
            stops[s] = stops.get(s, 0) + 1
    print("    what stopped the searches: " + ", ".join(f"{s} {c}" for s, c in sorted(stops.items(), key=lambda t: -t[1])))
    print("\n  Lift: thin-airfoil theory, good for the first few degrees only (at the limits the flow")
    print(f"  would separate). It leaves out compressibility: the simple correction would multiply it by "
          f"{fc['prandtl_glauert']:.2f} at Mach {fc['mach']:.2f},")
    print("  but that correction is unreliable this close to the speed where shocks form. CFD gives the real numbers.")
    print("  Ranges depend on the placeholder chord, skin, strain limit and need (README open questions).")

    name = f"{cfg.baseline.name}_flap_study" + ("_constlen" if v2 else "")
    outdir = Path(cfg.output.dir) / "flap_study"
    outdir.mkdir(parents=True, exist_ok=True)
    write_summary_csv(outdir / f"{name}_grid.csv", grid)
    write_summary_csv(outdir / f"{name}_limits.csv", limits)
    plots.flap_study_limits(limits, cfg, fc, f"Flap study - {name}", outdir / f"{name}_limits.png")
    plots.flap_study_map(grid, limits, cfg, f"Every shape tried - {name}", outdir / f"{name}_map.png")
    print(f"\nFiles written to {outdir} ({time.time() - t0:.0f} s)")
