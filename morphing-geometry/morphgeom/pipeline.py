"""Top level: load the baseline, make one shape with all its checks, write its files."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import export, plots
from .checks import run_all
from .morph import Morph, apply, equivalent_angles
from .naca import naca6_section
from .surfaces import Airfoil, Surfaces, build, loop, split_dat, stations


def load_baseline(cfg) -> Airfoil:
    """The baseline airfoil named in the config, as camber + thickness.

    cfg.baseline.source = 'naca6': thickness table + a = 1.0 mean line
                          'dat':   any .dat file, split approximately
    """
    src = cfg.baseline.source
    if src == "naca6":
        return naca6_section(cfg.baseline.thickness_file, cfg.baseline.cl_design, cfg.baseline.name)
    if src == "dat":
        return split_dat(cfg.baseline.dat_file, cfg.baseline.name)
    raise ValueError(f'Unknown baseline source "{src}". Use naca6 or dat.')


@dataclass
class Result:
    name: str             # file name stem, e.g. naca641412_TE+0.020_xh0.650_n3
    af: Airfoil
    morph: object         # the morph settings used
    mr: Morph             # the morph (morph.apply)
    base: Surfaces        # baseline surfaces at the output points
    shape: Surfaces       # morphed surfaces at the output points
    angles: dict          # equivalent deflection angles (deg)
    checks: dict          # every check value and pass/fail (checks.run_all)
    row: dict = None      # one flat summary row (for printing and summary.csv)


def make_shape(cfg, af: Airfoil = None) -> Result:
    """Baseline + one morphed shape + every check, from the settings.

    Pass af to reuse an already loaded baseline (faster in sweeps).
    cfg.morph.method picks the morph: 'camber' (v1) or 'constant_length' (v2).
    """
    if af is None:
        af = load_baseline(cfg)
    mp = cfg.morph
    x = stations(cfg.points.n_per_surface, cfg.points.spacing)
    mr = apply(af, mp, cfg)
    res = Result(name=export.shape_name(af.name, mp), af=af, morph=mp, mr=mr,
                 base=build(af.cam, af.thk, x), shape=mr.eval(x),
                 angles=equivalent_angles(mp), checks=run_all(af, mr, cfg))
    res.row = export.summary_row(res)
    return res


def write_outputs(res: Result, cfg, outdir, make_plots: bool = None) -> dict:
    """Write the coordinate files (and plots) for one shape.

    <name>.dat          Selig-style coordinates, chord = 1, 6 decimals (XFOIL, CFD meshing)
    <name>_cad_mm.txt   X Y Z in mm, Z = 0, one ordered loop (SolidWorks "Curve Through XYZ Points")
    <name>_*.png        overlay, camber lines, curvature (if make_plots)
    """
    if make_plots is None:
        make_plots = cfg.output.plots
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    xy = loop(res.shape, cfg.points.order)
    files = {"dat": outdir / f"{res.name}.dat"}
    export.write_selig(files["dat"], xy, res.name)
    if cfg.output.cad:
        files["cad"] = outdir / f"{res.name}_cad_mm.txt"
        export.write_cad_xyz(files["cad"], xy, cfg.structure.chord_m)
    if make_plots:
        files["overlay"] = outdir / f"{res.name}_overlay.png"
        files["camber"] = outdir / f"{res.name}_camber.png"
        files["curvature"] = outdir / f"{res.name}_curvature.png"
        plots.overlay(res, files["overlay"])
        plots.camber_lines(res, files["camber"])
        plots.curvature(res, cfg, files["curvature"])
    return files
