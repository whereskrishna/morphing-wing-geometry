"""Morphing wing shape generator (geometry team).

Takes a baseline airfoil plus a few morph settings, writes morphed airfoil
coordinate files, and checks each shape for buildability. See README.md.

    from morphgeom import load_config, make_shape, write_outputs
    cfg = load_config()
    res = make_shape(cfg)
"""
from .config import copy_config, load_config
from .pipeline import Result, load_baseline, make_shape, write_outputs

__all__ = ["load_config", "copy_config", "load_baseline", "make_shape", "write_outputs", "Result"]
