"""Read the settings into dot-accessible settings.

    cfg = load_config()                        # config/default.yaml
    cfg = load_config("config/sc20714.yaml")   # a case file: default.yaml + its changes
    cfg.morph.delta_TE = 0.03                  # change anything in a script

A case file only lists what differs from default.yaml; everything else is
taken from default.yaml. File paths in the YAML are relative to the project
folder.
"""
from __future__ import annotations

import copy
from pathlib import Path
from types import SimpleNamespace

import yaml

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / "config" / "default.yaml"


def _to_ns(d):
    if isinstance(d, dict):
        return SimpleNamespace(**{k: _to_ns(v) for k, v in d.items()})
    return d


def _merge(base: dict, changes: dict) -> dict:
    out = dict(base)
    for k, v in (changes or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _merge(out[k], v)
        else:
            if k not in out:
                raise KeyError(f'Unknown setting "{k}" in the case file (check the spelling against default.yaml).')
            out[k] = v
    return out


def load_config(path=None) -> SimpleNamespace:
    with open(DEFAULT) as f:
        d = yaml.safe_load(f)
    if path is not None:
        path = Path(path)
        if not path.is_absolute() and not path.exists():
            path = ROOT / path
        if path.resolve() != DEFAULT.resolve():
            with open(path) as f:
                d = _merge(d, yaml.safe_load(f))
    cfg = _to_ns(d)
    # make file paths absolute (they are written relative to the project folder)
    cfg.baseline.thickness_file = str(ROOT / cfg.baseline.thickness_file)
    if cfg.baseline.dat_file:
        cfg.baseline.dat_file = str(ROOT / cfg.baseline.dat_file)
    cfg.output.dir = str(ROOT / cfg.output.dir)
    return cfg


def copy_config(cfg) -> SimpleNamespace:
    """An independent copy, for changing settings shape by shape (sweeps)."""
    return copy.deepcopy(cfg)
