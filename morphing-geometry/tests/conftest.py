import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from morphgeom import load_baseline, load_config  # noqa: E402


@pytest.fixture(scope="session")
def root():
    return ROOT


@pytest.fixture(scope="session")
def cfg():
    """The NACA 64(1)-412 case: most tests check against its published tables."""
    return load_config("config/naca641412.yaml")


@pytest.fixture(scope="session")
def af(cfg):
    return load_baseline(cfg)


@pytest.fixture(scope="session")
def xref():
    """The 26 standard NACA table stations."""
    return np.concatenate([[0, 0.005, 0.0075, 0.0125, 0.025, 0.05, 0.075, 0.1], np.linspace(0.15, 1.0, 18)])
