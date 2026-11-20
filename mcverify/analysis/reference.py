"""Load measured / published reference beam data.

Directory layout (see reference_data/README.md)::

    reference_data/<beam>/pdd_<field>.csv                  depth_cm,dose
    reference_data/<beam>/profile_x_<field>_d<depth>.csv   position_cm,dose
    reference_data/<beam>/profile_y_<field>_d<depth>.csv   position_cm,dose
    reference_data/<beam>/output_factors.csv               field,output_factor

``<field>`` is the validation field name (e.g. ``jaw10cm``), ``<depth>`` is
written like ``1p5cm`` / ``10cm``. Doses may be in any unit: they are
re-normalised by the analysis.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


def _read_two_columns(path: Path) -> tuple[np.ndarray, np.ndarray]:
    a = np.genfromtxt(path, delimiter=",", names=True)
    names = a.dtype.names
    x, y = np.asarray(a[names[0]], float), np.asarray(a[names[1]], float)
    order = np.argsort(x)
    return x[order], y[order]


def depth_tag(depth_cm: float) -> str:
    return f"{depth_cm:g}".replace(".", "p") + "cm"


class ReferenceData:
    def __init__(self, root: str | Path):
        self.root = Path(root)

    def pdd(self, field: str):
        p = self.root / f"pdd_{field}.csv"
        return _read_two_columns(p) if p.exists() else None

    def profile(self, axis: str, field: str, depth_cm: float):
        p = self.root / f"profile_{axis.lower()}_{field}_d{depth_tag(depth_cm)}.csv"
        return _read_two_columns(p) if p.exists() else None

    def output_factors(self) -> dict[str, float]:
        p = self.root / "output_factors.csv"
        if not p.exists():
            return {}
        out = {}
        for line in p.read_text().splitlines()[1:]:
            if line.strip() and not line.startswith("#"):
                name, val = [s.strip() for s in line.split(",")[:2]]
                out[name] = float(val)
        return out
