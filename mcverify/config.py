"""Load and check the YAML configuration files."""
from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml

# Projected leaf widths at isocentre (cm), ordered from -Y to +Y (DICOM order).
LEAF_WIDTHS_CM = {
    "millennium120": [1.0] * 10 + [0.5] * 40 + [1.0] * 10,
    "hd120": [0.5] * 14 + [0.25] * 32 + [0.5] * 14,
}

_REQUIRED = {
    "": ["name", "sad_cm", "phase_space", "physics", "jaws", "mlc", "calibration"],
    "phase_space": ["chunk_dir"],
    "physics": ["modules", "production_cut_mm"],
    "jaws": ["upper", "lower"],
    "mlc": ["model", "travel_axis", "source_to_upstream_cm", "thickness_cm", "length_cm", "material"],
}


class ConfigError(ValueError):
    pass


def deep_update(base: dict, override: dict) -> dict:
    """Recursively merge ``override`` into a copy of ``base``."""
    out = deepcopy(base)
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_update(out[k], v)
        else:
            out[k] = deepcopy(v)
    return out


def load_yaml(path: str | Path) -> dict[str, Any]:
    with open(path) as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise ConfigError(f"{path}: top level must be a mapping")
    return data


def validate_machine(cfg: dict) -> dict:
    for section, keys in _REQUIRED.items():
        node = cfg if section == "" else cfg.get(section)
        if not isinstance(node, dict):
            raise ConfigError(f"missing section '{section}'")
        for k in keys:
            if k not in node:
                raise ConfigError(f"missing key '{section + '.' if section else ''}{k}'")
    model = cfg["mlc"]["model"]
    if model not in LEAF_WIDTHS_CM:
        raise ConfigError(f"unknown MLC model '{model}', choose from {sorted(LEAF_WIDTHS_CM)}")
    for name in ("upper", "lower"):
        ax = cfg["jaws"][name]["travel_axis"]
        if ax not in ("X", "Y"):
            raise ConfigError(f"jaws.{name}.travel_axis must be X or Y")
    if {cfg["jaws"]["upper"]["travel_axis"], cfg["jaws"]["lower"]["travel_axis"]} != {"X", "Y"}:
        raise ConfigError("one jaw pair must travel along X and the other along Y")
    if cfg["mlc"]["travel_axis"] not in ("X", "Y"):
        raise ConfigError("mlc.travel_axis must be X or Y")
    mat = cfg["mlc"]["material"]
    if isinstance(mat, dict):
        if len(mat["components"]) != len(mat["fractions"]):
            raise ConfigError("mlc.material components and fractions differ in length")
        if abs(sum(mat["fractions"]) - 1.0) > 1e-6:
            raise ConfigError("mlc.material fractions must sum to 1")
    return cfg


def load_machine(path: str | Path, overrides: dict | None = None) -> dict:
    cfg = load_yaml(path)
    if overrides:
        cfg = deep_update(cfg, overrides)
    return validate_machine(cfg)


def leaf_widths_cm(cfg: dict) -> list[float]:
    return list(LEAF_WIDTHS_CM[cfg["mlc"]["model"]])
