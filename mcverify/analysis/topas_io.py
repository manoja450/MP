"""Read TOPAS binned-scorer output (CSV and binary).

CSV layout written by OpenTOPAS (TsVBinnedScorer):
    # TOPAS Version: ...
    # X in 80 bins of 0.5 cm
    # Y in 1 bin  of 1 cm
    # Z in 300 bins of 0.1 cm
    # DoseToMedium ( Gy ) : Sum   Standard_Deviation
    i, j, k, sum, std          <- one row per voxel
Binary (.bin): float64 values written with k outermost and i innermost, each
voxel holding all reported columns in order. The header (.binheader) has the
same comment lines as the CSV.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

_BIN_RE = re.compile(r"^#\s*(\S+)\s+in\s+(\d+)\s+bins?\s+of\s+([-+0-9.eE]+)\s*(\S*)")
_QTY_RE = re.compile(r"^#\s*(\S+)\s*(?:\(\s*([^)]*?)\s*\))?\s*:\s*(.+)$")


@dataclass
class ScorerResult:
    quantity: str
    unit: str
    columns: list[str]
    shape: tuple[int, int, int]           # (Ni, Nj, Nk)
    bin_widths: list[float]
    data: dict[str, np.ndarray]           # column -> array of shape (Ni, Nj, Nk)

    def column(self, name: str) -> np.ndarray:
        return self.data[name]


def _parse_header(lines: list[str]) -> tuple[list, str, str, list[str]]:
    dims, quantity, unit, cols = [], "", "", []
    for line in lines:
        m = _BIN_RE.match(line)
        if m:
            dims.append((m.group(1), int(m.group(2)), float(m.group(3)), m.group(4)))
            continue
        m = _QTY_RE.match(line)
        if m and not line.startswith("# TOPAS") and "Results for scorer" not in line \
                and "Parameter File" not in line and "Scored in" not in line and "Binary file" not in line:
            quantity, unit = m.group(1), (m.group(2) or "")
            cols = m.group(3).split()
    return dims, quantity, unit, cols


def read_csv(path: str | Path) -> ScorerResult:
    path = Path(path)
    header, rows = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                header.append(line.rstrip("\n"))
            elif line.strip():
                rows.append(line)
    dims, quantity, unit, cols = _parse_header(header)
    if not cols:
        raise ValueError(f"{path}: could not find the quantity/columns header line")
    arr = np.loadtxt(rows, delimiter=",", ndmin=2) if rows else np.zeros((0, len(cols)))
    if dims:
        shape = tuple(d[1] for d in dims)
        while len(shape) < 3:
            shape = shape + (1,)
        idx = arr[:, :3].astype(int)
        vals = arr[:, 3:]
    else:  # a single, undivided component
        shape = (1, 1, 1)
        idx = np.zeros((arr.shape[0], 3), dtype=int)
        vals = arr
    if vals.shape[1] != len(cols):
        raise ValueError(f"{path}: {vals.shape[1]} value columns but header lists {cols}")
    data = {}
    for c, name in enumerate(cols):
        a = np.zeros(shape)
        a[idx[:, 0], idx[:, 1], idx[:, 2]] = vals[:, c]   # missing (sparsified) rows stay 0
        data[name] = a
    return ScorerResult(quantity, unit, cols, shape, [d[2] for d in dims], data)


def read_binary(header_path: str | Path) -> ScorerResult:
    header_path = Path(header_path)
    header = [l.rstrip("\n") for l in header_path.read_text().splitlines() if l.startswith("#")]
    dims, quantity, unit, cols = _parse_header(header)
    bin_path = header_path.with_suffix(".bin")
    raw = np.fromfile(bin_path, dtype="<f8")
    shape = tuple(d[1] for d in dims) or (1, 1, 1)
    ni, nj, nk = shape
    if raw.size != ni * nj * nk * len(cols):
        raise ValueError(f"{bin_path}: size {raw.size} does not match {shape} x {len(cols)} columns")
    arr = raw.reshape(nk, nj, ni, len(cols))           # k outermost, i innermost
    data = {name: np.transpose(arr[..., c], (2, 1, 0)).copy() for c, name in enumerate(cols)}
    return ScorerResult(quantity, unit, cols, shape, [d[2] for d in dims], data)


def read_scorer(path: str | Path) -> ScorerResult:
    path = Path(path)
    if path.suffix == ".csv":
        return read_csv(path)
    if path.suffix in (".binheader", ".bin"):
        return read_binary(path.with_suffix(".binheader"))
    raise ValueError(f"unsupported TOPAS output: {path}")


def to_1d(res: ScorerResult, axis: str | None, column: str = "Sum") -> np.ndarray:
    """Collapse a scorer that is binned along a single axis ('i','j','k')."""
    a = res.data[column]
    if axis is None:
        return np.array([a.sum()])
    keep = "ijk".index(axis)
    other = tuple(x for x in range(3) if x != keep)
    if any(a.shape[o] != 1 for o in other):
        raise ValueError(f"scorer has shape {a.shape}, expected it binned only along {axis}")
    return a.reshape(-1).copy() if a.size == a.shape[keep] else a.squeeze()
