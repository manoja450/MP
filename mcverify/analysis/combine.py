"""Combine independent TOPAS jobs of one field into a mean dose per history and
its statistical uncertainty (batch method: spread of the per-job means)."""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .topas_io import read_scorer, to_1d


@dataclass
class Combined:
    positions_cm: np.ndarray | None
    dose_per_history: np.ndarray          # Gy per history
    rel_unc: np.ndarray                   # 1 sigma, relative (batch method)
    n_jobs: int
    total_histories: int
    meta: dict = field(default_factory=dict)


def job_dirs(field_dir: str | Path) -> list[Path]:
    return sorted(p.parent for p in Path(field_dir).glob("job*/meta.json"))


def combine_scorer(field_dir: str | Path, scorer: str, require_all: bool = False) -> Combined:
    """Combine ``scorer`` over every finished job in ``field_dir``."""
    per_job, hist, meta0 = [], [], None
    missing = []
    for jd in job_dirs(field_dir):
        meta = json.loads((jd / "meta.json").read_text())
        sc = meta["scorers"][scorer]
        f = jd / sc["file"]
        if not f.exists():
            missing.append(str(f))
            continue
        res = read_scorer(f)
        per_job.append(to_1d(res, sc["axis"]) / meta["histories"])
        hist.append(meta["histories"])
        meta0 = meta0 or {"scorer": sc, **{k: meta[k] for k in ("field", "nominal_size_cm", "machine")}}
    if missing and require_all:
        raise FileNotFoundError(f"{len(missing)} outputs missing, e.g. {missing[0]}")
    if not per_job:
        raise FileNotFoundError(f"no finished jobs for scorer {scorer} in {field_dir}")
    a = np.vstack(per_job)
    w = np.asarray(hist, float)
    mean = (a * w[:, None]).sum(0) / w.sum()
    n = a.shape[0]
    if n > 1:
        sem = a.std(axis=0, ddof=1) / np.sqrt(n)
        rel = np.divide(sem, np.abs(mean), out=np.full_like(mean, np.inf), where=mean != 0)
    else:
        rel = np.full_like(mean, np.nan)
    pos = np.asarray(meta0["scorer"]["positions_cm"]) if meta0["scorer"].get("positions_cm") else None
    if pos is not None and np.any(np.diff(pos) < 0):        # PDD bins are stored deep -> shallow
        order = np.argsort(pos)
        pos, mean, rel = pos[order], mean[order], rel[order]
    return Combined(pos, mean, rel, n, int(w.sum()), meta0)
