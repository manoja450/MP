"""Synthetic inputs for the tests: a tiny IAEA phase space and fake TOPAS outputs."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.special import erf

HEADER_TEMPLATE = """$IAEA_INDEX:
   1000  // test file
$TITLE:
Synthetic TrueBeam-like phase space for unit tests
$FILE_TYPE:
   0
$CHECKSUM:
   {checksum}
$RECORD_CONTENTS:
    1     // X is stored ?
    1     // Y is stored ?
    0     // Z is stored ?
    1     // U is stored ?
    1     // V is stored ?
    1     // W is stored ?
    1     // Weight is stored ?
    0     // Extra floats stored ?
    0     // Extra longs stored ?
$RECORD_CONSTANT:
   26.7000   // Constant Z
$RECORD_LENGTH:
   {rl}
$BYTE_ORDER:
   1234
$ORIG_HISTORIES:
   {hist}
$PARTICLES:
   {n}
$PHOTONS:
   {g}
$ELECTRONS:
   {e}
$POSITRONS:
   {p}
"""


def write_fake_phsp(stem: Path, n: int = 1000, seed: int = 0, new_history_flags: bool = False) -> dict:
    """Write stem.IAEAheader + stem.IAEAphsp with n records (type, E, x, y, u, v, weight)."""
    rng = np.random.default_rng(seed)
    rl = 1 + 4 * 6     # type byte + E, x, y, u, v, weight (w is derived from u, v)
    rec = np.zeros((n, rl), dtype=np.uint8)
    types = rng.choice([1, 2, 3], size=n, p=[0.9, 0.08, 0.02]).astype(np.int8)
    types = -types  # all moving with negative cosZ sign? irrelevant for the tests
    rec[:, 0] = types.view(np.uint8)
    energy = rng.uniform(0.1, 6.0, n).astype("<f4")
    if new_history_flags:
        energy[::3] *= -1
    floats = np.column_stack([energy, rng.normal(0, 3, n), rng.normal(0, 3, n),
                              rng.normal(0, 0.05, n), rng.normal(0, 0.05, n), np.ones(n)]).astype("<f4")
    rec[:, 1:] = floats.view(np.uint8).reshape(n, 24)
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    (Path(str(stem) + ".IAEAphsp")).write_bytes(rec.tobytes())
    g, e, p = (int(np.sum(np.abs(types) == k)) for k in (1, 2, 3))
    (Path(str(stem) + ".IAEAheader")).write_text(
        HEADER_TEMPLATE.format(checksum=n * rl, rl=rl, hist=5 * n, n=n, g=g, e=e, p=p))
    return {"n": n, "photons": g, "electrons": e, "positrons": p, "rl": rl,
            "flags": int(np.sum(energy < 0))}


# ---------------------------------------------------------------------------
def write_topas_csv(path: Path, shape: tuple[int, int, int], bin_cm: tuple[float, float, float],
                    values: np.ndarray, std: np.ndarray | None = None) -> None:
    """Write a TOPAS-style binned CSV (Sum and Standard_Deviation columns)."""
    values = np.asarray(values, float).reshape(shape)
    std = np.zeros(shape) if std is None else np.asarray(std, float).reshape(shape)
    lines = ["# TOPAS Version: 4.x (synthetic)", "# Parameter File: run.txt",
             "# Results for scorer: test", "# Scored in component: Box"]
    for name, n, w in zip("XYZ", shape, bin_cm):
        lines.append(f"# {name} in {n} {'bin ' if n == 1 else 'bins'} of {w} cm")
    lines.append("# DoseToMedium ( Gy ) : Sum   Standard_Deviation   ")
    ni, nj, nk = shape
    for k in range(nk):
        for j in range(nj):
            for i in range(ni):
                lines.append(f"{i}, {j}, {k}, {float(values[i, j, k])!r}, {float(std[i, j, k])!r}")
    Path(path).write_text("\n".join(lines) + "\n")


def analytic_profile(x, half_field, sigma=0.3, horns=0.0):
    """Blurred top-hat with optional 'horns' (dose per history units)."""
    base = 0.5 * (erf((x + half_field) / (np.sqrt(2) * sigma)) - erf((x - half_field) / (np.sqrt(2) * sigma)))
    return base * (1 + horns * (x / max(half_field, 1e-9)) ** 2) + 0.002


def analytic_pdd(depth, dmax=1.5, mu=0.05):
    """Build-up then exponential fall-off with a smooth maximum near dmax."""
    k = 3.0 / dmax
    return (1 - np.exp(-k * depth)) * np.exp(-mu * depth)


def fill_fake_campaign(work: Path, rng_seed: int = 0, noise: float = 0.002,
                       field_scale: dict | None = None) -> None:
    """Write fake TOPAS outputs next to every meta.json produced by make_water_campaign."""
    rng = np.random.default_rng(rng_seed)
    field_scale = field_scale or {}
    for meta_p in sorted(Path(work).glob("*/job*/meta.json")):
        meta = json.loads(meta_p.read_text())
        fs = meta["field_setting"]
        scale = field_scale.get(meta["field"]["name"], 1.0)
        n_hist = meta["histories"]
        for name, sc in meta["scorers"].items():
            jd = meta_p.parent
            if sc["kind"] == "pdd":
                d = np.asarray(sc["positions_cm"])
                vals = analytic_pdd(d) * 1e-16 * n_hist * scale
                shape = (1, 1, len(d))
            elif sc["kind"].startswith("profile_"):
                x = np.asarray(sc["positions_cm"])
                if sc["kind"] == "profile_x":
                    lo, hi = fs["x1"], fs["x2"]
                else:
                    lo, hi = fs["y1"], fs["y2"]
                c, h = 0.5 * (lo + hi), 0.5 * (hi - lo)
                vals = analytic_profile(x - c, h) * 1e-16 * n_hist * scale
                shape = (len(x), 1, 1) if sc["axis"] == "i" else (1, len(x), 1)
            else:
                vals = np.array([1e-16 * n_hist * scale])
                shape = (1, 1, 1)
            vals = vals * (1 + noise * rng.standard_normal(vals.shape))
            write_topas_csv(jd / sc["file"], shape, (0.1, 0.1, 0.1), vals)
