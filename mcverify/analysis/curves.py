"""Beam-data metrics and 1D comparisons (PDDs and profiles).

Conventions
-----------
* depth / position arrays are in cm and must be increasing
* PDDs are normalised to 100% at d_max (found from a local quadratic fit,
  which is robust against Monte Carlo noise)
* profiles are normalised to 100% on the central axis (CAX), the standard
  choice for flattened beams; for FFF beams report the same quantities but
  read the field size from the 50% of CAX points (see Fogliata et al. 2012,
  renormalisation methods, if you need the formal FFF definitions)
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


# ---------------------------------------------------------------------------
# basic helpers
# ---------------------------------------------------------------------------
def _check(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if x.ndim != 1 or x.shape != y.shape:
        raise ValueError("x and y must be 1-D arrays of equal length")
    if np.any(np.diff(x) <= 0):
        raise ValueError("positions must be strictly increasing")
    return x, y


def crossing(x, y, level, side: str) -> float:
    """Position where y crosses ``level`` by linear interpolation.

    side='left' searches from the curve maximum towards smaller x,
    side='right' towards larger x. Returns nan if there is no crossing.
    """
    x, y = _check(x, y)
    i0 = int(np.argmax(y))
    rng = range(i0, 0, -1) if side == "left" else range(i0, len(x) - 1)
    for i in rng:
        a, b = (i - 1, i) if side == "left" else (i, i + 1)
        ya, yb = y[a], y[b]
        if (ya - level) * (yb - level) <= 0 and ya != yb:
            return float(x[a] + (level - ya) * (x[b] - x[a]) / (yb - ya))
    return float("nan")


def smooth(y, window: int = 5) -> np.ndarray:
    """Simple moving average (odd window) used only to locate maxima."""
    y = np.asarray(y, float)
    if window <= 1 or len(y) < window:
        return y.copy()
    k = np.ones(window) / window
    pad = window // 2
    return np.convolve(np.pad(y, pad, mode="edge"), k, mode="valid")


def peak_quadratic(x, y, half_width: float) -> tuple[float, float]:
    """Location and height of the maximum from a local polynomial fit around argmax.

    A cubic is used when enough points are available because depth-dose peaks
    are asymmetric (steep build-up, slow fall-off) and a parabola would bias
    d_max towards the slow side. Falls back to a parabola, then to argmax.
    """
    x, y = _check(x, y)
    i0 = int(np.argmax(smooth(y, 5)))
    sel = np.abs(x - x[i0]) <= half_width
    xs, ys = x[sel], y[sel]
    if xs.size < 3:
        return float(x[i0]), float(y[i0])
    lo, hi = xs.min(), xs.max()
    if xs.size >= 6:
        c3 = np.polyfit(xs - x[i0], ys, 3)
        roots = np.roots(np.polyder(c3))
        cands = [r.real + x[i0] for r in roots if abs(r.imag) < 1e-12
                 and np.polyval(np.polyder(c3, 2), r.real) < 0 and lo <= r.real + x[i0] <= hi]
        if cands:
            xm = min(cands, key=lambda v: abs(v - x[i0]))
            return float(xm), float(np.polyval(c3, xm - x[i0]))
    a, b, c = np.polyfit(xs, ys, 2)
    if a >= 0:
        return float(x[i0]), float(y[i0])
    xm = float(np.clip(-b / (2 * a), lo, hi))
    return xm, float(np.polyval([a, b, c], xm))


# ---------------------------------------------------------------------------
# PDD
# ---------------------------------------------------------------------------
@dataclass
class PDDMetrics:
    dmax_cm: float
    pdd10: float
    pdd20: float
    d20_d10: float


def normalise_pdd(depth, dose, fit_half_width_cm: float = 0.6) -> tuple[np.ndarray, float, float]:
    """Return (pdd in %, d_max, dose at d_max)."""
    depth, dose = _check(depth, dose)
    dmax, dpeak = peak_quadratic(depth, dose, fit_half_width_cm)
    return 100.0 * dose / dpeak, dmax, dpeak


def pdd_metrics(depth, dose) -> PDDMetrics:
    pdd, dmax, _ = normalise_pdd(depth, dose)
    p10 = float(np.interp(10.0, depth, pdd))
    p20 = float(np.interp(20.0, depth, pdd))
    return PDDMetrics(dmax, p10, p20, p20 / p10)


# ---------------------------------------------------------------------------
# Profiles
# ---------------------------------------------------------------------------
@dataclass
class ProfileMetrics:
    cax_value: float
    field_size_cm: float
    centre_cm: float
    penumbra_left_cm: float
    penumbra_right_cm: float
    flatness_pct: float
    symmetry_pct: float


def cax_value(x, y, half_width_cm: float = 0.3) -> float:
    """Mean dose within +-half_width of x = 0 (reduces MC noise)."""
    x, y = _check(x, y)
    sel = np.abs(x) <= half_width_cm
    if sel.sum() == 0:
        return float(np.interp(0.0, x, y))
    return float(np.mean(y[sel]))


def profile_metrics(x, y, cax_half_width_cm: float = 0.3) -> ProfileMetrics:
    x, y = _check(x, y)
    c = cax_value(x, y, cax_half_width_cm)
    if c <= 0:
        raise ValueError("central-axis dose is zero; is this an off-axis field?")
    yn = 100.0 * y / c
    ys = smooth(yn, 3)
    l50, r50 = crossing(x, ys, 50, "left"), crossing(x, ys, 50, "right")
    l80, r80 = crossing(x, ys, 80, "left"), crossing(x, ys, 80, "right")
    l20, r20 = crossing(x, ys, 20, "left"), crossing(x, ys, 20, "right")
    fs = r50 - l50
    centre = 0.5 * (l50 + r50)
    # flatness and symmetry inside the central 80% of the field
    inner = np.abs(x - centre) <= 0.4 * fs
    if inner.sum() >= 3:
        dmax_i, dmin_i = ys[inner].max(), ys[inner].min()
        flat = 100.0 * (dmax_i - dmin_i) / (dmax_i + dmin_i)
        xi = x[inner] - centre
        mirror = np.interp(-xi, x - centre, ys)
        sym = float(np.max(np.abs(ys[inner] - mirror) / ys[inner]) * 100.0)
    else:
        flat, sym = float("nan"), float("nan")
    return ProfileMetrics(c, fs, centre, l80 - l20, r20 - r80, float(flat), sym)


# ---------------------------------------------------------------------------
# 1-D gamma
# ---------------------------------------------------------------------------
def gamma_1d(x_ref, d_ref, x_eval, d_eval, dose_pct: float = 1.0, dta_mm: float = 1.0,
             threshold_pct: float = 10.0, global_norm: bool = True,
             norm_value: float | None = None, resample_factor: int = 10) -> np.ndarray:
    """Gamma index for every reference point (nan below the low-dose threshold).

    The evaluated curve is linearly interpolated onto a grid ``resample_factor``
    times finer than the DTA so that distance-to-agreement is resolved.
    Global gamma uses ``norm_value`` (default: max of the reference).
    """
    x_ref, d_ref = _check(x_ref, d_ref)
    x_eval, d_eval = _check(x_eval, d_eval)
    dta = dta_mm / 10.0
    norm = float(norm_value if norm_value is not None else d_ref.max())
    step = dta / resample_factor
    xf = np.arange(x_eval[0], x_eval[-1] + step / 2, step)
    df = np.interp(xf, x_eval, d_eval)
    gam = np.full(x_ref.shape, np.nan)
    search = 3.0 * dta
    for n, (xr, dr) in enumerate(zip(x_ref, d_ref)):
        if dr < threshold_pct / 100.0 * norm:
            continue
        if xr < x_eval[0] - search or xr > x_eval[-1] + search:
            continue
        sel = np.abs(xf - xr) <= search
        if not sel.any():
            continue
        dd = (dose_pct / 100.0) * (norm if global_norm else max(dr, 1e-12))
        g2 = ((xf[sel] - xr) / dta) ** 2 + ((df[sel] - dr) / dd) ** 2
        gam[n] = np.sqrt(g2.min())
    return gam


def pass_rate(gamma: np.ndarray) -> float:
    g = gamma[~np.isnan(gamma)]
    return float(100.0 * np.mean(g <= 1.0)) if g.size else float("nan")
