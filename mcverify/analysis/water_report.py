"""Analyse a finished water-tank validation campaign.

Produces, in ``out_dir``:
  summary_pdd.csv, summary_profiles.csv, output_factors.csv, checks.csv
  figures/<field>_pdd.png, figures/<field>_profiles.png
  REPORT.md  (human-readable summary with pass/fail against acceptance criteria)
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np

from .combine import combine_scorer, job_dirs
from .curves import gamma_1d, normalise_pdd, pass_rate, pdd_metrics, profile_metrics
from .reference import ReferenceData, depth_tag

CRITERIA = {
    "pdd_gamma": (1.0, 1.0),        # 1% / 1 mm, global
    "profile_gamma": (1.0, 1.0),
    "pdd_pass_rate_min": 95.0,
    "profile_pass_rate_min": 95.0,
    "of_tol_pct": {4.0: 1.0, 2.0: 2.0, 0.0: 3.0},   # field >= key -> tolerance
}


def _of_tolerance(size_cm: float) -> float:
    for k in sorted(CRITERIA["of_tol_pct"], reverse=True):
        if size_cm >= k:
            return CRITERIA["of_tol_pct"][k]
    return 3.0


def _fmt(v, nd=3):
    if v is None or (isinstance(v, float) and (math.isnan(v) or math.isinf(v))):
        return ""
    return f"{v:.{nd}f}" if isinstance(v, float) else str(v)


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: _fmt(v) if isinstance(v, float) else v for k, v in r.items()})


def _scorer_names(field_dir: Path) -> list[str]:
    jd = job_dirs(field_dir)
    if not jd:
        return []
    return list(json.loads((jd[0] / "meta.json").read_text())["scorers"])


def analyse_campaign(work_dir: str | Path, out_dir: str | Path, ref_dir: str | Path | None = None,
                     reference_field: str = "jaw10cm", make_plots: bool = True) -> dict:
    work_dir, out_dir = Path(work_dir), Path(out_dir)
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)
    ref = ReferenceData(ref_dir) if ref_dir else None
    ref_of = ref.output_factors() if ref else {}
    pdd_rows, prof_rows, of_rows, checks = [], [], [], []
    of_dose = {}
    fields = sorted(p for p in work_dir.iterdir() if p.is_dir() and job_dirs(p))

    for fdir in fields:
        name = fdir.name
        scorers = _scorer_names(fdir)
        try:
            pdd = combine_scorer(fdir, "PDD")
        except FileNotFoundError as e:
            checks.append({"field": name, "check": "outputs", "result": f"MISSING: {e}"})
            continue
        size = float(pdd.meta["nominal_size_cm"])
        kind = pdd.meta["field"]["kind"]

        # ---- PDD -------------------------------------------------------
        row = {"field": name, "kind": kind, "size_cm": size, "jobs": pdd.n_jobs,
               "histories": pdd.total_histories}
        if kind in ("jaw_square", "mlc_square") and pdd.dose_per_history.max() <= 0:
            checks.append({"field": name, "check": "PDD", "result": "FAIL - no dose scored in PDD column"})
        elif kind in ("jaw_square", "mlc_square"):
            m = pdd_metrics(pdd.positions_cm, pdd.dose_per_history)
            pct, dmax, dpk = normalise_pdd(pdd.positions_cm, pdd.dose_per_history)
            hi = pdd.dose_per_history > 0.5 * dpk
            row.update({"dmax_cm": m.dmax_cm, "pdd10": m.pdd10, "pdd20": m.pdd20, "d20_d10": m.d20_d10,
                        "median_rel_unc_pct": float(np.median(pdd.rel_unc[hi]) * 100)})
            r = ref.pdd(name) if ref else None
            if r is not None:
                rpct, rdmax, _ = normalise_pdd(*r)
                g = gamma_1d(r[0], rpct, pdd.positions_cm, pct, *CRITERIA["pdd_gamma"], norm_value=100.0)
                beyond = r[0] >= rdmax
                diff = np.interp(r[0][beyond], pdd.positions_cm, pct) - rpct[beyond]
                row.update({"ref_dmax_cm": rdmax, "gamma_1pct1mm_pass": pass_rate(g),
                            "max_abs_diff_beyond_dmax_pct": float(np.max(np.abs(diff))) if diff.size else float("nan")})
            if make_plots:
                _plot_pdd(out_dir / "figures" / f"{name}_pdd.png", name, pdd.positions_cm, pct, r)
            np.savetxt(out_dir / f"{name}_pdd.csv",
                       np.column_stack([pdd.positions_cm, pdd.dose_per_history, pdd.rel_unc, pct]),
                       delimiter=",", header="depth_cm,dose_per_history_Gy,rel_unc,pdd_pct", comments="")
        pdd_rows.append(row)

        # ---- profiles ----------------------------------------------------
        prof_curves = []
        for sname in [s for s in scorers if s.startswith("Prof")]:
            c = combine_scorer(fdir, sname)
            axis = c.meta["scorer"]["kind"][-1]
            depth = float(c.meta["scorer"]["depth_cm"])
            prow = {"field": name, "axis": axis, "depth_cm": depth}
            if kind in ("jaw_square", "mlc_square"):
                try:
                    pm = profile_metrics(c.positions_cm, c.dose_per_history)
                except ValueError as e:
                    checks.append({"field": name, "check": f"profile {sname}", "result": f"FAIL - {e}",
                                   "value": f"total dose {c.dose_per_history.sum():.3e} Gy/history"})
                    prof_rows.append(prow)
                    continue
                prow.update({"field_size_cm": pm.field_size_cm, "centre_cm": pm.centre_cm,
                             "penumbra_left_cm": pm.penumbra_left_cm, "penumbra_right_cm": pm.penumbra_right_cm,
                             "flatness_pct": pm.flatness_pct, "symmetry_pct": pm.symmetry_pct})
                norm = 100.0 * c.dose_per_history / pm.cax_value
                r = ref.profile(axis, name, depth) if ref else None
                if r is not None:
                    rp = profile_metrics(*r)
                    rnorm = 100.0 * r[1] / rp.cax_value
                    g = gamma_1d(r[0], rnorm, c.positions_cm, norm, *CRITERIA["profile_gamma"], norm_value=100.0)
                    prow.update({"ref_field_size_cm": rp.field_size_cm,
                                 "ref_penumbra_left_cm": rp.penumbra_left_cm,
                                 "gamma_1pct1mm_pass": pass_rate(g)})
                prof_curves.append((axis, depth, c.positions_cm, norm, r))
            elif kind == "asym_jaw":
                w = np.clip(c.dose_per_history, 0, None)
                centroid = float(np.sum(c.positions_cm * w) / np.sum(w)) if w.sum() > 0 else float("nan")
                prow["dose_centroid_cm"] = centroid
                fs = pdd.meta["field"]
                expect = 0.5 * ((fs["x1"] + fs["x2"]) if axis == "x" else (fs["y1"] + fs["y2"]))
                if abs(expect) > 0.5 and depth == 10.0:
                    ok = np.sign(centroid) == np.sign(expect)
                    checks.append({"field": name, "check": f"orientation {axis.upper()}",
                                   "result": "PASS" if ok else "FAIL - flip axis/rotation",
                                   "value": f"centroid {centroid:+.2f} cm, expected sign of {expect:+.2f}"})
            prof_rows.append(prow)
        if make_plots and prof_curves:
            _plot_profiles(out_dir / "figures" / f"{name}_profiles.png", name, prof_curves)

        # ---- output-factor voxel ------------------------------------------
        if "OF" in scorers:
            o = combine_scorer(fdir, "OF")
            of_dose[name] = (float(o.dose_per_history[0]), float(o.rel_unc[0]), size, kind)

    # ---- output factors, transmission ------------------------------------
    if reference_field in of_dose:
        d_ref, u_ref, _, _ = of_dose[reference_field]
        for name, (d, u, size, kind) in sorted(of_dose.items(), key=lambda t: t[1][2]):
            of = d / d_ref
            unc = of * math.sqrt(u ** 2 + u_ref ** 2) if name != reference_field else 0.0
            r = {"field": name, "kind": kind, "size_cm": size, "output_factor": of, "unc_1sigma": unc}
            if name in ref_of:
                diff = 100.0 * (of / ref_of[name] - 1.0)
                tol = _of_tolerance(size)
                r.update({"ref_output_factor": ref_of[name], "diff_pct": diff, "tolerance_pct": tol,
                          "result": "PASS" if abs(diff) <= tol else "FAIL"})
            of_rows.append(r)
            if kind == "mlc_closed":
                checks.append({"field": name, "check": "MLC transmission (closed / open 10x10 at 10 cm)",
                               "result": f"{100 * of:.2f} %", "value": f"+- {100 * unc:.2f} %"})
    else:
        checks.append({"field": reference_field, "check": "output factors",
                       "result": "SKIPPED - reference field not simulated"})

    _write_csv(out_dir / "summary_pdd.csv", pdd_rows)
    _write_csv(out_dir / "summary_profiles.csv", prof_rows)
    _write_csv(out_dir / "output_factors.csv", of_rows)
    _write_csv(out_dir / "checks.csv", checks)
    _write_report(out_dir / "REPORT.md", pdd_rows, prof_rows, of_rows, checks, ref_dir)
    return {"pdd": pdd_rows, "profiles": prof_rows, "output_factors": of_rows, "checks": checks}


# ---------------------------------------------------------------------------
def _plot_pdd(path, name, depth, pct, ref):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a, b) = plt.subplots(2, 1, figsize=(6, 5), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    a.plot(depth, pct, label="Monte Carlo", lw=1.5)
    if ref is not None:
        rp, _, _ = normalise_pdd(*ref)
        a.plot(ref[0], rp, "k.", ms=3, label="Reference")
        b.plot(ref[0], np.interp(ref[0], depth, pct) - rp, ".", ms=3)
    b.axhline(0, color="grey", lw=0.8)
    for lim in (-1, 1):
        b.axhline(lim, color="grey", ls=":", lw=0.8)
    a.set_ylabel("PDD (%)"); a.set_title(name); a.legend(frameon=False)
    b.set_xlabel("Depth (cm)"); b.set_ylabel("MC - ref (%)")
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def _plot_profiles(path, name, curves):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    axes = sorted({c[0] for c in curves})
    fig, ax = plt.subplots(1, len(axes), figsize=(6 * len(axes), 4), squeeze=False)
    for col, axis in enumerate(axes):
        a = ax[0, col]
        for (ax_, depth, x, y, ref) in curves:
            if ax_ != axis:
                continue
            line, = a.plot(x, y, lw=1.2, label=f"{depth:g} cm")
            if ref is not None:
                rp = profile_metrics(*ref)
                a.plot(ref[0], 100 * ref[1] / rp.cax_value, ".", ms=2, color=line.get_color())
        a.set_xlabel(f"{'Crossline X' if axis == 'x' else 'Inline Y'} (cm)")
        a.set_ylabel("Dose (% of CAX)"); a.set_title(f"{name}  (lines MC, dots reference)")
        a.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def _write_report(path, pdd_rows, prof_rows, of_rows, checks, ref_dir):
    L = ["# Water-tank validation report", ""]
    L.append(f"Reference data: `{ref_dir}`" if ref_dir else "Reference data: none supplied (MC-only metrics)")
    L += ["", "Acceptance criteria (project plan): PDD and profiles 1%/1 mm gamma >= 95%; "
          "output factors within 1% (>=4x4), 2% (2x2), 3% (1x1).", "", "## PDDs", "",
          "| field | jobs | histories | dmax (cm) | PDD10 | PDD20 | median unc (%) | gamma 1%/1mm | max diff >dmax (%) |",
          "|---|---|---|---|---|---|---|---|---|"]
    for r in pdd_rows:
        L.append(f"| {r['field']} | {r['jobs']} | {r['histories']:.3g} | {_fmt(r.get('dmax_cm'), 2)} | "
                 f"{_fmt(r.get('pdd10'), 1)} | {_fmt(r.get('pdd20'), 1)} | {_fmt(r.get('median_rel_unc_pct'), 2)} | "
                 f"{_fmt(r.get('gamma_1pct1mm_pass'), 1)} | {_fmt(r.get('max_abs_diff_beyond_dmax_pct'), 2)} |")
    L += ["", "## Profiles", "", "| field | axis | depth | field size | pen. L | pen. R | flatness | symmetry | gamma |",
          "|---|---|---|---|---|---|---|---|---|"]
    for r in prof_rows:
        L.append(f"| {r['field']} | {r['axis']} | {r['depth_cm']:g} | {_fmt(r.get('field_size_cm'), 2)} | "
                 f"{_fmt(r.get('penumbra_left_cm'), 2)} | {_fmt(r.get('penumbra_right_cm'), 2)} | "
                 f"{_fmt(r.get('flatness_pct'), 2)} | {_fmt(r.get('symmetry_pct'), 2)} | "
                 f"{_fmt(r.get('gamma_1pct1mm_pass'), 1)} |")
    L += ["", "## Output factors (relative to jaw10cm, 10 cm depth)", "",
          "| field | OF | +-1 sigma | reference | diff (%) | tol (%) | result |", "|---|---|---|---|---|---|---|"]
    for r in of_rows:
        L.append(f"| {r['field']} | {_fmt(r['output_factor'], 4)} | {_fmt(r['unc_1sigma'], 4)} | "
                 f"{_fmt(r.get('ref_output_factor'), 4)} | {_fmt(r.get('diff_pct'), 2)} | "
                 f"{_fmt(r.get('tolerance_pct'), 1)} | {r.get('result', '')} |")
    L += ["", "## Checks", ""] + [f"- **{c['field']}** {c['check']}: {c['result']} {c.get('value', '')}" for c in checks]
    Path(path).write_text("\n".join(L) + "\n")
