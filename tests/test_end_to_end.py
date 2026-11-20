"""Generate a campaign, fill it with synthetic 'TOPAS' outputs with known answers,
then run the full analysis and check the report."""
from pathlib import Path

import numpy as np
import pytest

from mcverify.analysis.water_report import analyse_campaign
from mcverify.config import load_machine, load_yaml
from mcverify.dicom.rtplan import read_rtplan, summarize
from mcverify.dicom.synthetic import make_arc_plan
from mcverify.topas.water import make_water_campaign
from tests.helpers import analytic_pdd, analytic_profile, fill_fake_campaign

ROOT = Path(__file__).resolve().parents[1]


def _campaign(tmp_path, fields_cm=(3, 10)):
    cfg = load_machine(ROOT / "config" / "truebeam_6x_millennium.yaml")
    vcfg = load_yaml(ROOT / "config" / "validation_water.yaml")
    vcfg["square_fields_cm"] = list(fields_cm)
    vcfg["jobs_per_field"] = 5
    chunks = [{"stem": f"/c/TB_c{k:03d}", "histories": 10000, "z_stored": False, "z_constant_cm": 26.7}
              for k in range(6)]
    make_water_campaign(cfg, vcfg, chunks, tmp_path / "work")
    return cfg, vcfg


def test_full_analysis_on_synthetic_outputs(tmp_path):
    _campaign(tmp_path)
    fill_fake_campaign(tmp_path / "work", field_scale={"jaw3cm": 0.85, "mlcClosed10cm": 0.015})
    # reference data identical to the truth -> everything must pass
    ref = tmp_path / "ref"
    ref.mkdir()
    d = np.arange(0.0, 30, 0.2)
    for name in ("jaw10cm", "jaw3cm"):
        np.savetxt(ref / f"pdd_{name}.csv", np.column_stack([d, analytic_pdd(d)]), delimiter=",",
                   header="depth_cm,dose", comments="")
    x = np.arange(-10, 10, 0.1)
    for dep in ("1p5cm", "5cm", "10cm", "20cm"):
        np.savetxt(ref / f"profile_x_jaw10cm_d{dep}.csv", np.column_stack([x, analytic_profile(x, 5.0)]),
                   delimiter=",", header="position_cm,dose", comments="")
    (ref / "output_factors.csv").write_text("field,output_factor\njaw10cm,1.0\njaw3cm,0.85\n")

    out = tmp_path / "results"
    res = analyse_campaign(tmp_path / "work", out, ref, make_plots=True)
    report = (out / "REPORT.md").read_text()
    assert "jaw10cm" in report and (out / "figures" / "jaw10cm_pdd.png").exists()

    pdd = {r["field"]: r for r in res["pdd"]}
    assert pdd["jaw10cm"]["gamma_1pct1mm_pass"] >= 95.0
    assert pdd["jaw10cm"]["jobs"] == 5
    prof = [r for r in res["profiles"] if r["field"] == "jaw10cm" and r["axis"] == "x"]
    assert all(abs(r["field_size_cm"] - 10.0) < 0.05 for r in prof)
    assert all(r["gamma_1pct1mm_pass"] >= 95.0 for r in prof)

    ofs = {r["field"]: r for r in res["output_factors"]}
    assert ofs["jaw3cm"]["output_factor"] == pytest.approx(0.85, rel=0.01)
    assert ofs["jaw3cm"]["result"] == "PASS"
    checks = {c["check"]: c for c in res["checks"]}
    assert checks["orientation Y"]["result"] == "PASS"
    assert checks["orientation X"]["result"] == "PASS"
    assert "1.5" in checks["MLC transmission (closed / open 10x10 at 10 cm)"]["result"]


def test_output_factor_failure_is_flagged(tmp_path):
    _campaign(tmp_path, fields_cm=(3, 10))
    fill_fake_campaign(tmp_path / "work", field_scale={"jaw3cm": 0.80})
    ref = tmp_path / "ref"
    ref.mkdir()
    (ref / "output_factors.csv").write_text("field,output_factor\njaw3cm,0.85\n")
    res = analyse_campaign(tmp_path / "work", tmp_path / "out", ref, make_plots=False)
    of = {r["field"]: r for r in res["output_factors"]}
    assert of["jaw3cm"]["result"] == "FAIL"


def test_rtplan_reader_on_synthetic_arc(tmp_path):
    p = make_arc_plan(tmp_path / "plan.dcm", mu=400, n_cp=9, start=181, stop=179, field_mm=30,
                      leaf_shift_mm_per_cp=1.0, fff=True)
    plan = read_rtplan(p)
    b = plan.beams[0]
    assert len(plan.beams) == 1 and b.mu == 400 and b.n_leaf_pairs == 60
    assert b.fluence_mode.startswith("NON_STANDARD")
    assert b.is_arc()
    assert b.segment_mu().sum() == pytest.approx(400)
    assert b.control_points[0].jaw_x_mm == (-25.0, 25.0)
    # jaws are only stored at CP0 and must be inherited
    assert b.control_points[5].jaw_y_mm == (-15.0, 15.0)
    # leaves drift +1 mm per control point
    i = 30
    assert b.control_points[4].mlc_pos_mm[i] - b.control_points[0].mlc_pos_mm[i] == pytest.approx(4.0)
    assert "ARC1" in summarize(plan)


def test_sparse_outputs_are_flagged_not_crashing(tmp_path):
    _campaign(tmp_path, fields_cm=(10,))
    fill_fake_campaign(tmp_path / "work")
    # make the PDD of jaw10cm mostly empty, like a run with far too few histories
    import json
    from tests.helpers import write_topas_csv
    for meta_p in (tmp_path / "work" / "jaw10cm").glob("job*/meta.json"):
        meta = json.loads(meta_p.read_text())
        n = len(meta["scorers"]["PDD"]["positions_cm"])
        vals = np.zeros(n)
        vals[::7] = 1e-12
        write_topas_csv(meta_p.parent / "PDD.csv", (1, 1, n), (0.1, 0.1, 0.2), vals)
    res = analyse_campaign(tmp_path / "work", tmp_path / "out", None, make_plots=False)
    assert any(c["check"] == "statistics" and c["result"].startswith("LOW") for c in res["checks"])
