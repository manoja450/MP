import json
import re
from pathlib import Path

import pytest

from mcverify.config import ConfigError, load_machine, load_yaml, validate_machine
from mcverify.topas.linac import FieldSetting, SourceSpec, linac, mlc_positions
from mcverify.topas.params import ParamWriter
from mcverify.topas.water import ValidationField, fields_from_config, make_water_campaign

ROOT = Path(__file__).resolve().parents[1]
MACHINE = ROOT / "config" / "truebeam_6x_millennium.yaml"
VALID = ROOT / "config" / "validation_water.yaml"
LINE = re.compile(r'^(s|d|dc|u|i|b|sv|dv|uv|iv):([A-Za-z0-9_/+\-]+) = (.+)$')


def params_of(text: str) -> dict:
    """Parse generated TOPAS text into {name: value-string}; fail on malformed lines."""
    out = {}
    for line in text.splitlines():
        if not line.strip() or line.startswith("#") or line.startswith("includeFile"):
            continue
        m = LINE.match(line)
        assert m, f"malformed TOPAS line: {line!r}"
        name = m.group(2)
        assert name not in out, f"parameter defined twice: {name}"
        out[name] = (m.group(1), m.group(3))
        if m.group(1) in ("sv", "dv", "uv", "iv"):
            toks = m.group(3).replace('" "', '"|"').split()
            n = int(toks[0])
            vals = re.findall(r'"[^"]*"', m.group(3)) if m.group(1) == "sv" else toks[1:]
            if m.group(1) == "dv":
                vals = vals[:-1]  # unit
            assert len(vals) == n, f"vector length mismatch in {name}"
    return out


def test_writer_formats():
    w = ParamWriter().d("Ge/A/HLX", 2, "cm").dv("Ge/A/V", [1, 2.5], "mm").sv("Ph/M", ["a", "b"]).b("X/Y", True)
    t = w.text()
    assert "d:Ge/A/HLX = 2. cm" in t
    assert "dv:Ge/A/V = 2 1. 2.5 mm" in t
    assert 'sv:Ph/M = 2 "a" "b"' in t
    assert 'b:X/Y = "True"' in t
    with pytest.raises(ValueError):
        ParamWriter().s("x", 'bad"quote')


def test_machine_config_validates():
    cfg = load_machine(MACHINE)
    assert cfg["mlc"]["model"] == "millennium120"
    bad = dict(cfg, jaws={"upper": dict(cfg["jaws"]["upper"]), "lower": dict(cfg["jaws"]["upper"])})
    with pytest.raises(ConfigError):
        validate_machine(bad)


def test_linac_block_is_wellformed():
    cfg = load_machine(MACHINE)
    src = SourceSpec(stem="/x/TB_c000", histories=1000, z_stored=False, z_constant_cm=26.7)
    p = params_of(linac(cfg, FieldSetting.square(10), src).text())
    assert p["Ge/PhspPlane/TransZ"][1] == "73.3 cm"
    assert p["So/Linac/PhaseSpaceMultipleUse"][1] == "0"
    assert p["Ge/JawY/Type"][1] == '"TsJaws"' and p["Ge/JawX/Type"][1] == '"TsJaws"'
    assert p["Ge/JawY/RotZ"][1] == "-90. deg" and p["Ge/JawX/RotZ"][1] == "0. deg"
    assert p["Ge/JawX/NegativeFieldSetting"] == ("dc", "-5. cm")
    assert p["Ge/MLC/Type"][1] == '"TsDivergingMLC"'
    assert p["Ge/MLC/LeafWidths"][1].startswith("60 ")
    # jaw centre = SAD - upstream distance - half thickness
    assert float(p["Ge/JawY/TransZ"][1].split()[0]) == pytest.approx(100 - 28.0 - 3.9)


def test_phase_space_without_constant_z_needs_config():
    cfg = load_machine(MACHINE)
    with pytest.raises(ValueError):
        linac(cfg, FieldSetting.square(10), SourceSpec(stem="x", histories=1, z_stored=False, z_constant_cm=None))


def test_mlc_offset_and_limits():
    cfg = load_machine(MACHINE, {"mlc": {"leaf_end_offset_mm": 0.5}})
    fs = FieldSetting(-5, 5, -5, 5, mlc_neg=[-2.0] * 60, mlc_pos=[3.0] * 60)
    neg, pos = mlc_positions(cfg, fs)
    assert neg[0] == pytest.approx(-2.05) and pos[0] == pytest.approx(3.05)
    with pytest.raises(ValueError):
        mlc_positions(cfg, FieldSetting(-5, 5, -5, 5, mlc_neg=[1.0] * 60, mlc_pos=[0.0] * 60))
    with pytest.raises(ValueError):
        mlc_positions(cfg, FieldSetting(-5, 5, -5, 5, mlc_neg=[0.0] * 10, mlc_pos=[0.0] * 10))


def test_mlc_square_field_shape():
    cfg = load_machine(MACHINE)
    fs = ValidationField("m", "mlc_square", size_cm=3).field_setting(cfg)
    open_pairs = [i for i, (a, b) in enumerate(zip(fs.mlc_neg, fs.mlc_pos)) if b > a]
    assert len(open_pairs) == 6          # 3 cm / 0.5 cm central leaves
    assert fs.x2 == pytest.approx(2.5)   # jaws 1 cm wider than the MLC field


def test_campaign_generation(tmp_path):
    cfg = load_machine(MACHINE)
    vcfg = load_yaml(VALID)
    vcfg["jobs_per_field"] = 3
    chunks = [{"stem": f"/c/TB_c{k:03d}", "histories": 5000, "z_stored": False, "z_constant_cm": 26.7}
              for k in range(5)]
    runs = make_water_campaign(cfg, vcfg, chunks, tmp_path)
    n_fields = len(fields_from_config(vcfg))
    assert len(runs) == 3 * n_fields
    seeds, chunk_use = set(), {}
    for r in runs:
        p = params_of(r.read_text())
        meta = json.loads((r.parent / "meta.json").read_text())
        assert int(p["Ts/Seed"][1]) not in seeds
        seeds.add(int(p["Ts/Seed"][1]))
        chunk_use.setdefault(meta["field"]["name"], []).append(meta["phase_space_chunk"])
        # each scorer's output file lives in the job directory
        for name, sc in meta["scorers"].items():
            assert p[f"Sc/{name}/OutputFile"][1].strip('"') == str(r.parent / name)
        # PDD bins: deepest first (TOPAS bin 0 is the most negative z)
        assert meta["scorers"]["PDD"]["positions_cm"][0] > meta["scorers"]["PDD"]["positions_cm"][-1]
    # jobs of one field never share a phase-space chunk
    assert all(len(v) == len(set(v)) for v in chunk_use.values())


def test_too_few_chunks_is_an_error(tmp_path):
    cfg = load_machine(MACHINE)
    vcfg = load_yaml(VALID)
    vcfg["jobs_per_field"] = 4
    with pytest.raises(ValueError):
        make_water_campaign(cfg, vcfg, [{"stem": "a", "histories": 1, "z_constant_cm": 26.7}] * 2, tmp_path)
