"""Water-tank validation runs (Paper 1, Aim 1): PDDs, profiles, output factors,
MLC transmission and MLC-defined fields.

Each job is a self-contained TOPAS parameter file plus a ``meta.json`` that
records how scorer bins map to depth / off-axis position, so the analysis code
never has to guess.

Scoring components live in separate *parallel worlds* (so the PDD column and the
profile strips may cross each other). Dose is scored as DoseToMedium; the mass
world material is water, so this equals dose to water.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ..config import leaf_widths_cm
from .linac import FieldSetting, SourceSpec, linac, run_control
from .params import ParamWriter


@dataclass
class ValidationField:
    name: str
    kind: str                    # jaw_square | mlc_square | asym_jaw | mlc_closed
    size_cm: float | None = None
    x1: float | None = None
    x2: float | None = None
    y1: float | None = None
    y2: float | None = None
    jaw_margin_cm: float = 1.0   # mlc_square: jaws this much wider than the MLC field
    closed_at_cm: float = -6.0   # mlc_closed: where the closed leaf pairs meet

    def field_setting(self, cfg: dict) -> FieldSetting:
        widths = np.array(leaf_widths_cm(cfg))
        edges = np.concatenate([[0.0], np.cumsum(widths)]) - widths.sum() / 2.0
        centres = 0.5 * (edges[:-1] + edges[1:])
        if self.kind == "jaw_square":
            return FieldSetting.square(self.size_cm)
        if self.kind == "asym_jaw":
            return FieldSetting(self.x1, self.x2, self.y1, self.y2)
        if self.kind == "mlc_square":
            h = self.size_cm / 2.0
            inside = np.abs(centres) < h
            neg = np.where(inside, -h, 0.0)
            pos = np.where(inside, h, 0.0)
            m = h + self.jaw_margin_cm
            return FieldSetting(-m, m, -h, h, mlc_neg=neg.tolist(), mlc_pos=pos.tolist())
        if self.kind == "mlc_closed":
            h = self.size_cm / 2.0
            c = self.closed_at_cm
            n = len(widths)
            return FieldSetting(-h, h, -h, h, mlc_neg=[c] * n, mlc_pos=[c] * n)
        raise ValueError(f"unknown field kind {self.kind}")

    def nominal_size_cm(self) -> float:
        if self.size_cm is not None:
            return float(self.size_cm)
        return float(min(self.x2 - self.x1, self.y2 - self.y1))


def fields_from_config(vcfg: dict) -> list[ValidationField]:
    out = [ValidationField(name=f"jaw{_tag(s)}", kind="jaw_square", size_cm=float(s))
           for s in vcfg.get("square_fields_cm", [])]
    for f in vcfg.get("extra_fields", []):
        out.append(ValidationField(**f))
    names = [f.name for f in out]
    if len(names) != len(set(names)):
        raise ValueError(f"duplicate field names: {names}")
    return out


def _tag(size: float) -> str:
    return f"{size:g}".replace(".", "p") + "cm"


# ---------------------------------------------------------------------------
def _resolution(vcfg: dict, size_cm: float) -> dict:
    p = vcfg["profiles"]
    small = size_cm <= p.get("small_field_limit_cm", 6.0)
    bin_mm = p["bin_mm_small"] if small else p["bin_mm_large"]
    thick_mm = min(p["thickness_mm"], max(1.0, 2.0 * size_cm))   # thin strips for tiny fields
    col_half = vcfg["pdd"].get("column_half_width_cm")
    if col_half in (None, "auto"):
        # 2x2 cm column for large fields (like a large ion chamber), shrinking for small fields
        col_half = min(1.0, 0.15 * size_cm)
    return {"bin_mm": bin_mm, "thick_mm": thick_mm, "col_half_cm": float(col_half)}


def _box(w: ParamWriter, name: str, parent: str, world: str, hl: tuple, bins: tuple, trans_z: float) -> None:
    w.s(f"Ge/{name}/Type", "TsBox")
    w.s(f"Ge/{name}/Parent", parent)
    # NO material: a parallel world with material switches on layered mass
    # geometry in OpenTOPAS. Dose uses the mass-world material (the water tank).
    w.b(f"Ge/{name}/IsParallel", True)
    w.s(f"Ge/{name}/ParallelWorldName", world)
    for ax, v in zip(("HLX", "HLY", "HLZ"), hl):
        w.d(f"Ge/{name}/{ax}", v, "cm")
    for ax, n in zip(("XBins", "YBins", "ZBins"), bins):
        w.i(f"Ge/{name}/{ax}", n)
    w.d(f"Ge/{name}/TransZ", trans_z, "cm")


def _scorer(w: ParamWriter, name: str, component: str, out_stem: Path) -> None:
    w.s(f"Sc/{name}/Quantity", "DoseToMedium")
    w.s(f"Sc/{name}/Component", component)
    w.s(f"Sc/{name}/OutputType", "csv")
    w.s(f"Sc/{name}/OutputFile", str(out_stem))
    w.s(f"Sc/{name}/IfOutputFileAlreadyExists", "Overwrite")
    w.sv(f"Sc/{name}/Report", ["Sum", "Standard_Deviation"])


def water_phantom_and_scorers(w: ParamWriter, cfg: dict, vcfg: dict, vf: ValidationField,
                              out_dir: Path) -> dict:
    """Write the tank and scorers; return scorer metadata (bin positions in cm)."""
    tank = vcfg["tank"]
    ssd = float(vcfg["ssd_cm"])
    sad = float(cfg["sad_cm"])
    hlz_tank = tank["depth_cm"] / 2.0
    surface_z = sad - ssd                      # z of water surface in IEC_F
    size = vf.nominal_size_cm()
    res = _resolution(vcfg, size)
    meta = {"ssd_cm": ssd, "resolution": res, "scorers": {}}

    w.section(f"Water tank (surface at SSD = {ssd} cm)")
    w.s("Ge/WaterTank/Type", "TsBox")
    w.s("Ge/WaterTank/Parent", "IEC_F")
    w.s("Ge/WaterTank/Material", "G4_WATER")
    w.d("Ge/WaterTank/HLX", tank["half_xy_cm"], "cm")
    w.d("Ge/WaterTank/HLY", tank["half_xy_cm"], "cm")
    w.d("Ge/WaterTank/HLZ", hlz_tank, "cm")
    w.d("Ge/WaterTank/TransZ", surface_z - hlz_tank, "cm")

    # depth d  <->  z (relative to tank centre) = hlz_tank - d
    w.section("PDD column (parallel world PDDWorld)")
    pdd = vcfg["pdd"]
    zmax = float(pdd["max_depth_cm"])
    nz = int(round(zmax * 10.0 / pdd["bin_mm"]))
    dz = zmax / nz
    _box(w, "PDDColumn", "WaterTank", "PDDWorld", (res["col_half_cm"], res["col_half_cm"], zmax / 2.0),
         (1, 1, nz), hlz_tank - zmax / 2.0)
    _scorer(w, "PDD", "PDDColumn", out_dir / "PDD")
    # TOPAS bin k=0 is at the most negative z, i.e. the DEEPEST bin
    depths = zmax - (np.arange(nz) + 0.5) * dz
    meta["scorers"]["PDD"] = {"file": "PDD.csv", "axis": "k", "kind": "pdd",
                              "positions_cm": depths.round(6).tolist()}

    prof = vcfg["profiles"]
    half_len = size / 2.0 + prof.get("margin_cm", 5.0)
    nb = int(round(2 * half_len * 10.0 / res["bin_mm"]))
    half_len = nb * res["bin_mm"] / 20.0
    xs = (-half_len + (np.arange(nb) + 0.5) * 2 * half_len / nb).round(6).tolist()
    t = res["thick_mm"] / 20.0
    for axis, world in (("X", "XProfWorld"), ("Y", "YProfWorld")):
        w.section(f"{'Crossline' if axis == 'X' else 'Inline'} profiles (parallel world {world})")
        for d in prof["depths_cm"]:
            name = f"Prof{axis}_d{_tag(float(d))}"
            hl = (half_len, t, t) if axis == "X" else (t, half_len, t)
            bins = (nb, 1, 1) if axis == "X" else (1, nb, 1)
            _box(w, name, "WaterTank", world, hl, bins, hlz_tank - float(d))
            _scorer(w, name, name, out_dir / name)
            meta["scorers"][name] = {"file": f"{name}.csv", "axis": "i" if axis == "X" else "j",
                                     "kind": f"profile_{axis.lower()}", "depth_cm": float(d),
                                     "positions_cm": xs}

    w.section("Output-factor voxel (parallel world OFWorld)")
    of = vcfg["output_factor"]
    h = min(of["voxel_mm"] / 20.0, 0.15 * size)
    _box(w, "OFVoxel", "WaterTank", "OFWorld", (h, h, h), (1, 1, 1), hlz_tank - float(of["depth_cm"]))
    _scorer(w, "OF", "OFVoxel", out_dir / "OF")
    meta["scorers"]["OF"] = {"file": "OF.csv", "axis": None, "kind": "point",
                             "depth_cm": float(of["depth_cm"]), "half_size_cm": h}
    return meta


def make_water_job(cfg: dict, vcfg: dict, vf: ValidationField, chunk: dict, job_index: int,
                   seed: int, job_dir: Path, histories: int | None = None) -> Path:
    """Write ``run.txt`` and ``meta.json`` for one job; return the run file path."""
    job_dir = Path(job_dir).resolve()
    job_dir.mkdir(parents=True, exist_ok=True)
    fs = vf.field_setting(cfg)
    src = SourceSpec(stem=chunk["stem"], histories=int(histories or chunk["histories"]),
                     z_stored=bool(chunk.get("z_stored", False)), z_constant_cm=chunk.get("z_constant_cm"))
    w = ParamWriter()
    w.comment(f"Auto-generated by mcverify -- machine {cfg['name']}, field {vf.name}, job {job_index}")
    w.comment("Do not edit by hand: change the YAML config and regenerate.")
    w.extend(linac(cfg, fs, src))
    meta = water_phantom_and_scorers(w, cfg, vcfg, vf, job_dir)
    run_control(w, seed=seed, threads=int(vcfg.get("threads_per_job", 1)))
    run = job_dir / "run.txt"
    w.write(run)
    meta.update({"machine": cfg["name"], "field": asdict(vf), "nominal_size_cm": vf.nominal_size_cm(),
                 "phsp_recycle": int(vcfg.get("phsp_recycle", 1)),
                 "job_index": job_index, "seed": seed, "histories": src.histories,
                 "phase_space_chunk": chunk["stem"],
                 "field_setting": asdict(fs)})
    (job_dir / "meta.json").write_text(json.dumps(meta, indent=1))
    return run


def make_water_campaign(cfg: dict, vcfg: dict, chunks: list[dict], work_dir: Path,
                        base_seed: int = 1000) -> list[Path]:
    """Create every job for every validation field. Returns the list of run files.

    Jobs cycle through the phase-space chunks. Within one field each job uses a
    different chunk, so jobs of the same field are statistically independent.
    """
    fields = fields_from_config(vcfg)
    recycle = int(vcfg.get("phsp_recycle", 1))
    if recycle < 1:
        raise ValueError("phsp_recycle must be >= 1")
    jobs_per_field = vcfg.get("jobs_per_field", 10)
    if jobs_per_field in ("all", None):
        jobs_per_field = len(chunks)
    jobs_per_field = int(jobs_per_field)
    if jobs_per_field > len(chunks):
        raise ValueError(f"jobs_per_field={jobs_per_field} but only {len(chunks)} phase-space chunks: "
                         "split the phase space into more chunks or reduce jobs_per_field")
    runs = []
    for fi, vf in enumerate(fields):
        for j in range(jobs_per_field):
            chunk = chunks[(fi * jobs_per_field + j) % len(chunks)] if vcfg.get("rotate_chunks", True) else chunks[j]
            seed = base_seed + 1000 * fi + j
            # recycling: TOPAS wraps to the start of the chunk at end of file, so
            # histories = recycle x chunk histories reads the chunk `recycle` times
            runs.append(make_water_job(cfg, vcfg, vf, chunk, j, seed, Path(work_dir) / vf.name / f"job{j:03d}",
                                       histories=int(chunk["histories"]) * recycle))
    return runs
