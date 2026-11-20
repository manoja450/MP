"""TOPAS blocks for the linac: world, IEC coordinate systems, physics, phase-space
source, jaws and MLC.

Coordinate systems follow the OpenTOPAS ``examples/MVLinac`` convention:
  IEC_F  fixed system, origin at isocentre, +z towards the source
  IEC_G  gantry      (rotates about y with the gantry angle)
  IEC_B  beam-limiting device (rotates about z with the collimator angle)
  IEC_S  patient support (rotates about z with the couch angle)
The X-ray target sits at z = +SAD in IEC_G. Jaw and MLC positions are given as
field settings projected to the isocentre plane (TsJaws / TsDivergingMLC).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

from ..config import leaf_widths_cm
from .params import Expr, ParamWriter

JAW_ROTZ_FOR_TRAVEL = {"X": 0.0, "Y": -90.0}   # from OpenTOPAS examples/MVLinac/Jaws.txt
MLC_ROTZ_FOR_TRAVEL = {"X": 0.0, "Y": -90.0}   # from OpenTOPAS examples/MVLinac/MLC.txt


@dataclass
class FieldSetting:
    """One beam aperture, all positions in cm projected to the isocentre plane.

    x1, y1 are the negative-side jaw positions (normally negative numbers).
    mlc_neg / mlc_pos: one value per leaf pair, DICOM order (from -Y to +Y).
    If they are None the leaves are retracted (jaw-defined field).
    """

    x1: float
    x2: float
    y1: float
    y2: float
    mlc_neg: Sequence[float] | None = None
    mlc_pos: Sequence[float] | None = None
    gantry_deg: float = 0.0
    collimator_deg: float = 0.0
    couch_deg: float = 0.0

    @classmethod
    def square(cls, size_cm: float, **kw) -> "FieldSetting":
        h = size_cm / 2.0
        return cls(-h, h, -h, h, **kw)

    def __post_init__(self):
        if self.x1 > self.x2 or self.y1 > self.y2:
            raise ValueError("jaw positions must satisfy x1 <= x2 and y1 <= y2")


@dataclass
class SourceSpec:
    """Which phase-space chunk to read and how many histories to run."""

    stem: str
    histories: int
    z_stored: bool = False
    z_constant_cm: float | None = None
    extra: dict = field(default_factory=dict)


def world_and_iec(w: ParamWriter, cfg: dict, fs: FieldSetting) -> None:
    w.section("World and IEC coordinate systems")
    w.s("Ge/World/Material", cfg["physics"].get("world_material", "G4_AIR"))
    for ax in ("HLX", "HLY", "HLZ"):
        w.d(f"Ge/World/{ax}", 1.5, "m")
    w.b("Ge/World/Invisible", True)
    w.d("Ge/SAD", cfg["sad_cm"], "cm")
    parents = {"IEC_F": "World", "IEC_G": "IEC_F", "IEC_B": "IEC_G", "IEC_S": "IEC_F"}
    for name, parent in parents.items():
        w.s(f"Ge/{name}/Parent", parent)
        w.s(f"Ge/{name}/Type", "Group")
        for t in ("TransX", "TransY", "TransZ"):
            w.d(f"Ge/{name}/{t}", 0.0, "cm")
    # IEC angles are applied in the opposite sense of TOPAS rotations (see MVLinac/Field.txt)
    w.d("Ge/IEC_G/RotY", 0.0 - fs.gantry_deg + 0.0, "deg")
    w.d("Ge/IEC_B/RotZ", 0.0 - fs.collimator_deg + 0.0, "deg")
    w.d("Ge/IEC_S/RotZ", 0.0 - fs.couch_deg + 0.0, "deg")


def physics(w: ParamWriter, cfg: dict) -> None:
    w.section("Physics")
    w.sv("Ph/Default/Modules", list(cfg["physics"]["modules"]))
    w.d("Ph/Default/CutForAllParticles", cfg["physics"]["production_cut_mm"], "mm")


def phase_space_source(w: ParamWriter, cfg: dict, src: SourceSpec, name: str = "Linac") -> None:
    ps = cfg["phase_space"]
    w.section("Phase-space source (particles just above the jaws)")
    plane = ps.get("plane_z_from_target_cm")
    w.s("Ge/PhspPlane/Parent", "IEC_G")
    w.s("Ge/PhspPlane/Type", "Group")
    if src.z_stored:
        # z in the file is measured from the target; after inversion z_world = SAD - z
        w.comment("Z stored per particle: plane group sits at the target (z = SAD)")
        w.d("Ge/PhspPlane/TransZ", cfg["sad_cm"], "cm")
    else:
        zc = plane if plane is not None else src.z_constant_cm
        if zc is None:
            raise ValueError("phase space has no stored Z: set phase_space.plane_z_from_target_cm "
                             "or make sure the header contains 'Constant Z'")
        w.comment(f"Z not stored: all particles lie on the plane {zc} cm below the target")
        w.d("Ge/PhspPlane/TransZ", cfg["sad_cm"] - float(zc), "cm")
    p = f"So/{name}"
    w.s(f"{p}/Type", "PhaseSpace")
    w.s(f"{p}/Component", "PhspPlane")
    w.s(f"{p}/PhaseSpaceFileName", src.stem)
    w.i(f"{p}/PhaseSpaceMultipleUse", 0)          # 0 => obey NumberOfHistoriesInRun
    w.i(f"{p}/NumberOfHistoriesInRun", src.histories)
    w.b(f"{p}/PhaseSpacePreCheck", False)
    w.b(f"{p}/PhaseSpaceIncludeEmptyHistories", bool(ps.get("include_empty_histories", False)))
    w.b(f"{p}/LimitedAssumePhotonIsNewHistory", bool(ps.get("limited_assume_photon_is_new_history", True)))
    w.b(f"{p}/LimitedAssumeFirstParticleIsNewHistory",
        bool(ps.get("limited_assume_first_particle_is_new_history", True)))
    w.b(f"{p}/IgnoreUnsupportedParticles", bool(ps.get("ignore_unsupported_particles", True)))
    for ax in ("X", "Y", "Z"):
        default = ax == "Z"
        w.b(f"{p}/PhaseSpaceInvert{ax}Axis", bool(ps.get(f"invert_{ax.lower()}_axis", default)))


def _one_jaw(w: ParamWriter, cfg: dict, key: str, lo: float, hi: float) -> None:
    j = cfg["jaws"][key]
    n = f"Jaw{j['travel_axis']}"
    w.comment(f"{key} jaw pair, travels along IEC {j['travel_axis']}")
    w.s(f"Ge/{n}/Parent", "IEC_B")
    w.s(f"Ge/{n}/Type", "TsJaws")
    w.s(f"Ge/{n}/Material", j.get("material", "G4_W"))
    w.d(f"Ge/{n}/LX", j["width_cm"], "cm")
    w.d(f"Ge/{n}/LY", j["length_cm"], "cm")
    w.d(f"Ge/{n}/LZ", j["thickness_cm"], "cm")
    w.d(f"Ge/{n}/SourceToUpstreamSurfaceDistance", j["source_to_upstream_cm"], "cm")
    w.d(f"Ge/{n}/SAD", Expr("Ge/SAD"), "cm")
    w.d(f"Ge/{n}/TransZ", cfg["sad_cm"] - j["source_to_upstream_cm"] - 0.5 * j["thickness_cm"], "cm")
    w.d(f"Ge/{n}/RotZ", JAW_ROTZ_FOR_TRAVEL[j["travel_axis"]], "deg")
    w.d(f"Ge/{n}/NegativeFieldSetting", lo, "cm", changeable=True)
    w.d(f"Ge/{n}/PositiveFieldSetting", hi, "cm", changeable=True)


def jaws(w: ParamWriter, cfg: dict, fs: FieldSetting) -> None:
    w.section("Jaws (TsJaws, focused, positions projected to isocentre)")
    for key in ("upper", "lower"):
        ax = cfg["jaws"][key]["travel_axis"]
        lo, hi = (fs.x1, fs.x2) if ax == "X" else (fs.y1, fs.y2)
        _one_jaw(w, cfg, key, lo, hi)


def mlc_material(w: ParamWriter, cfg: dict) -> str:
    mat = cfg["mlc"]["material"]
    if isinstance(mat, str):
        return mat
    n = mat["name"]
    w.sv(f"Ma/{n}/Components", list(mat["components"]))
    w.uv(f"Ma/{n}/Fractions", list(mat["fractions"]))
    w.d(f"Ma/{n}/Density", mat["density_g_cm3"], "g/cm3")
    w.s(f"Ma/{n}/DefaultColor", "grey")
    return n


def mlc_positions(cfg: dict, fs: FieldSetting) -> tuple[list[float], list[float]]:
    """Leaf openings (cm at isocentre) including the effective leaf-end offset."""
    widths = leaf_widths_cm(cfg)
    n = len(widths)
    max_open = float(cfg["mlc"].get("max_leaf_open_cm", 20.0))
    if fs.mlc_neg is None or fs.mlc_pos is None:
        neg, pos = [-max_open] * n, [max_open] * n
        return neg, pos
    if len(fs.mlc_neg) != n or len(fs.mlc_pos) != n:
        raise ValueError(f"expected {n} leaf positions per bank for {cfg['mlc']['model']}")
    off = float(cfg["mlc"].get("leaf_end_offset_mm", 0.0)) / 10.0
    neg, pos = [], []
    for a, b in zip(fs.mlc_neg, fs.mlc_pos):
        if a > b + 1e-9:
            raise ValueError(f"leaf pair overlaps: negative {a} > positive {b}")
        a2, b2 = a - off, b + off
        if a2 > b2:  # a negative offset must not make leaves cross
            mid = 0.5 * (a + b)
            a2 = b2 = mid
        neg.append(max(-max_open, a2))
        pos.append(min(max_open, b2))
    return neg, pos


def mlc(w: ParamWriter, cfg: dict, fs: FieldSetting) -> None:
    m = cfg["mlc"]
    w.section(f"MLC ({m['model']}, TsDivergingMLC, positions projected to isocentre)")
    mat = mlc_material(w, cfg)
    w.s("Ge/MLC/Parent", "IEC_B")
    w.s("Ge/MLC/Type", "TsDivergingMLC")
    w.s("Ge/MLC/Material", mat)
    w.d("Ge/MLC/SAD", Expr("Ge/SAD"), "cm")
    w.d("Ge/MLC/SourceToUpstreamSurfaceDistance", m["source_to_upstream_cm"], "cm")
    w.d("Ge/MLC/MaximumLeafOpen", m.get("max_leaf_open_cm", 20.0), "cm")
    w.d("Ge/MLC/Thickness", m["thickness_cm"], "cm")
    w.d("Ge/MLC/Length", m["length_cm"], "cm")
    w.s("Ge/MLC/LeafTravelAxis", "Xb")
    w.d("Ge/MLC/TransZ", cfg["sad_cm"] - m["source_to_upstream_cm"] - 0.5 * m["thickness_cm"], "cm")
    w.d("Ge/MLC/RotZ", MLC_ROTZ_FOR_TRAVEL[m["travel_axis"]], "deg")
    w.dv("Ge/MLC/LeafWidths", leaf_widths_cm(cfg), "cm")
    neg, pos = mlc_positions(cfg, fs)
    w.dv("Ge/MLC/NegativeFieldSetting", neg, "cm")
    w.dv("Ge/MLC/PositiveFieldSetting", pos, "cm")


def run_control(w: ParamWriter, seed: int, threads: int = 1, show_every: int = 1_000_000) -> None:
    w.section("Run control")
    w.i("Ts/Seed", seed)
    w.i("Ts/NumberOfThreads", threads)
    w.i("Ts/ShowHistoryCountAtInterval", show_every)
    w.b("Ts/ShowCPUTime", True)
    w.b("Ts/PauseBeforeQuit", False)
    w.b("Gr/Enable", False)
    w.b("Ts/UseQt", False)


def linac(cfg: dict, fs: FieldSetting, src: SourceSpec) -> ParamWriter:
    """Everything above the patient/phantom for one static field."""
    w = ParamWriter()
    world_and_iec(w, cfg, fs)
    physics(w, cfg)
    phase_space_source(w, cfg, src)
    jaws(w, cfg, fs)
    mlc(w, cfg, fs)
    return w
