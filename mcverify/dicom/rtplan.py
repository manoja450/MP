"""Read a DICOM RT Plan into simple Python objects.

Only what the dose recalculation needs is extracted: for every treatment beam
and control point the gantry / collimator / couch angles, jaw positions, MLC
leaf positions and the MU delivered in each segment.

DICOM conventions used (IEC 61217 beam-limiting-device coordinates):
* positions are in mm at the isocentre plane
* ``ASYMX``/``X`` jaws: [X1, X2];  ``ASYMY``/``Y`` jaws: [Y1, Y2]
* ``MLCX``: first N values = bank on the -X side (X1), next N = +X side (X2);
  leaves are ordered along +Y, matching ``LeafPositionBoundaries``
* values that do not change are only stored in the first control point, so
  later control points inherit them
* ``CumulativeMetersetWeight`` runs from 0 to ``FinalCumulativeMetersetWeight``
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pydicom


@dataclass
class ControlPoint:
    index: int
    gantry_deg: float
    collimator_deg: float
    couch_deg: float
    jaw_x_mm: tuple[float, float]
    jaw_y_mm: tuple[float, float]
    mlc_neg_mm: np.ndarray | None
    mlc_pos_mm: np.ndarray | None
    cumulative_weight: float
    isocenter_mm: tuple[float, float, float]


@dataclass
class Beam:
    number: int
    name: str
    beam_type: str                  # STATIC or DYNAMIC
    radiation_type: str
    energy_mv: float | None
    fluence_mode: str               # STANDARD or NON_STANDARD (FFF)
    sad_mm: float
    mu: float
    final_weight: float
    leaf_boundaries_mm: np.ndarray | None
    control_points: list[ControlPoint] = field(default_factory=list)

    @property
    def n_leaf_pairs(self) -> int:
        return 0 if self.leaf_boundaries_mm is None else len(self.leaf_boundaries_mm) - 1

    def segment_mu(self) -> np.ndarray:
        """MU delivered between control point i and i+1 (length = n_cp - 1)."""
        w = np.array([cp.cumulative_weight for cp in self.control_points])
        if self.final_weight <= 0:
            raise ValueError(f"beam {self.name}: FinalCumulativeMetersetWeight is {self.final_weight}")
        return np.diff(w) / self.final_weight * self.mu

    def is_arc(self) -> bool:
        g = [cp.gantry_deg for cp in self.control_points]
        return len(set(np.round(g, 3))) > 1


@dataclass
class Plan:
    path: str
    label: str
    patient_id: str
    beams: list[Beam]

    @property
    def total_mu(self) -> float:
        return float(sum(b.mu for b in self.beams))


def _limits(cp, keys):
    for dev in getattr(cp, "BeamLimitingDevicePositionSequence", []) or []:
        if dev.RTBeamLimitingDeviceType in keys:
            return [float(v) for v in dev.LeafJawPositions]
    return None


def read_rtplan(path: str | Path, treatment_beams_only: bool = True) -> Plan:
    ds = pydicom.dcmread(str(path))
    if getattr(ds, "Modality", "") != "RTPLAN":
        raise ValueError(f"{path} is not an RTPLAN (Modality={getattr(ds, 'Modality', None)})")

    mu = {}
    for fg in getattr(ds, "FractionGroupSequence", []) or []:
        for rb in getattr(fg, "ReferencedBeamSequence", []) or []:
            if hasattr(rb, "BeamMeterset"):
                mu[int(rb.ReferencedBeamNumber)] = float(rb.BeamMeterset)

    beams = []
    for b in ds.BeamSequence:
        if treatment_beams_only and getattr(b, "TreatmentDeliveryType", "TREATMENT") != "TREATMENT":
            continue
        number = int(b.BeamNumber)
        boundaries = None
        for bld in getattr(b, "BeamLimitingDeviceSequence", []) or []:
            if bld.RTBeamLimitingDeviceType in ("MLCX", "MLCY") and hasattr(bld, "LeafPositionBoundaries"):
                if bld.RTBeamLimitingDeviceType == "MLCY":
                    raise NotImplementedError("MLCY (leaves moving along Y) is not supported")
                boundaries = np.array([float(v) for v in bld.LeafPositionBoundaries])
        energy = None
        fluence = "STANDARD"
        state = {"gantry": 0.0, "coll": 0.0, "couch": 0.0, "jx": (-200.0, 200.0), "jy": (-200.0, 200.0),
                 "mlc": None, "iso": (0.0, 0.0, 0.0)}
        cps = []
        for k, cp in enumerate(b.ControlPointSequence):
            if hasattr(cp, "NominalBeamEnergy") and energy is None:
                energy = float(cp.NominalBeamEnergy)
            if hasattr(cp, "GantryAngle"):
                state["gantry"] = float(cp.GantryAngle)
            if hasattr(cp, "BeamLimitingDeviceAngle"):
                state["coll"] = float(cp.BeamLimitingDeviceAngle)
            if hasattr(cp, "PatientSupportAngle"):
                state["couch"] = float(cp.PatientSupportAngle)
            if hasattr(cp, "IsocenterPosition"):
                state["iso"] = tuple(float(v) for v in cp.IsocenterPosition)
            jx = _limits(cp, ("ASYMX", "X"))
            jy = _limits(cp, ("ASYMY", "Y"))
            if jx:
                state["jx"] = (jx[0], jx[1])
            if jy:
                state["jy"] = (jy[0], jy[1])
            leaves = _limits(cp, ("MLCX",))
            if leaves is not None:
                n = len(leaves) // 2
                state["mlc"] = (np.array(leaves[:n]), np.array(leaves[n:]))
            neg, pos = state["mlc"] if state["mlc"] is not None else (None, None)
            cps.append(ControlPoint(k, state["gantry"], state["coll"], state["couch"], state["jx"], state["jy"],
                                    None if neg is None else neg.copy(), None if pos is None else pos.copy(),
                                    float(cp.CumulativeMetersetWeight), state["iso"]))
        for pfm in getattr(b, "PrimaryFluenceModeSequence", []) or []:
            fluence = str(getattr(pfm, "FluenceMode", "STANDARD"))
            if fluence == "NON_STANDARD":
                fluence = f"NON_STANDARD:{getattr(pfm, 'FluenceModeID', '')}"
        beams.append(Beam(
            number=number, name=str(getattr(b, "BeamName", number)),
            beam_type=str(getattr(b, "BeamType", "STATIC")),
            radiation_type=str(getattr(b, "RadiationType", "PHOTON")),
            energy_mv=energy, fluence_mode=fluence,
            sad_mm=float(getattr(b, "SourceAxisDistance", 1000.0)),
            mu=mu.get(number, float("nan")),
            final_weight=float(getattr(b, "FinalCumulativeMetersetWeight", cps[-1].cumulative_weight if cps else 1)),
            leaf_boundaries_mm=boundaries, control_points=cps))
    return Plan(str(path), str(getattr(ds, "RTPlanLabel", "")), str(getattr(ds, "PatientID", "")), beams)


def summarize(plan: Plan) -> str:
    lines = [f"Plan '{plan.label}'  beams={len(plan.beams)}  total MU={plan.total_mu:.1f}"]
    for b in plan.beams:
        g = [cp.gantry_deg for cp in b.control_points]
        lines.append(f"  beam {b.number} '{b.name}': {b.beam_type} {b.energy_mv} MV {b.fluence_mode} "
                     f"MU={b.mu:.1f} CPs={len(b.control_points)} leaf pairs={b.n_leaf_pairs} "
                     f"gantry {g[0]:.1f}->{g[-1]:.1f} {'(arc)' if b.is_arc() else ''}")
    return "\n".join(lines)
