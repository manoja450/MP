"""Create small, valid-enough RTPLAN files for tests and phantom checks.

These are NOT clinical plans. They let us test the plan reader and, later, run
simple end-to-end checks (open field, static MLC field, a simple arc) before
touching patient data.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import ExplicitVRLittleEndian, generate_uid

RT_PLAN_STORAGE = "1.2.840.10008.5.1.4.1.1.481.5"

M120_BOUNDARIES_MM = np.concatenate([
    np.arange(-200, -100, 10), np.arange(-100, 100, 5), np.arange(100, 201, 10)]).astype(float)


def _device(kind: str, values) -> Dataset:
    d = Dataset()
    d.RTBeamLimitingDeviceType = kind
    d.LeafJawPositions = [float(v) for v in values]
    return d


def make_arc_plan(path: str | Path, mu: float = 400.0, n_cp: int = 9, start: float = 181.0,
                  stop: float = 179.0, direction: str = "CW", field_mm: float = 30.0,
                  leaf_shift_mm_per_cp: float = 0.0, energy: float = 6.0, fff: bool = False) -> Path:
    """One VMAT-like arc with a square MLC aperture that can drift per control point."""
    meta = FileMetaDataset()
    meta.MediaStorageSOPClassUID = RT_PLAN_STORAGE
    meta.MediaStorageSOPInstanceUID = generate_uid()
    meta.TransferSyntaxUID = ExplicitVRLittleEndian
    ds = FileDataset(str(path), {}, file_meta=meta, preamble=b"\0" * 128)
    ds.SOPClassUID = RT_PLAN_STORAGE
    ds.SOPInstanceUID = meta.MediaStorageSOPInstanceUID
    ds.Modality = "RTPLAN"
    ds.PatientID = "SYNTHETIC"
    ds.PatientName = "Synthetic^Phantom"
    ds.RTPlanLabel = "SYNTH_ARC"
    ds.StudyInstanceUID = generate_uid()
    ds.SeriesInstanceUID = generate_uid()

    fg = Dataset()
    fg.FractionGroupNumber = 1
    fg.NumberOfFractionsPlanned = 1
    rb = Dataset()
    rb.ReferencedBeamNumber = 1
    rb.BeamMeterset = float(mu)
    fg.ReferencedBeamSequence = Sequence([rb])
    ds.FractionGroupSequence = Sequence([fg])

    b = Dataset()
    b.BeamNumber = 1
    b.BeamName = "ARC1"
    b.BeamType = "DYNAMIC"
    b.RadiationType = "PHOTON"
    b.TreatmentDeliveryType = "TREATMENT"
    b.SourceAxisDistance = 1000.0
    b.NumberOfControlPoints = n_cp
    b.FinalCumulativeMetersetWeight = 1.0
    devs = []
    for kind, n in (("ASYMX", 1), ("ASYMY", 1)):
        d = Dataset()
        d.RTBeamLimitingDeviceType = kind
        d.NumberOfLeafJawPairs = n
        devs.append(d)
    mlc = Dataset()
    mlc.RTBeamLimitingDeviceType = "MLCX"
    mlc.NumberOfLeafJawPairs = len(M120_BOUNDARIES_MM) - 1
    mlc.LeafPositionBoundaries = M120_BOUNDARIES_MM.tolist()
    devs.append(mlc)
    b.BeamLimitingDeviceSequence = Sequence(devs)
    if fff:
        pfm = Dataset()
        pfm.FluenceMode = "NON_STANDARD"
        pfm.FluenceModeID = "FFF"
        b.PrimaryFluenceModeSequence = Sequence([pfm])

    span = ((stop - start) % 360.0) if direction == "CW" else -((start - stop) % 360.0)
    gantry = [(start + span * k / (n_cp - 1)) % 360.0 for k in range(n_cp)]
    centres = 0.5 * (M120_BOUNDARIES_MM[:-1] + M120_BOUNDARIES_MM[1:])
    h = field_mm / 2.0
    cps = []
    for k in range(n_cp):
        cp = Dataset()
        cp.ControlPointIndex = k
        cp.GantryAngle = float(gantry[k])
        cp.CumulativeMetersetWeight = k / (n_cp - 1)
        shift = k * leaf_shift_mm_per_cp
        inside = np.abs(centres) < h
        neg = np.where(inside, -h + shift, shift)
        pos = np.where(inside, h + shift, shift)
        items = [_device("MLCX", np.concatenate([neg, pos]))]
        if k == 0:
            cp.NominalBeamEnergy = energy
            cp.GantryRotationDirection = direction
            cp.BeamLimitingDeviceAngle = 0.0
            cp.PatientSupportAngle = 0.0
            cp.IsocenterPosition = [0.0, 0.0, 0.0]
            items = [_device("ASYMX", [-h - 10, h + 10]), _device("ASYMY", [-h, h])] + items
        else:
            cp.GantryRotationDirection = direction if k < n_cp - 1 else "NONE"
        cp.BeamLimitingDevicePositionSequence = Sequence(items)
        cps.append(cp)
    b.ControlPointSequence = Sequence(cps)
    ds.BeamSequence = Sequence([b])
    ds.save_as(str(path), enforce_file_format=True)
    return Path(path)
