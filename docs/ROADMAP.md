# Roadmap: Paper 1 (methods and validation)

Steps follow the project plan and the beginner's guide (Part E). ✅ means implemented and unit-tested; 🔲 means to do.

## Phase 1: Beam model in water (Aim 1)

| Step | What | Code | Status |
|---|---|---|---|
| 1 | Install OpenTOPAS; inspect and split phase spaces | `scripts/inspect_phsp.py`, `scripts/split_phsp.py` | ✅ code, 🔲 run on HPC |
| 2 | PDDs and profiles in water for jaw-defined fields 1–30 cm | `scripts/make_water_jobs.py`, `scripts/analyze_water.py` | ✅ code, 🔲 run |
| 3 | Jaws and MLC geometry; orientation test (`asymX`, `asymY`), MLC-defined fields | `config/*.yaml`, `topas/linac.py` | ✅ code, 🔲 replace the VERIFY geometry values |
| 4 | MLC tuning: transmission (`mlcClosed10cm`), field edges (`mlc10cm`, `mlc3cm`) → `leaf_end_offset_mm`, MLC density | config + rerun | 🔲 (sweeping-gap DLG comes with Step 7 time features) |
| 5 | Small-field output factors (1, 2, 3, 4 cm) vs TRS-483 and published data; source spot size sensitivity | `output_factors.csv` in the report | ✅ code, 🔲 run |
| 6 | Absolute calibration: dose/history at reference → histories per MU; backscatter correction | 🔲 `analysis/calibration.py` | 🔲 |

**Acceptance criteria:**
- PDD and profiles: 1%/1 mm gamma ≥ 95%, and ≤ 1% difference beyond d<sub>max</sub>.
- Output factors: ≤ 1% for fields ≥ 4×4 cm, ≤ 2% at 2×2, ≤ 3% at 1×1.

### Before trusting any result
- [ ] Replace every `VERIFY` value in `config/truebeam_6x_millennium.yaml` with the geometry from Varian's Monte Carlo data package.
- [ ] Smoke test passes (see `RUNNING_ON_HPC.md`).
- [ ] The orientation checks pass.
- [ ] Reference data is in place (`reference_data/README.md`).
- [ ] Production cut 0.05 mm vs 0.1 mm compared for 1×1 and 2×2 fields.

## Phase 2: Patient plans (Aim 2)

| Step | What | Status |
|---|---|---|
| 7 | RTPLAN → TOPAS: per-control-point gantry, collimator, jaws, MLC as `Step` time features; histories per control point ∝ MU (sequential time mode) | ✅ RTPLAN reader; 🔲 TOPAS generator |
| 7b | Check: open field → static MLC → one arc in a phantom vs the TPS | 🔲 (synthetic plans: `dicom/synthetic.py`) |
| 8 | CT → `TsDicomPatient` with Schneider HU conversion, couch, isocentre | 🔲 |
| 9 | Patient runs on SLURM; merge; export RTDOSE (dose-to-medium and dose-to-water) | 🔲 |
| 10 | 3D gamma (pymedphys) and DVH comparison vs the TPS; lung vs H&N statistics | 🔲 |

## Phase 3: Paper 2 (error injection)

Plan modifications (MLC, MU, gantry, collimator, CT calibration, beam model), re-simulation, and ROC analysis. Not started.
