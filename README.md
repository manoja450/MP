# MP: Monte Carlo independent dose verification for lung SBRT VMAT

This repository holds the code for a two-paper project:

1. **Paper 1 (methods):** an open-source TOPAS/Geant4 framework for independent Monte Carlo dose verification of VMAT plans, validated in water and applied to lung SBRT.
2. **Paper 2 (clinical result):** how sensitive Monte Carlo secondary checks are to delivery errors, comparing gamma analysis with DVH-based metrics.

The Python package `mcverify` does three jobs:
- It generates the TOPAS parameter files.
- It prepares the Varian phase-space files.
- It analyses the TOPAS output: PDDs, profiles, output factors, gamma, and later DVHs.

The simulations themselves run with **OpenTOPAS** on your HPC cluster.

> **Status:** Phase 1 (beam model and water validation, Steps 1–6) is implemented, and the RTPLAN reader is in place. Patient recalculation (Steps 7–10) and the error study (Paper 2) come next; see [`docs/ROADMAP.md`](docs/ROADMAP.md).
> The Python code is unit-tested, but **the generated TOPAS files have not yet been run in TOPAS.** The first task on the cluster is the smoke test in [`docs/RUNNING_ON_HPC.md`](docs/RUNNING_ON_HPC.md).

## Repository layout

```
config/                     machine model + validation campaign (YAML)
mcverify/
  phsp.py                   IAEA phase-space header reader and splitter
  config.py                 YAML loading/validation, MLC leaf widths
  topas/params.py           TOPAS parameter-file writer
  topas/linac.py            world, IEC frames, physics, phase-space source, jaws, MLC
  topas/water.py            water-tank validation jobs (PDD, profiles, output factors)
  analysis/topas_io.py      read TOPAS CSV/binary scorer output
  analysis/combine.py       combine jobs -> dose/history + batch uncertainty
  analysis/curves.py        PDD/profile metrics, 1-D gamma
  analysis/reference.py     load measured/published reference data
  analysis/water_report.py  full validation report (tables, figures, REPORT.md)
  dicom/rtplan.py           DICOM RT Plan reader (beams, control points, MLC, MU)
  dicom/synthetic.py        synthetic RTPLANs for tests and phantom checks
scripts/                    command-line entry points (see below)
slurm/run_topas_array.sbatch  one TOPAS job per SLURM array task
tests/                      pytest suite (synthetic data with known answers)
docs/                       roadmap, HPC guide, TG-268 reporting table
reference_data/             your reference beam data (not committed, see README there)
```

## Quick start

Python 3.10 or newer is required.

```bash
pip install -r requirements.txt
```

```bash
python -m pytest
```

On the cluster, after OpenTOPAS is installed (see `docs/RUNNING_ON_HPC.md`), the workflow is:

```bash
python scripts/inspect_phsp.py /data/phsp/TrueBeam_v2_6X_00
```

```bash
python scripts/split_phsp.py /data/phsp/TrueBeam_v2_6X_0* --out /scratch/$USER/phsp/TrueBeam_6X/chunks --chunks-per-file 10
```

```bash
python scripts/make_water_jobs.py --machine config/truebeam_6x_millennium.yaml --validation config/validation_water.yaml --work /scratch/$USER/mcv/water_6X
```

```bash
sbatch --array=1-N slurm/run_topas_array.sbatch /scratch/$USER/mcv/water_6X/jobs.txt
```

```bash
python scripts/analyze_water.py --work /scratch/$USER/mcv/water_6X --ref reference_data/6X --out results/water_6X
```

## Design decisions

These choices come from reading the OpenTOPAS source code.

| Decision | Reason |
|---|---|
| Phase spaces are **split into chunks**, one chunk per job | Every TOPAS job reads a phase space from the start of the file, so jobs that share a file would simulate identical particles. |
| `PhaseSpaceMultipleUse = 0` with an explicit `NumberOfHistoriesInRun` | If MultipleUse is greater than 0, TOPAS ignores the requested history count and runs the whole file. |
| `LimitedAssumePhotonIsNewHistory` and `LimitedAssumeFirstParticleIsNewHistory` | Varian IAEA files have no new-history flags (see OpenTOPAS `examples/PhaseSpace/ReadVarian.txt`). |
| The phase-space plane group sits at `SAD - Constant Z` and `PhaseSpaceInvertZAxis = True` | When Z is not stored, TOPAS puts particles at z = 0 of the source component. IAEA +z points downstream, while TOPAS IEC +z points to the source. |
| Jaws use `TsJaws`; the MLC uses `TsDivergingMLC` (focused leaf ends) | These are the built-in OpenTOPAS linac components (see `examples/MVLinac`). |
| MLC DLG and transmission are matched with an effective leaf-end offset and alloy density | `TsDivergingMLC` has no rounded leaf tips and no tongue-and-groove. We record this as a limitation in the paper; a detailed MLC could later be added as a TOPAS C++ extension. |
| Scorers sit in separate **parallel worlds** | This lets the PDD column, profile strips and output-factor voxel overlap. Dose uses the water material of the mass world. |
| Uncertainty is computed with the **batch method** (spread between independent jobs) | This is a standard, transparent method and is easy to report under TG-268. |

## Citation and licence

There is no licence file yet. Add one before making the repository public (MIT or BSD-3 are common for research code).
