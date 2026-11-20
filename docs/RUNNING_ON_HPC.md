# Running on the HPC cluster

## 1. Install OpenTOPAS

`scripts/install_opentopas.sh` follows the official guide (`OpenTOPAS_quickStart_Debian.md`, OpenTOPAS **v4.3.1** with **Geant4 11.4.3**), adapted for a cluster: no `sudo`, no graphics (Qt), and everything installed under `$HOME/Applications`.

Load a C++17 compiler and CMake 3.16 or newer. Module names differ per cluster; `module avail gcc` and `module avail cmake` list what exists.

```bash
module load gcc cmake
```

Download the sources and Geant4 data, about 5 GB. Run this on the **login node**, which has internet:

```bash
bash scripts/install_opentopas.sh download
```

Build, which takes 1–2 hours. You can run it on the login node, but many clusters prefer a short batch or interactive job for this:

```bash
JOBS=8 bash scripts/install_opentopas.sh build
```

Test with a 10,000-photon run in water:

```bash
bash scripts/install_opentopas.sh test
```

The build writes `$HOME/Applications/topas_env.sh`. `slurm/run_topas_array.sbatch` sources this file automatically. To use `topas` in your own terminal:

```bash
source $HOME/Applications/topas_env.sh
```

## 2. Python environment

```bash
python -m venv ~/venvs/mcv
```

```bash
source ~/venvs/mcv/bin/activate
```

```bash
pip install -r requirements.txt
```

```bash
python -m pytest
```

All tests should pass before you generate any jobs.

## 3. Phase-space files

1. Download the TrueBeam 6X (and later 6FFF) phase spaces from Varian. You need a myVarian account; check the licence terms, which may limit redistribution.
2. Inspect the header:

   ```bash
   python scripts/inspect_phsp.py /data/phsp/TrueBeam_v2_6X_00
   ```

   Check that `Z stored` is False and that `Constant Z` is the plane just above the jaws.

3. Split every file into chunks on fast scratch storage:

   ```bash
   python scripts/split_phsp.py /data/phsp/TrueBeam_v2_6X_0* --out /scratch/$USER/phsp/TrueBeam_6X/chunks --chunks-per-file 10
   ```

4. Put the chunk directory into `config/truebeam_6x_millennium.yaml` under `phase_space.chunk_dir`.

## 4. Smoke test (do this first)

Run a few short jobs to confirm that TOPAS accepts every generated parameter and that the geometry is oriented correctly.

```bash
python scripts/make_water_jobs.py --machine config/truebeam_6x_millennium.yaml --validation config/validation_water.yaml --work /scratch/$USER/mcv/smoke --only jaw10cm asymY asymX mlc10cm --jobs-per-field 1 --histories 200000
```

```bash
sbatch --array=1-4 slurm/run_topas_array.sbatch /scratch/$USER/mcv/smoke/jobs.txt
```

Check each job's `topas.log` for:

- **Unknown or unused parameter errors.** TOPAS stops on these. Report the exact message, and the generator will be fixed.
- **Overlap warnings** (`G4Exception ... GeomVol1002`).
- **The phase-space read.** Confirm the reported number of histories and that there are no "unsupported particle" aborts.
- **CPU time per history.** Use it to set `jobs_per_field` and the walltime.

Then run the analysis without reference data:

```bash
python scripts/analyze_water.py --work /scratch/$USER/mcv/smoke --out results/smoke
```

In `results/smoke/REPORT.md`, check two things:

- **Orientation checks** (`asymY`, `asymX`) must be PASS. If one fails, the field or phase space is mirrored. Flip `invert_x_axis` or `invert_y_axis` in the machine YAML, or check the jaw `RotZ`.
- **Field size of jaw10cm** should be close to 10 cm at 10 cm depth projected (about 11 cm at SSD 100 cm, 10 cm depth).

Optional: view the geometry once on a node with graphics. Set `b:Gr/Enable = "True"` and a small history count in one `run.txt`.

## 5. Full validation campaign

```bash
python scripts/make_water_jobs.py --machine config/truebeam_6x_millennium.yaml --validation config/validation_water.yaml --work /scratch/$USER/mcv/water_6X
```

```bash
sbatch --array=1-$(wc -l < /scratch/$USER/mcv/water_6X/jobs.txt) slurm/run_topas_array.sbatch /scratch/$USER/mcv/water_6X/jobs.txt
```

The script skips jobs that already have `done.flag`, so you can resubmit the same array after failures.

```bash
python scripts/analyze_water.py --work /scratch/$USER/mcv/water_6X --ref reference_data/6X --out results/water_6X
```

**Statistics target:** a median relative uncertainty of 0.5% or less for the PDD in the high-dose region. The report prints it. Increase `jobs_per_field` (and the number of chunks) until you reach it.
