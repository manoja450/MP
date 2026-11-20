#!/usr/bin/env python3
"""Generate the water-tank validation TOPAS jobs (Steps 2-5).

    python scripts/make_water_jobs.py --machine config/truebeam_6x_millennium.yaml \\
        --validation config/validation_water.yaml --work /scratch/$USER/mcv/water_6X

Writes WORK/<field>/jobNNN/{run.txt,meta.json} and WORK/jobs.txt (one run file
per line) for the SLURM array script.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mcverify.config import load_machine, load_yaml  # noqa: E402
from mcverify.phsp import load_chunk_manifests  # noqa: E402
from mcverify.topas.water import fields_from_config, make_water_campaign  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--machine", required=True)
ap.add_argument("--validation", required=True)
ap.add_argument("--work", required=True, help="output directory for jobs and results")
ap.add_argument("--only", nargs="*", help="only these field names (e.g. jaw10cm asymY)")
ap.add_argument("--jobs-per-field", help="override jobs_per_field (a number, or 'all' = one job per chunk)")
ap.add_argument("--recycle", type=int, help="override phsp_recycle (times each job reads its chunk)")
ap.add_argument("--histories", type=int, help="override histories per job (default: one pass through the chunk)")
ap.add_argument("--base-seed", type=int, default=1000)
a = ap.parse_args()

cfg = load_machine(a.machine)
vcfg = load_yaml(a.validation)
if a.jobs_per_field:
    vcfg["jobs_per_field"] = a.jobs_per_field if a.jobs_per_field == "all" else int(a.jobs_per_field)
if a.recycle:
    vcfg["phsp_recycle"] = a.recycle
if a.only:
    known = {f.name for f in fields_from_config(vcfg)}
    unknown = set(a.only) - known
    if unknown:
        sys.exit(f"unknown field(s) {sorted(unknown)}; known: {sorted(known)}")
    vcfg["square_fields_cm"] = [s for s in vcfg.get("square_fields_cm", [])
                                if f"jaw{s:g}".replace(".", "p") + "cm" in a.only]
    vcfg["extra_fields"] = [f for f in vcfg.get("extra_fields", []) if f["name"] in a.only]
chunks = load_chunk_manifests(cfg["phase_space"]["chunk_dir"])
if a.histories:
    for c in chunks:
        c["histories"] = a.histories
    vcfg["phsp_recycle"] = 1   # --histories is the exact count per job
work = Path(a.work)
runs = make_water_campaign(cfg, vcfg, chunks, work, base_seed=a.base_seed)
(work / "jobs.txt").write_text("\n".join(str(r) for r in runs) + "\n")
print(f"{len(runs)} jobs written under {work}")
print(f"Submit with:  sbatch --array=1-{len(runs)} slurm/run_topas_array.sbatch {work / 'jobs.txt'}")
