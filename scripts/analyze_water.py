#!/usr/bin/env python3
"""Analyse a finished water-tank campaign (Steps 2-5): combine jobs, compute
PDD/profile metrics, output factors, MLC transmission, orientation checks, and
compare with reference data if given.

    python scripts/analyze_water.py --work /scratch/$USER/mcv/water_6X \\
        --ref reference_data/6X --out results/water_6X
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mcverify.analysis.water_report import analyse_campaign  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--work", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--ref", help="reference data directory (optional)")
ap.add_argument("--reference-field", default="jaw10cm")
ap.add_argument("--no-plots", action="store_true")
a = ap.parse_args()
analyse_campaign(a.work, a.out, a.ref, a.reference_field, make_plots=not a.no_plots)
print(f"Report written to {Path(a.out) / 'REPORT.md'}")
