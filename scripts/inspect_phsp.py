#!/usr/bin/env python3
"""Print what is inside an IAEA phase-space header (Step 1).

    python scripts/inspect_phsp.py /data/phsp/TrueBeam_v2_6X_00
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mcverify.phsp import read_header  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("stems", nargs="+", help="phase-space path(s) WITHOUT extension")
for stem in ap.parse_args().stems:
    print(read_header(stem).summary(), "\n")
