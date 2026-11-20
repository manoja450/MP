#!/usr/bin/env python3
"""Print a summary of a DICOM RT Plan (beams, MU, control points, MLC).

    python scripts/inspect_rtplan.py RTPLAN.dcm
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mcverify.dicom.rtplan import read_rtplan, summarize  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("rtplan")
print(summarize(read_rtplan(ap.parse_args().rtplan)))
