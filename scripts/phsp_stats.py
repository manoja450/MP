#!/usr/bin/env python3
"""Read the first N particles of a phase space and print coordinate statistics.

Use it to confirm where the scoring plane is (Z) and which way particles move
before running TOPAS:
    python scripts/phsp_stats.py phsp/TrueBeam_6X/Varian_TrueBeam6MV_01
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mcverify.phsp import describe_sample, read_header, sample_particles  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("stem", help="phase-space path WITHOUT extension")
ap.add_argument("-n", type=int, default=1_000_000, help="particles to read (default 1e6)")
ap.add_argument("--header", action="store_true", help="also print the full header text")
a = ap.parse_args()
if a.header:
    from mcverify.phsp import find_header
    print(find_header(a.stem).read_text(errors="replace"))
print(read_header(a.stem).summary())
print(describe_sample(sample_particles(a.stem, a.n)))
