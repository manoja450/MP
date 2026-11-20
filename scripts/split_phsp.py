#!/usr/bin/env python3
"""Split IAEA phase-space files into independent chunks, one chunk per job (Step 1).

Every TOPAS job reads its phase space from the beginning, so jobs that share a
file would simulate identical particles. Splitting gives each job its own data.

    python scripts/split_phsp.py /data/phsp/TrueBeam_v2_6X_0{0..7} \\
        --out /scratch/$USER/phsp/TrueBeam_6X/chunks --chunks-per-file 10
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mcverify.phsp import split_phase_space  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("stems", nargs="+", help="phase-space path(s) WITHOUT extension (.IAEAheader/.IAEAphsp)")
ap.add_argument("--out", required=True, help="output directory for the chunks")
ap.add_argument("--chunks-per-file", type=int, default=10)
a = ap.parse_args()
total = 0
for stem in a.stems:
    chunks = split_phase_space(stem, a.out, a.chunks_per_file)
    total += len(chunks)
    print(f"{stem}: {len(chunks)} chunks, ~{chunks[0]['histories']:.3g} histories each "
          f"(photons {chunks[0]['photons']:.3g}, electrons {chunks[0]['electrons']:.3g})")
print(f"Wrote {total} chunks to {a.out}")
