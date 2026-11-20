#!/usr/bin/env python3
"""Quick look at raw TOPAS outputs of a campaign: total dose, peak value and
peak position for every scorer of every field. Use it when something looks wrong.

    python scripts/quicklook.py --work $HOME/MP/mcv_runs/smoke
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mcverify.analysis.topas_io import read_scorer, to_1d  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--work", required=True)
a = ap.parse_args()
for meta_p in sorted(Path(a.work).glob("*/job*/meta.json")):
    meta = json.loads(meta_p.read_text())
    fs = meta["field_setting"]
    print(f"\n== {meta_p.parent.parent.name}/{meta_p.parent.name}  histories={meta['histories']}  "
          f"jaws X[{fs['x1']:g},{fs['x2']:g}] Y[{fs['y1']:g},{fs['y2']:g}]  "
          f"MLC={'yes' if fs.get('mlc_neg') else 'retracted'}")
    for name, sc in meta["scorers"].items():
        f = meta_p.parent / sc["file"]
        if not f.exists():
            print(f"   {name:14s} MISSING {f.name}")
            continue
        d = to_1d(read_scorer(f), sc["axis"])
        pos = sc.get("positions_cm")
        if pos:
            pos = np.asarray(pos)
            k = int(np.argmax(d))
            nz = np.count_nonzero(d)
            centroid = float(np.sum(pos * d) / d.sum()) if d.sum() > 0 else float("nan")
            print(f"   {name:14s} sum={d.sum():.3e} max={d.max():.3e} at {pos[k]:+7.2f} cm  "
                  f"centroid={centroid:+7.2f} cm  nonzero bins={nz}/{d.size}")
        else:
            print(f"   {name:14s} value={d[0]:.3e}")
