"""Streaming check of the MH_{2,2} closed formula over the n=6 sweep.

For every row: recompute sum_v indeg(v)*outdeg(v) - #{(x,z): d(x,z)=2}
from the digraph6 string and compare with the stored MH_{2,2} rank.
"""

import glob
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.digraph import parse_digraph6
from src.mh2 import mh22_formula

ROOT = Path(__file__).resolve().parents[1]

bad = total = 0
for path in sorted(glob.glob(str(ROOT / "results" / "sweep_n6_shard*.jsonl"))):
    with open(path) as f:
        for line in f:
            r = json.loads(line)
            stored = next((x[2] for x in r["mh"] if x[0] == 2 and x[1] == 2), 0)
            n, arcs = parse_digraph6(r["g6"])
            if mh22_formula(n, arcs) != stored:
                bad += 1
                print("mismatch:", r["g6"])
            total += 1
print(f"checked {total} digraphs (n=6): {bad} mismatches")
sys.exit(1 if bad else 0)
