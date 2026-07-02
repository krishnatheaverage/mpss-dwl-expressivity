"""Exhaustive sweep: all digraphs on n vertices (via nauty geng|directg).

For every digraph up to isomorphism computes:
  - D-WL invariant (directed color refinement, canonical histograms)
  - magnitude homology profile: exact ranks MH_{k,l}, k <= KMAX, all l
  - path homology dims H_0..H_PMAX (GLMY, over Q)
  - reachability homology dims RH_0..RH_3

Writes one JSON line per digraph to results/sweep_n{n}.jsonl, then an
analysis pass (analyze.py) mines equivalence classes and witness pairs.

Every rank is exact (Fraction Gaussian elimination); no floats anywhere.
"""

import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.digraph import parse_digraph6
from src.dwl import dwl_colors
from src.maghom import magnitude_homology_ranks
from src.pathhom import path_homology_ranks
from src.reach import reachability_homology_ranks

KMAX = 3  # magnitude homology degrees 0..3
PMAX = 2  # path homology H_0..H_2


def invariants(n, arcs):
    mh = magnitude_homology_ranks(n, arcs, kmax=KMAX)
    ph = path_homology_ranks(n, arcs, pmax=PMAX)
    rh = reachability_homology_ranks(n, arcs)
    dwl = dwl_colors(n, arcs)
    return {
        "mh": sorted([k, l, r] for (k, l), r in mh.items()),
        "ph": ph,
        "rh": rh,
        "dwl": repr(dwl),
    }


def main(n):
    out_path = Path(__file__).resolve().parents[1] / "results" / f"sweep_n{n}.jsonl"
    proc = subprocess.run(
        f"geng -q {n} | directg -q", shell=True, capture_output=True, text=True
    )
    lines = proc.stdout.strip().split("\n")
    print(f"n={n}: {len(lines)} digraphs")
    t0 = time.time()
    with open(out_path, "w") as f:
        for i, line in enumerate(lines):
            g = parse_digraph6(line)
            inv = invariants(*g)
            inv["g6"] = line.strip()
            f.write(json.dumps(inv) + "\n")
            if (i + 1) % 500 == 0:
                el = time.time() - t0
                print(
                    f"  {i+1}/{len(lines)}  ({el:.0f}s, {el/(i+1)*1000:.0f} ms/graph)",
                    flush=True,
                )
    print(f"done in {time.time()-t0:.0f}s -> {out_path}")


if __name__ == "__main__":
    main(int(sys.argv[1]))
