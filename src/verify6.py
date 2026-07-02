"""Exactly re-verify every witness pair reported by analyze6.py.

Digest keys can only merge classes spuriously (never split), so every
claimed witness is recomputed from scratch here: full invariants via the
Fraction-based (slow, independent) pipeline and structural D-WL / DWL+
comparison. Any pair failing re-verification is reported and must be
discarded from the paper.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.digraph import INF, distance_matrix, parse_digraph6
from src.dwl import dwl_equivalent, dwl_plus_equivalent
from src.maghom import magnitude_homology_ranks
from src.pathhom import path_homology_ranks
from src.reach import reachability_homology_ranks

ROOT = Path(__file__).resolve().parents[1]


def invariants(g6):
    n, arcs = parse_digraph6(g6)
    D = distance_matrix(n, arcs)
    return {
        "mh": magnitude_homology_ranks(n, arcs, kmax=3),
        "ph": path_homology_ranks(n, arcs, pmax=2),
        "rh": reachability_homology_ranks(n, arcs),
        "dist": sorted(
            (10**6 if D[x][y] == INF else int(D[x][y]))
            for x in range(n)
            for y in range(n)
            if x != y
        ),
    }


def equal_under(key, a, b, ia, ib):
    ga, gb = parse_digraph6(a), parse_digraph6(b)
    if key == "dwl":
        return dwl_equivalent(ga, gb)
    if key == "dwlp":
        return dwl_plus_equivalent(ga, gb)
    if key == "dwl_dist":
        return dwl_equivalent(ga, gb) and ia["dist"] == ib["dist"]
    if key == "dwlp_dist":
        return dwl_plus_equivalent(ga, gb) and ia["dist"] == ib["dist"]
    if key == "mhph":
        return ia["mh"] == ib["mh"] and ia["ph"] == ib["ph"]
    if key == "hom":
        return all(ia[k] == ib[k] for k in ("mh", "ph", "rh"))
    return ia[key] == ib[key]


def main():
    data = json.loads((ROOT / "results" / "analysis_n6.json").read_text())
    ok = bad = 0
    for rel, pairs in data["witnesses"].items():
        ki, kj = rel.split("|")
        for a, b in pairs:
            ia, ib = invariants(a), invariants(b)
            same_i = equal_under(ki, a, b, ia, ib)
            diff_j = not equal_under(kj, a, b, ia, ib)
            if same_i and diff_j:
                ok += 1
                print(f"  ok  {ki}-equal,{kj}-diff: {a} vs {b}")
            else:
                bad += 1
                print(f"FAIL  {ki}-equal,{kj}-diff: {a} vs {b} "
                      f"(same_{ki}={same_i}, diff_{kj}={diff_j})")
    print(f"\n{ok} verified, {bad} failed")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
