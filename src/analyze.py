"""Mine sweep results for equivalence classes and witness pairs.

Questions answered per n:
  Q1  class counts: how many equivalence classes does each invariant
      induce on the set of iso-classes of digraphs?
  Q2  cross-refinement: for invariants I, J, how many I-classes split
      under J?  (I-equal but J-different = witness that J sees more)
  Q3  explicit witnesses:
      A: D-WL-equal, MH-different       (homology exceeds D-WL)
      B: D-WL-equal, PH-different
      C: MH-equal,  D-WL-different      (D-WL exceeds MH ranks)
      D: MH-equal,  PH-different        (E1 ranks do not determine E2)
      E: MH+PH+RH all equal, D-WL-different
      F: D-WL-equal and MH1-equal, MH{2,3}-different
         (magnitude homology beyond distance counts)
  MH profile = exact ranks MH_{k,l}, k <= 3 (all l);
  MH1 = the k=1 row only (equivalent to directed-distance multiset).
"""

import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.digraph import parse_digraph6


def arcs_str(g6):
    n, arcs = parse_digraph6(g6)
    return f"n={n} arcs={sorted(arcs)}"


def dist_multiset_key(g6):
    """Multiset of off-diagonal directed distances (inf included)."""
    from src.digraph import INF, distance_matrix

    n, arcs = parse_digraph6(g6)
    D = distance_matrix(n, arcs)
    vals = sorted(
        (10**6 if D[x][y] == INF else int(D[x][y]))
        for x in range(n)
        for y in range(n)
        if x != y
    )
    return json.dumps(vals)


def load(n):
    path = Path(__file__).resolve().parents[1] / "results" / f"sweep_n{n}.jsonl"
    rows = []
    with open(path) as f:
        for line in f:
            r = json.loads(line)
            r["mh_key"] = json.dumps(r["mh"])
            r["mh1_key"] = json.dumps([x for x in r["mh"] if x[0] == 1])
            r["ph_key"] = json.dumps(r["ph"])
            r["rh_key"] = json.dumps(r["rh"])
            r["dwl_key"] = r["dwl"]
            r["dist_key"] = dist_multiset_key(r["g6"])
            r["all_hom_key"] = r["mh_key"] + "|" + r["ph_key"] + "|" + r["rh_key"]
            rows.append(r)
    return rows


def classes(rows, key):
    d = defaultdict(list)
    for r in rows:
        d[r[key]].append(r)
    return d


def split_witnesses(rows, key_i, key_j, max_report=3):
    """I-classes that split under J; returns (count, example pairs)."""
    d = classes(rows, key_i)
    n_split = 0
    examples = []
    for group in d.values():
        if len(group) < 2:
            continue
        sub = defaultdict(list)
        for r in group:
            sub[r[key_j]].append(r)
        if len(sub) > 1:
            n_split += 1
            if len(examples) < max_report:
                reps = [v[0] for v in sub.values()][:2]
                examples.append((reps[0], reps[1]))
    return n_split, examples


def main(n):
    rows = load(n)
    total = len(rows)
    print(f"=== n={n}: {total} digraphs (up to iso) ===")
    for key, label in [
        ("dwl_key", "D-WL"),
        ("mh_key", "MH ranks (k<=3)"),
        ("mh1_key", "MH_(1,l) ranks row"),
        ("dist_key", "distance multiset"),
        ("ph_key", "PH (H_0..H_2)"),
        ("rh_key", "RH"),
        ("all_hom_key", "MH+PH+RH combined"),
    ]:
        print(f"  {label:32s}: {len(classes(rows, key))} classes")

    pairs = [
        ("dwl_key", "mh_key", "A: D-WL-equal, MH-different"),
        ("dwl_key", "ph_key", "B: D-WL-equal, PH-different"),
        ("mh_key", "dwl_key", "C: MH-equal, D-WL-different"),
        ("mh_key", "ph_key", "D: MH-equal, PH-different (E1 !=> E2)"),
        ("all_hom_key", "dwl_key", "E: MH+PH+RH-equal, D-WL-different"),
        ("mh1_key", "mh_key", "MH1-equal, MH{2,3}-different"),
        ("ph_key", "mh_key", "PH-equal, MH-different"),
        ("rh_key", "ph_key", "RH-equal, PH-different"),
        ("dist_key", "dwl_key", "dist-equal, D-WL-different"),
        ("dwl_key", "dist_key", "D-WL-equal, dist-different"),
    ]
    for ki, kj, label in pairs:
        cnt, ex = split_witnesses(rows, ki, kj)
        print(f"  {label:45s}: {cnt} classes split")
        for a, b in ex[:2]:
            print(f"      pair: {a['g6']}  vs  {b['g6']}")

    # F: D-WL-equal AND distance-multiset-equal, but MH different
    # (magnitude homology beyond both refinement and distance counting)
    d = classes(rows, "dwl_key")
    f_examples = []
    f_count = 0
    for group in d.values():
        if len(group) < 2:
            continue
        sub = defaultdict(list)
        for r in group:
            sub[r["dist_key"]].append(r)
        for g2 in sub.values():
            if len(g2) < 2:
                continue
            sub2 = defaultdict(list)
            for r in g2:
                sub2[r["mh_key"]].append(r)
            if len(sub2) > 1:
                f_count += 1
                if len(f_examples) < 3:
                    reps = [v[0] for v in sub2.values()][:2]
                    f_examples.append((reps[0], reps[1]))
    print(f"  {'F: D-WL & dist-equal, MH-different':45s}: {f_count} classes split")
    for a, b in f_examples:
        print(f"      pair: {a['g6']}  vs  {b['g6']}")
        print(f"        {arcs_str(a['g6'])}")
        print(f"        {arcs_str(b['g6'])}")

    # G: D-WL-equal AND dist-equal, PH-different
    g_examples = []
    g_count = 0
    for group in d.values():
        if len(group) < 2:
            continue
        sub = defaultdict(list)
        for r in group:
            sub[r["dist_key"]].append(r)
        for g2 in sub.values():
            if len(g2) < 2:
                continue
            sub2 = defaultdict(list)
            for r in g2:
                sub2[r["ph_key"]].append(r)
            if len(sub2) > 1:
                g_count += 1
                if len(g_examples) < 3:
                    reps = [v[0] for v in sub2.values()][:2]
                    g_examples.append((reps[0], reps[1]))
    print(f"  {'G: D-WL & dist-equal, PH-different':45s}: {g_count} classes split")
    for a, b in g_examples:
        print(f"      pair: {a['g6']}  vs  {b['g6']}")
        print(f"        {arcs_str(a['g6'])}")
        print(f"        {arcs_str(b['g6'])}")


if __name__ == "__main__":
    for n in [int(x) for x in sys.argv[1:]] or [3, 4]:
        main(n)
        print()
