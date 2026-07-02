"""Memory-light analysis of the n=6 sweep shards.

Two passes over results/sweep_n6_shard*.jsonl:
  pass 1: count classes per invariant and count I-classes split by J
          (via dict key -> Counter of J-keys, with compact digest keys)
  pass 2: extract representative witness pairs for each split relation.

Every witness pair printed here is re-verified EXACTLY (full invariant
recomputation, structural D-WL comparison) by verify_witnesses(), so
digest collisions cannot produce false claims.
"""

import glob
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

ROOT = Path(__file__).resolve().parents[1]


def rows():
    for path in sorted(glob.glob(str(ROOT / "results" / "sweep_n6_shard*.jsonl"))):
        with open(path) as f:
            for line in f:
                yield json.loads(line)


def keyed(r):
    return {
        "dwl": r["dwl"],
        "dwlp": r["dwlp"],
        "dist": r["dist"],
        "mh": json.dumps(r["mh"], separators=(",", ":")),
        "ph": json.dumps(r["ph"], separators=(",", ":")),
        "rh": json.dumps(r["rh"], separators=(",", ":")),
    }


RELATIONS = [
    ("dwl", "mh"), ("dwl", "ph"), ("dwl", "dist"),
    ("dist", "dwl"), ("mh", "dwl"), ("mh", "ph"),
    ("ph", "mh"), ("ph", "rh"),
    ("dwlp", "mh"), ("dwlp", "ph"),
    ("mhph", "rh"), ("hom", "dwl"),
    ("dwl_dist", "mh"), ("dwl_dist", "ph"),
    ("dwlp_dist", "mh"),
]


def combo_keys(k):
    k["mhph"] = k["mh"] + "|" + k["ph"]
    k["hom"] = k["mh"] + "|" + k["ph"] + "|" + k["rh"]
    k["dwl_dist"] = k["dwl"] + "|" + k["dist"]
    k["dwlp_dist"] = k["dwlp"] + "|" + k["dist"]
    return k


def main():
    # pass 1: class counts and split counts
    class_counts = defaultdict(set)
    groups = {rel: defaultdict(set) for rel in RELATIONS}
    total = 0
    for r in rows():
        total += 1
        k = combo_keys(keyed(r))
        for name in ("dwl", "dwlp", "dist", "mh", "ph", "rh", "hom"):
            class_counts[name].add(k[name])
        for ki, kj in RELATIONS:
            groups[(ki, kj)][k[ki]].add(k[kj])
    print(f"n=6: {total} digraphs")
    for name in ("dwl", "dwlp", "dist", "mh", "ph", "rh", "hom"):
        print(f"  classes[{name}] = {len(class_counts[name])}")
    split_relations = {}
    for ki, kj in RELATIONS:
        n_split = sum(1 for v in groups[(ki, kj)].values() if len(v) > 1)
        split_relations[(ki, kj)] = n_split
        print(f"  split {ki}-equal,{kj}-different: {n_split}")

    # pass 2: witnesses for the relations that matter
    want = [rel for rel in RELATIONS if split_relations[rel] > 0]
    seen = {rel: {} for rel in want}
    witnesses = {rel: [] for rel in want}
    for r in rows():
        k = combo_keys(keyed(r))
        for rel in want:
            ki, kj = rel
            if len(groups[rel][k[ki]]) > 1 and len(witnesses[rel]) < 3:
                bucket = seen[rel].setdefault(k[ki], {})
                if k[kj] not in bucket:
                    bucket[k[kj]] = r["g6"]
                    if len(bucket) == 2:
                        witnesses[rel].append(tuple(bucket.values()))
    print("\nwitness pairs (to be exactly re-verified):")
    for rel, pairs in witnesses.items():
        for a, b in pairs[:2]:
            print(f"  {rel[0]}-equal,{rel[1]}-diff: {a}  vs  {b}")
    out = ROOT / "results" / "analysis_n6.json"
    with open(out, "w") as f:
        json.dump(
            {
                "total": total,
                "classes": {k: len(v) for k, v in class_counts.items()},
                "splits": {f"{a}|{b}": c for (a, b), c in split_relations.items()},
                "witnesses": {f"{a}|{b}": v for (a, b), v in witnesses.items()},
            },
            f,
            indent=1,
        )
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
