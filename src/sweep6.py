"""Parallel exhaustive sweep for n=6 (1,540,944 digraphs).

Shards the digraph6 stream over worker processes; each worker writes
compact JSON lines. Colors are canonical digests (src/dwl.py), so keys
are comparable across shards. All ranks exact over Q, as everywhere.

Usage:
  python3 src/sweep6.py enumerate      # write results/all6.d6
  python3 src/sweep6.py sample 200     # timing estimate on a random sample
  python3 src/sweep6.py run [workers]  # full sweep -> results/sweep_n6_shard*.jsonl
"""

import json
import multiprocessing as mp
import random
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.digraph import INF, distance_matrix, parse_digraph6
from src.dwl import dwl_colors, dwl_plus_colors
from src.fastinv import fast_magnitude_ranks, fast_path_homology
from src.reach import reachability_homology_ranks

ROOT = Path(__file__).resolve().parents[1]
D6 = ROOT / "results" / "all6.d6"
KMAX = 3
PMAX = 2


def compact_row(line):
    g6 = line.strip()
    n, arcs = parse_digraph6(g6)
    mh = fast_magnitude_ranks(n, arcs, kmax=KMAX)
    ph = fast_path_homology(n, arcs, pmax=PMAX)
    rh = reachability_homology_ranks(n, arcs)
    dwl = dwl_colors(n, arcs)
    dwlp = dwl_plus_colors(n, arcs)
    D = distance_matrix(n, arcs)
    dist = sorted(
        (10**6 if D[x][y] == INF else int(D[x][y]))
        for x in range(n)
        for y in range(n)
        if x != y
    )
    import hashlib

    def hkey(obj):
        return hashlib.blake2b(repr(obj).encode(), digest_size=10).hexdigest()

    return json.dumps(
        {
            "g6": g6,
            "mh": sorted([k, l, r] for (k, l), r in mh.items()),
            "ph": ph,
            "rh": rh,
            "dwl": hkey(dwl),
            "dwlp": hkey(dwlp),
            "dist": hkey(dist),
        },
        separators=(",", ":"),
    )


def worker(args):
    shard_id, start, count = args
    out_path = ROOT / "results" / f"sweep_n6_shard{shard_id:02d}.jsonl"
    done = 0
    t0 = time.time()
    with open(D6) as f, open(out_path, "w") as out:
        for i, line in enumerate(f):
            if i < start:
                continue
            if i >= start + count:
                break
            out.write(compact_row(line) + "\n")
            done += 1
            if done % 20000 == 0:
                el = time.time() - t0
                print(
                    f"[shard {shard_id}] {done}/{count} ({el/done*1000:.0f} ms/graph)",
                    flush=True,
                )
    print(f"[shard {shard_id}] DONE {done} in {time.time()-t0:.0f}s", flush=True)
    return done


def main():
    cmd = sys.argv[1]
    if cmd == "enumerate":
        with open(D6, "w") as f:
            subprocess.run("geng -q 6 | directg -q", shell=True, stdout=f)
        n = sum(1 for _ in open(D6))
        print(f"wrote {D6}: {n} digraphs")
    elif cmd == "sample":
        k = int(sys.argv[2])
        lines = D6.read_text().splitlines()
        random.seed(0)
        sample = random.sample(lines, k)
        t0 = time.time()
        for line in sample:
            compact_row(line)
        per = (time.time() - t0) / k
        total_h = per * len(lines) / 3600
        print(f"{per*1000:.1f} ms/graph -> single-core {total_h:.1f} h "
              f"-> {total_h/9:.1f} h on 9 workers")
    elif cmd == "run":
        workers = int(sys.argv[2]) if len(sys.argv) > 2 else 9
        total = sum(1 for _ in open(D6))
        per = (total + workers - 1) // workers
        jobs = [(i, i * per, min(per, total - i * per)) for i in range(workers)]
        t0 = time.time()
        with mp.Pool(workers) as pool:
            counts = pool.map(worker, jobs)
        print(f"TOTAL {sum(counts)} digraphs in {(time.time()-t0)/3600:.2f} h")
    else:
        raise SystemExit(f"unknown command {cmd}")


if __name__ == "__main__":
    main()
