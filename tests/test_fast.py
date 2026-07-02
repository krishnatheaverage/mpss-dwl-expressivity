"""Cross-validate the fast sparse pipeline against the Fraction-based one.

The two implementations share no linear algebra code, so exact agreement
on thousands of digraphs is a strong correctness check for both.
"""

import random
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.digraph import parse_digraph6
from src.fastinv import fast_magnitude_ranks, fast_path_homology
from src.maghom import magnitude_homology_ranks
from src.pathhom import path_homology_ranks


def digraphs(n):
    out = subprocess.run(
        f"geng -q {n} | directg -q", shell=True, capture_output=True, text=True
    )
    return out.stdout.strip().split("\n")


def main():
    bad = 0
    for line in digraphs(4):
        g = parse_digraph6(line)
        if fast_magnitude_ranks(*g, kmax=3) != magnitude_homology_ranks(*g, kmax=3):
            print("MH mismatch:", line)
            bad += 1
        if fast_path_homology(*g, pmax=2) != path_homology_ranks(*g, pmax=2):
            print("PH mismatch:", line)
            bad += 1
    print(f"n=4: all 218 digraphs cross-validated, {bad} mismatches")

    lines = digraphs(5)
    random.seed(1)
    sample = random.sample(lines, 400)
    for line in sample:
        g = parse_digraph6(line)
        if fast_magnitude_ranks(*g, kmax=3) != magnitude_homology_ranks(*g, kmax=3):
            print("MH mismatch:", line)
            bad += 1
        if fast_path_homology(*g, pmax=2) != path_homology_ranks(*g, pmax=2):
            print("PH mismatch:", line)
            bad += 1
    print(f"n=5: 400-digraph sample cross-validated, {bad} total mismatches")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
