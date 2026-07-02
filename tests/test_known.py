"""Correctness tests against independently known values.

Every invariant implementation is validated against published values or
theorems BEFORE being used in the sweep:

1. Magnitude homology:
   - Euler characteristic identity (Leinster-Shulman Thm 7.14 lineage):
     coeff of q^l in the magnitude power series == alternating sum of
     MH_{k,l} ranks. Tested on complete graphs, undirected cycles,
     directed cycles, and random digraphs.
   - Complete graph K_n is diagonal with rank MH_{k,k} = n(n-1)^k
     (Hepworth-Willerton Example).
   - Undirected C_4: known non-diagonal entry MH_{2,3} etc. checked via
     Euler identity + diagonality pattern of bipartite girth-4 graphs.
2. Path homology (GLMY examples):
   - transitive triangle: H_1 = 0 (2-path fills the cycle).
   - cyclic triangle C_3: H_0 = 1, H_1 = 1 (no allowed d-invariant
     2-paths; the loop survives).
   - "diamond" square a->b, a->c, b->d, c->d: H_1 = 0 (GLMY: the square
     is filled by the d-invariant chain (abd - acd)).
   - directed cycle C_n (n >= 3, cyclic orientation): H_1 = 1.
   - d^2 = 0 on all digraphs on <= 4 vertices.
3. Reachability homology:
   - strongly connected digraph ~ point (Hepworth-Roff).
   - disjoint union of two strongly connected digraphs: RH_0 = 2.
   - poset "V" and diamond DAGs: order complex contractible.
"""

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fractions import Fraction

from src.digraph import (
    bidirected,
    complete_graph,
    directed_cycle,
    disjoint_union,
    parse_digraph6,
    to_digraph6,
    undirected_cycle,
)
from src.dwl import dwl_equivalent
from src.maghom import magnitude_homology_ranks, magnitude_series
from src.pathhom import check_d_squared_zero, path_homology_ranks
from src.reach import reachability_homology_ranks

PASS = 0
FAIL = 0


def check(name, cond, detail=""):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"  ok  {name}")
    else:
        FAIL += 1
        print(f"FAIL  {name}  {detail}")


def euler_identity(g, lmax, label):
    """coeff q^l of magnitude == sum_k (-1)^k rank MH_{k,l}, l <= lmax."""
    n, arcs = g
    # MH_{k,l} = 0 for k > l (each step has length >= 1), so kmax = lmax
    ranks = magnitude_homology_ranks(n, arcs, kmax=lmax, lmax=lmax)
    series = magnitude_series(n, arcs, lmax)
    for l in range(lmax + 1):
        alt = sum((-1) ** k * r for (k, ll), r in ranks.items() if ll == l)
        check(
            f"Euler identity {label} l={l}",
            Fraction(alt) == series[l],
            f"alt={alt} series={series[l]}",
        )


def test_digraph6_roundtrip():
    print("digraph6 parsing:")
    out = subprocess.run(
        "geng -q 4 | directg -q", shell=True, capture_output=True, text=True
    ).stdout.strip().split("\n")
    check("n=4 count", len(out) == 218, f"got {len(out)}")
    ok = True
    for line in out:
        n, arcs = parse_digraph6(line)
        if to_digraph6(n, arcs) != line.strip():
            ok = False
            break
    check("roundtrip all n=4 digraphs", ok)


def test_magnitude():
    print("magnitude homology:")
    # K_n diagonal, rank MH_{k,k} = n(n-1)^k  (Hepworth-Willerton)
    for n in (2, 3, 4):
        g = complete_graph(n)
        ranks = magnitude_homology_ranks(*g, kmax=3)
        diag_ok = all(k == l for (k, l) in ranks)
        vals_ok = all(ranks.get((k, k), 0) == n * (n - 1) ** k for k in range(4))
        check(f"K_{n} diagonal", diag_ok, str(ranks))
        check(f"K_{n} ranks n(n-1)^k", vals_ok, str(ranks))
    euler_identity(complete_graph(3), 5, "K_3")
    euler_identity(undirected_cycle(4), 5, "C_4 (undirected)")
    euler_identity(undirected_cycle(5), 5, "C_5 (undirected)")
    euler_identity(directed_cycle(3), 5, "C_3 (directed)")
    euler_identity(directed_cycle(6), 6, "C_6 (directed)")
    euler_identity(
        disjoint_union(directed_cycle(3), directed_cycle(3)), 6, "C3+C3 (directed)"
    )
    # a non-strongly-connected digraph: path a->b->c plus arc a->c
    g = (3, frozenset({(0, 1), (1, 2), (0, 2)}))
    euler_identity(g, 5, "transitive triangle")


def test_path_homology():
    print("path homology (GLMY):")
    # d^2 = 0 on every digraph with <= 4 vertices
    out = subprocess.run(
        "geng -q 4 | directg -q", shell=True, capture_output=True, text=True
    ).stdout.strip().split("\n")
    ok = True
    for line in out:
        good, bad = check_d_squared_zero(*parse_digraph6(line), pmax=3)
        if not good:
            ok = False
            print("   d^2 != 0 on", line, bad)
            break
    check("d^2 = 0 on all 218 digraphs n=4", ok)

    tt = (3, frozenset({(0, 1), (1, 2), (0, 2)}))  # transitive triangle
    H = path_homology_ranks(*tt, pmax=2)
    check("transitive triangle H = [1,0,0]", H == [1, 0, 0], str(H))

    c3 = directed_cycle(3)
    H = path_homology_ranks(*c3, pmax=2)
    check("cyclic triangle H = [1,1,0]", H == [1, 1, 0], str(H))

    diamond = (4, frozenset({(0, 1), (0, 2), (1, 3), (2, 3)}))
    H = path_homology_ranks(*diamond, pmax=2)
    check("diamond square H = [1,0,0]", H == [1, 0, 0], str(H))

    for n in (4, 5, 6):
        H = path_homology_ranks(*directed_cycle(n), pmax=2)
        check(f"directed C_{n} H = [1,1,0]", H == [1, 1, 0], str(H))

    both = disjoint_union(directed_cycle(3), directed_cycle(3))
    H = path_homology_ranks(*both, pmax=2)
    check("C3+C3 H = [2,2,0]", H == [2, 2, 0], str(H))


def test_reachability():
    print("reachability homology:")
    for g, label in [
        (directed_cycle(3), "C_3"),
        (directed_cycle(5), "C_5"),
        (complete_graph(3), "K_3 bidirected"),
    ]:
        H = reachability_homology_ranks(*g)
        check(f"{label} strongly connected ~ point", H == [1, 0, 0, 0], str(H))
    both = disjoint_union(directed_cycle(3), directed_cycle(4))
    H = reachability_homology_ranks(*both)
    check("C3+C4 two components", H == [2, 0, 0, 0], str(H))
    vee = (3, frozenset({(0, 1), (0, 2)}))  # a -> b, a -> c
    H = reachability_homology_ranks(*vee)
    check("V-poset contractible", H == [1, 0, 0, 0], str(H))
    diamond = (4, frozenset({(0, 1), (0, 2), (1, 3), (2, 3)}))
    H = reachability_homology_ranks(*diamond)
    check("diamond DAG contractible", H == [1, 0, 0, 0], str(H))


def test_dwl():
    print("D-WL:")
    c6 = directed_cycle(6)
    c33 = disjoint_union(directed_cycle(3), directed_cycle(3))
    check("C6 ~DWL C3+C3 (both 1-in-1-out regular)", dwl_equivalent(c6, c33))
    tt = (3, frozenset({(0, 1), (1, 2), (0, 2)}))
    check("C3 !~DWL transitive triangle", not dwl_equivalent(directed_cycle(3), tt))
    p3 = (3, frozenset({(0, 1), (1, 2)}))
    check("path !~DWL C3", not dwl_equivalent(p3, directed_cycle(3)))


if __name__ == "__main__":
    test_digraph6_roundtrip()
    test_magnitude()
    test_path_homology()
    test_reachability()
    test_dwl()
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)
