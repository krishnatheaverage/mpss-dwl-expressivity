"""Degree-2 magnitude homology of digraphs: closed formulas + decomposition.

The magnitude chain complex splits over ordered endpoint pairs (a,b)
(the differential never deletes endpoints). For k = 2 and grading l,
the (a,b)-summand is computed by the "midpoint graph" Gamma_l(a,b):

  vertices: V_l(a,b) = { y != a,b : d(a,y) + d(y,b) = l, both finite }
  for every 4-chain (a,x,y,b) of total length l (x != a,y; y != b):
    cond_A: d(a,x) + d(x,y) = d(a,y)   (deleting x is a face)
    cond_B: d(x,y) + d(y,b) = d(x,b)   (deleting y is a face)
    both        -> edge {x, y}
    exactly one -> pin on the surviving midpoint (y for A, x for B)

  rank MH_{2,l}(a,b) =
    c(Gamma) - 1                 if d(a,b) = l   (no pins occur)
    #components without a pin    otherwise.

Cf. Ivanov-Mukoseev (arXiv:2405.04748), Theorem 4.5, for an equivalent
basis description ("long paths modulo short congruence"); the component
form here is what the expressivity corollaries in the paper use.

Special case l = 2 (no 4-chains exist):
  MH_{2,2}(G) = sum_v indeg(v)*outdeg(v) - #{(x,z) : d(x,z) = 2}.

This module verifies both statements against the general implementation.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.digraph import INF, distance_matrix, in_neighbors, out_neighbors


def mh22_formula(n, arcs):
    """MH_{2,2} = sum indeg*outdeg - #{ordered pairs at distance 2}."""
    inn = in_neighbors(n, arcs)
    out = out_neighbors(n, arcs)
    walks = sum(len(inn[v]) * len(out[v]) for v in range(n))
    D = distance_matrix(n, arcs)
    d2 = sum(
        1 for x in range(n) for z in range(n) if x != z and D[x][z] == 2
    )
    return walks - d2


def midpoint_graph(n, D, a, b, l):
    """Vertices, edges, pinned set of Gamma_l(a,b)."""
    V = [
        y
        for y in range(n)
        if y != a and y != b and D[a][y] != INF and D[y][b] != INF
        and D[a][y] + D[y][b] == l
    ]
    Vset = set(V)
    edges = set()
    pinned = set()
    # 4-chains (a,x,y,b): x != a, x != y, y != b, total length l
    for x in range(n):
        if x == a or D[a][x] == INF:
            continue
        for y in range(n):
            if y == b or y == x or D[x][y] == INF or D[y][b] == INF:
                continue
            if D[a][x] + D[x][y] + D[y][b] != l:
                continue
            condA = D[a][x] + D[x][y] == D[a][y]  # face (a,y,b)
            condB = D[x][y] + D[y][b] == D[x][b]  # face (a,x,b)
            if condA and condB:
                # both faces valid: edge between midpoints y and x
                if x in Vset and y in Vset:
                    edges.add((min(x, y), max(x, y)))
            elif condA:
                if y in Vset:
                    pinned.add(y)
            elif condB:
                if x in Vset:
                    pinned.add(x)
    return V, edges, pinned


def mh2_pair_rank(n, D, a, b, l):
    """rank MH_{2,l}(a,b) via the midpoint graph."""
    V, edges, pinned = midpoint_graph(n, D, a, b, l)
    if not V:
        return 0
    # union-find over V
    parent = {v: v for v in V}

    def find(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v

    for u, v in edges:
        ru, rv = find(u), find(v)
        if ru != rv:
            parent[ru] = rv
    comps = {}
    for v in V:
        comps.setdefault(find(v), []).append(v)
    if D[a][b] == l:
        return len(comps) - 1
    return sum(1 for c in comps.values() if not any(v in pinned for v in c))


def mh2_via_decomposition(n, arcs, lmax):
    """dict {l: rank MH_{2,l}} computed pair-by-pair, l = 2..lmax."""
    D = distance_matrix(n, arcs)
    out = {}
    for l in range(2, lmax + 1):
        total = 0
        for a in range(n):
            for b in range(n):
                total += mh2_pair_rank(n, D, a, b, l)
        if total:
            out[l] = total
    return out
