"""Fast exact invariants for the large sweeps.

Same mathematical definitions as maghom.py / pathhom.py, computed via
sparse integer elimination (fastrank.sparse_rank), exact over Q.

Path homology avoids null-space bases entirely via the rank identity:
for the full face map T on span(A_p), split rows into non-allowed (N)
and allowed (D) faces; then
    dim Omega_p        = |A_p| - rank(N_p)
    rank(d|Omega_p)    = rank([N_p; D_p]) - rank(N_p)
    dim PH_p           = dim Omega_p - rank(d|Omega_p) - rank(d|Omega_{p+1})
(the middle line is the standard fact that the rank of D restricted to
ker N equals rank of the stacked matrix minus rank of N).

Cross-validated against the Fraction-based implementations on all
digraphs with n <= 4 and samples at n = 5 (tests/test_fast.py).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.digraph import INF, distance_matrix, out_neighbors
from src.fastrank import sparse_rank
from src.pathhom import _boundary_terms, allowed_paths


def fast_magnitude_ranks(n, arcs, kmax):
    """Exact ranks of MH_{k,l}, 0 <= k <= kmax, all l (sparse pipeline)."""
    D = distance_matrix(n, arcs)
    groups = {}
    level = {0: [(v,) for v in range(n)]}
    groups[(0, 0)] = level[0]
    for k in range(1, kmax + 2):
        new_level = {}
        for l, tups in level.items():
            for t in tups:
                row = D[t[-1]]
                last = t[-1]
                for v in range(n):
                    if v != last and row[v] != INF:
                        new_level.setdefault(l + row[v], []).append(t + (v,))
        for l, tups in new_level.items():
            groups[(k, l)] = tups
        level = new_level

    brank = {}
    for (k, l), tups in groups.items():
        if k == 0:
            continue
        below = groups.get((k - 1, l))
        if not below or not tups:
            brank[(k, l)] = 0
            continue
        index = {t: i for i, t in enumerate(below)}
        cols = []
        for t in tups:
            col = {}
            for i in range(1, k):
                a, b, c = t[i - 1], t[i], t[i + 1]
                if D[a][b] + D[b][c] == D[a][c]:
                    r = index.get(t[:i] + t[i + 1 :])
                    if r is not None:
                        col[r] = col.get(r, 0) + (-1) ** i
            cols.append({r: v for r, v in col.items() if v})
        brank[(k, l)] = sparse_rank(cols)

    ranks = {}
    for (k, l), tups in groups.items():
        if k > kmax:
            continue
        h = len(tups) - brank.get((k, l), 0) - brank.get((k + 1, l), 0)
        if h:
            ranks[(k, l)] = h
    return ranks


def _ph_cols(A, allowed_prev, p):
    """Columns of the face map on A_p, rows tagged ('n'|'d', face)."""
    cols_stacked = []
    cols_n = []
    for t in A[p]:
        cn, cs = {}, {}
        for s, face in _boundary_terms(t):
            tag = ("d", face) if face in allowed_prev else ("n", face)
            cs[tag] = cs.get(tag, 0) + s
            if tag[0] == "n":
                cn[tag] = cn.get(tag, 0) + s
        cols_stacked.append({k: v for k, v in cs.items() if v})
        cols_n.append({k: v for k, v in cn.items() if v})
    return cols_n, cols_stacked


def fast_path_homology(n, arcs, pmax):
    """Exact dims of GLMY path homology H_0..H_pmax (rank identities)."""
    A = allowed_paths(n, arcs, pmax + 1)
    allowed_sets = [set(ap) for ap in A]
    dim_omega = [len(A[0])]
    rank_d = [0] * (pmax + 2)  # rank of d: Omega_p -> Omega_{p-1}
    for p in range(1, pmax + 2):
        if not A[p]:
            dim_omega.append(0)
            rank_d[p] = 0
            continue
        cols_n, cols_stacked = _ph_cols(A, allowed_sets[p - 1], p)
        rn = sparse_rank([c for c in cols_n if c])
        rs = sparse_rank([c for c in cols_stacked if c])
        # columns that are entirely zero in N still count for dim Omega
        dim_omega.append(len(A[p]) - rn)
        rank_d[p] = rs - rn
    H = []
    for p in range(pmax + 1):
        H.append(dim_omega[p] - (rank_d[p] if p >= 1 else 0) - rank_d[p + 1])
    return H
