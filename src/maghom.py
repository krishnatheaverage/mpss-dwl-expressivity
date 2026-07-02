"""Magnitude homology of digraphs, with exact ranks over Q.

A digraph is a [0,inf]-enriched category via the directed shortest-path
quasimetric d (d(x,y) = inf over directed paths, inf if unreachable).
Following Hepworth-Willerton (graphs) / Leinster-Shulman (enriched):

  MC_{k,l} = free module on tuples (x_0,...,x_k), x_i != x_{i+1},
             all d(x_i, x_{i+1}) finite, and sum_i d(x_i,x_{i+1}) = l.
  d_i(x) drops x_i (0 < i < k) when d(x_{i-1},x_i) + d(x_i,x_{i+1})
             = d(x_{i-1},x_{i+1}); the differential is sum (-1)^i d_i.
  MH_{k,l} = H_k of the l-graded piece.

Also computes the magnitude power series of the digraph, so the
Leinster-Shulman Euler characteristic identity
  coeff of q^l in Mag(q)  ==  sum_k (-1)^k rank MH_{k,l}
can be used as an end-to-end correctness test.
"""

from fractions import Fraction
from itertools import product

from .digraph import INF, distance_matrix
from .linalg import rank


def _tuples_by_length(n, D, k):
    """All valid magnitude k-chains, grouped by total length l."""
    by_l = {}
    if k == 0:
        by_l[0] = [(v,) for v in range(n)]
        return by_l
    # extend (k-1)-tuples
    prev = _tuples_by_length(n, D, k - 1)
    for l, tups in prev.items():
        for t in tups:
            last = t[-1]
            for v in range(n):
                if v != last and D[last][v] is not INF and D[last][v] != INF:
                    d = D[last][v]
                    if d != INF:
                        by_l.setdefault(l + d, []).append(t + (v,))
    return by_l


def magnitude_chain_groups(n, arcs, kmax, lmax=None):
    """Return dict (k,l) -> list of tuples, for 0 <= k <= kmax.

    If lmax is given, tuples with total length > lmax are pruned during
    generation (each further step only increases the length, so this is
    exact for the graded pieces with l <= lmax).
    """
    D = distance_matrix(n, arcs)
    groups = {}
    level = {0: [(v,) for v in range(n)]}
    for l, tups in level.items():
        groups[(0, l)] = tups
    for k in range(1, kmax + 1):
        new_level = {}
        for l, tups in level.items():
            for t in tups:
                last = t[-1]
                row = D[last]
                for v in range(n):
                    if v != last and row[v] != INF:
                        nl = l + row[v]
                        if lmax is None or nl <= lmax:
                            new_level.setdefault(nl, []).append(t + (v,))
        for l, tups in new_level.items():
            groups[(k, l)] = tups
        level = new_level
    return groups, D


def boundary_matrix(tuples_k, tuples_km1, D):
    """Matrix of the magnitude differential MC_k -> MC_{k-1} (same l)."""
    index = {t: i for i, t in enumerate(tuples_km1)}
    rows = len(tuples_km1)
    cols = len(tuples_k)
    M = [[0] * cols for _ in range(rows)]
    for j, t in enumerate(tuples_k):
        k = len(t) - 1
        for i in range(1, k):
            a, b, c = t[i - 1], t[i], t[i + 1]
            if D[a][b] + D[b][c] == D[a][c]:
                s = t[:i] + t[i + 1 :]
                r = index.get(s)
                if r is not None:
                    M[r][j] += (-1) ** i
    return M


def magnitude_homology_ranks(n, arcs, kmax, lmax=None):
    """Exact ranks of MH_{k,l} for 0 <= k <= kmax (and l <= lmax if given).

    Needs chains up to kmax+1 (for the incoming boundary at k = kmax).
    Returns dict {(k,l): rank} with zero entries omitted.
    """
    groups, D = magnitude_chain_groups(n, arcs, kmax + 1, lmax=lmax)
    ls = sorted({l for (_, l) in groups})
    ranks = {}
    # precompute boundary ranks r(k,l) = rank of MC_{k,l} -> MC_{k-1,l}
    brank = {}
    for (k, l), tups in groups.items():
        if k == 0:
            continue
        below = groups.get((k - 1, l), [])
        if not below or not tups:
            brank[(k, l)] = 0
            continue
        M = boundary_matrix(tups, below, D)
        brank[(k, l)] = rank(M)
    for (k, l), tups in groups.items():
        if k > kmax:
            continue
        dim = len(tups)
        out_rank = brank.get((k, l), 0)  # MC_k -> MC_{k-1}
        in_rank = brank.get((k + 1, l), 0)  # MC_{k+1} -> MC_k
        h = dim - out_rank - in_rank
        if h:
            ranks[(k, l)] = h
    return ranks


def magnitude_series(n, arcs, lmax):
    """Power series of the magnitude of the digraph up to degree lmax.

    Z(q)_{xy} = q^{d(x,y)} (0 if unreachable); Mag(q) = sum entries Z^{-1}.
    Z = I + N with N strictly positive degrees, so Z^{-1} = sum (-N)^j.
    Coefficients returned as a list [c_0, ..., c_lmax] of Fractions.
    """
    D = distance_matrix(n, arcs)

    def poly_mul(p, q):
        r = [Fraction(0)] * (lmax + 1)
        for i, a in enumerate(p):
            if a:
                for j, b in enumerate(q):
                    if b and i + j <= lmax:
                        r[i + j] += a * b
        return r

    # N as matrix of polynomials (lists of coefficients)
    N = [[[Fraction(0)] * (lmax + 1) for _ in range(n)] for _ in range(n)]
    for x in range(n):
        for y in range(n):
            if x != y and D[x][y] != INF and D[x][y] <= lmax:
                N[x][y][int(D[x][y])] = Fraction(1)

    # accumulate sum over j of (-1)^j N^j, entrywise, truncated
    total = [[[Fraction(0)] * (lmax + 1) for _ in range(n)] for _ in range(n)]
    for x in range(n):
        total[x][x][0] = Fraction(1)
    power = None
    cur = [[list(N[x][y]) for y in range(n)] for x in range(n)]  # N^1
    sign = -1
    for _ in range(1, lmax + 1):
        for x in range(n):
            for y in range(n):
                for dgr in range(lmax + 1):
                    total[x][y][dgr] += sign * cur[x][y][dgr]
        # next power
        nxt = [[[Fraction(0)] * (lmax + 1) for _ in range(n)] for _ in range(n)]
        for x in range(n):
            for y in range(n):
                for z in range(n):
                    if any(cur[x][z]) and any(N[z][y]):
                        pm = poly_mul(cur[x][z], N[z][y])
                        row = nxt[x][y]
                        for dgr in range(lmax + 1):
                            row[dgr] += pm[dgr]
        cur = nxt
        sign = -sign
    coeffs = [Fraction(0)] * (lmax + 1)
    for x in range(n):
        for y in range(n):
            for dgr in range(lmax + 1):
                coeffs[dgr] += total[x][y][dgr]
    return coeffs
