"""GLMY path homology of digraphs over Q (Grigor'yan-Lin-Muranov-Yau).

Regular p-paths: vertex sequences with consecutive entries distinct.
Allowed p-paths: consecutive pairs are arcs.
Boundary: full alternating face map with irregular faces set to 0
(the standard regularized boundary; d^2 = 0, asserted in tests).
Omega_p = { v in span(A_p) : dv in span(A_{p-1}) }  (d-invariant paths).
H_p = homology of (Omega_*, d), exact ranks over Q.
"""

from fractions import Fraction

from .digraph import out_neighbors
from .linalg import nullspace, rank


def allowed_paths(n, arcs, pmax):
    """A[p] = list of allowed p-paths (as tuples), p = 0..pmax."""
    out = out_neighbors(n, arcs)
    A = [[(v,) for v in range(n)]]
    for p in range(1, pmax + 1):
        cur = []
        for t in A[-1]:
            for v in out[t[-1]]:
                cur.append(t + (v,))
        A.append(cur)
    return A


def _boundary_terms(t):
    """Full regularized boundary of a regular path: [(sign, face), ...]."""
    p = len(t) - 1
    terms = []
    for i in range(p + 1):
        face = t[:i] + t[i + 1 :]
        # regularization: drop faces with equal consecutive entries
        ok = all(face[j] != face[j + 1] for j in range(len(face) - 1))
        if ok:
            terms.append(((-1) ** i, face))
    return terms


def omega_bases(n, arcs, pmax):
    """Bases of Omega_p in A_p coordinates, p = 0..pmax.

    Returns (A, bases) where bases[p] is a list of vectors (Fractions)
    of length |A_p|.
    """
    A = allowed_paths(n, arcs, pmax)
    allowed_sets = [set(ap) for ap in A]
    bases = []
    for p in range(pmax + 1):
        if p == 0:
            bases.append(
                [[Fraction(int(i == j)) for j in range(len(A[0]))] for i in range(len(A[0]))]
            )
            continue
        # rows: non-allowed regular faces that appear; columns: A_p
        row_index = {}
        rows = []
        cols = len(A[p])
        entries = []  # (row, col, val)
        for j, t in enumerate(A[p]):
            for s, face in _boundary_terms(t):
                if face not in allowed_sets[p - 1]:
                    r = row_index.setdefault(face, len(row_index))
                    entries.append((r, j, s))
        if row_index:
            M = [[0] * cols for _ in range(len(row_index))]
            for r, j, s in entries:
                M[r][j] += s
            basis = nullspace(M, ncols=cols)
        else:
            basis = [
                [Fraction(int(i == j)) for j in range(cols)] for i in range(cols)
            ]
        bases.append(basis)
    return A, bases


def _boundary_on_omega(A, bases, p):
    """Matrix of d: Omega_p -> span(A_{p-1}), columns = Omega_p basis."""
    idx = {t: i for i, t in enumerate(A[p - 1])}
    rows = len(A[p - 1])
    cols = len(bases[p])
    M = [[Fraction(0)] * cols for _ in range(rows)]
    for j, vec in enumerate(bases[p]):
        for cidx, coef in enumerate(vec):
            if coef == 0:
                continue
            t = A[p][cidx]
            for s, face in _boundary_terms(t):
                r = idx.get(face)
                if r is not None:
                    M[r][j] += s * coef
                # faces outside A_{p-1} have total coefficient 0 by
                # construction of Omega_p; no need to track them
    return M


def path_homology_ranks(n, arcs, pmax):
    """Exact dims of GLMY path homology H_0..H_pmax over Q.

    Computes Omega up to pmax+1 (for the incoming boundary at pmax).
    """
    A, bases = omega_bases(n, arcs, pmax + 1)
    dims = [len(b) for b in bases]
    branks = [0] * (pmax + 2)  # branks[p] = rank d: Omega_p -> Omega_{p-1}
    for p in range(1, pmax + 2):
        if dims[p] == 0 or dims[p - 1] == 0:
            branks[p] = 0
            continue
        M = _boundary_on_omega(A, bases, p)
        branks[p] = rank(M)
    H = []
    for p in range(pmax + 1):
        out_rank = branks[p] if p >= 1 else 0
        in_rank = branks[p + 1]
        H.append(dims[p] - out_rank - in_rank)
    return H


def check_d_squared_zero(n, arcs, pmax=3):
    """Sanity: d^2 = 0 on span(A_p) with the regularized boundary."""
    A = allowed_paths(n, arcs, pmax)
    for p in range(2, pmax + 1):
        for t in A[p]:
            acc = {}
            for s1, f1 in _boundary_terms(t):
                for s2, f2 in _boundary_terms(f1):
                    acc[f2] = acc.get(f2, 0) + s1 * s2
            if any(v != 0 for v in acc.values()):
                return False, t
    return True, None
