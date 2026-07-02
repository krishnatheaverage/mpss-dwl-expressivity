"""Exact linear algebra over Q for small integer matrices.

All homology ranks in this project are computed with exact rational
arithmetic (Fraction-based Gaussian elimination); no floating point is
used anywhere, so every reported rank is exact.
"""

from fractions import Fraction


def rank(mat):
    """Exact rank over Q of a matrix given as list of rows (ints/Fractions)."""
    if not mat or not mat[0]:
        return 0
    m = [[Fraction(x) for x in row] for row in mat]
    nrows, ncols = len(m), len(m[0])
    r = 0
    for c in range(ncols):
        pivot = None
        for i in range(r, nrows):
            if m[i][c] != 0:
                pivot = i
                break
        if pivot is None:
            continue
        m[r], m[pivot] = m[pivot], m[r]
        pv = m[r][c]
        m[r] = [x / pv for x in m[r]]
        for i in range(nrows):
            if i != r and m[i][c] != 0:
                f = m[i][c]
                m[i] = [a - f * b for a, b in zip(m[i], m[r])]
        r += 1
        if r == nrows:
            break
    return r


def nullspace(mat, ncols=None):
    """Exact basis of the null space {v : mat @ v = 0} over Q.

    `mat` is a list of rows; if empty, `ncols` must give the ambient dim.
    Returns a list of basis vectors (lists of Fractions).
    """
    if not mat:
        assert ncols is not None
        return [[Fraction(int(i == j)) for j in range(ncols)] for i in range(ncols)]
    ncols = len(mat[0])
    m = [[Fraction(x) for x in row] for row in mat]
    nrows = len(m)
    pivots = []  # (row, col)
    r = 0
    for c in range(ncols):
        pivot = None
        for i in range(r, nrows):
            if m[i][c] != 0:
                pivot = i
                break
        if pivot is None:
            continue
        m[r], m[pivot] = m[pivot], m[r]
        pv = m[r][c]
        m[r] = [x / pv for x in m[r]]
        for i in range(nrows):
            if i != r and m[i][c] != 0:
                f = m[i][c]
                m[i] = [a - f * b for a, b in zip(m[i], m[r])]
        pivots.append((r, c))
        r += 1
        if r == nrows:
            break
    pivot_cols = {c for _, c in pivots}
    free_cols = [c for c in range(ncols) if c not in pivot_cols]
    basis = []
    for fc in free_cols:
        v = [Fraction(0)] * ncols
        v[fc] = Fraction(1)
        for pr, pc in pivots:
            v[pc] = -m[pr][fc]
        basis.append(v)
    return basis
