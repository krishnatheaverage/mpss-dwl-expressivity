"""Fast exact rank of sparse integer matrices over Q.

Columns are dicts {row: int}. Elimination is fraction-free: to eliminate
row r (pivot column p, pivot value a) from a column c with value b at r,
replace c by a*c - b*p, which scales c by a != 0 and so preserves rank.
Pivots prefer entries of absolute value 1, in which case entries stay
small; all arithmetic is exact Python integers. Used by the sweeps; the
Fraction-based dense routine in linalg.py is retained and the two are
cross-validated against each other in the test suite.
"""


def sparse_rank(cols):
    """Exact rank over Q; cols is a list of {row_index: int} dicts."""
    # active columns, with an index from row -> set of column ids
    cols = [dict(c) for c in cols if c]
    row_index = {}
    for j, c in enumerate(cols):
        for r in c:
            row_index.setdefault(r, set()).add(j)
    alive = set(range(len(cols)))
    rank = 0
    while alive:
        # choose pivot: prefer |value| == 1, then sparsest column
        best = None
        for j in alive:
            c = cols[j]
            for r, v in c.items():
                if v == 1 or v == -1:
                    best = (j, r, v)
                    break
            if best:
                break
        if best is None:
            # no unit entries anywhere: take any nonzero entry
            j = next(iter(alive))
            r, v = next(iter(cols[j].items()))
            best = (j, r, v)
        j, r, a = best
        piv = cols[j]
        alive.discard(j)
        rank += 1
        # eliminate row r from all other alive columns containing it
        for k in list(row_index.get(r, ())):
            if k == j or k not in alive:
                continue
            c = cols[k]
            b = c.get(r)
            if b is None:
                continue
            # c := a*c - b*piv
            if a != 1:
                for rr in c:
                    c[rr] *= a
            for rr, pv in piv.items():
                nv = c.get(rr, 0) - b * pv
                if nv:
                    c[rr] = nv
                    row_index.setdefault(rr, set()).add(k)
                else:
                    if rr in c:
                        del c[rr]
                        s = row_index.get(rr)
                        if s:
                            s.discard(k)
            if not c:
                alive.discard(k)
        # clear the pivot column from the index
        for rr in piv:
            s = row_index.get(rr)
            if s:
                s.discard(j)
    return rank


def rank_from_dense(mat):
    """Adapter: dense list-of-rows integer matrix -> sparse_rank."""
    if not mat or not mat[0]:
        return 0
    ncols = len(mat[0])
    cols = [{} for _ in range(ncols)]
    for r, row in enumerate(mat):
        for j, v in enumerate(row):
            if v:
                cols[j][r] = int(v)
    return sparse_rank(cols)
