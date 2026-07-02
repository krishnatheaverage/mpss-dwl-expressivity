"""Directed 1-WL (color refinement for digraphs), and the
reciprocity-aware refinement DWL+.

D-WL refines vertex colors by the multisets of in-neighbor colors and
out-neighbor colors separately (the standard directed color refinement,
the expressivity bound for directed message-passing GNNs). DWL+ refines
by the triple (in-only, out-only, reciprocal) neighbor-color multisets.

Colors are canonical ACROSS processes and runs: a color is the blake2b
digest of the structural signature, so color histograms from different
worker processes are directly comparable. (A digest collision could only
spuriously MERGE two classes, never separate equal ones; every witness
pair reported in the paper is re-verified with exact structural
comparison via `dwl_equivalent`, so collisions cannot create false
claims, and any claimed equality that matters is re-checked pairwise.)
"""

import hashlib

from .digraph import in_neighbors, out_neighbors

_INIT = "i"


def _h(sig):
    return hashlib.blake2b(sig.encode(), digest_size=12).hexdigest()


def dwl_colors(n, arcs, rounds=None):
    """Run D-WL refinement; return the per-round canonical histograms.

    The graph invariant is the tuple of sorted color histograms, one per
    round, run for `rounds` iterations (default n, which suffices for
    color-refinement stabilization since a partition of an n-set can
    strictly refine at most n-1 times).
    """
    if rounds is None:
        rounds = max(n, 1)
    inn = in_neighbors(n, arcs)
    out = out_neighbors(n, arcs)
    colors = [_INIT] * n
    history = []
    for _ in range(rounds):
        new = []
        for v in range(n):
            sig = (
                colors[v]
                + "|"
                + ",".join(sorted(colors[u] for u in inn[v]))
                + "|"
                + ",".join(sorted(colors[u] for u in out[v]))
            )
            new.append(_h(sig))
        colors = new
        hist = {}
        for c in colors:
            hist[c] = hist.get(c, 0) + 1
        history.append(tuple(sorted(hist.items())))
    return tuple(history)


def dwl_plus_colors(n, arcs, rounds=None):
    """Reciprocity-aware refinement DWL+: (in-only, out-only, reciprocal)."""
    if rounds is None:
        rounds = max(n, 1)
    inn = [set() for _ in range(n)]
    out = [set() for _ in range(n)]
    for u, v in arcs:
        out[u].add(v)
        inn[v].add(u)
    colors = [_INIT] * n
    history = []
    for _ in range(rounds):
        new = []
        for v in range(n):
            rec = inn[v] & out[v]
            sig = (
                colors[v]
                + "|"
                + ",".join(sorted(colors[u] for u in inn[v] - rec))
                + "|"
                + ",".join(sorted(colors[u] for u in out[v] - rec))
                + "|"
                + ",".join(sorted(colors[u] for u in rec))
            )
            new.append(_h(sig))
        colors = new
        hist = {}
        for c in colors:
            hist[c] = hist.get(c, 0) + 1
        history.append(tuple(sorted(hist.items())))
    return tuple(history)


def dwl_invariant(n, arcs):
    """Hashable D-WL graph invariant: per-round canonical histograms."""
    return (n, dwl_colors(n, arcs))


def dwl_equivalent(g1, g2):
    """Exact D-WL equivalence (same n; common round count)."""
    n1, a1 = g1
    n2, a2 = g2
    if n1 != n2:
        return False
    r = max(n1, n2, 1)
    return dwl_colors(n1, a1, rounds=r) == dwl_colors(n2, a2, rounds=r)


def dwl_plus_equivalent(g1, g2):
    n1, a1 = g1
    n2, a2 = g2
    if n1 != n2:
        return False
    r = max(n1, n2, 1)
    return dwl_plus_colors(n1, a1, rounds=r) == dwl_plus_colors(n2, a2, rounds=r)
