"""Directed 1-WL (color refinement for digraphs).

D-WL refines vertex colors by the multisets of in-neighbor colors and
out-neighbor colors separately (the standard directed color refinement,
the expressivity bound for directed message-passing GNNs).

Colors are globally canonical: a color is an interned string built
structurally, so color histograms are comparable ACROSS different
digraphs without running refinement on a disjoint union.
"""

from .digraph import in_neighbors, out_neighbors

_color_table = {}


def _intern(sig):
    c = _color_table.get(sig)
    if c is None:
        c = len(_color_table)
        _color_table[sig] = c
    return c


def dwl_colors(n, arcs, rounds=None):
    """Run D-WL refinement; return the per-round color histograms.

    The graph invariant is the tuple of sorted color histograms, one per
    round, run for `rounds` iterations (default n, which suffices for
    color-refinement stabilization).
    """
    if rounds is None:
        rounds = max(n, 1)
    inn = in_neighbors(n, arcs)
    out = out_neighbors(n, arcs)
    colors = [_intern(("init",))] * n
    history = []
    for _ in range(rounds):
        new = []
        for v in range(n):
            sig = (
                colors[v],
                tuple(sorted(colors[u] for u in inn[v])),
                tuple(sorted(colors[u] for u in out[v])),
            )
            new.append(_intern(sig))
        colors = new
        hist = {}
        for c in colors:
            hist[c] = hist.get(c, 0) + 1
        history.append(tuple(sorted(hist.items())))
    return tuple(history)


def dwl_invariant(n, arcs):
    """Hashable D-WL graph invariant: per-round canonical color histograms."""
    return (n, dwl_colors(n, arcs))


def dwl_equivalent(g1, g2):
    """Two digraphs are D-WL-equivalent iff their invariants agree.

    (Same n required; rounds = max(n1, n2) so histories are comparable.)
    """
    n1, a1 = g1
    n2, a2 = g2
    if n1 != n2:
        return False
    r = max(n1, n2, 1)
    return dwl_colors(n1, a1, rounds=r) == dwl_colors(n2, a2, rounds=r)
