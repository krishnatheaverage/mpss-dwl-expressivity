"""Reachability homology of a digraph (Hepworth-Roff).

RH_*(G) is the homology of the nerve of the reachability preorder.
The preorder category is equivalent to its condensation poset (strongly
connected components ordered by reachability), and equivalent categories
have homotopy-equivalent nerves, so RH_*(G) is the simplicial homology
of the order complex of the condensation poset. Exact ranks over Q.
"""

from .digraph import INF, distance_matrix
from .linalg import rank


def condensation_poset(n, arcs):
    """Strongly connected components + strict reachability order.

    Returns (num_components, set of strict relations (i, j) meaning
    component i strictly reaches component j).
    """
    D = distance_matrix(n, arcs)
    reach = [[D[x][y] != INF for y in range(n)] for x in range(n)]
    comp = [-1] * n
    nc = 0
    for v in range(n):
        if comp[v] == -1:
            for u in range(n):
                if reach[v][u] and reach[u][v]:
                    comp[u] = nc
            nc += 1
    rel = set()
    for x in range(n):
        for y in range(n):
            if comp[x] != comp[y] and reach[x][y]:
                rel.add((comp[x], comp[y]))
    return nc, rel


def order_complex_homology(nc, rel, dmax=3):
    """Simplicial homology (dims 0..dmax) of the order complex of a poset."""
    # chains: strictly increasing sequences in the partial order
    faces = {0: [(i,) for i in range(nc)]}
    for d in range(1, dmax + 2):
        cur = []
        for c in faces.get(d - 1, []):
            for j in range(nc):
                if (c[-1], j) in rel:
                    cur.append(c + (j,))
        faces[d] = cur
    branks = {}
    for d in range(1, dmax + 2):
        upper, lower = faces.get(d, []), faces.get(d - 1, [])
        if not upper or not lower:
            branks[d] = 0
            continue
        idx = {t: i for i, t in enumerate(lower)}
        M = [[0] * len(upper) for _ in range(len(lower))]
        for j, t in enumerate(upper):
            for i in range(len(t)):
                face = t[:i] + t[i + 1 :]
                M[idx[face]][j] += (-1) ** i
        branks[d] = rank(M)
    H = []
    for d in range(dmax + 1):
        dim = len(faces.get(d, []))
        H.append(dim - branks.get(d, 0) - branks.get(d + 1, 0))
    return H


def reachability_homology_ranks(n, arcs, dmax=3):
    nc, rel = condensation_poset(n, arcs)
    return order_complex_homology(nc, rel, dmax=dmax)
