"""Digraph representation, digraph6 parsing, directed distances.

A digraph on n vertices is stored as (n, arcs) with arcs a frozenset of
ordered pairs (u, v), u != v (no loops: directg output has none, and the
shortest-path quasimetric ignores loops anyway).
"""

from collections import deque

INF = float("inf")


def parse_digraph6(line):
    """Parse one line of nauty digraph6 format ('&' header) -> (n, arcs)."""
    s = line.strip()
    if not s:
        raise ValueError("empty digraph6 line")
    if s[0] != "&":
        raise ValueError(f"not digraph6: {s[:10]!r}")
    body = [ord(c) - 63 for c in s[1:]]
    if body[0] < 63:
        n = body[0]
        bits_start = 1
    else:
        # n >= 63: 126 then 3 bytes (or larger; not needed for n <= 62)
        n = (body[1] << 12) | (body[2] << 6) | body[3]
        bits_start = 4
    bits = []
    for byte in body[bits_start:]:
        for k in range(5, -1, -1):
            bits.append((byte >> k) & 1)
    arcs = set()
    idx = 0
    for u in range(n):
        for v in range(n):
            if idx < len(bits) and bits[idx]:
                if u != v:
                    arcs.add((u, v))
            idx += 1
    return n, frozenset(arcs)


def to_digraph6(n, arcs):
    """Encode (n, arcs) as a digraph6 line (n <= 62)."""
    assert n <= 62
    bits = []
    for u in range(n):
        for v in range(n):
            bits.append(1 if (u, v) in arcs else 0)
    while len(bits) % 6 != 0:
        bits.append(0)
    out = ["&", chr(n + 63)]
    for i in range(0, len(bits), 6):
        byte = 0
        for b in bits[i : i + 6]:
            byte = (byte << 1) | b
        out.append(chr(byte + 63))
    return "".join(out)


def out_neighbors(n, arcs):
    out = [[] for _ in range(n)]
    for u, v in arcs:
        out[u].append(v)
    return out


def in_neighbors(n, arcs):
    inn = [[] for _ in range(n)]
    for u, v in arcs:
        inn[v].append(u)
    return inn


def distance_matrix(n, arcs):
    """Directed shortest-path distances; INF where unreachable; d(x,x)=0."""
    out = out_neighbors(n, arcs)
    D = [[INF] * n for _ in range(n)]
    for s in range(n):
        D[s][s] = 0
        q = deque([s])
        while q:
            u = q.popleft()
            for v in out[u]:
                if D[s][v] is INF or D[s][v] > D[s][u] + 1:
                    D[s][v] = D[s][u] + 1
                    q.append(v)
    return D


def reverse(n, arcs):
    return n, frozenset((v, u) for (u, v) in arcs)


def disjoint_union(g1, g2):
    n1, a1 = g1
    n2, a2 = g2
    arcs = set(a1)
    for u, v in a2:
        arcs.add((u + n1, v + n1))
    return n1 + n2, frozenset(arcs)


def directed_cycle(n):
    return n, frozenset((i, (i + 1) % n) for i in range(n))


def bidirected(n, edges):
    """Undirected graph viewed as a digraph (each edge in both directions)."""
    arcs = set()
    for u, v in edges:
        arcs.add((u, v))
        arcs.add((v, u))
    return n, frozenset(arcs)


def complete_graph(n):
    return bidirected(n, [(i, j) for i in range(n) for j in range(i + 1, n)])


def undirected_cycle(n):
    return bidirected(n, [(i, (i + 1) % n) for i in range(n)])
