# Directed Weisfeiler–Leman vs. the magnitude–path spectral sequence

Exact-arithmetic comparison of directed color refinement (D-WL) with the
three homology theories organized by the magnitude–path spectral sequence
(MPSS) of a digraph: magnitude homology (E1), path homology (an axis of
E2), and reachability homology (the target).

Accompanies the manuscript in [`paper/`](paper/). All homology ranks are
computed exactly over Q (fraction-free Gaussian elimination); no floating
point is used anywhere.

## What's here

- `src/maghom.py` — magnitude homology `MH_{k,l}(G)` of a digraph (via the
  directed shortest-path quasimetric) and the magnitude power series.
- `src/pathhom.py` — GLMY path homology `PH_p(G)` over Q.
- `src/reach.py` — reachability homology `RH_*(G)` via the condensation
  poset's order complex.
- `src/dwl.py` — directed 1-WL refinement (D-WL) and the reciprocity-aware
  refinement DWL+.
- `src/digraph.py` — digraph6 I/O (nauty), directed distances, constructions.
- `src/sweep.py` — exhaustive sweep over all digraphs on `n` vertices.
- `src/analyze.py` — equivalence-class counts and witness-pair mining.
- `src/gen_tables.py` — regenerates `paper/tables_generated.tex`.
- `tests/test_known.py` — 69 checks against independently known values.
- `results/` — sweep outputs (`sweep_n{3,4,5}.jsonl`) and `NOTES.md`.

To our knowledge this is the first public implementation computing
magnitude homology of general digraphs, and the first computing all three
MPSS-related theories in one framework.

## Reproduce

Requirements: Python 3.10+ (standard library only) and nauty
(`geng`, `directg`) on PATH.

```
python3 tests/test_known.py          # 69/69 correctness checks
python3 src/sweep.py 3                # -> results/sweep_n3.jsonl
python3 src/sweep.py 4
python3 src/sweep.py 5                # ~4.5 min (9608 digraphs)
python3 src/analyze.py 3 4 5          # class counts + witness pairs
python3 src/gen_tables.py             # -> paper/tables_generated.tex
```

## Correctness

The implementation is validated before use, not asserted:

- **Magnitude homology** against its defining property — the graded Euler
  characteristic equals the magnitude power series (computed independently
  from the similarity matrix) — for all gradings `l <= 6` on seven
  digraphs/graphs; plus Hepworth–Willerton diagonality and ranks
  `n(n-1)^k` for complete graphs.
- **Path homology** against `d^2 = 0` exhaustively on all 218 digraphs
  with `n = 4`, and against the standard GLMY examples (directed cycles,
  filled/unfilled triangles and squares).
- **Reachability homology** against Hepworth–Roff contractibility of
  strongly connected digraphs and known small posets.

## Headline findings (all exhaustive on `n <= 5`, exact)

1. D-WL and `{MH, PH, RH}` are **incomparable**; minimal witnesses at
   `n = 4` (homology beats refinement: directed C4 vs. two digons) and
   `n = 3` (refinement beats all three: in-star vs. directed 2-path).
   Incomparability persists for DWL+, minimal at `n = 6`.
2. The **E1 rank function does not determine E2**: layered digraphs
   (every vertex a source or sink) with equal vertex/arc counts have
   identical magnitude homology in every bidegree but different path and
   reachability homology. Minimal at `n = 5`; none at `n <= 4`.
3. Magnitude homology separates D-WL-equivalent digraphs with **equal
   directed-distance multisets** (minimal at `n = 5`), e.g. the symmetric
   digraph of C5 vs. a non-circulant twin.

## License

MIT. Digraph enumeration uses nauty (McKay–Piperno).
