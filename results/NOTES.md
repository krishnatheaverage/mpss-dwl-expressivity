# Verified findings: MPSS invariants vs directed WL

All facts below are outputs of exact-rational computations in this repo,
cross-validated by tests/test_known.py (69/69: Leinster-Shulman Euler
characteristic identity vs the magnitude power series on 7 graphs, d^2=0
exhaustively on all 218 n=4 digraphs, GLMY examples, Hepworth-Roff
contractibility, Hepworth-Willerton K_n diagonality). Star witnesses
were additionally re-verified standalone with deeper degree ranges and
per-witness Euler-identity checks.

Conventions: all digraphs up to isomorphism from nauty (geng|directg);
MH = magnitude homology ranks MH_{k,l} over Q (sweep window k <= 3, all
l; witnesses re-checked to k <= 5, l <= 10); PH = GLMY path homology
dims H_0..H_2 over Q; RH = reachability homology dims RH_0..RH_3 over Q;
D-WL = directed color refinement (separate in/out multisets), canonical
per-round histograms; dist = multiset of off-diagonal directed distances
(inf allowed).

## Exhaustive equivalence-class counts (up to iso)

| n | digraphs | D-WL | MH(k<=3) | dist | PH | RH | MH+PH+RH |
|---|----------|------|----------|------|----|----|----------|
| 3 | 16       | 16   | 11       | 7    | 4  | 3  | 11       |
| 4 | 218      | 215  | 74       | 48   | 6  | 5  | 74       |
| 5 | 9608     | 9567 | 1415     | 237  | 11 | 8  | 1417     |

## T1. Incomparability of D-WL and the MPSS invariants

(a) Homology exceeds D-WL: minimal witnesses at n=4 (n=3 exhaustively
    has none; 2 classes split at n=4, 25 at n=5).
    Star witness: directed 4-cycle (&CGS_) vs two digons (&CGWO).
    Both 1-in-1-out regular => D-WL-equivalent. Separated by MH
    (bidegrees (2,4),(3,5) vs (2,2),(3,3)), by PH ([1,1,0] vs [2,0,0]),
    and by RH ([1,0,0,0] vs [2,0,0,0]).
(b) D-WL exceeds MH+PH+RH jointly: minimal witnesses at n=3
    (3 classes; 38 at n=4; 1197 at n=5).
    Star witness: in-star a->c<-b vs directed 2-path a->c->b.
    All three homology theories agree; D-WL separates in round 1.

## T2. E1 ranks do not determine E2 (MH ranks !=> PH)

- Exhaustive: at n <= 4, ZERO pairs with equal MH(k<=3) ranks have
  different PH. (Since full MH refines MH(k<=3)... more precisely the
  k<=3 statement implies the full-MH statement a fortiori on n<=4.)
- Minimal counterexample at n=5, and it is FULL-MH (all bidegrees):
    X = in-star K_{4->1} (&DACGO?, arcs (i,4), i=0..3)
    Y = directed K_{2,2} + isolated vertex (&DEK???)
  Neither has a composable pair of arcs => all magnitude chains vanish
  in degrees k >= 2 for both (one-line hand proof), and MH ranks are
  {(0,0):5, (1,1):4} identically for both.
  But PH(X) = [1,0,0], PH(Y) = [2,1,0] (the unfilled directed square),
  RH(X) = [1,0,0,0], RH(Y) = [2,1,0,0] (order complex of 2x2 bipartite
  poset ~ S^1).
  Machine-checked to k<=5, l<=10 + Euler identity per witness.
- Second n=5 pair (&DKCJN? vs &DKWMM?) agrees at k<=3 but SPLITS at
  bidegree (4,4) (ranks 4 vs 6): a truncation-only witness.

## T3. MH separates beyond D-WL and beyond distances

- At n <= 4: every homological separation of a D-WL class factors
  through the distance multiset (F,G = 0 exhaustively).
- At n=5: 3 classes are D-WL-equal AND dist-equal yet MH-different;
  1 of them PH-different too. Star witness:
    Circ(5;{2,3}) (&DKMME?, verified equal to the circulant)
    vs non-circulant twin (&DKLMI?).
  Both 2-in-2-out regular, same distance multiset (10x1, 10x2), both
  strongly connected (same RH [1,0,0,0]), D-WL-equivalent.
  MH differs at (2,3): 10 vs 6 and (3,4): 30 vs 26.
  PH differs qualitatively: [1,1,0] vs [1,0,1] (an H_1 hole vs an H_2
  hole). Machine-verified + per-witness Euler identity.

## T4. MH ranks do not determine the distance multiset

n=3 star witness (same pair as T1b): in-star ({1,1} plus inf's) vs
2-path ({1,1,2} plus inf's) have identical MH in all bidegrees (both
have no composable arcs). Chain-level dims see distances; homology
forgets them. Related structure lemma (proved by hand, matches all
computed data):
  MH_{1,1} = #arcs and MH_{1,l} = 0 for l >= 2, for EVERY digraph
  (any pair at distance l >= 2 has a geodesic midpoint, which bounds it).
So all fine structure of MH on digraphs lives in degrees k >= 2.

## T5. Within the MPSS family

- PH does not determine RH (1 class at n=4, 3 at n=5), but MH+PH
  jointly determine RH on ALL digraphs with n <= 5 (0 splits).
- PH-equal, MH-different: many (4 classes n=4, 8 n=5): E1 strictly
  refines E2 as expected.
- Combined MH+PH+RH classes = 1417 = 1415 (MH) + 2 (exactly the T2
  witnesses): internal consistency.

All proof obligations listed in earlier versions of these notes are
now discharged in the manuscript (paper/main.tex), which supersedes
this file; the tables there are regenerated from the sweep data by
src/gen_tables.py.
