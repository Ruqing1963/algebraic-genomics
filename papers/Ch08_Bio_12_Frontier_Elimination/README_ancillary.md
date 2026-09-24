# Ancillary files for *The double-stranded lattice sum by frontier elimination*

Ruqing Chen, 22 September 2026. Bio 12 of the applied strand; numbered XXV in the series.

## Contents

| file | what it is |
|---|---|
| `anccft25.tex` / `.pdf` | the note (8 pp.) |
| `ds_transfer_matrix.py` | the frontier elimination, the connectivity-class elimination, the matrix-tree connectivity lemma, and the two orderings; 7 checks |
| `ds_transfer_matrix_output.txt` | the transcript of one full run |
| `phix174_NC_001422.1.txt` | bacteriophage phiX174, RefSeq NC_001422.1, 5386 bp, plain sequence |

## Running it

```
python -u ds_transfer_matrix.py
```

**Standalone**: Python 3.8+ and the standard library. No `numpy`, no `sympy`, and no other
file of the series — `D_k` is rebuilt here from the definition rather than imported from
Bio 8, precisely so that check (D) is a real cross-check of that note's Table 1 and not a
tautology. About eighty-five minutes and 10 GB, almost all of it check (M), the connected
count at k = 10. `python -u ds_transfer_matrix.py --quick` skips (M) and runs in six minutes.

## The eight checks

| key | claim | mode |
|---|---|---|
| D | the double cover rebuilt from the definition reproduces `\|V\|`, `#fix_V`, `\|A\|`, `#fix_A`, `sum m` and `rank Lambda^-` of Bio 8's Table 1 at k = 10, 11, 12 | exact |
| T | the frontier elimination and the backtracking enumerator return the same box count, 10 and 528 | exact |
| C | on all 538 box points at k = 12, 11: `supp c` is connected **iff** the reduced Laplacian of the support is non-singular — Lemma 3 | exact |
| W | summing `N_seq` over the box returns 8 and 218 connected points and `#DS` = 10 and 864, reproducing Bio 8 | exact |
| O | coordinate width 6, 18, 62 against greedy width 6, 10, 26; at k = 10 the coordinate order does not finish | exact |
| K | **phiX174 at k = 10: 2,481,687,900 box points**, and the same total after pinning each of two orbits and summing | exact |
| N | the connectivity-class elimination returns Bio 8's 8 and 218 | exact |
| M | **phiX174 at k = 10: 909,103,210 box points with connected support** | exact |

"Exact" means exact integer arithmetic throughout: the elimination adds integers and never
divides, and the determinants of Lemma 3 are computed over **Q** with `fractions`.

## The result

Bio 8 left the k = 10 row of its Table 2 blank — `rank Lambda^- = 48`, and the backtracking
enumerator does not return. The box has

```
2,481,687,900  lattice points
```

and the elimination counts them in **8.4 s with 50,758 states**, because the cost is set by
the width of the frontier and not by the rank. Carrying the partition of the frontier as well
counts the ones with connected support: **909,103,210**, in 78 minutes and 9,749,212 states.
Both of Bio 8's empty counting columns are now filled. The box count is confirmed
independently by pinning one orbit variable to each of its values in turn and summing, for two
orbits with disjoint schedules.

| k | rank Lambda^- | box points | box states | connected | conn. states | % conn. | #DS |
|---|---|---|---|---|---|---|---|
| 13 | 0 | 1 | — | 1 | — | 100 | 1 |
| 12 | 4 | 10 | 4 | 8 | 8 | 80.0 | 10 |
| 11 | 12 | 528 | 20 | 218 | 136 | 41.3 | 864 |
| 10 | 48 | **2,481,687,900** | 50,758 | **909,103,210** | 9,749,212 | 36.6 | — |

The connected fraction falls with k, and the reason is visible: a point picks one arc from
each free `rho`-orbit, so at k = 10 the support is about 111 arcs on 125 vertices whatever the
point — sparse enough that a random choice usually falls apart. 36.6% do not.

## Two things the note settles, and one it does not

* **Connectivity was never a side condition.** A balanced digraph has strongly connected weak
  components, so its out-Laplacian has nullity equal to their number; deleting one row and
  column leaves it singular as soon as there are two. So `t(G_c) = 0` on a disconnected
  support by itself, and the connectivity predicate can be deleted from the lattice sum
  without changing it. This is the only form the "algebraic criterion for connected support"
  asked for in Bio 8 can take: the vanishing of a determinant whose entries are linear in `c`.
* **The genome's own coordinate is not the good elimination order.** Sweeping the chromosome
  keeps the frontier *bounded* — it is contained in the repeat classes straddling the
  coordinate — but at k = 10 that is 62 vertices wide and the elimination does not finish,
  while a greedy order on the constraint hypergraph is 26 wide and finishes in seconds. The
  reason is the interleaving of inverted repeats (Bio 1): pairwise crossing repeat intervals
  force every one of them to be open at the crossing point. Coordinate width is a measure of
  interleaving.
* **`#DS(phiX174, 10)` is still out of reach, and now for a stated reason.** The elimination
  sums any weight that factors into per-orbit and per-vertex terms (Prop. 5), and two of the
  three factors of `N_seq` do: `Loc(G_c)` is a product over vertices of a function of the
  out-degree, `prod c(a)!` a product over orbits. The third, `t(G_c)`, is a determinant of
  order up to `|V_c|` and does not factor over any bounded frontier. Point-by-point evaluation
  is cheap per point and hopeless in aggregate: 9.09e8 determinants of order up to 125. The
  wall is the arborescence factor, not the lattice — the index set is now known exactly, and
  it is the summand that cannot be carried. Bounds from Bio 8's Table 2:
  `8,610,708,632 <= #DS(phiX174,10) <= 1.368e36`.

## Traps recorded

* **A canonical class label is not a vertex index.** The connectivity state stores the
  partition of the frontier as a canonical label tuple; seeding the union-find with
  `par[q] = pt[q]` silently merges classes, because the vertex at index `pt[q]` generally
  carries a different label. Point each frontier vertex at the *first vertex carrying its
  label* instead. The symptom was a connected count of 0.
* **Finished components must be counted per class, not per vertex.** Iterating over the
  closing vertices and incrementing a counter once each makes every multi-vertex component
  look like several, so the "at most one finished component" rule rejects everything. Take the
  set of distinct roots first.
* The two together gave a clean, plausible, entirely wrong `0`; the check against Bio 8's
  published 8 and 218 is what caught them. A connectivity DP that returns 0 everywhere looks
  exactly like a correct DP on an instance with no connected points.
* **Frontier width must be read off the order, not off a run.** Reporting `max` of the widths
  observed before a run aborts understates the width of the order — at k = 10 the coordinate
  order shows 44 when abandoned at 1.5 million states and 50 when abandoned later, while its
  true width is 62. `frontier_widths()` computes it from the order alone.

## Cross-references

* the lattice sum, `Lambda^-`, its rank by the Hopf trace formula, and the k = 12 and k = 11
  values this note reproduces: Bio 8 (XXI);
* the double cover and its balance: Bio 6; its deck involution as an anti-automorphism: Bio 7;
* the interleaving of repeats, which is what the coordinate width measures: Bio 1.
