# Ancillary files for *The isoform decomposition polytope of a splicing graph, and a BEST identity for its integer points*

Ruqing Chen, 21 September 2026. Bio 11 of the applied strand; numbered XXIV in the series.

## Contents

| file | what it is |
|---|---|
| `anccft24.tex` / `.pdf` | the note |
| `isoform_polytope.py` | splicing DAGs, the path–arc matrix and its rank, the fibre polytope `Theta(w)` (integer points, vertices, minimum decompositions), the augmented multigraph and the BEST identity, the Dscam rank, a random search; 7 checks |
| `isoform_polytope_output.txt` | the transcript of one full run (about 3.5 minutes) |

## Running it

```
python -u isoform_polytope.py
```

Dependencies: `numpy` (modular rank of the Dscam matrix). Imports `pgl3_building.py` from
`../Paper 10` (`bareiss_det`) and `ds_lattice.py` from `../Bio 8` (the check-row printer), so
the folders must stay siblings.

## The seven checks

| key | claim | mode |
|---|---|---|
| R | the path–arc matrix of a splicing DAG (every arc on an s–t path) has rank `|E| − |V| + 2` | exact |
| D | `dim Theta(w) = #paths(supp w) − (|E_w| − |V_w| + 2)`; for events in series `prod k_i − sum(k_i − 1) − 1` | exact |
| U | deconvolution from junction counts is unique for a single alternative event and already of dimension 1 for two cassette exons in series | exact |
| S | Dscam (12 × 48 × 33 × 2): 38 016 isoforms, rank 92 (190 × 38016 matrix mod 2³¹−1), polytope dimension 37 924 | exact mod p |
| V | every minimum-support decomposition (real or integer) is a vertex of `Theta(w)`; vertices can be fractional for integer `w` (three parallel-arc pairs, `w ≡ 1`: four weights 1/2 beside an integer minimum of two paths) | exact |
| B | BEST identity: `N_seq(G_w + F return arcs) = sum over integer points of (F−1)! / prod theta_j!` on 4 chosen flows and 232 random ones | exact |
| I | real minimum vs integer minimum on 232 random small flows: never different; one instance with a fractional vertex | exact |

## Core conclusions

* **Dimension.** `dim Theta(w) = |P(w)| − rank`, rank = `|E_w| − |V_w| + 2`, because `w` is in the
  relative interior of the cone of its support paths. Unique iff the support paths are linearly
  independent.
* **Extreme rays ≠ vertices ≠ lattice points.** Rays of the flow cone are paths; vertices of
  `Theta(w)` are independent-support decompositions and contain every minimum decomposition;
  they need not be integral (a 3 × 3 path–arc submatrix of determinant 2 gives halves).
* **BEST identity.** Closing the flow with `F` return arcs, `t · Loc / (F! prod w(a)!)` equals the
  multinomial-weighted count of integer decompositions, i.e. the number of cyclic words in `F`
  isoforms with the given usage. A determinant for the weighted partition function; none for the
  plain point count.
* **Parsimony is not truth.** The tropomyosin-like flow was generated from 5 isoforms; its minimum
  decomposition has 4 paths.

## Conventions

`N_seq` is the per-sequence count of Bio 8 (`sum_T 1/d(T)`); a cyclic word of `F` isoforms can
be periodic (all `theta_j` sharing a common divisor), so the identity is stated and checked with
the `1/d` weighting, in exact rational arithmetic. `Loc(G) = prod_v (d+(v) − 1)!` over
vertices of positive out-degree; `t` = arborescences converging on a root (Bareiss determinant of
the reduced out-Laplacian). Integer points are enumerated by backtracking over the support
paths; vertices by solving every independent support of size ≤ rank exactly over `Fraction`.

## Cross-references

* the per-sequence normalisation and the closed-walk lattice sums: Bio 8 (XXI), Bio 9 (XXII);
* the arc-rooted single-genome convention: Bio 1.
