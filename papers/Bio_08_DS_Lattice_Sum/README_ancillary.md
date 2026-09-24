# Ancillary files for *Double-stranded genome assembly as a lattice sum, and the inversion-flip group*

Ruqing Chen, 20 September 2026. Bio 8 of the applied strand; numbered XXI in the series.

## Contents

| file | what it is |
|---|---|
| `anccft21.tex` / `.pdf` | the note |
| `ds_lattice.py` | the lattice sum: double cover with labelled arcs, `Lambda^-`, its rank, the box enumeration, `N_seq` at every point, the independent brute force; 9 checks |
| `ds_lattice_output.txt` | the transcript of one full run (about 2 s) |
| `ds_flips.py` | the inversion flip graph on all reconstructions; 4 checks |
| `ds_flips_output.txt` | its transcript (under 1 s) |
| `phix174_NC_001422.1.txt` | bacteriophage phiX174, RefSeq NC_001422.1, 5386 bp, plain sequence |

## Running it

```
python -u ds_lattice.py
python -u ds_flips.py
```

Dependencies: `numpy` (modular rank only). The scripts import `pgl3_building.py` from
`../Paper 10` (`bareiss_det`) and `ecoli_sandpile.py` from `../Bio 5` (`compacted_graph`,
`reduced_laplacian` for the single-strand and `|K(D_k)|` comparisons), so the folders must stay
siblings. `ds_flips.py` imports `ds_lattice.py`.

Environment knobs: `KS_ENUM=12,11,10` also attempts the box at k = 10 (rank 48; the
capacity-pruned backtracking produced no point in ten minutes on a loaded machine, so the
default is `12,11`); `BOX_CAP`, `BRUTE_CAP`, `POINT_FULL_MAX` bound the enumerations.

## The checks

`ds_lattice.py`

| key | claim | mode |
|---|---|---|
| D | `D_k` with labelled arcs is balanced; `rho` is a free involution on labelled arcs with reverse-complement contents; `m(rho a) = m(a)`; fixed contents have even multiplicity | exact |
| L | `rank Lambda^- = (|A| - |V| + 2 - #fixA - #fixV)/2` at k = 9, 10, 11, 12 (Hopf trace) | exact |
| B | every box point gives `c >= 0`, balanced, `c + rho c = m` (k = 12: 10 points; k = 11: 528) | exact |
| S | the spectra of `S` and `Sbar` are box points, and `N_seq` there equals the published single-strand count `|K(G_k)| Loc / prod m!` (2 at k = 12, 132/4 = 33 at k = 11) | exact |
| E | brute-force enumeration of closed walks with `usage + rho usage = m`, bucketed by spectrum, returns exactly the connected box points with exactly `N_seq` at each (10 and 864 reconstructions) | exact |
| N | at k = 11 the un-normalised sum `sum t(c) prod(d+ - 1)!` is 2756, the count is 864 | exact |

`ds_flips.py`

| key | claim | mode |
|---|---|---|
| A | every reconstruction is aperiodic as a cyclic content sequence (all weights `1/d(T)` are 1) | exact |
| F | every spectrum class is reached from `spec(S)` by flips alone; in fact the flip graph on all reconstructions is one component | exact |

Printed but not gated: the class-distance histograms, the same-spectrum flip distances
(4, 6, 8 only), the vertex pairs realising flips, and the reach when flips are restricted to
the inverted-repeat pairs of `S` (212 of 218 at k = 11).

## Conventions that matter

**Counting.** `N_seq(c) = t(c) prod_v (d_c^+(v) - 1)! / prod_a c(a)!` is the number of circular
sequences with spectrum `c`, each weighted by `1/d(T)` (rotational symmetry order); all
reconstructions here have `d(T) = 1`. The single-strand notes count in the arc-rooted convention
`t(m) Loc(m)` without the division; on one strand the divisor is a global constant, on the
lattice it is not.

**Contents.** After unitig contraction an arc is a maximal segment between consecutive repeated
k-mers on one strand; its content is the segment string. Contents with `m >= 2` are single
(k+1)-mers (a longer repeated segment would have repeated interior k-mers). `Lambda^-` lives on
contents, not on labelled arcs.

**Brute force.** Closed walks are generated rooted at their smallest content `a0` and weighted
`1/usage(a0)`; the walk length is exactly `sum(m)/2`. Bucket sums are asserted integral.

## Cross-references

* the double cover `D_k`, its balance and connectivity: Bio 6;
* the anti-automorphism theorem and the symmetric linking form that this note shows does not
  enter the count: Bio 7;
* the arc-rooted counting convention and phiX174 single-strand numbers: Bio 1;
* the p-fold decomposition without `rho` (polyploid phasing): Bio 9.
