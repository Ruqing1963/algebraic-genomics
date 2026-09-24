# Ancillary files for *Polyploid phasing as a p-fold decomposition of a balanced flow, and the homogeneous space of alleles*

Ruqing Chen, 20 September 2026. Bio 9 of the applied strand; numbered XXII in the series.

## Contents

| file | what it is |
|---|---|
| `anccft22.tex` / `.pdf` | the note |
| `polyploid_phasing.py` | the pooled compacted graph, the phasing lattice `Z_1 ∩ [0, m]`, `N_seq` at every point, the p-fold decomposition sums, the independent walk enumerations; 6 checks |
| `polyploid_phasing_output.txt` | the transcript of one full run (under 1 s) |
| `phix174_NC_001422.1.txt` | bacteriophage phiX174, RefSeq NC_001422.1, 5386 bp, plain sequence |

## Running it

```
python -u polyploid_phasing.py
```

Dependencies: none beyond the standard library, except that the script imports
`ds_lattice.py` from `../Bio 8` (`load_phix`, the check-row printer) and `pgl3_building.py`
from `../Paper 10` (`bareiss_det`), so the folders must stay siblings. Knobs: `BOX_CAP`,
`WALK_CAP`.

## The six checks

| key | claim | mode |
|---|---|---|
| P | the pooled compacted graph is balanced and gains exactly two branch vertices per substitution site away from repeats | exact |
| L | b unlinked biallelic sites give `2^b` equal-length ordered decompositions, `2^b N_seq(S)^2` labelled and `2^(b-1) N_seq(S)^2` unlabelled phasings (b = 1, 2, 3, 5); two sites within k of each other are one site (4 sites, 3 blocks → 8; `{600,606}` vs `{600}` → 2) | exact |
| W | p = 2: the walk enumeration of `T_1` (rooted at its smallest content, weight `1/usage(a0)`, usage ≤ m) reproduces `N_seq` on every spectrum, and `sum_{T_1} N_seq(m − spec T_1)` equals the lattice sum | exact |
| U | p = 2: unordered phasings enumerated as sets of explicit walk pairs `{T_1, T_2}` equal the multiset formula `prod_c binom(N(c)+r_c−1, r_c)` | exact |
| T | p = 3: ordered decompositions = `prod_sites binom(3, a_site)`, `3!` at a triallelic site (54, not `|Hom(pi_1, S_3)| = 216`); with two identical haplotypes the unlabelled count is `N binom(N+1, 2) = 6`, not `24/6` | exact |
| C | without the equal-length constraint the pool admits copy-number decompositions through the 13-mer repeat (22 of 24 for one site; 16 of 17 for `S+S`), including tandem-dimer circles of weight 1/2; the equal-length unlabelled count of `S+S` is `binom(N_seq+1, 2) = 3` | exact |

## Core conclusions

* **Labelled phasings = a lattice sum.** `#Lab = sum over ordered decompositions m = c_1+…+c_p
  (connected balanced box points of Z_1) of prod N_seq(c_i)`, with
  `N_seq(c) = t(c) prod(d+−1)! / prod c(a)! = sum_T 1/d(T)`. Same object as Bio 8 with the
  strand involution removed.
* **S_p acts on decompositions, not on the graph.** The map from the disjoint union of the
  haplotype graphs to the pool is a covering iff no k-mer is shared — i.e. never when there is
  anything to phase. No deck group, no monodromy.
* **The local fibre is a homogeneous space.** At a site with allele counts `(a_1, a_2, …)` the
  choice is an ordered partition of the labels: `S_p / prod S_{a_i}`, of size the multinomial
  coefficient. Unlinked sites on a repeat-free reference multiply; `|Hom(pi_1, S_p)| = (p!)^b`
  is wrong for p ≥ 3.
* **Linkage.** Sites inside one (k+1)-window are one multi-allelic site whose alleles are the
  distinct local haplotypes.
* **Unlabelled = multiset count**, `= #Lab / p!` exactly when all parts are pairwise distinct.
* **Free lengths expose copy-number ambiguity.** `S + S` has the same pooled spectrum as a
  2496-bp sub-circle plus an 8276-bp circle, and as two tandem dimers of complementary
  sub-circles (weight 1/2 each).

## Conventions that matter

**Compaction.** A k-mer is a branch vertex iff it has ≥ 2 distinct successors or ≥ 2 distinct
predecessors in the pool. Bio 5's "branch ⟺ repeated" is useless in a pool (every shared k-mer
is repeated). `N_seq` is intrinsic and does not depend on the compaction.

**Fractions.** `N_seq(c)` is a `Fraction`; it is a non-integer exactly when every sequence with
spectrum `c` is periodic (a tandem-dimer circle). The multiset formula is applied only to
decompositions without such parts; the script reports the rest separately.

**Equal length.** The substitution-only model is the restriction `|c_i| = N` for all parts;
`|c| = sum c(a) ell(a)` with `ell(a)` the number of (k+1)-mers of the content.

## Cross-references

* the lattice sum, `N_seq` and its `1/d(T)` weighting: Bio 8 (XXI);
* the phiX174 single-strand decomposition at k = 12 (`a X1 b Y1 a X2 b Y2`, `N_seq = 2`): Bio 1;
* the single-genome compaction criterion this note replaces: Bio 5.
