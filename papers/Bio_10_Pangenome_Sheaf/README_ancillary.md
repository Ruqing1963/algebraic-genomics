# Ancillary files for *The pan-genome sheaf sees nothing: the recombination group of a pan-genome graph, and its torsion*

Ruqing Chen, 21 September 2026. Bio 10 of the applied strand; numbered XXIII in the series.

## Contents

| file | what it is |
|---|---|
| `anccft23.tex` / `.pdf` | the note |
| `pangenome_sheaf.py` | the pan-genome sheaf and its Laplacian; the recombination group `R = Z_1(G)/sum Z_1(G_i)` by Smith form; the nerve of intersection components and `H_1`; the surface pools; the GenBank isolates; 5 checks |
| `pangenome_sheaf_output.txt` | the transcript of one full run (7 s) |
| `real_pangenome.py` | the same group on **complete E. coli chromosomes**: pooled branch vertices by lexsort-hashing with a brute-force audit, unitigs, the nerve of components, and a freeness certificate for `R`; 3 checks. **Standalone** — numpy only, no other file of the series |
| `real_pangenome_output.txt` | its transcript |
| `real_pangenome_p4.txt` .. `_p8.txt` | the four- to eight-strain transcripts |
| `RUN_LOCALLY.md` | how to run it on another machine, and what the flags do |
| `pangenome/*.fasta` | seven complete *E. coli* chromosomes from NCBI (O157:H7 Sakai, IAI39, K-12 DH10B, CFT073, S88, 536, UTI89); K-12 MG1655 comes from Bio 5 |
| `phix174_NC_001422.1.txt` | bacteriophage phiX174, RefSeq NC_001422.1, 5386 bp |
| `phix174_isolates.fasta` | the GenBank records returned by the E-utilities query `phiX174[All Fields] AND 5300:5500[SLEN] AND biomol_genomic[PROP]`, cached on first run |

## Running it

```
python -u pangenome_sheaf.py
```

Dependencies: `numpy`, `sympy` (Smith normal form over **Z**). Imports `Pool` and `mutate` from
`../Bio 9/polyploid_phasing.py`, the check-row printer from `../Bio 8/ds_lattice.py`, and
`pgl3_building.py` from `../Paper 10`; the folders must stay siblings. Without network access
check (D) is skipped (the cache file, if present, is used).

## The five checks

| key | claim | mode |
|---|---|---|
| F | the pan-genome sheaf (stalk `Z^{strains through v}`, projections) is `(+)_i Z_{G_i}`; `coker(sheaf Laplacian) = (+)_i (Z + K(|G_i|))`, verified by Smith form of the assembled Laplacian (`Z/38 + Z/38 + Z^2`; `(Z/84)^3 + Z^3`) | exact |
| B | two strains: `R` is free of rank `c(G_1 ∩ G_2) − 1` (8 pools; on a repeat-free reference with b bubbles, b − 1; on phiX174 the two-copy repeat glues shared stretches) | exact |
| S | six strains with the hemi-icosahedron (RP²) as nerve give `R = Z/2`; octahedron 0; 7-vertex torus `Z^2`; the nerve is verified to be exactly the designed triangulation | exact |
| T | twelve random three-strain pools: Smith form and nerve agree | exact |
| D | six phiX174 GenBank genomes, strand- and origin-normalised: `R = Z^38`, both ways | exact |

`real_pangenome.py`

| key | claim | mode |
|---|---|---|
| A | the hashed branch-vertex criterion agrees with brute force on random 20 kbp fragment pools — the candidate pass loses no true branch vertex | exact |
| V | the **F**<sub>p</sub> torsion detector reproduces Bio 10's exact Smith-form answers on the three surface pools, **including the Z/2 of the projective plane** — so a null result on real data is a result | exact |
| E | complete *E. coli* chromosomes at k = 31, each pool adding one strain: `Z^85`, `Z^29632`, `Z^66447`, `Z^87858`, `Z^104495`, `Z^116056`, `Z^117848`. All seven **certified free** by a unit-pivot reduction of `d_2` (every elementary divisor 1 ⟹ `im d_2` a direct summand ⟹ `R` free at every prime, so no prime-by-prime testing is needed at all). The two K-12 strains share 86 blocks and give `86 − 1`, which is Corollary 5 on real data; the eight-strain nerve has 744,598 one-cells and 1,683,865 two-cells, and the rank increments fall from 29,547 to 1,792 | exact |

Run it with `--strains N -k 31 --out FILE`; no argument at all runs the checks and all eight.
See `RUN_LOCALLY.md`. Times at k = 31 on one core, cold: 26 / 40 / 60 / 84 / 110 / 140 / 182 s
for 2 to 8 strains, of which the unit-pivot reduction is 0 to 28 s; most of the rest is the
branch-vertex pass over 40 Mbp.

## Core conclusions

* **The sheaf is trivial.** Projections never mix strains, so the sheaf is a sum of constant
  sheaves on the strain subgraphs; its cohomology and Laplacian cokernel are lists of per-strain
  quantities.
* **The recombination group.** `R = Z_1(G)/sum_i Z_1(G_i)` ≅ `H_1` of the nerve of intersection
  components (Mayer–Vietoris spectral sequence, two rows). Rank = number of independent
  strain-switching cycles.
* **Two strains: free**, rank `c(G_1 ∩ G_2) − 1`.
* **Three or more: torsion is possible**, and realised: a `Z/2` class is a mosaic cycle whose
  double, but not itself, is a sum of within-strain cycles. Invisible to every strain and every
  pair.
* **Real isolates:** six phiX174 genomes, `R = Z^38`, no torsion.
* **Eight complete chromosomes:** `R = Z^117848`, free at every prime with a certificate. The
  nerve of a real pan-genome is wide and sparse — each pair of divergent strains shares its
  core in tens of thousands of separate blocks — which is the opposite of the tight
  triangulated surface a torsion class needs.
* **The rank saturates.** Increments 29,547 / 36,815 / 21,411 / 16,637 / 11,561 / 1,792: each
  new chromosome after the third contributes less than the one before. `rank R` measures the
  diversity of a strain set, not its size. An observation on these eight genomes, not a law —
  the last four additions are all extraintestinal pathogenic strains.

## Traps recorded

* Fundamental cycles: with `parent[y] = (x, a, s)` and `s = +1` iff `a` is the arc `x → y`, the
  step `y → x` contributes `−s`. The opposite sign produced "rank sum Z_1(G_i) > rank Z_1(G)" and
  spurious `(Z/2)^22`; the code now asserts `∂x = 0` for every basis vector.
* **Pivot ordering decided a mathematical conclusion.** `unit_reduce` and `sparse_rank_mod`
  used to pick each pivot by `sorted(alive, key=len)` — an `O(n log n)` scan of the whole live
  column set *per pivot*, which at 278,000 columns was essentially the entire cost of the
  script (5-strain reduction: 14,080 s). Both now take the shortest column from a lazy heap:
  4 s for the same matrix, same rank. The better ordering also produces less fill-in, and at
  4 and 5 strains the entries then never leave `{0,±1}`, so the reduction *completes* where it
  had stalled — turning "no 2- or 3-torsion" into "free at every prime". A failed unit-pivot
  certificate says something about the elimination, not about the matrix.
* **A quadratic lookup hiding behind a generator expression.** `boundary_cols` took each
  component's representative with `next(v for v, cc in comp[s].items() if cc == c)` — one
  linear scan of the whole intersection per component. It reads like O(1) and is O(|G_s|); at
  six strains it outweighed everything else in the script. One pass now builds the table, and
  the 3-strain nerve step went 95 s → 1 s. The tell was that the "union-find" time grew far
  faster than the number of arcs.
* Random 40-mer segments at k = 12 share accidental 12-mers (≈ 10⁴ k-mers per pool); the surface
  pools use k = 16.
* One GenBank record is deposited rotated (4018 "substitutions" before normalisation); the
  isolates are rotated to the reference origin and, if needed, reverse-complemented before pooling.

## Cross-references

* the pooled graph, its compaction criterion and the phasing lattice: Bio 9 (XXII);
* the per-sequence count and the closed-walk lattice sums: Bio 8 (XXI);
* the isoform polytope (the same flow-decomposition object on a DAG): Bio 11 (XXIV).
