# Ancillary bundle: the bundle-chain decomposition of *Escherichia coli* K-12 MG1655

`ecoli_sandpile.py` reproduces every number below: 6 checks, 6 passing, **about 13 s**,
peak memory under 0.5 GB. Requires `numpy` and `sympy`, and imports `pgl3_building.py`
from the sibling `Paper 10` folder (exact p-adic Smith normal form); keep the folders as
siblings. Run `python ecoli_sandpile.py`. Archived output: `ecoli_sandpile_output.txt`.

**Data**, all cached in this folder so no network access is needed (the script re-fetches
the first two from the NCBI E-utilities API if they are missing):

* `NC_000913.3.fasta` — the complete genome, 4 641 652 bp, circular, pure A/C/G/T
* `NC_000913.3.ft.txt` — the NCBI feature table (4390 features parsed: rRNA, CDS,
  mobile_element, repeat_region), used to annotate the repeat families
* `phix174_NC_001422.1.txt` — the phage genome of the companion chapters, used by
  check (V) to validate this faster pipeline against their published numbers

## Why a 4.6 Mbp genome is tractable at all

In the convention of Bio 1–4 (vertices = distinct k-mers, one arc per (k+1)-mer
*occurrence*) every vertex satisfies d⁺(v) = d⁻(v) = mult(v). Hence

> a vertex survives unitig contraction **iff its k-mer is repeated**.

So the 4 641 652-vertex graph is never built. The script finds the repeated k-mers, and
the compacted graph is read off directly: one vertex per repeated k-mer, one arc per
consecutive pair of repeat occurrences around the circle. That is 24 k–34 k vertices
instead of 4.6 M. Bundle-chain contraction (Bio 3, Lemma 2) then runs in one pass, because
Bio 3's Lemma 6 says the chain steps form a partial injection, so the maximal chains are
disjoint and the contraction order is irrelevant.

## Checks

| key | content | arithmetic |
|-----|---------|-----------|
| (V) | **the fast pipeline against the published φX174 decomposition.** This code never builds G_k, contracts chains in one pass rather than iteratively, and takes the Smith form p-adically rather than by Euclidean elimination — three differences from the implementation the companion chapters used. Run on φX174 it returns, at k = 9, 10, 11, 12, the same branch counts (2, 10, 47, 146), chain counts (1, 2, 8, 28), skeleton sizes (1, 8, 37, 110), length-driven exponents (1, 2, 10, 36) and the same exact skeleton groups, including Z/4 + Z/1076338579 at k = 10 and all seventeen elementary divisors at k = 9. With a 4.6 Mbp target that admits no direct check, agreement between independent implementations on a genome that does is the only protection available | exact |
| (H) | k-mer multiplicities: positions are hashed to 64 bits and `np.unique` finds the repeated hashes; the actual k-mer **strings** at those positions then decide, so a hash collision can only add a spurious candidate, never a spurious repeat. 0 collisions occurred at any k | exact |
| (C) | the compacted graph is balanced with d⁺(v) = mult(v) (checked on up to 2000 vertices per k) | exact |
| (B) | maximal bundle chains found in one pass; the injectivity of the chain-step relation (Bio 3 Lemma 6) is asserted at run time on all 1541 chains | exact |
| (K) | skeleton group: exact invariant factors by p-adic Smith form (`pb.coker_structure`, Bareiss determinant + Z/pᵏ elimination) for skeletons of ≤ 200 vertices; log₁₀ of the order by float LU (`slogdet`) otherwise | exact / numeric |
| (R) | no rRNA-bearing chain ever has multiplicity 7; the largest is 5 (the + strand operons), and at every k the longest chain in the whole genome is the multiplicity-2 spine of the two co-oriented minus-strand operons | exact |

## Results

| k | branch vertices | bundle chains | skeleton | K(Sk) | log₁₀(length-driven) | log₁₀\|K(G_k)\| | log₁₀ #reconstructions |
|---|---|---|---|---|---|---|---|
| 150 | 23 933 | 107 | **115** | (Z/2)¹⁷+(Z/3)⁹+(Z/4)¹²+(Z/5)¹⁰+(Z/7)³+(Z/8)⁴+(Z/9)³+(Z/16)³+Z/27+(Z/32)²+Z/49+Z/21308303 | 10 439.18 | 10 488.89 | 20 735.14 |
| 100 | 26 631 | 166 | 175 | (Z/2)²⁹+(Z/3)¹⁶+(Z/4)¹⁵+(Z/5)¹⁵+(Z/7)²+(Z/8)⁴+(Z/9)²+(Z/16)⁶+Z/25+Z/27+(Z/32)²+Z/49+Z/81+Z/128+Z/512+Z/3448649 | 11 619.47 | 11 690.57 | 23 310.20 |
| 75 | 28 579 | 243 | 261 | order 10^99.80 | 12 463.28 | 12 563.08 | 25 132.01 |
| 51 | 31 011 | 342 | 366 | order 10^143.17 | 13 573.65 | 13 716.83 | 27 624.76 |
| 31 | 34 272 | 683 | 870 | order 10^393.75 | 14 956.79 | 15 350.54 | 31 748.76 |

`#reconstructions = |K(G_k)| · ∏_v (d⁺(v)−1)!` is the BEST-theorem count of circular
sequences with exactly this k-mer spectrum (as arc-rooted Eulerian circuits).

Two readings. First, **the whole irreducible ambiguity of a 4.6 Mbp genome is a
115-vertex digraph** at k = 150: 4 641 652 → 23 933 (unitig) → 115 (bundle chains).
Second, the split is extremely lopsided — at k = 150 the length-driven part is 10^10439
and the arrangement part only 10^49.7, i.e. **99.5 % of the exponent is removable by
longer k-mers and 0.5 % is not removable by any k**. The arrangement part nevertheless
grows by a factor 10^344 as k falls from 150 to 31.

## The rRNA operons: the brief's premise is wrong, and the theory says why

E. coli K-12 has seven rRNA operons, and the natural expectation is a repeat family of
multiplicity 7. **It does not exist, at any k.** The feature table places

* five operons on the + strand: 223771, 3941808, 4035531, 4166659, 4208147
* two on the − strand: 2726069, 3423423

The theory of Bio 1–4 is single-strand: two copies are a repeat of G_k only if they are
**co-oriented in S**. So the seven operons cannot form one family; they split by strand,
and the computation confirms it exactly:

* the **longest bundle chain in the whole genome**, at every k, is the multiplicity-**2**
  spine of the two minus-strand operons — s = 3881 at k = 150 (span 4030 bp), s = 4412 at
  k = 51 (span 4462 bp), annotated 16S + 23S + 5S ribosomal RNA;
* the five + strand operons give the multiplicity-**5** family, whose longest chain is the
  16S gene (s = 626, span 775 bp at k = 150; s = 745, span 775 bp at k = 31), with copies
  at 223997, 3942034, 4035757, 4166885, 4208373 — all five operons, confirmed 5/5;
* they are not identical to one another either, so the rRNA region also throws off chains
  of multiplicity **3 and 4**, one per subset of copies that happen to agree.

Multiplicity 7 does occur in the genome — but it is **IS5A**, not rRNA (5 chains, 88
k-mers at k = 31), as is multiplicity 8 and 9. The other long chains are all insertion
sequences: IS186A (×3, 1345 bp), IS2D (×3, 1331 bp), IS3B (×2, 1260 bp), IS3A (×3,
1255 bp), IS5U (×3, 1202 bp). At k = 31 the top multiplicities are 23 and 22 with spans of
32–35 bp: the REP/BIME extragenic palindromic elements, which vanish completely by k = 51.

## Scope

Single-strand only. A repeat whose copies lie on **opposite** strands is invisible here —
which is exactly why the seven operons appear as 5 + 2 rather than 7, and why the IS
element copies that sit in inverted orientation are not counted in the multiplicities
above. A bidirected theory would start from the signed Laplacian and does not exist
(Bio 2, Open (a)); this computation is a quantitative measure of how much it would change.

Exactness. Every structural number (branch vertices, chains, (r,s) multisets, skeleton
size, positions) is exact integer work. K(Sk) is exact at k = 150 and k = 100; at k = 75,
51, 31 only log₁₀ of its order is given, by float LU, because the Bareiss determinant of a
260–870 vertex skeleton is out of reach and a 100–400-digit determinant could not be
factored anyway. The two exact groups each carry one large prime factor we do not
interpret (21 308 303 at k = 150, 3 448 649 at k = 100); as in Bio 3 that is what the
sandpile group of a generic digraph looks like, and no topological reading is claimed.
