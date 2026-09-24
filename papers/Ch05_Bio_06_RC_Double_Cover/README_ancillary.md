# Ancillary bundle for "Double strands without a new lemma" (`rc_double_cover.tex`, Bio 6)

`rc_double_cover.py` reproduces every number: 6 checks, 6 passing, **about 33 s**.
Imports `ecoli_sandpile.py` from the sibling `Bio 5` folder (which imports
`pgl3_building.py` from `Paper 10`); keep the folders as siblings. Run
`python rc_double_cover.py`. Archived output: `rc_double_cover_output.txt`.
`phix174_NC_001422.1.txt` is here; the *E. coli* genome and feature table are read from
`../Bio 5/` (or re-fetched from NCBI if missing).

## The object, and why no new lemma is needed

`D_k(S)` has as arcs the (k+1)-mer occurrences in **S and in S̄**. Then

> `mult_D(v) = mult_S(v) + mult_S(v̄)` and `d⁺(v) = d⁻(v) = mult_D(v)`,

so **`D_k` is balanced**, and the unitig lemma, the bundle-chain lemma, the canonicity of
the skeleton and the splitting theorem all apply to it unchanged. The branch vertices are
the k-mers repeated *when both strands are counted*, so the fast pipeline of Bio 5 runs on
`D_k` with the positions of S̄ appended to those of S.

The task that prompted this note expected a *bidirected* bundle-chain lemma producing
summands `(Z/2r)^(s−1)`. None is needed, and `2r` is the wrong multiplicity: a family with
`a` forward and `b` reverse copies has multiplicity **`a + b`** in `D_k` (Theorem 3.1).

## Checks

| key | content | arithmetic |
|-----|---------|-----------|
| (B) | `D_k` is balanced and multiplicities add across strands, on 4 instances (φX174 k = 10, 11, 12; *E. coli* k = 150) | exact |
| (C) | `D_k` is connected **iff** S has an inverted repeat of length ≥ k, with the witness printed: φX174 at k = 12 is glued by `ATTAAGCTCATT` at 441 / its reverse complement at 1448; *E. coli* at k = 150 by a 150-mer at 19796 / 3584046 | exact |
| (R) | ρ(v) = v̄ is an involution with `A(ρu, ρv) = A(v, u)`, i.e. an isomorphism `D_k → D_k^op` | exact |
| (M) | the seven rRNA operons: in `D_k` they carry multiplicity-**7** chains, and the longest one (span 728 bp at k = 150, 100, 51) has its seven occurrences mapping one-to-one onto the seven operons, the two minus-strand ones appearing on the reverse strand of `D_k`. The r = 7 chains carry more k-mers than any other multiplicity over the operons (4154 / 5254 / 6660 at k = 150 / 100 / 51) | exact |
| (K) | `K(D_k)` exactly for φX174 at k = 10, 11, 12, beside `K(G_k)` | exact |
| (E) | the *E. coli* double-cover table | exact / numeric |

## Results

φX174 (all exact; the single-strand rows reproduce Bio 3, Table 1):

| k | graph | branch | chains | skeleton | K(Sk) | length-driven |
|---|---|---|---|---|---|---|
| 12 | G_k | 2 | 1 | 1 | 0 | Z/2 |
| 12 | **D_k** | 7 | 0 | 7 | **Z/3 + Z/11** | 0 |
| 11 | G_k | 10 | 2 | 8 | Z/3 + Z/11 | (Z/2)² |
| 11 | **D_k** | 30 | 7 | 23 | (Z/2)² + Z/49451 | (Z/2)⁷ |
| 10 | G_k | 47 | 8 | 37 | Z/4 + Z/1076338579 | (Z/2)¹⁰ |
| 10 | **D_k** | 125 | 23 | 97 | (Z/4)² + Z/79624452132275189516228869 | (Z/2)²⁸ |

**The phage's single global bit is a single-strand phenomenon.** `K(G_12) = Z/2` — the two
reconstructions Bio 2 exhibited — but `K(D_12) = Z/33`, and it is entirely arrangement:
`D_12` has no bundle chain at all.

*E. coli* K-12, double cover:

| k | branch | chains | skeleton | log₁₀ length-driven | log₁₀\|K(Sk)\| | log₁₀\|K(D_k)\| |
|---|---|---|---|---|---|---|
| 150 | 39 768 | 278 | 312 | 20 987.91 | 156.51 | 21 144.42 |
| 100 | 44 808 | 396 | 434 | 23 394.31 | 208.92 | 23 603.24 |
| 51 | 53 380 | 798 | 886 | 27 405.96 | 395.80 | 27 801.75 |

On one strand the same k give skeletons of 115, 175, 366 and arrangement exponents 49.71,
71.11, 143.17 (Bio 5). So the skeleton roughly triples and so does the irreducible
exponent: **inverted repeats are a large part of this genome's irreducible ambiguity, and
the single-strand theory sees none of them.**

## What this is *not*

`D_k` is the **double cover**, not the bidirected quotient, and the difference matters:

* `K(D_k)` does **not** count double-stranded reconstructions. An Eulerian circuit of
  `D_k` traverses all 2N arcs in one closed walk, which is not a genome; a double-stranded
  molecule is a decomposition of `D_k` into two ρ-conjugate circuits, and counting those
  is not solved here.
* No sandpile theory for the bidirected quotient is offered. The obstruction Bio 1 and
  Bio 2 recorded — reverse complementation reverses arc direction, so the quotient map is
  not a covering of digraphs and Reiner–Tseng does not apply — is untouched.
* ρ gives a perfect pairing on `K(D_k)` (Proposition 5.1, via Bowen–Franks duality and
  `Δ^op = Δᵀ`), but we do not determine whether it is symmetric or alternating, nor that
  it is independent of the sink. What is settled is the negative half of Bio 1's warning:
  ρ reverses arcs, so it does not act on `K(D_k)` at all and cannot split it into
  ±1-eigenspaces; no 2 is inverted anywhere.
