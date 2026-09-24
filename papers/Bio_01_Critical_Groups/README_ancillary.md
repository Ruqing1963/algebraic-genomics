# Ancillary bundle for "Assembly ambiguity is not determined by branch statistics: critical groups of de Bruijn graphs" — version 2 (`dbg_sandpile.tex`)

`dbg_sandpile.py` reproduces every computational claim of the paper: 30 checks, 30 passing,
about two minutes. Requires `sympy` (also reached through `pgl3_building.py` in the `Paper 10`
folder of the ANCFT series, which supplies the exact integer Smith normal form and Bareiss
determinant); the path is resolved relative to this file, so keep the two folders as
siblings. Run `python dbg_sandpile.py`. The archived output is `dbg_sandpile_output.txt`.

**Data.** `phix174_NC_001422.1.txt` is the 5386-bp circular genome of bacteriophage phiX174
(GenBank NC_001422.1), one uppercase A/C/G/T string, no header. It is read from beside the
script; no network access is needed to reproduce Table 1. All other sequences are generated
inside the script from fixed seeds. The two witness sequences of Theorem 8 are hard-coded
as `W_LO` and `W_HI`.

Checks, keyed as in the paper's battery table (Table 2):

| key | content | arithmetic |
|-----|---------|-----------|
| (T) | de Bruijn tower = line-digraph tower: DB(s,k+1) = L(DB(s,k)) by explicit construction; arborescence count s^(s^k - k - 1); the [ANCFT III] recursion at each step | exact |
| (L) | independent validation against Levine, JCTA 118 (2011) Thm 1.3: K(DB(2,n)) factor by factor, n = 2..5 | exact |
| (E) | BEST: #arc-rooted Eulerian circuits = \|K\| x prod (d-1)!, against exhaustive enumeration on instances with <= 14 arcs | exact |
| (R) | Theorem 5: one repeat of length L, multiplicity r -> (Z/r)^(L-k), 54 parameter settings; instances where the background makes its own repeat are detected by the degree multiset and skipped | exact |
| (I) | **Theorem 6 (new in v2)**: two repeats, disjoint = nested = (Z/2)^(2s-2), interleaved = (Z/2)^(2s-1). Verified three ways: (a) every clean random sequence instance matches the prediction; (b) the linear-time compaction produces exactly 2s branch vertices with connector multiplicities {1,1,1,1} (disjoint/nested) and {2,2} (interleaved); (c) symbolic Smith form of the explicit reduced Laplacians of the proof for s = 1..7 | exact |
| (M) | **Proposition 7 (new)**: three repeats of length L, multiplicity 2, in all 90 circular words with two of each letter. At s = 1 the compacted graph is the connector digraph D_w, so K(G_L) = K(D_w): 0 / Z/2 / Z/3 / (Z/2)^2 for no / one-pair / path / triangle crossings — the three non-trivial patterns all have F_2-rank 2, so the correction is NOT a chord-diagram rank (and can have odd order). At s = 2 the group is (Z/2)^3 + K(D_w). Groups are compared by elementary divisors, since (2,2,6) and (2,2,2,3) present the same group | exact |
| (N) | Theorem 8: recorded witnesses with identical (k, n_k, degree multiset, local factor), \|K\| = 2 vs 128; plus the randomized degree-class splitting rate | exact / sampled |
| (C) | collapse order k*: statistics over random sequences; then an EXHAUSTIVE search over binary sequences of length <= 10 for the minimal witnesses of k* = R, k* = R+1, k* < R | sampled / exact |
| (P) | Lemma 4: the linear-time unitig contraction leaves K unchanged (12/12, reductions up to 68 -> 2 vertices); self-loops leave the out-Laplacian unchanged (1, 2, 5 loops) | exact |
| (F) | **Table 1 (new in v2)**: the critical-group spectrum of phiX174 on the compacted graph, k = 10..13; R = 12 by binary search; k* = 13 = R+1 | exact |

Conventions. Sequences are circular; G_k(S) has one vertex per distinct k-mer and one arc
per (k+1)-mer occurrence, so out-degree = number of occurrences (arcs counted with
multiplicity, including parallel arcs inside a repeat). Reconstruction counts are
ARC-ROOTED (Remark 1 of the paper). The critical group is the cokernel of the transposed
out-Laplacian with the row and column of vertex 0 deleted.

Two implementation notes. `compact()` is linear in the number of arcs (one pass for
branch vertices, one walk per outgoing arc) and returns `(0, [])` for a single cycle;
it replaces the O(n^2) version of v1 and is what makes check (F) take milliseconds. The
exact-SNF routine is dense, so check (F) caps the compacted graph at 70 vertices and
reports the levels k <= 9 (146 to 1384 branch vertices) as skipped rather than computing
them slowly; a sparse or p-adic Smith form would reach them.

Provenance discipline. Rows marked "sampled" report statistics of a randomized search;
every other row is a deterministic exact-integer fact. Rows (L) and (E) are confrontations
with something outside the script — a published theorem and a brute-force enumeration.

