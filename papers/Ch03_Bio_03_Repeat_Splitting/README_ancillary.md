# Ancillary bundle for "The sandpile group of a de Bruijn graph splits along its repeats" (`repeat_splitting.tex`, Bio 3)

`repeat_splitting.py` reproduces every computational claim: 7 checks (six keys, (E) has
two rows), 7 passing, about 30 s. Requires `sympy`, reached through `pgl3_building.py` in the `Paper 10` folder of the
ANCFT series (exact integer Smith normal form); the path is resolved relative to this file,
so keep the folders as siblings. Run `python repeat_splitting.py`. The archived output is
`repeat_splitting_output.txt`. `phix174_NC_001422.1.txt` (5386 bp) is read from beside
the script; no network access is needed.

**The theorem.** In a balanced multidigraph, a bundle chain v_1 -> ... -> v_s (every v_t
sends all r of its out-arcs to v_{t+1}) contracts to a point at the exact cost of a direct
summand (Z/r)^(s-1). Iterating gives K(G) = (+)_chains (Z/r_c)^(len_c - 1) (+) K(skeleton).

Checks, keyed as in the paper's battery table:

| key | content | arithmetic |
|-----|---------|-----------|
| (A) | clean multi-repeat sequences: m in {2,3,4} repeats of multiplicities r_i in {2,3,4} and lengths 10..13, k in [min L - 3, min L]; K(G_k) computed directly equals (+)(Z/r_i)^(s_i-1) (+) K(D_w), where D_w is the connector digraph of the circular word of occurrences. 120 clean instances, 0 mismatches. Cleanliness is enforced by comparing the out-degree multiset with the one the planted repeats predict | exact |
| (B) | the lemma on arbitrary balanced digraphs: a random Eulerian digraph H (union of random closed walks, strongly connected), a vertex of out-degree r in 2..6 replaced by a bundle chain of length s in 2..4; K(G) equals (Z/r)^(s-1) (+) K(H). 60/60 | exact |
| (C) | the nine published values (m = 1, 2, 3 from Bio 1 v2) are the formula evaluated on D_w | exact |
| (D) | r_i = 2, all circular words with m letters each twice (m = 3, 4, 5: 30, 630, 22680 words): K(D_w) is NOT a function of the F2-rank of the unsigned chord intersection matrix (for m = 3 all three non-trivial patterns have rank 2; for m = 5 rank 2 carries nine groups), NOT of its integer Smith form ((1,1,0) carries Z/2 and Z/3). It is constant on isomorphism classes of the interlace graph for m <= 5 (4, 11, 34 classes, zero collisions) but NOT at m = 6: Merino-Moffatt-Noble's pair 123456132465 / 123654132564 (same interlace graph, Z/3+Z/6 vs Z/18) is recomputed from D_w, plus two further pairs. What determines the group is the SIGNED interlace matrix: K(D_w) = coker(A + I) (MMN Thm 7.3) on all 24543 words tested (m <= 5 exhaustive, 300 random each for m = 6..9) | exact |
| (E) | canonical skeleton: on 100 random balanced digraphs with planted chains, the maximal chains are pairwise disjoint and the chain partition, the (r, len) multiset and the skeleton's group are invariant under random relabelling; the skeleton is chain-free. Repeat graph: on 200 random circular words (m <= 5, r_i <= 4), K(RG_w) = (+)_i Z/r_i (+) K(D_w) where RG_w has r_i parallel arcs h_i -> t_i and one arc t_i -> h_j per adjacency | exact |
| (F) | phiX174: unitig-compact, then bundle-contract, at k = 12, 11, 10, 9, 8. The decomposition equals the directly computed group at k = 10, 11, 12; at k = 9 the compacted graph (146 branch vertices) is beyond the dense exact SNF but its skeleton (110) is not, extending the spectrum by one level; at k = 8 the skeleton (312) is still too large | exact |

Conventions. Groups are compared by ELEMENTARY DIVISORS (`elementary()`), never by raw
SNF diagonals — (2,2,6) and (2,2,2,3) are the same group. `bundle_contract()` finds a
vertex all of whose out-arcs go to one other vertex of the same out-degree, extends the
chain forward and backward maximally, contracts it to its head (dropping internal spine
arcs; an arc from tail back to head becomes a loop, harmless), and repeats until no chain
remains; it returns the list of (r, length) contracted. `compact()` is the r = 1 case and
is run first. The connector digraph keeps loops (they cancel in D_out - A).

Trap recorded. The clean filter in (A) rejects most random instances at small k because
random spacers produce chance repeats; k >= 9 and L >= 10 give a usable yield. A run with
k down to 6 produced only 8 clean instances out of 400 and tripped the minimum-count
threshold, which is why the parameters are what they are.
