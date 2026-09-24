# Ancillary bundle for "Multi-copy repeats: where the two-copy theory of assembly ambiguity stops" (`multicopy.tex`, Bio 4)

`multicopy.py` reproduces every computational claim: 4 checks, 4 passing, about 40 s.
It imports `repeat_splitting.py` from the sibling folder `Bio 3` (which in turn imports
`pgl3_building.py` from `Paper 10` for exact integer Smith normal forms), so keep the
folders as siblings. Run `python multicopy.py`. The archived output is
`multicopy_output.txt`. `phix174_NC_001422.1.txt` (5386 bp) is read from beside the
script; no network access is needed.

**The question.** Bio 3 splits the sandpile group of a de Bruijn graph into a
length-driven part and the group K(D_w) of the transition digraph of the circular word of
repeat occurrences. For two-copy repeats K(D_w) is the critical group of a bouquet,
coker(A + I) for Bouchet's signed interlace matrix (Merino-Moffatt-Noble 2025): a
presentation whose entries are pairwise, "topological" data. This note asks what survives
when repeats have three or more copies. Answer: the two-family case is a one-line theorem
(K = Z/c), and nothing pairwise survives beyond it.

Checks, keyed as in the paper's battery table:

| key | content | arithmetic |
|-----|---------|-----------|
| (A) | two repeat families of multiplicities r, s <= 5, all 119 circular words: K(D_w) = Z/c, c = number of A->B transitions (= number of A-blocks), and every c in 1..min(r, s) occurs | exact |
| (B) | census over all circular words for multiplicity vectors (3,3), (4,4), (3,3,3), (3,3,2), (2,2,3), (4,4,4), (3,3,3,3): number of words, of distinct transition digraphs D_w, of distinct groups, largest order, the non-cyclic groups; and the number of classes of pairwise restriction necklaces that carry more than one group. Witness: AABACCBBC (Z/3) and AABABCBCC (Z/4) have identical pairwise necklaces | exact |
| (C) | detaching a 3-in 3-out vertex into two 2-in 2-out vertices (nine ways): on 150 random balanced digraphs K(G) is a quotient of K(G') in only 140/1350 cases; witness K(G) = Z/4 + Z/5 whose detachments give Z/23, Z/25, Z/31, Z/32, Z/2+Z/13, ... | exact |
| (D) | phiX174 skeleton out-degree profile at k = 12..7: 2-in 2-out at k = 11 (8 vertices) and k = 10 (37), where coker(A + I) of a Hierholzer Euler circuit reproduces K(Sk) exactly; from k = 9 down (110, 312, 826 vertices) there are 6, 48, 207 vertices of out-degree >= 3, up to 6 | exact |

Conventions. Circular words are enumerated one per rotation class (first letter fixed,
then canonical rotation). "Pairwise restriction necklace" of (e, f) = the cyclic binary
word of the occurrences of e and f, up to rotation. Groups are compared by elementary
divisors. Detachment: the lone in-arc goes to the new vertex v'', the lone out-arc stays
at v', with one new arc v' -> v''.

Two dead ends recorded so they are not retried: (1) the linear (base-pointed) pairwise
restrictions determine the whole word for any m, so "K is a function of them" is vacuous;
only rotation classes are a real invariant, and those fail. (2) The Gram matrix
D^T (I - S) D of consecutive-occurrence differences is NOT Merino-Moffatt-Noble's A + I
even for r = 2 (AACBBC gives [[3,0,1],[0,3,-1],[1,-1,2]] against the identity), so that
route to an r >= 3 interlace matrix is closed.
