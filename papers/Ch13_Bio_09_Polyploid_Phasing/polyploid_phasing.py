# -*- coding: utf-8 -*-
"""
polyploid_phasing.py -- polyploid phasing as a p-fold decomposition of a balanced flow
                                                                    (Bio 9, exploratory)

THE MODEL.  p circular haplotypes H_1..H_p are sequenced together; what the k-mer spectrum
retains is the POOLED de Bruijn multigraph G_k(H_1 + ... + H_p), with arc multiplicities
m = spec(H_1) + ... + spec(H_p).  A (labelled) phasing is a p-tuple of circular sequences
(T_1..T_p) with spec(T_1) + ... + spec(T_p) = m; an unlabelled phasing is the multiset.  So

    #labelled phasings  =  sum over ordered decompositions m = c_1 + ... + c_p  of
                           prod_i N_seq(c_i),
    c_i >= 0 balanced with connected support (a point of the circulation lattice Z_1 in the
    box [0, m]),  N_seq(c) = t(c) prod(d^+ - 1)! / prod c(a)!  as in Bio 8.

This is the Bio 8 lattice sum with the involution rho removed and the single constraint
c + rho c = m replaced by a p-fold partition of m.  The symmetric group S_p permutes the
labels of a decomposition; it acts freely on the decompositions with pairwise distinct c_i,
and #unlabelled = #labelled / p! exactly there.  Where haplotypes coincide (homozygous
tracts, identical spectra) the count is a multiset count, prod_c binom(N(c) + r_c - 1, r_c).

WHAT THE BRIEF ASKED THAT IS NOT SO.  (i) S_p is not an automorphism of the pooled graph:
the pooled graph has no symmetry exchanging haplotypes (the haplotypes are different words),
and there is no p-sheeted covering whose deck group is S_p -- the natural map from the
disjoint union of the haplotype graphs to the pooled graph is not a covering (a shared k-mer
has p preimages of local degree 1 over a vertex of degree p).  S_p acts on the DECOMPOSITION
SET, not on the graph.  (ii) The phasing space is not H^1(G_k, S_p) = Hom(pi_1, S_p): at a
biallelic site carried by a haplotypes out of p, the local choice is an a-subset of the
labels, S_p/(S_a x S_{p-a}), a homogeneous space and not the group; only for p = 2 (or a site
with p distinct alleles) is it S_p itself.  For unlinked sites the labelled count is
prod_sites binom(p, a_site) (times the single-haplotype reconstruction counts), which the
script verifies, and which |Hom(pi_1, S_p)| = (p!)^b does not reproduce.

CHECKS (phiX174 haplotypes with synthetic SNPs, k = 12 unless stated)
 (P)  the pooled compacted graph: balanced, contents with multiplicity, 2 branch vertices
      per SNP away from repeats                                                         exact
 (W)  p = 2: brute-force walk enumeration of T_1 reproduces N_seq(c_1) on every spectrum
      and the labelled count sum_{T_1} N_seq(m - spec T_1) equals the lattice sum     exact
 (U)  p = 2: unordered phasings enumerated as multisets of walks (T_1, T_2) equal the
      multiset formula                                                                  exact
 (L)  linkage: b unlinked SNPs give 2^b equal-length decompositions, two SNPs within k of
      each other are one site                                                           exact
 (T)  p = 3: labelled decompositions = prod_sites binom(3, a_site) for unlinked sites,
      6 at a triallelic site; unlabelled = multiset count, not labelled/6               exact
 (C)  dropping the equal-length constraint admits copy-number decompositions through the
      repeat (one haplotype with 3 copies, the other with 1)                            exact

Run: python -u polyploid_phasing.py           (seconds)
"""

import sys, os, time, math, collections, itertools
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)
# _SIB_HELPER

def _sib(name):
    """the folder `name`, wherever it sits: beside this script, beside any ancestor of it, or
    one level inside any ancestor.  The Bio notes live in an `Algebraic Genomics` subtree
    while `Paper 10` stays in the project root, so no fixed level finds both."""
    here = os.path.dirname(os.path.abspath(__file__))
    seen = []
    for _ in range(4):
        here = os.path.dirname(here)
        seen.append(here)
        cand = os.path.join(here, name)
        if os.path.isdir(cand):
            return cand
    for base in seen:
        try:
            for sub in sorted(os.listdir(base)):
                cand = os.path.join(base, sub, name)
                if os.path.isdir(cand):
                    return cand
        except OSError:
            pass
    return os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), name)

sys.path.insert(0, _sib("lib"))
sys.path.insert(0, _sib("Ch07_Bio_08_DS_Lattice_Sum"))
import pgl3_building as pb
from ds_lattice import load_phix, row, ROWS

BOX_CAP = int(os.environ.get("BOX_CAP", "2000000"))
WALK_CAP = int(os.environ.get("WALK_CAP", "30000000"))
NEXT = {"A": "C", "C": "G", "G": "T", "T": "A"}


def det(M):
    return 1 if not M else pb.bareiss_det(M)


def mutate(S, sites):
    """sites: {pos: shift} -- replace S[pos] by the base `shift` steps further along A C G T"""
    L = list(S)
    for p, s in sites.items():
        b = L[p]
        for _ in range(s):
            b = NEXT[b]
        L[p] = b
    return "".join(L)


# --------------------------------------------------------------------- the pooled graph

class Pool:
    """unitig-compacted de Bruijn multigraph of several circular sequences.  A k-mer is a
    branch vertex iff it has >= 2 distinct successors or >= 2 distinct predecessors in the
    pool; arcs are the maximal unitigs between branch vertices, with content string, length
    ell = #(k+1)-mers, multiplicity m, and the list of (haplotype, position) occurrences."""

    def __init__(self, seqs, k):
        self.seqs, self.k = seqs, k
        succ = collections.defaultdict(set)
        pred = collections.defaultdict(set)
        for T in seqs:
            N = len(T)
            TT = T + T[:k + 1]
            for i in range(N):
                w, x = TT[i:i + k], TT[i + 1:i + k + 1]
                succ[w].add(x)
                pred[x].add(w)
        branch = {w for w in succ if len(succ[w]) >= 2 or len(pred[w]) >= 2}
        self.words = sorted(branch)
        self.vid = {w: i for i, w in enumerate(self.words)}
        self.nv = len(self.words)
        self.arcs = []            # labelled: dict(u, v, hap, pos, content)
        for hi, T in enumerate(seqs):
            N = len(T)
            TT = T + T + T[:k]
            R = [i for i in range(N) if TT[i:i + k] in branch]
            assert R, f"haplotype {hi} has no branch vertex at k={k}"
            for t in range(len(R)):
                p, q = R[t], R[(t + 1) % len(R)]
                if q <= p:
                    q += N
                content = TT[p:q + k]
                self.arcs.append(dict(u=self.vid[content[:k]], v=self.vid[content[-k:]],
                                      hap=hi, pos=p, content=content))
        self.contents = sorted({a["content"] for a in self.arcs})
        self.cid = {c: i for i, c in enumerate(self.contents)}
        self.nA = len(self.contents)
        self.m = [0] * self.nA
        for a in self.arcs:
            self.m[self.cid[a["content"]]] += 1
        self.tail = [self.vid[c[:k]] for c in self.contents]
        self.head = [self.vid[c[-k:]] for c in self.contents]
        self.ell = [len(c) - k for c in self.contents]
        self.out = collections.defaultdict(list)
        for a in range(self.nA):
            self.out[self.tail[a]].append(a)
        self.spec = []
        for hi in range(len(seqs)):
            c = [0] * self.nA
            for a in self.arcs:
                if a["hap"] == hi:
                    c[self.cid[a["content"]]] += 1
            self.spec.append(tuple(c))
        self.N = [len(T) for T in seqs]

    def balanced(self):
        o, i = collections.Counter(), collections.Counter()
        for a in range(self.nA):
            o[self.tail[a]] += self.m[a]
            i[self.head[a]] += self.m[a]
        return all(o[v] == i[v] for v in range(self.nv))

    def length(self, c):
        return sum(c[a] * self.ell[a] for a in range(self.nA))

    # --- lattice ------------------------------------------------------------------------
    def rank_Z1(self):
        par = list(range(self.nv))

        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]
                x = par[x]
            return x
        for a in range(self.nA):
            x, y = find(self.tail[a]), find(self.head[a])
            if x != y:
                par[x] = y
        b0 = len({find(v) for v in range(self.nv)})
        return self.nA - self.nv + b0

    def enumerate_box(self, cap):
        """all balanced integer c with 0 <= c <= m (the circulation lattice in the box)"""
        m, tail, head = self.m, self.tail, self.head
        by_vertex = collections.defaultdict(list)
        for a in range(self.nA):
            by_vertex[tail[a]].append(a)
            by_vertex[head[a]].append(a)
        order, seen = [], set()
        for v in range(self.nv):
            for a in by_vertex[v]:
                if a not in seen:
                    seen.add(a)
                    order.append(a)
        remaining = [0] * self.nv
        for a in order:
            remaining[tail[a]] += m[a]
            remaining[head[a]] += m[a]
        bal = [0] * self.nv
        cur = [0] * self.nA
        out = []

        class Overflow(Exception):
            pass

        def rec(i):
            if i == len(order):
                out.append(tuple(cur))
                if len(out) > cap:
                    raise Overflow
                return
            a = order[i]
            u, v = tail[a], head[a]
            remaining[u] -= m[a]
            remaining[v] -= m[a]
            for g in range(0, m[a] + 1):
                bal[u] += g
                bal[v] -= g
                if abs(bal[u]) <= remaining[u] and abs(bal[v]) <= remaining[v]:
                    cur[a] = g
                    rec(i + 1)
                bal[u] -= g
                bal[v] += g
            remaining[u] += m[a]
            remaining[v] += m[a]

        try:
            rec(0)
        except Overflow:
            return None
        return out

    def nseq(self, c):
        """(N_seq(c), number of weak components of the support); 0 if disconnected or empty"""
        outd = collections.Counter()
        for a in range(self.nA):
            if c[a]:
                outd[self.tail[a]] += c[a]
        supp = sorted(outd)
        if not supp:
            return 0, 0
        par = {v: v for v in supp}

        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]
                x = par[x]
            return x
        for a in range(self.nA):
            if c[a]:
                x, y = find(self.tail[a]), find(self.head[a])
                if x != y:
                    par[x] = y
        ncomp = len({find(v) for v in supp})
        if ncomp != 1:
            return 0, ncomp
        pos = {v: i for i, v in enumerate(supp[1:])}
        n = len(supp) - 1
        L = [[0] * n for _ in range(n)]
        for a in range(self.nA):
            if c[a]:
                u, v = self.tail[a], self.head[a]
                if u in pos:
                    L[pos[u]][pos[u]] += c[a]
                    if v in pos:
                        L[pos[u]][pos[v]] -= c[a]
        t = det(L)
        loc = 1
        for v in supp:
            loc *= math.factorial(outd[v] - 1)
        denom = 1
        for a in range(self.nA):
            denom *= math.factorial(c[a])
        ec = t * loc
        # sum over circular sequences T with spectrum c of 1/d(T), d(T) = rotational symmetry
        # order; a non-integer value means every such T is a tandem repeat (a dimer circle)
        return Fraction(ec, denom), 1

    # --- walks ----------------------------------------------------------------------------
    def walks(self, bound, max_len, cap, collect=None):
        """rooted closed walks with usage <= bound and length <= max_len, weighted
        1/usage(a0) at the minimal content a0: returns {spectrum: weight sum}"""
        nA, tail, head, out, ell = self.nA, self.tail, self.head, self.out, self.ell
        usage = [0] * nA
        path = []
        buckets = collections.defaultdict(Fraction)
        nodes = [0]

        class Overflow(Exception):
            pass

        def rec(cur, length, a0, start):
            nodes[0] += 1
            if nodes[0] > cap:
                raise Overflow
            if cur == start:
                buckets[tuple(usage)] += Fraction(1, usage[a0])
                if collect is not None:
                    collect.append(tuple(path))
            for a in out[cur]:
                if a < a0 or usage[a] >= bound[a] or length + ell[a] > max_len:
                    continue
                usage[a] += 1
                path.append(a)
                rec(head[a], length + ell[a], a0, start)
                path.pop()
                usage[a] -= 1

        try:
            for a0 in range(nA):
                if bound[a0] < 1 or ell[a0] > max_len:
                    continue
                usage[a0] = 1
                path.append(a0)
                rec(head[a0], ell[a0], a0, tail[a0])
                path.pop()
                usage[a0] = 0
        except Overflow:
            return None, nodes[0]
        return buckets, nodes[0]


def canon(w):
    L = len(w)
    ww = w + w
    return min(ww[i:i + L] for i in range(L))


def multiset_count(nseq_of, spectra):
    """number of multisets of sequences realising the multiset of spectra"""
    tot = 1
    for c, r in collections.Counter(spectra).items():
        n = nseq_of[c]
        if n.denominator != 1:
            return None                       # a periodic part: no multiset count
        tot *= math.comb(int(n) + r - 1, r)
    return tot


# ------------------------------------------------------------------------------ analysis

def analyse(name, seqs, k, equal_length=True, nested_check=True, box_cap=BOX_CAP):
    p = len(seqs)
    G = Pool(seqs, k)
    print(f"\n=== {name}: p = {p}, k = {k} ===")
    print(f"    branch vertices {G.nv}, contents {G.nA}, sum m = {sum(G.m)}, rank Z_1 = {G.rank_Z1()}, "
          f"balanced {G.balanced()}")
    pts = G.enumerate_box(box_cap)
    if pts is None:
        print(f"    box has more than {box_cap} points -- skipped")
        return None
    N_of = {}
    for c in pts:
        n, _ = G.nseq(c)
        if n:
            N_of[c] = n
    print(f"    box points {len(pts)}, with connected support {len(N_of)}")
    for hi in range(p):
        assert G.spec[hi] in N_of, "a haplotype's own spectrum is not a connected box point"
    res = dict(G=G, pts=len(pts), conn=len(N_of), N_of=N_of)

    lengths_ok = (lambda c: G.length(c) == G.N[0]) if equal_length else (lambda c: True)
    m = tuple(G.m)

    # ordered decompositions m = c_1 + ... + c_p into connected points
    decs = []
    if p == 2:
        for c1, n1 in N_of.items():
            c2 = tuple(m[a] - c1[a] for a in range(G.nA))
            if c2 in N_of and lengths_ok(c1) and lengths_ok(c2):
                decs.append((c1, c2))
    elif p == 3:
        keys = [c for c in N_of if lengths_ok(c)]
        for c1 in keys:
            for c2 in keys:
                if any(c1[a] + c2[a] > m[a] for a in range(G.nA)):
                    continue
                c3 = tuple(m[a] - c1[a] - c2[a] for a in range(G.nA))
                if c3 in N_of and lengths_ok(c3):
                    decs.append((c1, c2, c3))
    labelled = sum(math.prod(N_of[c] for c in d) for d in decs)
    unl_sets = {tuple(sorted(d)) for d in decs}
    periodic = sum(1 for d in decs if any(N_of[c].denominator != 1 for c in d))
    mc = [multiset_count(N_of, d) for d in unl_sets]
    unlabelled = sum(x for x in mc if x is not None)
    distinct = sum(1 for d in decs if len(set(d)) == p)
    if labelled.denominator == 1:
        labelled = int(labelled)
    print(f"    {'equal-length ' if equal_length else ''}ordered decompositions {len(decs)} "
          f"({distinct} with pairwise distinct parts); labelled phasings {labelled}; "
          f"unordered decompositions {len(unl_sets)}; unlabelled phasings {unlabelled}"
          + (f"  ({periodic} ordered decompositions have a periodic part -- a tandem-dimer circle --"
             f" and are counted with weight 1/d(T))" if periodic else ""))
    res.update(decs=decs, labelled=labelled, unlabelled=unlabelled, unl_sets=unl_sets, periodic=periodic)

    # length profile of the parts (copy-number decompositions when lengths are free)
    if not equal_length:
        prof = collections.Counter(tuple(sorted(G.length(c) for c in d)) for d in decs)
        print(f"    length profiles of the parts: " +
              ", ".join(f"{list(l)}x{n}" for l, n in sorted(prof.items())[:8]))
        res["profiles"] = prof

    # (W) p = 2: walk enumeration of T_1, independent of the box
    if p == 2:
        maxlen = G.N[0] if equal_length else sum(G.m[a] * G.ell[a] for a in range(G.nA))
        walks = []
        buckets, nodes = G.walks(G.m, maxlen, WALK_CAP, walks)
        if buckets is None:
            print(f"    walk enumeration exceeded {WALK_CAP} nodes")
        else:
            okN = all(N_of.get(c, 0) == v for c, v in buckets.items()) \
                and all(c in buckets for c in N_of if lengths_ok(c))
            lab2 = Fraction(0)
            for c1, v in buckets.items():
                c2 = tuple(m[a] - c1[a] for a in range(G.nA))
                if c2 in N_of and lengths_ok(c2):
                    lab2 += v * N_of[c2]
            print(f"    (W) {nodes} walk nodes, {len(buckets)} spectra; N_seq matches on every spectrum: {okN}; "
                  f"sum_T1 N_seq(m - spec T1) = {lab2} vs lattice {labelled}: {lab2 == labelled}")
            res["W"] = okN and lab2 == labelled
            # (U) unordered phasings as multisets of explicit walks
            if nested_check and labelled <= 20000:
                seqs1 = {}
                for w in walks:
                    seqs1[canon(w)] = None
                pairs = set()
                for w1 in seqs1:
                    c1 = [0] * G.nA
                    for a in w1:
                        c1[a] += 1
                    c2 = tuple(m[a] - c1[a] for a in range(G.nA))
                    if c2 not in N_of or not lengths_ok(c2):
                        continue
                    ws2 = []
                    b2, _ = G.walks(list(c2), G.length(c2), WALK_CAP, ws2)
                    for w2 in ws2:
                        cc = [0] * G.nA
                        for a in w2:
                            cc[a] += 1
                        if tuple(cc) == c2:
                            pairs.add(tuple(sorted((w1, canon(w2)))))
                if equal_length:
                    print(f"    (U) {len(seqs1)} distinct T_1, {len(pairs)} unordered pairs {{T_1, T_2}} "
                          f"vs multiset formula {unlabelled}: {len(pairs) == unlabelled}")
                    res["U"] = (len(pairs) == unlabelled)
                else:
                    per = sum(1 for d in unl_sets if multiset_count(N_of, d) is None)
                    print(f"    (U) {len(seqs1)} distinct T_1, {len(pairs)} unordered pairs {{T_1, T_2}} explicitly; "
                          f"the multiset formula gives {unlabelled} on the {len(unl_sets) - per} decompositions "
                          f"without a periodic part, and the remaining {len(pairs) - unlabelled} pairs each contain "
                          f"a tandem-dimer circle (weight 1/2 in the lattice sum)")
    return res


def main():
    t0 = time.time()
    S = load_phix()
    N = len(S)
    k = 12
    print("=== polyploid phasing as a p-fold decomposition of the pooled flow (phiX174) ===")
    print(f"    genome {N} bp; k = {k}; single-haplotype reconstruction count N_seq(S) = ?")
    G0 = Pool([S], k)
    nS, _ = G0.nseq(list(G0.m))
    nS = int(nS)
    print(f"    N_seq(S) = {nS}  (branch vertices {G0.nv}, contents {G0.nA})")

    # SNP sites, away from the repeat k-mers and > k apart
    rep_pos = set()
    for a in G0.arcs:
        rep_pos.update(range(a["pos"] - k, a["pos"] + k + 1))
    sites = [600, 1700, 2900, 3900, 4900]
    assert all(s not in rep_pos for s in sites), "a SNP site touches a repeat k-mer"

    okP = okL = okT = okW = okU = okC = True

    # ---------------------------------------------------------------- diploid, unlinked
    print("\n--- (L) diploid with b unlinked SNPs: 2^b equal-length decompositions ---")
    for b in (1, 2, 3, 5):
        H2 = mutate(S, {s: 1 for s in sites[:b]})
        r = analyse(f"diploid, {b} SNP(s) at {sites[:b]}", [S, H2], k)
        G = r["G"]
        okP &= (G.nv == G0.nv + 2 * b) and G.balanced()
        good = (len(r["decs"]) == 2 ** b and r["labelled"] == 2 ** b * nS * nS
                and r["unlabelled"] == 2 ** (b - 1) * nS * nS)
        okL &= good
        okW &= r.get("W", True)
        okU &= r.get("U", True)
        print(f"    prediction 2^{b} = {2**b} decompositions, labelled {2**b}*{nS}^2 = {2**b*nS*nS}, "
              f"unlabelled {2**(b-1)*nS*nS}: {'ok' if good else 'MISMATCH'}")

    # ---------------------------------------------------------------- diploid, linked pair
    print("\n--- (L) diploid: two SNPs within k of each other are one site ---")
    H2 = mutate(S, {600: 1, 606: 1, 1700: 1, 2900: 1})
    r = analyse("diploid, SNPs at 600, 606 (linked), 1700, 2900", [S, H2], k)
    good = (len(r["decs"]) == 2 ** 3 and r["labelled"] == 8 * nS * nS)
    okL &= good
    okW &= r.get("W", True)
    okU &= r.get("U", True)
    print(f"    4 SNPs but 3 linkage sites: prediction 2^3 = 8 decompositions: {'ok' if good else 'MISMATCH'}")
    # and the same pair with the second haplotype carrying only one of the two: still 2 alleles
    H2 = mutate(S, {600: 1, 606: 1})
    H3 = mutate(S, {600: 1})
    r = analyse("diploid, haplotypes differ at 600 and 606 vs 600 only (one linked site, 2 alleles)",
                [H3, H2], k)
    good = (len(r["decs"]) == 2)
    okL &= good
    print(f"    prediction 2 decompositions: {'ok' if good else 'MISMATCH'}")
    row("P", "pooled compacted graph is balanced and gains exactly two branch vertices per SNP away "
             "from repeats", "5 pools", "exact", okP)
    row("L", "b unlinked biallelic sites give 2^b equal-length decompositions and 2^b N_seq(S)^2 labelled "
             "phasings; two sites within k are one site", "6 pools", "exact", okL)
    row("W", "p=2: brute-force enumeration of T_1 reproduces N_seq on every spectrum and the labelled "
             "count", "walk enumeration", "exact", okW)
    row("U", "p=2: unordered phasings as multisets of explicit walk pairs equal the multiset formula",
        "nested enumeration", "exact", okU)

    # ---------------------------------------------------------------- triploid
    print("\n--- (T) triploid: the local fibre is S_3/(S_a x S_{3-a}), not S_3 ---")
    # sites: 600 minor in H2 only (a=1); 1700 in H2 and H3 (a=2); 2900 triallelic (3 alleles)
    H2 = mutate(S, {600: 1, 1700: 1, 2900: 1})
    H3 = mutate(S, {1700: 1, 2900: 2})
    r = analyse("triploid: site 600 (1 of 3), site 1700 (2 of 3), site 2900 triallelic", [S, H2, H3], k,
                nested_check=False)
    pred = 3 * 3 * 6
    good = (len(r["decs"]) == pred and r["labelled"] == pred * nS ** 3)
    okT &= good
    print(f"    prediction binom(3,1)*binom(3,2)*3! = {pred} ordered decompositions, labelled {pred}*{nS}^3 = "
          f"{pred*nS**3}: {'ok' if good else 'MISMATCH'};  |Hom(pi_1, S_3)| would be 6^3 = 216")
    print(f"    unlabelled {r['unlabelled']} vs labelled/3! = {r['labelled']/6:.2f}: "
          f"{'equal' if r['unlabelled']*6 == r['labelled'] else 'NOT equal (homozygous-type multisets)'}")
    H2 = mutate(S, {600: 1})
    H3 = S
    r = analyse("triploid: one site, haplotypes S, S', S (two identical)", [S, H2, H3], k, nested_check=False)
    good = (len(r["decs"]) == 3 and r["unlabelled"] == nS * math.comb(nS + 1, 2))
    okT &= good
    print(f"    prediction 3 ordered decompositions, unlabelled N_seq * binom(N_seq+1, 2) = "
          f"{nS*math.comb(nS+1,2)}: {'ok' if good else 'MISMATCH'}")
    row("T", "p=3: ordered decompositions = prod_sites binom(3, a_site), 3! at a triallelic site; the "
             "unlabelled count is the multiset count, not labelled/3!", "2 pools", "exact", okT)

    # ---------------------------------------------------------------- copy number
    print("\n--- (C) without the equal-length constraint: copy-number decompositions ---")
    H2 = mutate(S, {600: 1})
    r = analyse("diploid, 1 SNP, all lengths", [S, H2], k, equal_length=False)
    prof = r["profiles"]
    cnv = sum(n for l, n in prof.items() if l[0] != l[1])
    okC &= cnv > 0 and any(l[0] == l[1] == N for l in prof)
    print(f"    {cnv} ordered decompositions with unequal haplotype lengths (repeat copies moved between "
          f"haplotypes), beside the {sum(n for l, n in prof.items() if l[0] == l[1] == N)} equal-length ones")
    r = analyse("homozygous diploid S + S, all lengths", [S, S], k, equal_length=False, nested_check=False)
    prof = r["profiles"]
    eq = [d for d in r["decs"] if all(G0.length(c) == N for c in d)] if False else \
        [d for d in r["decs"] if all(r["G"].length(c) == N for c in d)]
    unl_eq = sum(multiset_count(r["N_of"], d) for d in {tuple(sorted(d)) for d in eq})
    print(f"    homozygous pool: {len(r['decs'])} ordered decompositions, of which {len(eq)} at equal length "
          f"with {unl_eq} unlabelled phasings = binom({nS}+1, 2) = {math.comb(nS + 1, 2)}: "
          f"{unl_eq == math.comb(nS + 1, 2)}")
    okC &= (unl_eq == math.comb(nS + 1, 2))
    row("C", "free lengths admit copy-number decompositions through the repeat: one haplotype carries "
             "both copies of the 13-mer repeat block, the other none", "2 pools", "exact", okC)

    print("\n=== battery summary ===")
    bad = [x for x in ROWS if not x[4]]
    print(f"  TOTAL {sum(1 for x in ROWS if x[4])}/{len(ROWS)} passed"
          + ("" if not bad else f"   FAILURES: {[(x[0], x[2]) for x in bad]}"))
    print(f"  elapsed {time.time()-t0:.0f}s")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
