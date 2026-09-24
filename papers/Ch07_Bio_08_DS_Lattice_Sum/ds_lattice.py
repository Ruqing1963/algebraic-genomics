# -*- coding: utf-8 -*-
"""
ds_lattice.py -- double-stranded reconstructions as a lattice sum      (Bio 8, exploratory)

THE QUESTION ("Solved, open, impossible", Open 5).  K(D_k) counts the Eulerian circuits of the
reverse-complement double cover D_k, and that is the wrong object: an Eulerian circuit of D_k
walks both strands in one closed walk of 2N arcs, whereas a double-stranded reconstruction is
a circular sequence T of N arcs with

        spec(T) + spec(Tbar) = spec(S) + spec(Sbar)        (multisets of (k+1)-mers),

i.e. one closed walk C through HALF the arcs of D_k whose rho-image walks the other half.

THE LATTICE.  Write m for the arc multiplicities of D_k (a function on contents, i.e. on
distinct (k+1)-mers, or after unitig contraction on distinct segments) and c for the
spectrum of T.  Then

        c + rho c = m,   c >= 0,   c balanced,   supp c connected,

and with f := 2c - m this reads   rho f = -f,   f == m (mod 2),   |f| <= m,   f balanced.
So f runs over the box points of Lambda^-, the (-1)-eigenlattice of rho on the circulation
lattice Z_1(D_k), and

    #DS(S,k)  =  sum_{f in Lambda^- cap box}  N_seq((m+f)/2),
    N_seq(c)  =  t(c) * prod_v (d_c^+(v) - 1)! / prod_a c(a)!        (BEST, per sequence)

where t(c) is the arborescence count of the multigraph with multiplicities c.  N_seq(c) = 0
exactly when supp c is disconnected.  The division by prod c(a)! is NOT optional: the labelled
BEST count t(c) prod (d^+ - 1)! is the number of arc-rooted circuits of the multigraph with
labelled parallel copies, and a sequence T with spectrum c has prod_a c(a)! such labellings.
On one strand this factor is a global constant (there is one lattice point); across the
lattice it varies from point to point, so the sum must be taken per sequence.

RANK.  rho reverses arcs, so on the cellular chain complex of D_k it induces (-P_A, P_V), a
chain map; the Hopf trace formula gives  tr(rho | Z_1) = #fixA + #fixV - tr(P_V | H_0)
(#fixA = rho-fixed contents, #fixV = rho-fixed vertices = palindromic repeated k-mers), hence

        rank Lambda^-  =  (|A_c| - |V| + 2 - #fixA - #fixV) / 2

whether D_k is connected (b_0 = 1, rho trivial on H_0) or a swapped pair (b_0 = 2, trace 0).
This is NOT the number of inverted repeats: for a disconnected D_k = G_k + Gbar_k the lattice
is all of Z_1(G_k) and only two of its box points are connected.  The brief's guess (and my
own in the planning note) that rank Lambda^- counts independent inversions is refuted at
k = 12 below, where one inverted-repeat pair meets a lattice of rank 4.

CHECKS
 (D)  D_k built with labelled arcs: balanced, rho a free involution on labelled arcs,
      m(rho a) = m(a) on contents, rho-fixed contents have even multiplicity            exact
 (L)  rank Lambda^- equals (|A_c| - |V| + 2 - #fixA - #fixV)/2 at k = 9..13              exact
 (B)  every box point gives c >= 0, balanced, c + rho c = m                              exact
 (S)  the spectra of S and Sbar are box points and N_seq there reproduces the published
      single-strand count |K(G_k)| Loc(G_k) / prod m!  (|K(G_12)| = 2, |K(G_11)| = 132)   exact
 (E)  independent brute force: every closed walk of D_k with c + rho c = m, enumerated by
      backtracking without reference to the lattice, bucketed by spectrum; bucket set and
      every bucket count agree with the lattice sum                                       exact
 (N)  the un-normalised sum  sum_f t(c_f) prod(d^+ - 1)!  disagrees with the brute force
      as soon as a content has multiplicity >= 2 (k = 11); the per-sequence sum agrees    exact

Run: python -u ds_lattice.py        (phiX174 NC_001422.1; a few seconds for k = 12, 11)
Imports Paper 10/pgl3_building.py (Bareiss determinant) and Bio 5/ecoli_sandpile.py.
"""

import sys, os, time, math, collections
from fractions import Fraction
import numpy as np

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
sys.path.insert(0, _sib("Ch15_Bio_05_Ecoli_Decomposition"))
import pgl3_building as pb
import ecoli_sandpile as es

COMP = str.maketrans("ACGT", "TGCA")
ROWS = []
BOX_CAP = int(os.environ.get("BOX_CAP", "3000000"))
BRUTE_CAP = int(os.environ.get("BRUTE_CAP", "60000000"))
POINT_FULL_MAX = int(os.environ.get("POINT_FULL_MAX", "40000"))
BOX_CAP_K10 = int(os.environ.get("BOX_CAP_K10", "300000"))
# k = 10 has rank 48 and 125 vertices: the capacity-pruned backtracking below spends its time in
# dead ends and produces no points in 10 minutes, so the box is enumerated at k = 12, 11 only
# (set KS_ENUM=12,11,10 to try).  The rank table still covers k = 9..13.
KS_ENUM = tuple(int(x) for x in os.environ.get("KS_ENUM", "12,11").split(","))
sys.setrecursionlimit(10000)


def row(key, desc, cases, mode, ok):
    ROWS.append((key, desc, cases, mode, bool(ok)))
    print(("  [ok]   " if ok else "  [FAIL] ") + f"({key}) {desc}   [{cases}] [{mode}]")
    return bool(ok)


def rc(s):
    return s.translate(COMP)[::-1]


def det(M):
    return 1 if not M else pb.bareiss_det(M)


def load_phix():
    for folder in (HERE, _sib("Ch15_Bio_05_Ecoli_Decomposition"), _sib("Ch01_Bio_01_Critical_Groups")):
        p = os.path.join(folder, "phix174_NC_001422.1.txt")
        if os.path.exists(p):
            return open(p).read().strip().upper()
    raise SystemExit("phiX174 sequence not found")


# ------------------------------------------------------------- the double cover, labelled

class Cover:
    """D_k(S) after unitig contraction, with every arc labelled by its occurrence
    (strand, position) and by its CONTENT (the segment string from the start of the tail
    k-mer to the end of the head k-mer).  Contents are the arcs of the content-level
    multigraph; m is their multiplicity."""

    def __init__(self, S, k):
        self.S, self.k, self.N = S, k, len(S)
        N = self.N
        seqs = [S, rc(S)]
        occ = collections.defaultdict(list)
        for si, T in enumerate(seqs):
            TT = T + T[:k]
            for i in range(N):
                occ[TT[i:i + k]].append((si, i))
        rep = {w: ps for w, ps in occ.items() if len(ps) >= 2}
        self.occ = rep
        self.words = sorted(rep)
        self.vid = {w: i for i, w in enumerate(self.words)}
        self.nv = len(self.words)
        self.mult = [len(rep[w]) for w in self.words]
        self.vfix = [self.vid[rc(w)] for w in self.words]          # rho on vertices
        self.n_fixV = sum(1 for v in range(self.nv) if self.vfix[v] == v)
        # labelled arcs
        self.arcs = []
        key = {}
        for si, T in enumerate(seqs):
            TT = T + T + T[:k]
            R = sorted(p for w in rep for (s, p) in rep[w] if s == si)
            for t in range(len(R)):
                p = R[t]
                q = R[(t + 1) % len(R)]
                if q <= p:
                    q += N
                ln = q - p + k
                content = TT[p:p + ln]
                a = dict(u=self.vid[content[:k]], v=self.vid[content[-k:]], strand=si, pos=p,
                         len=ln, content=content)
                key[(si, p)] = len(self.arcs)
                self.arcs.append(a)
        # rho on labelled arcs: (s, p, len) -> (1-s, N-p-len mod N, len)
        self.arho = []
        for a in self.arcs:
            j = key.get((1 - a["strand"], (N - a["pos"] - a["len"]) % N))
            self.arho.append(j)
        # contents
        self.contents = sorted({a["content"] for a in self.arcs})
        self.cid = {c: i for i, c in enumerate(self.contents)}
        self.nA = len(self.contents)
        self.m = [0] * self.nA
        for a in self.arcs:
            self.m[self.cid[a["content"]]] += 1
        self.tail = [self.vid[c[:k]] for c in self.contents]
        self.head = [self.vid[c[-k:]] for c in self.contents]
        self.crho = [self.cid.get(rc(c)) for c in self.contents]
        self.n_fixA = sum(1 for a in range(self.nA) if self.crho[a] == a)
        self.out = collections.defaultdict(list)
        for a in range(self.nA):
            self.out[self.tail[a]].append(a)
        # the spectra of S and Sbar as content vectors
        self.c_strand = [[0] * self.nA, [0] * self.nA]
        for a in self.arcs:
            self.c_strand[a["strand"]][self.cid[a["content"]]] += 1
        # orbit representatives (non-fixed contents only; fixed ones carry f = 0)
        self.reps = [a for a in range(self.nA) if self.crho[a] is not None and a < self.crho[a]]

    # --- sanity ------------------------------------------------------------------------
    def sanity(self):
        outm = collections.Counter()
        inm = collections.Counter()
        for a in range(self.nA):
            outm[self.tail[a]] += self.m[a]
            inm[self.head[a]] += self.m[a]
        balanced = all(outm[v] == inm[v] == self.mult[v] for v in range(self.nv))
        rho_ok = all(j is not None and j != i and self.arho[j] == i
                     and self.arcs[j]["content"] == rc(self.arcs[i]["content"])
                     for i, j in enumerate(self.arho))
        m_ok = all(self.crho[a] is not None and self.m[self.crho[a]] == self.m[a]
                   for a in range(self.nA))
        fix_even = all(self.m[a] % 2 == 0 for a in range(self.nA) if self.crho[a] == a)
        return balanced, rho_ok, m_ok, fix_even

    def components(self):
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
        return len({find(v) for v in range(self.nv)})

    # --- the lattice ---------------------------------------------------------------------
    def balance_matrix(self):
        """rows = vertices, columns = orbit representatives j: the coefficient of g_j in the
        net outflow of f at v, where f(rep) = g, f(rho rep) = -g"""
        B = [[0] * len(self.reps) for _ in range(self.nv)]
        for j, a in enumerate(self.reps):
            b = self.crho[a]
            B[self.tail[a]][j] += 1
            B[self.head[a]][j] -= 1
            B[self.tail[b]][j] -= 1
            B[self.head[b]][j] += 1
        return B

    def rank_lattice(self):
        B = self.balance_matrix()
        if not self.reps:
            return 0
        r = rank_mod(B, 2147483647)
        r2 = rank_mod(B, 2147483629)
        assert r == r2, "modular ranks disagree"
        return len(self.reps) - r

    def f_of_g(self, g):
        f = [0] * self.nA
        for j, a in enumerate(self.reps):
            f[a] = g[j]
            f[self.crho[a]] = -g[j]
        return f

    def enumerate_box(self, cap):
        """all f in Lambda^- with |f| <= m and f == m (mod 2), by backtracking over orbit
        representatives with the balance constraint pruned by the remaining capacity"""
        reps, m, tail, head, crho = self.reps, self.m, self.tail, self.head, self.crho
        by_vertex = collections.defaultdict(list)
        for j, a in enumerate(reps):
            b = crho[a]
            for v in (tail[a], head[a], tail[b], head[b]):
                by_vertex[v].append(j)
        order, seen = [], set()
        for v in range(self.nv):
            for j in by_vertex[v]:
                if j not in seen:
                    seen.add(j)
                    order.append(j)
        ends = {}
        remaining = [0] * self.nv
        for j in order:
            a = reps[j]
            b = crho[a]
            ends[j] = [(tail[a], 1), (head[a], -1), (tail[b], -1), (head[b], 1)]
            for v, _ in ends[j]:
                remaining[v] += m[a]
        bal = [0] * self.nv
        cur = [0] * len(reps)
        out = []

        class Overflow(Exception):
            pass

        def rec(i):
            if i == len(order):
                out.append(tuple(cur))
                if len(out) > cap:
                    raise Overflow
                return
            j = order[i]
            a = reps[j]
            e = ends[j]
            for v, _ in e:
                remaining[v] -= m[a]
            for g in range(-m[a], m[a] + 1, 2):
                for v, s in e:
                    bal[v] += s * g
                if all(abs(bal[v]) <= remaining[v] for v, _ in e):
                    cur[j] = g
                    rec(i + 1)
                for v, s in e:
                    bal[v] -= s * g
            for v, _ in e:
                remaining[v] += m[a]

        try:
            rec(0)
        except Overflow:
            return None
        return out

    # --- per point -----------------------------------------------------------------------
    def point(self, f, connectivity_only=False):
        """c = (m+f)/2 and its sequence count"""
        c = []
        for a in range(self.nA):
            assert (self.m[a] + f[a]) % 2 == 0 and abs(f[a]) <= self.m[a]
            c.append((self.m[a] + f[a]) // 2)
        outd = collections.Counter()
        ind = collections.Counter()
        for a in range(self.nA):
            if c[a]:
                outd[self.tail[a]] += c[a]
                ind[self.head[a]] += c[a]
        assert all(outd[v] == ind[v] for v in set(outd) | set(ind)), "c not balanced"
        assert all(c[a] + c[self.crho[a]] == self.m[a] for a in range(self.nA)), "c + rho c != m"
        supp = sorted(outd)
        # connectivity of the support
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
        if ncomp != 1 or connectivity_only:
            return c, 0, 0, 0, ncomp
        # t(c): reduced weighted out-Laplacian on the support
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
        assert t >= 0
        loc = 1
        for v in supp:
            loc *= math.factorial(outd[v] - 1)
        denom = 1
        for a in range(self.nA):
            denom *= math.factorial(c[a])
        ec = t * loc
        assert ec % denom == 0, "labelled count not divisible by prod c(a)!"
        return c, ec // denom, ec, t, ncomp

    def cross_strand(self):
        """the k-mers occurring on both strands: inverted-repeat pairs {w, rc w} (one line per
        pair) and palindromes w = rc w, with their positions in S (1-based)"""
        pairs, pals = [], []
        for w in self.words:
            ps = self.occ[w]
            if not (any(s == 0 for s, _ in ps) and any(s == 1 for s, _ in ps)):
                continue
            fwd = sorted(p + 1 for s, p in ps if s == 0)
            if rc(w) == w:
                pals.append((w, fwd))
            elif w < rc(w):
                rev = sorted(p + 1 for s, p in self.occ[rc(w)] if s == 0)
                pairs.append((w, fwd, rev))
        return pairs, pals

    def describe(self, c):
        """which strand-0 segments of S are dropped and which strand-1 segments taken, as
        genomic intervals of S (1-based)"""
        cS = self.c_strand[0]
        lost = [a for a in range(self.nA) if c[a] < cS[a]]
        gained = [a for a in range(self.nA) if c[a] > cS[a]]
        N = self.N

        def spans(ids, strand):
            ivs = []
            for a in ids:
                for arc in self.arcs:
                    if arc["strand"] == strand and self.cid[arc["content"]] == a:
                        p, ln = arc["pos"], arc["len"]
                        if strand == 1:
                            p = (N - p - ln) % N
                        ivs.append((p + 1, p + ln))
                        break
            return sorted(ivs)
        return spans(lost, 0), spans(gained, 1)

    # --- brute force ----------------------------------------------------------------------
    def brute(self, cap, walks=None):
        """every closed walk C of the content multigraph with usage(a) + usage(rho a) = m(a),
        i.e. every double-stranded reconstruction, by backtracking.  A walk is generated once
        per occurrence of its smallest content a0, so each is weighted 1/usage(a0); the sum is
        an integer bucket by bucket when the sequences are aperiodic (asserted).  With a list
        `walks`, every rooted walk is appended to it as a tuple of content ids."""
        nA, m, tail, head, crho, out = self.nA, self.m, self.tail, self.head, self.crho, self.out
        half = sum(m) // 2
        usage = [0] * nA
        path = []
        buckets = collections.defaultdict(Fraction)
        nodes = [0]

        class Overflow(Exception):
            pass

        def rec(cur, steps, a0, start):
            nodes[0] += 1
            if nodes[0] > cap:
                raise Overflow
            if steps == half:
                if cur == start:
                    buckets[tuple(usage)] += Fraction(1, usage[a0])
                    if walks is not None:
                        walks.append(tuple(path))
                return
            for a in out[cur]:
                if a < a0:
                    continue
                b = crho[a]
                if b == a:
                    if 2 * usage[a] >= m[a]:
                        continue
                elif usage[a] + usage[b] >= m[a]:
                    continue
                usage[a] += 1
                path.append(a)
                rec(head[a], steps + 1, a0, start)
                path.pop()
                usage[a] -= 1

        try:
            for a0 in range(nA):
                b = crho[a0]
                if (b == a0 and m[a0] < 2) or m[a0] < 1:
                    continue
                usage[a0] = 1
                path.append(a0)
                rec(head[a0], 1, a0, tail[a0])
                path.pop()
                usage[a0] = 0
        except Overflow:
            return None, nodes[0]
        return buckets, nodes[0]


def rank_mod(B, p):
    A = np.array(B, dtype=np.int64) % p
    nr, nc = A.shape
    r = 0
    for c in range(nc):
        if r == nr:
            break
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0:
            continue
        piv = r + int(nz[0])
        if piv != r:
            A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        fac = A[:, c].copy()
        fac[r] = 0
        A = (A - np.outer(fac, A[r])) % p
        r += 1
    return r


def single_strand(S, k):
    """|K(G_k)|, Loc(G_k), prod m! on the single-strand compacted graph (Bio 5 code path)"""
    G = es.compacted_graph(S, k)
    if G["nv"] == 0:
        return 1, 1, 1
    K = det([list(r) for r in zip(*es.reduced_laplacian(G["nv"], G["arcs"]))])
    loc = 1
    for ps in G["occ"].values():
        loc *= math.factorial(len(ps) - 1)
    # multiplicities of the (k+1)-mers of S (= of the compacted arcs by content)
    N = len(S)
    SS = S + S[:k + 1]
    cnt = collections.Counter(SS[i:i + k + 1] for i in range(N))
    pm = 1
    for x in cnt.values():
        pm *= math.factorial(x)
    return abs(K), loc, pm


def fmt(x):
    if x is None:
        return "-"
    if isinstance(x, int) and x >= 10 ** 12:
        return f"{x:.3e}"
    return str(x)


# ----------------------------------------------------------------------------------- main

def main():
    t0 = time.time()
    S = load_phix()
    print("=== double-stranded reconstructions of phiX174 as a lattice sum ===")
    print(f"    genome {len(S)} bp circular")

    # ------------------------------------------------------------------ (D), (L)
    print("\n=== (D) the double cover with labelled arcs; (L) rank of Lambda^- ===")
    print("      k   |V|  fixV   |A_c|  fixA   sum m   b0   rank Z_1   rank Lambda^-   (|A_c|-|V|+2-fixA-fixV)/2")
    covers = {}
    okD = okL = True
    for k in (13, 12, 11, 10, 9):
        D = Cover(S, k)
        covers[k] = D
        if D.nv == 0:
            print(f"     {k:2d}     0     0       0     0       0    -          -               0                          -"
                  "     (no repeated k-mer on either strand: D_k is the two disjoint cycles S, Sbar)")
            continue
        bal, rho_ok, m_ok, fix_even = D.sanity()
        okD &= bal and rho_ok and m_ok and fix_even
        b0 = D.components()
        rZ = D.nA - D.nv + b0
        rL = D.rank_lattice()
        pred2 = D.nA - D.nv + 2 - D.n_fixA - D.n_fixV
        good = (pred2 % 2 == 0) and (rL == pred2 // 2)
        okL &= good
        print(f"     {k:2d}  {D.nv:4d}  {D.n_fixV:4d}   {D.nA:5d}  {D.n_fixA:4d}   {sum(D.m):5d}   {b0:2d}   "
              f"{rZ:8d}   {rL:13d}   {pred2/2:26.1f}     {'ok' if good else 'MISMATCH'}"
              + ("" if (bal and rho_ok and m_ok and fix_even) else
                 f"   sanity: balanced {bal} rho {rho_ok} m {m_ok} fixed-even {fix_even}"))
    row("D", "D_k with labelled arcs is balanced; rho is a free involution on labelled arcs with "
             "rc contents; m(rho a) = m(a); rho-fixed contents have even multiplicity",
        "k = 9..12", "exact", okD)
    row("L", "rank Lambda^- = (|A_c| - |V| + 2 - #fixA - #fixV)/2, the Hopf trace of the "
             "arc-reversing involution -- not the number of inverted repeats",
        "k = 9..12", "exact", okL)
    print("\n    the k-mers shared by the two strands (what glues D_k):")
    for k in (12, 11):
        pairs, pals = covers[k].cross_strand()
        print(f"      k={k}: {len(pairs)} inverted-repeat pair(s), {len(pals)} palindrome(s); "
              f"rank Lambda^- = {covers[k].rank_lattice()}")
        for w, fwd, rev in pairs:
            print(f"         {w} at {fwd}   /   {rc(w)} at {rev}")
        for w, fwd in pals:
            print(f"         {w} = its own reverse complement, at {fwd}")

    # ------------------------------------------------------------ the lattice sum
    results = {}
    for k in KS_ENUM:
        D = covers[k]
        print(f"\n=== k = {k}: the box  Lambda^- cap {{|f| <= m, f == m mod 2}} ===")
        t1 = time.time()
        cap = BOX_CAP if k >= 11 else BOX_CAP_K10
        pts = D.enumerate_box(cap)
        if pts is None:
            print(f"    more than {cap} box points (rank {D.rank_lattice()}) -- enumeration abandoned"
                  f"   [{time.time()-t1:.0f}s]")
            results[k] = None
            continue
        print(f"    {len(pts)} box points   [{time.time()-t1:.1f}s]")
        if len(pts) > POINT_FULL_MAX:
            nconn = 0
            t1 = time.time()
            for g in pts:
                f = D.f_of_g(g)
                if D.point(f, connectivity_only=True)[4] == 1:
                    nconn += 1
            print(f"    connected support: {nconn} of {len(pts)} points   [{time.time()-t1:.0f}s]")
            print(f"    (more than {POINT_FULL_MAX} points: the determinant sum is not attempted here)")
            results[k] = dict(pts=len(pts), nconn=nconn, sum_seq=None, sum_unnorm=None, table=None)
            continue
        okB = True
        table = {}
        nconn = 0
        sum_seq = 0
        sum_unnorm = 0
        t1 = time.time()
        for g in pts:
            f = D.f_of_g(g)
            try:
                c, nseq, ec, t, ncomp = D.point(f)
            except AssertionError as e:
                okB = False
                print(f"    point {g}: {e}")
                continue
            table[tuple(c)] = (nseq, ec, t, ncomp, g)
            if ncomp == 1:
                nconn += 1
            sum_seq += nseq
            sum_unnorm += ec
        print(f"    every point: c >= 0, balanced, c + rho c = m   [{time.time()-t1:.1f}s]")
        print(f"    connected support (N_seq > 0): {nconn} of {len(pts)} points")
        print(f"    #DS = sum_f N_seq(c_f) = {sum_seq}      un-normalised sum_f t(c_f) prod(d+-1)! = {sum_unnorm}")
        row("B", f"k={k}: all {len(pts)} box points give c >= 0, balanced, c + rho c = m", f"{len(pts)} points",
            "exact", okB)

        # (S) the two single-strand points
        K1, loc1, pm1 = single_strand(S, k)
        okS = True
        for si, name in ((0, "S"), (1, "Sbar")):
            cS = tuple(D.c_strand[si])
            hit = table.get(cS)
            if hit is None:
                okS = False
                print(f"    spectrum of {name} is NOT a box point")
                continue
            nseq, ec, t, ncomp, g = hit
            good = (ec == K1 * loc1) and (nseq * pm1 == ec)
            okS &= good
            print(f"    spectrum of {name:4s}: t = {t} = |K(G_{k})| = {K1};  labelled count {ec} = |K| Loc = {K1*loc1};"
                  f"  N_seq = {nseq} = {ec}/{pm1}   {'ok' if good else 'MISMATCH'}")
        row("S", f"k={k}: the spectra of S and Sbar are box points and the lattice formula there returns the "
                 f"published single-strand count |K(G_{k})| Loc / prod m!", "2 points", "exact", okS)

        # the connected points, described
        conn = sorted(((v[0], c) for c, v in table.items() if v[3] == 1), key=lambda x: -x[0])
        print(f"    the connected points (spectrum classes), by N_seq:")
        shown = 0
        for nseq, c in conn:
            lost, gained = D.describe(list(c))
            if c == tuple(D.c_strand[0]):
                what = "= spec(S)"
            elif c == tuple(D.c_strand[1]):
                what = "= spec(Sbar)"
            elif k == 12:
                what = ("drops S" + "".join(f"[{a}..{b}]" for a, b in lost)
                        + "  takes rc(S" + "".join(f"[{a}..{b}]" for a, b in gained) + ")")
            else:
                what = (f"drops {len(lost)} segment(s) of S in {lost[0][0]}..{max(b for _, b in lost)}"
                        f", takes {len(gained)} of Sbar" if lost else f"takes {len(gained)} segment(s) of Sbar")
            print(f"      N_seq = {nseq:8d}   {what}")
            shown += 1
            if shown >= 14 and len(conn) > 16:
                print(f"      ... {len(conn) - shown} more")
                break
        hist = collections.Counter(v[0] for v in table.values() if v[3] == 1)
        print(f"    N_seq histogram over the connected points: "
              + ", ".join(f"{n}x{cnt}" for n, cnt in sorted(hist.items(), reverse=True)))
        results[k] = dict(pts=len(pts), nconn=nconn, sum_seq=sum_seq, sum_unnorm=sum_unnorm, table=table)

    # ------------------------------------------------------------------ (E) brute force
    for k in KS_ENUM:
        D = covers[k]
        if results.get(k) is None or results[k]["table"] is None:
            continue
        print(f"\n=== (E) k = {k}: brute force over closed walks with usage + rho usage = m ===")
        t1 = time.time()
        buckets, nodes = D.brute(BRUTE_CAP)
        if buckets is None:
            print(f"    more than {BRUTE_CAP} search nodes -- abandoned   [{time.time()-t1:.0f}s]")
            if k != 10:
                row("E", f"k={k}: brute force abandoned at {BRUTE_CAP} nodes", "-", "exact", False)
            results[k]["brute"] = None
            continue
        table = results[k]["table"]
        integral = all(v.denominator == 1 for v in buckets.values())
        total = sum(buckets.values())
        found = {c: int(v) for c, v in buckets.items()}
        predicted = {c: v[0] for c, v in table.items() if v[0] > 0}
        same_set = set(found) == set(predicted)
        same_counts = same_set and all(found[c] == predicted[c] for c in found)
        print(f"    {nodes} search nodes, {len(found)} spectra found, {total} reconstructions   [{time.time()-t1:.0f}s]")
        print(f"    all bucket weights integral: {integral}")
        print(f"    spectra found == lattice points with N_seq > 0: {same_set}"
              f"   ({len(found)} vs {len(predicted)})")
        print(f"    every bucket count == N_seq(c_f): {same_counts}")
        if not same_counts and same_set:
            bad = [(found[c], predicted[c]) for c in found if found[c] != predicted[c]][:5]
            print(f"    first mismatches (brute, lattice): {bad}")
        print(f"    brute-force total {total}  vs  lattice sum {results[k]['sum_seq']}"
              f"  vs  un-normalised sum {results[k]['sum_unnorm']}")
        okE = integral and same_counts and total == results[k]["sum_seq"]
        row("E", f"k={k}: the brute-force enumeration of double-stranded reconstructions returns exactly the "
                 f"lattice points with connected support, with the same count at every point",
            f"{len(found)} spectra, {total} reconstructions", "exact", okE)
        results[k]["brute"] = int(total)
        if k == 11:
            has_bundle = any(x >= 2 for x in D.m)
            okN = has_bundle and (results[k]["sum_unnorm"] != results[k]["sum_seq"])
            row("N", "k=11 has contents of multiplicity >= 2, and there the un-normalised BEST sum over the "
                     "lattice is NOT the reconstruction count; the per-sequence sum is",
                "1 value of k", "exact", okN)

    # ---------------------------------------------------------------- the comparison
    print("\n=== the comparison: what the blind double-cover group counts, and what is true ===")
    print("      k   |K(D_k)|     ec(D_k) = |K| Loc(D_k)    single-strand N_seq(S)+N_seq(Sbar)   "
          "box pts   connected   #DS (lattice)   #DS (brute)")
    for k in (12, 11, 10):
        D = covers[k]
        if k not in KS_ENUM and k != 10:
            continue
        M = es.reduced_laplacian(D.nv, [(D.tail[a], D.head[a]) for a in range(D.nA) for _ in range(D.m[a])])
        KD = abs(det([list(r) for r in zip(*M)]))
        locD = 1
        for v in range(D.nv):
            locD *= math.factorial(D.mult[v] - 1)
        K1, loc1, pm1 = single_strand(S, k)
        ss = 2 * (K1 * loc1 // pm1)
        r = results.get(k)
        print(f"     {k:2d}   {fmt(KD):>9s}   {fmt(KD*locD):>22s}   {fmt(ss):>32s}   "
              f"{fmt(r['pts'] if r else None):>7s}   {fmt(r['nconn'] if r else None):>9s}   "
              f"{fmt(r['sum_seq'] if r else None):>13s}   {fmt(r.get('brute') if r else None):>11s}")
    print("    (ec(D_k) counts Eulerian circuits of the double cover -- one walk through both strands;")
    print("     #DS counts circular sequences T with spec(T) + spec(Tbar) = spec(S) + spec(Sbar).)")

    print("\n=== battery summary ===")
    bad = [x for x in ROWS if not x[4]]
    print(f"  TOTAL {sum(1 for x in ROWS if x[4])}/{len(ROWS)} passed"
          + ("" if not bad else f"   FAILURES: {[(x[0], x[2]) for x in bad]}"))
    print(f"  elapsed {time.time()-t0:.0f}s")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
