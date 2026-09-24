# -*- coding: utf-8 -*-
r"""
bio17_verify.py -- inversions at inverted repeats are moves inside the double-stranded fibre.

Setting (Bio 8, 12, 13).  S is a circular sequence, k a word length, D_k(S) the labelled
reverse-complement double cover: vertices are the k-mers repeated in S + rc(S), arcs ("contents")
the stretches between consecutive repeated k-mers, rho the reverse complement.  A double-stranded
reconstruction is a circular T with the same double-stranded (k+1)-mer spectrum,
spec(T) + spec(rc T) = spec(S) + spec(rc S); #DS(S,k) counts them (T and rc T separately, up to
rotation), and Bio 8 writes it as a sum over the points f of the lattice Lambda^- in a box.

Claim tested here.  Write T = P w Z' rc(w) R with w a k-mer.  Replacing the segment
Z = w Z' rc(w) by rc(Z) -- an INVERSION FLANKED BY AN INVERTED REPEAT -- keeps the double-stranded
spectrum, because rc(Z) begins with w and ends with rc(w) (the junction words are unchanged) and the
words inside Z are traded for their reverse complements.  Together with the classical moves of the
single-stranded fibre (transpositions at repeated k-mers, Ukkonen / Pevzner) and global reverse
complementation, do these moves connect the whole of DS(S,k)?

Checks (all exact):
  (I) the inversion theorem at the string level: every inversion flanked by an inverted repeat of
      length k preserves the double-stranded (k+1)-spectrum -- on random sequences and on phiX174 at
      k = 10, 11, 12 -- while unflanked inversions of random segments almost never do;
  (C) three independent counts of DS(S,k) agree on random small sequences: brute-force enumeration of
      circular strings, Bio 8's lattice sum of N_seq, and the closure of S under the moves;
  (G) generation: the closure of S under {transpositions, IR-inversions, reverse complement} is all of
      DS(S,k) -- random cases and phiX174 at k = 12 (#DS = 10) and k = 11 (#DS = 864); the closures
      under the sub-families show what each move contributes;
  (L) the lattice: the inversion vectors (1 - rho) z, z a path from a vertex v to its reverse
      complement vbar, span a sublattice M of Lambda^-; its rank and index [Lambda^- : M], and the
      index of the cycle part (1 - rho) Z_1, on phiX174 at k = 10, 11, 12 and random covers;
  (H) the Tate dichotomy: on 600 random connected covers, [Lambda^- : (1-rho)Z_1] = 2 exactly when
      D_k has no palindromic vertex and no self-complementary content, and M = Lambda^- always;
  (W) phiX174 at k = 10 (#DS = 31 925 753 246 212): the inverted repeats of the genome itself, the
      lattice points one inversion reaches, and a long random walk of moves, every state checked.

    python -u bio17_verify.py            (about 5-15 minutes)
    python -u bio17_verify.py --quick    skips the long random walk of (W)

Needs ds_fermion.py of Bio 13 (imported from ../Bio 13) and the standard library.
"""
import os, sys, time, random, collections, itertools

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "Ch09_Bio_13_Fermionic_Partition"))
import ds_fermion as DF                                   # noqa: E402
from ds_fermion import Cover, rc, enumerate_box           # noqa: E402

ROWS = []


def row(key, desc, ok):
    ROWS.append((key, bool(ok)))
    print(("  [ok]   " if ok else "  [FAIL] ") + f"({key}) {desc}", flush=True)


# ------------------------------------------------------------------ strings

def cyc_kmers(s, K):
    ss = s + s[:K - 1]
    return collections.Counter(ss[i:i + K] for i in range(len(s)))


def ds_spec(s, K):
    c = cyc_kmers(s, K)
    c.update(cyc_kmers(rc(s), K))
    return c


def canon_str(s):
    return min(s[i:] + s[:i] for i in range(len(s)))


def ir_inversions(s, k):
    """all inversions of s flanked by an inverted repeat w ... rc(w) of length k (cyclic)"""
    N = len(s)
    ss = s + s[:k - 1]
    pos = collections.defaultdict(list)
    for i in range(N):
        pos[ss[i:i + k]].append(i)
    out = []
    for w, P in pos.items():
        Q = pos.get(rc(w))
        if not Q:
            continue
        for i in P:
            for j in Q:
                d = (j - i) % N
                if d == 0 or d + k > N:
                    continue
                t = s[i:] + s[:i]
                out.append(rc(t[:d + k]) + t[d + k:])
    return out


def brute_ds(S, k, nodecap=3_000_000):
    """DS(S,k) by exhaustive search: all circular strings with the double-stranded (k+1)-spectrum
    of S, returned as canonical rotations.  None if the search exceeds nodecap."""
    K, N = k + 1, len(S)
    D = ds_spec(S, K)
    key = lambda x: min(x, rc(x))
    need = collections.Counter()
    for x, c in D.items():
        need[key(x)] = c
    used = collections.Counter()
    out, s, nodes = set(), [], [0]

    def add(x, sg):
        kk = key(x)
        used[kk] += sg * (2 if x == rc(x) else 1)
        return used[kk] <= need[kk]

    def rec():
        nodes[0] += 1
        if nodes[0] > nodecap:
            return
        if len(s) == N:
            st = "".join(s)
            added, ok = [], True
            for p in range(N - K + 1, N):
                x = st[p:] + st[:K - (N - p)]
                ok &= add(x, 1)
                added.append(x)
            if ok and all(used[q] == need[q] for q in need):
                out.add(canon_str(st))
            for x in added:
                add(x, -1)
            return
        for ch in "ACGT":
            s.append(ch)
            if len(s) >= K:
                x = "".join(s[-K:])
                if add(x, 1):
                    rec()
                add(x, -1)
            else:
                rec()
            s.pop()
    rec()
    return None if nodes[0] > nodecap else out


# ------------------------------------------------------------------ circuits of contents

class Tok:
    """a double-stranded reconstruction as a cyclic sequence of contents of D_k"""

    def __init__(self, cov):
        self.cov, self.k = cov, cov.k
        self.tail, self.head, self.rho = cov.tail, cov.head, cov.crho
        self.vbar = [cov.vid[rc(w)] for w in cov.words]

    def circuit_S(self):
        return tuple(self.cov.cid[a["content"]] for a in self.cov.arcs if a["strand"] == 0)

    def valid(self, t):
        n = len(t)
        if any(self.head[t[i]] != self.tail[t[(i + 1) % n]] for i in range(n)):
            return False
        cnt = collections.Counter(t)
        m = self.cov.m
        for a in range(self.cov.nA):
            b = self.rho[a]
            if (2 * cnt[a] if b == a else cnt[a] + cnt[b]) != m[a]:
                return False
        return True

    @staticmethod
    def canon(t):
        return min(t[i:] + t[:i] for i in range(len(t)))

    def spell(self, t):
        return "".join(self.cov.contents[a][:-self.k] for a in t)

    def point(self, t):
        cnt = collections.Counter(t)
        return tuple(cnt[a] - cnt[self.rho[a]] for a in self.cov.reps)

    def moves(self, t, kinds):
        """descriptors of the moves available at t"""
        n = len(t)
        occ = collections.defaultdict(list)
        for i, a in enumerate(t):
            occ[self.tail[a]].append(i)
        out = []
        if "T" in kinds:
            for v, P in occ.items():
                for a_, b_, c_ in itertools.combinations(P, 3):
                    out.append(("T3", a_, b_, c_))
            vs = [v for v in occ if len(occ[v]) >= 2]
            for u, w in itertools.combinations(vs, 2):
                for a_, c_ in itertools.combinations(occ[u], 2):
                    for b_, d_ in itertools.combinations(occ[w], 2):
                        if (a_ < b_ < c_) != (a_ < d_ < c_):
                            out.append(("T4", a_, c_, b_, d_))
        if "I" in kinds:
            for v, P in occ.items():
                for i in P:
                    for j in occ.get(self.vbar[v], ()):
                        if j != i:
                            out.append(("I", i, j))
        if "R" in kinds:
            out.append(("R",))
        return out

    def apply(self, t, d):
        n = len(t)
        if d[0] == "T3":
            _, a_, b_, c_ = d
            r = t[a_:] + t[:a_]
            B, C = b_ - a_, c_ - a_
            return r[B:C] + r[:B] + r[C:]
        if d[0] == "T4":
            _, a_, c_, b_, d_ = d
            r = t[a_:] + t[:a_]
            x, y = (b_ - a_) % n, (d_ - a_) % n
            B, D = (x, y) if x < c_ - a_ else (y, x)
            C = c_ - a_
            return r[C:D] + r[B:C] + r[:B] + r[D:]
        if d[0] == "I":
            _, i, j = d
            r = t[i:] + t[:i]
            dd = (j - i) % n
            return tuple(self.rho[a] for a in reversed(r[:dd])) + r[dd:]
        return tuple(self.rho[a] for a in reversed(t))

    def closure(self, start, kinds, cap=200000):
        seen = {self.canon(start)}
        frontier = [start]
        bad = 0
        while frontier:
            nxt = []
            for t in frontier:
                for d in self.moves(t, kinds):
                    u = self.apply(t, d)
                    cu = self.canon(u)
                    if cu in seen:
                        continue
                    if not self.valid(u):
                        bad += 1
                        continue
                    seen.add(cu)
                    nxt.append(cu)
                    if len(seen) > cap:
                        return None, bad
            frontier = nxt
        return seen, bad


def lattice_count(cov):
    return sum(cov.n_seq(cov.point(g)) for g in enumerate_box(cov))


# ------------------------------------------------------------------ lattices

def diag_invariants(rows):
    """rank and Smith invariants of an integer matrix (list of rows), by unimodular elimination"""
    A = [list(r) for r in rows if any(r)]
    if not A:
        return 0, []
    m, n = len(A), len(A[0])
    diag, t = [], 0
    while t < m and t < n:
        piv = None
        for i in range(t, m):
            for j in range(t, n):
                if A[i][j] and (piv is None or abs(A[i][j]) < abs(A[piv[0]][piv[1]])):
                    piv = (i, j)
        if piv is None:
            break
        i, j = piv
        A[t], A[i] = A[i], A[t]
        for r in A:
            r[t], r[j] = r[j], r[t]
        while True:
            p = A[t][t]
            done = True
            for i in range(t + 1, m):
                if A[i][t]:
                    q = A[i][t] // p
                    A[i] = [x - q * y for x, y in zip(A[i], A[t])]
                    if A[i][t]:
                        done = False
            for j in range(t + 1, n):
                if A[t][j]:
                    q = A[t][j] // p
                    for r in A:
                        r[j] -= q * r[t]
                    if A[t][j]:
                        done = False
            if done:
                break
            # move the smallest remaining entry of row/column t to the pivot
            best = (abs(A[t][t]), t, t)
            for i in range(t + 1, m):
                if A[i][t] and abs(A[i][t]) < best[0]:
                    best = (abs(A[i][t]), i, t)
            for j in range(t + 1, n):
                if A[t][j] and abs(A[t][j]) < best[0]:
                    best = (abs(A[t][j]), t, j)
            _, i, j = best
            if i != t:
                A[t], A[i] = A[i], A[t]
            if j != t:
                for r in A:
                    r[t], r[j] = r[j], r[t]
        diag.append(abs(A[t][t]))
        t += 1
    # normalise the diagonal to Smith form: (a, b) -> (gcd, lcm)
    d = diag[:]
    changed = True
    while changed:
        changed = False
        for i in range(len(d)):
            for j in range(i + 1, len(d)):
                g = math_gcd(d[i], d[j])
                l = d[i] // g * d[j]
                if (d[i], d[j]) != (g, l):
                    d[i], d[j] = g, l
                    changed = True
    return len(d), [x for x in d if x != 1]


def math_gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def inversion_lattice(cov):
    """generators of M_cyc = (1-rho) Z_1 and M = M_cyc + <(1-rho) p_v>, p_v a path v -> vbar,
    in the orbit coordinates g_j = z(rep_j) - z(rho rep_j) of Lambda^-"""
    nv, nA = cov.nv, cov.nA
    tail, head, rho, reps = cov.tail, cov.head, cov.crho, cov.reps
    adj = collections.defaultdict(list)
    for a in range(nA):
        adj[tail[a]].append((head[a], a, 1))
        adj[head[a]].append((tail[a], a, -1))
    par, order, tree = {0: None}, [0], set()
    for x in order:
        for y, a, sg in adj[x]:
            if y not in par:
                par[y] = (x, a, sg)
                tree.add(a)
                order.append(y)
    assert len(par) == nv, "D_k is not connected"

    def up(x):
        z = collections.Counter()
        while par[x] is not None:
            p, a, sg = par[x]
            z[a] -= sg                      # the arc from x up to p: sg=+1 means p -> x
            x = p
        return z

    UP = {x: up(x) for x in range(nv)}

    def path(u, v):
        z = collections.Counter(UP[u])
        z.subtract(UP[v])
        return z

    def g_of(z):
        return [z[a] - z[rho[a]] for a in reps]

    cyc = []
    for a in range(nA):
        if a in tree:
            continue
        z = path(head[a], tail[a])
        z[a] += 1
        cyc.append(g_of(z))
    vbar = [cov.vid[rc(w)] for w in cov.words]
    anti = [g_of(path(v, vbar[v])) for v in range(nv)]
    var = cov.constraints()
    for g in cyc + anti:
        net = collections.Counter()
        for j, gj in enumerate(g):
            for v, c in var[j].items():
                net[v] += c * gj
        assert not any(net.values()), "a generator is not in Lambda^-"
    return cyc, anti


def lambda_rank(cov):
    var = cov.constraints()
    J = len(var)
    rows = []
    verts = sorted({v for d in var for v in d})
    for v in verts:
        rows.append([var[j].get(v, 0) for j in range(J)])
    r, _ = diag_invariants(rows)
    return J - r


# ------------------------------------------------------------------ main

def main():
    t0 = time.time()
    quick = "--quick" in sys.argv
    S = DF.load_phix()
    print("=== Bio 17: inversions at inverted repeats inside the double-stranded fibre ===")
    print(f"    phiX174 NC_001422.1, {len(S)} bp circular; python {sys.version.split()[0]}\n",
          flush=True)
    rnd = random.Random(20260924)

    # ---------------------------------------------------------- random small cases
    cases = []
    tries = 0
    while len(cases) < 24 and tries < 4000:
        tries += 1
        N = rnd.randrange(16, 34)
        k = rnd.randrange(3, 6)
        T = "".join(rnd.choice("ACGT") for _ in range(N))
        cov = Cover(T, k)
        if not cov.reps or cov.nv < 2:
            continue
        pts = enumerate_box(cov, cap=20000)
        if len(pts) > 20000:
            continue
        cnt = sum(cov.n_seq(cov.point(g)) for g in pts)
        if not (3 <= cnt <= 3000):
            continue
        # we want cases in which the strand choice matters: more than one lattice point is used
        used = sum(1 for g in pts if cov.n_seq(cov.point(g)))
        if used < 2:
            continue
        cases.append((T, k, cov, cnt))
    print(f"    {len(cases)} random test sequences (length 16-33, k = 3-5, 3 <= #DS <= 3000, at least "
          f"two lattice points carrying reconstructions)\n", flush=True)

    # (I) the inversion theorem
    print("=== (I) inversions flanked by inverted repeats preserve the double-stranded spectrum ===")
    okI, n_ir, n_ctrl, ctrl_keep = True, 0, 0, 0
    for T, k, cov, cnt in cases:
        D = ds_spec(T, k + 1)
        for u in ir_inversions(T, k):
            n_ir += 1
            okI &= ds_spec(u, k + 1) == D
        for _ in range(40):
            i, L = rnd.randrange(len(T)), rnd.randrange(2, len(T))
            t = T[i:] + T[:i]
            u = rc(t[:L]) + t[L:]
            n_ctrl += 1
            ctrl_keep += ds_spec(u, k + 1) == D
    print(f"    random sequences: {n_ir} IR-inversions, all preserve the spectrum: {okI}; "
          f"control: {ctrl_keep}/{n_ctrl} random unflanked inversions preserve it")
    for k in (12, 11, 10):
        D = ds_spec(S, k + 1)
        inv = ir_inversions(S, k)
        good = sum(1 for u in inv if ds_spec(u, k + 1) == D)
        okI &= good == len(inv)
        print(f"    phiX174 k={k}: {len(inv)} inversions flanked by inverted repeats of length {k}; "
              f"{good} preserve the double-stranded {k+1}-spectrum", flush=True)
    row("I", "an inversion flanked by an inverted repeat w ... rc(w) of length k preserves "
             "spec(T) + spec(rc T) at word length k", okI)

    # (C) three counts
    print("\n=== (C) DS(S,k): brute force = Bio 8 lattice sum = closure under the moves ===")
    okC, okG, nchecked = True, True, 0
    for idx, (T, k, cov, cnt) in enumerate(cases):
        B = brute_ds(T, k, nodecap=2_000_000)
        tk = Tok(cov)
        st = tk.circuit_S()
        assert tk.valid(st) and canon_str(tk.spell(st)) == canon_str(T)
        cl, bad = tk.closure(st, "TIR")
        cl_str = {canon_str(tk.spell(t)) for t in cl} if cl is not None else None
        if B is None:
            print(f"    N={len(T):2d} k={k}: brute force too large, skipped; lattice sum {cnt:>5}, "
                  f"closure {len(cl_str) if cl_str is not None else '?':>5}", flush=True)
            okC &= cl_str is not None and len(cl_str) == cnt
            continue
        nchecked += 1
        agree = cl_str is not None and len(B) == cnt and cl_str == B and bad == 0
        okC &= agree
        print(f"    N={len(T):2d} k={k}: brute force {len(B):>5}, lattice sum "
              f"{cnt:>5}, closure {len(cl_str) if cl_str is not None else '?':>5}"
              f"{'' if agree else '   MISMATCH'}", flush=True)
    okC &= nchecked >= 12
    print(f"    {nchecked} cases checked against exhaustive enumeration")
    row("C", "on every random case the exhaustive enumeration, Bio 8's lattice sum of N_seq and the "
             "closure of S under transpositions, IR-inversions and reverse complement give the same "
             "set (resp. number) of double-stranded reconstructions", okC)

    # (G) generation, and what each move family contributes
    print("\n=== (G) which moves generate DS(S,k) ===")
    fam = [("T", "transpositions"), ("TR", "transp. + rev. compl."), ("TI", "transp. + IR-inversions"),
           ("TIR", "all three")]
    tally = collections.Counter()
    for idx, (T, k, cov, cnt) in enumerate(cases):
        tk = Tok(cov)
        st = tk.circuit_S()
        sizes = {}
        for kinds, _ in fam:
            cl, bad = tk.closure(st, kinds)
            sizes[kinds] = len(cl)
            okG &= bad == 0
        okG &= sizes["TIR"] == cnt
        # transpositions alone = the single-stranded fibre = N_seq at the lattice point of S
        okG &= sizes["T"] == cov.n_seq(cov.point(tk.point(st)))
        for kinds, _ in fam:
            tally[kinds] += sizes[kinds] == cnt
    print("    closure equals all of DS(S,k), out of", len(cases), "random cases:")
    for kinds, name in fam:
        print(f"      {name:26s} {tally[kinds]:3d}")
    for k, want in ((12, 10), (11, 864)):
        cov = Cover(S, k)
        tk = Tok(cov)
        st = tk.circuit_S()
        assert tk.valid(st)
        line = []
        for kinds, name in fam:
            t1 = time.time()
            cl, bad = tk.closure(st, kinds)
            okG &= bad == 0
            line.append(f"{name} {len(cl)}")
            if kinds == "TIR":
                okG &= len(cl) == want
                pts = collections.Counter(tk.point(t) for t in cl)
                sample = list(cl)[:5]
                okG &= all(ds_spec(tk.spell(t), k + 1) == ds_spec(S, k + 1) for t in sample)
        print(f"    phiX174 k={k} (#DS = {want}): " + "; ".join(line)
              + f"; lattice points reached {len(pts)}", flush=True)
    row("G", "the closure of S under transpositions, inversions at inverted repeats and reverse "
             "complementation is the whole double-stranded fibre (random cases; phiX174 k = 11, 12)",
        okG)

    # (L) the inversion lattice
    print("\n=== (L) the lattice spanned by inversion vectors ===")
    okL = True
    print("      case                rank Lambda^-   rank/index (1-rho)Z_1     rank/index M")
    for label, cov in ([(f"phiX174 k={k}", Cover(S, k)) for k in (12, 11, 10)]
                       + [(f"random N={len(T)} k={k}", cov) for T, k, cov, _ in cases[:8]]):
        r = lambda_rank(cov)
        cyc, anti = inversion_lattice(cov)
        rc_, invc = diag_invariants(cyc)
        rm, invm = diag_invariants(cyc + anti)
        okL &= rm == r
        ic = 1
        for x in invc:
            ic *= x
        im = 1
        for x in invm:
            im *= x
        print(f"      {label:20s} {r:8d}      {rc_:5d} / {ic:<8d} {str(invc)[:22]:22s}"
              f" {rm:5d} / {im:<6d} {str(invm)[:20]}", flush=True)
    row("L", "the inversion vectors (1 - rho) z of paths z from v to vbar span a sublattice of full "
             "rank in Lambda^-; its index is reported", okL)

    # (H) the Tate dichotomy on random double covers
    print("\n=== (H) [Lambda^- : (1-rho)Z_1] = 2 exactly when there is no palindromic vertex and no "
          "self-complementary content ===")
    r5 = random.Random(5)
    stats, n, skipped, okH = {}, 0, 0, True
    while n < 600:
        N, k = r5.randrange(14, 60), r5.randrange(2, 8)
        T = "".join(r5.choice("ACGT") for _ in range(N))
        cov = Cover(T, k)
        if not cov.reps or cov.nv < 2:
            continue
        try:
            cyc, anti = inversion_lattice(cov)
        except AssertionError:                    # D_k disconnected
            skipped += 1
            continue
        n += 1
        r = lambda_rank(cov)
        rc_, ic = diag_invariants(cyc)
        rm, im = diag_invariants(cyc + anti)
        free = (cov.n_fixV == 0 and cov.n_fixA == 0)
        okH &= rc_ == r and rm == r and im == [] and ic == ([2] if free else [])
        key = ("no palindromic vertex" if cov.n_fixV == 0 else "palindromic vertex",
               "no self-compl. content" if cov.n_fixA == 0 else "self-compl. content")
        s = stats.setdefault(key, [0, 0, 0])
        s[0] += 1
        s[1] += ic != []
        s[2] += im != []
    for key, (c, a, b) in sorted(stats.items()):
        print(f"    {key[0]:22s} {key[1]:24s} {c:4d} covers: index of (1-rho)Z_1 is 2 in {a}, "
              f"index of M > 1 in {b}")
    print(f"    ({skipped} disconnected covers skipped)")
    row("H", "on 600 random connected double covers, [Lambda^- : (1-rho) Z_1] = 2 iff D_k has no "
             "palindromic vertex and no self-complementary content, and = 1 otherwise; the inversion "
             "lattice M is always all of Lambda^-", okH)

    # (W) phiX174 at k = 10
    print("\n=== (W) phiX174 at k = 10 ===")
    cov = Cover(S, 10)
    tk = Tok(cov)
    st = tk.circuit_S()
    okW = tk.valid(st)
    p0 = tk.point(st)
    one = [tk.apply(st, d) for d in tk.moves(st, "I")]
    okW &= all(tk.valid(u) for u in one)
    pts1 = {tk.point(u) for u in one}
    print(f"    {len(one)} inversions at inverted repeats of the genome itself; they reach "
          f"{len(pts1 - {p0})} other lattice points (of the box's 2,481,687,900)")
    if quick:
        print("    random walk skipped (--quick)")
    else:
        for kinds, steps in (("TI", 20000), ("TIR", 20000)):
            t = st
            seen_pts, seen_seq = {p0}, {tk.canon(st)}
            ok = True
            t1 = time.time()
            for s_ in range(steps):
                ds = tk.moves(t, kinds)
                u = tk.apply(t, rnd.choice(ds))
                if not tk.valid(u):
                    ok = False
                    break
                t = u
                seen_pts.add(tk.point(t))
                if s_ % 20 == 0:
                    seen_seq.add(tk.canon(t))
            okW &= ok
            print(f"    random walk, {steps} moves of kinds {kinds}: every state a valid "
                  f"reconstruction {ok}; {len(seen_pts)} distinct lattice points visited, "
                  f"{len(seen_seq)} distinct sequences among {steps // 20} sampled"
                  f"   [{time.time()-t1:.0f}s]", flush=True)
        # spell the final state and check its spectrum at the string level
        okW &= ds_spec(tk.spell(t), 11) == ds_spec(S, 11)
    row("W", "at k = 10 every state reached by the moves is a valid double-stranded reconstruction",
        okW)

    print("\n=== battery summary ===")
    bad = [k for k, ok in ROWS if not ok]
    print(f"  TOTAL {len(ROWS)-len(bad)}/{len(ROWS)} passed" + (f"   FAILURES: {bad}" if bad else ""))
    print(f"  elapsed {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
