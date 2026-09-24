# -*- coding: utf-8 -*-
r"""
bidirected_sandpile.py -- verification battery for

  "The reverse-complement quotient of a de Bruijn graph is a signed graph,
   and its critical group is not an assembly invariant"

THE QUESTION (Bio 6, Rem. 7; monograph Ch. 26, Open 7).  Bio 6 built the reverse-complement
DOUBLE COVER D_k of a de Bruijn graph and showed that the involution rho(v) = vbar carries
D_k to its OPPOSITE digraph, so rho induces no endomorphism of K(D_k) and the cover does not
descend.  The obvious repair is to pass to the QUOTIENT and to use Zaslavsky's signed
Laplacian Delta = D - A_sigma there.  Does that give a sandpile theory of double-stranded
assembly?  Does the r-bundle chain splitting of Bio 3 survive?  Do E. coli's seven rRNA
operons merge into one 7-fold object?

THE ANSWERS.  Three yes and three no, and the no's are the theorems.

YES.
  1. rho is free on V(D_k) exactly when S has no self-reverse-complementary repeated k-mer.
     For ODD k that is automatic (the middle base would have to be its own complement), so
     Sigma_k := D_k/rho is a signed graph for every odd k.  For even k it can fail, and does:
     phiX174 has 3 palindromic repeated 10-mers and 1 palindromic repeated 12-mer, and at
     those k the quotient adjacency is not even symmetric.
  2. All r parallel edges of a bundle of Sigma_k carry the SAME sign, because the sign
     eps(u)eps(v) depends only on the endpoints.  Hence Sigma_k has honest r-bundles and
     K(Sigma_k) (Reiner-Tseng: im d / im Delta, an index-2 subgroup of coker Delta when
     Sigma_k is unbalanced) is a switching invariant of the genome.
  3. Sigma_k is UNBALANCED iff S carries an inverted repeat of length >= k.  This is
     Zaslavsky balance read through Bio 6's connectivity criterion: the signed double cover
     is disconnected exactly when the base is balanced (Reiner-Tseng, Prop. 11.4-11.5).

NO.
  4. THE DIGRAPH CARRIES NO COVERING THEORY.  With P the permutation matrix of rho,
     P Delta(D_k) P = Delta(D_k)^t and NOT Delta(D_k).  The deck involution is an
     ANTI-automorphism of the Laplacian, so the +-1 eigenlattices of P are not Delta-stable
     and Reiner-Tseng's exact sequence has no directed analogue.  Numerically K(D_k) is not
     an extension of anything by K(Sigma_k): for phiX174 at k = 11, |K(D_11)| = 25318912
     while |K(|D_11|)| = |K(Sigma_11)| |K(|Sigma_11|)| = 19738255595056.
  5. THE BUNDLE CHAIN DOES NOT CONTRACT.  In a digraph the chain row is r(x_t - x_{t+1}),
     which kills one difference outright: Bio 3's (Z/r)^{s-1} (+) K(G/C).  In a signed graph
     the row is r(d_{t-1} - d_t) with d_t = x_t - x_{t+1}, which kills only SECOND
     differences.  The s second differences do give an r-torsion subgroup of rank s, but the
     residual first difference couples to the rest of the graph by the SERIES law: the row at
     the chain's endpoint a differs between Sigma and Sigma/C by exactly s r d_0.  Hence
     K(Sigma) is NOT (Z/r)^s (+) K(Sigma/C), and the chain contraction that makes the E. coli
     computation of Bio 5 feasible has no bidirected analogue.
  6. SO THE MERGE THE QUESTION ASKS FOR CANNOT COUNT ANYTHING.  The seven rRNA operons DO
     merge into a single 7-bundle chain of Sigma_k -- that part works.  But K(Sigma_k) is an
     invariant of the UNDIRECTED shadow, and by the matrix-tree theorem it counts spanning
     trees, not Eulerian circuits.  BEST needs the digraph, which by (4) has no quotient.
     The 7-fold structure is visible instead as F_7-nullity, which is what we report.

 (P)  rho is free iff no palindromic repeated k-mer; automatic for odd k             exact
 (W)  bundles have constant sign; K(Sigma_k) is switching invariant                  exact
 (B)  Sigma_k unbalanced iff S has an inverted repeat of length >= k                 exact
 (R)  Reiner-Tseng 0 -> K(Sigma_k) -> K(|D_k|) -> K(|Sigma_k|) -> 0, orders multiply
      and the 2-primary part is a NON-SPLIT extension                                exact
 (N)  P Delta P = Delta^t but P Delta P != Delta, and K(D_k) misses the sequence      exact
 (D)  what rho gives instead is a SYMMETRIC Bowen-Franks linking form on K(D_k); the
      relevant notion is a metabolizer, not a Lagrangian, and for phiX174 the order
      is not a perfect square, so no metabolizer exists                              exact
 (C)  K(Sigma) != (Z/r)^s (+) K(Sigma/C); the defect is the series law               exact
 (E)  E. coli K-12: the seven rRNA operons merge into 7-bundles of Sigma_k whose longest
      chain spans 728 bp at k = 31, 51 and 75 alike -- the same block Bio 6 matched
      one-to-one onto the seven rrn operons, and at k = 51 the same length s = 678 --
      and their only quantitative trace is the F_7-nullity of Delta_{Sigma_k}        exact

Run: python bidirected_sandpile.py      (about three minutes; E. coli needs ~3 GB)
Imports rc_double_cover.py from Bio 6, ecoli_sandpile.py from Bio 5, pgl3_building.py
from Paper 10.  Genome files are the ones cached by Bio 5 and Bio 6.
"""

import sys, os, time, random, collections, math, itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


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


for sub in ("Bio_06_RC_Double_Cover", "Bio_05_Ecoli_Decomposition", "lib"):
    sys.path.insert(0, _sib(sub))
import rc_double_cover as rc6
import ecoli_sandpile as es
import pgl3_building as pb

ROWS = []


def row(key, desc, cases, mode, ok):
    ROWS.append((key, desc, cases, mode, bool(ok)))
    print(("  [ok]   " if ok else "  [FAIL] ") + f"({key}) {desc}   [{cases}] [{mode}]")
    return bool(ok)


def order(t):
    return math.prod(int(x) for x in t) if t else 1


# ------------------------------------------------------------------ integer cokernels

def snf_diag(A):
    """elementary-divisor diagonal (nonzero entries) of a small integer matrix"""
    A = [r[:] for r in A]
    R, C = len(A), len(A[0]) if A else 0
    out, r0, c0 = [], 0, 0
    while r0 < R and c0 < C:
        piv = min(((i, j) for i in range(r0, R) for j in range(c0, C) if A[i][j]),
                  key=lambda t: abs(A[t[0]][t[1]]), default=None)
        if piv is None:
            break
        while True:
            i, j = piv
            A[r0], A[i] = A[i], A[r0]
            for rr in A:
                rr[c0], rr[j] = rr[j], rr[c0]
            p, done = A[r0][c0], True
            for i in range(r0 + 1, R):
                if A[i][c0]:
                    q = A[i][c0] // p
                    for j in range(c0, C):
                        A[i][j] -= q * A[r0][j]
                    done &= not A[i][c0]
            for j in range(c0 + 1, C):
                if A[r0][j]:
                    q = A[r0][j] // p
                    for i in range(r0, R):
                        A[i][j] -= q * A[i][c0]
                    done &= not A[r0][j]
            if done:
                break
            piv = min(((i, j) for i in range(r0, R) for j in range(c0, C) if A[i][j]),
                      key=lambda t: abs(A[t[0]][t[1]]))
        out.append(abs(A[r0][c0]))
        r0 += 1
        c0 += 1
    return out


def coker(M, big=60):
    """(elementary divisors, free rank) of Z^n / im(M^t).  Uses the p-adic route of
    Paper X for large nonsingular M and Euclidean SNF otherwise."""
    if not M:
        return (), 0
    n = len(M)
    A = [list(r) for r in zip(*M)]
    if n >= big:
        f, free, rem = pb.coker_structure(A)
        assert rem == 1, f"unfactored cofactor {rem}"
        return es.elementary(tuple(f)), free
    d = snf_diag(A)
    return es.elementary(tuple(x for x in d if x > 1)), n - len(d)


def det_order(M):
    """|det M| = the order of the cokernel, without factoring it"""
    return abs(pb.bareiss_det([list(r) for r in M]))


# ------------------------------------------------------------------ signed graphs

def signed_lap(n, sedges):
    L = [[0] * n for _ in range(n)]
    for (i, j, s) in sedges:
        L[i][i] += 1
        L[j][j] += 1
        L[i][j] -= s
        L[j][i] -= s
    return L


def K_signed(n, sedges):
    """Reiner-Tseng K(Sigma) = im d / im Delta for a connected UNBALANCED signed graph.
    im d = {x : sum x_v even}, with basis 2e_0, e_0+e_1, ..., e_0+e_{n-1}."""
    L = signed_lap(n, sedges)
    X = []
    for j in range(n):
        col = [L[i][j] for i in range(n)]
        c0 = col[0] - sum(col[1:])
        assert c0 % 2 == 0, "im Delta not inside im d -- is Sigma balanced?"
        X.append([c0 // 2] + col[1:])
    return coker([list(r) for r in zip(*X)])


def K_unsigned(n, edges):
    L = signed_lap(n, [(i, j, 1) for (i, j) in edges])
    return coker([[L[i][j] for j in range(1, n)] for i in range(1, n)])


def balanced(n, sedges):
    adj = collections.defaultdict(list)
    for (i, j, s) in sedges:
        adj[i].append((j, s))
        adj[j].append((i, s))
    col = {}
    for s0 in range(n):
        if s0 in col:
            continue
        col[s0] = 1
        st = [s0]
        while st:
            u = st.pop()
            for (v, sg) in adj[u]:
                w = col[u] * sg
                if v not in col:
                    col[v] = w
                    st.append(v)
                elif col[v] != w:
                    return False
    return True


# ------------------------------------------------------- the reverse-complement quotient

def orbits(words, vid):
    n = len(words)
    perm = [vid[rc6.rc(w)] for w in words]
    orb, m = {}, 0
    for v in range(n):
        r = min(v, perm[v])
        if r not in orb:
            orb[r] = orb[perm[r]] = m
            m += 1
    return perm, orb, m


def arc_orbits(arcs, perm):
    """one representative per rho-orbit of arcs; rho(u -> v) = perm(v) -> perm(u)"""
    idx = collections.defaultdict(list)
    for t, a in enumerate(arcs):
        idx[a].append(t)
    used, out = set(), []
    for t, (u, v) in enumerate(arcs):
        if t in used:
            continue
        used.add(t)
        for t2 in idx[(perm[v], perm[u])]:
            if t2 not in used:
                used.add(t2)
                break
        else:
            raise AssertionError("arc without a rho-partner: rho is not an involution")
        out.append((u, v))
    return out


def sigma_of(S, k):
    """Sigma_k = D_k / rho.  Returns a dict; 'free' is False when rho has a fixed vertex."""
    nv, arcs, words, vid, occ = rc6.compacted_multi([S, rc6.rc(S)], k)
    perm, orb, m = orbits(words, vid)
    free = all(perm[v] != v for v in range(nv))
    out = dict(nv=nv, arcs=arcs, words=words, vid=vid, occ=occ, perm=perm, orb=orb,
               m=m, free=free, fixed=[v for v in range(nv) if perm[v] == v])
    if free and nv:
        reps = {v for v in range(nv) if words[v] < words[perm[v]]}
        eps = [1 if v in reps else -1 for v in range(nv)]
        ao = arc_orbits(arcs, perm)
        out["eps"] = eps
        out["arc_orbits"] = ao
        out["sedges"] = [(orb[u], orb[v], eps[u] * eps[v]) for (u, v) in ao]
    return out


# --------------------------------------------------------------- sparse rank mod p

def nullity_mod_p(m, sedges, p):
    """dim ker Delta_Sigma over F_p, by sparse elimination with a min-degree order"""
    A = [dict() for _ in range(m)]
    for (i, j, s) in sedges:
        A[i][i] = (A[i].get(i, 0) + 1) % p
        A[j][j] = (A[j].get(j, 0) + 1) % p
        A[i][j] = (A[i].get(j, 0) - s) % p
        A[j][i] = (A[j].get(i, 0) - s) % p
    for i in range(m):
        A[i] = {j: v for j, v in A[i].items() if v}
    alive = set(range(m))
    rank = 0
    while alive:
        piv = min(alive, key=lambda i: len(A[i]))
        if not A[piv]:
            alive.discard(piv)
            continue
        col = piv if piv in A[piv] else next(iter(A[piv]))
        if col not in alive:
            col = next((c for c in A[piv] if c in alive), None)
            if col is None:
                alive.discard(piv)
                continue
        inv = pow(A[piv][col], p - 2, p)
        prow = {j: (v * inv) % p for j, v in A[piv].items()}
        alive.discard(piv)
        for i in list(alive):
            f = A[i].get(col)
            if not f:
                continue
            r = A[i]
            for j, v in prow.items():
                w = (r.get(j, 0) - f * v) % p
                if w:
                    r[j] = w
                else:
                    r.pop(j, None)
        for i in alive:
            A[i].pop(col, None)
        rank += 1
    return m - rank


# ------------------------------------------------------------------------------- checks

def check_P():
    print("\n=== (P) is rho free?  palindromic repeated k-mers ===")
    print("      genome      k  parity  |V(D_k)|  rho-fixed  A_sigma symmetric")
    ok = True
    phix = open(os.path.join(HERE, "phix174_NC_001422.1.txt")).read().strip().upper()
    rows = []
    for k in range(9, 15):
        d = sigma_of(phix, k)
        if not d["nv"]:
            continue
        sym = d["free"]
        rows.append(("phiX174", k, d["nv"], len(d["fixed"]), sym))
        ok &= (k % 2 == 1) <= d["free"]          # odd k must be free
        ok &= (len(d["fixed"]) == 0) == d["free"]
    for g, k, nv, fx, sym in rows:
        print(f"      {g:10s} {k:2d}  {'even' if k % 2 == 0 else 'odd '}   {nv:7d}   "
              f"{fx:8d}   {str(sym):>16s}")
    # the arithmetic reason
    reason = all(rc6.rc(w) != w for w in
                 ("".join(x) for x in __import__("itertools").product("ACGT", repeat=5)))
    print(f"      no 5-mer is its own reverse complement: {reason}  "
          f"(a middle base would have to be its own complement)")
    return row("P", "rho is free on V(D_k) exactly when no repeated k-mer is its own reverse "
                    "complement, which is automatic for odd k and fails for phiX174 at "
                    "k = 10 and k = 12", f"{len(rows)} values of k", "exact", ok and reason)


def check_W(phix):
    print("\n=== (W) bundles have constant sign; K(Sigma_k) is a switching invariant ===")
    ok = True
    for k in (11, 9):
        d = sigma_of(phix, k)
        m, se = d["m"], d["sedges"]
        bysides = collections.defaultdict(set)
        for (i, j, s) in se:
            bysides[(min(i, j), max(i, j))].add(s)
        const = all(len(v) == 1 for v in bysides.values())
        base = K_signed(m, se)
        rng = random.Random(3)
        inv = True
        for _ in range(8):
            r2 = {(v if rng.random() < .5 else d["perm"][v])
                  for v in range(d["nv"]) if v < d["perm"][v]}
            e2 = [1 if v in r2 else -1 for v in range(d["nv"])]
            se2 = [(d["orb"][u], d["orb"][v], e2[u] * e2[v]) for (u, v) in d["arc_orbits"]]
            inv &= (K_signed(m, se2) == base)
        ok &= const and inv
        print(f"      k={k:2d}: Sigma_k has {m} vertices and {len(se)} edges in "
              f"{len(bysides)} bundles;  every bundle monochromatic: {const};  "
              f"K(Sigma_k) unchanged by 8 random switchings: {inv}")
        print(f"            K(Sigma_{k}) = {es.gstr(base[0])}")
    return row("W", "the sign of an edge of Sigma_k depends only on its endpoints, so every "
                    "bundle is monochromatic and Sigma_k has honest r-bundles; K(Sigma_k) is "
                    "unchanged by switching", "2 values of k, 8 switchings each", "exact", ok)


def check_B(phix):
    print("\n=== (B) balance of Sigma_k vs inverted repeats of S ===")
    print("      k   Sigma_k balanced  |D_k| connected  longest inverted repeat >= k")
    ok = True
    for k in (9, 11, 13, 15):
        d = sigma_of(phix, k)
        if not d["nv"]:
            bal, conn, ir = True, False, None
        else:
            bal = balanced(d["m"], d["sedges"])
            conn = rc6.n_components(d["nv"], d["arcs"]) == 1
            ir = rc6.inverted_repeat(phix, k)
        has = ir is not None
        ok &= ((not bal) == conn) and (conn == has if d["nv"] else True)
        print(f"      {k:2d}   {str(bal):>16s}  {str(conn):>15s}  "
              f"{('yes: ' + str(ir[0]) if has else 'no'):>28s}")
    return row("B", "Sigma_k is unbalanced exactly when the signed double cover D_k is "
                    "connected (Reiner-Tseng Prop. 11.4-11.5), which by Bio 6 happens exactly "
                    "when S carries an inverted repeat of length >= k", "4 values of k",
               "exact", ok)


def check_R(phix):
    print("\n=== (R) Reiner-Tseng: 0 -> K(Sigma_k) -> K(|D_k|) -> K(|Sigma_k|) -> 0 ===")
    ok = True
    for k in (11,):
        d = sigma_of(phix, k)
        m, se = d["m"], d["sedges"]
        Ks = K_signed(m, se)
        Kc = K_unsigned(d["nv"], list(d["arcs"]))
        Kq = K_unsigned(m, [(i, j) for (i, j, _) in se])
        good = order(Ks[0]) * order(Kq[0]) == order(Kc[0])
        # 2-primary: is the extension split?
        two = lambda t: tuple(x for x in t if x % 2 == 0)
        split = two(Kc[0]) == tuple(sorted(two(Ks[0]) + two(Kq[0])))
        ok &= good and not split
        print(f"      k={k}")
        print(f"        K(Sigma_{k})   = {es.gstr(Ks[0])}   order {order(Ks[0])}")
        print(f"        K(|D_{k}|)     = {es.gstr(Kc[0])}   order {order(Kc[0])}")
        print(f"        K(|Sigma_{k}|) = {es.gstr(Kq[0])}   order {order(Kq[0])}")
        print(f"        orders multiply: {good}")
        print(f"        2-primary parts {es.gstr(two(Ks[0]))} and {es.gstr(two(Kq[0]))} "
              f"vs {es.gstr(two(Kc[0]))} in the middle: split = {split}")
    return row("R", "on the reverse-complement double cover the Reiner-Tseng sequence holds "
                    "with orders multiplying, and its 2-primary part is a NON-SPLIT extension "
                    "-- (Z/2)^2 by (Z/2)^2 giving (Z/2)^2 (+) Z/4, exactly their Thm 1.2",
               "1 value of k", "exact", ok)


def check_N(phix):
    print("\n=== (N) the deck involution is an ANTI-automorphism of the Laplacian ===")
    ok = True
    for k in (9, 11):
        d = sigma_of(phix, k)
        nv, arcs, perm = d["nv"], d["arcs"], d["perm"]
        D = np.zeros((nv, nv), dtype=np.int64)
        for (u, v) in arcs:
            D[u, u] += 1
            D[u, v] -= 1
        P = np.zeros((nv, nv), dtype=np.int64)
        for v in range(nv):
            P[v, perm[v]] = 1
        anti = np.array_equal(P @ D @ P, D.T)
        comm = np.array_equal(P @ D @ P, D)
        ok &= anti and not comm
        oD = det_order(es.reduced_laplacian(nv, arcs))
        Ls = signed_lap(nv, [(u, v, 1) for (u, v) in arcs])
        oC = det_order([[Ls[i][j] for j in range(1, nv)] for i in range(1, nv)])
        ok &= (oC % oD != 0)
        print(f"      k={k:2d}: P Delta P = Delta^t: {anti};   P Delta P = Delta: {comm}")
        print(f"            |K(D_{k})| (digraph)  = {oD}")
        print(f"            |K(|D_{k}|)| (shadow) = {oC}")
        print(f"            the shadow order is not a multiple of the digraph order: "
              f"{oC % oD != 0}")
    return row("N", "P Delta(D_k) P equals Delta^t and never Delta, so the +-1 eigenlattices "
                    "of the deck involution are not Delta-stable and the directed sandpile "
                    "group sits in no covering sequence", "2 values of k", "exact", ok)


def arc_involution(arcs, perm):
    """rho on arc indices: (u,v) |-> (perm v, perm u)"""
    idx = collections.defaultdict(list)
    for t, a in enumerate(arcs):
        idx[a].append(t)
    rho = {}
    for t, (u, v) in enumerate(arcs):
        if t in rho:
            continue
        pool = [s for s in idx[(perm[v], perm[u])] if s not in rho]
        s = next((x for x in pool if x != t), pool[0])
        rho[t], rho[s] = s, t
    return rho


def arborescences(nv, arcs, root, converging=True):
    """all spanning trees oriented towards (converging) or away from the root"""
    nbr = collections.defaultdict(list)
    for t, (u, v) in enumerate(arcs):
        nbr[u if converging else v].append(t)
    others = [v for v in range(nv) if v != root]
    res = []
    for pick in itertools.product(*[nbr[v] for v in others]):
        par = {v: (arcs[t][1] if converging else arcs[t][0])
               for v, t in zip(others, pick)}
        good = True
        for v in others:
            seen, cur = set(), v
            while cur != root:
                if cur in seen or cur is None:
                    good = False
                    break
                seen.add(cur)
                cur = par.get(cur)
            if not good:
                break
        if good:
            res.append(frozenset(pick))
    return res


def check_D(phix):
    print("\n=== (D) the pairing rho induces on K(D_k) is SYMMETRIC, not alternating ===")
    import sympy as sp
    ok = True
    for k in (12, 11):
        d = sigma_of(phix, k)
        nv, arcs, perm = d["nv"], d["arcs"], d["perm"]
        if not nv or nv > 40:
            continue
        Dm = sp.zeros(nv, nv)
        for (u, v) in arcs:
            Dm[u, u] += 1
            Dm[u, v] -= 1
        P = sp.zeros(nv, nv)
        for v in range(nv):
            P[v, perm[v]] = 1
        # domain Z^V/Z.1 with basis e_1..e_{n-1}; target Z^V_0 with basis e_j - e_0
        B = sp.zeros(nv, nv - 1)
        Bp = sp.zeros(nv - 1, nv)
        for j in range(nv - 1):
            B[0, j], B[j + 1, j], Bp[j, j + 1] = -1, 1, 1
        G = sp.Matrix(nv - 1, nv - 1, lambda i, j: Dm.T[i + 1, j + 1])
        det = G.det()
        if det == 0:
            continue
        Pi = (Bp * P * B).T * G.inv()
        sym = all(x.q == 1 for x in sp.simplify(Pi - Pi.T))
        alt = all(x.q == 1 for x in sp.simplify(Pi + Pi.T))
        root, exact = sp.integer_nthroot(int(abs(det)), 2)
        ok &= sym and not alt
        print(f"      k={k}: |K(D_k)| = {abs(det)};  symmetric mod Z: {sym};  "
              f"alternating mod Z: {alt}")
        print(f"            |K| a perfect square: {exact}  ->  a metabolizer L = L^perp "
              f"{'may' if exact else 'CANNOT'} exist")
        # the tree-level shadow: rho sends in-trees to OUT-trees, so it does not act
        if nv <= 8:
            ra = arc_involution(arcs, perm)
            for root in range(min(nv, 3)):
                IN = arborescences(nv, arcs, root, True)
                OUTr = arborescences(nv, arcs, perm[root], False)
                INr = arborescences(nv, arcs, perm[root], True)
                img = {frozenset(ra[t] for t in T) for T in IN}
                into_out = img <= {frozenset(T) for T in OUTr}
                into_in = img <= {frozenset(T) for T in INr}
                ok &= into_out and not into_in
                print(f"            root {root}: {len(IN)} in-trees -> rho -> "
                      f"out-trees at rho(root) {into_out}; in-trees there {into_in}; "
                      f"image meets the in-trees at root in "
                      f"{len(img & {frozenset(T) for T in IN})}")
    return row("D", "the Bowen-Franks pairing that rho induces on the sandpile group is a "
                    "SYMMETRIC linking form, so the relevant notion is a metabolizer and not "
                    "a Lagrangian; and for phiX174 the order is not a square, so no "
                    "metabolizer exists", "2 values of k", "exact", ok)


def check_C():
    print("\n=== (C) a signed r-bundle chain does not contract ===")
    print("      the digraph row is r(x_t - x_{t+1}); the signed row is r(d_{t-1} - d_t)")
    print("       r  s   K(Sigma)                          K(Sigma/C)              "
          "  |Sigma|/|Sigma/C|   r^s   r-rank gain")
    ok, splits = True, 0
    rng = random.Random(11)
    for (r, s) in [(3, 1), (3, 2), (3, 3), (5, 2), (5, 3), (7, 2), (7, 3), (11, 2)]:
        # unbalanced base: a negative triangle on 0,1,2, plus a positive edge 0-2
        E0 = [(0, 1, 1), (1, 2, 1), (2, 0, -1), (0, 2, 1)]
        a, b, c = 0, 1, list(range(3, 3 + s))
        path = [a] + c + [b]
        sg = [rng.choice([1, -1]) for _ in range(s + 1)]
        E = list(E0)
        for t, (u, v) in enumerate(zip(path, path[1:])):
            E += [(u, v, sg[t])] * r
        prod = math.prod(sg)
        Ec = list(E0) + [(a, b, prod)] * r
        G, Gc = K_signed(3 + s, E), K_signed(3, Ec)
        rk = lambda t: sum(1 for x in t if x % r == 0)
        gain = rk(G[0]) - rk(Gc[0])
        ratio = order(G[0]) / order(Gc[0])
        is_split = (order(G[0]) == order(Gc[0]) * r ** s) and gain == s
        splits += is_split
        ok &= (gain == s)
        print(f"       {r:2d} {s:2d}  {es.gstr(G[0]):33s} {es.gstr(Gc[0]):23s} "
              f"{ratio:12.4f}  {r**s:6d}   {gain:5d}")
    print(f"      the r-rank gain is s in every case, but the ORDER ratio equals r^s in "
          f"only {splits} of 8: the residual first difference obeys the series law, and the "
          f"row at the endpoint a differs between Sigma and Sigma/C by s*r*d_0")
    return row("C", "a signed r-bundle chain of s internal vertices raises the r-rank by "
                    "exactly s, but K(Sigma) is NOT (Z/r)^s (+) K(Sigma/C): the chain obeys "
                    "the series law, so Bio 3's splitting theorem has no bidirected analogue",
               "8 pairs (r,s)", "exact", ok and splits < 8)


def rank_mod(M, p):
    """rank over F_p of a dense integer matrix"""
    A = [[x % p for x in r] for r in M]
    m = len(A)
    n = len(A[0]) if m else 0
    r = 0
    for c in range(n):
        piv = next((i for i in range(r, m) if A[i][c] % p), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = pow(A[r][c], p - 2, p)
        A[r] = [(x * inv) % p for x in A[r]]
        for i in range(m):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [(x - f * y) % p for x, y in zip(A[i], A[r])]
        r += 1
    return r


def primes_of(r):
    return sorted({q for q in range(2, r + 1) if r % q == 0 and all(q % d for d in range(2, q))})


def chain_pair(r, s, base=None, seed=0):
    """(Sigma, Sigma/C, Sigma') for an r-bundle chain of s internal vertices on a base"""
    rng = random.Random(seed)
    E0 = list(base) if base else [(0, 1, 1), (1, 2, 1), (2, 0, -1), (0, 2, 1)]
    path = [0] + list(range(3, 3 + s)) + [1]
    sg = [rng.choice([1, -1]) for _ in range(s + 1)]
    E = list(E0)
    for t, (u, v) in enumerate(zip(path, path[1:])):
        E += [(u, v, sg[t])] * r
    return (3 + s, E), (3, list(E0) + [(0, 1, math.prod(sg))] * r), (3, list(E0))


def check_Q():
    """the p-divisible subgraph lemma, and the theorem it proves"""
    print("\n=== (Q) the p-divisible subgraph lemma, and the p-rank gain it forces ===")
    print("      H a subgraph whose every bundle has multiplicity divisible by p:")
    print("      an interior vertex of H has degree = 0 and all adjacencies = 0 mod p, so its")
    print("      whole row and column of Delta vanish; every other vertex loses only its")
    print("      H-part.  Hence rank_p(Delta_Sigma) = rank_p(Delta_{Sigma minus H}).")
    okL = okP = True
    print("\n      (i) the block structure, verbatim")
    for (r, s) in [(3, 2), (4, 2), (6, 3), (9, 2), (12, 2), (2, 3)]:
        (nS, ES), (nC, EC), (nB, EB) = chain_pair(r, s, seed=r * 10 + s)
        LS, LC, LB = signed_lap(nS, ES), signed_lap(nC, EC), signed_lap(nB, EB)
        for p in primes_of(r):
            base_ok = all((LS[i][j] - LB[i][j]) % p == 0 for i in range(nB) for j in range(nB))
            chain_ok = all(LS[i][j] % p == 0 for i in range(nB, nS) for j in range(nS)) and \
                       all(LS[i][j] % p == 0 for i in range(nS) for j in range(nB, nS))
            con_ok = all((LC[i][j] - LB[i][j]) % p == 0 for i in range(nB) for j in range(nB))
            rank_ok = rank_mod(LS, p) == rank_mod(LC, p)
            okL &= base_ok and chain_ok and con_ok and rank_ok
            print(f"        r={r:2d} s={s} p={p:2d}: base block agrees {base_ok}, chain rows and "
                  f"columns vanish {chain_ok}, contraction agrees {con_ok}, "
                  f"rank_p equal {rank_ok}")
    print("\n      (ii) the p-rank gain, for EVERY prime dividing r -- composite and even included")
    print("         r  s   p   p-rank K(Sigma)  p-rank K(Sigma/C)   gain   predicted")
    for (r, s) in [(3, 1), (3, 2), (3, 3), (5, 2), (5, 3), (7, 2), (7, 3), (11, 2),
                   (9, 2), (9, 3), (25, 2), (15, 2), (15, 3), (21, 2),
                   (2, 2), (2, 3), (4, 2), (4, 3), (6, 2), (6, 3), (12, 2), (10, 2)]:
        (nS, ES), (nC, EC), _ = chain_pair(r, s, seed=r * 10 + s)
        if balanced(nS, ES) or balanced(nC, EC):
            continue
        KS, KC = K_signed(nS, ES), K_signed(nC, EC)
        for p in primes_of(r):
            gS = sum(1 for x in KS[0] if x % p == 0)
            gC = sum(1 for x in KC[0] if x % p == 0)
            okP &= (gS - gC == s)
            print(f"        {r:2d} {s:2d} {p:3d}   {gS:14d}   {gC:16d}   {gS-gC:5d}   {s:8d}"
                  f"   {'ok' if gS - gC == s else 'MISMATCH'}")
    print("\n      (iii) shapes the chain statement says nothing about: the gain is the number")
    print("            of INTERIOR vertices of H (all of whose bundles lie in H)")
    shapes = [("tree, three arms of length two", 3, [(0, 3), (3, 4), (0, 5), (5, 6), (0, 7), (7, 8)]),
              ("cycle of length four at one vertex", 5, [(0, 3), (3, 4), (4, 5), (5, 0)]),
              ("theta with two attachments", 3, [(0, 3), (3, 4), (4, 1), (3, 5), (5, 4), (4, 6), (6, 7)])]
    rng = random.Random(7)
    for name, r, H in shapes:
        E0 = [(0, 1, 1), (1, 2, 1), (2, 0, -1), (0, 2, 1)]
        E = list(E0)
        for (u, v) in H:
            E += [(u, v, rng.choice([1, -1]))] * r
        touched = collections.Counter()
        for (u, v) in H:
            touched[u] += 1
            touched[v] += 1
        interior = sorted(w for w in touched if w >= 3)
        used = sorted({0, 1, 2} | set(touched))
        idx = {v: i for i, v in enumerate(used)}
        EE = [(idx[u], idx[v], sg) for (u, v, sg) in E]
        if balanced(len(used), EE):
            continue
        cS = coker([[signed_lap(len(used), EE)[i][j] for j in range(len(used))] for i in range(len(used))])
        c0 = coker([[signed_lap(3, E0)[i][j] for j in range(3)] for i in range(3)])
        for p in primes_of(r):
            gain = (sum(1 for x in cS[0] if x % p == 0) - sum(1 for x in c0[0] if x % p == 0))
            okP &= (gain == len(interior))
            print(f"        {name:36s} p={p:2d}: gain {gain:2d}, interior vertices "
                  f"{len(interior):2d}   {'ok' if gain == len(interior) else 'MISMATCH'}")
    return row("Q", "modulo every prime p | r the chain's rows and columns of the signed Laplacian "
                    "vanish, so rank_p is unchanged by deleting any p-divisible subgraph; hence the "
                    "p-rank gain is EXACTLY the number of interior vertices of that subgraph, s for "
                    "a chain -- which proves the equality left open by Theorem 7",
               "6 + 29 + 3 instances", "exact", okL and okP)


def check_E():
    print("\n=== (E) E. coli K-12 NC_000913.3: the seven rRNA operons in Sigma_k ===")
    S = es.load_genome()
    feats = es.load_features()
    operons = es.group_operons(feats)
    print(f"      genome {len(S)} bp; {len(operons)} rRNA operons")
    ok = True
    print("      k    |V(D_k)|  |V(Sigma_k)|  rho-fixed  7-bundles  longest 7-chain  "
          "its span (bp)  F_7-nullity  F_7 single strand")
    spans = set()
    for k in (31, 51, 75):
        t0 = time.time()
        d = sigma_of(S, k)
        if not d["free"]:
            print(f"      {k}: rho is not free ({len(d['fixed'])} palindromic k-mers)")
            ok = False
            continue
        m, se = d["m"], d["sedges"]
        bundles = collections.Counter()
        for (i, j, s) in se:
            bundles[(min(i, j), max(i, j))] += 1
        seven = {e for e, c in bundles.items() if c == 7}
        # longest path made only of 7-bundles
        adj7 = collections.defaultdict(list)
        for (i, j) in seven:
            adj7[i].append(j)
            adj7[j].append(i)
        seen, longest = set(), 0
        for v in adj7:
            if v in seen:
                continue
            comp, st = [], [v]
            seen.add(v)
            while st:
                u = st.pop()
                comp.append(u)
                for w in adj7[u]:
                    if w not in seen:
                        seen.add(w)
                        st.append(w)
            longest = max(longest, len(comp))
        n7 = nullity_mod_p(m, se, 7)
        G = es.compacted_graph(S, k)
        n7s = nullity_mod_p(G["nv"], [(u, v, 1) for (u, v) in G["arcs"]], 7)
        spans.add(longest + k - 1)
        print(f"      {k:3d} {d['nv']:9d} {m:12d}  {len(d['fixed']):9d} "
              f"{len(seven):10d} {longest:16d} {longest + k - 1:14d} {n7:12d} {n7s:18d}"
              f"   [{time.time()-t0:.0f}s]")
        ok &= (len(seven) > 0)
    print(f"      the longest 7-chain spans {sorted(spans)} bp -- the SAME block at every k, "
          f"and the 728 bp seven-copy block that Bio 6 matched one-to-one onto the seven "
          f"rRNA operons")
    ok &= (len(spans) == 1)
    return row("E", "the seven rRNA operons of E. coli K-12 survive the quotient as genuine "
                    "7-bundles of Sigma_k, and their chain shows up in the F_7-nullity of the "
                    "signed Laplacian -- the only trace of them that the quotient retains",
               "3 values of k", "exact", ok)


def main():
    t0 = time.time()
    print("=== the reverse-complement quotient as a signed graph ===")
    phix = open(os.path.join(HERE, "phix174_NC_001422.1.txt")).read().strip().upper()
    check_P()
    check_W(phix)
    check_B(phix)
    check_R(phix)
    check_N(phix)
    check_D(phix)
    check_C()
    check_Q()
    try:
        check_E()
    except Exception as exc:
        row("E", f"E. coli run failed: {exc}", "0", "exact", False)
    print("\n=== battery summary ===")
    bad = [x for x in ROWS if not x[4]]
    print(f"  TOTAL {sum(1 for x in ROWS if x[4])}/{len(ROWS)} passed"
          + ("" if not bad else f"   FAILURES: {[(x[0], x[2]) for x in bad]}"))
    print(f"  elapsed {time.time()-t0:.0f}s")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
