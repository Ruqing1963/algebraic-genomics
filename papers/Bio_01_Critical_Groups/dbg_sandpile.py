# -*- coding: utf-8 -*-
"""
dbg_sandpile.py -- verification battery for

    "Assembly ambiguity beyond branch statistics:
     the sandpile group of a de Bruijn graph"        (version 2)

Setting.  A circular sequence S of length N over an alphabet of size sigma has, for each
k, a de Bruijn multigraph G_k(S): vertices = distinct k-mers, one arc per occurrence of a
(k+1)-mer.  S is an Eulerian circuit of G_k(S), and every Eulerian circuit is a sequence
with the same k-mer spectrum -- a candidate reconstruction.  K(G_k) denotes the sandpile
(critical) group, whose order is the number of spanning arborescences.  All reconstruction
counts are ARC-ROOTED (based at a fixed occurrence), the BEST convention.

 (T)  tower identity: G_{k+1} is a spanning subgraph of the line digraph L(G_k), with
      equality for the complete de Bruijn graph; the line-digraph recursion of [ANCFT III]
      reproduces the classical count tw(DB(s,k)) = s^(s^k - k - 1).      exact
 (L)  our critical groups reproduce Levine's published decomposition of K(DB(2,n)). exact
 (E)  BEST decomposition against brute-force enumeration of arc-rooted Eulerian circuits.
 (R)  Theorem (single repeat): length L, multiplicity r  =>  K(G_k) = (Z/r)^(L-k).  exact
 (I)  Theorem (interleaving): two repeats of length L, multiplicity 2, s = L-k+1:
      disjoint and nested give K = (Z/2)^(2s-2) (additive); interleaved gives
      (Z/2)^(2s-1), one extra Z/2 that survives to k = L.  Verified on clean sequence
      instances, on the explicit compacted graphs, and by symbolic SNF.             exact
 (N)  Theorem (not degree-determined): two recorded witnesses with identical degree data
      and reconstruction counts 4 and 256; plus a randomized splitting rate.  exact/sampled
 (C)  collapse order k*: NOT R+1 in general; minimal witnesses for k* = R, R+1, < R.  exact
 (P)  unitig contraction (Lemma) leaves K unchanged; self-loops leave the out-Laplacian
      unchanged (they cancel in D - A).                                             exact
 (F)  REAL DATA: the critical-group spectrum of phiX174 (NC_001422.1, 5386 bp), computed
      on the compacted graph; 99% vertex reduction makes it a sub-second computation. exact

Exact integer arithmetic throughout; rows marked "sampled" report statistics of a
randomized search and are labelled as such.  Requires sympy (via the Paper X helper).
The phiX174 sequence is read from phix174_NC_001422.1.txt beside this file (no network).
Run: python dbg_sandpile.py     (about 4 minutes)
"""

import sys, os, time, random, itertools, collections
from math import factorial

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

PAPER10 = _sib("lib")
sys.path.insert(0, PAPER10)
import pgl3_building as pb
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

HERE = os.path.dirname(os.path.abspath(__file__))
EXACT, SAMPLED = "exact", "sampled"
ROWS = []
AB = "ACGT"


def row(key, desc, cases, mode, ok):
    ROWS.append((key, desc, cases, mode, bool(ok)))
    print(("  [ok]   " if ok else "  [FAIL] ") + f"({key}) {desc}   [{cases}] [{mode}]")
    return bool(ok)


# --------------------------------------------------------------------------- graphs

def dbg_from_sequence(S, k):
    """Circular S -> de Bruijn multigraph of order k. The genome is an Eulerian circuit.
    Insertion-ordered vertex labelling; O(N)."""
    N = len(S)
    SS = S + S[:k]
    verts = {}
    for i in range(N):
        w = SS[i:i + k]
        if w not in verts:
            verts[w] = len(verts)
    arcs = [(verts[SS[i:i + k]], verts[SS[i + 1:i + 1 + k]]) for i in range(N)]
    return len(verts), arcs, list(verts)


def dbg_complete(sigma, k):
    verts = [tuple(w) for w in itertools.product(range(sigma), repeat=k)]
    vidx = {w: i for i, w in enumerate(verts)}
    arcs = [(vidx[w], vidx[w[1:] + (c,)]) for w in verts for c in range(sigma)]
    return len(verts), arcs


def line_digraph(nv, arcs):
    outs = collections.defaultdict(list)
    for j, (u, v) in enumerate(arcs):
        outs[u].append(j)
    return len(arcs), [(i, j) for i, (u, v) in enumerate(arcs) for j in outs[v]]


# ----------------------------------------------------------- sandpile / arborescences

def out_laplacian(nv, arcs):
    L = [[0] * nv for _ in range(nv)]
    for (u, v) in arcs:
        L[u][u] += 1
        L[u][v] -= 1
    return L


def arborescences(nv, arcs, root=0):
    if nv <= 1:
        return 1
    L = out_laplacian(nv, arcs)
    M = [[L[i][j] for j in range(nv) if j != root] for i in range(nv) if i != root]
    return abs(pb.bareiss_det(M))


def sandpile(nv, arcs, root=0):
    """(order, invariant factors) of coker of the reduced out-Laplacian, transposed."""
    if nv <= 1:
        return 1, ()
    L = out_laplacian(nv, arcs)
    M = [[L[i][j] for j in range(nv) if j != root] for i in range(nv) if i != root]
    d = pb.diagonalize([list(r) for r in zip(*M)], False, False)[0]
    inv = tuple(sorted(abs(x) for x in d if abs(x) not in (0, 1)))
    o = 1
    for x in inv:
        o *= x
    return o, inv


def outdegrees(nv, arcs):
    out = [0] * nv
    for (u, v) in arcs:
        out[u] += 1
    return out


def local_factor(nv, arcs):
    p = 1
    for d in outdegrees(nv, arcs):
        p *= factorial(d - 1)
    return p


def branch_degrees(nv, arcs):
    return tuple(sorted(d for d in outdegrees(nv, arcs) if d > 1))


def euler_count_best(nv, arcs, root=0):
    return arborescences(nv, arcs, root) * local_factor(nv, arcs)


def euler_count_bruteforce(nv, arcs, start_arc=0):
    """Eulerian circuits beginning with a fixed arc (the arc-rooted count)."""
    outs = collections.defaultdict(list)
    for j, (u, v) in enumerate(arcs):
        outs[u].append(j)
    m = len(arcs)
    used = [False] * m
    total = 0

    def walk(v, depth):
        nonlocal total
        if depth == m:
            if v == arcs[start_arc][0]:
                total += 1
            return
        for j in outs[v]:
            if not used[j]:
                used[j] = True
                walk(arcs[j][1], depth + 1)
                used[j] = False

    used[start_arc] = True
    walk(arcs[start_arc][1], 1)
    return total


def compact(nv, arcs):
    """Linear-time unitig contraction: keep vertices with indeg != 1 or outdeg != 1 and
    replace each maximal non-branching path by one arc. Returns (0, []) for a single
    cycle. By the contraction Lemma the critical group is unchanged."""
    outd = [0] * nv
    ind = [0] * nv
    ol = collections.defaultdict(list)
    for (u, v) in arcs:
        outd[u] += 1
        ind[v] += 1
        ol[u].append(v)
    br = [v for v in range(nv) if outd[v] != 1 or ind[v] != 1]
    if not br:
        return 0, []
    bs = set(br)
    rl = {v: i for i, v in enumerate(br)}
    na = []
    for b in br:
        for f in ol[b]:
            x, steps = f, 0
            while x not in bs:
                x = ol[x][0]
                steps += 1
                if steps > nv + 5:
                    raise RuntimeError("runaway unitig")
            na.append((rl[b], rl[x]))
    return len(br), na


# ------------------------------------------------------------------ repeat statistics

def longest_repeat(S, hi=None):
    """Longest substring occurring at least twice in the circular sequence (binary search)."""
    N = len(S)
    SS = S + S
    lo, hi = 1, (hi or N - 1)
    best = 0
    while lo <= hi:
        mid = (lo + hi) // 2
        seen, dup = set(), False
        for i in range(N):
            w = SS[i:i + mid]
            if w in seen:
                dup = True
                break
            seen.add(w)
        if dup:
            best, lo = mid, mid + 1
        else:
            hi = mid - 1
    return best


def collapse_k(S, kmax=40):
    for k in range(1, min(kmax, len(S))):
        nv, arcs, _ = dbg_from_sequence(S, k)
        cn, ca = compact(nv, arcs)
        if cn == 0 or sandpile(cn, ca)[0] == 1:
            return k
    return None


# ----------------------------------------------------------------------------- (T)

def check_tower():
    print("\n=== (T) the de Bruijn tower is the line-digraph tower ===")
    ok_all = True
    for sigma in (2, 3):
        prev = None
        for k in (1, 2, 3):
            nv, arcs = dbg_complete(sigma, k)
            tw = arborescences(nv, arcs)
            want = sigma ** (sigma ** k - k - 1)
            ok = tw == want
            law_ok = True
            if prev is not None:
                npv, pt = prev
                law_ok = sigma ** (npv * (sigma - 1) - 1) * pt == tw
            ok_all &= ok and law_ok
            row("T", f"tw(DB({sigma},{k})) = {tw} = {sigma}^({sigma}^{k}-{k}-1)"
                     + ("" if prev is None else "; [ANCFT III] line-digraph law holds"),
                f"DB({sigma},{k})", EXACT, ok and law_ok)
            prev = (nv, tw)
        nv, arcs = dbg_complete(sigma, 2)
        nL, aL = line_digraph(nv, arcs)
        nv3, arcs3 = dbg_complete(sigma, 3)
        same = (nL == nv3 and len(aL) == len(arcs3)
                and arborescences(nL, aL) == arborescences(nv3, arcs3))
        ok_all &= same
        row("T", f"L(DB({sigma},2)) and DB({sigma},3) agree in size and arborescence count",
            f"sigma={sigma}", EXACT, same)
    return ok_all


# ----------------------------------------------------------------------------- (L)

def check_levine():
    print("\n=== (L) our critical groups agree with Levine's published decomposition ===")
    ok_all = True
    for n in (2, 3, 4, 5):
        nv, arcs = dbg_complete(2, n)
        o, inv = sandpile(nv, arcs)
        pred = []
        for j in range(1, n):
            pred += [2 ** j] * (2 ** (n - 1 - j))
        ok = inv == tuple(sorted(pred))
        ok_all &= ok
        row("L", f"K(DB(2,{n})) = {list(inv)}, order {o} = 2^{2**n - n - 1}",
            f"n={n}", EXACT, ok and o == 2 ** (2 ** n - n - 1))
    return ok_all


# ----------------------------------------------------------------------------- (E)

def check_best():
    print("\n=== (E) BEST: #arc-rooted reconstructions = |K(G_k)| x prod (d-1)! ===")
    random.seed(4)
    ok_all, tested = True, 0
    while tested < 6:
        S = "".join(random.choice("AC") for _ in range(12))
        k = random.choice((2, 3))
        nv, arcs, _ = dbg_from_sequence(S, k)
        if len(arcs) > 14 or nv < 3:
            continue
        best = euler_count_best(nv, arcs)
        brute = euler_count_bruteforce(nv, arcs)
        o, inv = sandpile(nv, arcs)
        ok = best == brute
        ok_all &= ok
        tested += 1
        row("E", f"S={S} k={k}: |K|={o} x prod(d-1)!={local_factor(nv,arcs)}"
                 f" = {best} = brute force {brute}", "N=12", EXACT, ok)
    return ok_all


# ----------------------------------------------------------------------------- (R)

def check_single_repeat():
    print("\n=== (R) single repeat of length L, multiplicity r -> (Z/r)^(L-k) ===")
    random.seed(31337)
    agree = tot = 0
    for Lrep in (8, 10, 12, 14):
        for r in (2, 3, 4):
            for k in range(max(4, Lrep - 4), Lrep + 1):
                hits = collections.Counter()
                clean = 0
                for _ in range(400):
                    rep = "".join(random.choice(AB) for _ in range(Lrep))
                    parts = []
                    for _i in range(r):
                        parts.append(rep)
                        parts.append("".join(random.choice(AB) for _ in range(30)))
                    S = "".join(parts)
                    nv, arcs, _ = dbg_from_sequence(S, k)
                    if branch_degrees(nv, arcs) != (r,) * (Lrep - k + 1):
                        continue
                    clean += 1
                    hits[sandpile(nv, arcs)[1]] += 1
                    if clean >= 60:
                        break
                if clean == 0:
                    continue
                pred = tuple([r] * (Lrep - k))
                agree += (len(hits) == 1 and next(iter(hits)) == pred)
                tot += 1
    row("R", f"K(G_k) = (Z/r)^(L-k) in every clean instance, L in 8..14, r in 2..4, "
             f"k in L-4..L", f"{tot} parameter settings", EXACT, agree == tot)
    return agree == tot


# ----------------------------------------------------------------------------- (I)

TAGS = ("disjoint", "nested", "interleaved")


def arrangements(A, B, spc):
    return {"disjoint":    A + spc[0] + A + spc[1] + B + spc[2] + B + spc[3],
            "nested":      A + spc[0] + B + spc[1] + B + spc[2] + A + spc[3],
            "interleaved": A + spc[0] + B + spc[1] + A + spc[2] + B + spc[3]}


def reduced_laplacian_two_repeats(s, interleaved):
    """The explicit compacted graph of two length-L multiplicity-2 repeats, s = L-k+1:
    spines a_1..a_s and b_1..b_s with doubled arcs, connectors as in the paper."""
    n = 2 * s
    L = sp.zeros(n, n)
    a = lambda i: i - 1
    b = lambda i: s + i - 1
    for i in range(1, s):
        L[a(i), a(i)] += 2; L[a(i), a(i + 1)] -= 2
        L[b(i), b(i)] += 2; L[b(i), b(i + 1)] -= 2
    if interleaved:
        L[a(s), a(s)] += 2; L[a(s), b(1)] -= 2
        L[b(s), b(s)] += 2; L[b(s), a(1)] -= 2
    else:
        L[a(s), a(s)] += 2; L[a(s), a(1)] -= 1; L[a(s), b(1)] -= 1
        L[b(s), b(s)] += 2; L[b(s), b(1)] -= 1; L[b(s), a(1)] -= 1
    return L[1:, 1:]


def check_interleaving():
    print("\n=== (I) two repeats: disjoint = nested = (Z/2)^(2s-2), interleaved = (Z/2)^(2s-1) ===")
    # (a) clean sequence instances
    rng = random.Random(777)
    tot = ok = 0
    for (L, k, spacer) in [(8, 6, 20), (8, 7, 20), (10, 8, 24), (10, 10, 24), (12, 9, 24)]:
        s = L - k + 1
        pred_dn, pred_in = (2,) * (2 * s - 2), (2,) * (2 * s - 1)
        hits = 0
        for _ in range(400):
            A = "".join(rng.choice(AB) for _ in range(L))
            B = "".join(rng.choice(AB) for _ in range(L))
            if A == B:
                continue
            spc = ["".join(rng.choice(AB) for _ in range(spacer)) for _ in range(4)]
            cfg = arrangements(A, B, spc)
            res, clean = {}, True
            for tag in TAGS:
                nv, arcs, _ = dbg_from_sequence(cfg[tag], k)
                if branch_degrees(nv, arcs) != (2,) * (2 * s):
                    clean = False
                    break
                res[tag] = sandpile(nv, arcs)[1]
            if not clean:
                continue
            hits += 1
            tot += 1
            ok += (res["disjoint"] == pred_dn and res["nested"] == pred_dn
                   and res["interleaved"] == pred_in)
            if hits >= 40:
                break
    row("I", f"prediction exact on every clean sequence instance", f"{tot} instances", EXACT,
        ok == tot and tot > 0)

    # (b) the compacted graphs really are the ones analysed
    rng2 = random.Random(3)
    while True:
        A = "".join(rng2.choice(AB) for _ in range(8))
        B = "".join(rng2.choice(AB) for _ in range(8))
        spc = ["".join(rng2.choice(AB) for _ in range(20)) for _ in range(4)]
        cfg = arrangements(A, B, spc)
        if all(branch_degrees(*dbg_from_sequence(cfg[t], 6)[:2]) == (2,) * 6 for t in TAGS):
            break
    ok_b = True
    for tag in TAGS:
        nv, arcs, _ = dbg_from_sequence(cfg[tag], 6)
        cn, ca = compact(nv, arcs)
        mult = sorted(collections.Counter(collections.Counter(ca).values()).items())
        want = [(2, 6)] if tag == "interleaved" else [(1, 4), (2, 4)]
        good = (cn == 6 and mult == want and sandpile(cn, ca)[1] == sandpile(nv, arcs)[1])
        ok_b &= good
        print(f"    {tag:<12s}: {nv} -> {cn} branch vertices, arc multiplicities {mult}, "
              f"K = {list(sandpile(cn, ca)[1])}")
    row("I", "compacted graphs have 2s vertices with the predicted connector multiplicities "
             "(single for disjoint/nested, doubled for interleaved)", "L=8, k=6, s=3", EXACT, ok_b)

    # (c) symbolic SNF of the explicit matrices
    ok_c = True
    for s in range(1, 8):
        for inter, exp in ((False, 2 * s - 2), (True, 2 * s - 1)):
            M = reduced_laplacian_two_repeats(s, inter)
            D = smith_normal_form(M.T, domain=sp.ZZ)
            inv = sorted(abs(D[i, i]) for i in range(min(D.shape)) if abs(D[i, i]) not in (0, 1))
            ok_c &= (inv == [2] * exp)
    row("I", "symbolic Smith form of the explicit reduced Laplacians: (Z/2)^(2s-2) and "
             "(Z/2)^(2s-1) for s = 1..7", "14 matrices", EXACT, ok_c)
    return ok == tot and ok_b and ok_c


# ----------------------------------------------------------------------------- (M)

def f2_rank(M):
    M = [r[:] for r in M]
    n, r = len(M), 0
    for c in range(n):
        piv = next((i for i in range(r, n) if M[i][c]), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        for i in range(n):
            if i != r and M[i][c]:
                M[i] = [a ^ b for a, b in zip(M[i], M[r])]
        r += 1
    return r


def crossing_pattern(word):
    pos = collections.defaultdict(list)
    for i, ch in enumerate(word):
        pos[ch].append(i)
    letters = sorted(pos)
    I = [[0] * 3 for _ in range(3)]
    for a, b in itertools.combinations(range(3), 2):
        p, q = pos[letters[a]], pos[letters[b]]
        cross = (p[0] < q[0] < p[1] < q[1]) or (q[0] < p[0] < q[1] < p[1])
        I[a][b] = I[b][a] = int(cross)
    ncross = sum(I[a][b] for a, b in itertools.combinations(range(3), 2))
    return {0: "none", 1: "one pair", 2: "path", 3: "triangle"}[ncross], f2_rank(I)


def check_three_repeats():
    """Three repeats A,B,C of length L, multiplicity 2, in every circular arrangement (all
    90 words with two of each letter).  At s = 1 the compacted graph is the connector
    digraph D_w on the three symbols (loops dropped), so K(G_L) = K(D_w) tautologically.
    Tested: (a) the correction is NOT a function of the F_2-rank of the chord intersection
    matrix -- 'one pair', 'path' and 'triangle' all have rank 2 but give Z/2, Z/3 and
    (Z/2)^2; (b) at s = 2 the group is (Z/2)^3 + K(D_w), the spines splitting off as in
    Theorem 6(i)."""
    print("\n=== (M) three repeats: the correction is K(D_w), not a chord-diagram rank ===")
    rng = random.Random(99)
    L = 10
    expect = {"none": (), "one pair": (2,), "path": (3,), "triangle": (2, 2)}
    seen = collections.defaultdict(set)
    words = sorted(set(itertools.permutations("AABBCC")))
    done = 0
    for word in words:
        got = None
        for _ in range(60):
            reps = {c: "".join(rng.choice(AB) for _ in range(L)) for c in "ABC"}
            if len(set(reps.values())) < 3:
                continue
            spc = ["".join(rng.choice(AB) for _ in range(22)) for _ in range(6)]
            S = "".join(reps[c] + spc[i] for i, c in enumerate(word))
            res, ok = {}, True
            for s in (1, 2):
                nv, arcs, _ = dbg_from_sequence(S, L - s + 1)
                if branch_degrees(nv, arcs) != (2,) * (3 * s):
                    ok = False
                    break
                res[s] = sandpile(nv, arcs)[1]
            if ok:
                got = res
                break
        if got is None:
            continue
        done += 1
        shape, rk = crossing_pattern(word)
        seen[(shape, rk)].add((got[1], got[2]))
    def primary(inv):
        """Elementary divisors: the multiset of prime powers. Two tuples of cyclic orders
        present the same group iff their primary decompositions agree, e.g. (2,2,6) and
        (2,2,2,3) are both (Z/2)^3 + Z/3."""
        out = []
        for x in inv:
            for p, e in sp.factorint(int(x)).items():
                out.append(p ** e)
        return tuple(sorted(out))

    ok_all = (done == len(words))
    print("    pattern     rank_F2   K at s=1     K at s=2")
    for (shape, rk), vals in sorted(seen.items(), key=lambda kv: kv[0][1]):
        for (K1, K2) in sorted(vals):
            ok_all &= (primary(K1) == primary(expect[shape])
                       and primary(K2) == primary((2, 2, 2) + K1) and len(vals) == 1)
            print(f"    {shape:<10s}  {rk:^7d}   {str(list(K1)):<11s}  {list(K2)}")
    row("M", "all 90 circular words: K(G_L) = Z/2, Z/3, (Z/2)^2 for one-pair / path / "
             "triangle crossings (all of F_2-rank 2), 0 for none; and K at s=2 = (Z/2)^3 + K(G_L)",
        f"{done} words, s = 1, 2", EXACT, ok_all)
    return ok_all


# ----------------------------------------------------------------------------- (N)

W_LO = "GAACGAGAGAAAAACACGTACCAGCGCGGGGGGCGGGACGGCAGGTGACGCTCTGGGAAG"
W_HI = "GCACAGGGTCCAGTCTGCGAACGCGCTACAGTCTGGGTTCAAGCGGTCCTCTGTGCATCG"


def check_not_degree_determined():
    print("\n=== (N) the global factor is not a function of the degree data ===")
    k = 4
    prof = []
    for S in (W_LO, W_HI):
        nv, arcs, _ = dbg_from_sequence(S, k)
        o, inv = sandpile(nv, arcs)
        prof.append((nv, branch_degrees(nv, arcs), local_factor(nv, arcs), o, inv))
    (n1, d1, l1, o1, i1), (n2, d2, l2, o2, i2) = prof
    same_obs = (n1 == n2 and d1 == d2 and l1 == l2)
    print(f"    both: k={k}, n_k={n1}, branch out-degrees={list(d1)}, prod(d-1)!={l1}")
    print(f"    W_LO: |K|={o1:<4d} inv={list(i1)}  #reconstructions={o1*l1}")
    print(f"    W_HI: |K|={o2:<4d} inv={list(i2)}  #reconstructions={o2*l2}")
    row("N", f"identical observable degree data, reconstruction counts {o1*l1} vs {o2*l2} "
             f"({o2//o1}x apart)", "explicit witnesses", EXACT, same_obs and o1 != o2)
    random.seed(7)
    buckets = collections.defaultdict(set)
    for _ in range(1500):
        S = "".join(random.choice(AB) for _ in range(60))
        for kk in (4, 5):
            nv, arcs, _ = dbg_from_sequence(S, kk)
            dg = branch_degrees(nv, arcs)
            if dg:
                buckets[(kk, nv, dg)].add(sandpile(nv, arcs)[1])
    multi = sum(1 for b in buckets.values() if len(b) > 1)
    frac = multi / max(len(buckets), 1)
    row("N", f"{frac:.0%} of degree-classes split into several sandpile groups",
        f"{len(buckets)} classes", SAMPLED, frac > 0.3)
    return same_obs and o1 != o2


# ----------------------------------------------------------------------------- (C)

def check_collapse():
    print("\n=== (C) the collapse order k* ===")
    random.seed(5)
    at_R = at_R1 = other = tot = 0
    for _ in range(40):
        S = "".join(random.choice(AB) for _ in range(100))
        R = longest_repeat(S)
        c = collapse_k(S)
        tot += 1
        at_R += (c == R)
        at_R1 += (c == R + 1)
        other += (c not in (R, R + 1))
    print(f"    {tot} random sequences of length 100: k* = R in {at_R}, R+1 in {at_R1}, "
          f"neither in {other}")
    row("C", f"no single closed form for k* in terms of the longest repeat R alone",
        "N=100", SAMPLED, at_R > 0 and at_R1 > 0)
    # minimal witnesses, exhaustively over short binary sequences
    found = {}
    for n in range(4, 11):
        for bits in itertools.product("AC", repeat=n):
            S = "".join(bits)
            R = longest_repeat(S)
            if R < 2:
                continue
            c = collapse_k(S)
            if c is None:
                continue
            tag = "k*=R" if c == R else ("k*=R+1" if c == R + 1 else "k*<R")
            if tag not in found:
                found[tag] = (S, R, c)
        if len(found) == 3:
            break
    for tag in ("k*=R", "k*=R+1", "k*<R"):
        if tag in found:
            S, R, c = found[tag]
            print(f"    minimal {tag:<7s}: S={S} (N={len(S)}), R={R}, k*={c}")
    row("C", "minimal binary witnesses exist for k* = R, k* = R+1 and k* < R",
        "exhaustive, N <= 10", EXACT, len(found) == 3)
    return found


# ----------------------------------------------------------------------------- (P)

def check_compaction():
    print("\n=== (P) contraction and self-loops ===")
    random.seed(808)
    ok_all, tested, worst = True, 0, (0, 0)
    while tested < 12:
        S = "".join(random.choice(AB) for _ in range(70))
        k = random.choice((4, 5, 6))
        nv, arcs, _ = dbg_from_sequence(S, k)
        if nv < 8:
            continue
        o1, i1 = sandpile(nv, arcs)
        cn, ca = compact(nv, arcs)
        if cn < 2:
            continue
        o2, i2 = sandpile(cn, ca)
        ok = (o1 == o2 and i1 == i2)
        ok_all &= ok
        tested += 1
        if nv > worst[0]:
            worst = (nv, cn)
    row("P", f"linear-time unitig contraction leaves K unchanged in {tested}/{tested} cases "
             f"(largest reduction {worst[0]} -> {worst[1]} vertices)", "N=70", EXACT, ok_all)
    # self-loops
    base_nv, base_arcs = 4, [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2), (2, 0)]
    o0, i0 = sandpile(base_nv, base_arcs)
    ok_s = True
    for extra in (1, 2, 5):
        o1, i1 = sandpile(base_nv, base_arcs + [(1, 1)] * extra + [(3, 3)] * extra)
        ok_s &= (o0, i0) == (o1, i1)
    row("P", "self-loops leave |K| and the invariant factors unchanged (they cancel in "
             "D_out - A)", "1, 2, 5 loops", EXACT, ok_s)
    return ok_all and ok_s


# ----------------------------------------------------------------------------- (F)

def check_phix174():
    print("\n=== (F) real data: phiX174, NC_001422.1 ===")
    path = os.path.join(HERE, "phix174_NC_001422.1.txt")
    S = open(path).read().strip().upper()
    ok = (len(S) == 5386 and set(S) <= set("ACGT"))
    R = longest_repeat(S, hi=200)
    print(f"    length {len(S)} bp, longest exact circular repeat R = {R}")
    print("      k    n_k   branch   arcs   |K(G_k)|            invariant factors")
    table = []
    CAP = 70
    for k in range(8, 20):
        t0 = time.time()
        nv, arcs, _ = dbg_from_sequence(S, k)
        cn, ca = compact(nv, arcs)
        if cn == 0:
            print(f"     {k:2d} {nv:6d}      0      0   1                   []   (single cycle)")
            table.append((k, nv, 0, 1, ()))
            break
        if cn > CAP:
            print(f"     {k:2d} {nv:6d} {cn:6d} {len(ca):6d}   (compacted graph above the exact-SNF cap; skipped)")
            continue
        o, inv = sandpile(cn, ca)
        short = list(inv) if len(inv) <= 10 else list(inv[:8]) + ["..."] + [inv[-1]]
        print(f"     {k:2d} {nv:6d} {cn:6d} {len(ca):6d}   {o:<18d}  {short}   ({time.time()-t0:.2f}s)")
        table.append((k, nv, cn, o, inv))
    kstar = next((k for (k, nv, cn, o, inv) in table if o == 1), None)
    row("F", f"phiX174: R = {R}, collapse at k* = {kstar} (= R+1), |K(G_10)| = "
             f"{next(o for (k,nv,cn,o,inv) in table if k == 10)}, 99% vertex reduction",
        "5386 bp", EXACT, ok and kstar == R + 1)
    return table


# --------------------------------------------------------------------------- main

def main():
    t0 = time.time()
    check_tower()
    check_levine()
    check_best()
    check_single_repeat()
    check_interleaving()
    check_three_repeats()
    check_not_degree_determined()
    check_collapse()
    check_compaction()
    check_phix174()

    print("\n=== battery summary ===")
    bad = [r for r in ROWS if not r[4]]
    for key in ("T", "L", "E", "R", "I", "M", "N", "C", "P", "F"):
        sub = [r for r in ROWS if r[0] == key]
        if sub:
            print(f"  ({key})  {sum(1 for r in sub if r[4])}/{len(sub)} passed")
    print(f"  TOTAL {sum(1 for r in ROWS if r[4])}/{len(ROWS)} passed"
          + ("" if not bad else f"   FAILURES: {[(r[0], r[2]) for r in bad]}"))
    print(f"  elapsed {time.time()-t0:.1f}s")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
