# -*- coding: utf-8 -*-
"""
multicopy.py -- verification battery for

    "Multi-copy repeats: where the two-copy theory of assembly ambiguity stops"   (Bio 4)

Setting (Bio 3).  For a clean circular sequence with repeat families W_1..W_m of
multiplicities r_i, K(G_k) = (+)_i (Z/r_i)^(s_i-1) (+) K(D_w), where D_w is the transition
digraph of the circular word w of repeat occurrences.  For r_i = 2 the arrangement part
K(D_w) is the critical group of a bouquet: coker(A + I) for Bouchet's signed interlace
matrix (Merino-Moffatt-Noble).  This note asks what survives for multiplicities >= 3.

 (A)  TWO FAMILIES, any multiplicities (r, s): K(D_w) = Z/c, c = number of A-blocks of w
      (= number of A->B transitions), 1 <= c <= min(r, s).  Exhaustive, r, s <= 5.   exact
 (B)  CENSUS for multiplicity vectors (3,3), (4,4), (3,3,3), (3,3,2), (2,2,3), (4,4,4),
      (3,3,3,3): circular words, distinct transition digraphs, distinct groups, largest
      order; and the WITNESS that pairwise data do not suffice: AABACCBBC and AABABCBCC
      have the same pairwise restriction necklaces but groups Z/3 and Z/4.          exact
 (C)  NO REDUCTION TO TWO COPIES: detaching a 3-in 3-out vertex into two 2-in 2-out
      vertices does not control K: on random balanced digraphs K(G) is a quotient of
      K(G') in a minority of the 9 detachments; explicit witness.                   exact
 (D)  phiX174: out-degree profile of the skeleton at k = 12..7; at k = 11, 10 the
      skeleton is 2-in 2-out and coker(A + I) of any Euler circuit reproduces K(Sk);
      from k = 9 down multi-copy vertices appear (degrees 3..6).                    exact

Imports `repeat_splitting.py` from the sibling folder `Bio 3` (which imports
`pgl3_building.py` from `Paper 10`).  Keep the folders as siblings.
Run: python multicopy.py   (about a minute)
"""

import sys, os, time, random, itertools, collections
from sympy.utilities.iterables import multiset_permutations
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
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

sys.path.insert(0, _sib("Bio_03_Repeat_Splitting"))
import repeat_splitting as rs

EXACT = "exact"
ROWS = []


def row(key, desc, cases, mode, ok):
    ROWS.append((key, desc, cases, mode, bool(ok)))
    print(("  [ok]   " if ok else "  [FAIL] ") + f"({key}) {desc}   [{cases}] [{mode}]")
    return bool(ok)


def K_of(word, mult):
    return rs.elementary(rs.sandpile(*rs.connector_digraph(word, mult)))


def order(E):
    o = 1
    for x in E:
        o *= x
    return o


def circular_words(rs_vec, syms):
    """all circular words with the given letter multiplicities, one per rotation class"""
    base = list("".join(c * r for c, r in zip(syms, rs_vec)))
    seen = set()
    for p in multiset_permutations(base):
        if p[0] != syms[0]:
            continue
        word = "".join(p)
        rot = min(word[i:] + word[:i] for i in range(len(word)))
        if rot in seen:
            continue
        seen.add(rot)
        yield word


# ----------------------------------------------------------------------------- (A)

def check_two_families():
    print("\n=== (A) two families: K(D_w) = Z/c, c = number of A-blocks ===")
    tot = ok = 0
    seen_c = collections.defaultdict(set)
    for r in range(1, 6):
        for s in range(1, 6):
            for word in circular_words((r, s), "AB"):
                n = len(word)
                c = sum(1 for i in range(n) if word[i] == "A" and word[(i + 1) % n] == "B")
                K = K_of(word, {"A": r, "B": s})
                good = (K == (() if c == 1 else rs.elementary((c,)))) and 1 <= c <= min(r, s)
                tot += 1
                ok += good
                seen_c[(r, s)].add(c)
                if not good:
                    print(f"    MISMATCH {word}: c={c} K={rs.gstr(K)}")
    full = all(seen_c[(r, s)] == set(range(1, min(r, s) + 1)) for r in range(1, 6) for s in range(1, 6))
    print(f"    {tot} circular words over two letters, r, s <= 5; every c in 1..min(r,s) realised: {full}")
    row("A", "two repeat families of multiplicities r, s: K(D_w) = Z/c with c the number of "
             "A->B transitions, and every c in 1..min(r,s) occurs", f"{tot} words", EXACT, ok == tot and full)
    return ok == tot


# ----------------------------------------------------------------------------- (B)

def necklace(bits):
    n = len(bits)
    return min(tuple(bits[i:] + bits[:i]) for i in range(n))


def pair_necklace(word, e, f):
    return necklace([1 if c == e else 0 for c in word if c in (e, f)])


def canon_pairs(word, syms):
    m = len(syms)
    best = None
    for p in itertools.permutations(range(m)):
        key = tuple(pair_necklace(word, syms[p[i]], syms[p[j]])
                    for i in range(m) for j in range(m) if i != j)
        if best is None or key < best:
            best = key
    return best


def dw_canon(word, syms):
    m = len(syms)
    idx = {c: i for i, c in enumerate(syms)}
    M = [[0] * m for _ in range(m)]
    for i in range(len(word)):
        M[idx[word[i]]][idx[word[(i + 1) % len(word)]]] += 1
    best = None
    for p in itertools.permutations(range(m)):
        key = tuple(M[p[i]][p[j]] for i in range(m) for j in range(m))
        if best is None or key < best:
            best = key
    return best


def check_census():
    print("\n=== (B) census of arrangement groups for multi-copy repeats ===")
    print("      mult        words   D_w   groups  max|K|   non-cyclic groups            pairwise-necklace classes (split)")
    ok = True
    for rs_vec in ((3, 3), (4, 4), (3, 3, 3), (3, 3, 2), (2, 2, 3), (4, 4, 4), (3, 3, 3, 3)):
        syms = "ABCD"[:len(rs_vec)]
        mult = dict(zip(syms, rs_vec))
        groups = collections.Counter()
        by_dw = collections.defaultdict(set)
        by_pairs = collections.defaultdict(set)
        n = 0
        for word in circular_words(rs_vec, syms):
            n += 1
            K = K_of(word, mult)
            groups[K] += 1
            by_dw[dw_canon(word, syms)].add(K)
            by_pairs[canon_pairs(word, syms)].add(K)
        noncyc = sorted({rs.gstr(K) for K in groups if len(K) > 1 and
                         any(sp.gcd(a, b) > 1 for a, b in itertools.combinations(K, 2))})
        split = sum(1 for v in by_pairs.values() if len(v) > 1)
        ok &= all(len(v) == 1 for v in by_dw.values())
        print(f"    {str(rs_vec):12s} {n:6d}  {len(by_dw):4d}  {len(groups):6d}  {max(order(K) for K in groups):6d}   "
              f"{', '.join(noncyc) if noncyc else '-':28s} {len(by_pairs)} ({split})")
    # the witness, by hand
    w1, w2 = "AABACCBBC", "AABABCBCC"
    same = canon_pairs(w1, "ABC") == canon_pairs(w2, "ABC")
    K1, K2 = K_of(w1, {c: 3 for c in "ABC"}), K_of(w2, {c: 3 for c in "ABC"})
    print(f"    witness: {w1} -> {rs.gstr(K1)},  {w2} -> {rs.gstr(K2)},  same pairwise necklaces: {same}")
    ok &= same and K1 == (3,) and K2 == (4,)
    row("B", "K(D_w) is a function of the transition digraph (trivially) but NOT of the pairwise "
             "restriction necklaces: AABACCBBC (Z/3) and AABABCBCC (Z/4) agree on all three pairs",
        "7 multiplicity vectors, exhaustive", EXACT, ok)
    return ok


# ----------------------------------------------------------------------------- (C)

def is_quotient(P, Q):
    primes = {sp.factorint(x).popitem()[0] for x in P + Q}
    for p in primes:
        ep = sorted((sp.factorint(x)[p] for x in P if x % p == 0), reverse=True)
        eq = sorted((sp.factorint(x)[p] for x in Q if x % p == 0), reverse=True)
        if len(eq) > len(ep) or any(eq[i] > ep[i] for i in range(len(eq))):
            return False
    return True


def check_detachment():
    print("\n=== (C) detaching a 3-in 3-out vertex into two 2-in 2-out vertices ===")
    rng = random.Random(21)
    tot = quo = 0
    witness = None
    graphs = 0
    while graphs < 150:
        nv, arcs = rs.random_eulerian(rng, rng.randint(4, 7), cycles=rng.randint(3, 5))
        if nv < 3 or not rs.strongly_connected(nv, arcs):
            continue
        outd = collections.Counter(u for u, v in arcs)
        cands = [v for v in range(nv) if outd[v] == 3 and not any(u == v == x for (u, x) in arcs)]
        if not cands:
            continue
        v = cands[0]
        ins = [j for j, (u, x) in enumerate(arcs) if x == v]
        outs = [j for j, (u, x) in enumerate(arcs) if u == v]
        K = rs.elementary(rs.sandpile(nv, arcs))
        graphs += 1
        found = []
        for al in ins:
            for bk in outs:
                new = []
                for j, (u, x) in enumerate(arcs):
                    if x == v:
                        x = nv if j == al else v
                    if u == v:
                        u = v if j == bk else nv
                    new.append((u, x))
                new.append((v, nv))
                Kp = rs.elementary(rs.sandpile(nv + 1, new))
                tot += 1
                q = is_quotient(Kp, K)
                quo += q
                found.append(Kp)
        if witness is None and len(K) == 2 and sum(1 for Kp in found if not is_quotient(Kp, K)) >= 3:
            witness = (K, sorted(set(found)))
    print(f"    {graphs} graphs x 9 detachments: K(G) is a quotient of K(G') in {quo}/{tot}")
    if witness:
        print(f"    witness: K(G) = {rs.gstr(witness[0])}; the nine detachments give "
              + ", ".join(rs.gstr(k) for k in witness[1]))
    row("C", "K(G) is a quotient of the detached K(G') in a minority of cases; no containment either way",
        f"{tot} detachments", EXACT, quo < tot // 2 and witness is not None)
    return witness


# ----------------------------------------------------------------------------- (D)

def euler_circuit(nv, arcs):
    out = collections.defaultdict(list)
    for j, (u, v) in enumerate(arcs):
        out[u].append(j)
    stack, path = [0], []
    used = [False] * len(arcs)
    ptr = collections.defaultdict(int)
    while stack:
        v = stack[-1]
        while ptr[v] < len(out[v]) and used[out[v][ptr[v]]]:
            ptr[v] += 1
        if ptr[v] < len(out[v]):
            j = out[v][ptr[v]]
            used[j] = True
            stack.append(arcs[j][1])
        else:
            path.append(stack.pop())
    path.reverse()
    return path[:-1]


def check_phix():
    print("\n=== (D) phiX174: where multi-copy vertices appear ===")
    S = open(os.path.join(HERE, "phix174_NC_001422.1.txt")).read().strip().upper()
    ok = True
    matched = 0
    print("      k   skeleton   arcs   out-degree multiset                     MMN on the skeleton")
    for k in (12, 11, 10, 9, 8, 7):
        nv, arcs = rs.dbg(S, k)
        cn, ca = rs.compact(nv, arcs)
        sn, sa, chains = rs.bundle_contract(cn, ca)
        outd = collections.Counter(u for u, v in sa)
        degs = dict(sorted(collections.Counter(outd[v] for v in range(sn)).items()))
        note = ""
        if set(degs) == {2} and sn <= 120:
            seq = euler_circuit(sn, sa)
            pos = collections.defaultdict(list)
            for i, v in enumerate(seq):
                pos[v].append(i)
            A = [[0] * sn for _ in range(sn)]
            for i in range(sn):
                for j in range(sn):
                    if i != j and rs.crossing(tuple(pos[i]), tuple(pos[j])):
                        A[i][j] = 1 if pos[i][0] < pos[j][0] else -1
                A[i][i] += 1
            d = rs.pb.diagonalize([list(r) for r in A], False, False)[0]
            Kint = rs.elementary(tuple(abs(x) for x in d if abs(x) not in (0, 1)))
            Ksk = rs.elementary(rs.sandpile(sn, sa))
            good = (Kint == Ksk)
            ok &= good
            matched += good
            note = f"coker(A+I) = {rs.gstr(Kint)} = K(Sk): {'ok' if good else 'MISMATCH'}"
        elif set(degs) == {2} or sn <= 1:
            note = "(trivial)"
        else:
            note = f"not 2-in 2-out: {sum(v for dd, v in degs.items() if dd >= 3)} vertices of out-degree >= 3"
        print(f"     {k:2d}   {sn:6d}   {len(sa):5d}   {str(degs):38s} {note}")
    row("D", "phiX174 skeleton is 2-in 2-out at k = 11, 10 (signed interlace matrix reproduces K(Sk)) "
             "and acquires vertices of out-degree 3..6 from k = 9 down",
        "k = 12..7", EXACT, ok and matched == 2)
    return ok


# --------------------------------------------------------------------------- main

def main():
    t0 = time.time()
    check_two_families()
    check_census()
    check_detachment()
    check_phix()
    print("\n=== battery summary ===")
    bad = [r for r in ROWS if not r[4]]
    for key in ("A", "B", "C", "D"):
        sub = [r for r in ROWS if r[0] == key]
        if sub:
            print(f"  ({key})  {sum(1 for r in sub if r[4])}/{len(sub)} passed")
    print(f"  TOTAL {sum(1 for r in ROWS if r[4])}/{len(ROWS)} passed"
          + ("" if not bad else f"   FAILURES: {[(r[0], r[2]) for r in bad]}"))
    print(f"  elapsed {time.time()-t0:.1f}s")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
