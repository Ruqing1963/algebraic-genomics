# -*- coding: utf-8 -*-
"""
repeat_splitting.py -- verification battery for

    "The sandpile group of a de Bruijn graph splits along its repeats"   (Bio 3)

THEOREM (bundle-chain contraction).  Let G be a balanced multidigraph and v_1 -> ... -> v_s
a chain in which every v_t (t < s) has out-degree r with ALL r arcs going to v_{t+1}, and
v_s has out-degree r.  Contract the chain to one vertex (in-arcs of v_1, out-arcs of v_s)
to get G'.  Then   K(G)  =  (Z/r)^(s-1)  (+)  K(G').
Proof: generators x_v of coker(Delta^T); put z = x_{v_1}, y_t = x_{v_t} - x_{v_{t+1}}.
Rows of v_t (t<s) become r*y_t.  Arcs INTO the chain land on v_1 = z.  Arcs OUT of the
chain leave v_s = z - sum y_u with coefficient r, and r*y_u = 0, so those rows become
the rows of G'.  Unimodular.  Iterating to a chain-free SKELETON gives

    K(G_k(S))  =  (+)_chains (Z/r_c)^(len_c - 1)  (+)  K(skeleton),

and for a sequence whose repeats are "clean" the skeleton is the connector digraph D_w of
the circular word of repeat occurrences (multiplicities kept).

 (A)  the clean case, randomised: m in {2,3,4} repeats, multiplicities r_i in {2,3,4},
      lengths L_i, several k:  K(G_k) = (+)_i (Z/r_i)^(s_i-1) (+) K(D_w).          exact
 (B)  the lemma on arbitrary balanced digraphs: random Eulerian digraph H, plant a bundle
      chain of multiplicity r and length s at a vertex; K(G) = (Z/r)^(s-1) (+) K(H). exact
 (C)  the special cases already published: m=1 (Thm 5 of Bio 1), m=2 (Thm 6), m=3
      (Prop 8) are recovered by the formula.                                        exact
 (F)  phiX174: bundle-contract the compacted graph at k = 9,10,11,12; list the chains
      (r, length), the skeleton size and K(skeleton); verify the decomposition against
      K(G_k) computed directly where feasible, and EXTEND the spectrum to k = 9.     exact

Groups are compared by elementary divisors.  Exact integer arithmetic throughout.
Run: python repeat_splitting.py   (a few minutes)
"""

import sys, os, time, random, itertools, collections

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

HERE = os.path.dirname(os.path.abspath(__file__))
EXACT = "exact"
ROWS = []
AB = "ACGT"


def row(key, desc, cases, mode, ok):
    ROWS.append((key, desc, cases, mode, bool(ok)))
    print(("  [ok]   " if ok else "  [FAIL] ") + f"({key}) {desc}   [{cases}] [{mode}]")
    return bool(ok)


# --------------------------------------------------------------------------- basics

def dbg(S, k):
    N = len(S)
    SS = S + S[:k]
    verts = {}
    for i in range(N):
        w = SS[i:i + k]
        if w not in verts:
            verts[w] = len(verts)
    return len(verts), [(verts[SS[i:i + k]], verts[SS[i + 1:i + 1 + k]]) for i in range(N)]


def sandpile(nv, arcs, root=0):
    """invariant factors (sorted) of coker of the transposed reduced out-Laplacian"""
    if nv <= 1:
        return ()
    L = [[0] * nv for _ in range(nv)]
    for (u, v) in arcs:
        L[u][u] += 1
        L[u][v] -= 1
    M = [[L[i][j] for j in range(nv) if j != root] for i in range(nv) if i != root]
    d = pb.diagonalize([list(r) for r in zip(*M)], False, False)[0]
    return tuple(sorted(abs(x) for x in d if abs(x) not in (0, 1)))


def elementary(inv):
    """multiset of prime powers; the canonical form for comparing finite abelian groups"""
    out = []
    for x in inv:
        for p, e in sp.factorint(int(x)).items():
            out.append(p ** e)
    return tuple(sorted(out))


def order(inv):
    o = 1
    for x in inv:
        o *= x
    return o


def gstr(inv):
    if not inv:
        return "0"
    c = collections.Counter(elementary(inv))
    return " + ".join(f"(Z/{p})^{e}" if e > 1 else f"Z/{p}" for p, e in sorted(c.items()))


def compact(nv, arcs):
    """unitig (r = 1) contraction, linear time"""
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
            x = f
            while x not in bs:
                x = ol[x][0]
            na.append((rl[b], rl[x]))
    return len(br), na


def bundle_contract(nv, arcs):
    """Contract every maximal bundle chain.  Returns (skeleton nv, skeleton arcs,
    list of (r, length) of the chains contracted).  A vertex v is a bundle head onto w if
    all of its out-arcs go to the single vertex w != v (multiplicity r = outdeg v) and
    outdeg w == r (balanced => all in-arcs of w come from v)."""
    arcs = list(arcs)
    chains = []
    changed = True
    while changed:
        changed = False
        outd = collections.Counter(u for u, v in arcs)
        targets = collections.defaultdict(set)
        for (u, v) in arcs:
            targets[u].add(v)
        # find a vertex whose arcs all go to one other vertex with the same out-degree
        for v in range(nv):
            if outd[v] == 0 or len(targets[v]) != 1:
                continue
            w = next(iter(targets[v]))
            if w == v or outd[w] != outd[v]:
                continue
            r = outd[v]
            # extend the chain forward as far as possible
            chain = [v, w]
            while True:
                x = chain[-1]
                if len(targets[x]) == 1:
                    y = next(iter(targets[x]))
                    if y not in chain and outd[y] == r:
                        chain.append(y)
                        continue
                break
            # extend backward
            preds = collections.defaultdict(set)
            for (u, x) in arcs:
                preds[x].add(u)
            while True:
                x = chain[0]
                if len(preds[x]) == 1:
                    u = next(iter(preds[x]))
                    if u not in chain and outd[u] == r and targets[u] == {x}:
                        chain.insert(0, u)
                        continue
                break
            s = len(chain)
            # contract: all chain vertices -> chain[0]; drop internal arcs
            head, tail = chain[0], chain[-1]
            cs = set(chain)
            new = []
            for (u, x) in arcs:
                if u in cs and x in cs:
                    continue                      # internal spine arc
                nu = head if u in cs else u
                nx = head if x in cs else x
                new.append((nu, nx))
            arcs = new
            chains.append((r, s))
            changed = True
            break
    # relabel surviving vertices
    alive = sorted({u for u, v in arcs} | {v for u, v in arcs})
    if not alive:
        return 1, [], chains
    rl = {v: i for i, v in enumerate(alive)}
    return len(alive), [(rl[u], rl[v]) for (u, v) in arcs], chains


def predicted(chains, skel_inv):
    parts = list(skel_inv)
    for r, s in chains:
        parts += [r] * (s - 1)
    return elementary(tuple(parts))


def connector_digraph(word, mult):
    """D_w: vertices = symbols, one arc per cyclic adjacency (loops kept; harmless)."""
    syms = sorted(set(word))
    idx = {c: i for i, c in enumerate(syms)}
    arcs = [(idx[word[i]], idx[word[(i + 1) % len(word)]]) for i in range(len(word))]
    return len(syms), arcs


# ----------------------------------------------------------------------------- (A)

def clean_instance(rng, m, rs, Ls, k, spacer=24):
    """random sequence with m repeats of lengths Ls and multiplicities rs; returns
    (S, word) or None if a chance repeat spoils the clean hypothesis at this k."""
    syms = "ABCDEFG"[:m]
    reps = {}
    for c, L in zip(syms, Ls):
        reps[c] = "".join(rng.choice(AB) for _ in range(L))
    if len(set(reps.values())) < m:
        return None
    word = list("".join(c * r for c, r in zip(syms, rs)))
    rng.shuffle(word)
    S = "".join(reps[c] + "".join(rng.choice(AB) for _ in range(spacer)) for c in word)
    nv, arcs = dbg(S, k)
    outd = collections.Counter(u for u, v in arcs)
    want = collections.Counter()
    for c, r, L in zip(syms, rs, Ls):
        s = L - k + 1
        if s >= 1:
            want[r] += s
    have = collections.Counter(d for d in outd.values() if d > 1)
    if have != want:
        return None
    return S, "".join(word)


def check_clean():
    print("\n=== (A) clean multi-repeat sequences: K(G_k) = (+)(Z/r_i)^(s_i-1) (+) K(D_w) ===")
    rng = random.Random(2026)
    tot = ok = 0
    shown = 0
    for trial in range(4000):
        m = rng.choice((2, 3, 3, 4))
        rs = [rng.choice((2, 2, 3, 4)) for _ in range(m)]
        Ls = [rng.choice((10, 11, 12, 13)) for _ in range(m)]
        k = rng.randint(max(9, min(Ls) - 3), min(Ls))
        inst = clean_instance(rng, m, rs, Ls, k)
        if inst is None:
            continue
        S, word = inst
        nv, arcs = dbg(S, k)
        K = sandpile(nv, arcs)
        dn, da = connector_digraph(word, rs)
        Kd = sandpile(dn, da)
        chains = [(r, L - k + 1) for r, L in zip(rs, Ls)]
        pred = predicted(chains, Kd)
        good = (elementary(K) == pred)
        tot += 1
        ok += good
        if not good:
            print(f"    MISMATCH word={word} r={rs} L={Ls} k={k}: K={gstr(K)} vs pred={gstr(tuple(pred))}")
        if shown < 6:
            shown += 1
            print(f"    word={word:<10s} r={rs} L={Ls} k={k}: K = {gstr(K):<28s} "
                  f"pred (+)(Z/r)^(s-1) + K(D_w)={gstr(Kd):<12s} {'ok' if good else 'MISMATCH'}")
        if tot >= 120:
            break
    row("A", "K(G_k) = (+)(Z/r_i)^(s_i-1) (+) K(D_w) on every clean instance, m in 2..4, "
             "r_i in 2..4", f"{tot} instances", EXACT, ok == tot and tot >= 40)
    return ok == tot


# ----------------------------------------------------------------------------- (B)

def random_eulerian(rng, n, cycles=4, maxlen=6):
    """a random balanced multidigraph: union of random closed walks (may have loops/parallels)"""
    arcs = []
    for _ in range(cycles):
        L = rng.randint(2, maxlen)
        walk = [rng.randrange(n) for _ in range(L)]
        for i in range(L):
            arcs.append((walk[i], walk[(i + 1) % L]))
    # keep only the vertices used, relabel
    used = sorted({u for u, v in arcs} | {v for u, v in arcs})
    rl = {v: i for i, v in enumerate(used)}
    return len(used), [(rl[u], rl[v]) for (u, v) in arcs]


def strongly_connected(nv, arcs):
    adj = collections.defaultdict(list)
    radj = collections.defaultdict(list)
    for (u, v) in arcs:
        adj[u].append(v)
        radj[v].append(u)
    def reach(a, start):
        seen, st = {start}, [start]
        while st:
            x = st.pop()
            for y in a[x]:
                if y not in seen:
                    seen.add(y)
                    st.append(y)
        return seen
    return len(reach(adj, 0)) == nv and len(reach(radj, 0)) == nv


def plant_chain(nv, arcs, v, r, s):
    """replace vertex v (out-degree must be r) by a chain of s vertices with multiplicity r"""
    outs = [(u, x) for (u, x) in arcs if u == v]
    ins = [(u, x) for (u, x) in arcs if x == v and u != v]
    loops = [(u, x) for (u, x) in arcs if u == v and x == v]
    if len(outs) != r or loops:
        return None
    rest = [(u, x) for (u, x) in arcs if u != v and x != v]
    ids = [v] + list(range(nv, nv + s - 1))       # v_1 = v, v_2..v_s new
    new = rest + [(u, ids[0]) for (u, x) in ins]
    for t in range(s - 1):
        new += [(ids[t], ids[t + 1])] * r
    new += [(ids[-1], x) for (u, x) in outs]
    return nv + s - 1, new


def check_lemma():
    print("\n=== (B) the bundle-chain lemma on random balanced digraphs ===")
    rng = random.Random(7)
    tot = ok = 0
    while tot < 60:
        nv, arcs = random_eulerian(rng, rng.randint(4, 7))
        if nv < 3 or not strongly_connected(nv, arcs):
            continue
        outd = collections.Counter(u for u, v in arcs)
        cands = [v for v in range(nv) if outd[v] >= 2 and not any(u == v == x for (u, x) in arcs)]
        if not cands:
            continue
        v = rng.choice(cands)
        r = outd[v]
        s = rng.randint(2, 4)
        res = plant_chain(nv, arcs, v, r, s)
        if res is None:
            continue
        nv2, arcs2 = res
        K_big = sandpile(nv2, arcs2)
        K_small = sandpile(nv, arcs)
        pred = elementary(tuple(list(K_small) + [r] * (s - 1)))
        good = (elementary(K_big) == pred)
        tot += 1
        ok += good
        if not good:
            print(f"    MISMATCH: n={nv} r={r} s={s}: K(G)={gstr(K_big)} vs K(H)+(Z/r)^(s-1)={gstr(K_small)}+...")
    row("B", f"K(G) = (Z/r)^(s-1) (+) K(G') for a planted bundle chain, r in 2..6, s in 2..4",
        f"{tot} random balanced digraphs", EXACT, ok == tot)
    return ok == tot


# ----------------------------------------------------------------------------- (C)

def check_special_cases():
    print("\n=== (C) the published special cases are instances ===")
    ok = True
    # m = 1: D_w is one vertex with r loops -> K trivial -> (Z/r)^(s-1)
    for r, s in ((2, 3), (3, 4), (4, 2)):
        dn, da = connector_digraph("A" * r, [r])
        ok &= (predicted([(r, s)], sandpile(dn, da)) == elementary(tuple([r] * (s - 1))))
    # m = 2, r = 2: AABB -> 0 ; ABAB -> Z/2
    ok &= (sandpile(*connector_digraph("AABB", [2, 2])) == ())
    ok &= (elementary(sandpile(*connector_digraph("ABAB", [2, 2]))) == (2,))
    # m = 3: AABBCC -> 0 ; AABCBC -> Z/2 ; ABACBC -> Z/3 ; ABCABC -> (Z/2)^2
    exp = {"AABBCC": (), "AABCBC": (2,), "ABACBC": (3,), "ABCABC": (2, 2)}
    for w, e in exp.items():
        ok &= (elementary(sandpile(*connector_digraph(w, [2, 2, 2]))) == e)
    row("C", "m = 1 (Thm 5), m = 2 (Thm 6) and m = 3 (Prop 8) of the companion paper are "
             "the formula evaluated on D_w", "9 cases", EXACT, ok)
    return ok


# ----------------------------------------------------------------------------- (F)

def check_phix():
    print("\n=== (F) phiX174: bundle-contract the compacted graph and read the skeleton ===")
    S = open(os.path.join(HERE, "phix174_NC_001422.1.txt")).read().strip().upper()
    print("      k   compacted   chains (r,len)                       skeleton   K(skeleton)          "
          "K(G_k) from decomposition           direct")
    ok_all = True
    ext = {}
    for k in (12, 11, 10, 9, 8):
        nv, arcs = dbg(S, k)
        cn, ca = compact(nv, arcs)
        if cn == 0:
            print(f"     {k:2d}   single cycle")
            continue
        sn, sa, chains = bundle_contract(cn, ca)
        Ks = sandpile(sn, sa) if sn <= 120 else None
        if Ks is None:
            print(f"     {k:2d}   {cn:6d}     {len(chains)} chains{'':30s} {sn:6d}   (skeleton too large)")
            continue
        pred = predicted(chains, Ks)
        direct = sandpile(cn, ca) if cn <= 70 else None
        nch = len(chains)
        contrib = sum(s - 1 for r, s in chains)
        rset = sorted({r for r, s in chains})
        chain_str = f"{nch} chains, r in {rset}, sum(len-1) = {contrib}"
        d_str = gstr(direct) if direct is not None else "(not computed)"
        good = (direct is None) or (elementary(direct) == pred)
        ok_all &= good
        ext[k] = (cn, chains, sn, Ks, pred)
        print(f"     {k:2d}   {cn:6d}     {chain_str:<36s} {sn:6d}   {gstr(Ks):<20s} "
              f"{gstr(tuple(pred)):<34s} {d_str}  {'ok' if good else 'MISMATCH'}")
    row("F", "phiX174 k=10..12: K(G_k) computed directly equals the bundle-chain decomposition; "
             "k = 9 (146 branch vertices) becomes computable after contraction",
        "compacted graphs", EXACT, ok_all and 9 in ext)
    return ext


# ----------------------------------------------------------------------------- (D)

def chords(word):
    pos = collections.defaultdict(list)
    for i, c in enumerate(word):
        pos[c].append(i)
    return {c: tuple(p) for c, p in pos.items()}


def crossing(a, b):
    (a1, a2), (b1, b2) = a, b
    return (a1 < b1 < a2) != (a1 < b2 < a2)


def intersection_matrix(word):
    ch = chords(word)
    syms = sorted(ch)
    return [[1 if (i != j and crossing(ch[s], ch[t])) else 0 for j, t in enumerate(syms)]
            for i, s in enumerate(syms)]


def rank_mod2(M):
    rows = [sum(x << j for j, x in enumerate(r)) for r in M]
    rank = 0
    for bit in range(len(M)):
        piv = next((i for i in range(rank, len(rows)) if rows[i] >> bit & 1), None)
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        for i in range(len(rows)):
            if i != rank and rows[i] >> bit & 1:
                rows[i] ^= rows[rank]
        rank += 1
    return rank


def graph_canon(M):
    n = len(M)
    best = None
    for p in itertools.permutations(range(n)):
        key = tuple(M[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or key < best:
            best = key
    return best


def signed_interlace(word):
    """Bouchet's signed interlace matrix, read from the linear word: A[e][f] = +1 if the
    occurrences run e f e f, -1 if f e f e, 0 if not interlaced (skew-symmetric)."""
    ch = chords(word)
    syms = sorted(ch)
    A = [[0] * len(syms) for _ in syms]
    for i, s in enumerate(syms):
        for j, t in enumerate(syms):
            if i == j or not crossing(ch[s], ch[t]):
                continue
            A[i][j] = 1 if ch[s][0] < ch[t][0] else -1
    return A


def coker_A_plus_I(word):
    A = signed_interlace(word)
    for i in range(len(A)):
        A[i][i] += 1
    d = pb.diagonalize([list(r) for r in A], False, False)[0]
    return elementary(tuple(abs(x) for x in d if abs(x) not in (0, 1)))


def check_chord_invariants(mmax=5):
    """r_i = 2: on which invariant of the chord diagram does K(D_w) depend?  Exhaustive
    over all circular words with m letters each twice (rotation fixed by word[0] = 'A').
    Then: K(D_w) = coker(A + I) for the SIGNED interlace matrix (Merino-Moffatt-Noble,
    Thm 7.3, the bouquet case), and the unsigned interlace graph does NOT determine the
    group at m = 6 (their Remark 6.7)."""
    print("\n=== (D) two-copy repeats: K(D_w) against invariants of the chord diagram ===")
    ok = True
    # Merino-Moffatt-Noble's pair (1..6 -> a..f) and two further pairs found by search
    pairs = [("abcdefacbdfe", "abcfedacbefd"), ("abacdcefbdef", "abacdebcfdfe"),
             ("abacdefbcdfe", "abacdefbcedf")]
    for w1, w2 in pairs:
        H1, H2 = graph_canon(intersection_matrix(w1)), graph_canon(intersection_matrix(w2))
        K1 = elementary(sandpile(*connector_digraph(w1, {c: 2 for c in w1})))
        K2 = elementary(sandpile(*connector_digraph(w2, {c: 2 for c in w2})))
        same_H = (H1 == H2)
        print(f"    m=6 pair {w1} / {w2}: same interlace graph = {same_H}; "
              f"K = {gstr(K1)} / {gstr(K2)}  ({'DIFFERENT' if K1 != K2 else 'same'})")
        ok &= same_H and K1 != K2
    ok &= (elementary(sandpile(*connector_digraph("abcdefacbdfe", {c: 2 for c in "abcdef"})))
           == (2, 3, 3))                                                    # Z/3 + Z/6
    ok &= (elementary(sandpile(*connector_digraph("abcfedacbefd", {c: 2 for c in "abcdef"})))
           == (2, 9))                                                       # Z/18
    # signed matrix: K(D_w) = coker(A + I), exhaustive m <= 5 and random m = 6..9
    rng = random.Random(3)
    tot = good = 0
    for m in range(2, 10):
        letters = "abcdefghi"[:m]
        if m <= 5:
            words = ["".join(w) for w in set(itertools.permutations(letters * 2)) if w[0] == "a"]
        else:
            words = []
            for _ in range(300):
                w = list(letters * 2)
                rng.shuffle(w)
                words.append("".join(w))
        for word in words:
            K = elementary(sandpile(*connector_digraph(word, {c: 2 for c in letters})))
            tot += 1
            good += (K == coker_A_plus_I(word))
    print(f"    K(D_w) = coker(signed interlace matrix + I) on {good}/{tot} words "
          f"(all m <= 5 exhaustive, 300 random each for m = 6..9)")
    ok &= (good == tot)
    for m in range(3, mmax + 1):
        letters = "ABCDE"[:m]
        by_rank = collections.defaultdict(set)
        by_snf = collections.defaultdict(set)
        by_graph = collections.defaultdict(set)
        n = 0
        for w in set(itertools.permutations(letters * 2)):
            if w[0] != "A":
                continue
            word = "".join(w)
            n += 1
            M = intersection_matrix(word)
            K = elementary(sandpile(*connector_digraph(word, {c: 2 for c in letters})))
            r2 = rank_mod2(M)
            snf = tuple(sorted(abs(x) for x in pb.diagonalize([list(r) for r in M], False, False)[0]))
            by_rank[r2].add(K)
            by_snf[snf].add(K)
            by_graph[graph_canon(M)].add(K)
        rank_coll = {k: sorted(v) for k, v in by_rank.items() if len(v) > 1}
        snf_coll = {k: sorted(v) for k, v in by_snf.items() if len(v) > 1}
        graph_coll = {k: sorted(v) for k, v in by_graph.items() if len(v) > 1}
        print(f"    m={m}: {n} words, {len(by_graph)} intersection graphs")
        print(f"      groups sharing one F2-rank:       {rank_coll}")
        print(f"      groups sharing one integer SNF:   {snf_coll}")
        print(f"      intersection-graph classes with more than one group: {len(graph_coll)}")
        for k, v in graph_coll.items():
            print(f"        {k}: {v}")
        if m == 3:
            ok &= (2 in rank_coll and len(rank_coll[2]) == 3)     # Z/2, Z/3, (Z/2)^2 all rank 2
            ok &= any(len(v) > 1 for v in snf_coll.values())     # (1,1,0) carries Z/2 and Z/3
        ok &= not graph_coll
    row("D", f"r_i = 2: K(D_w) is NOT a function of the F2-rank nor of the integer Smith form of "
             f"the unsigned intersection matrix; it is constant on interlace-graph classes for "
             f"m <= {mmax} but NOT at m = 6 (Merino-Moffatt-Noble Rem. 6.7: Z/3+Z/6 vs Z/18); "
             f"it equals coker(A + I) for the signed interlace matrix (their Thm 7.3)",
        f"m = 3..{mmax} exhaustive; 3 pairs at m = 6; {tot} words for A + I", EXACT, ok)
    return ok


# ----------------------------------------------------------------------------- (E)

def maximal_chains(nv, arcs):
    """the bundle-chain steps v -> w (all out-arcs of v go to w != v, outdeg w = outdeg v)
    form a partial injection on vertices; its maximal paths are the maximal chains."""
    outd = collections.Counter(u for u, v in arcs)
    targets = collections.defaultdict(set)
    for (u, v) in arcs:
        targets[u].add(v)
    step = {}
    for v in range(nv):
        if outd[v] and len(targets[v]) == 1:
            w = next(iter(targets[v]))
            if w != v and outd[w] == outd[v]:
                step[v] = w
    pred = {w: v for v, w in step.items()}
    assert len(pred) == len(step), "a chain step target had two chain predecessors"
    chains = []
    seen = set()
    for v in range(nv):
        if v in seen or v in pred:
            continue                    # start only at heads (no predecessor)
        c = [v]
        while c[-1] in step and step[c[-1]] not in c:
            c.append(step[c[-1]])
        seen |= set(c)
        if len(c) >= 2:
            chains.append(tuple(c))
    # closed chains (a cycle of steps with no head) -- the whole component is one chain
    rest = [v for v in range(nv) if v not in seen and v in step]
    while rest:
        v = rest[0]
        c = [v]
        while step[c[-1]] not in c:
            c.append(step[c[-1]])
        seen |= set(c)
        chains.append(tuple(c))
        rest = [v for v in rest if v not in seen]
    return chains


def repeat_graph(word, mult):
    """the repeat graph of a clean circular sequence: vertices h_i, t_i per repeat family,
    r_i parallel arcs h_i -> t_i, and one arc t_i -> h_j per cyclic adjacency i j in w."""
    syms = sorted(set(word))
    h = {c: 2 * i for i, c in enumerate(syms)}
    t = {c: 2 * i + 1 for i, c in enumerate(syms)}
    arcs = []
    for c in syms:
        arcs += [(h[c], t[c])] * mult[c]
    arcs += [(t[word[i]], h[word[(i + 1) % len(word)]]) for i in range(len(word))]
    return 2 * len(syms), arcs


def check_confluence():
    print("\n=== (E) the skeleton is canonical; the repeat graph is the s_i = 2 case ===")
    rng = random.Random(11)
    tot = ok = 0
    while tot < 100:
        nv, arcs = random_eulerian(rng, rng.randint(4, 8))
        if nv < 3 or not strongly_connected(nv, arcs):
            continue
        for _ in range(rng.randint(1, 3)):
            outd = collections.Counter(u for u, v in arcs)
            cands = [v for v in range(nv) if outd[v] >= 2 and not any(u == v == x for (u, x) in arcs)]
            if not cands:
                break
            v = rng.choice(cands)
            res = plant_chain(nv, arcs, v, outd[v], rng.randint(2, 3))
            if res is not None:
                nv, arcs = res
        ch = maximal_chains(nv, arcs)
        verts = [v for c in ch for v in c]
        disjoint = len(verts) == len(set(verts))
        # random relabelling: the partition into maximal chains, the (r, len) multiset and
        # the skeleton's group must not change
        perm = list(range(nv))
        rng.shuffle(perm)
        arcs2 = [(perm[u], perm[v]) for (u, v) in arcs]
        ch2 = maximal_chains(nv, arcs2)
        same_partition = ({frozenset(perm[v] for v in c) for c in ch} == {frozenset(c) for c in ch2})
        sn, sa, chains = bundle_contract(nv, arcs)
        sn2, sa2, chains2 = bundle_contract(nv, arcs2)
        same_contraction = (sorted(chains) == sorted(chains2) == sorted(
            (collections.Counter(u for u, v in arcs)[c[0]], len(c)) for c in ch)
            and sn == sn2 and len(sa) == len(sa2) and elementary(sandpile(sn, sa)) == elementary(sandpile(sn2, sa2)))
        no_chain_left = not maximal_chains(sn, sa) or sn == 1
        good = disjoint and same_partition and same_contraction and no_chain_left
        tot += 1
        ok += good
        if not good:
            print(f"    FAIL n={nv}: disjoint={disjoint} partition={same_partition} "
                  f"contraction={same_contraction} chain-free={no_chain_left}")
    row("E", "maximal bundle chains are pairwise disjoint; their partition, the (r, len) multiset "
             "and the skeleton's group are invariant under relabelling; the skeleton is chain-free",
        f"{tot} random balanced digraphs with planted chains", EXACT, ok == tot)
    # repeat graph: K(RG) = (+)_i Z/r_i (+) K(D_w)
    tot2 = ok2 = 0
    for _ in range(200):
        m = rng.randint(2, 5)
        syms = "ABCDE"[:m]
        mult = {c: rng.choice((2, 2, 3, 4)) for c in syms}
        word = list("".join(c * mult[c] for c in syms))
        rng.shuffle(word)
        word = "".join(word)
        Krg = elementary(sandpile(*repeat_graph(word, mult)))
        Kd = sandpile(*connector_digraph(word, mult))
        pred = elementary(tuple(list(Kd) + [mult[c] for c in syms]))
        tot2 += 1
        ok2 += (Krg == pred)
        if Krg != pred:
            print(f"    MISMATCH word={word} mult={mult}: K(RG)={gstr(Krg)} pred={gstr(pred)}")
    row("E", "the repeat graph (repeats as multi-arcs, unique segments as arcs) has "
             "K(RG) = (+)_i Z/r_i (+) K(D_w)", f"{tot2} random circular words, m <= 5, r_i <= 4",
        EXACT, ok2 == tot2)
    return ok == tot and ok2 == tot2


# --------------------------------------------------------------------------- main

def main():
    t0 = time.time()
    check_clean()
    check_lemma()
    check_special_cases()
    check_chord_invariants()
    check_confluence()
    check_phix()
    print("\n=== battery summary ===")
    bad = [r for r in ROWS if not r[4]]
    for key in ("A", "B", "C", "D", "E", "F"):
        sub = [r for r in ROWS if r[0] == key]
        if sub:
            print(f"  ({key})  {sum(1 for r in sub if r[4])}/{len(sub)} passed")
    print(f"  TOTAL {sum(1 for r in ROWS if r[4])}/{len(ROWS)} passed"
          + ("" if not bad else f"   FAILURES: {[(r[0], r[2]) for r in bad]}"))
    print(f"  elapsed {time.time()-t0:.1f}s")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
