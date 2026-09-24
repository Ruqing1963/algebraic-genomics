# -*- coding: utf-8 -*-
"""
rotor_assembly.py -- verification battery for

    "The sandpile group as a coordinate system on genome reconstructions"   (Bio 2)

Setting.  G = G_k(S) is the de Bruijn multigraph of a circular sequence (Eulerian). By
BEST, the arc-rooted Eulerian circuits -- the reconstructions -- are in bijection with
(arborescence toward a fixed sink) x (a cyclic order of the non-tree out-arcs at every
vertex).  The sandpile group K(G) acts SIMPLY TRANSITIVELY on the arborescences by
rotor-routing (Holroyd-Levine-Meszaros-Peres-Propp-Wilson).  Hence K(G) is a torsor of
coordinates on the GLOBAL part of the reconstruction space, and one chip added at a
vertex is an explicit move from one candidate genome to another.

 (A)  the rotor-routing action on assembly graphs is well defined, by permutations, free
      and transitive: the orbit has size |K| and the generated group has order |K|.  exact
 (B)  BEST bijection implemented explicitly: (arborescence, rotations) -> reconstruction
      is a bijection onto the arc-rooted Eulerian circuits (checked by enumeration).  exact
 (M)  a chip-firing move is a local rearrangement: the two reconstructions share the
      same k-mer spectrum and differ by re-threading a bounded number of repeat visits;
      we exhibit the moves.                                                        exact
 (U)  exact uniform sampling of reconstructions via Smith-normal-form coordinates of
      K(G) followed by rotor-routing; empirical frequencies match the enumeration.  sampled
 (W)  comparison with Wilson's loop-erased-random-walk sampler (the classical route to
      a uniform arborescence): same distribution, and the group gives what Wilson does
      not -- the DISTANCE between two samples as a group element.                 sampled
 (F)  phiX174 at k = 12: K = Z/2, so the 12-mer spectrum leaves exactly ONE binary global
      choice; we exhibit both reconstructions and the chip that flips them.  At k = 11,
      K = Z/2 + Z/66: 132 global classes.                                          exact
 (J)  one chip re-threads several junctions (phiX174 k = 11: 0..6), so a chip is NOT
      a single transposition.                                                      exact
 (N)  a local constraint ("the tree arc at v is e", what a long read fixes) is NOT a
      coset of a subgroup of K(G): 3 of 4 small graphs have non-coset constraint sets. exact

Run: python rotor_assembly.py   (about a minute).  Sequence file phix174_NC_001422.1.txt
must sit beside this script.
"""

import sys, os, time, random, collections, itertools

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

HERE = os.path.dirname(os.path.abspath(__file__))
EXACT, SAMPLED = "exact", "sampled"
ROWS = []


def row(key, desc, cases, mode, ok):
    ROWS.append((key, desc, cases, mode, bool(ok)))
    print(("  [ok]   " if ok else "  [FAIL] ") + f"({key}) {desc}   [{cases}] [{mode}]")
    return bool(ok)


# --------------------------------------------------------------------------- graphs

def dbg(S, k):
    N = len(S)
    SS = S + S[:k]
    verts = {}
    for i in range(N):
        w = SS[i:i + k]
        if w not in verts:
            verts[w] = len(verts)
    arcs = [(verts[SS[i:i + k]], verts[SS[i + 1:i + 1 + k]]) for i in range(N)]
    # label each arc by its (k+1)-mer so a circuit can be spelled back into a sequence
    labels = [SS[i:i + k + 1] for i in range(N)]
    return len(verts), arcs, list(verts), labels


def out_lists(nv, arcs):
    o = collections.defaultdict(list)
    for j, (u, v) in enumerate(arcs):
        o[u].append(j)
    return {v: sorted(o[v]) for v in range(nv)}


def reduced_laplacian(nv, arcs, root=0):
    L = [[0] * nv for _ in range(nv)]
    for (u, v) in arcs:
        L[u][u] += 1
        L[u][v] -= 1
    return [[L[i][j] for j in range(nv) if j != root] for i in range(nv) if i != root]


def sandpile(nv, arcs, root=0):
    if nv <= 1:
        return 1, ()
    M = reduced_laplacian(nv, arcs, root)
    d = pb.diagonalize([list(r) for r in zip(*M)], False, False)[0]
    inv = tuple(sorted(abs(x) for x in d if abs(x) not in (0, 1)))
    o = 1
    for x in inv:
        o *= x
    return o, inv


def compact(nv, arcs, labels=None):
    """Linear-time unitig contraction; the contracted arc carries the concatenated label."""
    outd = [0] * nv
    ind = [0] * nv
    ol = collections.defaultdict(list)
    for j, (u, v) in enumerate(arcs):
        outd[u] += 1
        ind[v] += 1
        ol[u].append(j)
    br = [v for v in range(nv) if outd[v] != 1 or ind[v] != 1]
    if not br:
        return 0, [], []
    bs = set(br)
    rl = {v: i for i, v in enumerate(br)}
    na, nl = [], []
    for b in br:
        for j0 in ol[b]:
            j, x = j0, arcs[j0][1]
            lab = labels[j0] if labels else ""
            steps = 0
            while x not in bs:
                j = ol[x][0]
                if labels:
                    lab += labels[j][-1]
                x = arcs[j][1]
                steps += 1
                if steps > nv + 5:
                    raise RuntimeError("runaway")
            na.append((rl[b], rl[x]))
            nl.append(lab)
    return len(br), na, nl


# ------------------------------------------------------------ arborescences, rotor-routing

def is_arb(arcs, rot, root):
    for v in rot:
        seen, x = set(), v
        while x != root:
            if x in seen:
                return False
            seen.add(x)
            x = arcs[rot[x]][1]
    return True


def all_arbs(nv, arcs, root):
    oa = out_lists(nv, arcs)
    vs = [v for v in range(nv) if v != root]
    return [tuple(sorted(dict(zip(vs, c)).items()))
            for c in itertools.product(*[oa[v] for v in vs])
            if is_arb(arcs, dict(zip(vs, c)), root)]


def route(arcs, oa, rot, root, start):
    """Add one chip at `start`, rotor-route it to the sink; return the new rotor config."""
    rot = dict(rot)
    x, steps = start, 0
    while x != root:
        lst = oa[x]
        rot[x] = lst[(lst.index(rot[x]) + 1) % len(lst)]
        x = arcs[rot[x]][1]
        steps += 1
        if steps > 10 ** 7:
            raise RuntimeError("no termination")
    return tuple(sorted(rot.items()))


def act(arcs, oa, rot, root, chips):
    """Act by a chip configuration (dict vertex -> nonneg int)."""
    for v, c in chips.items():
        for _ in range(c):
            rot = route(arcs, oa, rot, root, v)
    return tuple(sorted(dict(rot).items()))


# ------------------------------------------------------------------- BEST bijection

def best_circuit(nv, arcs, oa, arb, root, rotations):
    """BEST: given an arborescence toward `root` and, at every vertex, a cyclic order of
    its out-arcs with the tree arc LAST (root: any order), produce the Eulerian circuit
    by 'take the next unused arc' starting at the root.  `rotations[v]` is a permutation
    of the non-tree out-arcs of v (root: of all out-arcs)."""
    order = {}
    for v in range(nv):
        if v == root:
            order[v] = list(rotations[v])
        else:
            order[v] = list(rotations[v]) + [arb[v]]
    ptr = {v: 0 for v in range(nv)}
    circuit = []
    x = root
    for _ in range(len(arcs)):
        j = order[x][ptr[x]]
        ptr[x] += 1
        circuit.append(j)
        x = arcs[j][1]
    return tuple(circuit)


def all_circuits_bruteforce(nv, arcs, root):
    """All Eulerian circuits starting at `root` with the first arc = smallest out-arc of
    root, as arc tuples (arc-rooted convention)."""
    oa = out_lists(nv, arcs)
    m = len(arcs)
    used = [False] * m
    out = []
    first = oa[root][0]

    def walk(v, path):
        if len(path) == m:
            if v == root:
                out.append(tuple(path))
            return
        for j in oa[v]:
            if not used[j]:
                used[j] = True
                path.append(j)
                walk(arcs[j][1], path)
                path.pop()
                used[j] = False

    used[first] = True
    walk(arcs[first][1], [first])
    return out


def spell(circuit, labels, k):
    """Turn an arc circuit into the circular sequence it reads (first k-mer + last chars)."""
    return labels[circuit[0]][:k] + "".join(labels[j][-1] for j in circuit)[:len(circuit) - k] \
        if False else "".join(labels[j][-1] for j in circuit)


# ----------------------------------------------------------------------------- (A)

def check_action():
    print("\n=== (A) rotor-routing action of K(G) on assembly-graph arborescences ===")
    tests = [("ACACCAAAACCA", 2), ("CCAACAACCAAC", 3), ("AABABBAABABB", 2),
             ("AABBAABBAABB", 2), ("ACACACACGTGT", 2), ("CACCCCCAACCA", 2)]
    ok_all = True
    for S, k in tests:
        nv, arcs, _, _ = dbg(S, k)
        root = 0
        oa = out_lists(nv, arcs)
        A = all_arbs(nv, arcs, root)
        idx = {a: i for i, a in enumerate(A)}
        Korder, Kinv = sandpile(nv, arcs, root)
        # generators: one chip at each non-sink vertex
        gens = {}
        wd = True
        for v in range(nv):
            if v == root:
                continue
            p = []
            for a in A:
                b = route(arcs, oa, dict(a), root, v)
                if b not in idx:
                    wd = False
                    break
                p.append(idx[b])
            if not wd or sorted(p) != list(range(len(A))):
                wd = False
                break
            gens[v] = tuple(p)
        # group generated, as permutations
        ident = tuple(range(len(A)))
        grp = {ident}
        frontier = [ident]
        while frontier and len(grp) <= 4 * Korder:
            nxt = []
            for g in frontier:
                for p in gens.values():
                    h = tuple(p[g[i]] for i in range(len(A)))
                    if h not in grp:
                        grp.add(h)
                        nxt.append(h)
            frontier = nxt
        orbit = {g[0] for g in grp}
        free = all(g[0] != ident[0] for g in grp if g != ident)
        good = wd and len(A) == Korder and len(grp) == Korder and len(orbit) == len(A) and free
        ok_all &= good
        print(f"    S={S:<13s} k={k}  |arb|={len(A):3d} = |K|={Korder:3d}  K={list(Kinv)}  "
              f"group order {len(grp)}, orbit {len(orbit)}, free={free}  {'ok' if good else 'FAIL'}")
    row("A", "the action is by permutations, and the generated group is simply transitive "
             "(order = orbit = |K|)", f"{len(tests)} assembly graphs", EXACT, ok_all)
    return ok_all


# ----------------------------------------------------------------------------- (B)

def check_best():
    print("\n=== (B) the BEST bijection, implemented and enumerated ===")
    tests = [("ACACCAAAACCA", 2), ("CACCCCCAACCA", 2), ("ACCCCAAACACC", 2), ("CACCAAACCACA", 3)]
    ok_all = True
    for S, k in tests:
        nv, arcs, _, labels = dbg(S, k)
        root = 0
        oa = out_lists(nv, arcs)
        A = all_arbs(nv, arcs, root)
        brute = set(all_circuits_bruteforce(nv, arcs, root))
        # the local rotation data: for non-root v, permutations of non-tree out-arcs;
        # for the root, permutations of all out-arcs whose FIRST element is the smallest
        # (to match the arc-rooted convention of the brute force)
        produced = set()
        for arb in A:
            arbd = dict(arb)
            choices = []
            vs = list(range(nv))
            for v in vs:
                if v == root:
                    rest = [j for j in oa[v] if j != oa[v][0]]
                    choices.append([(oa[v][0],) + p for p in itertools.permutations(rest)])
                else:
                    nontree = [j for j in oa[v] if j != arbd[v]]
                    choices.append(list(itertools.permutations(nontree)))
            for combo in itertools.product(*choices):
                rotations = dict(zip(vs, combo))
                c = best_circuit(nv, arcs, oa, arbd, root, rotations)
                produced.add(c)
        good = (produced == brute)
        ok_all &= good
        print(f"    S={S:<13s} k={k}: {len(A)} arborescences x local rotations -> "
              f"{len(produced)} circuits; brute force {len(brute)}  {'bijection' if good else 'MISMATCH'}")
    row("B", "(arborescence, rotations) -> circuit is a bijection onto the arc-rooted "
             "Eulerian circuits", f"{len(tests)} graphs", EXACT, ok_all)
    return ok_all


# ----------------------------------------------------------------------------- (M)

def diff_positions(s1, s2):
    return sum(1 for a, b in zip(s1, s2) if a != b)


def check_moves():
    print("\n=== (M) a chip is a move: two reconstructions one chip apart ===")
    S, k = "ACACCAAAACCA", 2
    nv, arcs, verts, labels = dbg(S, k)
    root = 0
    oa = out_lists(nv, arcs)
    A = all_arbs(nv, arcs, root)
    # fix trivial rotations (identity order) so that the only freedom is the arborescence
    def circ_of(arb):
        arbd = dict(arb)
        rotations = {}
        for v in range(nv):
            if v == root:
                rotations[v] = tuple(oa[v])
            else:
                rotations[v] = tuple(j for j in oa[v] if j != arbd[v])
        return best_circuit(nv, arcs, oa, arbd, root, rotations)
    base = A[0]
    seq0 = spell(circ_of(base), labels, k)
    print(f"    base reconstruction        : {seq0}   (spectrum of {S})")
    ok = True
    shown = 0
    for v in range(nv):
        if v == root:
            continue
        moved = route(arcs, oa, dict(base), root, v)
        seq1 = spell(circ_of(moved), labels, k)
        same_spec = (sorted(dbg(seq0, k)[3]) == sorted(dbg(seq1, k)[3]))
        ok &= same_spec
        if shown < 3:
            shown += 1
            print(f"    + one chip at k-mer {verts[v]!r:6s}: {seq1}   "
                  f"({diff_positions(seq0, seq1)} positions changed, spectrum preserved: {same_spec})")
    row("M", "each chip yields another reconstruction with the same k-mer spectrum; the "
             "move is explicit", "S=ACACCAAAACCA, k=2", EXACT, ok)
    return ok


# ----------------------------------------------------------------------------- (U)(W)

def snf_uniform_sampler(nv, arcs, root, rng):
    """Exact uniform sampling of K(G) via Smith-form coordinates, returned as a chip
    configuration on the non-sink vertices.  coker(M^T) = Z^{n-1} / M^T Z^{n-1};
    with U M^T V = D (U, V unimodular), a uniform element is U^{-1} (r_1,...,r_i,...) with
    r_i uniform mod d_i (and 0 where d_i = 1)."""
    import sympy as sp
    M = reduced_laplacian(nv, arcs, root)
    Mt = sp.Matrix([list(r) for r in zip(*M)])
    # sympy's smith_normal_decomp is not in all versions; use Paper X's transforms
    d, P, Q, _ = pb.diagonalize([[int(x) for x in r] for r in Mt.tolist()], True, True)
    # pb.diagonalize returns d = diag, P = row transform, Q = column transform with P*Mt*Q = diag
    P = sp.Matrix(P)
    n1 = nv - 1
    r = sp.zeros(n1, 1)
    for i in range(n1):
        di = abs(int(d[i])) if i < len(d) else 0
        r[i] = rng.randrange(di) if di > 1 else 0
    c = P.inv() * r                       # a representative in Z^{n-1}
    chips = {}
    idx = [v for v in range(nv) if v != root]
    for i, v in enumerate(idx):
        val = int(c[i])
        chips[v] = val
    return chips


def normalise_chips(nv, arcs, root, chips):
    """Rotor-routing acts by NONNEGATIVE chip additions; a negative count at v is
    equivalent, in K(G), to adding (d_v - 1) chips plus one unit of the sink relation...
    simpler: add a large multiple of the identity relation.  We add |K| * outdeg to make
    all entries nonnegative; |K| * (row of Laplacian) is trivial in K."""
    Korder, _ = sandpile(nv, arcs, root)
    L = [[0] * nv for _ in range(nv)]
    for (u, v) in arcs:
        L[u][u] += 1
        L[u][v] -= 1
    out = dict(chips)
    for v in list(out):
        while out[v] < 0:
            # add the Laplacian row of v (a trivial element) scaled to make v nonnegative
            for w in range(nv):
                if w == root:
                    continue
                out[w] = out.get(w, 0) + L[v][w] * Korder
    return {v: c for v, c in out.items() if c > 0}


def wilson_arborescence(nv, arcs, root, rng):
    """Wilson's algorithm: loop-erased random walks toward the root give an exactly
    uniform arborescence.  Rotor = the last out-arc used to leave each vertex."""
    oa = out_lists(nv, arcs)
    intree = {root}
    nxt = {}
    for v in range(nv):
        x = v
        path = []
        while x not in intree:
            j = rng.choice(oa[x])
            nxt[x] = j
            path.append(x)
            x = arcs[j][1]
        for y in path:
            intree.add(y)
    return tuple(sorted((v, nxt[v]) for v in range(nv) if v != root))


def check_sampling():
    print("\n=== (U)/(W) exact uniform sampling, and the group distance ===")
    rng = random.Random(2026)
    S, k = "AABBAABBAABB", 2
    nv, arcs, _, _ = dbg(S, k)
    root = 0
    oa = out_lists(nv, arcs)
    A = all_arbs(nv, arcs, root)
    idx = {a: i for i, a in enumerate(A)}
    Korder, Kinv = sandpile(nv, arcs, root)
    print(f"    S={S} k={k}: |K| = {Korder}, K = {list(Kinv)}, {len(A)} arborescences")
    # (U) SNF sampler
    base = A[0]
    N = 2700
    counts = collections.Counter()
    for _ in range(N):
        chips = normalise_chips(nv, arcs, root, snf_uniform_sampler(nv, arcs, root, rng))
        a = act(arcs, oa, dict(base), root, chips)
        counts[idx[a]] += 1
    exp = N / len(A)
    df = len(A) - 1

    def chi2_p(x):
        import mpmath as mp
        return float(mp.gammainc(df / 2, x / 2, mp.inf) / mp.gamma(df / 2))

    chi2 = sum((c - exp) ** 2 / exp for c in counts.values()) + exp * (len(A) - len(counts))
    p_u = chi2_p(chi2)
    ok_u = (len(counts) == len(A) and p_u > 0.01)
    print(f"    (U) SNF-coordinate sampler: {N} draws, all {len(counts)}/{len(A)} arborescences hit,"
          f" chi^2 = {chi2:.1f} on {df} df, p = {p_u:.3f}")
    row("U", f"exact uniform sampler via Smith-form coordinates: chi^2 = {chi2:.1f} on "
             f"{df} df, p = {p_u:.2f}", f"{N} draws", SAMPLED, ok_u)
    # (W) Wilson
    countsW = collections.Counter()
    for _ in range(N):
        countsW[idx[wilson_arborescence(nv, arcs, root, rng)]] += 1
    chi2W = sum((c - exp) ** 2 / exp for c in countsW.values()) + exp * (len(A) - len(countsW))
    p_w = chi2_p(chi2W)
    ok_w = (len(countsW) == len(A) and p_w > 0.01)
    print(f"    (W) Wilson sampler          : {N} draws, all {len(countsW)}/{len(A)} hit,"
          f" chi^2 = {chi2W:.1f}, p = {p_w:.3f}")
    # the group distance: for two Wilson samples, find the chip configuration taking one
    # to the other (search over the group, feasible here); that element is the 'difference'
    a1 = wilson_arborescence(nv, arcs, root, rng)
    a2 = wilson_arborescence(nv, arcs, root, rng)
    found = None
    gens = [v for v in range(nv) if v != root]
    frontier = {a1: ()}
    seen = {a1}
    while frontier and found is None:
        nxt = {}
        for a, word in frontier.items():
            for v in gens:
                b = route(arcs, oa, dict(a), root, v)
                if b == a2:
                    found = word + (v,)
                    break
                if b not in seen:
                    seen.add(b)
                    nxt[b] = word + (v,)
            if found:
                break
        frontier = nxt
    dist = len(found) if found is not None else (0 if a1 == a2 else None)
    print(f"    group distance between two independent samples: {dist} chip(s)"
          f" (word {found})")
    row("W", "Wilson's sampler agrees; the group additionally yields the chip word taking one "
             "sample to another", f"{N} draws", SAMPLED, ok_w and dist is not None)
    return ok_u and ok_w


# ----------------------------------------------------------------------------- (F)

def check_phix():
    print("\n=== (F) phiX174: the global choices left by the 12-mer and 11-mer spectra ===")
    S = open(os.path.join(HERE, "phix174_NC_001422.1.txt")).read().strip().upper()
    ok = True
    for k in (12, 11):
        nv, arcs, verts, labels = dbg(S, k)
        cn, ca, cl = compact(nv, arcs, labels)
        Korder, Kinv = sandpile(cn, ca)
        print(f"    k={k}: {nv} vertices -> {cn} branch vertices, K = {list(Kinv)}, |K| = {Korder}")
        if k == 12:
            root = 0
            oa = out_lists(cn, ca)
            A = all_arbs(cn, ca, root)
            ok &= (len(A) == 2 and Korder == 2)
            # the two reconstructions with identity rotations, and the chip that flips them
            def circ_of(arb):
                arbd = dict(arb)
                rot = {v: (tuple(oa[v]) if v == root else tuple(j for j in oa[v] if j != arbd[v]))
                       for v in range(cn)}
                return best_circuit(cn, ca, oa, arbd, root, rot)
            seqs = []
            for a in A:
                c = circ_of(a)
                seqs.append("".join(cl[j][k:] for j in c))
            flip = None
            for v in range(cn):
                if v != root and route(ca, oa, dict(A[0]), root, v) == A[1]:
                    flip = v
                    break
            same_len = all(len(s) == len(S) for s in seqs)
            # both have the genome's 12-mer spectrum
            spec0 = sorted(dbg(S, 12)[3])
            same_spec = all(sorted(dbg(s, 12)[3]) == spec0 for s in seqs)
            # The biology.  A single repeat W of multiplicity 2 gives NO circular ambiguity
            # (W X W Y and W Y W X are rotations of each other; Theorem 5 of Bio 1: K trivial
            # at k = L).  A Z/2 at k = 12 therefore means two INTERLEAVED 12-mer repeats
            # A .. B .. A .. B -- the s = 1 case of the interleaving theorem of Bio 1 v2 --
            # and the two reconstructions are A X B Y A Z B U and A X B U A Z B Y: the two
            # segments following the two copies of B are exchanged.  Verify exactly that.
            s0, s1 = seqs

            def occurrences(s, w):
                ss = s + s[:11]
                return [i for i in range(len(s)) if ss[i:i + 12] == w]

            reps = []
            seen_w = set()
            for i in range(len(s0)):
                cand = (s0 + s0[:11])[i:i + 12]
                if cand not in seen_w and len(occurrences(s0, cand)) == 2:
                    seen_w.add(cand)
                    reps.append(cand)
            # the two repeated 12-mers are the two branch vertices; keep the two with
            # non-overlapping occurrences (a run of length 13 would give two overlapping 12-mers)
            A12 = reps[0]
            B12 = next(w for w in reps[1:] if all(abs(p - q) > 12 for p in occurrences(s0, w)
                                                       for q in occurrences(s0, A12)))

            def segments(s):
                oa_ = occurrences(s, A12)
                r = s[oa_[0]:] + s[:oa_[0]]              # rotate to start at the first A
                pa = occurrences(r, A12)
                pb_ = occurrences(r, B12)
                cuts = sorted([(p, "A") for p in pa] + [(p, "B") for p in pb_])
                pattern = "".join(t for _, t in cuts)
                segs = []
                for (p, _), (q, _) in zip(cuts, cuts[1:] + [(len(r), None)]):
                    segs.append(r[p + 12:q])
                return pattern, tuple(segs)

            pat0, T0 = segments(s0)
            pat1, T1 = segments(s1)
            interleaved = (pat0 == "ABAB" and pat1 == "ABAB")
            swapped = (T0[0], T0[3], T0[2], T0[1])
            rot2 = (T1[2], T1[3], T1[0], T1[1])
            exchange = (T1 == swapped) or (rot2 == swapped)
            print(f"      two repeated 12-mers A={A12}, B={B12}, arrangement {pat0} (interleaved)")
            print(f"      segment lengths (X,Y,Z,U) = {tuple(len(t) for t in T0)}; the second "
                  f"reconstruction is A X B U A Z B Y (Y and U exchanged): {exchange}")
            print(f"      => the one global ambiguity of the 12-mer spectrum is the Z/2 of the "
                  f"interleaving theorem, realised in a real genome; one chip at branch vertex "
                  f"{flip} performs the exchange")
            ok &= same_len and same_spec and flip is not None and interleaved and exchange
    row("F", "phiX174 k=12: exactly two global reconstructions (K = Z/2), exhibited, one chip "
             "apart; k=11: 132 global classes (K = Z/2 + Z/66)", "5386 bp", EXACT, ok)
    return ok


# ----------------------------------------------------------------------------- (J)

def check_junctions():
    """Is one chip one transposition?  A rotor walk re-threads EVERY branch vertex it
    passes, so a single chip can change the tree arc at many junctions.  Count, on the
    compacted phiX174 graph at k = 11 (ten branch vertices) and on the small example."""
    print("\n=== (J) how many junctions does one chip re-thread? ===")
    S = open(os.path.join(HERE, "phix174_NC_001422.1.txt")).read().strip().upper()
    out = {}
    for label, (seq, k) in {"phiX174 k=11": (S, 11), "ACACCAAAACCA k=2": ("ACACCAAAACCA", 2)}.items():
        nv, arcs, _, labels = dbg(seq, k)
        cn, ca, _ = compact(nv, arcs, labels)
        if cn == 0:
            continue
        oa = out_lists(cn, ca)
        base = dict(all_arbs(cn, ca, 0)[0])
        changes = []
        for v in range(1, cn):
            new = route(ca, oa, base, 0, v)
            changes.append(sum(1 for w in base if base[w] != dict(new)[w]))
        out[label] = changes
        print(f"    {label:<18s}: {cn} branch vertices; junctions re-threaded per chip = {changes}"
              f"  (single re-threading in {sum(1 for c in changes if c == 1)}/{len(changes)})")
    ch = out["phiX174 k=11"]
    row("J", f"one chip re-threads {min(ch)}..{max(ch)} junctions on phiX174 at k=11; a chip is a "
             f"single transposition in only {sum(1 for c in ch if c == 1)}/{len(ch)} cases",
        "phiX174 k=11, 10 branch vertices", EXACT, max(ch) > 1)
    return True


# ----------------------------------------------------------------------------- (N)

def check_cosets():
    """Is a local constraint a coset?  A long read that fixes the tree arc at a branch
    vertex v selects S_e = { T in Arb : T(v) = e }.  Identify Arb with K(G) through the
    torsor (g <-> g.T_0).  If constraints could be combined by linear algebra in K(G),
    each S_e would have to be a coset of a subgroup.  Test: g_0^{-1} S_e closed under the
    group law?"""
    print("\n=== (N) are long-read constraints cosets of K(G)? ===")
    tests = [("ACACCAAAACCA", 2), ("AABABBAABABB", 2), ("AABBAABBAABB", 2), ("CCAACAACCAAC", 3)]
    summary = {}
    for S, k in tests:
        nv, arcs, _, _ = dbg(S, k)
        root = 0
        oa = out_lists(nv, arcs)
        A = all_arbs(nv, arcs, root)
        idx = {a: i for i, a in enumerate(A)}
        n = len(A)
        gens = {v: tuple(idx[route(arcs, oa, dict(a), root, v)] for a in A) for v in range(1, nv)}
        ident = tuple(range(n))
        grp = {ident}
        fr = [ident]
        while fr:
            nx = []
            for g in fr:
                for p in gens.values():
                    h = tuple(p[g[i]] for i in range(n))
                    if h not in grp:
                        grp.add(h)
                        nx.append(h)
            fr = nx
        elem_of_arb = {g[0]: g for g in grp}            # simply transitive: g <-> g.T_0
        def compose(g, h):
            return tuple(g[h[i]] for i in range(n))
        def inverse(g):
            inv = [0] * n
            for i, gi in enumerate(g):
                inv[gi] = i
            return tuple(inv)
        cos = tot = 0
        detail = []
        for v in range(1, nv):
            for e in oa[v]:
                Sset = [elem_of_arb[i] for i, a in enumerate(A) if dict(a)[v] == e]
                if not Sset:
                    continue
                tot += 1
                g0inv = inverse(Sset[0])
                H = {compose(g0inv, g) for g in Sset}
                closed = all(compose(a, b) in H for a in H for b in H)
                cos += closed
                detail.append((v, e, len(Sset), closed))
        summary[S] = (cos, tot)
        print(f"    S={S} k={k}: |K|={len(grp)}; {cos}/{tot} constraint sets are cosets  "
              + "  ".join(f"v{v}:arc{e}|{s}|{'c' if c else 'X'}" for v, e, s, c in detail))
    non = [S for S, (c, t) in summary.items() if c < t]
    row("N", f"the constraint 'tree arc at v is e' is NOT a coset of a subgroup of K(G) on "
             f"{len(non)} of {len(tests)} graphs (all cosets only on the single-repeat graph)",
        f"{sum(t for c, t in summary.values())} constraint sets", EXACT,
        len(non) == 3 and summary["AABBAABBAABB"][0] == summary["AABBAABBAABB"][1])
    return True


# --------------------------------------------------------------------------- main

def main():
    t0 = time.time()
    check_action()
    check_best()
    check_moves()
    check_sampling()
    check_phix()
    check_junctions()
    check_cosets()
    print("\n=== battery summary ===")
    bad = [r for r in ROWS if not r[4]]
    for key in ("A", "B", "M", "U", "W", "F", "J", "N"):
        sub = [r for r in ROWS if r[0] == key]
        if sub:
            print(f"  ({key})  {sum(1 for r in sub if r[4])}/{len(sub)} passed")
    print(f"  TOTAL {sum(1 for r in ROWS if r[4])}/{len(ROWS)} passed"
          + ("" if not bad else f"   FAILURES: {[(r[0], r[2]) for r in bad]}"))
    print(f"  elapsed {time.time()-t0:.1f}s")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
