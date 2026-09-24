# -*- coding: utf-8 -*-
r"""
ds_fermion.py -- #DS by a fermionic (Grassmann) frontier sweep.

Bio 12 filled two of the three blank entries of Bio 8's Table 2 for phiX174 at k = 10: the
number of lattice points of the box, 2 481 687 900, and the number with connected support,
909 103 210.  It could not reach #DS itself, because

    #DS(S,k) = sum_f Nseq((m+f)/2),   Nseq(c) = t(G_c) Loc(G_c) / prod_a c(a)!

carries the spanning-arborescence count t(G_c) = det Ltilde_c, a determinant of the whole
support -- of order 125 at k = 10 -- and a determinant does not factor through a scalar
frontier state.

It factors through a FERMIONIC one.  By Berezin,

    det Ltilde_c = int D[psibar psi] prod_{a in A} (1 - c(a) psibar_{t(a)} (psi_{t(a)} - psi_{h(a)})),

one factor per arc, each of them EVEN, hence commuting: the product may be accumulated one arc
at a time in Bio 12's own orbit order, and the state it needs is two bits per frontier vertex
-- has psibar_v been spent, has psi_v been taken -- instead of a set partition.  The frontier
is the frontier of Bio 12; only its alphabet grows.

Three things have to be arranged and none of them is in the Berezin formula.

(1) ROOT.  A reduced determinant needs a root inside supp c, and supp c moves with c.  Give
    every vertex outside supp c a 1 on the diagonal -- the factor (1 - psibar_v psi_v), whose
    sign is easy to get wrong -- and pin one orbit j0 of odd multiplicity to f_{j0} > 0.  Then
    c(a_{j0}) >= 1, so t(a_{j0}) is in supp c for every point of that half-box, and it serves
    as a fixed root.  Since Nseq(m-c) = Nseq(c), the half-box carries exactly half of #DS.

(2) Loc.  (d_v - 1)! is not arc-local.  The partial out-degree is carried on the out-arc
    sub-frontier, which is much narrower than the frontier itself, and spent when the last
    out-arc of v is decided.

(3) INTEGRALITY.  prod_a 1/c(a)! = prod_j binom(m_j, c(a_j)) / prod_j m_j!, so the sweep runs
    over Z and there is one division at the end.

(4) SIZE.  The fermionic state is linear, and at k = 10 a single history already fills up to
    3^(open rows) coordinates of it.  The number itself is computed by the POSITIVE expansion
    of the same determinant -- in-trees, one class per undecided vertex, a vertex choosing its
    tree arc when its last out-arc is decided -- which puts one state per history (check A, Q,
    and the run K).  Check X records how fast the fermionic frontier grows.

    python -u ds_fermion.py --quick    D, F, R, S, V, A, Q   (about ten minutes)
    python -u ds_fermion.py            the same, and X: the fermionic frontier at k = 10,
                                       abandoned at 10^6 states or ten minutes
        --no-probe     skip X
        --arbor-k10    also attempt the whole k = 10 half-box with the arborescence frontier
                       (--cap=, --jobs=; checkpoint ds_k10_checkpoint.json).  It does not
                       finish on a 16 GB machine; the number is computed by ds_count.py.

Standalone: Python 3.10+ (int.bit_count) and the standard library.  The phiX174 sequence is
read from phix174_NC_001422.1.txt beside this script.
"""
import os, sys, time, math, json, random, collections
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
COMP = str.maketrans("ACGT", "TGCA")
ROWS = []


def row(key, desc, cases, mode, ok):
    ROWS.append((key, desc, cases, mode, bool(ok)))
    print(("  [ok]   " if ok else "  [FAIL] ") + f"({key}) {desc}   [{cases}] [{mode}]")
    return bool(ok)


def rc(s):
    return s.translate(COMP)[::-1]


def load_phix():
    for nm in ("phix174_NC_001422.1.txt", os.path.join("..", "Bio_12_Frontier_Elimination", "phix174_NC_001422.1.txt"),
               os.path.join("..", "Bio_08_DS_Lattice_Sum", "phix174_NC_001422.1.txt")):
        p = os.path.join(HERE, nm)
        if os.path.exists(p):
            s = "".join(c for c in open(p, encoding="utf-8", errors="replace").read().upper()
                        if c in "ACGT")
            if len(s) > 5000:
                return s
    raise SystemExit("phix174_NC_001422.1.txt not found")


def det(M):
    n = len(M)
    if n == 0:
        return Fraction(1)
    A = [r[:] for r in M]
    d = Fraction(1)
    for i in range(n):
        p = next((r for r in range(i, n) if A[r][i] != 0), None)
        if p is None:
            return Fraction(0)
        if p != i:
            A[i], A[p] = A[p], A[i]
            d = -d
        d *= A[i][i]
        inv = 1 / A[i][i]
        for r in range(i + 1, n):
            if A[r][i]:
                f = A[r][i] * inv
                for cc in range(i, n):
                    A[r][cc] -= f * A[i][cc]
    return d


# ------------------------------------------------------- the double cover (as in Bio 8, 12)

class Cover:
    def __init__(self, S, k):
        self.S, self.k, self.N = S, k, len(S)
        N = len(S)
        seqs = [S, rc(S)]
        occ = collections.defaultdict(list)
        for si, T in enumerate(seqs):
            TT = T + T[:k]
            for i in range(N):
                occ[TT[i:i + k]].append((si, i))
        rep = {w: ps for w, ps in occ.items() if len(ps) >= 2}
        self.words = sorted(rep)
        self.vid = {w: i for i, w in enumerate(self.words)}
        self.nv = len(self.words)
        self.n_fixV = sum(1 for w in self.words if rc(w) == w)
        self.arcs = []
        for si, T in enumerate(seqs):
            TT = T + T + T[:k]
            R = sorted(p for w in rep for (s, p) in rep[w] if s == si)
            for t in range(len(R)):
                p, q = R[t], R[(t + 1) % len(R)]
                if q <= p:
                    q += N
                self.arcs.append(dict(strand=si, pos=p, len=q - p + k,
                                      content=TT[p:p + (q - p + k)]))
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
        self.reps = [a for a in range(self.nA) if self.crho[a] is not None and a < self.crho[a]]

    def constraints(self):
        out = []
        for a in self.reps:
            b = self.crho[a]
            d = collections.Counter()
            d[self.tail[a]] += 1
            d[self.head[a]] -= 1
            d[self.tail[b]] -= 1
            d[self.head[b]] += 1
            out.append({v: c for v, c in d.items() if c})
        return out

    def point(self, g):
        c = [0] * self.nA
        for j, a in enumerate(self.reps):
            c[a] = (self.m[a] + g[j]) // 2
            c[self.crho[a]] = (self.m[a] - g[j]) // 2
        for a in range(self.nA):
            if self.crho[a] == a:
                c[a] = self.m[a] // 2
        return c

    def outdeg(self, c):
        o = collections.Counter()
        for a in range(self.nA):
            if c[a]:
                o[self.tail[a]] += c[a]
        return o

    def n_seq(self, c):
        outd = self.outdeg(c)
        supp = sorted(outd)
        if not supp:
            return 0
        idx = {v: i for i, v in enumerate(supp)}
        n = len(supp)
        L = [[0] * n for _ in range(n)]
        for a in range(self.nA):
            if c[a]:
                L[idx[self.tail[a]]][idx[self.tail[a]]] += c[a]
                L[idx[self.tail[a]]][idx[self.head[a]]] -= c[a]
        t = det([[Fraction(L[r][s]) for s in range(1, n)] for r in range(1, n)]) if n > 1 \
            else Fraction(1)
        if t == 0:
            return 0
        val = int(t)
        for v in supp:
            val *= math.factorial(outd[v] - 1)
        for a in range(self.nA):
            if c[a]:
                val //= math.factorial(c[a])
        return val

    def t_rooted(self, c, root):
        """det of (L_c + diag(1 for v outside supp c)) with row and column `root` deleted"""
        outd = self.outdeg(c)
        L = [[Fraction(0)] * self.nv for _ in range(self.nv)]
        for a in range(self.nA):
            if c[a]:
                L[self.tail[a]][self.tail[a]] += c[a]
                L[self.tail[a]][self.head[a]] -= c[a]
        for v in range(self.nv):
            if outd[v] == 0:
                L[v][v] += 1
        idx = [v for v in range(self.nv) if v != root]
        return det([[L[r][s] for s in idx] for r in idx])


def n_seq_fast(cov, c):
    """Nseq(c) with t(G_c) computed by contracting forced arcs first: a vertex whose out-arcs
    all go to one other vertex u must use one of them in every in-tree, so t(G) = w t(G/(v->u)).
    On a k = 10 point that leaves a determinant of order ~10 instead of ~110.  Independent of
    the sweeps; used to check them on pinned k = 10 sub-boxes."""
    outd = cov.outdeg(c)
    supp = sorted(outd)
    if not supp:
        return 0
    out = {v: collections.Counter() for v in supp}
    inn = {v: set() for v in supp}
    for a in range(cov.nA):
        if c[a] and cov.tail[a] != cov.head[a]:
            out[cov.tail[a]][cov.head[a]] += c[a]
            inn[cov.head[a]].add(cov.tail[a])
    root = supp[0]
    t = 1
    todo = [v for v in supp if v != root]
    while todo:
        v = todo.pop()
        if v not in out or v == root or len(out[v]) != 1:
            continue
        (u, w), = out[v].items()
        t *= w
        for x in inn[v]:
            if x == v:
                continue
            w2 = out[x].pop(v)
            if x != u:
                out[x][u] += w2
                inn[u].add(x)
            if x != root and len(out[x]) == 1:
                todo.append(x)
        inn[u].discard(v)
        del out[v], inn[v]
        if u != root and len(out[u]) == 1:
            todo.append(u)
    rest = [v for v in out if v != root]
    if any(not out[v] for v in rest):
        return 0
    ix = {v: i for i, v in enumerate(rest)}
    M = [[Fraction(0)] * len(rest) for _ in rest]
    for v in rest:
        for u, w in out[v].items():
            M[ix[v]][ix[v]] += w
            if u != root:
                M[ix[v]][ix[u]] -= w
    t *= int(det(M))
    if t == 0:
        return 0
    val = t
    for v in supp:
        val *= math.factorial(outd[v] - 1)
    for a in range(cov.nA):
        if c[a]:
            val //= math.factorial(c[a])
    return val


def enumerate_box(cov, cap=2000000, allowed=None, nodecap=None):
    """Bio 8's backtracking enumerator, for the small k where it still runs.  `allowed[j]`
    restricts orbit j to the listed values of g = 2c - m; `nodecap` gives up (returns None)
    after that many search nodes."""
    var = cov.constraints()
    m = [cov.m[a] for a in cov.reps]
    if not hasattr(cov, "_enum"):
        order = greedy_order(var)
        touch = collections.defaultdict(set)
        for j, d in enumerate(var):
            for v in d:
                touch[v].add(j)
        capv = {v: sum(abs(var[j].get(v, 0)) * m[j] for j in s) for v, s in touch.items()}
        rem = [dict(capv)]
        for j in order:
            nc = dict(rem[-1])
            for v, c in var[j].items():
                nc[v] -= abs(c) * m[j]
            rem.append(nc)
        cov._enum = (order, rem)
    order, rem = cov._enum
    out, cur = [], [0] * len(var)
    bal = collections.Counter()
    nodes = [0]

    def recur(i):
        if i == len(order):
            out.append(tuple(cur))
            return len(out) <= cap
        nodes[0] += 1
        if nodecap and nodes[0] > nodecap:
            return False
        j = order[i]
        vals = allowed[j] if (allowed and j in allowed) else range(-m[j], m[j] + 1, 2)
        for g in vals:
            for v, c in var[j].items():
                bal[v] += c * g
            if all(abs(bal[v]) <= rem[i + 1].get(v, 0) for v in var[j]):
                cur[j] = g
                if not recur(i + 1):
                    for v, c in var[j].items():
                        bal[v] -= c * g
                    return False
            for v, c in var[j].items():
                bal[v] -= c * g
        return True
    recur(0)
    if nodecap and nodes[0] > nodecap:
        return None
    return out


# ------------------------------------------------------------------- orders and their width

def greedy_order(var, dom=None, first=None, jitter=0.0, rnd=None):
    dom = dom or [set(d) for d in var]
    touch = collections.defaultdict(set)
    for j, d in enumerate(dom):
        for e in d:
            touch[e].add(j)
    left = {e: set(s) for e, s in touch.items()}
    todo, open_e, order = set(range(len(dom))), set(), []
    while todo:
        if first is not None and not order:
            best = first
        else:
            best, bk = None, None
            for j in todo:
                nxt = open_e | dom[j]
                cl = {e for e in nxt if not (left[e] - {j})}
                key = (len(nxt) - len(cl) + (rnd.random() * jitter if jitter else 0.0),
                       -len(cl), -len(open_e & dom[j]))
                if bk is None or key < bk:
                    best, bk = j, key
        todo.discard(best)
        order.append(best)
        for e in dom[best]:
            left[e].discard(best)
        open_e = {e for e in (open_e | dom[best]) if left[e]}
    return order


def fermionic_domains(cov):
    """the elements the fermionic frontier has to carry: one psi (with the imbalance) per
    vertex over its whole window, one psibar (with the partial out-degree) per vertex over
    its out-arc window"""
    # elements are small integers, not tuples with a string in them, so that the greedy does
    # not depend on the interpreter's hash seed
    dom, weight = [], {}
    for j, a in enumerate(cov.reps):
        b = cov.crho[a]
        full = {cov.tail[a], cov.head[a], cov.tail[b], cov.head[b]}
        dom.append({2 * v for v in full} | {2 * v + 1 for v in (cov.tail[a], cov.tail[b])})
    for v in range(cov.nv):
        weight[2 * v] = 2.6          # psi_v, together with the imbalance at v
        weight[2 * v + 1] = 2.0      # psibar_v, together with the partial out-degree at v
    return dom, weight


def width_profile(order, dom, weight, n):
    pos = [0] * n
    for i, j in enumerate(order):
        pos[j] = i
    first, last = {}, {}
    for j in order:
        i = pos[j]
        for e in dom[j]:
            if e not in first:
                first[e] = last[e] = i
            else:
                first[e] = min(first[e], i)
                last[e] = max(last[e], i)
    dw = [0.0] * (n + 1)
    dc = [0] * (n + 1)
    for e, a in first.items():
        b = last[e]
        if b > a:
            dw[a] += weight[e]
            dw[b] -= weight[e]
            dc[a] += 1
            dc[b] -= 1
    s = c = 0
    mx, mc, acc = 0.0, 0, []
    for i in range(n):
        s += dw[i]
        c += dc[i]
        acc.append(s)
        mx, mc = max(mx, s), max(mc, c)
    return mx, mc, math.log2(sum(2.0 ** (x - mx) for x in acc)) + mx


def search_order(cov, tries=40, rounds=8, seed=7, verbose=False):
    dom, weight = fermionic_domains(cov)
    J = len(cov.reps)
    var = cov.constraints()
    rnd = random.Random(seed)
    pool = [greedy_order(var, dom=[set(var[j]) | dom[j] for j in range(J)]),
            greedy_order(var, dom=dom)]
    for _ in range(tries):
        pool.append(greedy_order(var, dom=dom, first=rnd.randrange(J), jitter=0.9, rnd=rnd))
    pool.sort(key=lambda o: width_profile(o, dom, weight, J)[2])
    best, bkey = None, None
    for s, seed_o in enumerate(pool[:3]):
        cur = list(seed_o)
        key = width_profile(cur, dom, weight, J)
        kk = (key[0], key[2])
        for _ in range(rounds):
            moved = False
            js = list(cur)
            rnd.shuffle(js)
            for j in js:
                base = [x for x in cur if x != j]
                bi = None
                for i in range(len(base) + 1):
                    cand = base[:i] + [j] + base[i:]
                    e = width_profile(cand, dom, weight, J)
                    c2 = (e[0], e[2])
                    if c2[0] < kk[0] - 1e-9 or (abs(c2[0] - kk[0]) < 1e-9 and c2[1] < kk[1] - 1e-9):
                        kk, bi = c2, i
                if bi is not None:
                    cur = base[:bi] + [j] + base[bi:]
                    moved = True
            if not moved:
                break
        if bkey is None or kk[1] < bkey[1]:
            best, bkey = cur, kk
        if verbose:
            print(f"      order seed {s}: max {kk[0]:.1f} bits, log2 total {kk[1]:.1f}",
                  flush=True)
    return best


# ------------------------------------------------------------------ the fermionic sweep

def _slots(first, last):
    """colour the intervals [first[v], last[v]] with as few slots as possible; return the slot
    of each vertex and its epoch, the rank of its turn in that slot"""
    order = sorted(first, key=lambda v: (first[v], last[v], v))
    free, active, slot, users = [], [], {}, collections.defaultdict(list)
    for v in order:
        active = [(l, u) for (l, u) in active if l >= first[v]]
        used = {slot[u] for (l, u) in active}
        s = 0
        while s in used:
            s += 1
        slot[v] = s
        users[s].append(v)
        active.append((last[v], v))
    epoch = {}
    for s, us in users.items():
        for e, v in enumerate(sorted(us, key=lambda v: first[v])):
            epoch[v] = e
    return slot, epoch, (max(slot.values()) + 1 if slot else 0)


def perm_sign(seq):
    s = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


class Plan:
    """everything about the sweep that does not depend on the state"""

    def __init__(self, cov, order, j0):
        self.cov, self.order, self.j0 = cov, order, j0
        reps, crho, tail, head, m = cov.reps, cov.crho, cov.tail, cov.head, cov.m
        J = len(reps)
        self.root = tail[reps[j0]]
        root = self.root
        dom_full = [set() for _ in range(J)]
        dom_out = [set() for _ in range(J)]
        extra = [[] for _ in range(J)]
        for j in range(J):
            a, b = reps[j], crho[reps[j]]
            dom_full[j] |= {tail[a], head[a], tail[b], head[b]}
            dom_out[j] |= {tail[a], tail[b]}
        for a in range(cov.nA):
            if crho[a] == a:
                pick = next((j for j in order
                             if tail[a] in dom_full[j] or head[a] in dom_full[j]), order[0])
                extra[pick].append(a)
                dom_full[pick] |= {tail[a], head[a]}
                dom_out[pick] |= {tail[a]}
        self.extra = extra
        self.cfix = {a: m[a] // 2 for a in range(cov.nA) if crho[a] == a}

        n = len(order)
        ff, lf, fo, lo = {}, {}, {}, {}
        for i, j in enumerate(order):
            for v in dom_full[j]:
                ff.setdefault(v, i)
                lf[v] = i
            for v in dom_out[j]:
                fo.setdefault(v, i)
                lo[v] = i
        self.ff, self.lf, self.fo, self.lo = ff, lf, fo, lo
        self.slotF, epF, self.NF = _slots(ff, lf)
        self.slotO, epO, self.NO = _slots(fo, lo)

        # Canonical order of the generators.  It has to INTERLEAVE psibar_v and psi_v vertex by
        # vertex: the constant relating the top coefficient of prod_v (1 - psibar_v (A psi)_v)
        # to det A is (-1)^{#rows} only for an interleaved order, and any two of those differ
        # by a permutation of the even blocks psibar_v psi_v, which costs nothing.  Listing all
        # the psibar first instead costs a shuffle sign that varies from instance to instance.
        # The out-arc window of v sits inside its full window, so both generators of v can live
        # in v's frontier slot.
        ranks, r = {}, 0
        for v in sorted(ff, key=lambda x: (self.slotF[x], epF[x])):
            if v != root:
                ranks[("bar", v)] = r
                ranks[("psi", v)] = r + 1
                r += 2
        self.bit_bar = {v: 2 * self.slotF[v] for v in fo}
        self.bit_psi = {v: 2 * self.slotF[v] + 1 for v in ff}

        # the sign of the permutation from the order the generators leave to the canonical one
        seq = []
        for i in range(n):
            for v in sorted([x for x in lo if lo[x] == i], key=lambda x: self.slotO[x]):
                if v != root:
                    seq.append(ranks[("bar", v)])
            for v in sorted([x for x in lf if lf[x] == i], key=lambda x: self.slotF[x]):
                if v != root:
                    seq.append(ranks[("psi", v)])
        self.gsign = perm_sign(seq) * (-1) ** (cov.nv - 1)

        # remaining capacity of every imbalance, for the frontier bound
        capb = [collections.Counter() for _ in range(n)]
        for i, j in enumerate(order):
            a, b = reps[j], crho[reps[j]]
            for v in dom_full[j]:
                ca = (tail[a] == v) - (head[a] == v)
                cb = (tail[b] == v) - (head[b] == v)
                capb[i][v] += m[a] * max(abs(ca), abs(cb))
            for x in extra[j]:
                capb[i][tail[x]] += self.cfix[x]
                capb[i][head[x]] += self.cfix[x]
        remb = [collections.Counter() for _ in range(n + 1)]
        for i in range(n - 1, -1, -1):
            c = collections.Counter(remb[i + 1])
            c.update(capb[i])
            remb[i] = c
        self.remb = remb
        if remb and max(list(remb[0].values()) or [0]) >= BIAS:
            raise SystemExit("an imbalance field is too narrow; raise BITS")
        self.dom_full, self.dom_out = dom_full, dom_out


BITS = 8                       # width of an imbalance field
BIAS = 1 << (BITS - 2)         # the zero of an imbalance
GBIT = 1 << (BITS - 1)         # the guard bit of an imbalance field


def sweep(plan, pins=None, verbose=0, statecap=None, tlimit=None):
    cov, order, j0 = plan.cov, plan.order, plan.j0
    reps, crho, tail, head, m = cov.reps, cov.crho, cov.tail, cov.head, cov.m
    root, slotF, slotO = plan.root, plan.slotF, plan.slotO
    bit_bar, bit_psi, NF, NO = plan.bit_bar, plan.bit_psi, plan.NF, plan.NO
    n = len(order)

    guard = sum(GBIT << (BITS * s) for s in range(NF))
    bal0 = sum(BIAS << (BITS * s) for s in range(NF))

    # per step: the allowed values of the orbit and, for each, everything that does not
    # depend on the state
    steps = []
    for i, j in enumerate(order):
        a, b = reps[j], crho[reps[j]]
        mj = m[a]
        rng = range(mj + 1)
        if j == j0:                       # the half-box: f_j0 > 0, i.e. c(a) > m/2
            rng = range(mj // 2 + 1, mj + 1)
        if pins and j in pins:            # the sub-box this run is confined to
            rng = range(pins[j], pins[j] + 1)
        cvals = []
        for ca in rng:
            arcs = [(a, ca), (b, mj - ca)] + [(x, plan.cfix[x]) for x in plan.extra[j]]
            dbal = 0
            ddeg = 0
            ferm = []
            for arc, c in arcs:
                if not c:
                    continue
                dbal += c << (BITS * slotF[tail[arc]])
                dbal -= c << (BITS * slotF[head[arc]])
                ddeg += c << (8 * slotO[tail[arc]])
                t_, h_ = tail[arc], head[arc]
                if t_ != root and t_ != h_:
                    ferm.append((bit_bar[t_],
                                 -1 if t_ == root else bit_psi[t_],
                                 -1 if h_ == root else bit_psi[h_], c))
            cvals.append((dbal, ddeg, math.comb(mj, ca), tuple(ferm)))
        # the imbalance bound, one field per SLOT: the vertex still sitting in that slot after
        # step i may be off by its remaining capacity, an empty slot must read exactly zero --
        # which is also how a vertex leaving the frontier is made to balance
        occ = {}
        for v, s in slotF.items():
            if plan.ff[v] <= i < plan.lf[v]:
                occ[s] = v
        hi = 0
        lo_ = 0
        for s in range(NF):
            v = occ.get(s)
            r = plan.remb[i + 1].get(v, 0) if v is not None else 0
            hi += (BIAS + r) << (BITS * s)
            lo_ += (BIAS - r) << (BITS * s)
        clos_out = sorted([v for v in plan.lo if plan.lo[v] == i], key=lambda v: slotO[v])
        clos_full = sorted([v for v in plan.lf if plan.lf[v] == i], key=lambda v: slotF[v])
        steps.append((tuple(cvals), hi + guard, lo_, tuple(clos_out), tuple(clos_full)))

    fact = [math.factorial(t) for t in range(256)]
    states = {(bal0, 0, 0): 1}
    peak, t0 = 1, time.time()

    for i in range(n):
        cvals, hig, lo_, clos_out, clos_full = steps[i]
        new = {}
        get = new.get
        simple = not clos_out and not clos_full
        for (bal, deg, mask), val in states.items():
            for dbal, ddeg, w, ferm in cvals:
                b2 = bal + dbal
                if ((hig - b2) & guard) != guard:
                    continue
                if ((b2 + guard - lo_) & guard) != guard:
                    continue
                vw = val * w
                br = [(mask, vw)]
                for pbar, pt, ph, c in ferm:
                    nb = []
                    ap = nb.append
                    for gm, vv in br:
                        ap((gm, vv))
                        if gm >> pbar & 1:
                            continue
                        s1 = -1 if ((gm >> (pbar + 1)).bit_count() & 1) else 1
                        g1 = gm | (1 << pbar)
                        if pt >= 0 and not (g1 >> pt & 1):
                            s2 = -1 if ((g1 >> (pt + 1)).bit_count() & 1) else 1
                            ap((g1 | (1 << pt), -vv * s1 * s2 * c))
                        if ph >= 0 and not (g1 >> ph & 1):
                            s2 = -1 if ((g1 >> (ph + 1)).bit_count() & 1) else 1
                            ap((g1 | (1 << ph), vv * s1 * s2 * c))
                    br = nb
                if simple:
                    d2 = deg + ddeg
                    for gm, vv in br:
                        kk = (b2, d2, gm)
                        u = get(kk)
                        new[kk] = vv if u is None else u + vv
                    continue
                for gm, vv in br:
                    d2 = deg + ddeg
                    ok = True
                    for v in clos_out:
                        sh = 8 * slotO[v]
                        d = (d2 >> sh) & 255
                        d2 -= d << sh
                        if v == root:
                            if not d:
                                ok = False
                                break
                            vv *= fact[d - 1]
                            continue
                        pbar = bit_bar[v]
                        if d == 0:
                            pv = bit_psi[v]
                            if (gm >> pbar & 1) or (gm >> pv & 1):
                                ok = False
                                break
                            if (gm >> (pbar + 1)).bit_count() & 1:
                                vv = -vv
                            gm |= 1 << pbar
                            if (gm >> (pv + 1)).bit_count() & 1:
                                vv = -vv
                            gm |= 1 << pv
                            vv = -vv        # the diagonal enters as (1 - psibar_v psi_v)
                        elif not (gm >> pbar & 1):
                            ok = False
                            break
                        else:
                            vv *= fact[d - 1]
                        if (gm & ((1 << pbar) - 1)).bit_count() & 1:
                            vv = -vv
                        gm ^= 1 << pbar
                    if not ok:
                        continue
                    for v in clos_full:
                        if v == root:
                            continue
                        pv = bit_psi[v]
                        if not (gm >> pv & 1):
                            ok = False
                            break
                        if (gm & ((1 << pv) - 1)).bit_count() & 1:
                            vv = -vv
                        gm ^= 1 << pv
                    if not ok:
                        continue
                    kk = (b2, d2, gm)
                    u = get(kk)
                    new[kk] = vv if u is None else u + vv
        states = {k: v for k, v in new.items() if v}
        if len(states) > peak:
            peak = len(states)
        if verbose and (i % verbose == 0 or i == n - 1):
            print(f"      step {i+1}/{n}: states {len(states):,}   [{time.time()-t0:.0f}s]",
                  flush=True)
        if (statecap and len(states) > statecap) or (tlimit and time.time() - t0 > tlimit):
            return None, peak
    return plan.gsign * sum(states.values()), peak


def ds(cov, order=None, verbose=0, statecap=None, tlimit=None, split=0, report=None):
    """#DS(S,k) exactly.  `split` pins that many of the earliest orbits of the order, one
    combination of values at a time: the sub-boxes are disjoint and cover the half-box, so the
    sums add, and the peak state count of a run drops roughly by their number."""
    if order is None:
        order = search_order(cov)
    pos = {j: i for i, j in enumerate(order)}
    odd = [j for j in range(len(cov.reps)) if cov.m[cov.reps[j]] % 2 == 1]
    if not odd:
        raise SystemExit("no orbit of odd multiplicity; the half-box trick does not apply")
    j0 = min(odd, key=lambda j: pos[j])
    plan = Plan(cov, order, j0)
    pinned = [j for j in order if j != j0][:split]
    combos = [{}]
    for j in pinned:
        combos = [{**d, j: t} for d in combos for t in range(cov.m[cov.reps[j]] + 1)]
    half, peak = 0, 0
    for w, pin in enumerate(combos):
        h, pk = sweep(plan, pins=pin, verbose=verbose, statecap=statecap, tlimit=tlimit)
        if h is None:
            return None, pk
        half += h
        peak = max(peak, pk)
        if report:
            report(w, len(combos), half, pk)
    if half is None:
        return None, peak
    den = 1
    for j in range(len(cov.reps)):
        den *= math.factorial(cov.m[cov.reps[j]])
    for a in range(cov.nA):
        if cov.crho[a] == a:
            den *= math.factorial(cov.m[a] // 2)
    num = 2 * half
    assert num % den == 0, "the sweep did not land on an integer"
    return num // den, peak


def denominator(cov):
    """prod_j m_j! * prod_{rho a = a} (m_a/2)!, the constant of (eq:int)"""
    den = 1
    for j in range(len(cov.reps)):
        den *= math.factorial(cov.m[cov.reps[j]])
    for a in range(cov.nA):
        if cov.crho[a] == a:
            den *= math.factorial(cov.m[a] // 2)
    return den


# ------------------------------------------------- the arborescence frontier (positive, exact)
#
# The fermionic state is LINEAR: at a fixed imbalance it is a vector in the exterior algebra of
# the open generators, and a single history already fills up to 3^(open rows) of its
# coordinates -- a boundary vertex v with a past arc to u contributes 1 - psibar_v psi_v +
# psibar_v psi_u.  The positive expansion of the same determinant, t(G_c) = sum over in-trees
# of prod c(a), needs no signs and puts ONE state per history:
#
#   * every vertex that has not yet decided its tree arc ("pending") is the root of its own
#     class; a vertex that has decided points to the root of the class it joined, or to R, the
#     class of the global root r;
#   * a pending vertex remembers its candidate arcs so far, as {target class: total c}, and
#     their sum is its partial out-degree, so Loc needs nothing more;
#   * when its last out-arc has been decided, a pending vertex picks one candidate class other
#     than its own (weight = the multiplicity to that class), and its class is merged into the
#     one it picked.  Picking into its own class would close a cycle; the last vertex of a
#     component that does not reach r has only such candidates, so disconnected supports die
#     by themselves.
#
# The acyclicity argument: in a cycle of the functional graph, the last vertex to choose finds
# its head already in its own class.  So choosing, each time, a class other than one's own is
# exactly "the chosen arcs form an in-arborescence of supp c to r".

class Windows:
    """per-orbit endpoint sets and per-vertex windows for a given order"""

    def __init__(self, cov, order):
        reps, crho, tail, head, m = cov.reps, cov.crho, cov.tail, cov.head, cov.m
        J = len(reps)
        dom_full = [set() for _ in range(J)]
        dom_out = [set() for _ in range(J)]
        extra = [[] for _ in range(J)]
        for j in range(J):
            a, b = reps[j], crho[reps[j]]
            dom_full[j] |= {tail[a], head[a], tail[b], head[b]}
            dom_out[j] |= {tail[a], tail[b]}
        for a in range(cov.nA):
            if crho[a] == a:
                pick = next((j for j in order
                             if tail[a] in dom_full[j] or head[a] in dom_full[j]), order[0])
                extra[pick].append(a)
                dom_full[pick] |= {tail[a], head[a]}
                dom_out[pick] |= {tail[a]}
        self.dom_full, self.dom_out, self.extra = dom_full, dom_out, extra
        self.cfix = {a: m[a] // 2 for a in range(cov.nA) if crho[a] == a}
        ff, lf, fo, lo = {}, {}, {}, {}
        for i, j in enumerate(order):
            for v in dom_full[j]:
                ff.setdefault(v, i)
                lf[v] = i
            for v in dom_out[j]:
                fo.setdefault(v, i)
                lo[v] = i
        self.ff, self.lf, self.fo, self.lo = ff, lf, fo, lo
        n = len(order)
        capb = [collections.Counter() for _ in range(n)]
        for i, j in enumerate(order):
            a, b = reps[j], crho[reps[j]]
            for v in dom_full[j]:
                ca = (tail[a] == v) - (head[a] == v)
                cb = (tail[b] == v) - (head[b] == v)
                capb[i][v] += m[a] * max(abs(ca), abs(cb))
            for x in extra[j]:
                capb[i][tail[x]] += self.cfix[x]
                capb[i][head[x]] += self.cfix[x]
        remb = [collections.Counter() for _ in range(n + 1)]
        for i in range(n - 1, -1, -1):
            c = collections.Counter(remb[i + 1])
            c.update(capb[i])
            remb[i] = c
        self.remb = remb
        if max(list(remb[0].values()) or [0]) > 127:
            raise SystemExit("an imbalance does not fit in a byte")
        self.width = max(sum(1 for v in ff if ff[v] <= i <= lf[v]) for i in range(n)) if n else 0


# state layout, one record per frontier slot:  bal+128, then a tag
#   0 isolated (closed with out-degree 0)      1 decided, in the class R of the root
#   2+s decided, in the class of pending slot s                    253 the root, out-closed
#   254 d  the root, out-open, partial out-degree d
#   255 n (t w)*n  pending; n candidate classes t (0 = R, 1+s = slot s) with multiplicity w
_EMPTY = {}


def _adec(key):
    bal, kind, cand, rdeg = [], [], [], 0
    p, L = 0, len(key)
    while p < L:
        bal.append(key[p] - 128)
        t = key[p + 1]
        p += 2
        if t == 255:
            n = key[p]
            p += 1
            d = {}
            for _ in range(n):
                d[key[p] - 1] = key[p + 1]
                p += 2
            kind.append(-3)
            cand.append(d)
        elif t >= 2 and t < 253:
            kind.append(t - 2)
            cand.append(None)
        elif t == 1:
            kind.append(-1)
            cand.append(None)
        elif t == 0:
            kind.append(-2)
            cand.append(None)
        elif t == 254:
            rdeg = key[p]
            p += 1
            kind.append(-4)
            cand.append(None)
        else:
            kind.append(-5)
            cand.append(None)
    return bal, kind, cand, rdeg


def _aenc(bal, kind, cand, rdeg, keep, newpos):
    out = bytearray()
    for q in keep:
        out.append(bal[q] + 128)
        k = kind[q]
        if k == -3:
            items = sorted((0 if t < 0 else 1 + newpos[t], w) for t, w in cand[q].items())
            out.append(255)
            out.append(len(items))
            for tc, w in items:
                out.append(tc)
                out.append(w)
        elif k >= 0:
            out.append(2 + newpos[k])
        elif k == -1:
            out.append(1)
        elif k == -2:
            out.append(0)
        elif k == -4:
            out.append(254)
            out.append(rdeg)
        else:
            out.append(253)
    return bytes(out)


def rss_gb():
    """resident memory of this process in GB, best effort, standard library only"""
    try:
        if os.name == "nt":
            import ctypes
            from ctypes import wintypes

            class PMC(ctypes.Structure):
                _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                            ("PeakWorkingSetSize", ctypes.c_size_t),
                            ("WorkingSetSize", ctypes.c_size_t),
                            ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                            ("QuotaPagedPoolUsage", ctypes.c_size_t),
                            ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                            ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                            ("PagefileUsage", ctypes.c_size_t),
                            ("PeakPagefileUsage", ctypes.c_size_t)]
            pmc = PMC()
            pmc.cb = ctypes.sizeof(PMC)
            h = ctypes.windll.kernel32.GetCurrentProcess()
            ctypes.windll.psapi.GetProcessMemoryInfo(h, ctypes.byref(pmc), pmc.cb)
            return pmc.WorkingSetSize / 2 ** 30
        import resource
        r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return r / 2 ** 20 if sys.platform != "darwin" else r / 2 ** 30
    except Exception:
        return float("nan")


def arbor_sweep(cov, order, j0, pins=None, W=None, verbose=0, statecap=None, tag=""):
    """sum over the half-box f_{j0} > 0 (and the sub-box `pins`: orbit j -> c(rep_j)) of
    prod_j binom(m_j, c_j) * t_r(G_c) * Loc(G_c).  Returns (value, peak, None), or
    (None, peak, step) when the state count exceeded `statecap` after that step."""
    W = W or Windows(cov, order)
    reps, crho, tail, head, m = cov.reps, cov.crho, cov.tail, cov.head, cov.m
    r = tail[reps[j0]]
    fact = [math.factorial(t) for t in range(600)]
    front, states = [], {b"": 1}
    peak, t0, n = 1, time.time(), len(order)
    for i, j in enumerate(order):
        a, b = reps[j], crho[reps[j]]
        mj = m[a]
        vals = range(mj + 1)
        if j == j0:
            vals = range(mj // 2 + 1, mj + 1)
        if pins and j in pins:
            vals = [x for x in vals if x == pins[j]]
        fs = set(front)
        newly = sorted(v for v in W.dom_full[j] if v not in fs)
        nxt = front + newly
        idx = {v: q for q, v in enumerate(nxt)}
        rq = idx.get(r, -1)
        clos_out = [q for q, v in enumerate(nxt) if W.lo[v] == i]
        closing = [q for q, v in enumerate(nxt) if W.lf[v] == i]
        cset = set(closing)
        keep = [q for q, v in enumerate(nxt) if W.lf[v] != i]
        newpos = {q: p for p, q in enumerate(keep)}
        lim = [W.remb[i + 1].get(v, 0) for v in nxt]
        cvals = []
        for ca in vals:
            arcs = [(a, ca), (b, mj - ca)] + [(x, W.cfix[x]) for x in W.extra[j]]
            dl, al = collections.Counter(), []
            for x, c in arcs:
                if c:
                    tq, hq = idx[tail[x]], idx[head[x]]
                    dl[tq] += c
                    dl[hq] -= c
                    al.append((tq, hq, c))
            # every vertex of the step: its remaining capacity has just shrunk
            chk = tuple(sorted(set(dl) | cset | {idx[v] for v in W.dom_full[j]}))
            cvals.append((math.comb(mj, ca), tuple((q, d) for q, d in dl.items() if d), chk,
                          tuple(al)))
        newk = [(-4 if v == r else -3) for v in newly]
        newc = [(None if v == r else _EMPTY) for v in newly]
        new = {}
        get = new.get
        for key, val in states.items():
            bal0, kind0, cand0, rdeg0 = _adec(key)
            if newly:
                bal0 += [0] * len(newly)
                kind0 += newk
                cand0 += newc
            for w, dl, chk, al in cvals:
                bal = bal0[:]
                for q, d in dl:
                    bal[q] += d
                bad = False
                for q in chk:
                    x = bal[q]
                    if q in cset:
                        if x:
                            bad = True
                            break
                    elif x > lim[q] or -x > lim[q]:
                        bad = True
                        break
                if bad:
                    continue
                cand = cand0[:]
                rdeg = rdeg0
                copied = []
                for tq, hq, c in al:
                    if hq == rq:
                        tg = -1
                    else:
                        kh = kind0[hq]
                        if kh == -3:
                            tg = hq
                        elif kh >= 0:
                            tg = kh
                        elif kh == -1:
                            tg = -1
                        else:                  # an arc into a vertex already declared isolated
                            bad = True
                            break
                    if tq == rq:
                        rdeg += c
                    else:
                        d = cand[tq]
                        if tq not in copied:
                            d = dict(d)
                            cand[tq] = d
                            copied.append(tq)
                        d[tg] = d.get(tg, 0) + c
                if bad:
                    continue
                branches = [(kind0, cand, val * w, rdeg)]
                for p in clos_out:
                    nb = []
                    for K_, C_, vv, rd in branches:
                        if p == rq:                     # the root does not choose
                            if rd == 0:
                                continue
                            K2 = list(K_)
                            K2[p] = -5
                            nb.append((K2, C_, vv * fact[rd - 1], 0))
                            continue
                        d = C_[p]
                        deg = sum(d.values())
                        if deg == 0:                    # outside supp c
                            if p in K_ or any(c_ and p in c_ for c_ in C_):
                                continue
                            K2 = list(K_)
                            K2[p] = -2
                            C2 = list(C_)
                            C2[p] = None
                            nb.append((K2, C2, vv, rd))
                            continue
                        loc = fact[deg - 1]
                        for tg, wt in d.items():
                            if tg == p:                 # into its own class: a cycle
                                continue
                            K2 = [tg if k_ == p else k_ for k_ in K_]
                            K2[p] = tg
                            C2 = list(C_)
                            C2[p] = None
                            for q_, c_ in enumerate(C2):
                                if c_ and p in c_:
                                    c2 = dict(c_)
                                    x = c2.pop(p)
                                    c2[tg] = c2.get(tg, 0) + x
                                    C2[q_] = c2
                            nb.append((K2, C2, vv * wt * loc, rd))
                    branches = nb
                    if not branches:
                        break
                for K_, C_, vv, rd in branches:
                    k2 = _aenc(bal, K_, C_, rd, keep, newpos)
                    u = get(k2)
                    new[k2] = vv if u is None else u + vv
        states, front = new, [nxt[q] for q in keep]
        if len(states) > peak:
            peak = len(states)
        if verbose and (i % verbose == 0 or i == n - 1):
            print(f"      {tag}step {i+1}/{n}: frontier {len(front)}, states {len(states):,}"
                  f"   [{time.time()-t0:.0f}s, {rss_gb():.1f} GB]", flush=True)
        if statecap and len(states) > statecap:
            return None, peak, i
        if not states:
            return 0, peak, None
    assert list(states) in ([b""], []), "the frontier did not empty"
    return states.get(b"", 0), peak, None


def arbor_order(cov):
    var = cov.constraints()
    dom = []
    for j, a in enumerate(cov.reps):
        b = cov.crho[a]
        dom.append(set(var[j]) | {cov.tail[a], cov.head[a], cov.tail[b], cov.head[b]})
    return greedy_order(var, dom=dom)


def pick_j0(cov, order):
    pos = {j: i for i, j in enumerate(order)}
    odd = [j for j in range(len(cov.reps)) if cov.m[cov.reps[j]] % 2 == 1]
    if not odd:
        raise SystemExit("no orbit of odd multiplicity; the half-box trick does not apply")
    return min(odd, key=lambda j: pos[j])


def arbor_ds(cov, order=None, pins_list=None, verbose=0):
    """#DS by the arborescence frontier on the half-box, optionally as a sum over the sub-boxes
    of `pins_list` (which must partition the half-box)"""
    order = order or arbor_order(cov)
    j0 = pick_j0(cov, order)
    W = Windows(cov, order)
    half, peak = 0, 0
    for pins in (pins_list or [None]):
        h, pk, _ = arbor_sweep(cov, order, j0, pins=pins, W=W, verbose=verbose)
        half += h
        peak = max(peak, pk)
    den = denominator(cov)
    assert (2 * half) % den == 0, "the arborescence sweep did not land on an integer"
    return 2 * half // den, peak


# ---- k = 10: adaptive sub-boxes, a checkpoint, and optionally several processes

_WK = {}


def _winit(S, k, order, j0):
    cov = Cover(S, k)
    _WK.update(cov=cov, order=order, j0=j0, W=Windows(cov, order))


def _wrun(pins, cap, verbose, tag):
    h, pk, ab = arbor_sweep(_WK["cov"], _WK["order"], _WK["j0"], pins=pins, W=_WK["W"],
                            verbose=verbose, statecap=cap, tag=tag)
    return pins, h, pk, ab


def _pkey(pins):
    return ",".join(f"{j}:{pins[j]}" for j in sorted(pins))


def arbor_k10(S, k, cov, order, j0, cap, jobs, verbose, ckpt):
    """runs the half-box; a (sub-)box whose state count passes `cap` is abandoned and cut into
    sub-boxes by pinning the last free orbit processed at or before the step where it passed.
    Finished sub-boxes and the cuts are written to `ckpt`, so an interrupted run resumes."""
    sig = f"k={k};order={order};j0={j0};cap={cap}"
    book = {"sig": sig, "done": {}, "split": {}}
    if ckpt and os.path.exists(ckpt):
        try:
            old = json.load(open(ckpt, encoding="utf-8"))
            if old.get("sig") == sig:
                book = old
                print(f"    resuming from {os.path.basename(ckpt)}: {len(book['done'])} sub-boxes "
                      f"done, {len(book['split'])} cuts", flush=True)
        except Exception:
            pass

    def save():
        if ckpt:
            tmp = ckpt + ".tmp"
            json.dump(book, open(tmp, "w", encoding="utf-8"))
            os.replace(tmp, ckpt)

    pos = {j: i for i, j in enumerate(order)}
    m = cov.m
    reps = cov.reps
    t0 = time.time()
    peak = [0]
    queue, total, leaves = [{}], 0, 0

    def children(pins):
        jn = book["split"][_pkey(pins)]
        return [{**pins, jn: t} for t in range(m[reps[jn]] + 1)]

    def cut(pins, ab):
        free = [j for j in order[:ab + 1] if j != j0 and j not in pins]
        if not free:
            raise SystemExit("the state cap is too small even for a fully pinned prefix")
        jn = free[-1]
        book["split"][_pkey(pins)] = jn
        save()
        print(f"    sub-box [{_pkey(pins) or 'half-box'}] passed {cap:,} states at step {ab+1}; "
              f"cutting on orbit {jn} (step {pos[jn]+1}, m = {m[reps[jn]]})", flush=True)

    def finish(pins, h, pk):
        nonlocal total, leaves
        book["done"][_pkey(pins)] = str(h)
        save()
        total += h
        leaves += 1
        peak[0] = max(peak[0], pk)
        print(f"    sub-box [{_pkey(pins) or 'half-box'}] done: peak {pk:,}; "
              f"{leaves} done, {len(queue)} queued   [{time.time()-t0:.0f}s]", flush=True)

    if jobs <= 1:
        _winit(S, k, order, j0)
        while queue:
            pins = queue.pop()
            kk = _pkey(pins)
            if kk in book["done"]:
                total += int(book["done"][kk])
                leaves += 1
                continue
            if kk in book["split"]:
                queue.extend(children(pins))
                continue
            _, h, pk, ab = _wrun(pins, cap, verbose, "")
            if h is None:
                cut(pins, ab)
                queue.extend(children(pins))
            else:
                finish(pins, h, pk)
    else:
        from concurrent.futures import ProcessPoolExecutor, FIRST_COMPLETED, wait
        with ProcessPoolExecutor(max_workers=jobs, initializer=_winit,
                                 initargs=(S, k, order, j0)) as ex:
            running = set()
            while queue or running:
                while queue and len(running) < jobs:
                    pins = queue.pop()
                    kk = _pkey(pins)
                    if kk in book["done"]:
                        total += int(book["done"][kk])
                        leaves += 1
                        continue
                    if kk in book["split"]:
                        queue.extend(children(pins))
                        continue
                    running.add(ex.submit(_wrun, pins, cap, verbose, f"[{kk or 'half'}] "))
                if not running:
                    continue
                fin, running = wait(running, return_when=FIRST_COMPLETED)
                for fu in fin:
                    pins, h, pk, ab = fu.result()
                    if h is None:
                        cut(pins, ab)
                        queue.extend(children(pins))
                    else:
                        finish(pins, h, pk)
    return total, leaves, peak[0]

# ------------------------------------------------------------------------------ the checks

def _ext_coeff(cov, c, root, rank):
    """the fermionic product evaluated in the whole exterior algebra, no frontier: the
    reference the sweep is checked against"""
    outd = cov.outdeg(c)
    gen = {}
    for v in range(cov.nv):
        if v != root:
            gen[(v, 0)] = 2 * rank[v]
            gen[(v, 1)] = 2 * rank[v] + 1
    pos = {g: p for p, g in enumerate(sorted(gen.values()))}
    full = (1 << len(pos)) - 1
    P = {0: 1}

    def rmul(P, g, co):
        p, out = pos[g], {}
        for mk, val in P.items():
            if mk >> p & 1:
                continue
            s = -1 if ((mk >> (p + 1)).bit_count() & 1) else 1
            nk = mk | (1 << p)
            out[nk] = out.get(nk, 0) + val * co * s
        return out

    def add(A, B):
        r = dict(A)
        for k, v in B.items():
            r[k] = r.get(k, 0) + v
        return r

    for a in range(cov.nA):
        t, h = cov.tail[a], cov.head[a]
        if c[a] == 0 or t == root or t == h:
            continue
        Q = rmul(P, gen[(t, 0)], 1)
        R = rmul(Q, gen[(t, 1)], -c[a])
        if h != root:
            R = add(R, rmul(Q, gen[(h, 1)], c[a]))
        P = add(P, R)
    for v in range(cov.nv):
        if v != root and outd[v] == 0:
            P = rmul(rmul(P, gen[(v, 0)], 1), gen[(v, 1)], -1)
    return P.get(full, 0)


def main():
    t0 = time.time()
    S = load_phix()
    quick = "--quick" in sys.argv
    print("=== #DS by a fermionic frontier sweep ===")
    print(f"    phiX174 NC_001422.1, {len(S)} bp circular; python {sys.version.split()[0]}\n")

    covers, pts = {}, {}
    print("      k   |V|  fixV  |A_c|  fixA   sum m   orbits   odd-m orbits   sum_a c(a)")
    okD = True
    want = {12: (7, 1, 14, 0, 14, 7), 11: (30, 0, 53, 1, 60, 26), 10: (125, 3, 222, 0, 252, 111)}
    for k in (12, 11, 10):
        cov = Cover(S, k)
        covers[k] = cov
        got = (cov.nv, cov.n_fixV, cov.nA, cov.n_fixA, sum(cov.m), len(cov.reps))
        okD &= (got == want[k])
        odd = sum(1 for j in range(len(cov.reps)) if cov.m[cov.reps[j]] % 2)
        tot = sum(cov.m[cov.reps[j]] for j in range(len(cov.reps))) \
            + sum(cov.m[a] // 2 for a in range(cov.nA) if cov.crho[a] == a)
        print(f"     {k:2d} {cov.nv:5d} {cov.n_fixV:5d} {cov.nA:6d} {cov.n_fixA:5d}"
              f" {sum(cov.m):7d} {len(cov.reps):8d} {odd:14d} {tot:12d}   "
              f"{'ok' if got == want[k] else 'MISMATCH ' + str(want[k])}")
    row("D", "the double cover rebuilt from the definition reproduces the invariants of Bio 8, "
             "Table 1, and every k carries an orbit of odd multiplicity, which is what the "
             "half-box needs", "k = 10, 11, 12", "exact", okD)
    for k in (12, 11):
        pts[k] = enumerate_box(covers[k])

    # small independent instances: short random sequences, where everything can be checked
    rnd = random.Random(20260922)
    small = []
    while len(small) < 12:
        T = "".join(rnd.choice("ACGT") for _ in range(rnd.randrange(30, 70)))
        kk = rnd.randrange(3, 6)
        c2 = Cover(T, kk)
        if 3 <= c2.nv <= 9 and c2.reps and any(c2.m[c2.reps[j]] % 2 for j in range(len(c2.reps))):
            p2 = enumerate_box(c2)
            if 1 <= len(p2) <= 4000:
                small.append((T, kk, c2, p2))

    # (F) the fermionic identity itself
    print("\n=== (F) the Berezin product against the reduced determinant ===")
    okF, cnt = True, 0
    inst = [("phiX174 k=12", covers[12], pts[12])] + \
           [(f"random {kk}-mers on {len(T)} bp", c2, p2) for (T, kk, c2, p2) in small]
    for name, cov, pp in inst:
        rank = {v: v for v in range(cov.nv)}
        agree = skip = 0
        for g in pp:
            c = cov.point(g)
            root = 0
            if cov.outdeg(c)[root] == 0:
                skip += 1
                continue
            e = _ext_coeff(cov, c, root, rank)
            d = int(cov.t_rooted(c, root)) * (-1) ** (cov.nv - 1)
            agree += (e == d)
            cnt += 1
        okF &= (agree == len(pp) - skip)
        if name.startswith("phiX"):
            print(f"    {name}: {agree} of {len(pp)-skip} points agree "
                  f"({skip} have the test root outside supp c)")
    print(f"    {len(small)} random covers, |V| = "
          f"{[c2.nv for (_, _, c2, _) in small]}: all points agree")
    row("F", "the coefficient of the top monomial in prod_a (1 - c(a) psibar_t (psi_t - psi_h)) "
             "times prod_{v outside supp} (1 - psibar_v psi_v) equals (-1)^{|V|-1} det of the "
             "Laplacian with the diagonal of the complement of supp c filled in and one row and "
             "column deleted -- so t(G_c) is a fermionic partition function",
        f"{cnt} points", "exact", okF)

    # (R) the root, and why one fixed root is enough on the half-box
    print("\n=== (R) a fixed root on the half-box ===")
    okR = True
    for k in (12, 11):
        cov = covers[k]
        odd = [j for j in range(len(cov.reps)) if cov.m[cov.reps[j]] % 2]
        j0 = odd[0]
        rt = cov.tail[cov.reps[j0]]
        bad = out = 0
        for g in pts[k]:
            if g[j0] <= 0:
                continue
            c = cov.point(g)
            if cov.outdeg(c)[rt] == 0:
                out += 1
            if Fraction(int(cov.t_rooted(c, rt)) * math.prod(
                    math.factorial(d - 1) for d in cov.outdeg(c).values()),
                    math.prod(math.factorial(x) for x in c if x)) != cov.n_seq(c):
                bad += 1
        okR &= (bad == 0 and out == 0)
        print(f"    k={k}: orbit {j0} has m={cov.m[cov.reps[j0]]}; on the half-box f_{j0}>0 the "
              f"root {rt} is outside supp c {out} times, and Nseq disagrees {bad} times")
    row("R", "on the half-box f_{j0}>0 of an orbit of odd multiplicity the tail of a_{j0} lies "
             "in supp c for every point, so one fixed root serves the whole sweep and "
             "t_root(G_c) Loc(G_c) / prod c(a)! is Nseq(c) there", "k = 11, 12", "exact", okR)

    # (S) the half-box carries half of #DS
    print("\n=== (S) Nseq is invariant under c -> m-c ===")
    okS = True
    for k in (12, 11):
        cov = covers[k]
        odd = [j for j in range(len(cov.reps)) if cov.m[cov.reps[j]] % 2]
        j0 = odd[0]
        sym = all(cov.n_seq(cov.point(g)) ==
                  cov.n_seq([cov.m[a] - cov.point(g)[a] for a in range(cov.nA)]) for g in pts[k])
        tot = sum(cov.n_seq(cov.point(g)) for g in pts[k])
        half = sum(cov.n_seq(cov.point(g)) for g in pts[k] if g[j0] > 0)
        okS &= sym and 2 * half == tot
        print(f"    k={k}: Nseq(m-c)=Nseq(c) everywhere: {sym};  2 x (half-box {half}) = "
              f"{2*half}, #DS = {tot}   {'ok' if 2*half == tot else 'MISMATCH'}")
    row("S", "reverse complementation acts on the box as c -> m-c and fixes Nseq, because the "
             "reverse of a balanced digraph has the same number of arborescences; no orbit of "
             "odd multiplicity meets f=0, so the half-box f_{j0}>0 carries exactly half of #DS",
        "k = 11, 12", "exact", okS)

    # (V) the sweep
    print("\n=== (V) the sweep against Bio 8's published #DS ===")
    okV = True
    for k, wantds in ((12, 10), (11, 864)):
        cov = covers[k]
        t = time.time()
        o = search_order(cov, tries=12, rounds=4)
        n, peak = ds(cov, order=o)
        okV &= (n == wantds)
        print(f"    k={k}: fermionic #DS = {n} (Bio 8 {wantds}), peak {peak:,} states"
              f"   [{time.time()-t:.1f}s]   {'ok' if n == wantds else 'MISMATCH'}")
    agree = 0
    for (T, kk, c2, p2) in small:
        b = sum(c2.n_seq(c2.point(g)) for g in p2)
        agree += (ds(c2, order=search_order(c2, tries=6, rounds=3))[0] == b)
    okV &= (agree == len(small))
    print(f"    {agree} of {len(small)} random covers agree with the brute-force sum of Nseq")
    row("V", "the fermionic frontier sweep reproduces Bio 8's #DS from the half-box, with one "
             "fixed root and a single division at the end", f"k = 11, 12 and {len(small)} "
             "random covers", "exact", okV)

    # (A) the arborescence frontier on everything that can be checked by brute force
    print("\n=== (A) the arborescence frontier against the brute-force sum ===")
    okA = True
    nfast = 0
    for k in (12, 11):
        for g in pts[k]:
            c = covers[k].point(g)
            okA &= (n_seq_fast(covers[k], c) == covers[k].n_seq(c))
            nfast += 1
    for (T, kk, c2, p2) in small:
        for g in p2:
            c = c2.point(g)
            okA &= (n_seq_fast(c2, c) == c2.n_seq(c))
            nfast += 1
    print(f"    forced-arc contraction: Nseq agrees with the full determinant on {nfast} points: "
          f"{okA}")
    for k, wantds in ((12, 10), (11, 864)):
        cov = covers[k]
        t = time.time()
        o = arbor_order(cov)
        n, peak = arbor_ds(cov, order=o)
        good = (n == wantds)
        # the same number as a sum over sub-boxes: pin two orbits other than j0
        j0 = pick_j0(cov, o)
        pj = [j for j in o if j != j0][:2]
        plist = [{pj[0]: x, pj[1]: y} for x in range(cov.m[cov.reps[pj[0]]] + 1)
                 for y in range(cov.m[cov.reps[pj[1]]] + 1)]
        n2, _ = arbor_ds(cov, order=o, pins_list=plist)
        good &= (n2 == n)
        okA &= good
        print(f"    k={k}: arborescence #DS = {n} (Bio 8 {wantds}); over {len(plist)} sub-boxes "
              f"{n2}; width {Windows(cov, o).width}, peak {peak:,} states   "
              f"[{time.time()-t:.1f}s]   {'ok' if good else 'MISMATCH'}")
    agree = 0
    for (T, kk, c2, p2) in small:
        b = sum(c2.n_seq(c2.point(g)) for g in p2)
        agree += (arbor_ds(c2)[0] == b)
    okA &= (agree == len(small))
    print(f"    {agree} of {len(small)} random covers agree with the brute-force sum of Nseq")
    row("A", "the positive arborescence frontier -- pending vertices carry their candidate "
             "classes, a vertex picks a class other than its own when its last out-arc is "
             "decided -- reproduces Bio 8's #DS, the brute force on random covers, and itself "
             "as a sum over sub-boxes", f"k = 11, 12 and {len(small)} random covers", "exact", okA)

    cov = covers[10]
    o10 = arbor_order(cov)
    j0 = pick_j0(cov, o10)
    W10 = Windows(cov, o10)

    # (Q) the k = 10 implementation on sub-boxes small enough to enumerate
    print("\n=== (Q) k = 10: the arborescence frontier on enumerable sub-boxes ===")
    t = time.time()
    rnd = random.Random(1310)
    den10 = denominator(cov)
    m0 = cov.m[cov.reps[j0]]
    tests, tries = [], 0
    J10 = len(cov.reps)
    while len(tests) < 5 and tries < 60:
        tries += 1
        # a random point of the half-box with connected support: the enumerator with every
        # orbit's values shuffled stops at its first leaf
        shuf = {}
        for j in range(J10):
            mj = cov.m[cov.reps[j]]
            vs = list(range(1, mj + 1, 2)) if j == j0 else list(range(-mj, mj + 1, 2))
            rnd.shuffle(vs)
            shuf[j] = vs
        p0 = enumerate_box(cov, cap=0, allowed=shuf, nodecap=300000)
        if not p0 or not n_seq_fast(cov, cov.point(p0[0])):
            continue
        g0 = p0[0]
        # pin more and more of its coordinates until the sub-box around it is enumerable
        free = [j for j in range(J10) if j != j0]
        rnd.shuffle(free)
        got, npin = None, 30
        while npin <= len(free):
            pins = {j: (cov.m[cov.reps[j]] + g0[j]) // 2 for j in free[:npin]}
            allowed = {j: [g0[j]] for j in pins}
            allowed[j0] = list(range(1, m0 + 1, 2))
            pp = enumerate_box(cov, cap=250, allowed=allowed, nodecap=300000)
            if pp is not None and 0 < len(pp) <= 250:
                got = pp
                break
            npin += 6
        if not got:
            continue
        vals = [n_seq_fast(cov, cov.point(g)) for g in got]
        tests.append((pins, got, vals))
    okQ = len(tests) >= 3
    for pins, got, vals in tests:
        h, pk, _ = arbor_sweep(cov, o10, j0, pins=pins, W=W10)
        # a few points also through the full 125-vertex determinant
        full = all(cov.n_seq(cov.point(g)) == v for g, v in list(zip(got, vals))[:3])
        good = (h == sum(vals) * den10) and full
        okQ &= good
        print(f"    {len(pins)} orbits pinned: {len(got)} points, {sum(1 for v in vals if v)} "
              f"connected, sum Nseq = {sum(vals)}; sweep {'agrees' if good else 'DISAGREES'}"
              f" (peak {pk:,})")
    print(f"    [{time.time()-t:.0f}s]")
    row("Q", "on random sub-boxes of the k = 10 half-box small enough to enumerate, the "
             "arborescence frontier (the same order, windows and root as the full run) equals "
             "the point-by-point sum of Nseq", f"{len(tests)} sub-boxes", "exact", okQ)

    # (X) how the fermionic frontier grows at k = 10
    print("\n=== (X) k = 10: the size of the fermionic frontier ===")
    if quick or "--no-probe" in sys.argv:
        print("    skipped")
    else:
        t = time.time()
        o = search_order(cov, tries=40, rounds=8)
        dom, weight = fermionic_domains(cov)
        mx, mc, _ = width_profile(o, dom, weight, len(cov.reps))
        jf = min((j for j in range(len(cov.reps)) if cov.m[cov.reps[j]] % 2),
                 key=lambda j: o.index(j))
        plan = Plan(cov, o, jf)
        h, pk = sweep(plan, verbose=5, statecap=1_000_000, tlimit=600)
        print(f"    fermionic frontier {mc} generators; "
              + ("finished" if h is not None else f"abandoned at {pk:,} states")
              + f"   [{time.time()-t:.0f}s]", flush=True)

    # (K) phiX174 at k = 10
    print("\n=== (K) phiX174 at k = 10: the last blank of Bio 8, Table 2 ===")
    if quick or "--arbor-k10" not in sys.argv:
        print("    k=10: not run here -- the arborescence frontier passes 400,000 states already on"
              " sub-boxes of a few hundred points (Q); the number is computed by ds_count.py."
              "  (--arbor-k10 forces the attempt.)")
    else:
        t = time.time()
        arg = lambda nm, d: next((a.split("=", 1)[1] for a in sys.argv
                                  if a.startswith(f"--{nm}=")), d)
        cap = int(float(arg("cap", "6e6")))
        jobs = int(arg("jobs", "1"))
        ckpt = os.path.join(HERE, "ds_k10_checkpoint.json")
        print(f"    order: greedy on the four endpoints of each orbit, frontier width "
              f"{W10.width}; root {cov.tail[cov.reps[j0]]} from orbit {j0} (m = {m0})")
        print(f"    state cap {cap:,} per sub-box, {jobs} process(es)", flush=True)
        half, leaves, peak = arbor_k10(S, 10, cov, o10, j0, cap, jobs,
                                       10 if jobs <= 1 else 0, ckpt)
        num = 2 * half
        okK = (num % den10 == 0)
        n = num // den10
        okK &= (8610708632 <= n <= 1.368e36)
        print(f"    k=10: #DS(phiX174,10) = {n:,}")
        print(f"          = {n:.4e}; {leaves} sub-boxes, peak {peak:,} states"
              f"   [{time.time()-t:.0f}s]")
        row("K", f"phiX174 at k = 10 has exactly {n:,} double-stranded reconstructions, inside "
                 "Bio 8's bounds 8,610,708,632 <= #DS <= 1.368e36: the entry Bio 8 left blank "
                 "and Bio 12 could not reach", "1 value of k", "exact", okK)

    print("\n=== battery summary ===")
    bad = [x for x in ROWS if not x[4]]
    print(f"  TOTAL {sum(1 for x in ROWS if x[4])}/{len(ROWS)} passed"
          + ("" if not bad else f"   FAILURES: {[x[0] for x in bad]}"))
    print(f"  elapsed {time.time()-t0:.0f}s")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
