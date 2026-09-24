# -*- coding: utf-8 -*-
r"""
ds_transfer_matrix.py -- the double-stranded lattice sum by frontier elimination.

Bio 8 wrote the number of double-stranded reconstructions of a circular sequence S as a sum
over the lattice points of a box,

    #DS(S,k) = sum over f in Lambda^- cap {|f| <= m, f = m mod 2, supp((m+f)/2) connected}
                   of  N_seq((m+f)/2),        N_seq(c) = t(G_c) Loc(G_c) / prod_a c(a)! ,

and enumerated the box by backtracking.  That works at rank(Lambda^-) = 4 and 12 and dies at
rank 48: phiX174 at k = 10 is the blank row of Bio 8's Table 2.

This note replaces the enumeration by a FRONTIER ELIMINATION.  The constraint "f is a
circulation" is one linear equation per vertex, each involving only the four endpoints of an
orbit {a, rho a}.  Processing the orbits in some order and carrying, as state, the partial
imbalance at the vertices that are still open, turns the count into a product of sparse
transfer operators.  Its cost is governed by the width of the frontier, not by the rank.

    python -u ds_transfer_matrix.py

Standalone: Python 3.8+ and the standard library.  No numpy, no sympy, no other file of the
series; the phiX174 sequence is read from phix174_NC_001422.1.txt beside this script.
"""
import os, sys, time, math, collections
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
    for nm in ("phix174_NC_001422.1.txt", os.path.join("..", "Ch07_Bio_08_DS_Lattice_Sum", "phix174_NC_001422.1.txt")):
        p = os.path.join(HERE, nm)
        if os.path.exists(p):
            txt = open(p, encoding="utf-8", errors="replace").read()
            s = "".join(c for c in txt.upper() if c in "ACGT")
            if len(s) > 5000:
                return s
    raise SystemExit("phix174_NC_001422.1.txt not found")


# --------------------------------------------------------------- the double cover, labelled

class Cover:
    """D_k(S): vertices = k-mers repeated in the two-strand pool, labelled arcs = the maximal
    stretches between consecutive such k-mers along each strand, rho = reverse complement.
    Rebuilt here from the definition rather than imported, and checked in (D) against the
    invariants Bio 8 published."""

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
                content = TT[p:p + (q - p + k)]
                self.arcs.append(dict(strand=si, pos=p, len=q - p + k, content=content))
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
        # the spectra of S and Sbar, as content vectors
        self.c_strand = [[0] * self.nA, [0] * self.nA]
        for a in self.arcs:
            self.c_strand[a["strand"]][self.cid[a["content"]]] += 1

    def constraints(self):
        """orbit j -> {vertex: coefficient of g_j in the net outflow}, where f(rep) = g_j and
        f(rho rep) = -g_j.  At most four vertices, and fewer when endpoints coincide."""
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

    def rank_lattice(self):
        var = self.constraints()
        if not var:
            return 0
        p = 2147483647
        cols = [{v: c % p for v, c in d.items()} for d in var]
        rows = collections.defaultdict(set)
        for j, c in enumerate(cols):
            for r in c:
                rows[r].add(j)
        alive, rank = set(range(len(cols))), 0
        while alive:
            j = min(alive, key=lambda x: len(cols[x]))
            c = cols[j]
            alive.discard(j)
            if not c:
                continue
            r = min(c, key=lambda x: len(rows[x]))
            inv = pow(c[r], p - 2, p)
            rank += 1
            for x in [y for y in rows[r] if y in alive]:
                f = (cols[x][r] * inv) % p
                for rr, vv in c.items():
                    nv = (cols[x].get(rr, 0) - f * vv) % p
                    if nv:
                        if rr not in cols[x]:
                            rows[rr].add(x)
                        cols[x][rr] = nv
                    elif rr in cols[x]:
                        del cols[x][rr]
                        rows[rr].discard(x)
            for rr in list(c):
                rows[rr].discard(j)
        return len(self.reps) - rank

    # --- the point attached to a lattice vector -----------------------------------------
    def point(self, g):
        c = [0] * self.nA
        for j, a in enumerate(self.reps):
            c[a] = (self.m[a] + g[j]) // 2
            c[self.crho[a]] = (self.m[a] - g[j]) // 2
        for a in range(self.nA):
            if self.crho[a] == a:
                c[a] = self.m[a] // 2
        return c

    def n_seq(self, c):
        """t(G_c) Loc(G_c) / prod c(a)!  -- zero when the support is disconnected, because
        the reduced Laplacian of a disconnected graph is singular (Lemma 3)"""
        outd = collections.Counter()
        for a in range(self.nA):
            if c[a]:
                outd[self.tail[a]] += c[a]
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


def det(M):
    """exact fraction-free-enough determinant; M is a list of rows of Fractions"""
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


# ------------------------------------------------------------------ frontier elimination

def coord_order(cov):
    """the orbits in the order of their position along the reference genome"""
    pos = {}
    for a in cov.arcs:
        c = cov.cid[a["content"]]
        p = a["pos"] if a["strand"] == 0 else (cov.N - a["pos"] - a["len"]) % cov.N
        pos[c] = min(pos.get(c, 1 << 60), p)
    return sorted(range(len(cov.reps)), key=lambda j: pos.get(cov.reps[j], 1 << 60))


def greedy_order(var, extra=None):
    """the order that greedily keeps the frontier narrow"""
    dom = [set(var[j]) | (set(extra[j]) if extra else set()) for j in range(len(var))]
    touch = collections.defaultdict(set)
    for j, d in enumerate(dom):
        for v in d:
            touch[v].add(j)
    left = {v: set(s) for v, s in touch.items()}
    todo, open_v, order = set(range(len(var))), set(), []
    while todo:
        best, bk = None, None
        for j in todo:
            nxt = open_v | dom[j]
            cl = {v for v in nxt if not (left[v] - {j})}
            key = (len(nxt) - len(cl), -len(cl))
            if bk is None or key < bk:
                best, bk = j, key
        todo.discard(best)
        order.append(best)
        for v in dom[best]:
            left[v].discard(best)
        open_v = {v for v in (open_v | dom[best]) if left[v]}
    return order


def frontier_widths(var, order, dom=None):
    """the frontier of an order, independently of running the elimination: a vertex is open
    from the first orbit touching it to the last."""
    dom = dom or [set(d) for d in var]
    touch = collections.defaultdict(set)
    for j, d in enumerate(dom):
        for v in d:
            touch[v].add(j)
    done, open_v, widths = set(), set(), []
    for j in order:
        done.add(j)
        open_v = {v for v in (open_v | dom[j]) if touch[v] - done}
        widths.append(len(open_v))
    return widths


def _schedule(var, m, order, dom=None):
    dom = dom or [set(d) for d in var]
    touch = collections.defaultdict(set)
    for j, d in enumerate(dom):
        for v in d:
            touch[v].add(j)
    capv = {v: sum(abs(var[j].get(v, 0)) * m[j] for j in s) for v, s in touch.items()}
    rem = [dict(capv)]
    for j in order:
        nc = dict(rem[-1])
        for v, c in var[j].items():
            nc[v] -= abs(c) * m[j]
        rem.append(nc)
    return touch, rem


def count_box(cov, order, cap=40_000_000, weights=None):
    """|Lambda^- cap box| by frontier elimination.  `weights[j]` optionally gives, for each
    value g of orbit j, a multiplier, so that any product-form weight can be summed."""
    var = cov.constraints()
    m = [cov.m[a] for a in cov.reps]
    touch, rem = _schedule(var, m, order)
    states, frontier, peak, done, widths = {(): 1}, [], 0, set(), []
    for i, j in enumerate(order):
        done.add(j)
        newly = [v for v in var[j] if v not in frontier]
        nxt = frontier + newly
        closing = [v for v in nxt if not (touch[v] - done)]
        keep = [v for v in nxt if v not in closing]
        vidx = {v: p for p, v in enumerate(nxt)}
        kidx, cidx = [vidx[v] for v in keep], [vidx[v] for v in closing]
        delta = [(vidx[v], c) for v, c in var[j].items()]
        lim = [rem[i + 1][v] for v in keep]
        new = collections.defaultdict(int)
        pad = len(newly)
        for st, n in states.items():
            base = list(st) + [0] * pad
            for g in range(-m[j], m[j] + 1, 2):
                b = base[:]
                for idx, c in delta:
                    b[idx] += c * g
                if any(b[t] for t in cidx):
                    continue
                if any(abs(b[t]) > lim[q] for q, t in enumerate(kidx)):
                    continue
                w = weights[j][g] if weights else 1
                if w:
                    new[tuple(b[t] for t in kidx)] += n * w
        states, frontier = dict(new), keep
        peak = max(peak, len(states))
        widths.append(len(frontier))
        if peak > cap:
            return None, peak, widths
    return states.get((), 0), peak, widths


def count_connected(cov, cap=5_000_000, verbose=False, tlimit=None):
    """the same elimination carrying, in addition, the partition of the frontier into the
    connected classes of the chosen support.  A class with no vertex left on the frontier is
    a finished component; the support is connected iff exactly one ever appears.  Classes
    that carry no chosen arc are not components of the support and are ignored."""
    var = cov.constraints()
    m = [cov.m[a] for a in cov.reps]
    tail, head, crho, reps = cov.tail, cov.head, cov.crho, cov.reps
    ends = [[tail[reps[j]], head[reps[j]], tail[crho[reps[j]]], head[crho[reps[j]]]]
            for j in range(len(reps))]
    dom = [set(var[j]) | set(ends[j]) for j in range(len(reps))]
    order = greedy_order(var, extra=ends)
    touch, rem = _schedule(var, m, order, dom)
    # the state is packed into one bytes object -- fin, then the balances (offset by 128), the
    # class labels and the support flags, each as long as the frontier.  At the sizes reached
    # here a tuple-of-tuples key costs several times more memory than the count itself.
    def pack(bal, part, flags, fin):
        return bytes([fin]) + bytes(x + 128 for x in bal) + bytes(part) + bytes(flags)

    def unpack(key):
        n = (len(key) - 1) // 3
        return ([x - 128 for x in key[1:1 + n]], key[1 + n:1 + 2 * n],
                key[1 + 2 * n:], key[0])

    states = {pack((), (), (), 0): 1}
    frontier, peak, done, t0 = [], 0, set(), time.time()
    for i, j in enumerate(order):
        done.add(j)
        a, b = reps[j], crho[reps[j]]
        newly = [v for v in dom[j] if v not in frontier]
        nxt = frontier + newly
        closing = [v for v in nxt if not (touch[v] - done)]
        keep = [v for v in nxt if v not in closing]
        vidx = {v: p for p, v in enumerate(nxt)}
        kidx, cidx = [vidx[v] for v in keep], [vidx[v] for v in closing]
        delta = [(vidx[v], c) for v, c in var[j].items()]
        lim = [rem[i + 1][v] for v in keep]
        nf, last = len(frontier), (i == len(order) - 1)
        new = collections.defaultdict(int)
        for key, n in states.items():
            st, pt, fl, fin = unpack(key)
            base = list(st) + [0] * len(newly)
            for g in range(-m[j], m[j] + 1, 2):
                bal = base[:]
                for idx, c in delta:
                    bal[idx] += c * g
                if any(bal[t] for t in cidx):
                    continue
                if any(abs(bal[t]) > lim[q] for q, t in enumerate(kidx)):
                    continue
                # union-find over nxt, seeded with the incoming partition.  A canonical class
                # label is NOT a vertex index: point each frontier vertex at the first vertex
                # carrying its label.
                par = list(range(len(nxt)))
                first = {}
                for q in range(nf):
                    if pt[q] not in first:
                        first[pt[q]] = q
                for q in range(nf):
                    par[q] = first[pt[q]]
                hasarc = [False] * len(nxt)
                for q in range(nf):
                    hasarc[q] = bool(fl[q])

                def find(x, P=par):
                    while P[x] != x:
                        P[x] = P[P[x]]
                        x = P[x]
                    return x
                if (m[j] + g) // 2:
                    x, y = find(vidx[tail[a]]), find(vidx[head[a]])
                    if x != y:
                        par[x] = y
                    hasarc[vidx[tail[a]]] = hasarc[vidx[head[a]]] = True
                if (m[j] - g) // 2:
                    x, y = find(vidx[tail[b]]), find(vidx[head[b]])
                    if x != y:
                        par[x] = y
                    hasarc[vidx[tail[b]]] = hasarc[vidx[head[b]]] = True
                # a class is part of the support iff one of its vertices meets a chosen arc
                arcroot = {find(q) for q in range(len(nxt)) if hasarc[q]}
                live = {find(t) for t in kidx}
                # a component is finished once no vertex of its class is still on the frontier;
                # count each CLASS once, not once per vertex
                closed = {find(t) for t in cidx} - live
                fin2 = fin + len(closed & arcroot)
                if fin2 > 1 or (last and fin2 != 1):
                    continue
                lab, part, flags = {}, [], []
                for t in kidx:
                    r = find(t)
                    if r not in lab:
                        lab[r] = len(lab)
                    part.append(lab[r])
                    flags.append(1 if r in arcroot else 0)   # per position, so that the three
                new[pack([bal[t] for t in kidx], part, flags, fin2)] += n   # arrays match
        states, frontier = dict(new), keep
        peak = max(peak, len(states))
        if verbose and (i % 10 == 0 or last):
            print(f"        orbit {i+1}/{len(order)}: frontier {len(frontier)}, "
                  f"states {len(states):,}   [{time.time()-t0:.0f}s]", flush=True)
        if len(states) > cap or (tlimit and time.time() - t0 > tlimit):
            return None, peak
    return sum(states.values()), peak


def enumerate_box(cov, cap=200000):
    """the backtracking enumerator of Bio 8, for the small k where it still runs"""
    var = cov.constraints()
    m = [cov.m[a] for a in cov.reps]
    order = greedy_order(var)
    touch, rem = _schedule(var, m, order)
    out, cur = [], [0] * len(var)
    bal = collections.Counter()

    def rec(i):
        if i == len(order):
            out.append(tuple(cur))
            return len(out) <= cap
        j = order[i]
        for g in range(-m[j], m[j] + 1, 2):
            for v, c in var[j].items():
                bal[v] += c * g
            ok = all(abs(bal[v]) <= rem[i + 1].get(v, 0) for v in var[j])
            if ok:
                cur[j] = g
                if not rec(i + 1):
                    for v, c in var[j].items():
                        bal[v] -= c * g
                    return False
            for v, c in var[j].items():
                bal[v] -= c * g
        return True
    rec(0)
    return out


# ------------------------------------------------------------------------------ the checks

def main():
    t0 = time.time()
    S = load_phix()
    print("=== the double-stranded lattice sum by frontier elimination ===")
    print(f"    phiX174 NC_001422.1, {len(S)} bp circular; python {sys.version.split()[0]}\n")

    covers = {}
    print("      k   |V|  fixV  |A_c|  fixA   sum m   orbits   rank Lambda^-")
    okD = True
    want = {12: (7, 1, 14, 0, 14, 4), 11: (30, 0, 53, 1, 60, 12), 10: (125, 3, 222, 0, 252, 48)}
    for k in (12, 11, 10):
        cov = Cover(S, k)
        covers[k] = cov
        r = cov.rank_lattice()
        got = (cov.nv, cov.n_fixV, cov.nA, cov.n_fixA, sum(cov.m), r)
        okD &= (got == want[k])
        print(f"     {k:2d} {cov.nv:5d} {cov.n_fixV:5d} {cov.nA:6d} {cov.n_fixA:5d} {sum(cov.m):7d}"
              f" {len(cov.reps):8d} {r:15d}   {'ok' if got == want[k] else 'MISMATCH ' + str(want[k])}")
    row("D", "the double cover rebuilt from the definition reproduces the invariants of Bio 8, "
             "Table 1: |V|, fixed vertices, contents, fixed contents, sum m and rank Lambda^-",
        "k = 10, 11, 12", "exact", okD)

    # (T) the elimination against the enumerator
    print("\n=== (T) frontier elimination against the backtracking enumerator ===")
    okT = True
    for k, want_n in ((12, 10), (11, 528)):
        cov = covers[k]
        n, peak, w = count_box(cov, greedy_order(cov.constraints()))
        e = len(enumerate_box(cov))
        good = (n == e == want_n)
        okT &= good
        print(f"    k={k}: elimination {n}, enumeration {e}, Bio 8 {want_n}, "
              f"peak states {peak}, max frontier {max(w)}   {'ok' if good else 'MISMATCH'}")
    row("T", "the frontier elimination returns the same box-point count as Bio 8's enumerator",
        "k = 11, 12", "exact", okT)

    # (C) connectivity is the vanishing of the matrix-tree determinant
    print("\n=== (C) connectivity of the support = non-vanishing of t(G_c) ===")
    okC, tot = True, 0
    for k in (12, 11):
        cov = covers[k]
        agree = dis = 0
        for g in enumerate_box(cov):
            c = cov.point(g)
            supp = sorted({cov.tail[a] for a in range(cov.nA) if c[a]}
                          | {cov.head[a] for a in range(cov.nA) if c[a]})
            par = {v: v for v in supp}

            def find(x):
                while par[x] != x:
                    par[x] = par[par[x]]
                    x = par[x]
                return x
            for a in range(cov.nA):
                if c[a]:
                    x, y = find(cov.tail[a]), find(cov.head[a])
                    if x != y:
                        par[x] = y
            conn = len({find(v) for v in supp}) == 1
            if conn == (cov.n_seq(c) > 0):
                agree += 1
            else:
                dis += 1
        tot += agree + dis
        okC &= (dis == 0)
        print(f"    k={k}: {agree} agree, {dis} disagree, of {agree+dis} box points")
    row("C", "supp(c) is connected exactly when the reduced Laplacian of the support is "
             "non-singular, so the connectivity filter of Bio 8 is an optimisation and not a "
             "side condition -- Lemma 3", f"{tot} points", "exact", okC)

    # (W) the weighted sum reproduces #DS
    print("\n=== (W) the lattice sum reproduces Bio 8's #DS ===")
    okW = True
    for k, want_ds, want_conn in ((12, 10, 8), (11, 864, 218)):
        cov = covers[k]
        pts = enumerate_box(cov)
        vals = [cov.n_seq(cov.point(g)) for g in pts]
        ds = sum(vals)
        cn = sum(1 for v in vals if v)
        good = (ds == want_ds and cn == want_conn)
        okW &= good
        print(f"    k={k}: connected points {cn} (Bio 8 {want_conn}), "
              f"#DS = {ds} (Bio 8 {want_ds})   {'ok' if good else 'MISMATCH'}")
    row("W", "summing N_seq over the box reproduces the published #DS and the published number "
             "of connected points", "k = 11, 12", "exact", okW)

    # (O) which order makes the frontier narrow
    print("\n=== (O) the genomic coordinate order is not the good order ===")
    lines = {}
    for k in (12, 11, 10):
        cov = covers[k]
        var = cov.constraints()
        for tag, o in (("coordinate", coord_order(cov)), ("greedy", greedy_order(var))):
            width = max(frontier_widths(var, o))      # a property of the order alone
            t = time.time()
            cap = 40_000_000 if tag == "greedy" else 1_500_000
            n, peak, _ = count_box(cov, o, cap=cap)
            lines[(k, tag)] = (n, peak, width)
            print(f"    k={k:2d} {tag:11s}: width {width:3d}, "
                  f"{('points %d' % n) if n is not None else 'ABORTED':>22},"
                  f" peak states {peak:>9,}   [{time.time()-t:.1f}s]", flush=True)
    narrower = all(lines[(k, "greedy")][2] <= lines[(k, "coordinate")][2] for k in (12, 11, 10))
    narrower &= lines[(10, "coordinate")][0] is None and lines[(10, "greedy")][0] is not None
    row("O", "ordering the orbits along the genome gives a frontier of width "
             f"{lines[(10,'coordinate')][2]} at k=10 and the elimination does not finish, while "
             f"the greedy order has width {lines[(10,'greedy')][2]} and finishes in seconds: the "
             "genome's own coordinate is not the good elimination order, because inverted "
             "repeats interleave", "k = 10, 11, 12", "exact", narrower)

    # (K) the blank row of Bio 8's Table 2
    print("\n=== (K) phiX174 at k = 10: the blank row of Bio 8, Table 2 ===")
    cov = covers[10]
    t = time.time()
    n10, peak10, w10 = count_box(cov, greedy_order(cov.constraints()))
    print(f"    rank Lambda^- = 48, orbits = {len(cov.reps)}")
    print(f"    box points = {n10:,}   peak states {peak10:,}, max frontier {max(w10)}"
          f"   [{time.time()-t:.1f}s]")
    # an independent check: pin one coordinate and sum over its values
    o = greedy_order(cov.constraints())
    parts = []
    for j in (0, 7):
        mj = cov.m[cov.reps[j]]
        s = 0
        for g in range(-mj, mj + 1, 2):
            wts = [None] * len(cov.reps)
            wts[j] = collections.defaultdict(int)
            wts[j][g] = 1
            wt = [{gg: (1 if (jj != j or gg == g) else 0)
                   for gg in range(-cov.m[cov.reps[jj]], cov.m[cov.reps[jj]] + 1, 2)}
                  for jj in range(len(cov.reps))]
            nn, _, _ = count_box(cov, o, weights=wt)
            s += nn
        parts.append((j, s))
        print(f"    pinning orbit {j} and summing over its values: {s:,}"
              f"   {'ok' if s == n10 else 'MISMATCH'}")
    okK = (n10 == 2481687900) and all(s == n10 for _, s in parts)
    row("K", f"phiX174 at k = 10 has exactly {n10:,} box points, a row Bio 8 left blank; the "
             "count is confirmed by pinning a coordinate and summing over its values",
        "1 value of k", "exact", okK)

    # (N) the connected points
    print("\n=== (N) connected points by frontier elimination with connectivity classes ===")
    okN = True
    for k, want_conn in ((12, 8), (11, 218)):
        n, peak = count_connected(covers[k])
        good = (n == want_conn)
        okN &= good
        print(f"    k={k}: connected = {n} (Bio 8 {want_conn}), peak states {peak:,}"
              f"   {'ok' if good else 'MISMATCH'}")
    row("N", "carrying the partition of the frontier into connected classes counts the "
             "connected points directly, reproducing Bio 8's 8 and 218",
        "k = 11, 12", "exact", okN)
    t = time.time()
    if "--quick" in sys.argv:
        print("    k=10: skipped (--quick); it takes about eighty minutes and 10 GB")
    else:
        n, peak = count_connected(covers[10], cap=40_000_000, verbose=True)
        print(f"    k=10: connected = {n:,} of {n10:,} box points "
              f"({100.0*n/n10:.1f}%)   peak {peak:,} states   [{time.time()-t:.0f}s]")
        row("M", f"phiX174 at k = 10 has exactly {n:,} box points with connected support, the "
                 "second blank column of Bio 8's Table 2; the connected fraction falls from "
                 "80% at k=12 through 41.3% at k=11 to 36.6% here",
            "1 value of k", "exact", n == 909103210)

    print("\n=== battery summary ===")
    bad = [x for x in ROWS if not x[4]]
    print(f"  TOTAL {sum(1 for x in ROWS if x[4])}/{len(ROWS)} passed"
          + ("" if not bad else f"   FAILURES: {[x[0] for x in bad]}"))
    print(f"  elapsed {time.time()-t0:.0f}s")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
