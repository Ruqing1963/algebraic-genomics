# -*- coding: utf-8 -*-
r"""
ds_count.py -- #DS(phiX174, 10) by compiled enumeration of the half-box.

ds_fermion.py proves the pieces: the Berezin identity, the fixed root on the half-box, the
halving by c -> m-c, the integrality, and an exact positive frontier (arborescence classes)
that reproduces #DS at k = 11, 12 and on enumerable k = 10 sub-boxes.  At k = 10 both
frontiers are too wide for Python.  This script instead VISITS every point of the half-box
f_{j0} > 0 -- 1 240 843 950 of them, half of Bio 12's 2 481 687 900 -- with the same
imbalance-bounded backtracking as Bio 8, compiled with numba, and evaluates Nseq(c) at each:

  * connectivity of supp c by union-find (Bio 12, Lemma 3: t = 0 otherwise);
  * t(G_c) by contracting forced arcs -- a vertex whose out-arcs all go to one vertex u must
    use one of them in every in-tree, so t(G) = w t(G / (v->u)) -- and a determinant of the
    small remainder;
  * everything modulo five primes below 2^31, Loc and 1/prod c(a)! from tables; the half-sum
    is recovered by the Chinese remainder theorem, which is exact because it is below Bio 8's
    upper bound 1.368e36 < the product of the primes.

Checks, all exact:
  (P) half-box point count and connected count at k = 10 against Bio 12 (halved: the
      involution c -> m-c has no fixed point on the box and preserves connectivity);
  (W) #DS = 10 and 864 at k = 12, 11, and brute force on random covers;
  (X) on k = 10 sub-boxes, against the exact Python arborescence frontier of ds_fermion.py.

    python -u ds_count.py            everything
    python -u ds_count.py --quick    everything but the full k = 10 run
    --threads=N  numba threads (default: all)

Needs numpy and numba (Anaconda has both), and ds_fermion.py beside it.
"""
import os, sys, time, math, json, random, collections

import numpy as np
import numba
from numba import njit, prange

import ds_fermion as F

HERE = os.path.dirname(os.path.abspath(__file__))
ROWS = []


def row(key, desc, cases, ok):
    ROWS.append((key, bool(ok)))
    print(("  [ok]   " if ok else "  [FAIL] ") + f"({key}) {desc}   [{cases}] [exact]", flush=True)


# ------------------------------------------------------------------ modular arithmetic

def _is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


PRIMES = []
_q = (1 << 31) - 1
while len(PRIMES) < 5:
    if _is_prime(_q):
        PRIMES.append(_q)
    _q -= 2
MAXF = 600


def tables():
    P = np.array(PRIMES, dtype=np.int64)
    fac = np.zeros((len(PRIMES), MAXF), dtype=np.int64)
    ifac = np.zeros((len(PRIMES), MAXF), dtype=np.int64)
    for k, p in enumerate(PRIMES):
        f = 1
        for n in range(MAXF):
            if n:
                f = f * n % p
            fac[k, n] = f
            ifac[k, n] = pow(f, p - 2, p)
    return P, fac, ifac


def crt(res):
    x, M = 0, 1
    for r, p in zip(res, PRIMES):
        t = ((int(r) - x) * pow(M, -1, p)) % p
        x += M * t
        M *= p
    return x, M


# ------------------------------------------------------------------ the compiled kernel

@njit(cache=True)
def _find(uf, x):
    while uf[x] != x:
        uf[x] = uf[uf[x]]
        x = uf[x]
    return x


@njit(cache=True)
def _powmod(a, e, p):
    r = 1
    a %= p
    while e > 0:
        if e & 1:
            r = r * a % p
        a = a * a % p
        e >>= 1
    return r


@njit(cache=True)
def _leaf(c, nA, nv, tail, head, root, P, fac, ifac,
          deg, uf, outN, outW, no, inN, ni, alive, stack, M, idx, res):
    """1 and res[k] = Nseq(c) mod P[k] if supp c is connected; 0 if not; -1 on overflow"""
    for v in range(nv):
        deg[v] = 0
        uf[v] = v
    for a in range(nA):
        if c[a] > 0:
            deg[tail[a]] += c[a]
            x = _find(uf, tail[a])
            y = _find(uf, head[a])
            if x != y:
                uf[x] = y
    rr = _find(uf, root)
    for v in range(nv):
        if deg[v] > 0 and _find(uf, v) != rr:
            return 0
    K = outN.shape[1]
    for v in range(nv):
        no[v] = 0
        ni[v] = 0
        alive[v] = 1 if deg[v] > 0 else 0
    for a in range(nA):
        w = c[a]
        if w == 0:
            continue
        t = tail[a]
        h = head[a]
        if t == h:
            continue
        f = -1
        for q in range(no[t]):
            if outN[t, q] == h:
                f = q
                break
        if f >= 0:
            outW[t, f] += w
        else:
            if no[t] >= K or ni[h] >= K:
                return -1
            outN[t, no[t]] = h
            outW[t, no[t]] = w
            no[t] += 1
            inN[h, ni[h]] = t
            ni[h] += 1
    nk = P.shape[0]
    for k in range(nk):
        res[k] = 1
    sp = 0
    for v in range(nv):
        if alive[v] == 1 and v != root and no[v] == 1:
            stack[sp] = v
            sp += 1
    while sp > 0:
        sp -= 1
        v = stack[sp]
        if alive[v] == 0 or no[v] != 1:
            continue
        u = outN[v, 0]
        w = outW[v, 0]
        for k in range(nk):
            res[k] = res[k] * (w % P[k]) % P[k]
        for q in range(ni[u]):
            if inN[u, q] == v:
                ni[u] -= 1
                inN[u, q] = inN[u, ni[u]]
                break
        no[v] = 0
        for qq in range(ni[v]):
            x = inN[v, qq]
            wx = 0
            for q in range(no[x]):
                if outN[x, q] == v:
                    wx = outW[x, q]
                    no[x] -= 1
                    outN[x, q] = outN[x, no[x]]
                    outW[x, q] = outW[x, no[x]]
                    break
            if x != u:
                f = -1
                for q in range(no[x]):
                    if outN[x, q] == u:
                        f = q
                        break
                if f >= 0:
                    outW[x, f] += wx
                else:
                    if no[x] >= K or ni[u] >= K:
                        return -1
                    outN[x, no[x]] = u
                    outW[x, no[x]] = wx
                    no[x] += 1
                    inN[u, ni[u]] = x
                    ni[u] += 1
            if x != root and alive[x] == 1 and no[x] == 1:
                if sp >= stack.shape[0]:
                    return -1
                stack[sp] = x
                sp += 1
        ni[v] = 0
        alive[v] = 0
        if u != root and no[u] == 1:
            if sp >= stack.shape[0]:
                return -1
            stack[sp] = u
            sp += 1
    # what is left: a small reduced Laplacian
    s = 0
    for v in range(nv):
        if alive[v] == 1 and v != root:
            if no[v] == 0:
                return 0
            idx[v] = s
            s += 1
    if s > M.shape[0]:
        return -1
    for k in range(nk):
        p = P[k]
        for i in range(s):
            for j in range(s):
                M[i, j] = 0
        for v in range(nv):
            if alive[v] == 1 and v != root:
                i = idx[v]
                for q in range(no[v]):
                    u = outN[v, q]
                    w = outW[v, q] % p
                    M[i, i] = (M[i, i] + w) % p
                    if u != root:
                        M[i, idx[u]] = (M[i, idx[u]] - w) % p
        d = 1
        for col in range(s):
            piv = -1
            for r in range(col, s):
                if M[r, col] != 0:
                    piv = r
                    break
            if piv < 0:
                d = 0
                break
            if piv != col:
                for j in range(s):
                    tmp = M[col, j]
                    M[col, j] = M[piv, j]
                    M[piv, j] = tmp
                d = (p - d) % p
            d = d * M[col, col] % p
            inv = _powmod(M[col, col], p - 2, p)
            for r in range(col + 1, s):
                if M[r, col] != 0:
                    fct = M[r, col] * inv % p
                    for j in range(col, s):
                        M[r, j] = (M[r, j] - fct * M[col, j]) % p
        res[k] = res[k] * d % p
    for k in range(nk):
        for v in range(nv):
            if deg[v] > 0:
                res[k] = res[k] * fac[k, deg[v] - 1] % P[k]
        for a in range(nA):
            if c[a] > 0:
                res[k] = res[k] * ifac[k, c[a]] % P[k]
    return 1


@njit(cache=True)
def _run_one(pref, D, J, oa, ob, om, olo, ohi, tail, head, nA, nv, rem, bal0, c0, root,
             P, fac, ifac, out):
    """all points below one prefix: out = [points, connected, overflow, sum mod P[0..]]"""
    nk = P.shape[0]
    bal = bal0.copy()
    c = c0.copy()
    cur = np.zeros(J + 1, dtype=np.int64)
    applied = np.zeros(J + 1, dtype=np.int64)
    K = 64
    deg = np.zeros(nv, dtype=np.int64)
    uf = np.zeros(nv, dtype=np.int64)
    outN = np.zeros((nv, K), dtype=np.int64)
    outW = np.zeros((nv, K), dtype=np.int64)
    no = np.zeros(nv, dtype=np.int64)
    inN = np.zeros((nv, K), dtype=np.int64)
    ni = np.zeros(nv, dtype=np.int64)
    alive = np.zeros(nv, dtype=np.int64)
    stack = np.zeros(nv * K, dtype=np.int64)
    M = np.zeros((96, 96), dtype=np.int64)
    idx = np.zeros(nv, dtype=np.int64)
    res = np.zeros(nk, dtype=np.int64)
    for i in range(out.shape[0]):
        out[i] = 0
    for i in range(D):
        x = pref[i]
        a = oa[i]
        b = ob[i]
        bal[tail[a]] += x
        bal[head[a]] -= x
        bal[tail[b]] += om[i] - x
        bal[head[b]] -= om[i] - x
        c[a] = x
        c[b] = om[i] - x
    for i in range(D):
        for e in (tail[oa[i]], head[oa[i]], tail[ob[i]], head[ob[i]]):
            if abs(bal[e]) > rem[D, e]:
                return
    if D == J:
        st = _leaf(c, nA, nv, tail, head, root, P, fac, ifac,
                   deg, uf, outN, outW, no, inN, ni, alive, stack, M, idx, res)
        out[0] += 1
        if st == 1:
            out[1] += 1
            for k in range(nk):
                out[3 + k] = (out[3 + k] + res[k]) % P[k]
        elif st < 0:
            out[2] += 1
        return
    lev = D
    cur[lev] = olo[lev] - 1
    applied[lev] = 0
    while lev >= D:
        a = oa[lev]
        b = ob[lev]
        m = om[lev]
        if applied[lev] == 1:
            x = cur[lev]
            bal[tail[a]] -= x
            bal[head[a]] += x
            bal[tail[b]] -= m - x
            bal[head[b]] += m - x
            applied[lev] = 0
        cur[lev] += 1
        if cur[lev] > ohi[lev]:
            lev -= 1
            continue
        x = cur[lev]
        bal[tail[a]] += x
        bal[head[a]] -= x
        bal[tail[b]] += m - x
        bal[head[b]] -= m - x
        applied[lev] = 1
        c[a] = x
        c[b] = m - x
        r1 = rem[lev + 1]
        if (abs(bal[tail[a]]) > r1[tail[a]] or abs(bal[head[a]]) > r1[head[a]]
                or abs(bal[tail[b]]) > r1[tail[b]] or abs(bal[head[b]]) > r1[head[b]]):
            continue
        if lev == J - 1:
            st = _leaf(c, nA, nv, tail, head, root, P, fac, ifac,
                       deg, uf, outN, outW, no, inN, ni, alive, stack, M, idx, res)
            out[0] += 1
            if st == 1:
                out[1] += 1
                for k in range(nk):
                    out[3 + k] = (out[3 + k] + res[k]) % P[k]
            elif st < 0:
                out[2] += 1
            continue
        lev += 1
        cur[lev] = olo[lev] - 1
        applied[lev] = 0


@njit(parallel=True, cache=True)
def _run_batch(prefs, D, J, oa, ob, om, olo, ohi, tail, head, nA, nv, rem, bal0, c0, root,
               P, fac, ifac, out):
    for i in prange(prefs.shape[0]):
        _run_one(prefs[i], D, J, oa, ob, om, olo, ohi, tail, head, nA, nv, rem, bal0, c0,
                 root, P, fac, ifac, out[i])


# ------------------------------------------------------------------ the problem, as arrays

class Box:
    """the half-box f_{j0} > 0 of a cover (optionally a sub-box: pins j -> c(rep_j)) laid out
    for the kernel, in Bio 12's greedy order"""

    def __init__(self, cov, j0, pins=None):
        self.cov, self.j0 = cov, j0
        var = cov.constraints()
        order = F.greedy_order(var)
        J = len(order)
        reps, crho, tail, head, m = cov.reps, cov.crho, cov.tail, cov.head, cov.m
        self.J = J
        self.oa = np.array([reps[j] for j in order], dtype=np.int64)
        self.ob = np.array([crho[reps[j]] for j in order], dtype=np.int64)
        self.om = np.array([m[reps[j]] for j in order], dtype=np.int64)
        lo, hi = [], []
        for j in order:
            mj = m[reps[j]]
            l, h = 0, mj
            if j == j0:
                l = mj // 2 + 1
            if pins and j in pins:
                l, h = max(l, pins[j]), min(h, pins[j])
            lo.append(l)
            hi.append(h)
        self.olo = np.array(lo + [0], dtype=np.int64)
        self.ohi = np.array(hi + [-1], dtype=np.int64)
        self.tail = np.array(tail, dtype=np.int64)
        self.head = np.array(head, dtype=np.int64)
        nv = cov.nv
        rem = np.zeros((J + 1, nv), dtype=np.int64)
        for i in range(J - 1, -1, -1):
            rem[i] = rem[i + 1]
            a, b = reps[order[i]], crho[reps[order[i]]]
            for v in {tail[a], head[a], tail[b], head[b]}:
                ca = (tail[a] == v) - (head[a] == v)
                cb = (tail[b] == v) - (head[b] == v)
                rem[i, v] += m[a] * max(abs(ca), abs(cb))
        self.rem = rem
        bal0 = np.zeros(nv, dtype=np.int64)
        c0 = np.zeros(cov.nA, dtype=np.int64)
        for a in range(cov.nA):
            if crho[a] == a:
                c0[a] = m[a] // 2
                bal0[tail[a]] += m[a] // 2
                bal0[head[a]] -= m[a] // 2
        self.bal0, self.c0 = bal0, c0
        self.root = tail[reps[j0]]
        self.nA, self.nv = cov.nA, nv

    def prefixes(self, target):
        """all feasible assignments of the first D levels, D the least with >= target of them"""
        J = self.J
        for D in range(0, J + 1):
            out = []
            bal = [int(x) for x in self.bal0]
            cur = []
            tail, head = self.tail, self.head

            def rec(i):
                if i == D:
                    out.append(list(cur))
                    return
                a, b, m = int(self.oa[i]), int(self.ob[i]), int(self.om[i])
                ends = (int(tail[a]), int(head[a]), int(tail[b]), int(head[b]))
                for x in range(int(self.olo[i]), int(self.ohi[i]) + 1):
                    bal[ends[0]] += x
                    bal[ends[1]] -= x
                    bal[ends[2]] += m - x
                    bal[ends[3]] -= m - x
                    if all(abs(bal[e]) <= self.rem[i + 1, e] for e in ends):
                        cur.append(x)
                        rec(i + 1)
                        cur.pop()
                    bal[ends[0]] -= x
                    bal[ends[1]] += x
                    bal[ends[2]] -= m - x
                    bal[ends[3]] += m - x
            rec(0)
            if len(out) >= target or D == J:
                arr = np.zeros((len(out), max(D, 1)), dtype=np.int64)
                for p, pr in enumerate(out):
                    arr[p, :D] = pr
                return D, arr


def run_box(box, target=4000, batches=1, verbose=False, ckpt=None, sig=""):
    """(points, connected, overflow, half-sum of Nseq) over the box"""
    P, fac, ifac = TABLES
    D, prefs = box.prefixes(target)
    rng = np.random.default_rng(20260922)
    perm = rng.permutation(len(prefs))
    prefs = prefs[perm]
    tot = np.zeros(3 + len(PRIMES), dtype=object)
    tot[:] = 0
    book = {"sig": sig, "done": {}}
    if ckpt and os.path.exists(ckpt):
        try:
            old = json.load(open(ckpt, encoding="utf-8"))
            if old.get("sig") == sig:
                book = old
                print(f"    resuming: {len(book['done'])} of {batches} batches already done",
                      flush=True)
        except Exception:
            pass
    chunks = np.array_split(np.arange(len(prefs)), batches)
    t0 = time.time()
    tdone, ncomp = 0, 0
    for bi, ch in enumerate(chunks):
        if str(bi) in book["done"]:
            vals = [int(x) for x in book["done"][str(bi)]]
        else:
            t1 = time.time()
            out = np.zeros((len(ch), 3 + len(PRIMES)), dtype=np.int64)
            _run_batch(prefs[ch], D, box.J, box.oa, box.ob, box.om, box.olo, box.ohi,
                       box.tail, box.head, box.nA, box.nv, box.rem, box.bal0, box.c0,
                       box.root, P, fac, ifac, out)
            vals = [int(out[:, 0].sum()), int(out[:, 1].sum()), int(out[:, 2].sum())]
            for k, p in enumerate(PRIMES):
                vals.append(int(sum(int(x) for x in out[:, 3 + k]) % p))
            tdone += time.time() - t1
            ncomp += 1
            book["done"][str(bi)] = [str(x) for x in vals]
            if ckpt:
                tmp = ckpt + ".tmp"
                json.dump(book, open(tmp, "w", encoding="utf-8"))
                os.replace(tmp, ckpt)
        for i in range(3):
            tot[i] += vals[i]
        for k, p in enumerate(PRIMES):
            tot[3 + k] = (tot[3 + k] + vals[3 + k]) % p
        if verbose:
            el = time.time() - t0
            left = batches - bi - 1
            eta = tdone / ncomp * left if ncomp else 0
            print(f"      batch {bi+1}/{batches}: {int(tot[0]):,} points, "
                  f"{int(tot[1]):,} connected   [{el:.0f}s, eta ~{eta/60:.0f} min]", flush=True)
    half, _ = crt([tot[3 + k] for k in range(len(PRIMES))])
    return int(tot[0]), int(tot[1]), int(tot[2]), half, D, len(prefs)


TABLES = tables()


def main():
    t0 = time.time()
    quick = "--quick" in sys.argv
    thr = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--threads=")), None)
    if thr:
        numba.set_num_threads(int(thr))
    S = F.load_phix()
    print("=== #DS(phiX174, 10) by compiled enumeration of the half-box ===")
    print(f"    python {sys.version.split()[0]}, numba {numba.__version__}, "
          f"{numba.get_num_threads()} threads; primes {PRIMES}\n", flush=True)
    covers = {k: F.Cover(S, k) for k in (12, 11, 10)}

    # (W) small k and random covers
    print("=== (W) against Bio 8 and against brute force ===")
    okW = True
    for k, want, wpts, wconn in ((12, 10, 10, 8), (11, 864, 528, 218)):
        cov = covers[k]
        j0 = F.pick_j0(cov, F.greedy_order(cov.constraints()))
        npts, ncon, nof, half, D, npf = run_box(Box(cov, j0), target=8)
        good = (2 * half == want and 2 * npts == wpts and 2 * ncon == wconn and nof == 0)
        okW &= good
        print(f"    k={k}: half-box {npts} points ({ncon} connected), #DS = {2*half} "
              f"(Bio 8 {want})   {'ok' if good else 'MISMATCH'}")
    rnd = random.Random(20260922)
    nr = agree = 0
    while nr < 12:
        T = "".join(rnd.choice("ACGT") for _ in range(rnd.randrange(30, 70)))
        kk = rnd.randrange(3, 6)
        c2 = F.Cover(T, kk)
        if not (3 <= c2.nv <= 9 and c2.reps
                and any(c2.m[c2.reps[j]] % 2 for j in range(len(c2.reps)))):
            continue
        p2 = F.enumerate_box(c2)
        if not (1 <= len(p2) <= 4000):
            continue
        nr += 1
        brute = sum(c2.n_seq(c2.point(g)) for g in p2)
        j0 = F.pick_j0(c2, F.greedy_order(c2.constraints()))
        _, _, nof, half, _, _ = run_box(Box(c2, j0), target=8)
        agree += (2 * half == brute and nof == 0)
    okW &= (agree == nr)
    print(f"    {agree} of {nr} random covers agree with the brute-force sum of Nseq")
    row("W", "the compiled enumeration with forced-arc contraction and five-prime arithmetic "
             "reproduces Bio 8's #DS, point and connected counts at k = 11, 12, and the brute "
             "force on random covers", f"k = 11, 12 and {nr} random covers", okW)

    # (X) k = 10 sub-boxes against the exact Python arborescence frontier
    print("\n=== (X) k = 10 sub-boxes: compiled enumeration against the arborescence frontier ===")
    cov = covers[10]
    o_arb = F.arbor_order(cov)
    j0 = F.pick_j0(cov, o_arb)
    W = F.Windows(cov, o_arb)
    den = F.denominator(cov)
    rnd = random.Random(1310)
    J10 = len(cov.reps)
    okX, nx, t = True, 0, time.time()
    for npin in (40, 34, 30, 30, 26):
        g0 = None
        for _ in range(40):
            shuf = {}
            for j in range(J10):
                mj = cov.m[cov.reps[j]]
                vs = list(range(1, mj + 1, 2)) if j == j0 else list(range(-mj, mj + 1, 2))
                rnd.shuffle(vs)
                shuf[j] = vs
            p0 = F.enumerate_box(cov, cap=0, allowed=shuf, nodecap=300000)
            if p0 and F.n_seq_fast(cov, cov.point(p0[0])):
                g0 = p0[0]
                break
        if g0 is None:
            continue
        free = [j for j in range(J10) if j != j0]
        rnd.shuffle(free)
        pins = {j: (cov.m[cov.reps[j]] + g0[j]) // 2 for j in free[:npin]}
        h, pk, ab = F.arbor_sweep(cov, o_arb, j0, pins=pins, W=W, statecap=400000)
        if h is None:
            print(f"    {npin} orbits pinned: arborescence frontier above 400,000 states, skipped")
            continue
        npts, ncon, nof, half, _, _ = run_box(Box(cov, j0, pins), target=64)
        good = (h == half * den and nof == 0)
        okX &= good
        nx += 1
        print(f"    {npin} orbits pinned: {npts:,} points, {ncon:,} connected, sum Nseq = {half:,}"
              f"; frontier (peak {pk:,}) {'agrees' if good else 'DISAGREES'}", flush=True)
    okX &= nx >= 3
    print(f"    [{time.time()-t:.0f}s]")
    row("X", "on k = 10 sub-boxes the compiled enumeration and the exact arborescence frontier "
             "of ds_fermion.py -- two unrelated algorithms -- give the same sum of Nseq",
        f"{nx} sub-boxes", okX)

    # (K) the full half-box
    print("\n=== (K) phiX174 at k = 10: the whole half-box ===")
    if quick:
        print("    skipped (--quick)")
    else:
        t = time.time()
        box = Box(cov, j0)
        ck = os.path.join(HERE, "ds_count_k10_checkpoint.json")
        npts, ncon, nof, half, D, npf = run_box(
            box, target=20000, batches=100, verbose=True, ckpt=ck,
            sig=f"k=10;j0={j0};primes={PRIMES}")
        n = 2 * half
        M4 = PRIMES[0] * PRIMES[1] * PRIMES[2] * PRIMES[3] * PRIMES[4]
        print(f"    {npf:,} prefixes of depth {D}; root {box.root}, orbit {j0}")
        print(f"    half-box points    {npts:,}   (Bio 12: 2,481,687,900 / 2 = 1,240,843,950)")
        print(f"    connected          {ncon:,}   (Bio 12:   909,103,210 / 2 =   454,551,605)")
        print(f"    overflows          {nof}")
        print(f"    #DS(phiX174, 10) = {n:,}")
        print(f"                     = {float(n):.6e}   [{time.time()-t:.0f}s]")
        okP = (npts == 1240843950 and ncon == 454551605 and nof == 0)
        row("P", "the half-box has 1,240,843,950 points of which 454,551,605 have connected "
                 "support, exactly half of Bio 12's two counts", "k = 10", okP)
        okK = okP and (8610708632 <= n <= 1.368e36) and half < M4 // 2
        row("K", f"#DS(phiX174, 10) = {n:,}, inside Bio 8's bounds 8,610,708,632 <= #DS <= "
                 "1.368e36", "k = 10", okK)

    print("\n=== battery summary ===")
    bad = [k for k, ok in ROWS if not ok]
    print(f"  TOTAL {len(ROWS)-len(bad)}/{len(ROWS)} passed" + (f"   FAILURES: {bad}" if bad else ""))
    print(f"  elapsed {time.time()-t0:.0f}s")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
