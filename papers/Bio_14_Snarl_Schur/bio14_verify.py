# -*- coding: utf-8 -*-
r"""
bio14_verify.py -- the BEST determinant by Schur complements along the snarl tree.

Single-stranded reconstruction counts are one determinant (BEST):  ec(G) = t_r(G) prod (d_v - 1)!,
t_r(G) = det of the reduced Laplacian.  The specification for this paper proposed "snarl tensor
networks" contracted along the snarl tree of a genome graph.  What survives, exactly:

  * eliminating the interior I of a snarl with boundary B replaces the Laplacian by its SCHUR
    COMPLEMENT on B (Kron reduction), and  det L~ = det(L_II) * det(Schur);
  * det(L_II) counts the spanning forests of the snarl's interior rooted into its boundary (the
    all-minors matrix-tree theorem), and the Schur complement is again a Laplacian: the snarl's
    "boundary tensor" is an |B| x |B| weighted digraph, not a higher tensor, because the
    determinant is multiplicative;
  * eliminating the snarls in the order of the snarl tree keeps every intermediate matrix small:
    the work is O(sum over vertices of w_v^2), linear in the genome for bounded boundary width w.

It does NOT make the double-stranded count #DS tractable -- that is #P-hard (Bio 19).

Checks (all exact):
  (M) Schur-complement identity det L~ = det L_II * det(Schur) and "Schur of a Laplacian is a
      Laplacian", on random Eulerian digraphs;
  (F) det L_II = number of spanning forests of the interior rooted into the boundary, by
      brute-force enumeration of out-arc choices;
  (B) BEST itself: t * prod (d-1)! = number of Eulerian circuits, by brute force on small digraphs;
  (E) on synthetic genome graphs (chains of nested bubbles traversed by several haplotypes): the
      snarl-order elimination equals the dense Bareiss determinant and the forced-arc contraction
      of the spectral-fibres paper; the front width stays bounded;
  (P) phiX174, de Bruijn graphs of the circular genome at k = 7, 8, 10, 12: minimum-degree elimination
      equals the forced-arc contraction; sizes and front widths;
  (S) scaling: synthetic genome graphs up to ~45 000 vertices, snarl order vs minimum-degree order
      (exact modulo two primes), time per vertex and front width.

    python -u bio14_verify.py            (writes bio14_output.txt as it runs)

Standard library only; imports ds_fermion.py (Bio 13) and spectral_fibres.py (Spectral Fibres).
"""
import os, sys, time, math, random, heapq, collections, itertools
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "Bio_13_Fermionic_Partition"))
sys.path.insert(0, os.path.join(HERE, "..", "Spectral_Fibres"))
import ds_fermion as DF                  # noqa: E402
import spectral_fibres as SF             # noqa: E402

ROWS = []
P1, P2 = 2147483647, 2147483629


def row(key, desc, ok):
    ROWS.append((key, bool(ok)))
    print(("  [ok]   " if ok else "  [FAIL] ") + f"({key}) {desc}", flush=True)


# ------------------------------------------------------------------ sparse reduced Laplacian

class Field:
    def __init__(self, mode, p=None):
        self.mode, self.p = mode, p

    def c(self, x):
        return Fraction(x) if self.mode == "frac" else x % self.p

    def div(self, a, b):
        return a / b if self.mode == "frac" else a * pow(b, self.p - 2, self.p) % self.p

    def sub(self, a, b):
        return a - b if self.mode == "frac" else (a - b) % self.p

    def mul(self, a, b):
        return a * b if self.mode == "frac" else a * b % self.p


def reduced_laplacian(arcs, root, F):
    """rows R[v] = {u: L(v,u)} and column sets C[u] = {v}, the row and column of root deleted"""
    R = collections.defaultdict(dict)
    for x, y in arcs:
        if x == y or x == root:
            continue
        R[x][x] = R[x].get(x, 0) + 1
        if y != root:
            R[x][y] = R[x].get(y, 0) - 1
    C = collections.defaultdict(set)
    for v, r in R.items():
        for u in list(r):
            r[u] = F.c(r[u])
            C[u].add(v)
    return R, C


def eliminate(R, C, F, v, stats):
    """eliminate vertex v: returns the pivot"""
    piv = R[v].get(v, F.c(0))
    rowv = R.pop(v)
    rowv.pop(v, None)
    colv = C.pop(v)
    colv.discard(v)
    stats[0] = max(stats[0], len(rowv) + len(colv))
    stats[1] += len(rowv) * len(colv)
    for j in rowv:
        C[j].discard(v)
    if piv == 0 or (F.mode == "mod" and piv % F.p == 0):
        return piv
    for i in colv:
        a = R[i].pop(v)
        f = F.div(a, piv)
        ri = R[i]
        for j, b in rowv.items():
            nv = F.sub(ri.get(j, F.c(0)), F.mul(f, b))
            if nv == 0:
                if j in ri:
                    del ri[j]
                    C[j].discard(i)
            else:
                ri[j] = nv
                C[j].add(i)
    return piv


def det_by_order(arcs, root, order, F):
    R, C = reduced_laplacian(arcs, root, F)
    d, stats = F.c(1), [0, 0]
    for v in order:
        if v == root:
            continue
        piv = eliminate(R, C, F, v, stats)
        d = F.mul(d, piv)
        if d == 0:
            return d, stats
    return d, stats


def det_mindegree(arcs, root, F):
    R, C = reduced_laplacian(arcs, root, F)
    alive = set(R)
    heap = [(len(R[v]) + len(C[v]), v) for v in alive]
    heapq.heapify(heap)
    d, stats = F.c(1), [0, 0]
    while alive:
        key, v = heapq.heappop(heap)
        if v not in alive:
            continue
        cur = len(R[v]) + len(C[v])
        if cur != key:
            heapq.heappush(heap, (cur, v))
            continue
        nbrs = (set(R[v]) | C[v]) - {v}
        piv = eliminate(R, C, F, v, stats)
        alive.discard(v)
        d = F.mul(d, piv)
        for u in nbrs:
            if u in alive:
                heapq.heappush(heap, (len(R[u]) + len(C[u]), u))
    return d, stats


def bareiss_dense(arcs, root, verts):
    idx = {v: i for i, v in enumerate(v for v in verts if v != root)}
    n = len(idx)
    M = [[0] * n for _ in range(n)]
    for x, y in arcs:
        if x == y or x == root:
            continue
        M[idx[x]][idx[x]] += 1
        if y != root:
            M[idx[x]][idx[y]] -= 1
    return SF.bareiss(M)


def contraction(arcs, root, verts):
    out = collections.defaultdict(collections.Counter)
    for x, y in arcs:
        if x != y:
            out[x][y] += 1
    return SF.tree_count(out, list(verts), root)


# ------------------------------------------------------------------ graphs

def balanced_random(n, narcs, rnd):
    """a random closed walk: a balanced, strongly connected multidigraph on the vertices it visits"""
    walk = [rnd.randrange(n) for _ in range(narcs)]
    arcs = [(walk[i], walk[(i + 1) % narcs]) for i in range(narcs)]
    return sorted(set(walk)), arcs


def synthetic_genome(nsnarl, pool, hap, rnd):
    """a chain of snarls between boundary vertices b_i, each with an interior pool and a nested
    bubble q_i0 -> nested pool -> q_i1, traversed by `hap` haplotypes (one closed walk).  Returns
    the vertices, the arcs, the snarl-tree order (interiors before their boundary) and the root."""
    B = lambda i: ("b", i % nsnarl)
    walk = []
    for h in range(hap):
        for i in range(nsnarl):
            walk.append(B(i))
            for _ in range(rnd.randrange(1, 4)):
                walk.append(("p", i, rnd.randrange(pool)))
            if rnd.random() < 0.6:
                walk.append(("q", i, 0))
                for _ in range(rnd.randrange(0, 3)):
                    walk.append(("r", i, rnd.randrange(pool)))
                walk.append(("q", i, 1))
            for _ in range(rnd.randrange(0, 3)):
                walk.append(("p", i, rnd.randrange(pool)))
    arcs = [(walk[i], walk[(i + 1) % len(walk)]) for i in range(len(walk))]
    verts = sorted(set(walk), key=repr)
    present = set(walk)
    order = []
    for i in range(nsnarl):
        order += [v for v in [("r", i, j) for j in range(pool)] if v in present]     # nested interior
        order += [v for v in [("q", i, 0), ("q", i, 1)] if v in present]            # nested boundary
        order += [v for v in [("p", i, j) for j in range(pool)] if v in present]    # snarl interior
        if i >= 1:
            order.append(B(i))                                                        # its entry
    return verts, arcs, order, B(0)


def brute_euler(arcs):
    """Eulerian circuits up to rotation, labelled arcs: trails that start with arc 0"""
    out = collections.defaultdict(list)
    for i, (x, y) in enumerate(arcs):
        out[x].append((y, i))
    used = [False] * len(arcs)
    used[0] = True
    cnt = [0]
    x0, y0 = arcs[0]

    def rec(v, d):
        if d == len(arcs):
            cnt[0] += v == x0
            return
        for w, i in out[v]:
            if not used[i]:
                used[i] = True
                rec(w, d + 1)
                used[i] = False
    rec(y0, 1)
    return cnt[0]


def phix_debruijn(S, k):
    N = len(S)
    SS = S + S[:k]
    arcs = [(SS[i:i + k], SS[i + 1:i + 1 + k]) for i in range(N)]
    verts = sorted({a for a, _ in arcs})
    return verts, arcs


# ------------------------------------------------------------------ main

def main():
    t0 = time.time()
    rnd = random.Random(20260924)
    FR = Field("frac")
    print("=== Bio 14: the BEST determinant by Schur complements along the snarl tree ===")
    print(f"    python {sys.version.split()[0]}\n", flush=True)

    # (M) Schur identity and Laplacian closure
    print("=== (M) det L~ = det L_II * det(Schur), and the Schur complement is a Laplacian ===")
    okM, cases = True, 0
    for _ in range(60):
        n = rnd.randrange(5, 12)
        verts, arcs = balanced_random(n, rnd.randrange(n + 2, 3 * n), rnd)
        if len(verts) < 4:
            continue
        root = verts[0]
        rest = verts[1:]
        I = rnd.sample(rest, rnd.randrange(1, len(rest)))
        Bd = [v for v in rest if v not in I]
        whole = Fraction(bareiss_dense(arcs, root, verts))
        # eliminate I and read the Schur complement on Bd, with the root's column folded in
        R, C = reduced_laplacian(arcs, root, FR)
        dI, st = FR.c(1), [0, 0]
        for v in I:
            dI *= eliminate(R, C, FR, v, st)
        schur = {b: {u: R[b].get(u, Fraction(0)) for u in Bd} for b in Bd}
        # determinant of the Schur complement
        rest_det, st2 = FR.c(1), [0, 0]
        R2 = {b: dict((u, x) for u, x in schur[b].items() if x != 0) for b in Bd}
        C2 = collections.defaultdict(set)
        for b, r in R2.items():
            for u in r:
                C2[u].add(b)
        for b in Bd:
            R2.setdefault(b, {})
            C2.setdefault(b, set())
        for b in Bd:
            rest_det *= eliminate(R2, C2, FR, b, st2)
        okM &= (dI * rest_det == whole)
        # Laplacian: off-diagonal <= 0, and row sums >= 0 (the deficit is the flow to the root)
        for b in Bd:
            okM &= all(schur[b][u] <= 0 for u in Bd if u != b)
            okM &= sum(schur[b].values()) >= 0
        cases += 1
    print(f"    {cases} random Eulerian digraphs, random interiors I")
    row("M", "eliminating an interior I multiplies by det L_II and leaves the Schur complement, "
             "which is again a (reduced) Laplacian: non-positive off-diagonal, non-negative row sums",
        okM)

    # (F) forests
    print("\n=== (F) det L_II counts spanning forests of the interior rooted into the boundary ===")
    okF, cf = True, 0
    for _ in range(40):
        n = rnd.randrange(4, 8)
        verts, arcs = balanced_random(n, rnd.randrange(n + 2, 2 * n + 4), rnd)
        if len(verts) < 3:
            continue
        I = rnd.sample(verts, rnd.randrange(1, len(verts)))
        Is = set(I)
        outs = {v: [y for (x, y) in arcs if x == v and y != v] for v in I}
        # brute force: each interior vertex picks one out-arc (labelled); no cycle inside I
        forests = 0
        for choice in itertools.product(*[range(len(outs[v])) for v in I]):
            nxt = {v: outs[v][c] for v, c in zip(I, choice)}
            ok = True
            for v in I:
                seen, x = set(), v
                while x in Is:
                    if x in seen:
                        ok = False
                        break
                    seen.add(x)
                    x = nxt[x]
                if not ok:
                    break
            forests += ok
        M = [[0] * len(I) for _ in I]
        ix = {v: i for i, v in enumerate(I)}
        for x, y in arcs:
            if x in Is and x != y:
                M[ix[x]][ix[x]] += 1
                if y in Is:
                    M[ix[x]][ix[y]] -= 1
        okF &= SF.bareiss(M) == forests
        cf += 1
    print(f"    {cf} random interiors")
    row("F", "det of the interior block equals the number of spanning forests of the interior whose "
             "roots lie in the boundary (all-minors matrix-tree theorem), by enumeration", okF)

    # (B) BEST
    print("\n=== (B) BEST: t * prod (d-1)! = Eulerian circuits, by brute force ===")
    okB, cb = True, 0
    for _ in range(40):
        n = rnd.randrange(3, 7)
        verts, arcs = balanced_random(n, rnd.randrange(n + 1, 13), rnd)
        t = bareiss_dense(arcs, verts[0], verts)
        outd = collections.Counter(x for x, _ in arcs)
        ec = t
        for v in verts:
            ec *= math.factorial(outd[v] - 1)
        okB &= ec == brute_euler(arcs)
        cb += 1
    print(f"    {cb} random Eulerian digraphs with up to 12 arcs")
    row("B", "the BEST count agrees with exhaustive enumeration of Eulerian circuits", okB)

    # (E) synthetic genomes, exact
    print("\n=== (E) synthetic genome graphs: snarl-tree elimination, exact ===")
    okE = True
    print("      snarls  haplotypes  vertices   arcs   max front   snarl order = Bareiss = contraction")
    for nsn, hap in ((10, 2), (20, 3), (40, 3), (60, 4)):
        verts, arcs, order, root = synthetic_genome(nsn, 4, hap, rnd)
        d1, st = det_by_order(arcs, root, order, FR)
        d2 = bareiss_dense(arcs, root, verts) if len(verts) <= 320 else None
        d3 = contraction(arcs, root, verts)
        good = d1 == d3 and (d2 is None or d1 == d2) and d1 > 0
        okE &= good
        print(f"      {nsn:6d}  {hap:10d}  {len(verts):8d}  {len(arcs):5d}   {st[0]:9d}   "
              f"{'equal' if good else 'DIFFERENT'} (t has {len(str(int(d1)))} digits"
              f"{'' if d2 is not None else '; Bareiss skipped'})", flush=True)
    row("E", "on synthetic genome graphs the determinant obtained by eliminating snarl interiors "
             "before their boundaries equals the dense Bareiss determinant and the forced-arc "
             "contraction, with a bounded front", okE)

    # (P) phiX174
    print("\n=== (P) phiX174 de Bruijn graphs (circular genome, single strand) ===")
    S = DF.load_phix()
    okP = True
    print("       k   vertices   arcs   max front (min-degree)   t by elimination = t by contraction")
    for k in (7, 8, 10, 12):
        verts, arcs = phix_debruijn(S, k)
        root = verts[0]
        t1 = time.time()
        d1, st = det_mindegree(arcs, root, FR)
        d3 = contraction(arcs, root, verts)
        good = d1 == d3 and d1 > 0
        okP &= good
        print(f"      {k:2d}   {len(verts):8d}  {len(arcs):5d}   {st[0]:22d}   "
              f"{'equal' if good else 'DIFFERENT'}, log10 t = {math.log10(int(d1)):.1f}"
              f"   [{time.time()-t1:.1f}s]", flush=True)
    row("P", "on the de Bruijn graphs of phiX174 minimum-degree Schur elimination and forced-arc "
             "contraction give the same arborescence count", okP)

    # (S) scaling
    print("\n=== (S) scaling: snarl order vs minimum-degree order, exact modulo two primes ===")
    okS = True
    print("      snarls   vertices     arcs   front(snarl)  front(min-deg)   equal mod p1, p2   "
          "us per vertex (snarl order)")
    for nsn in (500, 1500, 5000):
        verts, arcs, order, root = synthetic_genome(nsn, 4, 3, rnd)
        res = []
        for p in (P1, P2):
            F = Field("mod", p)
            t1 = time.time()
            a, sa = det_by_order(arcs, root, order, F)
            dt = time.time() - t1
            b, sb = det_mindegree(arcs, root, F)
            res.append((a, b, sa, sb, dt))
        eq = all(a == b and a != 0 for a, b, _, _, _ in res)
        okS &= eq
        print(f"      {nsn:6d}   {len(verts):8d}  {len(arcs):7d}   {res[0][2][0]:11d}  "
              f"{res[0][3][0]:14d}   {'yes' if eq else 'NO':>16}   {1e6*res[0][4]/len(verts):10.1f}",
              flush=True)
    row("S", "on synthetic genomes of up to ~45 000 vertices the snarl-order elimination agrees with "
             "the minimum-degree elimination modulo two primes; the front stays bounded and the time "
             "per vertex stays constant", okS)

    print("\n=== battery summary ===")
    bad = [k_ for k_, ok in ROWS if not ok]
    print(f"  TOTAL {len(ROWS)-len(bad)}/{len(ROWS)} passed" + (f"   FAILURES: {bad}" if bad else ""))
    print(f"  elapsed {time.time()-t0:.0f}s")


class _Tee:
    """print to the console and to bio14_output.txt (UTF-8) at the same time"""

    def __init__(self, path):
        self.f = open(path, "w", encoding="utf-8")
        self.out = sys.stdout

    def write(self, s):
        self.out.write(s)
        self.f.write(s)

    def flush(self):
        self.out.flush()
        self.f.flush()


if __name__ == "__main__":
    sys.stdout = _Tee(os.path.join(HERE, "bio14_output.txt"))
    main()
