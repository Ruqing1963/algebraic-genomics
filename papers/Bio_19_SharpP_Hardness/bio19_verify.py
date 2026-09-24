# -*- coding: utf-8 -*-
r"""
bio19_verify.py -- #DS is #P-hard: undirected Eulerian circuits inside genomes.

The reduction.  Let G be a connected loopless multigraph with all degrees 2 or 4 and at least three
edges, and let eps(G) be its number of Eulerian circuits (cyclic edge sequences, up to rotation and
reversal).  Counting eps(G) is #P-complete (Brightwell-Winkler 2005), already for 4-regular graphs
(Ge-Stefankovic 2012).  From G we build, in polynomial time, a circular DNA sequence S and an even k:

  * each vertex x of G becomes a distinct PALINDROMIC k-mer P_x = u_x rc(u_x);
  * S is spelled along one Eulerian circuit of G: P_x0 Sigma_e1 P_x1 Sigma_e2 P_x2 ...  with a
    spacer Sigma_e for each edge whose first and last letters are ASSIGNED so that the deg(x) <= 4
    letters following the occurrences of P_x in S + rc(S) are pairwise different (otherwise the
    (k+1)-mers P_x c repeat); the interior of each spacer is random, and the result is checked.

Then the double cover D_k(S) has exactly the P_x as vertices (all palindromic), one rho-orbit
{P_x Sigma_e P_y, P_y rc(Sigma_e) P_x} of multiplicity m = 1 per edge, and its quotient by rho is G.
A lattice point of Bio 8's box chooses one arc per orbit, i.e. an ORIENTATION of G; it is in the box iff
the orientation is Eulerian, and it carries N_seq = ec(orientation) reconstructions (BEST).  Hence

    #DS(S, k) = sum over Eulerian orientations O of ec(O) = 2 eps(G),

every directed traversal of every undirected Eulerian circuit being counted once.  A polynomial-time
(or "meshless", closed-form) evaluation of #DS would count Eulerian circuits of 4-regular graphs.
As by-products: the number of box points (Bio 12's first column) is the number of Eulerian
orientations of G, itself #P-complete (Mihail-Winkler 1996); and in this class Bio 17's generation
conjecture is Kotzig's theorem (inversions at palindromic vertices = kappa-transformations).

Checks (all exact):
  (E) the identity sum_O ec(O) = 2 eps(G) on random Eulerian multigraphs and 4-regular graphs,
      eps by brute-force enumeration of trails, ec by the BEST theorem;
  (R) the realization: D_k(S) built from the definition (ds_fermion.Cover) has vertex set {P_x},
      all palindromic, one orbit of multiplicity 1 per edge, no self-complementary content, and
      quotient graph exactly G;
  (D) #DS(S,k) by Bio 8's lattice sum equals 2 eps(G), and the number of box points equals the
      number of Eulerian orientations;
  (K) the closure of S under Bio 17's moves is all of DS(S,k) (Kotzig);
  (S) the size of the reduction: |S| is linear in |E(G)|, while #DS grows exponentially.

    python -u bio19_verify.py

Standard library only; imports ds_fermion.py (Bio 13) and bio17_verify.py (Bio 17).
"""
import os, sys, time, math, random, collections
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "Bio_13_Fermionic_Partition"))
sys.path.insert(0, os.path.join(HERE, "..", "Bio_17_Inversions_Fibre"))
import ds_fermion as DF                               # noqa: E402
from ds_fermion import Cover, rc, enumerate_box        # noqa: E402
import bio17_verify as B17                             # noqa: E402

ROWS = []


def row(key, desc, ok):
    ROWS.append((key, bool(ok)))
    print(("  [ok]   " if ok else "  [FAIL] ") + f"({key}) {desc}", flush=True)


# ------------------------------------------------------------------ graphs

def connected(n, edges):
    adj = collections.defaultdict(set)
    for x, y in edges:
        adj[x].add(y)
        adj[y].add(x)
    seen, st = {0}, [0]
    while st:
        v = st.pop()
        for w in adj[v]:
            if w not in seen:
                seen.add(w)
                st.append(w)
    return len(seen) == n


def random_4regular(n, rnd, simple=False):
    while True:
        stubs = [v for v in range(n) for _ in range(4)]
        rnd.shuffle(stubs)
        edges = [(stubs[2 * i], stubs[2 * i + 1]) for i in range(2 * n)]
        if any(x == y for x, y in edges):
            continue
        if simple and len({tuple(sorted(e)) for e in edges}) < len(edges):
            continue
        if connected(n, edges):
            return edges


def random_eulerian(n, rnd, ncyc):
    """a union of random cycles (even degrees), loopless and connected"""
    while True:
        edges = []
        for _ in range(ncyc):
            L = rnd.randrange(2, n + 1)
            cyc = rnd.sample(range(n), L)
            for i in range(L):
                edges.append((cyc[i], cyc[(i + 1) % L]))
        deg = collections.Counter(v for e in edges for v in e)
        if len(edges) >= 3 and all(x != y for x, y in edges) and connected(n, edges) \
                and len(deg) == n and max(deg.values()) <= 4:
            return edges


def eps_bruteforce(n, edges):
    """undirected Eulerian circuits up to rotation and reversal = the trails that start with edge 0
    traversed as x0 -> y0 and return to x0"""
    adj = collections.defaultdict(list)
    for i, (x, y) in enumerate(edges):
        adj[x].append((y, i))
        adj[y].append((x, i))
    E = len(edges)
    used = [False] * E
    x0, y0 = edges[0]
    used[0] = True
    cnt = [0]

    def rec(v, d):
        if d == E:
            cnt[0] += v == x0
            return
        for w, i in adj[v]:
            if not used[i]:
                used[i] = True
                rec(w, d + 1)
                used[i] = False
    rec(y0, 1)
    return cnt[0]


def ec_best(n, arcs):
    """Eulerian circuits (up to rotation) of a balanced connected digraph with labelled arcs, BEST"""
    outd = collections.Counter(x for x, _ in arcs)
    L = [[Fraction(0)] * n for _ in range(n)]
    for x, y in arcs:
        if x != y:
            L[x][x] += 1
            L[x][y] -= 1
    t = DF.det([[L[r][s] for s in range(1, n)] for r in range(1, n)]) if n > 1 else Fraction(1)
    val = int(t)
    for v in range(n):
        val *= math.factorial(outd[v] - 1)
    return val


def sum_orientations(n, edges):
    E = len(edges)
    tot, neo = 0, 0
    for bits in range(1 << E):
        arcs = [(x, y) if not (bits >> i) & 1 else (y, x) for i, (x, y) in enumerate(edges)]
        bal = collections.Counter()
        for x, y in arcs:
            bal[x] += 1
            bal[y] -= 1
        if any(bal.values()):
            continue
        neo += 1
        tot += ec_best(n, arcs)
    return tot, neo


def euler_circuit(n, edges):
    """one undirected Eulerian circuit (Hierholzer): list of (x, edge index, y)"""
    adj = collections.defaultdict(list)
    for i, (x, y) in enumerate(edges):
        adj[x].append((y, i))
        adj[y].append((x, i))
    used = [False] * len(edges)
    stack, out = [(0, None, None)], []
    while stack:
        v, e, u = stack[-1]
        while adj[v] and used[adj[v][-1][1]]:
            adj[v].pop()
        if adj[v]:
            w, i = adj[v].pop()
            used[i] = True
            stack.append((w, i, v))
        else:
            stack.pop()
            if e is not None:
                out.append((u, e, v))
    out.reverse()
    return out


# ------------------------------------------------------------------ the realization

def realize(n, edges, k, Ls, rnd, tries=200):
    """a circular sequence S whose double cover D_k(S) has quotient G (see the docstring);
    G must have maximum degree 4"""
    assert k % 2 == 0
    deg = collections.Counter(v for e in edges for v in e)
    assert max(deg.values()) <= 4, "the realization needs maximum degree 4"
    W = euler_circuit(n, edges)
    assert len(W) == len(edges) and all(W[i][2] == W[(i + 1) % len(W)][0] for i in range(len(W)))
    for _ in range(tries):
        P = []
        while len(P) < n:
            u = "".join(rnd.choice("ACGT") for _ in range(k // 2))
            p = u + rc(u)
            if p not in P:
                P.append(p)
        # The spacer of traversal i runs from P_{x_i} to P_{x_{i+1}}.  In S + rc(S) the palindrome
        # P_x occurs deg(x) times; the letters that follow it are the first letters of the outgoing
        # spacers and the complements of the last letters of the incoming ones.  They must be
        # pairwise different, or the (k+1)-mers P_x c repeat and create extra vertices.  With
        # deg(x) <= 4 they are made a subset of {A, C, G, T} deterministically: at the j-th visit
        # of x the outgoing spacer starts with "AC"[j] and the incoming one ends with "CA"[j]
        # (complement "GT"[j]).
        visit = collections.Counter()
        first, last = [None] * len(W), [None] * len(W)
        for i, (x, e, y) in enumerate(W):
            j = visit[x]
            visit[x] += 1
            first[i] = "AC"[j]
            last[(i - 1) % len(W)] = "CA"[j]
        sig = []
        for i in range(len(W)):
            while True:
                s = first[i] + "".join(rnd.choice("ACGT") for _ in range(Ls - 2)) + last[i]
                if s != rc(s):
                    sig.append(s)
                    break
        S = "".join(P[x] + sig[i] for i, (x, e, y) in enumerate(W))
        cov = Cover(S, k)
        ok = (sorted(cov.words) == sorted(P) and cov.n_fixV == n and len(cov.reps) == len(edges)
              and cov.n_fixA == 0 and all(cov.m[a] == 1 for a in range(cov.nA)))
        if not ok:
            continue
        idx = {p: i for i, p in enumerate(P)}
        quo = sorted(tuple(sorted((idx[cov.words[cov.tail[a]]], idx[cov.words[cov.head[a]]])))
                     for a in cov.reps)
        if quo == sorted(tuple(sorted(e)) for e in edges):
            return S, cov
    return None, None


def ds_lattice(cov):
    pts = enumerate_box(cov)
    return sum(cov.n_seq(cov.point(g)) for g in pts), len(pts)


# ------------------------------------------------------------------ main

def main():
    t0 = time.time()
    rnd = random.Random(20260924)
    print("=== Bio 19: #DS is #P-hard -- undirected Eulerian circuits inside genomes ===")
    print(f"    python {sys.version.split()[0]}\n", flush=True)

    graphs = []
    for n in (3, 4, 5, 6):
        for _ in range(3):
            graphs.append(("4-regular", n, random_4regular(n, rnd)))
    for n in (5, 6, 7):
        graphs.append(("4-regular simple", n, random_4regular(n, rnd, simple=True)))
    for _ in range(8):
        n = rnd.randrange(3, 7)
        g = random_eulerian(n, rnd, rnd.randrange(2, 4))
        if len(g) <= 14:
            graphs.append(("Eulerian", n, g))

    # (E) the combinatorial identity
    print("=== (E) sum over Eulerian orientations of ec(O) = 2 eps(G) ===")
    okE, data = True, []
    for kind, n, g in graphs:
        eps = eps_bruteforce(n, g)
        s, neo = sum_orientations(n, g)
        okE &= s == 2 * eps
        data.append((kind, n, g, eps, neo))
        print(f"    {kind:17s} n={n} |E|={len(g):2d}: eps(G) = {eps:6d}, Eulerian orientations {neo:5d}, "
              f"sum ec = {s:7d}   {'ok' if s == 2*eps else 'MISMATCH'}", flush=True)
    row("E", "summing BEST's ec over the Eulerian orientations of G counts each undirected Eulerian "
             "circuit twice (once per direction of traversal)", okE)

    # (R), (D), (K) the genome
    print("\n=== (R) realization, (D) #DS by Bio 8's lattice sum, (K) closure under Bio 17's moves ===")
    okR, okD, okK = True, True, True
    k, Ls = 12, 18
    for kind, n, g, eps, neo in data:
        S, cov = realize(n, g, k, Ls, rnd)
        if S is None:
            okR = False
            print(f"    {kind} n={n}: realization failed")
            continue
        t1 = time.time()
        ds, npts = ds_lattice(cov)
        okD &= ds == 2 * eps and npts == neo
        tk = B17.Tok(cov)
        st = tk.circuit_S()
        cl, bad = tk.closure(st, "IR")
        okK &= cl is not None and len(cl) == ds and bad == 0
        print(f"    {kind:17s} n={n} |E|={len(g):2d}: |S| = {len(S):4d} bp, k = {k}, D_k: {cov.nv} "
              f"palindromic vertices, {len(cov.reps)} orbits (m = 1); box points {npts:5d} "
              f"(= {neo}), #DS = {ds:7d} (= 2 x {eps}); closure {len(cl) if cl else '?':>7}"
              f"   [{time.time()-t1:.1f}s]", flush=True)
    row("R", "from every test graph G a circular sequence S is built whose double cover D_k(S) has "
             "exactly the palindromic words P_x as vertices, one orbit of multiplicity 1 per edge, no "
             "self-complementary content, and quotient G", okR)
    row("D", "#DS(S,k) computed by Bio 8's lattice sum equals 2 eps(G), and the box points of Bio 12 "
             "are exactly the Eulerian orientations of G", okD)
    row("K", "the closure of S under inversions at (palindromic) inverted repeats and reverse "
             "complementation is all of DS(S,k): Bio 17's conjecture holds on this class, by "
             "Kotzig's theorem", okK)

    # (S) size of the reduction
    print("\n=== (S) the reduction is linear, the count is exponential ===")
    okS = True
    print("      n  |E|   |S| (bp)   Eulerian orientations         #DS = 2 eps(G)     [s]")
    prev = None
    for n in (4, 5, 6, 7, 8):
        g = random_4regular(n, rnd, simple=(n >= 5))
        S, cov = realize(n, g, k, Ls, rnd)
        t1 = time.time()
        ds, npts = ds_lattice(cov)
        okS &= len(S) == len(g) * (k + Ls)
        print(f"     {n:2d}  {len(g):3d}   {len(S):7d}   {npts:21d}   {ds:22d}   {time.time()-t1:5.1f}",
              flush=True)
    row("S", "|S| = |E(G)| (k + |Sigma|) grows linearly with the graph", okS)

    print("\n=== battery summary ===")
    bad = [k_ for k_, ok in ROWS if not ok]
    print(f"  TOTAL {len(ROWS)-len(bad)}/{len(ROWS)} passed" + (f"   FAILURES: {bad}" if bad else ""))
    print(f"  elapsed {time.time()-t0:.0f}s")


class _Tee:
    """print to the console and to bio19_output.txt (UTF-8) at the same time"""

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
    sys.stdout = _Tee(os.path.join(HERE, "bio19_output.txt"))
    main()
