# -*- coding: utf-8 -*-
"""
isoform_polytope.py -- the isoform decomposition polytope of a splicing DAG   (Bio 11, exploratory)

THE OBJECT.  A splicing graph is a DAG G with source s (transcription start) and sink t
(polyadenylation site); an isoform is an s-t path; junction/exon read counts give an integer
flow w >= 0 on the arcs with Kirchhoff conservation at the interior vertices and value F.
The decompositions of w into isoforms form the fibre polytope

      Theta(w) = { theta in R^{P(w)}_{>=0} :  sum_j theta_j 1_{P_j} = w },

P(w) the s-t paths inside supp(w).  It is a polytope (every theta_j <= F).

WHAT IS TRUE.
 (R)  The path-arc matrix of the paths inside supp(w) has rank |E_w| - |V_w| + 2: the flow
      space of a DAG in which every arc lies on an s-t path.
 (D)  dim Theta(w) = |P(w)| - (|E_w| - |V_w| + 2), because w is in the relative interior of the
      cone spanned by P(w).  Deconvolution is unique iff this is 0.  For a series of
      alternative events with k_1..k_r branches: dim = prod k_i - sum (k_i - 1) - 1, which is 0
      iff at most one event has k_i >= 2.  Dscam (12, 48, 33, 2): 38016 isoforms, rank 92,
      dim 37924.
 (V)  Every minimum-support decomposition (MFD, real or integer weights) is a vertex of
      Theta(w) -- a dependent support could be shrunk -- but vertices of Theta(w) need not be
      integral for integer w: three parallel-arc pairs in series give a vertex with four
      weights 1/2 while the MFD has two integer paths.  So "MFD = integer points that are
      convex combinations of extreme rays" (the brief) is not the right statement; extreme
      rays of the flow cone are the paths, vertices of Theta(w) are the independent-support
      decompositions, and integrality is a separate matter.
 (B)  BEST identity.  Add F return arcs t -> s to the multigraph G_w; its per-sequence count
      N_seq (arborescences x prod (d+ - 1)! / prod w(a)! / F!) equals
            sum_{theta in Theta(w) cap Z^P}  (F-1)! / prod_j theta_j!,
      the number of cyclic arrangements of F isoforms with multiplicities theta (weighted
      1/d as in Bio 8/9).  A determinant thus computes the multinomial-weighted partition
      function of the integer decompositions; the plain number of integer points does not
      have such a formula.
 (I)  Real MFD vs integer MFD on random small instances: reported (searched, see log).

Run: python -u isoform_polytope.py           (about a minute; the Dscam rank dominates)
"""

import sys, os, time, math, itertools, random, collections
from fractions import Fraction
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)
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

sys.path.insert(0, _sib("lib"))
sys.path.insert(0, _sib("Bio_08_DS_Lattice_Sum"))
import pgl3_building as pb
from ds_lattice import row, ROWS


def det(M):
    return 1 if not M else pb.bareiss_det(M)


# ------------------------------------------------------------------------------- DAGs

class DAG:
    def __init__(self, name, nv, arcs, s, t):
        """arcs: list of (u, v); vertices 0..nv-1 in topological order"""
        self.name, self.nv, self.arcs, self.s, self.t = name, nv, list(arcs), s, t
        self.nA = len(self.arcs)
        self.out = collections.defaultdict(list)
        for i, (u, v) in enumerate(self.arcs):
            assert u < v
            self.out[u].append(i)
        self._paths = None

    def paths(self):
        if self._paths is None:
            res = []
            stack = [(self.s, [])]
            while stack:
                v, p = stack.pop()
                if v == self.t:
                    res.append(tuple(p))
                    continue
                for a in self.out[v]:
                    stack.append((self.arcs[a][1], p + [a]))
            self._paths = sorted(res)
        return self._paths

    def n_paths(self):
        """by dynamic programming, for graphs too large to list"""
        cnt = [0] * self.nv
        cnt[self.t] = 1
        for v in range(self.nv - 1, -1, -1):
            if v != self.t:
                cnt[v] = sum(cnt[self.arcs[a][1]] for a in self.out[v])
        return cnt[self.s]

    def vec(self, path):
        x = [0] * self.nA
        for a in path:
            x[a] = 1
        return x

    def flow(self, theta):
        w = [0] * self.nA
        for p, th in theta.items():
            for a in p:
                w[a] += th
        return w


def series(name, blocks):
    """series composition of alternative events; a block is the list of branch lengths
    (number of arcs) between consecutive junction vertices"""
    arcs = []
    nv = 1
    cur = 0
    for block in blocks:
        nxt = None
        internal = []
        for ln in block:
            chain = [cur]
            for _ in range(ln - 1):
                chain.append(nv)
                nv += 1
            internal.append(chain)
        nxt = nv
        nv += 1
        for chain in internal:
            chain.append(nxt)
            for i in range(len(chain) - 1):
                arcs.append((chain[i], chain[i + 1]))
        cur = nxt
    return DAG(name, nv, sorted(arcs), 0, cur)


def parallel_pairs(name, n):
    """n consecutive vertices joined by two parallel arcs each: 2^n paths"""
    arcs = []
    for i in range(n):
        arcs += [(i, i + 1), (i, i + 1)]
    return DAG(name, n + 1, arcs, 0, n)


# ------------------------------------------------------------------- linear algebra

def rank_frac(rows):
    A = [[Fraction(x) for x in r] for r in rows]
    m = len(A)
    n = len(A[0]) if m else 0
    r = 0
    for c in range(n):
        piv = next((i for i in range(r, m) if A[i][c] != 0), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = 1 / A[r][c]
        A[r] = [x * inv for x in A[r]]
        for i in range(m):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        r += 1
    return r


def rank_mod(B, p=2147483647):
    A = np.array(B, dtype=np.int64) % p
    nr, nc = A.shape
    r = 0
    for c in range(nc):
        if r == nr:
            break
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0:
            continue
        piv = r + int(nz[0])
        if piv != r:
            A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        fac = A[:, c].copy()
        fac[r] = 0
        A = (A - np.outer(fac, A[r])) % p
        r += 1
    return r


def solve_exact(cols, w):
    """theta with sum_j theta_j cols[j] = w, if the columns are independent and the system is
    consistent; else None"""
    m, r = len(w), len(cols)
    A = [[Fraction(cols[j][i]) for j in range(r)] + [Fraction(w[i])] for i in range(m)]
    piv_cols = []
    rr = 0
    for c in range(r):
        piv = next((i for i in range(rr, m) if A[i][c] != 0), None)
        if piv is None:
            return None
        A[rr], A[piv] = A[piv], A[rr]
        inv = 1 / A[rr][c]
        A[rr] = [x * inv for x in A[rr]]
        for i in range(m):
            if i != rr and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[rr])]
        piv_cols.append(c)
        rr += 1
    for i in range(rr, m):
        if A[i][r] != 0:
            return None
    return [A[i][r] for i in range(r)]


# ------------------------------------------------------------------- the polytope

def support_paths(G, w):
    return [p for p in G.paths() if all(w[a] > 0 for a in p)]


def integer_points(G, w, paths, cap=2000000):
    """all integer theta >= 0 on `paths` with sum theta_j 1_{P_j} = w"""
    out = []
    rem = list(w)
    theta = [0] * len(paths)
    nodes = [0]

    def rec(j):
        nodes[0] += 1
        if nodes[0] > cap:
            raise RuntimeError("cap")
        if j == len(paths):
            if all(x == 0 for x in rem):
                out.append(tuple(theta))
            return
        p = paths[j]
        mx = min(rem[a] for a in p)
        for th in range(mx, -1, -1):
            for a in p:
                rem[a] -= th
            theta[j] = th
            rec(j + 1)
            for a in p:
                rem[a] += th
        theta[j] = 0

    rec(0)
    return out


def vertices(G, w, paths, r):
    """vertices of Theta(w): feasible theta with linearly independent support of size <= r"""
    cols = [G.vec(p) for p in paths]
    verts = set()
    for k in range(1, r + 1):
        for S in itertools.combinations(range(len(paths)), k):
            th = solve_exact([cols[j] for j in S], w)
            if th is None or any(x < 0 for x in th):
                continue
            if any(x == 0 for x in th):
                continue                       # not a basic solution with this support
            v = [Fraction(0)] * len(paths)
            for j, x in zip(S, th):
                v[j] = x
            verts.add(tuple(v))
    return sorted(verts)


def nseq_augmented(G, w, F):
    """N_seq of G_w plus F return arcs t -> s"""
    L = collections.defaultdict(int)
    outd = collections.Counter()
    for a, (u, v) in enumerate(G.arcs):
        if w[a]:
            L[(u, u)] += w[a]
            L[(u, v)] -= w[a]
            outd[u] += w[a]
    L[(G.t, G.t)] += F
    L[(G.t, G.s)] -= F
    outd[G.t] += F
    supp = sorted(outd)
    pos = {v: i for i, v in enumerate(supp[1:])}
    n = len(supp) - 1
    M = [[0] * n for _ in range(n)]
    for (u, v), c in L.items():
        if u in pos and v in pos:
            M[pos[u]][pos[v]] += c
        elif u in pos and u == v:
            M[pos[u]][pos[u]] += c
    t = det(M)
    loc = 1
    for v in supp:
        loc *= math.factorial(outd[v] - 1)
    denom = math.factorial(F)
    for a in range(G.nA):
        denom *= math.factorial(w[a])
    return Fraction(t * loc, denom), t


def analyse(G, w=None, theta_true=None, small=True, verbose=True):
    P = G.paths() if small else None
    nP = len(P) if small else G.n_paths()
    m, n = G.nA, G.nv
    if verbose:
        print(f"\n=== {G.name}: |V| = {n}, |E| = {m}, s-t paths {nP}, |E|-|V|+2 = {m - n + 2} ===")
    res = dict(nP=nP, mn2=m - n + 2)
    if small:
        cols = [G.vec(p) for p in P]
        r = rank_frac([list(c) for c in zip(*cols)]) if cols else 0
        res["rank"] = r
        res["dim"] = nP - r
        if verbose:
            print(f"    rank of the path-arc matrix {r} (= |E|-|V|+2: {r == m - n + 2});"
                  f" dim Theta(w) for full-support w: {nP - r}"
                  f"{'  -> deconvolution UNIQUE' if nP - r == 0 else ''}")
    else:
        # rank modulo a large prime on the full path matrix, built column by column
        return res
    if w is None:
        return res
    F = sum(w[a] for a in G.out[G.s])
    Pw = support_paths(G, w)
    res["F"] = F
    pts = integer_points(G, w, Pw)
    res["npts"] = len(pts)
    rw = rank_frac([list(c) for c in zip(*[G.vec(p) for p in Pw])])
    verts = vertices(G, w, Pw, rw)
    frac = [v for v in verts if any(x.denominator != 1 for x in v)]
    res["nverts"], res["nfrac"] = len(verts), len(frac)
    mfd_int = min(sum(1 for x in th if x) for th in pts)
    mfd_real = min(sum(1 for x in v if x) for v in verts)
    res["mfd_int"], res["mfd_real"] = mfd_int, mfd_real
    # every integer minimum-support decomposition is a vertex
    mins = [tuple(Fraction(x) for x in th) for th in pts if sum(1 for x in th if x) == mfd_int]
    vset = set(verts)
    res["mins_are_vertices"] = all(mv in vset for mv in mins)
    # BEST identity
    lhs, t = nseq_augmented(G, w, F)
    rhs = sum(Fraction(math.factorial(F - 1), math.prod(math.factorial(x) for x in th)) for th in pts)
    res["best"] = (lhs == rhs)
    if verbose:
        print(f"    w = {w}  (F = {F}); paths in supp(w): {len(Pw)}; rank there {rw}; dim Theta(w) = {len(Pw) - rw}")
        print(f"    integer points {len(pts)}; vertices {len(verts)} of which fractional {len(frac)};"
              f" MFD real {mfd_real}, integer {mfd_int}; all integer minima are vertices: {res['mins_are_vertices']}")
        if frac:
            v = frac[0]
            print(f"    a fractional vertex: " + ", ".join(f"{x}" for x in v))
        print(f"    BEST identity: N_seq(G_w + F returns) = {lhs}  vs  sum_theta (F-1)!/prod theta! = {rhs}"
              f"   {'ok' if lhs == rhs else 'MISMATCH'}   (arborescences {t})")
    return res


# ----------------------------------------------------------------------------- main

def main():
    t0 = time.time()
    random.seed(11)
    print("=== the isoform decomposition polytope of a splicing DAG ===")

    okR = okD = okV = okB = okU = True
    models = [
        ("cassette exon (skipping)", series("cassette exon", [[2, 1]])),
        ("mutually exclusive exon pair", series("mutually exclusive exons", [[2, 2]])),
        ("alternative 5' splice site", series("alternative 5' site", [[1, 1]])),
        ("two cassettes in series", series("two cassettes", [[2, 1], [2, 1]])),
        ("three cassettes in series", series("three cassettes", [[2, 1], [2, 1], [2, 1]])),
        ("tropomyosin-like: MXE, cassette, MXE, cassette, 4-way MXE",
         series("tropomyosin-like", [[2, 2], [2, 1], [2, 2], [2, 1], [2, 2, 2, 2]])),
        ("three parallel-arc pairs in series", parallel_pairs("parallel pairs x3", 3)),
    ]

    # ------------------------------------------------------------ (R), (D): rank and dimension
    print("\n--- (R)/(D) rank of the path-arc matrix and the dimension of Theta(w) ---")
    dims = {}
    for label, G in models:
        r = analyse(G, verbose=True)
        okR &= (r["rank"] == r["mn2"])
        dims[label] = r["dim"]
    # series formula: dim = prod k_i - sum (k_i - 1) - 1
    pred = {"cassette exon (skipping)": 2 - 1 - 1, "mutually exclusive exon pair": 0,
            "alternative 5' splice site": 0, "two cassettes in series": 4 - 2 - 1,
            "three cassettes in series": 8 - 3 - 1,
            "tropomyosin-like: MXE, cassette, MXE, cassette, 4-way MXE": 2 * 2 * 2 * 2 * 4 - (1 + 1 + 1 + 1 + 3) - 1,
            "three parallel-arc pairs in series": 8 - 3 - 1}
    okD &= all(dims[k] == v for k, v in pred.items())
    okU &= (dims["cassette exon (skipping)"] == 0 and dims["mutually exclusive exon pair"] == 0
            and dims["alternative 5' splice site"] == 0 and dims["two cassettes in series"] == 1)
    row("R", "the path-arc matrix of a splicing DAG has rank |E| - |V| + 2 (the flow space)",
        f"{len(models)} graphs", "exact", okR)
    row("D", "dim Theta(w) = #paths - (|E|-|V|+2); for events in series prod k_i - sum(k_i - 1) - 1",
        f"{len(models)} graphs", "exact", okD)
    row("U", "deconvolution from junction counts is unique for a single alternative event and "
             "already non-unique (dim 1) for two cassette exons in series", "4 graphs", "exact", okU)

    # ------------------------------------------------------------ Dscam
    print("\n--- Dscam: four blocks of mutually exclusive exons, 12 x 48 x 33 x 2 ---")
    D = series("Dscam", [[2] * 12, [2] * 48, [2] * 33, [2] * 2])
    nP = D.n_paths()
    print(f"    |V| = {D.nv}, |E| = {D.nA}, isoforms {nP}, |E|-|V|+2 = {D.nA - D.nv + 2}")
    t1 = time.time()
    # path matrix column by column (38016 columns) via the product structure
    blocks = [[2] * 12, [2] * 48, [2] * 33, [2] * 2]
    P = D.paths()
    assert len(P) == nP
    M = np.zeros((D.nA, nP), dtype=np.int64)
    for j, p in enumerate(P):
        M[list(p), j] = 1
    rD = rank_mod(M)
    print(f"    rank of the {D.nA} x {nP} path-arc matrix modulo 2^31-1: {rD}   [{time.time()-t1:.0f}s]")
    print(f"    dim of the isoform polytope for a full-support flow: {nP - rD}")
    okR &= (rD == D.nA - D.nv + 2)
    row("S", f"Dscam: {nP} isoforms, junction counts determine a {nP - rD}-dimensional polytope of "
             f"abundances (rank {rD} = |E|-|V|+2)", "1 graph", "exact mod p", rD == D.nA - D.nv + 2
        and nP - rD == 38016 - 92)

    # ------------------------------------------------------------ (V), (B): flows
    print("\n--- (V)/(B) integer points, vertices, MFD and the BEST identity ---")
    # two cassettes: the classical unidentifiable pair
    G = models[3][1]
    P = G.paths()
    w = G.flow({P[0]: 2, P[3]: 3, P[1]: 1})          # inclusion-inclusion, skip-skip, mixed
    r = analyse(G, w)
    okV &= r["mins_are_vertices"]
    okB &= r["best"]
    # three parallel pairs: the fractional vertex
    G = models[6][1]
    w = [1] * G.nA
    r = analyse(G, w)
    okV &= r["mins_are_vertices"] and r["nfrac"] > 0 and r["mfd_int"] == 2 and r["mfd_real"] == 2
    okB &= r["best"]
    # tropomyosin-like with a sparse true isoform set
    G = models[5][1]
    P = G.paths()
    chosen = random.sample(P, 5)
    w = G.flow({p: random.randint(1, 3) for p in chosen})
    r = analyse(G, w)
    okV &= r["mins_are_vertices"]
    okB &= r["best"]
    # three cassettes, random flows
    G = models[4][1]
    P = G.paths()
    for trial in range(3):
        chosen = random.sample(P, random.randint(2, 4))
        w = G.flow({p: random.randint(1, 3) for p in chosen})
        r = analyse(G, w, verbose=(trial == 0))
        okV &= r["mins_are_vertices"]
        okB &= r["best"]
    row("V", "every minimum-support decomposition is a vertex of Theta(w); vertices can be fractional "
             "for integer w (three parallel pairs: four weights 1/2 beside an integer MFD of two paths)",
        "6 flows", "exact", okV)
    row("B", "BEST identity: N_seq of G_w with F return arcs equals sum over integer decompositions of "
             "(F-1)!/prod theta_j!", "6 flows", "exact", okB)

    # ------------------------------------------------------------ (I): real vs integer MFD
    print("\n--- (I) real MFD vs integer MFD, random small DAGs ---")
    found = None
    n_inst = n_frac = 0
    for inst in range(400):
        nv = random.randint(4, 6)
        arcs = []
        for u in range(nv - 1):
            for v in range(u + 1, nv):
                for _ in range(random.choice([0, 1, 1, 2]) if v - u <= 2 else random.choice([0, 0, 1])):
                    arcs.append((u, v))
        try:
            G = DAG(f"random {inst}", nv, sorted(arcs), 0, nv - 1)
        except AssertionError:
            continue
        P = G.paths()
        if not (3 <= len(P) <= 12):
            continue
        chosen = random.sample(P, random.randint(2, min(4, len(P))))
        w = G.flow({p: random.randint(1, 3) for p in chosen})
        try:
            r = analyse(G, w, verbose=False)
        except RuntimeError:
            continue
        n_inst += 1
        n_frac += (r["nfrac"] > 0)
        okV &= r["mins_are_vertices"]
        okB &= r["best"]
        if r["mfd_real"] < r["mfd_int"] and found is None:
            found = (G, w, r)
    print(f"    {n_inst} instances; {n_frac} with a fractional vertex; BEST identity and vertex property held "
          f"on all: {okB and okV}")
    if found:
        G, w, r = found
        print(f"    instance with real MFD {r['mfd_real']} < integer MFD {r['mfd_int']}: arcs {G.arcs}, w = {w}")
    else:
        print(f"    no instance with real MFD < integer MFD among these (integer weights were never worse)")
    row("I", "on random small flows the minimum-support real and integer decompositions were compared "
             "(result in the log); the vertex property and the BEST identity held throughout",
        f"{n_inst} instances", "exact", okB and okV)

    print("\n=== battery summary ===")
    bad = [x for x in ROWS if not x[4]]
    print(f"  TOTAL {sum(1 for x in ROWS if x[4])}/{len(ROWS)} passed"
          + ("" if not bad else f"   FAILURES: {[(x[0], x[2]) for x in bad]}"))
    print(f"  elapsed {time.time()-t0:.0f}s")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
