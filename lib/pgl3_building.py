"""
Paper X battery: A~_2 complexes, Hodge and Hecke Laplacians, the cubic companion, signed shadows
(pgl3_building.py).  Exact integer arithmetic (Python ints, own Smith diagonalization); numpy only for
floating spectra (checks marked "numeric").

Complexes:
  Thin (q = 1): triangulated tori R^2/L with the A~_2 tiling, L = 3Z^2 (9 vertices), 6Z^2 (36 vertices).
  Thick (q = 2): finite quotients of an A~_2 building of order 2 built from a triangle presentation over
  PG(2,2) compatible with the Singer cycle lambda: x -> x+1 (found by exhaustive search, axioms re-verified),
  as abelian voltage covers of the 3-vertex type quotient X_3 with H_1(X_3) = Z/2+Z/2+Z/14:
  Y_21 (Z/7), Y_24 ((Z/2)^3), Y_42 (Z/14), Y_168 (full).  Vertex links are verified to be Heawood graphs.

Checks (letters as in the paper's battery table):
 (P)  triangle-presentation axioms; simplicial; Heawood links; regularity q^2+q+1 / q+1; flag condition
 (D)  d_1 d_2 = 0; Betti numbers; Euler characteristic n - E + F
 (H)  Smith forms of Delta_0, Delta_1, Delta_2: coker = Z^{b_i} + Jac_i;  DKM critical group K(Y) = Z_1/im(d_2 d_2^t)
 (L)  Hecke-Laplacian Delta^H = (q^3-1)I + A_1 - qA_2: two-sided kernel Z1, equal cofactors,
      |Jac^H| = (1/n) prod_j |P_j(1)| over nontrivial joint eigenvalues (numeric), Jac^H != Jac_0
 (C)  I - C^(3) unimodularly equivalent to (-Delta^H) + I + I via explicit U, V
 (T)  Tr C^k, k <= 6, directly and via the Newton recursion; Tr C = Tr C^2 = 0, Tr C^3 = n(q-1)^2(q+1);
      Tr C^k = 0 for 3 not dividing k;  net traces tr_m >= 0 for m <= M_0 and a bound beyond
 (R)  spectrum of C^(3) (numeric): trivial 1, q, q^2 simple; max nontrivial modulus; Ramanujan test
 (S)  signed shadow of Delta^H and of Delta_1: explicit U (I - B^t) V = M + (-I); strong connectivity,
      double loops, B non-negative; Smith form check for the Delta^H shadow
"""
import itertools, math, sys
import numpy as np
import sympy as sp

# ============================================================ exact integer linear algebra
def diagonalize(M, transforms=False, coltransforms=False):
    """P M Q = D (diagonal, not divisibility-ordered).  Returns (d, P) with P unimodular (row ops);
    with coltransforms=True returns (d, P, Qinv) where Qinv is the inverse of the column transform."""
    A = [list(map(int, r)) for r in M]; m = len(A); n = len(A[0]) if m else 0
    P = [[int(i == j) for j in range(m)] for i in range(m)] if transforms else None
    Qi = [[int(i == j) for j in range(n)] for i in range(n)] if coltransforms else None
    Qc = [[int(i == j) for j in range(n)] for i in range(n)] if coltransforms else None   # Q itself (columns)
    def rowop(r, i, q):
        A[r] = [a - q*b for a, b in zip(A[r], A[i])]
        if transforms: P[r] = [a - q*b for a, b in zip(P[r], P[i])]
    def rowswap(r, i):
        A[r], A[i] = A[i], A[r]
        if transforms: P[r], P[i] = P[i], P[r]
    def colop(c, j, q):
        for t in range(m): A[t][c] -= q*A[t][j]
        if coltransforms:
            Qi[j] = [a + q*b for a, b in zip(Qi[j], Qi[c])]     # Q[:,c] -= q Q[:,j]  <=>  Qinv[j,:] += q Qinv[c,:]
            for t in range(n): Qc[t][c] -= q*Qc[t][j]
    def colswap(c, j):
        for t in range(m): A[t][c], A[t][j] = A[t][j], A[t][c]
        if coltransforms:
            Qi[c], Qi[j] = Qi[j], Qi[c]
            for t in range(n): Qc[t][c], Qc[t][j] = Qc[t][j], Qc[t][c]
    t = 0
    while t < min(m, n):
        piv = None
        for r in range(t, m):
            for c in range(t, n):
                if A[r][c] != 0 and (piv is None or abs(A[r][c]) < abs(A[piv[0]][piv[1]])): piv = (r, c)
        if piv is None: break
        rowswap(t, piv[0]); colswap(t, piv[1])
        while True:
            changed = False
            for r in range(t+1, m):
                if A[r][t] != 0:
                    q = A[r][t] // A[t][t]; rowop(r, t, q)
                    if A[r][t] != 0: rowswap(t, r); changed = True; break
            if changed: continue
            for c in range(t+1, n):
                if A[t][c] != 0:
                    q = A[t][c] // A[t][t]; colop(c, t, q)
                    if A[t][c] != 0: colswap(t, c); changed = True; break
            if not changed: break
        t += 1
    d = [abs(A[i][i]) for i in range(min(m, n))]
    if coltransforms: return d, P, Qi, Qc
    return d, P

def invariant_factors(d, ncols):
    primes = {}
    for x in d:
        if x in (0, 1): continue
        for p, e in sp.factorint(x).items(): primes.setdefault(p, []).append(e)
    k = max((len(v) for v in primes.values()), default=0)
    inv = [1]*k
    for p, es in primes.items():
        for i, e in enumerate(sorted(es, reverse=True)): inv[k-1-i] *= p**e
    inv = sorted(x for x in inv if x > 1)
    free = sum(1 for x in d if x == 0) + (ncols - len(d))
    return inv, free

def group_str(inv, free):
    from collections import Counter
    parts = ([f"Z^{free}"] if free > 1 else ["Z"] if free == 1 else [])
    for v, mult in sorted(Counter(inv).items()):
        parts.append(f"(Z/{v})^{mult}" if mult > 1 else f"Z/{v}")
    return " + ".join(parts) if parts else "0"

def torsion_order(d):
    o = 1
    for x in d:
        if x not in (0, 1): o *= x
    return o

# ---------------------------------------------------------------- p-adic Smith forms (no coefficient explosion)
def padic_valuations(M, p, k):
    """Smith form of M over Z/p^k with valuation-minimal pivots. Returns the list of pivot valuations
    (k stands for an entry that is 0 mod p^k). Exact for the p-primary part of coker M whenever no
    nonzero invariant factor is divisible by p^k."""
    m = p**k; A = [[x % m for x in row] for row in M]; N = len(A); nc = len(A[0]) if N else 0
    def val(x):
        v = 0
        while x % p == 0: x //= p; v += 1
        return v
    vals = []; t = 0
    while t < min(N, nc):
        best = None
        for i in range(t, N):
            row = A[i]
            for j in range(t, nc):
                x = row[j]
                if x:
                    vv = val(x)
                    if best is None or vv < best[0]: best = (vv, i, j)
                    if vv == 0: break
            if best is not None and best[0] == 0: break
        if best is None:
            vals.extend([k]*(min(N, nc) - t)); break
        vv, i, j = best
        A[t], A[i] = A[i], A[t]
        if j != t:
            for row in A: row[t], row[j] = row[j], row[t]
        pv = p**vv; u = A[t][t] // pv; uinv = pow(u, -1, m)
        A[t] = [(x*uinv) % m for x in A[t]]
        rt = A[t]
        for i in range(t+1, N):
            f = A[i][t] // pv
            if f:
                ri = A[i]; A[i] = [(a - f*b) % m for a, b in zip(ri, rt)]
        for j in range(t+1, nc): rt[j] = 0
        vals.append(vv); t += 1
    return vals

BIGP = (2147483647, 2147483629)
def exact_rank(M):
    return max(sum(1 for v in padic_valuations(M, p, 1) if v == 0) for p in BIGP)

def coker_structure(M, primes=None, kmax=64):
    """(invariant factors > 1, free rank, complete?) of coker_Z M.  For nonsingular M the primes come from the
    Bareiss determinant (small factors by trial division; a remaining cofactor is reported); otherwise `primes`
    must be supplied and completeness is left to the caller."""
    N = len(M); nc = len(M[0]); r = exact_rank(M); free = N - r; zeros_expected = min(N, nc) - r
    remainder = 1
    if primes is None:
        assert r == N == nc, "singular matrix: supply candidate primes"
        D = abs(bareiss_det(M))
        f = sp.factorint(D, limit=10**6)
        primes = sorted(f.keys()); remainder = 1
        for p_, e in f.items():
            if p_ > 10**6 and not sp.isprime(p_): remainder *= p_**e
    inv = {}
    for p in primes:
        if remainder > 1 and remainder % p == 0: continue
        k = 12
        while True:
            vals = padic_valuations(M, p, k)
            zeros = sum(1 for v in vals if v == k)
            if zeros == zeros_expected or k >= kmax: break
            k *= 2
        for v in vals:
            if 0 < v < k: inv.setdefault(p, []).append(v)
    # assemble invariant factors from the p-primary decomposition
    depth = max((len(v) for v in inv.values()), default=0)
    factors = [1]*depth
    for p, es in inv.items():
        for i, e in enumerate(sorted(es, reverse=True)): factors[depth-1-i] *= p**e
    factors = sorted(x for x in factors if x > 1)
    return factors, free, remainder

def mat(M): return np.array(M, dtype=object)
def eye(n): return [[int(i == j) for j in range(n)] for i in range(n)]
def zeros(m, n): return [[0]*n for _ in range(m)]
def T(M): return [list(r) for r in zip(*M)]
def mul(A, B):
    return (mat(A) @ mat(B)).tolist()
def block(rows):
    return [sum((list(b[i]) for b in row), []) for row in rows for i in range(len(row[0]))]

def bareiss_det(M):
    A = [list(map(int, r)) for r in M]; N = len(A); sign = 1; prev = 1
    for k in range(N-1):
        if A[k][k] == 0:
            sw = next((i for i in range(k+1, N) if A[i][k] != 0), None)
            if sw is None: return 0
            A[k], A[sw] = A[sw], A[k]; sign = -sign
        for i in range(k+1, N):
            for j in range(k+1, N):
                A[i][j] = (A[i][j]*A[k][k] - A[i][k]*A[k][j]) // prev
        prev = A[k][k]
    return sign*A[N-1][N-1]

# ============================================================ complexes
class Complex:
    """Pure 2-dimensional Delta/simplicial complex with a Z/3 vertex typing.
    verts: list; typ: list of types; edges: list of (a,b) with typ[b] = typ[a]+1; tris: list of (a,b,c) types t,t+1,t+2."""
    def __init__(self, name, q, typ, edges, tris):
        self.name, self.q, self.typ = name, q, typ
        self.n = len(typ); self.edges = edges; self.tris = tris
        self.eidx = {e: k for k, e in enumerate(edges)}
        E, F, n = len(edges), len(tris), self.n
        self.d1 = zeros(n, E); self.d2 = zeros(E, F)
        for k, (a, b) in enumerate(edges): self.d1[b][k] += 1; self.d1[a][k] -= 1
        for k, (a, b, c) in enumerate(tris):      # oriented boundary: (b,c) - (a,c) + (a,b), with our typed orientations
            self.d2[self.eidx[(a, b)]][k] += 1
            self.d2[self.eidx[(b, c)]][k] += 1
            self.d2[self.eidx[(c, a)]][k] += 1    # edge c->a (type t+2 -> t) is the typed orientation of {a,c}
        A1 = zeros(n, n)
        for (a, b) in edges: A1[a][b] += 1
        self.A1 = A1; self.A2 = T(A1)

def thin_torus(m):
    """A~_2 tiling of R^2 modulo mZ^2 (3 | m), directions d1=(1,0), d2=(0,1), d3=(-1,-1) raise type."""
    assert m % 3 == 0
    verts = [(a, b) for a in range(m) for b in range(m)]; vidx = {v: k for k, v in enumerate(verts)}
    typ = [(a + b) % 3 for (a, b) in verts]
    D = [(1, 0), (0, 1), (-1, -1)]
    add = lambda v, d: ((v[0]+d[0]) % m, (v[1]+d[1]) % m)
    edges = sorted({(vidx[v], vidx[add(v, d)]) for v in verts for d in D})
    tris = set()
    for v in verts:
        for i in range(3):
            for j in range(3):
                if i != j:
                    a, b, c = vidx[v], vidx[add(v, D[i])], vidx[add(add(v, D[i]), D[j])]
                    t = (a, b, c)
                    while typ[t[0]] != 0: t = (t[1], t[2], t[0])     # canonical rotation: first vertex of type 0
                    tris.add(t)
    return Complex(f"T_{m}Z2 (n={m*m})", 1, typ, edges, sorted(tris))

# --- triangle presentation over PG(2,2), Singer cycle lambda(x) = line x+1 in the labelling below
PTS = list(range(7))
LINES = [frozenset(s) for s in [(0,1,3),(1,2,4),(2,3,5),(3,4,6),(4,5,0),(5,6,1),(6,0,2)]]
LAM = (1, 2, 3, 4, 5, 6, 0)
ORBITS = [(0,1,3),(0,2,6),(0,4,5),(1,2,4),(1,5,6),(2,3,5),(3,4,6)]
TRIPLES = {r for o in ORBITS for r in ((o[0],o[1],o[2]), (o[1],o[2],o[0]), (o[2],o[0],o[1]))}

def check_presentation():
    lamline = {x: LINES[LAM[x]] for x in PTS}
    A = all((y, z, x) in TRIPLES for (x, y, z) in TRIPLES)
    B = all(((y in lamline[x]) == any((x, y, z) in TRIPLES for z in PTS)) for x in PTS for y in PTS)
    Cc = all(len([z for z in PTS if (x, y, z) in TRIPLES]) <= 1 for x in PTS for y in PTS)
    nofix = all(x not in lamline[x] for x in PTS)
    return A and B and Cc and nofix and len(TRIPLES) == 21

def thick_cover(moduli):
    """Abelian voltage cover of the 3-vertex type quotient X_3 with group ⊕ Z/m_i (m_i | (2,2,14))."""
    edges3 = [(i, x) for i in range(3) for x in PTS]; eidx = {e: k for k, e in enumerate(edges3)}
    tris3 = [(i, o) for i in range(3) for o in ORBITS]
    E = len(edges3)
    d1 = zeros(3, E)
    for k, (i, x) in enumerate(edges3): d1[(i+1) % 3][k] += 1; d1[i][k] -= 1
    d2 = zeros(E, len(tris3))
    for k, (i, (x, y, z)) in enumerate(tris3):
        d2[eidx[(i, x)]][k] += 1; d2[eidx[((i+1) % 3, y)]][k] += 1; d2[eidx[((i+2) % 3, z)]][k] += 1
    # H_1(X_3) = ker d1 / im d2 with explicit coordinates
    D1 = sp.Matrix(d1); ns = D1.nullspace()
    K = sp.Matrix.hstack(*[v*sp.ilcm(*[sp.fraction(x)[1] for x in v]) for v in ns])
    d, P = diagonalize([[int(K[i, j]) for j in range(K.cols)] for i in range(K.rows)], True)
    Pinv = sp.Matrix(P).inv(); r = K.cols; Ksat = Pinv[:, :r]
    KsatPinv = (Ksat.T*Ksat).inv()*Ksat.T
    Y = KsatPinv*sp.Matrix(d2)
    assert all(x.q == 1 for x in Y)
    dd, Q = diagonalize([[int(Y[i, j]) for j in range(Y.cols)] for i in range(Y.rows)], True); Qm = sp.Matrix(Q)
    tor_idx = [i for i, x in enumerate(dd) if x != 1]; orders = [dd[i] for i in tor_idx]
    assert orders == [2, 2, 14] and all(x != 0 for x in dd)
    def h1class(c):
        w = Qm*(KsatPinv*sp.Matrix(c))
        return tuple(int(w[i]) % orders[k] for k, i in enumerate(tor_idx))
    tree = {(0, 0), (1, 0)}
    path = {0: [0]*E, 1: [0]*E, 2: [0]*E}
    path[1][eidx[(0, 0)]] = 1; path[2][eidx[(0, 0)]] = 1; path[2][eidx[(1, 0)]] = 1
    proj = lambda c: tuple(ci % mi for ci, mi in zip(c, moduli) if mi > 1)
    mods = [mi for mi in moduli if mi > 1]
    volt = {}
    for e in edges3:
        i, x = e
        if e in tree: volt[e] = proj((0, 0, 0)); continue
        c = [0]*E; c[eidx[e]] = 1
        volt[e] = proj(h1class([path[i][k] + c[k] - path[(i+1) % 3][k] for k in range(E)]))
    A = sorted({proj(c) for c in itertools.product(range(2), range(2), range(14))})
    add = lambda g, h: tuple((a+b) % m for a, b, m in zip(g, h, mods))
    verts = [(i, g) for i in range(3) for g in A]; vidx = {v: k for k, v in enumerate(verts)}
    typ = [i for (i, g) in verts]
    edges = sorted({(vidx[(i, g)], vidx[((i+1) % 3, add(g, volt[(i, x)]))]) for (i, x) in edges3 for g in A})
    tris = set()
    for (i, (x, y, z)) in tris3:
        for g in A:
            g1 = add(g, volt[(i, x)]); g2 = add(g1, volt[((i+1) % 3, y)]); g3 = add(g2, volt[((i+2) % 3, z)])
            assert g3 == g
            tris.add((vidx[(i, g)], vidx[((i+1) % 3, g1)], vidx[((i+2) % 3, g2)]))
    name = "Y_" + str(len(verts)) + " (" + "x".join(f"Z/{m}" for m in mods) + ")"
    return Complex(name, 2, typ, edges, sorted(tris))

def heawood_like(vs, es):
    if len(vs) != 14 or len(es) != 21: return False
    adj = {v: set() for v in vs}
    for a, b in es:
        if a == b or b in adj[a]: return False
        adj[a].add(b); adj[b].add(a)
    if any(len(adj[v]) != 3 for v in adj): return False
    girth = 99
    for s in adj:
        dist = {s: 0}; par = {s: None}; Qu = [s]
        for u in Qu:
            for w in adj[u]:
                if w not in dist: dist[w] = dist[u]+1; par[w] = u; Qu.append(w)
                elif par[u] != w: girth = min(girth, dist[u]+dist[w]+1)
        if len(dist) != 14: return False
    return girth == 6

def local_checks(X):
    n, q = X.n, X.q
    und = {(min(a, b), max(a, b)) for (a, b) in X.edges}
    simplicial = len(und) == len(X.edges) and all(a != b for (a, b) in X.edges) and \
                 len({tuple(sorted(t)) for t in X.tris}) == len(X.tris)
    nbrs = {v: set() for v in range(n)}
    for a, b in und: nbrs[a].add(b); nbrs[b].add(a)
    linkedges = {v: set() for v in range(n)}
    for t in X.tris:
        for v in t:
            o = [u for u in t if u != v]; linkedges[v].add((min(o), max(o)))
    if q == 2:
        links = all(heawood_like(sorted(nbrs[v]), sorted(linkedges[v])) for v in range(n))
    else:   # thin: link is a hexagon
        links = all(len(nbrs[v]) == 6 and len(linkedges[v]) == 6 for v in range(n))
    out_reg = all(sum(X.A1[v]) == q*q+q+1 for v in range(n)) and all(sum(X.A2[v]) == q*q+q+1 for v in range(n))
    tri_per_edge = {}
    for (a, b, c) in X.tris:
        for e in ((a, b), (b, c), (c, a)): tri_per_edge[e] = tri_per_edge.get(e, 0) + 1
    edge_reg = all(tri_per_edge.get(e, 0) == q+1 for e in X.edges)
    A1 = mat(X.A1); trA3 = int(np.trace(A1 @ A1 @ A1)); flag = (trA3 == 3*len(X.tris))
    return simplicial, links, out_reg, edge_reg, flag

# ============================================================ Hodge theory
SMALL_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]

def torsion_general(M):
    """coker structure for a possibly singular integer matrix: p-adic parts for the small primes;
    the caller validates completeness against a pseudo-determinant."""
    inv, free, _ = coker_structure(M, primes=SMALL_PRIMES)
    return inv, free

def hodge(X, full=True):
    d1, d2 = X.d1, X.d2; n, E, F = X.n, len(X.edges), len(X.tris)
    d1d2 = mul(d1, d2); assert all(x == 0 for r in d1d2 for x in r), "d1 d2 != 0"
    D0 = mul(d1, T(d1)); out = {"D0": D0}
    red0 = [[D0[i][j] for j in range(1, n)] for i in range(1, n)]
    inv0, f0, rem0 = coker_structure(red0); out["coker0"] = (inv0, 1); out["rem0"] = rem0
    r1 = exact_rank(d1); r2 = exact_rank(d2)
    out["betti"] = (n - r1, E - r1 - r2, F - r2)
    if full:
        D1 = [[a + b for a, b in zip(ra, rb)] for ra, rb in zip(mul(T(d1), d1), mul(d2, T(d2)))]
        D2 = mul(T(d2), d2)
        out["D1"] = D1; b1 = E - r1 - r2
        if b1 == 0:
            inv1, f1, rem1 = coker_structure(D1); out["coker1"] = (inv1, 0); out["rem1"] = rem1
        else:
            out["coker1"] = torsion_general(D1); out["rem1"] = 1
        inv2, f2 = torsion_general(D2); out["coker2"] = (inv2, f2)
        # completeness: |Tor coker M| = pdet(M) / det Gram(ker_Z M) for symmetric M with saturated kernel lattice
        ev1 = np.linalg.eigvalsh(np.array(D1, dtype=float)); ev2 = np.linalg.eigvalsh(np.array(D2, dtype=float))
        out["logpdet"] = (float(np.sum(np.log(ev1[ev1 > 1e-8]))), float(np.sum(np.log(ev2[ev2 > 1e-8]))))
        def gram_log(K):          # K: list of basis vectors (rows)
            if not K: return 0.0
            G = mul(K, T(K)); return math.log(abs(bareiss_det(G)))
        # ker Delta_2 = ker d2 ; ker Delta_1 = ker d1 ∩ ker d2^t  (saturated bases from column transforms of 0/±1 matrices)
        dd2, _, _, Qc2 = diagonalize(d2, False, True); rk2 = sum(1 for x in dd2 if x != 0)
        K2 = [[Qc2[t][c] for t in range(F)] for c in range(rk2, F)]               # basis vectors of ker_Z d2 (as rows)
        assert all(x == 0 for r in mul(d2, T(K2)) for x in r)
        out["gram2"] = gram_log(K2)
        if b1 > 0:
            S = d1 + T(d2); ddS, _, _, QcS = diagonalize(S, False, True); rkS = sum(1 for x in ddS if x != 0)
            K1 = [[QcS[t][c] for t in range(E)] for c in range(rkS, E)]
            assert all(x == 0 for r in mul(S, T(K1)) for x in r)
            out["gram1"] = gram_log(K1)
        else:
            out["gram1"] = 0.0
        # DKM critical group K(Y) = ker d1 / im(d2 d2^t): P d1 Q = D, the last E-rank columns of Q span ker_Z d1,
        # and coordinates of v in that basis are the last E-rank rows of Qinv v.
        dd, _, Qi, _ = diagonalize(d1, False, True)
        rank = sum(1 for x in dd if x != 0)
        Kcoord = [Qi[i] for i in range(rank, E)]                       # (E-rank) x E
        L = mul(Kcoord, mul(d2, T(d2)))                                  # coordinates of the columns of d2 d2^t
        top = mul([Qi[i] for i in range(rank)], mul(d2, T(d2)))
        assert all(x == 0 for r in top for x in r), "im(d2 d2^t) not in ker d1"
        if len(L) == len(L[0]) and exact_rank(L) == len(L):
            invL, fL, remL = coker_structure(L); out["DKM"] = (invL, 0); out["remDKM"] = remL
        else:
            out["DKM"] = torsion_general(L); out["remDKM"] = 1
    return out

# ============================================================ Hecke-Laplacian and companion
def hecke(X):
    n, q = X.n, X.q; A1, A2 = X.A1, X.A2
    DH = [[(q**3 - 1)*(i == j) + A1[i][j] - q*A2[i][j] for j in range(n)] for i in range(n)]
    ker_right = all(sum(r) == 0 for r in DH); ker_left = all(sum(c) == 0 for c in T(DH))
    if q > 1:
        redH = [[DH[i][j] for j in range(1, n)] for i in range(1, n)]
        invH, _, remH = coker_structure(redH); fH = 1
        kappaH = abs(bareiss_det(redH))
    else:
        invH, fH = torsion_general(DH); kappaH = 0
    # equal cofactors: reduced determinants at three different sinks
    def reduced(M, s): return [[M[i][j] for j in range(n) if j != s] for i in range(n) if i != s]
    cof = {abs(bareiss_det(reduced(DH, s))) for s in (0, n//2, n-1)}
    commute = mul(A1, A2) == mul(A2, A1)
    # joint spectrum (numeric): A1 normal (A1 A2 = A2 A1 with A2 = A1^t)
    ev = np.linalg.eigvals(np.array(A1, dtype=float))
    triv = np.argmin(np.abs(ev - (q*q+q+1)))
    prodL = 1.0
    for k, a in enumerate(ev):
        if k == triv: continue
        b = np.conj(a)
        prodL *= abs(1 - a + q*b - q**3)
    return DH, (invH, fH), kappaH, ker_right and ker_left, cof, commute, prodL/n

def companion(X):
    n, q = X.n, X.q; A1, A2 = X.A1, X.A2; I = eye(n); Z = zeros(n, n)
    C = block([[A1, [[-q*x for x in r] for r in A2], [[q**3*x for x in r] for r in I]], [I, Z, Z], [Z, I, Z]])
    return C

def traces(X, kmax=30, kdirect=6):
    n, q = X.n, X.q; A1, A2 = mat(X.A1), mat(X.A2); I = mat(eye(n))
    P = {}
    P[1] = A1; P[2] = A1 @ P[1] - 2*q*A2; P[3] = A1 @ P[2] - q*(A2 @ P[1]) + 3*q**3*I
    for k in range(4, kmax+1): P[k] = A1 @ P[k-1] - q*(A2 @ P[k-2]) + q**3*P[k-3]
    tr = {k: int(np.trace(P[k])) for k in range(1, kmax+1)}
    C = np.array(companion(X), dtype=object); Ck = C.copy(); direct = {}
    for k in range(1, kdirect+1):
        direct[k] = int(np.trace(Ck)); Ck = Ck @ C
    net = {m: sum(sp.mobius(m//d)*tr[d] for d in sp.divisors(m)) for m in range(1, kmax+1)}
    return tr, direct, net

def spectrum(X):
    """eigenvalues of C^(3); the nine trivial ones are {1,q,q^2} x {1,w,w^2} (grading)."""
    q = X.q; C = np.array(companion(X), dtype=float); ev = np.linalg.eigvals(C)
    w = np.exp(2j*np.pi/3); triv = []
    for t in (1, q, q*q):
        for k in range(3):
            i = min((i for i in range(len(ev)) if i not in triv), key=lambda i: abs(ev[i] - t*w**k)); triv.append(i)
    nontriv = np.array([ev[i] for i in range(len(ev)) if i not in triv])
    return ev, nontriv

def reduced_net_traces(tr, kmax3):
    """Power sums of the reduced spectrum Lambda_Q (one representative per w-orbit, cubed): p_d = Tr C^{3d}/3."""
    p = {d: tr[3*d] // 3 for d in range(1, kmax3+1)}
    ok = all(tr[3*d] % 3 == 0 for d in range(1, kmax3+1))
    net = {m: sum(sp.mobius(m//d)*p[d] for d in sp.divisors(m)) for m in range(1, kmax3+1)}
    return p, net, ok

def trace_bound_reduced(X, M, kmax):
    """smallest m0 with the lower bound for the reduced net traces positive for all m0 <= m <= kmax."""
    n, q = X.n, X.q; N = n - 1; Q6, Q3, M3 = q**6, q**3, M**3
    def lower(m):
        main = Q6**m + Q3**m + 1 - N*M3**m
        return main - sum(Q6**d + Q3**d + 1 + N*M3**d for d in range(1, m//2 + 1))
    return next((m for m in range(1, kmax+1) if all(lower(mm) > 0 for mm in range(m, kmax+1))), None)

# ============================================================ signed shadow
def signed_shadow(M):
    """Non-negative B with U (I - B^t) V = M + (-I_aux) + (-I_N); requires M_ee >= 1."""
    N = len(M)
    pos = [(e, f) for e in range(N) for f in range(N) if e != f and M[e][f] > 0]
    P_ = len(pos); size = N + P_ + N
    B = zeros(size, size)
    for e in range(N):
        for f in range(N):
            if e != f and M[e][f] < 0: B[f][e] += -M[e][f]           # arcs f -> e
    for k, (e, f) in enumerate(pos):
        w = N + k
        B[f][w] += M[e][f]; B[w][e] += 1; B[w][w] += 2               # f -> w (M_ef times), w -> e, two loops
    for e in range(N):
        v = N + P_ + e
        assert M[e][e] >= 1
        B[e][v] += 1; B[v][e] += M[e][e] - 1; B[v][v] += 2
    IB = np.eye(size, dtype=np.int64) - np.array(B, dtype=np.int64).T
    R = IB[:N, :N]; X = IB[:N, N:]; Y = IB[N:, :N]; Daux = IB[N:, N:]
    # U(I-B^t)V = (R + XY) + (-I) with U = [[I,X],[0,I]], V = [[I,0],[Y,I]] holds iff Daux = -I; check R+XY = M
    ident = np.array_equal(Daux, -np.eye(size - N, dtype=np.int64)) and np.array_equal(R + X @ Y, np.array(M, dtype=np.int64))
    nonneg = all(x >= 0 for r in B for x in r)
    # strong connectivity
    adj = [[j for j in range(size) if B[i][j] > 0] for i in range(size)]
    def reach(s):
        seen = {s}; st = [s]
        while st:
            u = st.pop()
            for w in adj[u]:
                if w not in seen: seen.add(w); st.append(w)
        return seen
    radj = [[i for i in range(size) if B[i][j] > 0] for j in range(size)]
    def rreach(s):
        seen = {s}; st = [s]
        while st:
            u = st.pop()
            for w in radj[u]:
                if w not in seen: seen.add(w); st.append(w)
        return seen
    strong = len(reach(0)) == size and len(rreach(0)) == size
    double_loops = sum(1 for i in range(size) if B[i][i] == 2)
    return B, ident, nonneg, strong, double_loops, P_

# ============================================================ main
def report(X, full=True, shadow1=True):
    n, q, E, F = X.n, X.q, len(X.edges), len(X.tris)
    print(f"===== {X.name}: q={q}, n={n}, E={E}, F={F}, chi={n-E+F} =====")
    simp, links, out_reg, edge_reg, flag = local_checks(X)
    print(f"  (P) simplicial={simp}; links={'Heawood' if q==2 else 'hexagon'}:{links}; {q*q+q+1} type+1 and type-1 neighbours:{out_reg}; "
          f"{q+1} triangles per edge:{edge_reg}; flag (Tr A1^3 = 3F):{flag}")
    H = hodge(X, full)
    print(f"  (D) d1 d2 = 0; Betti (b0,b1,b2) = {H['betti']}; predicted b2 = chi-1+b1 = {n-E+F-1+H['betti'][1]}")
    print(f"  (H) coker Delta_0 = {group_str(*H['coker0'])}" + (f"  [unfactored cofactor {H['rem0']}]" if H['rem0'] > 1 else ""))
    if full:
        inv1, f1 = H['coker1']; inv2, f2 = H['coker2']
        t1 = math.prod(inv1); t2 = math.prod(inv2)
        l1 = math.log(t1) if t1 > 1 else 0.0; l2 = math.log(t2) if t2 > 1 else 0.0
        print(f"      coker Delta_1 = {group_str(inv1, f1)}" + (f"  [unfactored cofactor {H['rem1']}]" if H['rem1'] > 1 else "")
              + f";  |Tor| = {t1};  log|Tor| + log Gram(ker) = {l1 + H['gram1']:.6f} vs log pdet = {H['logpdet'][0]:.6f}: complete {abs(l1 + H['gram1'] - H['logpdet'][0]) < 1e-6}")
        print(f"      coker Delta_2 = {group_str(inv2, f2)};  |Tor| = {t2};  log|Tor| + log Gram(ker) = {l2 + H['gram2']:.6f} vs log pdet = {H['logpdet'][1]:.6f}: complete {abs(l2 + H['gram2'] - H['logpdet'][1]) < 1e-6}")
        print(f"      DKM K(Y) = Z_1/im(d2 d2^t) = {group_str(*H['DKM'])}" + (f"  [unfactored cofactor {H['remDKM']}]" if H['remDKM'] > 1 else ""))
    DH, cokH, kappaH, twosided, cof, commute, prodL = hecke(X)
    print(f"  (L) Delta^H: two-sided kernel Z1:{twosided}; A1A2=A2A1:{commute}; coker Delta^H = {group_str(*cokH)}; "
          f"|Jac^H|={kappaH}; equal cofactors:{len(cof)==1 and cof=={kappaH}}; (1/n)prod|P_j(1)| (numeric)={prodL:.6f}")
    # (C) unimodular reduction of I - C
    C = companion(X); I = eye(n); Z = zeros(n, n)
    IC = [[(i == j) - C[i][j] for j in range(len(C))] for i in range(len(C))]
    U = block([[I, [[-(q*X.A2[i][j] - q**3*(i == j)) for j in range(n)] for i in range(n)], [[q**3*(i == j) for j in range(n)] for i in range(n)]], [Z, I, Z], [Z, Z, I]])
    V = block([[I, Z, Z], [I, I, Z], [I, I, I]])
    red = mul(mul(U, IC), V)
    Pone = [[(i == j) - X.A1[i][j] + q*X.A2[i][j] - q**3*(i == j) for j in range(n)] for i in range(n)]
    target = block([[Pone, Z, Z], [Z, I, Z], [Z, Z, I]])
    negDH = [[-x for x in r] for r in DH]
    print(f"  (C) U(I-C)V = P(1) + I + I: {red == target};  P(1) = -Delta^H: {Pone == negDH}")
    # (T) traces
    tr, direct, net = traces(X)
    print(f"  (T) Tr C^k (direct, k<=6) = {[direct[k] for k in range(1,7)]}; via Newton recursion = {[tr[k] for k in range(1,7)]}; agree:{all(direct[k]==tr[k] for k in direct)}")
    trA3 = int(np.trace(mat(X.A1) @ mat(X.A1) @ mat(X.A1)))
    print(f"      Tr C^3 = {tr[3]} = Tr A1^3 - 3q|E| + 3nq^3 = {trA3 - 3*q*E + 3*n*q**3}; flag formula n(q-1)^2(q+1) = {n*(q-1)**2*(q+1)}: {tr[3]==n*(q-1)**2*(q+1)} "
          f"(flag condition holds: {flag}; Tr A1^3 = {trA3}, 3F = {3*F});  Tr C^k = 0 for 3 not| k (k<=30): {all(tr[k]==0 for k in tr if k%3)}")
    print(f"      Tr C^(3j), j=1..10: {[tr[3*j] for j in range(1,11)]}")
    print(f"      net traces tr_m, m=1..30: {[net[m] for m in range(1,31)]};  all >= 0: {all(v >= 0 for v in net.values())}")
    ev, nontriv = spectrum(X)
    Mmax = float(np.max(np.abs(nontriv))) if len(nontriv) else 0.0; Mmin = float(np.min(np.abs(nontriv))) if len(nontriv) else 0.0
    pQ, netQ, ok3 = reduced_net_traces(tr, 10)
    m0 = trace_bound_reduced(X, Mmax, 10) if (q > 1 and Mmax < q*q) else None
    print(f"  (R) genuinely nontrivial |lambda| in [{Mmin:.6f}, {Mmax:.6f}] (numeric; the nine trivial ones removed); "
          f"Perron gap to q^2={q*q}: {Mmax < q*q - 1e-9}; Ramanujan (all |lambda|=q): {abs(Mmax-q)<1e-6 and abs(Mmin-q)<1e-6}")
    print(f"      reduced spectrum Lambda_Q: p_d = Tr C^(3d)/3 integral:{ok3}; power sums {[pQ[d] for d in range(1,6)]}; "
          f"net traces m=1..10: {[netQ[m] for m in range(1,11)]}; all >= 0: {all(v >= 0 for v in netQ.values())}; provably > 0 for m >= {m0}")
    # (S) shadows
    if q > 1:
        B, ident, nonneg, strong, dl, P_ = signed_shadow(DH)
        snf_msg = ""
        if len(B) <= 400:
            IBt = [[(i == j) - B[j][i] for j in range(len(B))] for i in range(len(B))]
            prs = sorted(set(SMALL_PRIMES) | set(sp.factorint(kappaH, limit=10**7).keys()))
            inv_sh, f_sh, _ = coker_structure(IBt, primes=prs)
            snf_msg = f"; coker(I-B^t) = {group_str(inv_sh, f_sh)} = coker Delta^H: {(inv_sh, f_sh)==cokH}"
        print(f"  (S) shadow of Delta^H: size {len(B)} ({P_} auxiliaries); U(I-B^t)V = M+(-I):{ident}; B>=0:{nonneg}; strongly connected:{strong}; "
              f"double loops:{dl}{snf_msg}")
    else:
        print(f"  (S) q=1: Delta^H = A_1 - A_2 is antisymmetric with zero diagonal; no Hecke shadow (see paper)")
    if full and shadow1:
        B1, ident1, nonneg1, strong1, dl1, P1 = signed_shadow(H["D1"])
        print(f"      shadow of Delta_1: size {len(B1)} ({P1} auxiliaries); U(I-B^t)V = Delta_1+(-I):{ident1}; B>=0:{nonneg1}; strongly connected:{strong1}; double loops:{dl1}")
    print()

if __name__ == "__main__":
    print("triangle presentation axioms (A),(B),(C), no absolute points, |T|=21:", check_presentation(), "\n")
    report(thin_torus(3))
    report(thin_torus(6))
    report(thick_cover((1, 1, 7)))
    report(thick_cover((2, 2, 2)))
    report(thick_cover((1, 1, 14)), full=True, shadow1=False)
    report(thick_cover((2, 2, 7)), full=False)
    report(thick_cover((2, 2, 14)), full=False)
    print("done.")
