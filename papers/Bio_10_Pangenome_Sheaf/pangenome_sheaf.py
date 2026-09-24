# -*- coding: utf-8 -*-
"""
pangenome_sheaf.py -- the pan-genome sheaf sees nothing; the recombination group does
                                                                    (Bio 10, exploratory)

THE OBJECTS.  A pan-genome graph is the pooled de Bruijn multigraph G of strains H_1..H_p
together with the strain walks (the GFA "P-lines"): the decomposition m = sum spec(H_i) is
KNOWN.  G_i is the subgraph traversed by strain i.

 (F)  The pan-genome cellular sheaf: stalk Z^{S(v)} at a vertex (S(v) = strains through v),
      Z^{S(a)} at an arc, restriction = projection.  It is the direct sum of the constant
      sheaves Z_{G_i}, so H^0 = Z^p, H^1 = (+)_i H^1(G_i), and the integer cokernel of the
      sheaf Laplacian is (+)_i (Z (+) K(|G_i|)).  Verified; there is no interaction term.

 (R)  The recombination group  R = Z_1(G) / sum_i Z_1(G_i):  the circulations of the pooled
      graph modulo those that live inside single strains, i.e. the cycles that exist only by
      switching strain.  Mayer-Vietoris for the cover {G_i} of the 1-complex G gives
             R  ~=  H_1(C),   C_q = (+)_{|sigma| = q+1} H_0(G_sigma),
      the first homology of the nerve of intersection components (Cech differential).  So:
        - p = 2:  R is free of rank c(G_1 cap G_2) - 1   (components of the shared subgraph
          minus one; a chain of b bubbles gives b - 1);
        - p >= 3: R can have TORSION.  Six strains whose intersection pattern is the
          6-vertex triangulation of the projective plane give R = Z/2; the octahedron
          (sphere) gives 0 and the 7-vertex torus gives Z^2.  The torsion is a genuinely
          pan-genomic integer invariant: no single strain and no pair sees it.
      Both sides are computed independently: R by a Smith form of the inclusion of the
      strain cycle lattices into Z^A, H_1(C) from the intersection components.

 (D)  Real isolates of phiX174 fetched from NCBI (if the network allows): the nerve, R, and
      the check R = H_1(C).

Run: python -u pangenome_sheaf.py           (seconds; plus the download on first use)
Imports polyploid_phasing.Pool from ../Bio 9 and pgl3_building from ../Paper 10.
"""

import sys, os, time, math, random, collections, itertools, urllib.request, urllib.parse
import numpy as np
import sympy
from sympy.matrices.normalforms import smith_normal_form

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
sys.path.insert(0, _sib("Bio_09_Polyploid_Phasing"))
import pgl3_building as pb
from ds_lattice import load_phix, row, ROWS
from polyploid_phasing import Pool, mutate

random.seed(10)
BASES = "ACGT"


# ------------------------------------------------------------------ integer linear algebra

def snf_invariants(M):
    """invariant factors (> 1) and rank of an integer matrix (list of rows)"""
    if not M or not M[0]:
        return [], 0
    A = sympy.Matrix(M)
    r = A.rank()
    if r == 0:
        return [], 0
    S = smith_normal_form(A, domain=sympy.ZZ)
    d = [abs(int(S[i, i])) for i in range(min(S.shape)) if S[i, i] != 0]
    return [x for x in d if x != 1], r


def group_str(inv, free):
    parts = [f"Z/{d}" for d in inv]
    if free:
        parts.append(f"Z^{free}" if free > 1 else "Z")
    return " + ".join(parts) if parts else "0"


# ------------------------------------------------------------------ cycle lattices

def cycle_basis(nA, tail, head, contents):
    """fundamental cycles (integer vectors in Z^nA) of the sub-multigraph on the given
    contents, one per non-tree content, via an undirected spanning forest"""
    contents = sorted(contents)
    adj = collections.defaultdict(list)          # v -> [(content, other end, sign)]
    for a in contents:
        u, v = tail[a], head[a]
        adj[u].append((a, v, +1))
        adj[v].append((a, u, -1))
    parent = {}
    tree = set()
    for root in adj:
        if root in parent:
            continue
        parent[root] = (None, None, 0)
        stack = [root]
        while stack:
            x = stack.pop()
            for a, y, s in adj[x]:
                if y not in parent:
                    parent[y] = (x, a, s)             # y reached from x along a with sign s
                    tree.add(a)
                    stack.append(y)
    # signed walk from a vertex up to the root: parent[y] = (x, a, s) with s = +1 iff a is
    # the arc x -> y, so the step y -> x traverses a against its direction when s = +1
    def to_root(v):
        vec = collections.Counter()
        while parent[v][0] is not None:
            x, a, s = parent[v]
            vec[a] -= s
            v = x
        return vec
    basis = []
    for a in contents:
        if a in tree:
            continue
        u, v = tail[a], head[a]
        vec = collections.Counter({a: 1})
        pu, pv = to_root(u), to_root(v)
        # cycle: a (u -> v), then the walk v -> root, then the walk root -> u (= -(u -> root))
        for c, s in pv.items():
            vec[c] += s
        for c, s in pu.items():
            vec[c] -= s
        x = [0] * nA
        for c, s in vec.items():
            x[c] = s
        # certificate: a circulation
        bal = collections.Counter()
        for c, s in vec.items():
            bal[head[c]] += s
            bal[tail[c]] -= s
        assert all(b == 0 for b in bal.values()), "fundamental cycle is not a circulation"
        if any(x):
            basis.append(x)
    return basis


def recombination_group(G):
    """R = Z_1(G) / sum_i Z_1(G_i) via a Smith form of the sublattice generators"""
    all_contents = range(G.nA)
    full = cycle_basis(G.nA, G.tail, G.head, all_contents)
    rank_Z1 = len(full)
    gens = []
    p = len(G.seqs)
    for i in range(p):
        ci = {G.cid[a["content"]] for a in G.arcs if a["hap"] == i}
        gens += cycle_basis(G.nA, G.tail, G.head, ci)
    if not gens:
        return [], rank_Z1, rank_Z1, 0
    inv, r = snf_invariants([list(col) for col in zip(*gens)])   # nA x #gens
    return inv, rank_Z1 - r, rank_Z1, r


# ------------------------------------------------------------------ the nerve of components

def strain_subgraphs(G):
    p = len(G.seqs)
    V, A = [], []
    for i in range(p):
        ci = {G.cid[a["content"]] for a in G.arcs if a["hap"] == i}
        vi = {G.tail[a] for a in ci} | {G.head[a] for a in ci}
        V.append(vi)
        A.append(ci)
    return V, A


def components(vertices, arcs, tail, head):
    """connected components of the subgraph (vertices, arcs): {vertex: component id}"""
    par = {v: v for v in vertices}

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for a in arcs:
        x, y = find(tail[a]), find(head[a])
        if x != y:
            par[x] = y
    roots = sorted({find(v) for v in vertices})
    rid = {r: i for i, r in enumerate(roots)}
    return {v: rid[find(v)] for v in vertices}, len(roots)


def nerve_complex(G, verbose=False):
    """the chain complex C_q = (+)_{|sigma|=q+1} H_0(G_sigma); returns H_1 invariants, rank,
    and the component counts per subset"""
    V, A = strain_subgraphs(G)
    p = len(V)
    comp = {}          # sigma -> (vertex->component, count)
    counts = {}
    for q in range(p):
        for sigma in itertools.combinations(range(p), q + 1):
            vs = set.intersection(*(V[i] for i in sigma))
            as_ = set.intersection(*(A[i] for i in sigma))
            if not vs:
                continue
            comp[sigma] = components(vs, as_, G.tail, G.head)
            counts[sigma] = comp[sigma][1]
    # chain groups
    cells = {q: [(s, c) for s, (cm, n) in comp.items() if len(s) == q + 1 for c in range(n)]
             for q in range(p)}
    index = {q: {cell: i for i, cell in enumerate(cells[q])} for q in cells}

    def boundary(q):
        rows, cols = len(cells[q - 1]), len(cells[q])
        M = [[0] * cols for _ in range(rows)]
        for j, (s, c) in enumerate(cells[q]):
            # a representative vertex of component c of G_s
            rep = next(v for v, cc in comp[s][0].items() if cc == c)
            for pos in range(len(s)):
                face = s[:pos] + s[pos + 1:]
                cf = comp[face][0][rep]
                M[index[q - 1][(face, cf)]][j] += (-1) ** pos
        return M
    # H_1 = ker d_1 / im d_2 ; torsion(H_1) = torsion(C_1 / im d_2), rank = dim ker d_1 - rank d_2
    n1 = len(cells[1])
    d1 = boundary(1) if n1 else []
    r1 = sympy.Matrix(d1).rank() if n1 else 0
    if p >= 3 and cells[2]:
        d2 = boundary(2)
        inv, r2 = snf_invariants(d2)
    else:
        inv, r2 = [], 0
    rank_H1 = n1 - r1 - r2
    if verbose:
        print("      intersection components: " + ", ".join(
            f"{''.join(str(i + 1) for i in s)}:{n}" for s, n in sorted(counts.items(), key=lambda kv: (len(kv[0]), kv[0]))))
    return inv, rank_H1, counts


# ------------------------------------------------------------------ the sheaf

def sheaf_check(G):
    """the pan-genome sheaf F = (+)_i Z_{G_i}: coker of its Laplacian vs (+)_i (Z + K(|G_i|))"""
    V, A = strain_subgraphs(G)
    ok = True
    total = []
    for i in range(len(V)):
        vs, as_ = sorted(V[i]), sorted(A[i])
        pos = {v: j for j, v in enumerate(vs)}
        n = len(vs)
        L = [[0] * n for _ in range(n)]
        for a in as_:
            u, v = pos[G.tail[a]], pos[G.head[a]]
            if u == v:
                continue
            L[u][u] += 1
            L[v][v] += 1
            L[u][v] -= 1
            L[v][u] -= 1
        inv, r = snf_invariants(L)
        total.append((inv, n - r))
    # the sheaf Laplacian is block diagonal in the strain decomposition: assemble and compare
    blocks = []
    for i in range(len(V)):
        vs, as_ = sorted(V[i]), sorted(A[i])
        pos = {v: j for j, v in enumerate(vs)}
        n = len(vs)
        L = [[0] * n for _ in range(n)]
        for a in as_:
            u, v = pos[G.tail[a]], pos[G.head[a]]
            if u == v:
                continue
            L[u][u] += 1
            L[v][v] += 1
            L[u][v] -= 1
            L[v][u] -= 1
        blocks.append(L)
    N = sum(len(b) for b in blocks)
    Lf = [[0] * N for _ in range(N)]
    off = 0
    for b in blocks:
        for i, r_ in enumerate(b):
            for j, x in enumerate(r_):
                Lf[off + i][off + j] = x
        off += len(b)
    inv_f, r_f = snf_invariants(Lf)
    inv_sum = sorted(x for inv, _ in total for x in inv)
    free_sum = sum(f for _, f in total)
    ok = (sorted(inv_f) == inv_sum) and (N - r_f == free_sum)
    return ok, inv_f, N - r_f


# ------------------------------------------------------------------ synthetic pools

def rand_seq(L):
    return "".join(random.choice(BASES) for _ in range(L))


def surface_pool(name, faces, nverts, seglen=40):
    """strains = vertices of a triangulated surface; a segment F per face shared by its
    three strains, a segment E per edge shared by its two strains; strain i reads, for each
    incident edge ij with faces ijk, ijl:   F_ijk E_ij F_ijl S  F_ijl E_ij F_ijk S'
    (both junction directions, so G_i cap G_j is connected), with unique spacers S, S'."""
    edges = sorted({tuple(sorted(e)) for f in faces for e in itertools.combinations(f, 2)})
    F = {tuple(sorted(f)): rand_seq(seglen) for f in faces}
    E = {e: rand_seq(seglen) for e in edges}
    strains = []
    for i in range(nverts):
        s = ""
        for e in edges:
            if i not in e:
                continue
            fs = [f for f in F if set(e) <= set(f)]
            assert len(fs) == 2, (e, fs)
            k_, l_ = fs
            s += F[k_] + E[e] + F[l_] + rand_seq(seglen) + F[l_] + E[e] + F[k_] + rand_seq(seglen)
        strains.append(s)
    return name, strains


HEMI_ICOSAHEDRON = [(0, 1, 2), (0, 2, 3), (0, 3, 4), (0, 4, 5), (0, 5, 1),
                    (1, 2, 4), (2, 3, 5), (3, 4, 1), (4, 5, 2), (5, 1, 3)]
OCTAHEDRON = [(0, 1, 2), (0, 2, 3), (0, 3, 4), (0, 4, 1), (5, 1, 2), (5, 2, 3), (5, 3, 4), (5, 4, 1)]
TORUS7 = [(0, 1, 3), (1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 0), (5, 6, 1), (6, 0, 2),
          (0, 3, 2), (1, 4, 3), (2, 5, 4), (3, 6, 5), (4, 0, 6), (5, 1, 0), (6, 2, 1)]


# ------------------------------------------------------------------ real isolates

def fetch_isolates(n=6):
    cache = os.path.join(HERE, "phix174_isolates.fasta")
    if not os.path.exists(cache):
        try:
            term = urllib.parse.quote('phiX174[All Fields] AND 5300:5500[SLEN] AND biomol_genomic[PROP]')
            url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&retmax=40"
                   f"&term={term}")
            xml = urllib.request.urlopen(url, timeout=60).read().decode()
            ids = [x.split("</Id>")[0] for x in xml.split("<Id>")[1:]]
            if not ids:
                return []
            url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta"
                   f"&retmode=text&id={','.join(ids)}")
            txt = urllib.request.urlopen(url, timeout=120).read().decode()
            open(cache, "w").write(txt)
        except Exception as e:
            print(f"    (network unavailable: {e})")
            return []
    txt = open(cache).read()
    recs = []
    for block in txt.split(">")[1:]:
        lines = block.splitlines()
        name = lines[0].split()[0]
        seq = "".join(lines[1:]).upper()
        if set(seq) <= set("ACGT") and 5300 <= len(seq) <= 5500:
            recs.append((name, seq))
    ref = load_phix()
    K = 12
    COMP = str.maketrans("ACGT", "TGCA")

    def kset(s):
        ss = s + s[:K]
        return {ss[i:i + K] for i in range(len(s))}
    refk = kset(ref)
    out, seen = [], set()
    for name, seq in recs:
        rc = seq.translate(COMP)[::-1]
        if len(kset(rc) & refk) > len(kset(seq) & refk):
            seq = rc                                   # deposited on the other strand
        # rotate to the reference origin when the reference's first 12-mer occurs once
        pos = [i for i in range(len(seq)) if (seq + seq[:K])[i:i + K] == ref[:K]]
        if len(pos) == 1:
            seq = seq[pos[0]:] + seq[:pos[0]]
        share = len(kset(seq) & refk) / len(refk)
        if share < 0.5 or seq in seen:
            continue
        seen.add(seq)
        out.append((name, seq, share))
        if len(out) >= n:
            break
    return out


# ------------------------------------------------------------------ main

def describe(name, G, expect_inv=None, expect_rank=None, verbose_nerve=True):
    inv, rk, rZ1, rS = recombination_group(G)
    inv_c, rk_c, counts = nerve_complex(G, verbose=verbose_nerve)
    agree = (sorted(inv) == sorted(inv_c)) and rk == rk_c
    print(f"    {name}: |V| = {G.nv}, |A| = {G.nA}, rank Z_1(G) = {rZ1}, rank sum Z_1(G_i) = {rS};"
          f"  R = {group_str(inv, rk)}   H_1(nerve) = {group_str(inv_c, rk_c)}   {'agree' if agree else 'DISAGREE'}")
    good = agree
    if expect_inv is not None:
        good &= sorted(inv) == sorted(expect_inv)
    if expect_rank is not None:
        good &= rk == expect_rank
    return good, inv, rk


def main():
    t0 = time.time()
    S = load_phix()
    N = len(S)
    k = 12
    print("=== the pan-genome sheaf and the recombination group ===")

    # ---------------------------------------------------------------- (F) the sheaf
    print("\n--- (F) the pan-genome sheaf is a sum of constant sheaves ---")
    okF = True
    pools = [("diploid, 3 SNPs", [S, mutate(S, {600: 1, 1700: 1, 2900: 1})]),
             ("three strains, 2+2+1 SNPs", [S, mutate(S, {600: 1, 1700: 1}), mutate(S, {1700: 1, 2900: 1, 3900: 1})])]
    for name, seqs in pools:
        G = Pool(seqs, k)
        ok, inv, free = sheaf_check(G)
        okF &= ok
        print(f"    {name}: coker(sheaf Laplacian) = {group_str(inv, free)} = (+)_i (Z + K(|G_i|)): {ok}")
    row("F", "the pan-genome sheaf (stalk Z^{strains through v}, projections) is (+)_i Z_{G_i}: its "
             "Laplacian cokernel is the direct sum of the strains' undirected critical groups plus Z^p; "
             "no interaction term", "2 pools", "exact", okF)

    # ---------------------------------------------------------------- (B) two strains: bubbles
    print("\n--- (B) two strains: R is free of rank c(G_1 cap G_2) - 1 ---")
    okB = True
    print("    (the reference has a two-copy repeat a..b, so shared stretches containing copies of a or b")
    print("     merge into one component of G_1 cap G_2: the rank is c(G_1 cap G_2) - 1, not #bubbles - 1)")
    for b in (1, 2, 3, 5):
        sites = [600, 1700, 2900, 3900, 4900][:b]
        G = Pool([S, mutate(S, {s: 1 for s in sites})], k)
        _, _, counts = nerve_complex(G)
        c12 = counts.get((0, 1), 0)
        good, inv, rk = describe(f"{b} SNP(s), c(G_1 cap G_2) = {c12}", G, expect_inv=[], expect_rank=c12 - 1,
                                 verbose_nerve=False)
        okB &= good
    G = Pool([S, mutate(S, {600: 1, 606: 1, 1700: 1})], k)
    _, _, counts = nerve_complex(G)
    c12 = counts.get((0, 1), 0)
    good, inv, rk = describe(f"SNPs 600, 606 (linked), 1700, c(G_1 cap G_2) = {c12}", G, expect_inv=[],
                             expect_rank=c12 - 1, verbose_nerve=False)
    okB &= good
    # a repeat-free reference: rank exactly #bubbles - 1
    R0 = rand_seq(3000)
    for b in (2, 3, 4):
        sites = [400 * (j + 1) for j in range(b)]
        G = Pool([R0, mutate(R0, {s: 1 for s in sites})], k)
        good, inv, rk = describe(f"repeat-free reference, {b} SNPs", G, expect_inv=[], expect_rank=b - 1,
                                 verbose_nerve=False)
        okB &= good
    row("B", "two strains: R = Z_1(G)/(Z_1(G_1) + Z_1(G_2)) is free of rank c(G_1 cap G_2) - 1; on a "
             "repeat-free reference with b bubbles that is b - 1", "8 pools", "exact", okB)

    # ---------------------------------------------------------------- (S) surfaces
    print("\n--- (S) three or more strains: R = H_1(nerve of intersection components), with torsion ---")
    okS = True
    ks = 16          # random 40-mers: at k = 12 accidental shared k-mers between segments are likely
    for name, faces, nv, exp_inv, exp_rk in (("octahedron (sphere)", OCTAHEDRON, 6, [], 0),
                                              ("7-vertex torus", TORUS7, 7, [], 2),
                                              ("hemi-icosahedron (projective plane)", HEMI_ICOSAHEDRON, 6, [2], 0)):
        nm, strains = surface_pool(name, faces, nv)
        G = Pool(strains, ks)
        _, _, counts = nerve_complex(G)
        simplices = {tuple(sorted(f[:j] + f[j + 1:])) for f in faces for j in range(3)} \
            | {tuple(sorted(f)) for f in faces} | {(i,) for i in range(nv)}
        designed = all(counts.get(s, 0) == (1 if s in simplices else 0)
                       for q in range(nv) for s in itertools.combinations(range(nv), q + 1))
        print(f"    nerve equals the designed triangulation (one component per simplex, none elsewhere): {designed}")
        good, inv, rk = describe(f"{name}, {nv} strains, k = {ks}", G, expect_inv=exp_inv, expect_rank=exp_rk)
        okS &= good and designed
    row("S", "six strains whose intersection nerve is the projective plane have R = Z/2; the sphere gives 0 "
             "and the torus Z^2; R equals H_1 of the nerve of intersection components in every case",
        "3 pools", "exact", okS)

    # ---------------------------------------------------------------- (T) random three-strain pools
    print("\n--- (T) random three-strain pools from a shared segment library ---")
    okT = True
    tors = 0
    for trial in range(12):
        lib = [rand_seq(30) for _ in range(6)]
        strains = []
        for i in range(3):
            order = [random.randrange(6) for _ in range(random.randint(4, 7))]
            strains.append("".join(lib[j] + rand_seq(random.choice([0, 0, 12])) for j in order))
        try:
            G = Pool(strains, k)
        except AssertionError:
            continue
        good, inv, rk = describe(f"trial {trial}", G, verbose_nerve=False)
        okT &= good
        tors += bool(inv)
    print(f"    torsion appeared in {tors} of the random pools")
    row("T", "R computed by Smith form agrees with H_1(nerve) on random three-strain pools", "12 pools",
        "exact", okT)

    # ---------------------------------------------------------------- (D) real isolates
    print("\n--- (D) phiX174 isolates from NCBI ---")
    iso = fetch_isolates(6)
    if len(iso) >= 2:
        names = [n for n, _, _ in iso]
        print(f"    {len(iso)} isolates (strand- and origin-normalised to NC_001422.1): "
              + ", ".join(f"{n} ({len(s)} bp, {sh:.3f} of the reference 12-mers)" for n, s, sh in iso))
        seqs = [s for _, s, _ in iso]
        diff = [[sum(1 for x, y in zip(a, b) if x != y) if len(a) == len(b) else -1 for b in seqs] for a in seqs]
        print("    pairwise substitution counts after normalisation (-1: different lengths): " + str(diff))
        G = Pool(seqs, k)
        good, inv, rk = describe(f"{len(iso)} isolates at k = {k}", G)
        row("D", f"{len(iso)} phiX174 isolates: R computed both ways agrees; R = {group_str(inv, rk)}",
            "1 pool", "exact", good)
    else:
        print("    no isolates available (offline); skipped")

    print("\n=== battery summary ===")
    bad = [x for x in ROWS if not x[4]]
    print(f"  TOTAL {sum(1 for x in ROWS if x[4])}/{len(ROWS)} passed"
          + ("" if not bad else f"   FAILURES: {[(x[0], x[2]) for x in bad]}"))
    print(f"  elapsed {time.time()-t0:.0f}s")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
