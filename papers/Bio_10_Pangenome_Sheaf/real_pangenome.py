# -*- coding: utf-8 -*-
r"""
real_pangenome.py -- the recombination group of a REAL bacterial pan-genome.

Bio 10 defined, for a pan-genome graph G pooled from strains H_1..H_p,

        R  =  Z_1(G) / sum_i Z_1(G_i)            (the recombination group)

and proved R = H_1(N), N the nerve of intersection components:
        C_q = (+)_{|sigma| = q+1} H_0(G_sigma),   Cech differential.
Torsion in R was realised only synthetically, by six strains whose nerve is the six-vertex
triangulation of the projective plane.  Open problem 1 of that note asks whether a real
pan-genome has any.  The nerve needs only CONNECTED COMPONENTS of the shared subgraphs, not
the cycle lattices, so it is computable at genome scale; that is what this script does.

PIPELINE (up to eight complete E. coli chromosomes, ~40 Mbp)
  1. load, and orient each genome to the strand of K-12 MG1655 by k-mer overlap
     (a circular genome's k-mer set is rotation-invariant, so only the strand matters);
  2. branch vertices of the POOLED graph: a k-mer with >= 2 distinct successors or >= 2
     distinct predecessors in the pool.  64-bit hashing finds candidates, exact string
     comparison decides, so no collision survives  (the criterion is Bio 9's, not Bio 5's:
     in a pool every shared k-mer is "repeated" and Bio 5's criterion contracts nothing);
  3. unitigs between branch vertices = the contents; two occurrences in different strains
     are the same content iff the strings agree;
  4. for every nonempty sigma the components of G_sigma = {a : sigma subset S(a)}, by
     union-find; this is the nerve N;
  5. H_1(N): the rational rank from sparse ranks of the Cech boundaries, and freeness from a
     unit-pivot reduction of d_2 over Z -- no Smith normal form at any point.

STANDALONE.  Needs Python 3.8+ and numpy, nothing else: no sympy, and no other file of the
series.  The genomes are fetched from NCBI into ./pangenome/ on first use.

    python -u real_pangenome.py --strains 5 -k 31 --primes 2,3 --out real_pangenome_p5.txt

`--cache` stores the reduced boundary matrices, so re-running with further `--primes` skips
the reduction, which is the expensive part.  `--help` lists the rest.
"""
import sys, os, time, heapq, argparse, collections, itertools, pickle, random, urllib.request
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "pangenome")
COMP = str.maketrans("ACGT", "TGCA")
MASK = (1 << 64) - 1
BIG = 2147483647            # 2^31 - 1, the prime the "rational" ranks are computed modulo
ROWS = []

STRAINS = [
    ("NC_000913.3", "K-12 MG1655"),
    ("NC_010473.1", "K-12 DH10B"),
    ("NC_002695.2", "O157:H7 Sakai"),
    ("NC_004431.1", "CFT073"),
    ("NC_011750.1", "IAI39"),
    ("NC_011742.1", "S88"),
    ("NC_008253.1", "536"),
    ("NC_007946.1", "UTI89"),
]

ARGS = None                 # set in main()


def row(key, desc, cases, mode, ok):
    ROWS.append((key, desc, cases, mode, bool(ok)))
    print(("  [ok]   " if ok else "  [FAIL] ") + f"({key}) {desc}   [{cases}] [{mode}]")
    return bool(ok)


def rc(s):
    return s.translate(COMP)[::-1]


def _sib(name):
    """the folder `name`, wherever it sits: beside this script, beside any ancestor of it, or
    one level inside any ancestor.  Only used to notice a genome already downloaded by another
    note of the series; a standalone copy of this script simply fetches its own."""
    here = os.path.abspath(HERE)
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
    return os.path.join(os.path.dirname(HERE), name)


def load(acc):
    """the chromosome `acc`, from ./pangenome/, or from Bio 5 if this script sits in the
    series tree, or from NCBI."""
    for cand in (os.path.join(DATA, f"{acc}.fasta"), os.path.join(_sib("Bio_05_Ecoli_Decomposition"), f"{acc}.fasta")):
        if os.path.exists(cand):
            txt = open(cand, encoding="utf-8", errors="replace").read()
            break
    else:
        url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta"
               f"&retmode=text&id={acc}")
        print(f"    fetching {acc} from NCBI ...", flush=True)
        os.makedirs(DATA, exist_ok=True)
        txt = urllib.request.urlopen(url, timeout=300).read().decode("utf-8", "replace")
        open(os.path.join(DATA, f"{acc}.fasta"), "w", encoding="utf-8").write(txt)
    seq = "".join(l.strip() for l in txt.splitlines() if not l.startswith(">")).upper()
    return "".join(c if c in "ACGT" else "A" for c in seq)


# ------------------------------------------------------------------ the pooled graph

def kmer_hashes(S, k):
    """hash of the k-mer starting at every position of the circular sequence S"""
    N = len(S)
    SS = S + S[:k]
    return np.fromiter((hash(SS[i:i + k]) & MASK for i in range(N)), dtype=np.uint64, count=N)


def multi_base_hashes(h, side):
    """the hash values h that carry at least two distinct `side` bases.  Done by lexsort,
    with no arithmetic on h: an earlier version formed h*4 + base, which OVERFLOWS uint64
    and silently dropped every k-mer whose hash exceeded 2^62."""
    order = np.lexsort((side, h))
    hs, bs = h[order], side[order]
    first = np.ones(len(hs), dtype=bool)
    if len(hs) > 1:
        first[1:] = (hs[1:] != hs[:-1]) | (bs[1:] != bs[:-1])
    hd = hs[first]                       # one entry per distinct (hash, base) pair
    start = np.ones(len(hd), dtype=bool)
    if len(hd) > 1:
        start[1:] = hd[1:] != hd[:-1]
    idx = np.nonzero(start)[0]
    cnt = np.diff(np.append(idx, len(hd)))
    return hd[idx][cnt >= 2]


def branch_vertices(seqs, k, audit=False):
    """k-mers with >= 2 distinct successors or >= 2 distinct predecessors in the pool.
    Candidates by hashing on (k-mer, next base) and (k-mer, previous base); the surviving
    k-mers are then re-derived as exact strings and re-tested, so a hash collision can only
    propose a candidate that the string pass discards, and no true branch vertex is lost."""
    B = {"A": 0, "C": 1, "G": 2, "T": 3}
    hs, nxt, prv = [], [], []
    for S in seqs:
        N = len(S)
        SS = S + S[:k + 1]
        hs.append(kmer_hashes(S, k))
        nxt.append(np.frombuffer(bytes(B[SS[i + k]] for i in range(N)), dtype=np.uint8))
        prv.append(np.frombuffer(bytes(B[S[i - 1]] for i in range(N)), dtype=np.uint8))
    h = np.concatenate(hs)
    cand = np.zeros(len(h), dtype=bool)
    for arrs in (nxt, prv):
        side = np.concatenate(arrs)
        cand |= np.isin(h, multi_base_hashes(h, side))
    del h
    words = set()
    off = 0
    for S in seqs:
        N = len(S)
        SS = S + S[:k]
        for i in np.nonzero(cand[off:off + N])[0].tolist():
            words.add(SS[i:i + k])
        off += N
    succ = collections.defaultdict(set)
    pred = collections.defaultdict(set)
    for S in seqs:
        N = len(S)
        SS = S + S[:k + 1]
        for i in range(N):
            w = SS[i:i + k]
            if w in words:
                succ[w].add(SS[i + k])
                pred[w].add(S[i - 1])
    branch = {w for w in words if len(succ[w]) >= 2 or len(pred[w]) >= 2}
    if audit:                       # brute force, for the short sequences of the audit check
        bs, bp = collections.defaultdict(set), collections.defaultdict(set)
        for S in seqs:
            N = len(S)
            SS = S + S[:k + 1]
            for i in range(N):
                bs[SS[i:i + k]].add(SS[i + k])
                bp[SS[i:i + k]].add(S[i - 1])
        exact = {w for w in bs if len(bs[w]) >= 2 or len(bp[w]) >= 2}
        assert branch == exact, (len(branch), len(exact), len(branch ^ exact))
    return branch


def unitigs(seqs, k, branch):
    """contents of the compacted pooled graph: the maximal segments between consecutive
    branch k-mers along each strain, as strings; returns
      contents  : sorted list of distinct content strings
      strain_of : content index -> frozenset of strain indices
      ends      : content index -> (tail k-mer, head k-mer)
      vstrain   : branch k-mer -> frozenset of strain indices"""
    cont = collections.defaultdict(set)
    vstrain = collections.defaultdict(set)
    for gi, S in enumerate(seqs):
        N = len(S)
        SS = S + S + S[:k]
        pos = [i for i in range(N) if SS[i:i + k] in branch]
        assert pos, f"strain {gi} has no branch vertex at k={k}"
        for t in range(len(pos)):
            p, q = pos[t], pos[(t + 1) % len(pos)]
            if q <= p:
                q += N
            cont[SS[p:q + k]].add(gi)
            vstrain[SS[p:p + k]].add(gi)
    contents = sorted(cont)
    return (contents,
            [frozenset(cont[c]) for c in contents],
            [(c[:k], c[-k:]) for c in contents],
            {w: frozenset(v) for w, v in vstrain.items()})


# ------------------------------------------------------------------ sparse elimination

def components(vert_ids, arcs):
    """union-find over the given arcs, restricted to the given vertex ids"""
    par = {v: v for v in vert_ids}

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for (u, v) in arcs:
        a, b = find(u), find(v)
        if a != b:
            par[a] = b
    comp = {}
    for v in vert_ids:
        comp[v] = find(v)
    roots = sorted(set(comp.values()))
    rid = {r: i for i, r in enumerate(roots)}
    return {v: rid[c] for v, c in comp.items()}, len(roots)


def _index_rows(cols):
    rows = collections.defaultdict(set)
    for j, c in enumerate(cols):
        for r in c:
            rows[r].add(j)
    return rows


def _progress(tag, done, total, t0, every=20000):
    if ARGS is not None and ARGS.quiet:
        return
    if done and done % every == 0:
        el = time.time() - t0
        rate = done / el if el else 0
        print(f"        {tag}: {done:,} / {total:,} columns   [{el:.0f}s, {rate:,.0f}/s]",
              flush=True)


def sparse_rank_mod(cols, p, tag="rank"):
    """rank over F_p of a sparse matrix given as a list of {row: value} dicts, by
    Markowitz-ordered elimination.  Columns are taken shortest-first from a lazy heap; the
    earlier version re-sorted the whole live set at every pivot, which is O(n log n) per
    pivot and was the dominant cost of the whole script at genome scale."""
    cols = [{r: v % p for r, v in c.items() if v % p} for c in cols]
    cols = [c for c in cols if c]
    rows = _index_rows(cols)
    alive = set(range(len(cols)))
    heap = [(len(c), j) for j, c in enumerate(cols)]
    heapq.heapify(heap)
    rank, t0 = 0, time.time()
    while heap:
        n, j = heapq.heappop(heap)
        if j not in alive or len(cols[j]) != n:
            continue                                   # a stale heap entry
        c = cols[j]
        if not c:
            alive.discard(j)
            continue
        r = min(c, key=lambda x: len(rows[x]))
        inv = pow(c[r], p - 2, p)
        rank += 1
        alive.discard(j)
        _progress(tag, rank, len(cols), t0)
        for x in [y for y in rows[r] if y in alive]:
            f = (cols[x][r] * inv) % p
            cx = cols[x]
            for rr, vv in c.items():
                nv = (cx.get(rr, 0) - f * vv) % p
                if nv:
                    if rr not in cx:
                        rows[rr].add(x)
                    cx[rr] = nv
                elif rr in cx:
                    del cx[rr]
                    rows[rr].discard(x)
            heapq.heappush(heap, (len(cx), x))
        for rr in list(c):
            rows[rr].discard(j)
    return rank


def unit_reduce(cols, cap=10 ** 12, tag="unit-pivot"):
    """Integer elimination of a sparse matrix (list of {row: value}) using ONLY +-1 pivots,
    shortest column first.  Unit pivots are unimodular row/column operations, so if every
    column is eliminated this way the matrix has all elementary divisors equal to 1.

    Returns (pivots, residual), residual being the columns that could not be reduced because
    no +-1 entry was left.  residual == [] certifies:  all elementary divisors are 1, hence
    im(d) is a DIRECT SUMMAND, hence coker(d) -- and every submodule of it -- is free.

    A column with no unit entry is set aside and not re-examined until some later pivot
    modifies it, at which point it re-enters the heap; so the columns still alive when the
    heap empties are exactly those that never regained a unit entry."""
    cols = [dict(c) for c in cols]
    rows = _index_rows(cols)
    alive = set(j for j, c in enumerate(cols) if c)
    heap = [(len(cols[j]), j) for j in alive]
    heapq.heapify(heap)
    pivots, t0 = 0, time.time()
    while heap:
        n, j = heapq.heappop(heap)
        if j not in alive or len(cols[j]) != n:
            continue                                   # a stale heap entry
        c = cols[j]
        units = [r for r, v in c.items() if v == 1 or v == -1]
        if not units:
            continue                                   # set aside; re-enters if modified
        r = min(units, key=lambda x: len(rows[x]))
        piv = c[r]
        pivots += 1
        alive.discard(j)
        _progress(tag, pivots, len(cols), t0)
        for x in [y for y in rows[r] if y in alive]:
            cx = cols[x]
            f = cx[r] // piv
            for rr, vv in c.items():
                nv = cx.get(rr, 0) - f * vv
                if nv:
                    if abs(nv) > cap:
                        raise OverflowError("entry blow-up in unit_reduce")
                    if rr not in cx:
                        rows[rr].add(x)
                    cx[rr] = nv
                elif rr in cx:
                    del cx[rr]
                    rows[rr].discard(x)
            heapq.heappush(heap, (len(cx), x))
        for rr in list(c):
            rows[rr].discard(j)
    return pivots, [cols[j] for j in sorted(alive) if cols[j]]


# ------------------------------------------------------------------ the nerve

def nerve(p, vid_strain, arc_ends, arc_strain, primes=(2, 3), verbose=True, cache=None):
    """C_0 -> C_1 -> C_2 of the nerve of components, and dim H_1 over Q and over F_p.

    Torsion is detected without any Smith normal form: H_0(N) is free, so universal
    coefficients give dim_{F_p} H_1(N;F_p) = rank_Q H_1(N) + (p-torsion rank of H_1(N)).
    Any excess of the F_p dimension over the rational rank is therefore p-torsion.  This is
    exactly the cheap torsion detector asked for in Bio 10's open problem 3."""
    comp, counts = {}, {}
    for q in range(min(p, 3)):
        for sigma in itertools.combinations(range(p), q + 1):
            ss = frozenset(sigma)
            vs = [v for v, S in vid_strain.items() if ss <= S]
            if not vs:
                continue
            ar = [arc_ends[a] for a in range(len(arc_ends)) if ss <= arc_strain[a]]
            comp[sigma], counts[sigma] = components(vs, ar)
    cells = {q: [(s, c) for s in sorted(comp) if len(s) == q + 1 for c in range(counts[s])]
             for q in range(min(p, 3))}
    index = {q: {cell: i for i, cell in enumerate(cells[q])} for q in cells}
    # one representative vertex per component, in a single pass over each comp[s].  Looking
    # the representative up per cell instead -- next(v for v, cc in comp[s].items() if cc==c)
    # -- is a linear scan of comp[s] for every one of its components, i.e. quadratic in the
    # size of the intersection, and at six strains that alone outweighed the whole reduction.
    reps = {}
    for s, cc in comp.items():
        r = {}
        for v, c in cc.items():
            if c not in r:
                r[c] = v
        reps[s] = r
    if verbose:
        print("      nerve cells: " + ", ".join(f"q={q}: {len(cells[q]):,}" for q in sorted(cells)))

    def boundary_cols(q):
        out = []
        for (s, c) in cells[q]:
            rep = reps[s][c]
            col = collections.defaultdict(int)
            for pos in range(len(s)):
                face = s[:pos] + s[pos + 1:]
                col[index[q - 1][(face, comp[face][rep])]] += (-1) ** pos
            out.append({r: v for r, v in col.items() if v})
        return out

    n1 = len(cells.get(1, []))
    if n1 == 0:
        return {}, 0, counts, (0, 0, True)
    d1 = boundary_cols(1)
    d2 = boundary_cols(2) if cells.get(2) else []

    # the reduction is the expensive step; cache it so further primes are free
    state = None
    if cache and os.path.exists(cache):
        with open(cache, "rb") as f:
            state = pickle.load(f)
        if state.get("shape") != (n1, len(d1), len(d2)):
            print(f"      cache {os.path.basename(cache)} is for a different instance -- ignored")
            state = None
        elif verbose:
            print(f"      reduction restored from {os.path.basename(cache)}"
                  f" ({state['elapsed']:.0f}s when it was computed)")
    if state is None:
        t = time.time()
        r1Q = sparse_rank_mod(d1, BIG, tag="d_1 rank")
        piv, stuck = unit_reduce(d2)
        state = {"shape": (n1, len(d1), len(d2)), "r1Q": r1Q, "piv": piv, "stuck": stuck,
                 "elapsed": time.time() - t}
        if cache:
            with open(cache, "wb") as f:
                pickle.dump(state, f, protocol=4)
            if verbose:
                print(f"      reduction cached in {os.path.basename(cache)}")
    r1Q, piv, stuck = state["r1Q"], state["piv"], state["stuck"]
    unimodular = not stuck
    r2Q = piv if unimodular else piv + sparse_rank_mod(stuck, BIG, tag="residual rank")
    rkQ = n1 - r1Q - r2Q
    if verbose:
        print(f"      rank d_1 = {r1Q}, rank d_2 = {r2Q:,} of {len(d2):,} columns;"
              f" unit-pivot reduction {'COMPLETE' if unimodular else f'left {len(stuck):,} columns'}"
              f"   [{state['elapsed']:.0f}s]")
    dims = {}
    if not unimodular:
        # the unit-pivot phase is a sequence of unimodular column operations over Z, so it is
        # valid modulo every prime and its pivots stay units there: rank_{F_p}(d_2) is
        # piv + rank_{F_p}(stuck), and only the residual block has to be redone per prime.
        # (Doing sparse_rank_mod on the whole of d_2 once per prime costs as much again as
        # the reduction itself, for no extra information.)
        for q in primes:
            t = time.time()
            dims[q] = n1 - sparse_rank_mod(d1, q, tag=f"d_1 mod {q}") \
                - (piv + sparse_rank_mod(stuck, q, tag=f"residual mod {q}"))
            if verbose:
                print(f"      dim H_1(N; F_{q}) = {dims[q]:,}   ({q}-torsion rank"
                      f" {dims[q] - rkQ})   [{time.time()-t:.0f}s]", flush=True)
    return dims, rkQ, counts, (r1Q, r2Q, unimodular)


# ------------------------------------------------------------------ synthetic pools

def surface_pool(faces, nverts, seglen=40, seed=1):
    """strains = vertices of a triangulated surface; a segment F per face shared by its
    three strains, a segment E per edge shared by its two strains; strain i reads, for each
    incident edge ij with faces ijk, ijl:   F_ijk E_ij F_ijl S  F_ijl E_ij F_ijk S'
    (both junction directions, so G_i cap G_j is connected), with unique spacers S, S'.
    Bio 10's construction, reproduced here so that this script stands alone."""
    rng = random.Random(seed)

    def rand_seq(L):
        return "".join(rng.choice("ACGT") for _ in range(L))
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
    return strains


HEMI_ICOSAHEDRON = [(0, 1, 2), (0, 2, 3), (0, 3, 4), (0, 4, 5), (0, 5, 1),
                    (1, 2, 4), (2, 3, 5), (3, 4, 1), (4, 5, 2), (5, 1, 3)]
OCTAHEDRON = [(0, 1, 2), (0, 2, 3), (0, 3, 4), (0, 4, 1),
              (5, 1, 2), (5, 2, 3), (5, 3, 4), (5, 4, 1)]
TORUS7 = [(0, 1, 3), (1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 0), (5, 6, 1), (6, 0, 2),
          (0, 3, 2), (1, 4, 3), (2, 5, 4), (3, 6, 5), (4, 0, 6), (5, 1, 0), (6, 2, 1)]


def pooled_nerve_input(seqs, k):
    """this script's own pipeline, from sequences to the arguments of nerve()"""
    branch = branch_vertices(seqs, k)
    contents, cstrain, ends, vstrain = unitigs(seqs, k, branch)
    vid = {w: i for i, w in enumerate(sorted(vstrain))}
    return ({vid[w]: S for w, S in vstrain.items()},
            [(vid[a], vid[b]) for (a, b) in ends], cstrain, len(branch), len(contents))


def validate_detector():
    """the F_p detector against the three synthetic pools of Bio 10, whose answers are known
    exactly by Smith normal form: octahedron (S^2) -> 0, seven-vertex torus -> Z^2,
    hemi-icosahedron (RP^2) -> Z/2.  If the detector cannot see the Z/2 it is worthless.
    The pools go through the SAME branch/unitig pipeline as the real chromosomes, so this
    validates the code path that is actually used below."""
    print("\n=== (V) the F_p torsion detector against the synthetic pools of Bio 10 ===")
    ok = True
    for name, faces, nv, exp_rank, exp_tors in (
            ("octahedron (S^2)", OCTAHEDRON, 6, 0, {}),
            ("7-vertex torus", TORUS7, 7, 2, {}),
            ("hemi-icosahedron (RP^2)", HEMI_ICOSAHEDRON, 6, 0, {2: 1})):
        seqs = surface_pool(faces, nv)
        vid_strain, arc_ends, cstrain, _, _ = pooled_nerve_input(seqs, 16)
        dims, rkQ, _, _ = nerve(nv, vid_strain, arc_ends, cstrain, primes=(2, 3), verbose=False)
        tors = {q: d - rkQ for q, d in dims.items() if d > rkQ}
        good = (rkQ == exp_rank) and (tors == exp_tors)
        ok &= good
        print(f"    {name:26s} rank {rkQ} (expected {exp_rank}), torsion {tors or '{}'} "
              f"(expected {exp_tors or '{}'})   {'ok' if good else 'MISMATCH'}")
    return row("V", "the F_p detector reproduces the exact Smith-form answers of Bio 10 on the "
                    "three surface pools, including the Z/2 of the projective plane -- so a "
                    "null result on real data means something", "3 pools", "exact", ok)


def audit_branch_criterion():
    rng = random.Random(10)
    frag = load(STRAINS[0][0])
    ok = True
    for _ in range(6):
        kk = rng.choice([11, 15, 21])
        cuts = [frag[rng.randrange(0, len(frag) - 60000):][:20000] for _ in range(3)]
        try:
            branch_vertices(cuts, kk, audit=True)
        except AssertionError as e:
            ok = False
            print(f"    audit FAILED at k={kk}: {e}")
    return row("A", "the hashed branch-vertex criterion agrees with brute force on random 20 kbp "
                    "fragment pools -- the candidate pass loses no true branch vertex",
               "6 pools", "exact", ok)


# ----------------------------------------------------------------------------- main

class Tee(object):
    def __init__(self, path):
        self.f = open(path, "w", encoding="utf-8")
        self.out = sys.stdout

    def write(self, s):
        self.out.write(s)
        self.f.write(s)

    def flush(self):
        self.out.flush()
        self.f.flush()


def main():
    global ARGS
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[2],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--strains", "-n", type=int, default=int(os.environ.get("NMAX", 8)),
                    help=f"how many of the {len(STRAINS)} chromosomes to pool (default 8)")
    ap.add_argument("-k", "--kmer", default=os.environ.get("KS", "31"),
                    help="word length(s), comma-separated (default 31)")
    ap.add_argument("--primes", default=os.environ.get("PRIMES", "2,3"),
                    help="primes tested when the unit-pivot certificate does not complete "
                         "(default 2,3; 2 is the prime the projective-plane construction "
                         "produces, so it is the one to spend time on first)")
    ap.add_argument("--out", default=None, help="also write the transcript to this file")
    ap.add_argument("--cache", action="store_true",
                    help="store the reduced boundary matrices under ./cache/, so re-running "
                         "with further --primes skips the reduction")
    ap.add_argument("--skip-checks", action="store_true",
                    help="go straight to the chromosomes, skipping checks (A) and (V)")
    ap.add_argument("--quiet", action="store_true", help="no progress lines during elimination")
    ARGS = ap.parse_args()
    KS = [int(x) for x in ARGS.kmer.split(",")]
    PRIMES = tuple(int(x) for x in ARGS.primes.split(",") if x.strip())
    if ARGS.out:
        sys.stdout = Tee(ARGS.out)

    t0 = time.time()
    print("=== the recombination group of a real Escherichia coli pan-genome ===")
    print(f"    strains {ARGS.strains}, k = {ARGS.kmer}, primes {PRIMES}, "
          f"python {sys.version.split()[0]}, numpy {np.__version__}")
    if not ARGS.skip_checks:
        audit_branch_criterion()
        validate_detector()
    strains = STRAINS[:ARGS.strains]
    seqs, names = [], []
    ref = load(strains[0][0])
    K0 = 31
    refk = {ref[i:i + K0] for i in range(0, len(ref) - K0, 7)}
    print()
    for acc, nm in strains:
        S = load(acc)
        fwd = sum(1 for i in range(0, len(S) - K0, 101) if S[i:i + K0] in refk)
        r = rc(S)
        rev = sum(1 for i in range(0, len(r) - K0, 101) if r[i:i + K0] in refk)
        if rev > fwd:
            S = r
        seqs.append(S)
        names.append(nm)
        print(f"    {acc:14s} {nm:16s} {len(S):>9,} bp   strand "
              f"{'reverse-complemented' if rev > fwd else 'as deposited'}"
              f"   (fwd {fwd} / rev {rev} sampled 31-mers shared with the reference)")
    p = len(seqs)

    for k in KS:
        print(f"\n=== k = {k} ===")
        t1 = time.time()
        branch = branch_vertices(seqs, k)
        print(f"    branch vertices of the pooled graph: {len(branch):,}   [{time.time()-t1:.0f}s]")
        t1 = time.time()
        contents, cstrain, ends, vstrain = unitigs(seqs, k, branch)
        print(f"    contents (unitigs): {len(contents):,}   [{time.time()-t1:.0f}s]")
        vid = {w: i for i, w in enumerate(sorted(vstrain))}
        arc_ends = [(vid[a], vid[b]) for (a, b) in ends]
        vid_strain = {vid[w]: S for w, S in vstrain.items()}
        core_v = sum(1 for S in vid_strain.values() if len(S) == p)
        core_a = sum(1 for S in cstrain if len(S) == p)
        print(f"    core (in all {p} strains): {core_v:,} vertices, {core_a:,} contents;"
              f" strain-specific contents: {sum(1 for S in cstrain if len(S) == 1):,}")
        cache = None
        if ARGS.cache:
            os.makedirs(os.path.join(HERE, "cache"), exist_ok=True)
            cache = os.path.join(HERE, "cache", f"nerve_p{p}_k{k}.pkl")
        t1 = time.time()
        dims, rkQ, counts, (r1, r2, unimod) = nerve(p, vid_strain, arc_ends, cstrain,
                                                    primes=PRIMES, cache=cache)
        tors = {q: d - rkQ for q, d in dims.items() if d > rkQ}
        tested = ", ".join(str(q) for q in sorted(dims))
        verdict = ("R is FREE of rank {:,} -- certified, every elementary divisor of d_2 is 1"
                   .format(rkQ) if unimod else
                   (f"rational rank {rkQ:,} (mod {BIG}, an upper bound), "
                    + (f"torsion ranks {tors}" if tors else
                       f"no torsion at p = {tested} -- other primes untested")))
        print(f"    R = H_1(nerve): {verdict}   [{time.time()-t1:.0f}s]")
        pair = {s: c for s, c in counts.items() if len(s) == 2}
        print(f"    intersection components: singletons {[counts[(i,)] for i in range(p)]},"
              f" pairs min/median/max "
              f"{min(pair.values())}/{sorted(pair.values())[len(pair)//2]}/{max(pair.values())}"
              + (f", triples max {max(c for s, c in counts.items() if len(s) == 3)}"
                 if p >= 3 else ""))
        row(f"E{k}", f"k={k}: the recombination group of {p} complete E. coli chromosomes is "
                     + (f"free of rank {rkQ:,}, certified by a unit-pivot reduction of d_2"
                        if unimod else
                        f"of rational rank at most {rkQ:,} with torsion {tors}" if tors else
                        f"of rational rank at most {rkQ:,}, with no p-torsion for p = " + tested),
            f"{p} strains", "exact", True)

    print("\n=== summary ===")
    bad = [x for x in ROWS if not x[4]]
    print(f"  {sum(1 for x in ROWS if x[4])}/{len(ROWS)} passed"
          + ("" if not bad else f"   FAILURES: {[(x[0], x[2]) for x in bad]}"))
    print(f"  elapsed {time.time()-t0:.0f}s")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
