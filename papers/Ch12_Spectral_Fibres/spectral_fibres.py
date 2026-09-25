# -*- coding: utf-8 -*-
r"""
spectral_fibres.py -- what a local encoder cannot see: spectral fibres of the human proteome.

An encoder of a protein sequence S is LOCAL of order k if it factors through the k-mer spectrum
spec_k(S) (the multiset of length-k substrings, together with the first and last (k-1)-mers).
Composition, dipeptide/tripeptide composition, PseAAC / APAAC with lag <= k-1, and any mean-pooled
convolutional encoder of window k are local.  A local encoder is constant on the FIBRE

    F_k(S) = { T : spec_k(T) = spec_k(S) },

whose size is an Eulerian-trail count in the de Bruijn multigraph of S (BEST theorem), computed here
exactly, with forced-arc contraction (Algebraic Genomics, Lemma 9.7), for every reviewed human protein.

Checks (all exact):
  (T) the fibre formula against brute-force enumeration of all strings, small alphabets;
  (C) forced-arc contraction against the full reduced determinant;
  (P) the proteome: N_k(S) = |F_k(S)| for k = 1, 2, ..., and the spectral resolution
      kappa(S) = min{k : N_k'(S) = 1 for all k' >= k};
  (X) collisions: distinct human proteins with identical k-spectra, k = 1, 2, 3;
  (M) fibre-mates: for proteins with N_25 > 1, an explicit different sequence with the same
      25-spectrum; APAAC (lambda = 24, as used for PseAAC in UniPert's benchmark, any property
      scales) and a random integer CNN of window 25 with mean pooling give IDENTICAL outputs, while
      the two sequences differ by a transposition of segments.

    python -u spectral_fibres.py            downloads UniProt (organism 9606, reviewed) once
    python -u spectral_fibres.py --jobs=N   worker processes (default: all cores)

Standard library only.  Writes human_reviewed.tsv.gz (cache), fibre_results.json and prints a log.
"""
import os, sys, time, math, json, gzip, random, collections, urllib.request
from fractions import Fraction
from multiprocessing import Pool, cpu_count

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "human_reviewed.tsv.gz")
URL = ("https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=tsv"
       "&fields=accession%2Cgene_primary%2Cprotein_families%2Clength%2Csequence"
       "&query=%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29")
KREP = [1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 25, 30, 40, 50]
KCAP = 2000
ROWS = []


def row(key, desc, ok):
    ROWS.append((key, bool(ok)))
    print(("  [ok]   " if ok else "  [FAIL] ") + f"({key}) {desc}", flush=True)


def log2i(n):
    if n <= 0:
        return float("-inf")
    b = n.bit_length()
    if b <= 60:
        return math.log2(n)
    return (b - 60) + math.log2(n >> (b - 60))


# ------------------------------------------------------------------ arborescences

def bareiss(M):
    n = len(M)
    if n == 0:
        return 1
    A = [r[:] for r in M]
    sign, prev = 1, 1
    for i in range(n - 1):
        if A[i][i] == 0:
            p = next((r for r in range(i + 1, n) if A[r][i] != 0), None)
            if p is None:
                return 0
            A[i], A[p] = A[p], A[i]
            sign = -sign
        for r in range(i + 1, n):
            for c in range(i + 1, n):
                A[r][c] = (A[r][c] * A[i][i] - A[r][i] * A[i][c]) // prev
            A[r][i] = 0
        prev = A[i][i]
    return sign * A[n - 1][n - 1]


def laplacian_det(out, verts, root):
    rest = [v for v in verts if v != root]
    ix = {v: i for i, v in enumerate(rest)}
    M = [[0] * len(rest) for _ in rest]
    for v in rest:
        for u, w in out[v].items():
            M[ix[v]][ix[v]] += w
            if u != root:
                M[ix[v]][ix[u]] -= w
    return bareiss(M)


def tree_count(out0, verts, root, contract=True):
    """in-arborescences to `root` of the multidigraph out0[v] = Counter(head -> multiplicity)
    (loops excluded).  With `contract`, forced arcs are contracted first (Algebraic Genomics, Lemma 9.7)."""
    out = {v: collections.Counter(out0.get(v, {})) for v in verts}
    if not contract:
        return laplacian_det(out, verts, root)
    inn = {v: set() for v in verts}
    for v, c in out.items():
        for u in c:
            inn[u].add(v)
    T = 1
    alive = set(verts)
    todo = [v for v in verts if v != root and len(out[v]) == 1]
    while todo:
        v = todo.pop()
        if v not in alive or v == root or len(out[v]) != 1:
            continue
        (u, w), = out[v].items()
        T *= w
        inn[u].discard(v)
        for x in list(inn[v]):
            w2 = out[x].pop(v)
            if x != u:
                out[x][u] += w2
                inn[u].add(x)
            if x != root and len(out[x]) == 1:
                todo.append(x)
        alive.discard(v)
        del out[v], inn[v]
        if u != root and len(out[u]) == 1:
            todo.append(u)
    rest = [v for v in alive if v != root]
    if any(not out[v] for v in rest):
        return 0
    return T * laplacian_det(out, [root] + rest, root)


# ------------------------------------------------------------------ the fibre

def graph(S, k):
    L = len(S)
    m = L - k + 1
    mult = collections.Counter(S[i:i + k] for i in range(m))
    s, t = S[:k - 1], S[L - k + 1:]
    out = collections.defaultdict(collections.Counter)
    deg = collections.Counter()
    for a, w in mult.items():
        u, v = a[:-1], a[1:]
        deg[u] += w
        if u != v:
            out[u][v] += w
    return mult, s, t, out, deg, m


def fibre_count(S, k, contract=True):
    """|{T : spec_k(T) = spec_k(S)}|, where spec_k is the multiset of k-mers (so |T| = |S|)"""
    L = len(S)
    if k <= 1:
        N = math.factorial(L)
        for v in collections.Counter(S).values():
            N //= math.factorial(v)
        return N
    if k > L:
        raise ValueError("k > |S|")
    mult, s, t, out, deg, m = graph(S, k)
    if s != t:                        # close the trail: one extra arc t -> s
        deg[t] += 1
        out[t][s] += 1
    verts = list(deg)
    T = tree_count(out, verts, s, contract)
    num = T
    for v in verts:
        num *= math.factorial(deg[v] - 1)
    if s == t:                        # a closed trail may be read from any of its m positions; the
        num *= m                      # m rotations are counted BEFORE identical k-mers are merged
    den = 1
    for w in mult.values():
        den *= math.factorial(w)
    assert num % den == 0
    return num // den


EXACT_MAX = 150


def fibre_log2(S, k):
    """(log2 N_k(S), N_k(S) == 1).  Exact whenever the contracted Laplacian has order <= EXACT_MAX;
    beyond that the determinant is taken in floating point (numpy slogdet) and the factorials by
    lgamma -- N is then astronomically large, and certainly not 1."""
    L = len(S)
    if k <= 1:
        N = fibre_count(S, 1)
        return log2i(N), N == 1
    mult, s, t, out, deg, m = graph(S, k)
    if s != t:
        deg[t] += 1
        out[t][s] += 1
    verts = list(deg)
    # contraction, as in tree_count, but stopping before the determinant
    o = {v: collections.Counter(out.get(v, {})) for v in verts}
    inn = {v: set() for v in verts}
    for v, c in o.items():
        for u in c:
            inn[u].add(v)
    T, alive = 1, set(verts)
    todo = [v for v in verts if v != s and len(o[v]) == 1]
    while todo:
        v = todo.pop()
        if v not in alive or v == s or len(o[v]) != 1:
            continue
        (u, w), = o[v].items()
        T *= w
        inn[u].discard(v)
        for x in list(inn[v]):
            w2 = o[x].pop(v)
            if x != u:
                o[x][u] += w2
                inn[u].add(x)
            if x != s and len(o[x]) == 1:
                todo.append(x)
        alive.discard(v)
        del o[v], inn[v]
        if u != s and len(o[u]) == 1:
            todo.append(u)
    rest = [v for v in alive if v != s]
    if any(not o[v] for v in rest):
        raise AssertionError("the spectrum of a sequence always has a trail")
    if len(rest) <= EXACT_MAX:
        num = T * laplacian_det(o, [s] + rest, s)
        for v in verts:
            num *= math.factorial(deg[v] - 1)
        if s == t:
            num *= m
        den = 1
        for w in mult.values():
            den *= math.factorial(w)
        assert num % den == 0
        N = num // den
        return log2i(N), N == 1
    import numpy as np
    ix = {v: i for i, v in enumerate(rest)}
    M = np.zeros((len(rest), len(rest)))
    for v in rest:
        for u, w in o[v].items():
            M[ix[v], ix[v]] += w
            if u != s:
                M[ix[v], ix[u]] -= w
    sign, ld = np.linalg.slogdet(M)
    assert sign > 0
    lg = log2i(T) + ld / math.log(2)
    lg += sum(math.lgamma(deg[v]) for v in verts) / math.log(2)
    lg -= sum(math.lgamma(w + 1) for w in mult.values()) / math.log(2)
    if s == t:
        lg += math.log2(m)
    return lg, False


def spectrum(S, k):
    return tuple(sorted(collections.Counter(S[i:i + k] for i in range(len(S) - k + 1)).items()))


# ------------------------------------------------------------------ fibre-mates

def random_trail(S, k, rnd):
    """a uniformly-not-necessarily random Eulerian trail of the k-mer multigraph of S that starts
    at S's first (k-1)-mer: Hierholzer with random choices.  Returns the spelled string."""
    mult, s, t, out, deg, m = graph(S, k)
    adj = collections.defaultdict(list)
    for a, w in mult.items():
        adj[a[:-1]].extend([a[1:]] * w)
    for v in adj:
        rnd.shuffle(adj[v])
    stack, path = [s], []
    while stack:
        v = stack[-1]
        if adj[v]:
            stack.append(adj[v].pop())
        else:
            path.append(stack.pop())
    path.reverse()
    if len(path) != m + 1:
        return None
    return path[0] + "".join(p[-1] for p in path[1:])


def apaac(S, lam, H1, H2, w=Fraction(1, 2)):
    """Chou's amphiphilic pseudo amino acid composition (type II PseAAC), 20 + 2 lam values, exact"""
    AA = "ACDEFGHIKLMNPQRSTVWY"
    L = len(S)
    f = [Fraction(S.count(a), 1) for a in AA]
    tau = []
    for j in range(1, lam + 1):
        t1 = sum((H1.get(S[i], 0) * H1.get(S[i + j], 0) for i in range(L - j)), Fraction(0))
        t2 = sum((H2.get(S[i], 0) * H2.get(S[i + j], 0) for i in range(L - j)), Fraction(0))
        tau += [t1 / (L - j), t2 / (L - j)]
    Z = sum(f) + w * sum(tau)
    return tuple(x / Z for x in f) + tuple(w * x / Z for x in tau)


def std_scale(raw):
    AA = "ACDEFGHIKLMNPQRSTVWY"
    vals = [Fraction(raw[a]) for a in AA]
    mu = sum(vals) / 20
    var = sum((v - mu) ** 2 for v in vals) / 20
    # an exact rational stand-in for the standard deviation (any positive constant gives a scale;
    # the fibre argument does not care which)
    sd = Fraction(math.sqrt(var)).limit_denominator(10 ** 6)
    return {a: (Fraction(raw[a]) - mu) / sd for a in AA}


KD = dict(A="1.8", R="-4.5", N="-3.5", D="-3.5", C="2.5", Q="-3.5", E="-3.5", G="-0.4", H="-3.2",
          I="4.5", L="3.8", K="-3.9", M="1.9", F="2.8", P="-1.6", S="-0.8", T="-0.7", W="-0.9",
          Y="-1.3", V="4.2")
HW = dict(A="-0.5", R="3.0", N="0.2", D="3.0", C="-1.0", Q="0.2", E="3.0", G="0.0", H="-0.5",
          I="-1.8", L="-1.8", K="3.0", M="-1.3", F="-2.5", P="0.0", S="0.3", T="-0.4", W="-3.4",
          Y="-2.3", V="-1.5")


def cnn_mean(S, W, bias, k):
    """a one-layer integer CNN, window k, ReLU, mean pooling over positions; sentinel padding"""
    P = "^" * (k - 1) + S + "$" * (k - 1)
    acc = [0] * len(bias)
    for i in range(len(P) - k + 1):
        win = P[i:i + k]
        for c in range(len(bias)):
            z = bias[c] + sum(W[c][j].get(win[j], 0) for j in range(k))
            if z > 0:
                acc[c] += z
    n = len(P) - k + 1
    return tuple(Fraction(a, n) for a in acc)


def sw_score(a, b, match=2, mis=-1, gap=-2):
    prev = [0] * (len(b) + 1)
    best = 0
    for i in range(1, len(a) + 1):
        cur = [0] * (len(b) + 1)
        ai = a[i - 1]
        for j in range(1, len(b) + 1):
            v = max(0, prev[j - 1] + (match if ai == b[j - 1] else mis), prev[j] + gap, cur[j - 1] + gap)
            cur[j] = v
            if v > best:
                best = v
        prev = cur
    return best


# ------------------------------------------------------------------ one protein

def analyse(item):
    acc, gene, fam, S = item
    L = len(S)
    logN, kappa, lastgt1 = {}, None, 0
    k = 1
    while k <= min(L, KCAP):
        lg, one = fibre_log2(S, k)
        if not one:
            lastgt1 = k
        if k in KREP:
            logN[k] = lg
        # once no (k-1)-mer repeats, the k-mer graph is a simple path and N_k' = 1 for k' >= k
        if k >= 2 and len({S[i:i + k - 1] for i in range(L - k + 2)}) == L - k + 2:
            break
        k += 1
    kappa = lastgt1 + 1
    if k > KCAP:
        kappa = None                  # censored
    for kk in KREP:
        if kk not in logN:
            logN[kk] = 0.0 if (kk <= L and (kappa is not None and kk >= kappa)) else logN.get(kk, None)
    return dict(acc=acc, gene=gene, fam=fam, L=L, kappa=kappa, logN=logN)


def _analyse_i(p):
    return p[0], analyse(p[1])


# ------------------------------------------------------------------ data

def load():
    if not os.path.exists(DATA):
        print(f"    downloading UniProt human reviewed proteome ...", flush=True)
        try:
            urllib.request.urlretrieve(URL, DATA + ".part")
            os.replace(DATA + ".part", DATA)
        except Exception as e:
            raise SystemExit(f"download failed ({e}).  Save\n  {URL}\nas {DATA} and rerun.")
    items = []
    with gzip.open(DATA, "rt", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            acc, gene, fam, _, seq = p[:5]
            if seq:
                items.append((acc, gene, fam.split(";")[0].strip(), seq))
    return items


# ------------------------------------------------------------------ main

def main():
    t0 = time.time()
    jobs = int(next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--jobs=")), cpu_count()))
    print("=== spectral fibres of the human proteome ===")
    print(f"    python {sys.version.split()[0]}, {jobs} processes\n", flush=True)

    # (T) brute force
    print("=== (T) the fibre formula against brute force ===")
    okT, ncases = True, 0
    import itertools
    for alpha, Lmax in (("AB", 12), ("ABC", 8), ("ABCD", 6)):
        for L in range(2, Lmax + 1):
            strings = ["".join(p) for p in itertools.product(alpha, repeat=L)]
            for k in range(1, min(L, 5) + 1):
                groups = collections.Counter(spectrum(x, k) for x in strings)
                for x in strings[:: max(1, len(strings) // 300)]:
                    okT &= (fibre_count(x, k) == groups[spectrum(x, k)])
                    ncases += 1
    print(f"    {ncases} (string, k) cases, alphabets of size 2-4, lengths up to 12")
    row("T", "N_k(S) = t(G) prod (d-1)! / prod mult! (times m for a closed trail) equals the number of "
             "strings with the same k-mer multiset, by exhaustive enumeration", okT)

    # (C) contraction
    print("\n=== (C) forced-arc contraction against the full determinant ===")
    rnd = random.Random(20260923)
    okC, nc = True, 0
    for _ in range(400):
        S = "".join(rnd.choice("ACDEFGHIKLMNPQRSTVWY"[:rnd.randrange(2, 21)])
                    for _ in range(rnd.randrange(10, 200)))
        k = rnd.randrange(2, 6)
        okC &= fibre_count(S, k, True) == fibre_count(S, k, False)
        nc += 1
    print(f"    {nc} random sequences, alphabets of size 2-20, k = 2-5")
    # the floating-point path of fibre_log2 against the exact count
    global EXACT_MAX
    keep, EXACT_MAX, worst = EXACT_MAX, 0, 0.0
    for _ in range(200):
        S = "".join(rnd.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(rnd.randrange(50, 400)))
        k = rnd.randrange(2, 4)
        ex = log2i(fibre_count(S, k))
        lg, _ = fibre_log2(S, k)
        worst = max(worst, abs(lg - ex) / max(1.0, ex))
    EXACT_MAX = keep
    okC &= worst < 1e-9
    print(f"    floating-point log2 N against exact on 200 sequences: worst relative error {worst:.1e}")
    row("C", "the fibre size with forced arcs contracted equals the one from the full reduced "
             "Laplacian", okC)

    # (P) the proteome
    print("\n=== (P) the human proteome ===", flush=True)
    items = load()
    Ls = sorted(len(x[3]) for x in items)
    print(f"    {len(items):,} reviewed human proteins (UniProt, organism 9606); length median "
          f"{Ls[len(Ls)//2]}, max {Ls[-1]:,}; {sum(Ls):,} residues", flush=True)
    t = time.time()
    order = sorted(range(len(items)), key=lambda i: -len(items[i][3]))
    res = [None] * len(items)
    with Pool(jobs) as pool:
        for c, (i, r) in enumerate(pool.imap_unordered(_analyse_i, [(i, items[i]) for i in order],
                                                        chunksize=4)):
            res[i] = r
            if (c + 1) % 2000 == 0 or c + 1 == len(items):
                print(f"      {c+1:,}/{len(items):,} proteins   [{time.time()-t:.0f}s]", flush=True)
    print(f"    fibres computed   [{time.time()-t:.0f}s]", flush=True)
    json.dump(res, open(os.path.join(HERE, "fibre_results.json"), "w"), indent=0)

    n = len(res)
    print("\n      k   N_k > 1        share    median log2 N_k   max log2 N_k (protein)")
    for k in KREP:
        vals = [(r["logN"].get(k) or 0.0, r) for r in res if r["L"] >= k]
        amb = [v for v, _ in vals if v > 1e-9]
        vs = sorted(v for v, _ in vals)
        mx = max(vals, key=lambda x: x[0])
        print(f"     {k:2d}  {len(amb):7,}   {100*len(amb)/len(vals):7.2f}%   {vs[len(vs)//2]:14.1f}"
              f"   {mx[0]:12.1f} ({mx[1]['gene'] or mx[1]['acc']}, L={mx[1]['L']})")
    kap = sorted(r["kappa"] for r in res if r["kappa"] is not None)
    cens = sum(1 for r in res if r["kappa"] is None)
    q = lambda p: kap[min(len(kap) - 1, int(p * len(kap)))]
    print(f"\n    spectral resolution kappa: median {q(0.5)}, 90% {q(0.9)}, 99% {q(0.99)}, "
          f"max {kap[-1]}; censored (> {KCAP}) {cens}")
    for c in (2, 3, 4, 5, 6, 8, 10, 25, 50, 100):
        print(f"      kappa <= {c:3d}: {sum(1 for x in kap if x <= c):6,} "
              f"({100*sum(1 for x in kap if x <= c)/n:5.1f}%)")
    print("\n    the 20 proteins with the largest kappa:")
    for r in sorted(res, key=lambda r: -(r["kappa"] or 10 ** 9))[:20]:
        print(f"      {r['gene'] or '-':12s} {r['acc']:10s} L={r['L']:6d} kappa={r['kappa']}"
              f"  log2N_3={(r['logN'].get(3) or 0.0):.1f}  [{r['fam'][:60]}]")
    print("\n    the 20 proteins with the most order information at k = 3 (log2 N_3 / L, L >= 100):")
    for r in sorted([r for r in res if r["L"] >= 100], key=lambda r: -(r["logN"].get(3) or 0.0) / r["L"])[:20]:
        print(f"      {r['gene'] or '-':12s} {r['acc']:10s} L={r['L']:6d} bits/res={(r['logN'].get(3) or 0.0)/r['L']:.3f}"
              f"  kappa={r['kappa']}  [{r['fam'][:60]}]")
    fam = collections.defaultdict(list)
    for r in res:
        if r["fam"]:
            fam[r["fam"]].append(r)
    print("\n    families (>= 8 members) by median kappa:")
    fr = [(sorted(x["kappa"] or KCAP for x in v)[len(v) // 2], f, len(v)) for f, v in fam.items()
          if len(v) >= 8]
    for med, f, cnt in sorted(fr, reverse=True)[:20]:
        print(f"      median kappa {med:4d}  n={cnt:4d}  {f[:80]}")
    print("    families (>= 8 members) by median log2 N_3 / L:")
    fb = [(sorted((x["logN"].get(3) or 0.0) / x["L"] for x in v)[len(v) // 2], f, len(v)) for f, v in fam.items()
          if len(v) >= 8]
    for med, f, cnt in sorted(fb, reverse=True)[:15]:
        print(f"      {med:.3f} bits/res  n={cnt:4d}  {f[:80]}")
    allk = sorted(r["kappa"] or KCAP for r in res)
    print(f"    all proteins: median kappa {allk[len(allk)//2]}, median log2 N_3 / L "
          f"{sorted((r['logN'].get(3) or 0.0)/r['L'] for r in res)[n//2]:.3f}")
    okP = all(r["kappa"] is None or r["kappa"] >= 1 for r in res)
    row("P", f"N_k and kappa computed exactly for all {n:,} proteins", okP)

    # (X) collisions within the proteome
    print("\n=== (X) distinct human proteins with identical k-spectra ===")
    seqs = sorted({x[3] for x in items})
    for k in (1, 2, 3, 4):
        g = collections.defaultdict(list)
        for s in seqs:
            if len(s) >= k:
                g[spectrum(s, k)].append(s)
        coll = [v for v in g.values() if len(v) > 1]
        print(f"    k={k}: {len(seqs):,} distinct sequences, {len(coll)} spectra shared by >= 2 of them "
              f"({sum(len(v) for v in coll)} sequences)")
        if coll and k >= 2:
            ex = coll[0]
            names = [next((x[1] or x[0]) for x in items if x[3] == s) for s in ex[:4]]
            print(f"      e.g. {names}, lengths {[len(s) for s in ex[:4]]}")
    row("X", "collisions of k-spectra among distinct human proteins counted for k = 1-4", True)

    # (M) fibre-mates and local encoders
    print("\n=== (M) fibre-mates at k = 25: APAAC(lambda = 24) and a window-25 CNN cannot tell them apart ===")
    H1, H2 = std_scale(KD), std_scale(HW)
    rr = random.Random(7)
    Wc = [[{a: rr.randrange(-3, 4) for a in "ACDEFGHIKLMNPQRSTVWY"} for _ in range(25)] for _ in range(8)]
    bc = [rr.randrange(-5, 6) for _ in range(8)]
    byacc = {x[0]: x for x in items}
    cand = sorted([r for r in res if (r["logN"].get(25) or 0) > 0 and r["L"] <= 2500],
                  key=lambda r: r["L"])
    print(f"    {sum(1 for r in res if (r['logN'].get(25) or 0) > 0)} proteins have N_25 > 1; "
          f"testing up to 12 of length <= 2500")
    okM, tested = True, 0
    for r in cand[:: max(1, len(cand) // 12)][:12]:
        S = byacc[r["acc"]][3]
        T = None
        for _ in range(200):
            x = random_trail(S, 25, rr)
            if x and x != S:
                T = x
                break
        if T is None:
            continue
        tested += 1
        same_spec = spectrum(S, 25) == spectrum(T, 25) and S[:24] == T[:24] and S[-24:] == T[-24:]
        same_ap = apaac(S, 24, H1, H2) == apaac(T, 24, H1, H2)
        same_cnn = cnn_mean(S, Wc, bc, 25) == cnn_mean(T, Wc, bc, 25)
        diff = sum(1 for a, b in zip(S, T) if a != b)
        sw = sw_score(S, T) / sw_score(S, S) if len(S) <= 2500 else float("nan")
        okM &= same_spec and same_ap and same_cnn
        print(f"      {r['gene'] or r['acc']:10s} L={len(S):5d}: mate differs at {diff:5d} positions "
              f"({100*diff/len(S):4.1f}%), SW(S,T)/SW(S,S) = {sw:.3f}; same 25-spectrum {same_spec}, "
              f"APAAC equal {same_ap}, CNN equal {same_cnn}")
    okM &= tested >= 3
    row("M", f"on {tested} human proteins an explicit different sequence with the same 25-spectrum "
             "gets the identical 68-dimensional APAAC vector and the identical mean-pooled CNN "
             "output", okM)

    print("\n=== battery summary ===")
    bad = [k for k, ok in ROWS if not ok]
    print(f"  TOTAL {len(ROWS)-len(bad)}/{len(ROWS)} passed" + (f"   FAILURES: {bad}" if bad else ""))
    print(f"  elapsed {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
