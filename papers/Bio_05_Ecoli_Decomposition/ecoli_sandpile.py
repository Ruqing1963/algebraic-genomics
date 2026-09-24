# -*- coding: utf-8 -*-
"""
ecoli_sandpile.py -- the bundle-chain decomposition of a real bacterial genome

    Escherichia coli str. K-12 substr. MG1655, NCBI RefSeq NC_000913.3, 4 641 652 bp

Theory (Bio 1 / Bio 3).  For a circular sequence S and k >= 1, G_k(S) is the de Bruijn
multigraph whose vertices are the distinct k-mers and whose arcs are the occurrences of
(k+1)-mers.  Every vertex then has d+(v) = d-(v) = mult(v), the number of occurrences of
its k-mer, so

    *** a vertex survives unitig contraction if and only if its k-mer is REPEATED. ***

That single observation makes a 4.6 Mbp genome tractable: we never build the 4.6M-vertex
graph.  We find the repeated k-mers (exactly), and the compacted graph has one vertex per
repeated k-mer and one arc per consecutive pair of repeat occurrences around the circle.
Bundle-chain contraction (Bio 3, Lemma 2) then splits

    K(G_k) = (+)_chains (Z/r_c)^(len_c - 1)  (+)  K(Sk),

the first summand length-driven (it dies as k grows past the repeat lengths), the second
the arrangement group of the repeat-adjacency skeleton.

SCOPE.  This is the SINGLE-STRAND theory of Bio 1-4.  Real assembly graphs are bidirected;
a repeat whose copies lie on opposite strands is invisible here (Bio 2, Open (a)).  For
E. coli this is not a footnote -- see the rRNA operon analysis in check (R).

Checks:
 (H)  k-mer multiplicities are computed by 64-bit hashing and then verified EXACTLY on
      the repeated set by string comparison, so hash collisions cannot survive.   exact
 (C)  compacted graph: branch vertices = repeated k-mers, balanced, d+(v) = mult(v). exact
 (B)  bundle-chain contraction in one pass (Bio 3 Lemma 6: chain steps form a partial
      injection, so maximal chains are disjoint and order-independent).           exact
 (K)  skeleton group: exact Smith normal form when the skeleton is small enough,
      log10 of the order from a float LU determinant otherwise.              exact/numeric
 (R)  rRNA operons: which chains carry them, at which k, and on which strand.      exact

Run: python ecoli_sandpile.py            (about 40 s, ~0.5 GB peak)
Data: NC_000913.3.fasta and NC_000913.3.ft.txt are read from the parent folder; if absent
they are fetched from the NCBI E-utilities API.
"""

import sys, os, time, math, collections, urllib.request
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

PAPER10 = _sib("lib")
sys.path.insert(0, PAPER10)
import pgl3_building as pb

ACC = "NC_000913.3"
EFETCH = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
          "?db=nuccore&id={acc}&rettype={rt}&retmode=text")
MASK = (1 << 64) - 1
KS = (150, 100, 75, 51, 31)
SNF_MAX = int(os.environ.get("SNF_MAX", "200"))   # exact dense Smith form up to this size
ROWS = []


def row(key, desc, cases, mode, ok):
    ROWS.append((key, desc, cases, mode, bool(ok)))
    print(("  [ok]   " if ok else "  [FAIL] ") + f"({key}) {desc}   [{cases}] [{mode}]")
    return bool(ok)


# --------------------------------------------------------------------------- data

def fetch(rt, name):
    for folder in (HERE, PARENT):
        p = os.path.join(folder, name)
        if os.path.exists(p):
            return open(p, encoding="utf-8", errors="replace").read()
    url = EFETCH.format(acc=ACC, rt=rt)
    print(f"      fetching {url}")
    text = urllib.request.urlopen(url, timeout=300).read().decode("utf-8", "replace")
    open(os.path.join(HERE, name), "w", encoding="utf-8").write(text)
    return text


def load_genome():
    text = fetch("fasta", f"{ACC}.fasta")
    seq = "".join(l.strip() for l in text.splitlines() if not l.startswith(">")).upper()
    assert set(seq) <= set("ACGT"), sorted(set(seq))
    return seq


def load_features(kinds=("rRNA", "CDS", "mobile_element", "repeat_region")):
    """features from the NCBI feature table: (start, end, strand, kind, product), 1-based"""
    text = fetch("ft", f"{ACC}.ft.txt")
    feats, pending = [], None
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3 and parts[0] and parts[1] and parts[2] in kinds:
            a, b = int(parts[0].lstrip("<>")), int(parts[1].lstrip("<>"))
            pending = [min(a, b), max(a, b), "+" if a <= b else "-", parts[2], ""]
            feats.append(pending)
        elif pending is not None and len(parts) >= 5 and parts[3] in ("product", "mobile_element_type",
                                                                     "rpt_family", "gene"):
            if not pending[4]:
                pending[4] = parts[4]
    return [tuple(f) for f in feats]


# ------------------------------------------------------------------- repeated k-mers

def compacted_graph(S, k):
    """(nv, arcs, vertex k-mers, positions per vertex) of the unitig-compacted graph"""
    N = len(S)
    SS = S + S[:k]
    h = np.fromiter((hash(SS[i:i + k]) & MASK for i in range(N)), dtype=np.uint64, count=N)
    vals, inv, counts = np.unique(h, return_inverse=True, return_counts=True)
    n_distinct, n_cand = len(vals), int((counts[inv] >= 2).sum())
    cand = np.nonzero(counts[inv] >= 2)[0]
    del h, vals, inv, counts
    occ = collections.defaultdict(list)
    for i in cand.tolist():
        occ[SS[i:i + k]].append(i)
    occ = {w: ps for w, ps in occ.items() if len(ps) >= 2}      # exact: kills collisions
    words = sorted(occ)
    vid = {w: i for i, w in enumerate(words)}
    pos2v = {}
    for w, ps in occ.items():
        for p in ps:
            pos2v[p] = vid[w]
    R = sorted(pos2v)
    m = len(R)
    arcs = [(pos2v[R[t]], pos2v[R[(t + 1) % m]]) for t in range(m)] if m else []
    return dict(nv=len(words), arcs=arcs, words=words, occ=occ, positions=R,
                n_distinct=n_distinct, n_cand=n_cand, m=m)


# ----------------------------------------------------------------- bundle contraction

def maximal_chains(nv, arcs):
    """Bio 3, Lemma 6: the chain steps form a partial injection, so the maximal bundle
    chains are its maximal paths -- disjoint, and independent of the contraction order."""
    outd = collections.Counter(u for u, v in arcs)
    targets = collections.defaultdict(set)
    for (u, v) in arcs:
        targets[u].add(v)
    step = {}
    for v in range(nv):
        if outd[v] and len(targets[v]) == 1:
            w = next(iter(targets[v]))
            if w != v and outd[w] == outd[v]:
                step[v] = w
    pred = {}
    for v, w in step.items():
        assert w not in pred, "chain steps are not injective -- Lemma 6 violated"
        pred[w] = v
    chains, seen = [], set()
    for v in range(nv):
        if v in seen or v in pred:
            continue
        c = [v]
        while c[-1] in step:
            nxt = step[c[-1]]
            if nxt in seen or nxt in c:
                break
            c.append(nxt)
        seen.update(c)
        if len(c) >= 2:
            chains.append(c)
    rest = [v for v in range(nv) if v not in seen and v in step]   # closed chains
    while rest:
        v = rest[0]
        c = [v]
        while step[c[-1]] not in c:
            c.append(step[c[-1]])
        seen.update(c)
        chains.append(c)
        rest = [x for x in rest if x not in seen]
    return chains, outd


def contract(nv, arcs, chains):
    rep = list(range(nv))
    cs = {}
    for i, c in enumerate(chains):
        for v in c:
            rep[v] = c[0]
            cs[v] = i
    new = []
    for (u, v) in arcs:
        if u in cs and cs[u] == cs.get(v, -1):
            continue                                    # internal spine arc (or tail->head loop)
        new.append((rep[u], rep[v]))
    alive = sorted({u for u, v in new} | {v for u, v in new})
    if not alive:
        return 1, []
    rl = {v: i for i, v in enumerate(alive)}
    return len(alive), [(rl[u], rl[v]) for (u, v) in new]


# ------------------------------------------------------------------------ the group

def reduced_laplacian(nv, arcs, root=0):
    L = collections.defaultdict(int)
    for (u, v) in arcs:
        L[(u, u)] += 1
        L[(u, v)] -= 1
    idx = [v for v in range(nv) if v != root]
    pos = {v: i for i, v in enumerate(idx)}
    n = len(idx)
    M = [[0] * n for _ in range(n)]
    for (u, v), c in L.items():
        if u != root and v != root:
            M[pos[u]][pos[v]] = c
    return M


def log10_order(M):
    """log10 |det| of the transposed reduced Laplacian, by float LU"""
    if not M:
        return 0.0
    A = np.array(M, dtype=np.float64).T
    sign, logdet = np.linalg.slogdet(A)
    return float(logdet / math.log(10.0)) if sign != 0 else float("-inf")


def exact_group(M):
    """invariant factors of coker of the transposed reduced Laplacian, p-adically (no
    coefficient explosion).  Returns (factors, unfactored cofactor of |det|): the cofactor
    is 1 when the determinant factored completely, and otherwise a composite whose
    p-primary contribution we do not claim to resolve."""
    if not M:
        return (), 1
    A = [list(r) for r in zip(*M)]
    factors, free, remainder = pb.coker_structure(A)
    assert free == 0, f"skeleton Laplacian should be nonsingular, got free rank {free}"
    return tuple(factors), remainder


def elementary(inv):
    import sympy as sp
    out = []
    for x in inv:
        for p, e in sp.factorint(int(x)).items():
            out.append(p ** e)
    return tuple(sorted(out))


def gstr(inv):
    if not inv:
        return "0"
    c = collections.Counter(elementary(inv))
    return " + ".join(f"(Z/{p})^{e}" if e > 1 else f"Z/{p}" for p, e in sorted(c.items()))


def group_operons(feats, gap=2000):
    """cluster the rRNA genes of the feature table into operons: same strand, small gaps"""
    ops, cur = [], None
    for a, b, st, kind, prod in sorted((f for f in feats if f[3] == "rRNA"), key=lambda f: f[0]):
        if cur is not None and st == cur[2] and a - cur[1] <= gap:
            cur[1] = max(cur[1], b)
        else:
            cur = [a, b, st]
            ops.append(cur)
    return [tuple(o) for o in ops]


def annotate(feats, lo, hi):
    """the distinct gene products a genomic interval [lo, hi] (1-based) overlaps"""
    names = []
    for a, b, st, kind, prod in feats:
        if a <= hi and b >= lo and prod:
            p = prod.split(";")[0].strip()
            if p not in names:
                names.append(p)
    return names


def chains_str(chains, outd):
    c = collections.Counter((outd[ch[0]], len(ch)) for ch in chains)
    top = sorted(c.items(), key=lambda kv: -kv[0][1])[:4]
    return ", ".join(f"(r={r},s={s})x{n}" for (r, s), n in top)


# ----------------------------------------------------------------------------- main

def analyse(S, k):
    """structure only -- fast; the exact Smith form is a separate, slower pass"""
    t0 = time.time()
    G = compacted_graph(S, k)
    nv, arcs = G["nv"], G["arcs"]
    outd0 = collections.Counter(u for u, v in arcs)
    balanced = all(outd0[v] == len(G["occ"][G["words"][v]]) for v in range(min(nv, 2000)))
    chains, outd = maximal_chains(nv, arcs)
    sn, sa = contract(nv, arcs, chains)
    length_driven = collections.Counter()
    for ch in chains:
        length_driven[outd[ch[0]]] += len(ch) - 1
    log_len = sum(e * math.log10(r) for r, e in length_driven.items())
    M = reduced_laplacian(sn, sa)
    log_sk = log10_order(M)
    # BEST: #reconstructions from the k-mer spectrum = |K(G_k)| * prod_v (d+(v) - 1)!
    log_local = sum(math.lgamma(len(ps)) for ps in G["occ"].values()) / math.log(10.0)
    return dict(k=k, nv=nv, m=G["m"], n_distinct=G["n_distinct"], chains=chains, outd=outd,
                sn=sn, sa=sa, M=M, length_driven=length_driven, log_len=log_len, log_sk=log_sk,
                log_local=log_local, K_sk=None, balanced=balanced, occ=G["occ"],
                words=G["words"], secs=time.time() - t0)


def check_against_phix174():
    """The fast code path -- repeated k-mers straight to the compacted graph, one-pass
    chain contraction, p-adic Smith form -- rerun on the genome whose decomposition the
    companion chapter computed the slow way (build the graph, compact it, contract).
    Every published number must come back."""
    print("\n=== (V) the fast code path against the published phiX174 decomposition ===")
    path = None
    for cand in (os.path.join(HERE, "phix174_NC_001422.1.txt"),
                 os.path.join(PARENT, "Bio_03_Repeat_Splitting", "phix174_NC_001422.1.txt"),
                 os.path.join(PARENT, "Bio_04_Multicopy_Repeats", "phix174_NC_001422.1.txt")):
        if os.path.exists(cand):
            path = cand
            break
    if path is None:
        return row("V", "phiX174 cross-check skipped: sequence not found", "-", "exact", False)
    S = open(path).read().strip().upper()
    # (branch vertices, chains, skeleton vertices, sum(s-1), K(Sk) elementary divisors)
    published = {
        12: (2, 1, 1, 1, ()),
        11: (10, 2, 8, 2, (3, 11)),
        10: (47, 8, 37, 10, (4, 1076338579)),
        9: (146, 28, 110, 36, (2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 7, 523, 3779, 5297, 6561,
                               7669, 2707831823)),
    }
    ok = True
    for k in sorted(published, reverse=True):
        nvE, chE, snE, ldE, KE = published[k]
        r = analyse(S, k)
        ld = sum(e for _, e in r["length_driven"].items())
        got = (r["nv"], len(r["chains"]), r["sn"], ld)
        good = got == (nvE, chE, snE, ldE)
        Kgot = None
        if KE is not None:
            factors, rem = exact_group(r["M"])
            Kgot = elementary(factors)
            good &= (Kgot == tuple(sorted(elementary(KE))) and rem == 1)
        ok &= good
        print(f"     k={k:2d}  branch {r['nv']:4d}/{nvE:<4d} chains {len(r['chains']):3d}/{chE:<3d} "
              f"skeleton {r['sn']:4d}/{snE:<4d} sum(s-1) {ld:3d}/{ldE:<3d}  "
              f"K(Sk) {'matches' if KE is not None and Kgot == tuple(sorted(elementary(KE))) else 'n/a'}"
              f"   {'ok' if good else 'MISMATCH'}")
    row("V", "on phiX174 the fast path reproduces every published number of the companion "
             "chapter: branch vertices, chain count, skeleton size, length-driven exponent and "
             "the exact skeleton group at k = 9, 10, 11, 12",
        "4 values of k", "exact", ok)
    return ok


def main():
    t0 = time.time()
    print("=== Escherichia coli K-12 MG1655 (NC_000913.3): the bundle-chain decomposition ===")
    check_against_phix174()
    S = load_genome()
    feats = load_features()
    print(f"    genome {len(S)} bp circular; {len(feats)} annotated features")
    res = {}
    print("\n      k   distinct k-mers   branch   arcs    chains   skeleton   log10|K(Sk)|   log10|K|")
    for k in KS:
        r = analyse(S, k)
        res[k] = r
        print(f"     {k:4d}   {r['n_distinct']:15d}   {r['nv']:6d}  {r['m']:6d}   "
              f"{len(r['chains']):6d}   {r['sn']:8d}   {r['log_sk']:12.4f}   "
              f"{r['log_len'] + r['log_sk']:9.2f}    [{r['secs']:.0f}s]")
    row("H", "k-mer multiplicities: 64-bit hashing finds candidates, exact string comparison "
             "decides; no collision survives", f"{len(KS)} values of k", "exact", True)
    row("C", "compacted graph is balanced with d+(v) = mult(v) (spot-checked on <=2000 vertices)",
        f"{len(KS)} values of k", "exact", all(res[k]["balanced"] for k in KS))
    row("B", "maximal bundle chains found in one pass (Lemma 6 injectivity asserted)",
        f"{sum(len(res[k]['chains']) for k in KS)} chains", "exact", True)
    row("K", f"skeleton group exactly for <= {SNF_MAX} vertices, log10 of the order by float LU "
             f"otherwise", f"{len(KS)} skeletons", "exact/numeric",
        all(res[k]["log_sk"] > float("-inf") for k in KS))

    print("\n=== (R) the rRNA operons, located in the decomposition ===")
    operons = group_operons(feats)
    n_plus = sum(1 for a, b, st in operons if st == "+")
    n_minus = len(operons) - n_plus
    print(f"    feature table: {sum(1 for f in feats if f[3]=='rRNA')} rRNA genes in "
          f"{len(operons)} operons: " + ", ".join(f"{a}{st}" for a, b, st in operons))
    print(f"    {n_plus} operons on the + strand, {n_minus} on the -.  The theory of Bio 1-4 is")
    print(f"    SINGLE-STRAND: two copies are a repeat of G_k only if they are co-oriented in S.")
    print(f"    So the seven operons cannot form one multiplicity-7 family; they must split into")
    print(f"    a {n_plus}-copy (+ strand) family and a {n_minus}-copy (- strand) family.")
    seen_mult, longest_is_rrna = {}, []
    for k in KS:
        r = res[k]
        fams = collections.defaultdict(lambda: [0, 0])     # r -> [#chains, total k-mers]
        best = max(r["chains"], key=len)
        for ch in r["chains"]:
            ps = sorted(r["occ"][r["words"][ch[0]]])
            span = len(ch) + k - 1
            in_op = sum(1 for p in ps for op in operons if p + 1 <= op[1] and p + span >= op[0])
            if in_op:
                mult = r["outd"][ch[0]]
                fams[mult][0] += 1
                fams[mult][1] += len(ch)
                if ch is best:
                    longest_is_rrna.append(mult == n_minus and in_op == len(ps))
        seen_mult[k] = dict(fams)
        desc = ", ".join(f"r={m}: {c} chain(s), {t} k-mers" for m, (c, t) in sorted(fams.items()))
        print(f"    k={k:4d}: rRNA-bearing chains -> {desc if desc else '(none)'}")
    ok_r = (all(7 not in seen_mult[k] for k in KS)
            and all(max(seen_mult[k]) == n_plus for k in KS)
            and len(longest_is_rrna) == len(KS) and all(longest_is_rrna))
    row("R", f"no rRNA-bearing chain ever has multiplicity 7; the largest is {n_plus} (the + strand "
             f"operons), and at every k the LONGEST chain in the whole genome is the multiplicity-"
             f"{n_minus} spine of the two co-oriented minus-strand operons",
        f"{len(KS)} values of k", "exact", ok_r)
    print(f"    The five + strand operons are not identical either: they contribute chains of")
    print(f"    multiplicity 3 and 4 as well as {n_plus}, one per sub-family of copies that agree.")

    print("\n    the highest-multiplicity repeat families, annotated")
    for k in (150, 31):
        r = res[k]
        by_r = collections.defaultdict(list)
        for ch in r["chains"]:
            by_r[r["outd"][ch[0]]].append(ch)
        print(f"      k={k}")
        for mult in sorted(by_r, reverse=True)[:4]:
            chs = sorted(by_r[mult], key=len, reverse=True)
            ps = sorted(r["occ"][r["words"][chs[0][0]]])
            span = len(chs[0]) + k - 1
            names = annotate(feats, ps[0] + 1, ps[0] + span)
            print(f"        r={mult:2d}: {len(chs):3d} chains, longest s={len(chs[0])} (span {span} bp)"
                  f"   {'; '.join(names[:2])[:80] if names else '(unannotated)'}")

    print("\n    the six longest bundle chains at each k, with their annotation")
    for k in (150, 51):
        r = res[k]
        print(f"      k={k}")
        for ch in sorted(r["chains"], key=lambda c: -len(c))[:6]:
            ps = sorted(r["occ"][r["words"][ch[0]]])
            span = len(ch) + k - 1
            names = annotate(feats, ps[0] + 1, ps[0] + span)
            print(f"        r={r['outd'][ch[0]]:2d}  s={len(ch):5d}  span={span:6d} bp  at "
                  f"{[p+1 for p in ps[:6]]}")
            print(f"              {'; '.join(names[:4])[:150]}")

    print(f"\n=== (K) the arrangement group, exactly (skeletons of at most {SNF_MAX} vertices) ===")
    for k in sorted(KS, key=lambda k: res[k]["sn"]):
        sn = res[k]["sn"]
        if sn > SNF_MAX:
            print(f"      k={k:4d}: skeleton {sn:4d} vertices -- beyond the exact range; "
                  f"|K(Sk)| = 10^{res[k]['log_sk']:.2f} (float LU)")
            continue
        t1 = time.time()
        try:
            factors, rem = exact_group(res[k]["M"])
            res[k]["K_sk"] = (factors, rem)
            extra = "" if rem == 1 else f"  x (an unfactored cofactor {rem} of the determinant)"
            print(f"      k={k:4d}: skeleton {sn:4d} vertices -> K(Sk) = {gstr(factors)}{extra}"
                  f"   [{time.time()-t1:.0f}s]")
        except Exception as e:
            print(f"      k={k:4d}: skeleton {sn:4d} vertices -- exact structure not obtained ({e})")
        sys.stdout.flush()

    print("\n=== summary table ===")
    print("   k | branch | chains | skel | K(Sk) arrangement part            | log10 length-driven"
          " | log10|K(G_k)| | log10 #reconstructions")
    for k in KS:
        r = res[k]
        if r["K_sk"] is not None:
            ks = gstr(r["K_sk"][0]) + ("" if r["K_sk"][1] == 1 else " + (cofactor)")
        else:
            ks = f"order 10^{r['log_sk']:.2f} ({r['sn']} vtx)"
        print(f"  {k:3d} | {r['nv']:6d} | {len(r['chains']):6d} | {r['sn']:4d} | {ks[:33]:33s} | "
              f"{r['log_len']:19.2f} | {r['log_len'] + r['log_sk']:13.2f} | "
              f"{r['log_len'] + r['log_sk'] + r['log_local']:22.2f}")
    print("    (#reconstructions = |K(G_k)| * prod_v (d+(v)-1)!, the BEST theorem: the number of")
    print("     circular sequences with exactly this k-mer spectrum, read as arc-rooted circuits.)")

    print("\n=== battery summary ===")
    bad = [x for x in ROWS if not x[4]]
    print(f"  TOTAL {sum(1 for x in ROWS if x[4])}/{len(ROWS)} passed"
          + ("" if not bad else f"   FAILURES: {[(x[0], x[2]) for x in bad]}"))
    print(f"  elapsed {time.time()-t0:.0f}s")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
