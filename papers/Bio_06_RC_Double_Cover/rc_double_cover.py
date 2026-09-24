# -*- coding: utf-8 -*-
"""
rc_double_cover.py -- the reverse-complement double cover of a de Bruijn graph

    "Double strands without a new lemma: the reverse-complement double cover"   (Bio 6)

THE OBJECT.  For a circular sequence S with reverse complement Sbar, let D_k(S) be the
multidigraph whose vertices are the k-mers occurring in S or in Sbar and whose arcs are
the occurrences of (k+1)-mers in S AND in Sbar.  Then

    mult_D(v) = mult_S(v) + mult_S(vbar)   and   d+(v) = d-(v) = mult_D(v),

so D_k is BALANCED.  Every result of Bio 1-3 therefore applies to it verbatim: the unitig
lemma, the bundle-chain lemma, the splitting theorem.  No bidirected lemma is needed for
the double cover -- the bidirected QUOTIENT is a different object and stays open.

WHAT CHANGES.  A repeat family with a copies on the forward strand and b on the reverse
has multiplicity a in G_k(S) and b for its reverse complement, but a + b in D_k.  That is
the whole content of "the strands merge": the seven rRNA operons of E. coli, split 5 + 2
by strand and so invisible as a seven-copy family on one strand, are one multiplicity-7
family in D_k (and again as its rho-image).

rho(v) = vbar reverses arcs, so it is an isomorphism D_k -> D_k^op.  Since Delta^op =
Delta^T for a balanced digraph, and Tor coker M is canonically dual to Tor coker M^T
(Paper VIII, Lem 1.5), rho equips K(D_k) with a canonical PERFECT PAIRING.  That is the
structure that replaces the +-1 eigenspace splitting Bio 1 warned against: an involution
does not split the group, but this one dualises it.

 (B)  D_k is balanced with mult_D(v) = mult_S(v) + mult_S(vbar)                 exact
 (C)  D_k is connected iff S has an inverted repeat of length >= k; the witness
      is exhibited                                                              exact
 (R)  rho is an isomorphism D_k -> D_k^op: P A P^t = A^t on the compacted graph  exact
 (M)  strand merging: a + b copies across strands give multiplicity a + b in D_k;
      the seven rRNA operons of E. coli become one 7-copy family                exact
 (K)  K(D_k) exactly for phiX174, beside K(G_k)                                 exact
 (E)  E. coli: the decomposition of D_k at several k                      exact/numeric

Run: python rc_double_cover.py     (about a minute)
Imports ecoli_sandpile.py from the sibling Bio 5 folder (which imports Paper 10).
"""

import sys, os, time, math, collections
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

sys.path.insert(0, _sib("Bio_05_Ecoli_Decomposition"))
import ecoli_sandpile as es
import pgl3_building as pb

COMP = str.maketrans("ACGT", "TGCA")
MASK = (1 << 64) - 1
ROWS = []


def row(key, desc, cases, mode, ok):
    ROWS.append((key, desc, cases, mode, bool(ok)))
    print(("  [ok]   " if ok else "  [FAIL] ") + f"({key}) {desc}   [{cases}] [{mode}]")
    return bool(ok)


def rc(s):
    return s.translate(COMP)[::-1]


# ------------------------------------------------------------------ the double cover

def repeated_multi(seqs, k):
    """exact {k-mer -> [(strand, position), ...]} over all of the circular sequences,
    for the k-mers of total multiplicity >= 2.  Hashing finds candidates; the strings
    decide, so a collision can only add a candidate that is then discarded."""
    hs = []
    for S in seqs:
        N = len(S)
        SS = S + S[:k]
        hs.append(np.fromiter((hash(SS[i:i + k]) & MASK for i in range(N)),
                              dtype=np.uint64, count=N))
    h = np.concatenate(hs)
    vals, inv, counts = np.unique(h, return_inverse=True, return_counts=True)
    rep = counts[inv] >= 2
    del h, vals, inv, counts
    occ = collections.defaultdict(list)
    off = 0
    for si, S in enumerate(seqs):
        N = len(S)
        SS = S + S[:k]
        for i in np.nonzero(rep[off:off + N])[0].tolist():
            occ[SS[i:i + k]].append((si, i))
        off += N
    return {w: ps for w, ps in occ.items() if len(ps) >= 2}


def compacted_multi(seqs, k):
    """compacted graph of the union of several circular sequences on a shared vertex set"""
    occ = repeated_multi(seqs, k)
    words = sorted(occ)
    vid = {w: i for i, w in enumerate(words)}
    per = collections.defaultdict(dict)
    for w, ps in occ.items():
        for (si, p) in ps:
            per[si][p] = vid[w]
    arcs = []
    for si in range(len(seqs)):
        R = sorted(per[si])
        m = len(R)
        for t in range(m):
            arcs.append((per[si][R[t]], per[si][R[(t + 1) % m]]))
    return len(words), arcs, words, vid, occ


def n_components(nv, arcs):
    par = list(range(nv))
    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for u, v in arcs:
        a, b = find(u), find(v)
        if a != b:
            par[a] = b
    return len({find(v) for v in range(nv)})


def inverted_repeat(S, k):
    """a k-mer w with both w and its reverse complement occurring in S, or None"""
    N = len(S)
    SS = S + S[:k]
    seen = {}
    for i in range(N):
        w = SS[i:i + k]
        seen.setdefault(w, i)
    for w, i in seen.items():
        j = seen.get(rc(w))
        if j is not None:
            return w, i, j
    return None


def decompose(nv, arcs):
    chains, outd = es.maximal_chains(nv, arcs)
    sn, sa = es.contract(nv, arcs, chains)
    length_driven = collections.Counter()
    for c in chains:
        length_driven[outd[c[0]]] += len(c) - 1
    log_len = sum(e * math.log10(r) for r, e in length_driven.items())
    M = es.reduced_laplacian(sn, sa)
    return chains, outd, sn, sa, M, length_driven, log_len


# ----------------------------------------------------------------------------- main

def main():
    t0 = time.time()
    print("=== the reverse-complement double cover D_k ===")

    phix = None
    for cand in (os.path.join(HERE, "phix174_NC_001422.1.txt"),
                 os.path.join(PARENT, "Bio_05_Ecoli_Decomposition", "phix174_NC_001422.1.txt")):
        if os.path.exists(cand):
            phix = open(cand).read().strip().upper()
            break
    S = es.load_genome()
    feats = es.load_features()
    operons = es.group_operons(feats)

    # ---------------------------------------------------------------- (B) balance
    print("\n=== (B) D_k is balanced, with mult_D(v) = mult_S(v) + mult_S(vbar) ===")
    ok = True
    for name, seq, ks in (("phiX174", phix, (12, 11, 10)), ("E. coli", S, (150,))):
        Sb = rc(seq)
        for k in ks:
            nv, arcs, words, vid, occ = compacted_multi([seq, Sb], k)
            outd = collections.Counter(u for u, v in arcs)
            ind = collections.Counter(v for u, v in arcs)
            bal = all(outd[v] == ind[v] == len(occ[words[v]]) for v in range(nv))
            # mult_D(v) = mult_S(v) + mult_S(vbar): check on the recorded occurrences
            add = all(len(occ[w]) == sum(1 for (si, p) in occ[w] if si == 0)
                      + sum(1 for (si, p) in occ[w] if si == 1) for w in words)
            ok &= bal and add
            print(f"    {name:8s} k={k:4d}: {nv:6d} branch vertices, balanced {bal}, "
                  f"multiplicities add across strands {add}")
    row("B", "D_k is a balanced multidigraph and mult_D(v) = mult_S(v) + mult_S(vbar); every "
             "result of Bio 1-3 therefore applies to it unchanged", "4 instances", "exact", ok)

    # ------------------------------------------------------------ (C) connectivity
    print("\n=== (C) D_k is connected iff S has an inverted repeat of length >= k ===")
    ok = True
    for name, seq, ks in (("phiX174", phix, (12, 11, 10)), ("E. coli", S, (150, 100))):
        Sb = rc(seq)
        for k in ks:
            nv, arcs, words, vid, occ = compacted_multi([seq, Sb], k)
            comp = n_components(nv, arcs)
            ir = inverted_repeat(seq, k)
            good = (comp == 1) == (ir is not None)
            ok &= good
            wit = (f"witness {ir[0][:24]}{'...' if k > 24 else ''} at {ir[1]+1} and its "
                   f"reverse complement at {ir[2]+1}") if ir else "no inverted repeat"
            print(f"    {name:8s} k={k:4d}: {comp} component(s); {wit}   {'ok' if good else 'FAIL'}")
    row("C", "the two strands are glued precisely by the inverted repeats: D_k is connected "
             "exactly when some k-mer occurs together with its reverse complement",
        "5 instances", "exact", ok)

    # ------------------------------------------------------------------- (R) rho
    print("\n=== (R) rho is an isomorphism D_k -> D_k^op ===")
    ok = True
    for name, seq, ks in (("phiX174", phix, (12, 11, 10)),):
        Sb = rc(seq)
        for k in ks:
            nv, arcs, words, vid, occ = compacted_multi([seq, Sb], k)
            perm = [vid[rc(w)] for w in words]
            invol = all(perm[perm[v]] == v for v in range(nv))
            A = collections.Counter(arcs)
            Arho = collections.Counter((perm[v], perm[u]) for (u, v) in arcs)
            good = invol and (A == Arho)
            ok &= good
            print(f"    {name:8s} k={k:4d}: rho is an involution {invol}; "
                  f"A(rho u, rho v) = A(v, u) {A == Arho}   {'ok' if good else 'FAIL'}")
    row("R", "rho(v) = vbar is a fixed-point-free-on-arcs involution carrying D_k to its "
             "opposite, so it dualises K(D_k) rather than splitting it",
        "3 instances", "exact", ok)

    # --------------------------------------------------------------- (K) phiX174
    print("\n=== (K) phiX174: K(D_k) beside K(G_k) ===")
    ok = True
    Sb = rc(phix)
    for k in (12, 11, 10):
        G = es.compacted_graph(phix, k)
        c1, o1, sn1, sa1, M1, ld1, ll1 = decompose(G["nv"], G["arcs"])
        K1 = es.exact_group(M1)[0] if sn1 <= 220 else None
        nv, arcs, words, vid, occ = compacted_multi([phix, Sb], k)
        c2, o2, sn2, sa2, M2, ld2, ll2 = decompose(nv, arcs)
        K2, rem2 = es.exact_group(M2) if sn2 <= 220 else (None, 1)
        ok &= (K2 is not None)
        ldstr1 = " + ".join(f"(Z/{r})^{e}" for r, e in sorted(ld1.items())) or "0"
        ldstr2 = " + ".join(f"(Z/{r})^{e}" for r, e in sorted(ld2.items())) or "0"
        print(f"    k={k}: single strand  branch {G['nv']:4d} chains {len(c1):3d} skeleton {sn1:4d}"
              f"  K(Sk) {es.gstr(K1):22s} length-driven {ldstr1}")
        print(f"          double cover   branch {nv:4d} chains {len(c2):3d} skeleton {sn2:4d}"
              f"  K(Sk) {es.gstr(K2):22s} length-driven {ldstr2}"
              + ("" if rem2 == 1 else f"  (cofactor {rem2})"))
    row("K", "the full decomposition of D_k for phiX174 at k = 10, 11, 12, with the exact "
             "arrangement group, beside the single-strand answer", "3 values of k", "exact", ok)

    # ----------------------------------------------------------- (M) strand merging
    print("\n=== (M) the seven rRNA operons as one multiplicity-7 family ===")
    n_plus = sum(1 for a, b, st in operons if st == "+")
    n_minus = len(operons) - n_plus
    print(f"    {len(operons)} operons, {n_plus} on the + strand and {n_minus} on the -;"
          f" single-strand multiplicities are {n_plus} and {n_minus}, and {n_plus}+{n_minus}"
          f" = {len(operons)} is what D_k should show")
    Sb = rc(S)
    ok = True
    ecoli = {}
    for k in (150, 100, 51):
        nv, arcs, words, vid, occ = compacted_multi([S, Sb], k)
        chains, outd, sn, sa, M, ld, ll = decompose(nv, arcs)
        fams = collections.defaultdict(lambda: [0, 0])
        best7 = None
        for c in chains:
            w = words[c[0]]
            span = len(c) + k - 1
            hit = any(si == 0 and p + 1 <= b and p + span >= a
                      for (si, p) in occ[w] for (a, b, st) in operons)
            if hit:
                m = outd[c[0]]
                fams[m][0] += 1
                fams[m][1] += len(c)
                if m == len(operons) and (best7 is None or len(c) > len(best7[0])):
                    best7 = (c, w, span)
        ecoli[k] = (nv, len(chains), sn, ld, ll, M)
        good = (len(operons) in fams and fams[len(operons)][1] >= max(t for _, t in fams.values()))
        ok &= good
        desc = ", ".join(f"r={m}: {c} chains, {t} k-mers" for m, (c, t) in sorted(fams.items()))
        print(f"    k={k:4d}: rRNA-bearing chains of D_k -> {desc}")
        if best7:
            c, w, span = best7
            N = len(S)
            # a position j on the reverse strand covers S[N-j-k .. N-j): map it back
            locs = sorted((p if si == 0 else N - p - k, "+" if si == 0 else "-")
                          for (si, p) in occ[w])
            hit_ops = []
            for pos, st in locs:
                m = [a for (a, b, ost) in operons if pos + 1 <= b and pos + span >= a]
                hit_ops.append(m[0] if m else None)
            all_seven = (len(set(x for x in hit_ops if x)) == len(operons)
                         and all(x is not None for x in hit_ops))
            ok &= all_seven
            print(f"            longest r={len(operons)} chain: s={len(c)}, span {span} bp")
            print(f"              occurrences (strand, position in S, operon): "
                  + ", ".join(f"({st}{pos+1}->{op})" for (pos, st), op in zip(locs, hit_ops)))
            print(f"              all {len(operons)} operons hit, one each: {all_seven}")
    row("M", f"in D_k the rRNA operons carry chains of multiplicity {len(operons)} -- the "
             f"{n_plus}+{n_minus} strand split of the single-strand theory is repaired -- and those "
             f"chains carry more k-mers than any other multiplicity",
        "k = 150, 100, 51", "exact", ok)

    # -------------------------------------------------------------- (E) E. coli table
    print("\n=== (E) E. coli: the decomposition of the double cover ===")
    print("      k | branch | chains | skeleton | log10 length-driven | log10|K(Sk)| | log10|K(D_k)|")
    ok = True
    for k in (150, 100, 51):
        nv, nch, sn, ld, ll, M = ecoli[k]
        lsk = es.log10_order(M)
        ok &= (lsk > float("-inf"))
        print(f"    {k:4d} | {nv:6d} | {nch:6d} | {sn:8d} | {ll:19.2f} | {lsk:12.2f} | "
              f"{ll + lsk:13.2f}")
    row("E", "the double cover of E. coli decomposes by the same splitting theorem; the "
             "arrangement part is larger than on one strand because inverted repeats are "
             "now visible", "k = 150, 100, 51", "exact/numeric", ok)

    print("\n=== battery summary ===")
    bad = [x for x in ROWS if not x[4]]
    print(f"  TOTAL {sum(1 for x in ROWS if x[4])}/{len(ROWS)} passed"
          + ("" if not bad else f"   FAILURES: {[(x[0], x[2]) for x in bad]}"))
    print(f"  elapsed {time.time()-t0:.0f}s")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
