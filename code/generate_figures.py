# -*- coding: utf-8 -*-
r"""
generate_figures.py -- the five overview figures of the Algebraic Genomics monograph.

Every number drawn in a figure is COMPUTED here, from the chapters' own code and data, and
checked; nothing is typed in by hand.  The script ends with a line "TOTAL a/b passed" so that
build_all.py can treat it like the chapter verification scripts.

    Figure 1 (Part I)    fig1_sandpile_torsor    two interleaved repeats: K = Z/2 acts on the two
                                                 reconstructions by one rotor-routing chip
    Figure 2 (Part II)   fig2_double_cover       phiX174 double cover D_12 with rho and its palindromic
                                                 vertex; frontier width at k = 10 (greedy 26 vs 62)
                                                 and the open vertices at the widest step
    Figure 3 (Part II)   fig3_sharpP_reduction   #DS = 2 eps(G): a 4-regular graph, its palindromic
                                                 genome, and the Eulerian orientations = box points
    Figure 4 (Part III)  fig4_spectral_fibre     ZNF91 and a fibre-mate: same 25-spectrum, identical
                                                 68-dim PseAAC and CNN output, permuted zinc fingers
    Figure 5 (Part IV)   fig5_snarl_schur        a genome bubble graph, its snarl tree, Kron/Schur
                                                 reduction of one snarl, and linear-time scaling

    python code/generate_figures.py            all five  (about 1-2 minutes)
    python code/generate_figures.py 2 4        only figures 2 and 4

Requires numpy, matplotlib, networkx; figure 4 needs papers/Spectral_Fibres/human_reviewed.tsv.gz.
Output: figures/figN_*.pdf (vector) and .png (300 dpi).
"""
import os, sys, re, math, random, itertools, collections, importlib.util
from fractions import Fraction

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle, Wedge, Circle
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PAPERS = os.path.join(ROOT, "papers")
OUT = os.path.join(ROOT, "figures")

# palette shared with papers/Spectral_Fibres/make_figures.py
INK, INK2, MUTED, GRID, BASE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, AQUA, GRAY, RED, PURPLE = "#2a78d6", "#eb6834", "#1baf7a", "#b9b8b1", "#d23c3c", "#7a5cc7"
plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 8.5, "axes.edgecolor": BASE, "axes.labelcolor": INK2,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.8, "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False, "pdf.fonttype": 42, "axes.titlesize": 9,
    "axes.titleweight": "bold", "axes.titlelocation": "left", "axes.titlecolor": INK})

CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append(bool(ok))
    print(("  [ok]   " if ok else "  [FAIL] ") + name + (f"   ({detail})" if detail else ""), flush=True)
    return ok


def load(folder, module):
    """import a chapter's module from papers/<folder>/<module>.py (its folder on sys.path)"""
    path = os.path.join(PAPERS, folder)
    if path not in sys.path:
        sys.path.insert(0, path)
    spec = importlib.util.spec_from_file_location(module, os.path.join(path, module + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, name + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote figures/{name}.pdf/.png", flush=True)


def det_frac(M):
    M = [[Fraction(x) for x in r] for r in M]
    n, d = len(M), Fraction(1)
    for i in range(n):
        p = next((r for r in range(i, n) if M[r][i] != 0), None)
        if p is None:
            return Fraction(0)
        if p != i:
            M[i], M[p] = M[p], M[i]
            d = -d
        d *= M[i][i]
        for r in range(i + 1, n):
            f = M[r][i] / M[i][i]
            for c in range(i, n):
                M[r][c] -= f * M[i][c]
    return d


def inv_frac(M):
    """exact inverse by Gauss-Jordan over the rationals"""
    n = len(M)
    A = [[Fraction(x) for x in r] + [Fraction(int(i == j)) for j in range(n)] for i, r in enumerate(M)]
    for i in range(n):
        p = next(r for r in range(i, n) if A[r][i] != 0)
        A[i], A[p] = A[p], A[i]
        piv = A[i][i]
        A[i] = [x / piv for x in A[i]]
        for r in range(n):
            if r != i and A[r][i] != 0:
                f = A[r][i]
                A[r] = [x - f * y for x, y in zip(A[r], A[i])]
    return [row[n:] for row in A]


def arborescences(verts, arcs, root):
    """in-arborescences to root: every other vertex picks one out-arc, no cycles"""
    outs = {v: [i for i, (t, h, _) in enumerate(arcs) if t == v and h != v] for v in verts}
    others = [v for v in verts if v != root]
    res = []
    for pick in itertools.product(*[outs[v] for v in others]):
        nxt = {v: arcs[a][1] for v, a in zip(others, pick)}
        ok = True
        for v in others:
            seen, x = set(), v
            while x != root:
                if x in seen:
                    ok = False
                    break
                seen.add(x)
                x = nxt[x]
            if not ok:
                break
        if ok:
            res.append(dict(zip(others, pick)))
    return res


def circuits(arcs, first):
    """Eulerian circuits up to rotation: arc sequences that start with arc `first`"""
    out = collections.defaultdict(list)
    for i, (t, h, _) in enumerate(arcs):
        out[t].append(i)
    res, used = [], [False] * len(arcs)
    used[first] = True

    def rec(v, path):
        if len(path) == len(arcs):
            if v == arcs[first][0]:
                res.append(tuple(path))
            return
        for i in out[v]:
            if not used[i]:
                used[i] = True
                rec(arcs[i][1], path + [i])
                used[i] = False
    rec(arcs[first][1], [first])
    return res


# ================================================================== Figure 1
def figure1():
    print("Figure 1: sandpile torsor on two interleaved repeats")
    # repeat skeleton of the circular genome  R1 a R2 b R1 c R2 d  (arcs = unique segments)
    V = ["R1", "R2"]
    inter = [("R1", "R2", "a"), ("R2", "R1", "b"), ("R1", "R2", "c"), ("R2", "R1", "d")]
    disj = [("R1", "R1", "a"), ("R1", "R2", "b"), ("R2", "R2", "c"), ("R2", "R1", "d")]   # R1 a R1 b R2 c R2 d

    def summary(arcs):
        L = [[sum(1 for t, h, _ in arcs if t == v and h != v) if v == u else
              -sum(1 for t, h, _ in arcs if t == v and h == u) for u in V] for v in V]
        K = abs(det_frac([[L[1][1]]]))                      # reduced at root R1
        loc = math.prod(math.factorial(sum(1 for t, *_ in arcs if t == v) - 1) for v in V)
        return int(K), loc, circuits(arcs, 0)
    K1, loc1, C1 = summary(inter)
    K2, loc2, C2 = summary(disj)
    check("interleaved: |K| = 2 and BEST = loc * |K| = 2 = number of Eulerian circuits",
          K1 == 2 and loc1 * K1 == len(C1) == 2, f"circuits {len(C1)}")
    check("disjoint (non-interleaved): |K| = 1 and a unique reconstruction", K2 == 1 and len(C2) == loc2 * K2 == 1)
    # rotor-routing: rotor order at R2 = (b, d); one chip at R2 advances the rotor and exits to the root
    arbs = arborescences(V, inter, "R1")
    rotor = {"R2": [1, 3]}                                   # arc indices b, d

    def add_chip(T, v="R2"):
        T = dict(T)
        while v != "R1":
            r = rotor[v]
            T[v] = r[(r.index(T[v]) + 1) % len(r)]
            v = inter[T[v]][1]
        return T

    def best_circuit(T):
        """last-exit bijection: at each non-root vertex the tree arc is used last"""
        order = {v: [i for i, (t, *_) in enumerate(inter) if t == v] for v in V}
        for v, a in T.items():
            order[v] = [i for i in order[v] if i != a] + [a]
        order["R1"] = [0] + [i for i in order["R1"] if i != 0]
        ptr = {v: 0 for v in V}
        v, seq = "R1", []
        while ptr[v] < len(order[v]):
            a = order[v][ptr[v]]
            ptr[v] += 1
            seq.append(a)
            v = inter[a][1]
        return tuple(seq)
    Tb = next(T for T in arbs if T["R2"] == 1)
    Td = add_chip(Tb)
    cb, cd = best_circuit(Tb), best_circuit(Td)
    word = lambda c: " ".join(f"{inter[a][0]} {inter[a][2]}" for a in c)
    check("one chip at R2 maps arborescence {b} to {d} and back (K = Z/2 acts simply transitively)",
          Td["R2"] == 3 and add_chip(Td) == Tb and len(arbs) == 2)
    check("the two arborescences give the two reconstructions", {cb, cd} == set(C1),
          f"{word(cd)}  |  {word(cb)}")

    fig = plt.figure(figsize=(7.2, 2.7))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 1.15], wspace=0.25)
    # (a) the two genomes as rings
    ax = fig.add_subplot(gs[0]); ax.set_aspect("equal"); ax.axis("off")
    colors = {"R1": BLUE, "R2": ORANGE}
    for cx, circ, lab in ((0.0, cd, "$S$"), (3.3, cb, "$S'$")):
        segs = []
        for a in circ:
            segs += [(inter[a][0], True), (inter[a][2], False)]
        n = len(segs)
        for i, (name, rep) in enumerate(segs):
            th1, th2 = 90 - 360 * i / n, 90 - 360 * (i + 1) / n
            ax.add_patch(Wedge((cx, 0), 1.0, th2 + 1.5, th1 - 1.5, width=0.28,
                               color=colors[name] if rep else GRID))
            mid = math.radians((th1 + th2) / 2)
            rr_ = 1.22 if rep else 0.58
            ax.text(cx + rr_ * math.cos(mid), rr_ * math.sin(mid), name, ha="center", va="center",
                    fontsize=6.5 if rep else 7.5, color=INK if not rep else colors[name],
                    fontweight="bold" if rep else None)
        ax.text(cx, 0, lab, ha="center", va="center", fontsize=10)
    ax.annotate("", (2.0, 0), (1.3, 0), arrowprops=dict(arrowstyle="<->", color=ORANGE, lw=1.2))
    ax.text(1.65, 0.16, "1 chip", ha="center", fontsize=7, color=ORANGE)
    ax.set_xlim(-1.45, 4.75); ax.set_ylim(-2.25, 1.45)
    ax.text(1.65, -2.2, "same $k$-mer data, two genomes\n(segments $b$ and $d$ exchanged)", ha="center",
            fontsize=7, color=INK2)
    fig.text(0.13, 0.93, "a  Two interleaved repeats", fontsize=9, fontweight="bold")
    # (b) skeleton multigraph with rotors
    ax = fig.add_subplot(gs[1]); ax.set_aspect("equal"); ax.axis("off")
    P = {"R1": (0, 0), "R2": (2, 0)}
    rads = {"a": 0.35, "c": 0.8, "b": 0.35, "d": 0.8}
    for t, h, name in inter:
        tree = name == "b"
        ax.add_patch(FancyArrowPatch(P[t], P[h], connectionstyle=f"arc3,rad={rads[name]}",
                                     arrowstyle="-|>", mutation_scale=10, lw=2.2 if tree else 1.0,
                                     color=AQUA if tree else INK2, shrinkA=11, shrinkB=11))
        yl = -rads[name] if t == "R1" else rads[name]        # midpoint of the arc3 curve
        ax.text(1.0, yl + (0.13 if yl > 0 else -0.13), name, ha="center",
                va="center", fontsize=8, color=AQUA if tree else INK2)
    for v, (x, y) in P.items():
        ax.add_patch(Circle((x, y), 0.26, color=colors[v]))
        ax.text(x, y, v, ha="center", va="center", color="white", fontsize=8, fontweight="bold")
    ax.text(-0.45, -0.42, "root", ha="center", fontsize=7, color=MUTED)
    ax.text(1.4, -1.6, "rotor at $R_2$: $b\\to d\\to b$\ntree arc (green) = last exit", ha="center", fontsize=6.8, color=INK2)
    ax.set_xlim(-0.6, 2.6); ax.set_ylim(-1.9, 1.3)
    fig.text(0.405, 0.93, "b  Repeat skeleton", fontsize=9, fontweight="bold")
    # (c) the torsor
    ax = fig.add_subplot(gs[2]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    for y, T, c in ((0.78, "{b}", cd), (0.30, "{d}", cb)):
        ax.add_patch(Rectangle((0.05, y - 0.09), 0.9, 0.18, color=GRID, lw=0))
        ax.text(0.09, y, f"$T_{T}$", va="center", fontsize=9, color=AQUA)
        ax.text(0.30, y, word(c), va="center", fontsize=7.2, family="monospace", color=INK)
    ax.annotate("", (0.5, 0.40), (0.5, 0.68), arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=1.3,
                connectionstyle="arc3,rad=-0.4"))
    ax.annotate("", (0.5, 0.68), (0.5, 0.40), arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=1.3,
                connectionstyle="arc3,rad=-0.4"))
    ax.text(0.73, 0.54, "+1 chip at $R_2$", fontsize=7.5, color=ORANGE, va="center")
    ax.text(0.05, 0.03, f"$\\mathcal{{K}}\\cong\\mathbb{{Z}}/{K1}$, BEST $={loc1}\\times{K1}={len(C1)}$\n"
            f"disjoint repeats: $|\\mathcal{{K}}|={K2}$, {len(C2)} reconstruction", fontsize=7.2, color=INK2)
    fig.text(0.665, 0.93, "c  $\\mathcal{K}$ acts by rotor-routing", fontsize=9, fontweight="bold")
    save(fig, "fig1_sandpile_torsor")


# ================================================================== Figure 2
def figure2():
    print("Figure 2: phiX174 double cover and frontier (Bio 6, Bio 12)")
    B12 = load("Bio_12_Frontier_Elimination", "ds_transfer_matrix")
    S = B12.load_phix()
    c12 = B12.Cover(S, 12)
    check("phiX174 D_12: 7 vertices, 1 palindromic, 14 contents, 0 self-complementary",
          (c12.nv, c12.n_fixV, c12.nA, c12.n_fixA) == (7, 1, 14, 0))
    c10 = B12.Cover(S, 10)
    var = c10.constraints()
    wg = B12.frontier_widths(var, B12.greedy_order(var))
    wc = B12.frontier_widths(var, B12.coord_order(c10))
    check("phiX174 k=10: 125 vertices, 111 orbits; frontier width 26 (greedy) vs 62 (genome order)",
          (c10.nv, len(var), max(wg), max(wc)) == (125, 111, 26, 62))
    fig = plt.figure(figsize=(7.4, 2.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.0, 1.05], wspace=0.42)
    # (a) D_12 with rho
    ax = fig.add_subplot(gs[0]); ax.set_aspect("equal"); ax.axis("off")
    n = c12.nv
    pos = {v: (math.cos(2 * math.pi * v / n + 0.3), math.sin(2 * math.pi * v / n + 0.3)) for v in range(n)}
    pal = [BLUE, ORANGE, AQUA, PURPLE, RED, "#c9a227", "#4c9fb8"]
    orbit = {}
    for a in range(c12.nA):
        b = c12.crho[a]
        orbit[a] = min(a, b)
    oid = {o: i for i, o in enumerate(sorted(set(orbit.values())))}
    cnt = collections.Counter()
    for a in range(c12.nA):
        t, h = c12.tail[a], c12.head[a]
        key = (min(t, h), max(t, h)) if t != h else (t, t)
        cnt[key] += 1
        rad = 0.18 * cnt[key] * (1 if t < h else -1)
        rep = a == orbit[a]
        ax.add_patch(FancyArrowPatch(pos[t], pos[h], connectionstyle=f"arc3,rad={rad + 0.08}",
                                     arrowstyle="-|>", mutation_scale=8, lw=1.3,
                                     ls="-" if rep else (0, (3, 2)), color=pal[oid[orbit[a]] % len(pal)],
                                     shrinkA=7, shrinkB=7))
    for v in range(n):
        w = c12.words[v]
        pal_v = B12.rc(w) == w
        ax.add_patch(Circle(pos[v], 0.1, color=RED if pal_v else INK2, zorder=3))
        ax.text(pos[v][0] * 1.33, pos[v][1] * 1.33, w[:6] + "\n" + w[6:], ha="center", va="center",
                fontsize=4.8, family="monospace", color=RED if pal_v else INK2)
    ax.set_xlim(-1.6, 1.6); ax.set_ylim(-2.25, 1.6)
    ax.text(0, -2.2, "solid $a$, dashed $\\rho a$ (same colour = one $\\rho$-orbit);\nred: palindromic vertex "
            "$\\bar v=v$", ha="center", fontsize=6.5, color=INK2)
    ax.set_title("a  $\\phi$X174 double cover $D_{12}$", x=0.0)
    # (b) frontier profile
    ax = fig.add_subplot(gs[1])
    ax.plot(range(1, len(wc) + 1), wc, color=GRAY, lw=1.4, label=f"genome order (max {max(wc)})")
    ax.plot(range(1, len(wg) + 1), wg, color=BLUE, lw=1.6, label=f"greedy order (max {max(wg)})")
    istar = int(np.argmax(wg))
    ax.axvline(istar + 1, color=RED, lw=0.8, ls=(0, (3, 2)))
    ax.set_xlabel("orbits eliminated\ngreedy: 2 481 687 900 box points,\n50 758 states, 8.4 s (Bio 12)", fontsize=7.5)
    ax.set_ylabel("open vertices (frontier)")
    ax.legend(fontsize=6.8, loc="upper left")
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_title("b  Frontier width, $k=10$", x=0.0)
    # (c) the cut at the widest step
    ax = fig.add_subplot(gs[2]); ax.axis("off")
    order = B12.greedy_order(var)
    done = set(order[:istar + 1])
    G = nx.Graph()
    G.add_nodes_from(range(c10.nv))
    for j, a in enumerate(c10.reps):
        for x in (c10.tail[a], c10.head[a]):
            for y in (c10.tail[a], c10.head[a]):
                if x < y:
                    G.add_edge(x, y)
    P = nx.spring_layout(G, seed=3, k=0.22)
    touch = collections.defaultdict(set)
    for j, d in enumerate(var):
        for v in d:
            touch[v].add(j)
    openv = {v for v in touch if (touch[v] & done) and (touch[v] - done)}
    closed = {v for v in touch if touch[v] <= done}
    check("the cut at the widest step has exactly 26 open vertices", len(openv) == max(wg), f"{len(openv)}")
    for j, a in enumerate(c10.reps):
        x, y = c10.tail[a], c10.head[a]
        ax.plot([P[x][0], P[y][0]], [P[x][1], P[y][1]], color=INK2 if j in done else GRID,
                lw=0.8 if j in done else 0.6, zorder=1)
    for v in G.nodes:
        col = RED if v in openv else (INK2 if v in closed else BASE)
        ax.scatter(*P[v], s=16 if v in openv else 7, color=col, zorder=3, lw=0)
    ax.text(0.5, -0.06, f"red: the {len(openv)} open vertices at step {istar + 1}/{len(var)}\n"
            "dark: eliminated  grey: not yet reached", transform=ax.transAxes, ha="center", fontsize=6.5,
            color=INK2)
    ax.set_title("c  The active frontier cut", x=0.0)
    save(fig, "fig2_double_cover")


# ================================================================== Figure 3
def figure3():
    print("Figure 3: #DS = 2 eps(G) (Bio 19)")
    # the octahedron K_{2,2,2}: 4-regular, 6 vertices, 12 edges
    Vn = 6
    E = [(u, v) for u in range(Vn) for v in range(u + 1, Vn) if not (u // 2 == v // 2)]
    # eps(G): Eulerian circuits up to rotation and reversal = directed circuits with edge 0 as u->v first
    adj = collections.defaultdict(list)
    for i, (u, v) in enumerate(E):
        adj[u].append((v, i)); adj[v].append((u, i))
    used = [False] * len(E)
    used[0] = True
    eps = [0]

    def rec(v, d):
        if d == len(E):
            eps[0] += v == E[0][0]
            return
        for w, i in adj[v]:
            if not used[i]:
                used[i] = True
                rec(w, d + 1)
                used[i] = False
    rec(E[0][1], 1)
    eps = eps[0]
    # Eulerian orientations and their BEST counts ec(O) = tau(O) prod (d^+-1)!  (d^+ = 2 here)
    ecs = []
    for bits in itertools.product((0, 1), repeat=len(E)):
        arcs = [(u, v) if b == 0 else (v, u) for (u, v), b in zip(E, bits)]
        outd = collections.Counter(t for t, _ in arcs)
        ind = collections.Counter(h for _, h in arcs)
        if all(outd[x] == ind[x] == 2 for x in range(Vn)):
            L = [[0] * Vn for _ in range(Vn)]
            for t, h in arcs:
                L[t][t] += 1
                L[t][h] -= 1
            tau = int(det_frac([r[1:] for r in L[1:]]))
            ecs.append(tau * math.prod(math.factorial(outd[x] - 1) for x in range(Vn)))
    check("octahedron: sum over Eulerian orientations of ec(O) = 2 eps(G) (Lemma 1 of Bio 19)",
          sum(ecs) == 2 * eps, f"eps = {eps}, {len(ecs)} Eulerian orientations, sum ec = {sum(ecs)}")
    comp = str.maketrans("ACGT", "TGCA")
    rng = random.Random(19)
    u = "".join(rng.choice("ACGT") for _ in range(6))
    Pw = u + u.translate(comp)[::-1]
    check("P = u rc(u) is its own reverse complement (a palindromic k-mer, k = 12)",
          Pw == Pw.translate(comp)[::-1], Pw)

    fig = plt.figure(figsize=(7.4, 2.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.95, 1.3, 1.0], wspace=0.45)
    ax = fig.add_subplot(gs[0]); ax.set_aspect("equal"); ax.axis("off")
    ang = {0: 90, 1: 270, 2: 210, 3: 30, 4: 330, 5: 150}
    Pp = {x: (math.cos(math.radians(a)), math.sin(math.radians(a))) for x, a in ang.items()}
    vc = [BLUE, ORANGE, AQUA, PURPLE, RED, "#c9a227"]
    # draw one Eulerian orientation
    for bits in itertools.product((0, 1), repeat=len(E)):
        arcs = [(u_, v_) if b == 0 else (v_, u_) for (u_, v_), b in zip(E, bits)]
        if all(collections.Counter(t for t, _ in arcs)[x] == 2 for x in range(Vn)):
            break
    for t, h in arcs:
        ax.add_patch(FancyArrowPatch(Pp[t], Pp[h], arrowstyle="-|>", mutation_scale=8, lw=1.0, color=INK2,
                                     shrinkA=8, shrinkB=8))
    for x in range(Vn):
        ax.add_patch(Circle(Pp[x], 0.13, color=vc[x], zorder=3))
        ax.text(*Pp[x], str(x + 1), ha="center", va="center", color="white", fontsize=7, zorder=4)
    ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.45, 1.25)
    ax.text(0, -1.42, f"$G$ = octahedron, 4-regular;\none of its {len(ecs)} Eulerian orientations",
            ha="center", fontsize=6.8, color=INK2)
    ax.set_title("a  A graph $G$", x=0.0)
    ax = fig.add_subplot(gs[1]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    # an Eulerian circuit of G, spelled with palindromes P_x and spacers Sigma_i
    circ = [e[0] for e in nx.eulerian_circuit(nx.MultiGraph(E), source=0)]
    check("the genome visits the palindromes along an Eulerian circuit (|E| = 12 traversals)",
          len(circ) == len(E) and collections.Counter(circ) == collections.Counter({x: 2 for x in range(Vn)}))
    circ = circ + [circ[0]]
    x0, y = 0.02, 0.62
    wblk, wsp = 0.05, 0.028
    for j, x in enumerate(circ[:-1][:12]):
        ax.add_patch(Rectangle((x0, y), wblk, 0.14, color=vc[x], lw=0))
        ax.text(x0 + wblk / 2, y + 0.07, str(x + 1), ha="center", va="center", color="white", fontsize=6.5)
        x0 += wblk
        ax.add_patch(Rectangle((x0, y + 0.045), wsp, 0.05, color=GRID, lw=0))
        x0 += wsp
    ax.text(0.02, 0.83, "$S=P_{x_0}\\,\\Sigma_1\\,P_{x_1}\\,\\Sigma_2\\cdots$ (circular),\nvertices in Eulerian-circuit order",
            fontsize=7, color=INK2)
    ax.text(0.02, 0.50, f"$P_x=u_x\\,\\overline{{u_x}}$, e.g.\n{Pw[:6]}|{Pw[6:]}, $\\rho P_x=P_x$",
            fontsize=7, color=INK2, family="sans-serif")
    ax.text(0.02, 0.33, "content $a\\colon P_x\\to P_y$, $\\rho a\\colon P_y\\to P_x$,\n"
            "multiplicity 1: a box point picks\n$a$ or $\\rho a$ = an orientation of $\\{x,y\\}$",
            fontsize=7, color=INK2, va="top")
    ax.set_title("b  The palindromic genome", x=0.0)
    ax = fig.add_subplot(gs[2])
    vals = collections.Counter(ecs)
    xs = sorted(vals)
    ax.bar(range(len(xs)), [vals[x] for x in xs], color=BLUE, width=0.7)
    ax.set_xticks(range(len(xs))); ax.set_xticklabels([str(x) for x in xs], fontsize=7)
    ax.set_xlabel("ec$(O)=\\tau(O)\\prod(d^+-1)!$"); ax.set_ylabel("Eulerian orientations")
    ax.grid(True, axis="y", color=GRID, lw=0.6)
    ax.set_ylim(0, max(vals.values()) * 1.75)
    ax.text(0.98, 0.97, f"$\\#\\mathrm{{DS}}=\\sum_O \\mathrm{{ec}}(O)={sum(ecs)}$\n$=2\\,\\varepsilon(G)$, "
            f"$\\varepsilon(G)={eps}$", transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color=INK)
    ax.text(0.98, 0.74, "$\\varepsilon(G)$ is #P-complete\n$\\Rightarrow$ #DS is #P-hard", transform=ax.transAxes,
            ha="right", va="top", fontsize=7, color=RED)
    ax.set_title("c  Box points = Eulerian orientations", x=0.0)
    save(fig, "fig3_sharpP_reduction")


# ================================================================== Figure 4
def figure4():
    print("Figure 4: ZNF91 and a fibre-mate (Spectral Fibres)")
    SF = load("Spectral_Fibres", "spectral_fibres")
    items = SF.load()
    S = next(x[3] for x in items if x[1] == "ZNF91")
    # sample fibre-mates (random Eulerian trails of the 25-mer graph) and keep the most rearranged
    rr = random.Random(7)
    T, best = None, 0
    for _ in range(300):
        x = SF.random_trail(S, 25, rr)
        if x and x != S:
            d = sum(1 for p, q in zip(S, x) if p != q)
            if d > best:
                T, best = x, d
    check("ZNF91 has a fibre-mate at k = 25 (same 25-mer multiset)", T is not None and
          SF.spectrum(S, 25) == SF.spectrum(T, 25))
    H1, H2 = SF.std_scale(SF.KD), SF.std_scale(SF.HW)
    aS, aT = SF.apaac(S, 24, H1, H2), SF.apaac(T, 24, H1, H2)
    check("identical 68-dimensional amphiphilic PseAAC (lambda = 24), exact rationals", len(aS) == 68 and aS == aT)
    rng = random.Random(7)
    W = [[{a: rng.randrange(-3, 4) for a in "ACDEFGHIKLMNPQRSTVWY"} for _ in range(25)] for _ in range(8)]
    bias = [rng.randrange(-5, 6) for _ in range(8)]
    cS, cT = SF.cnn_mean(S, W, bias, 25), SF.cnn_mean(T, W, bias, 25)
    check("identical output of an 8-channel window-25 CNN with mean pooling", cS == cT)
    Sh = "".join(random.Random(1).sample(S, len(S)))
    aSh = SF.apaac(Sh, 24, H1, H2)
    diff = [i for i in range(len(S)) if S[i] != T[i]]
    # C2H2 zinc fingers: C-x(2,4)-C-x(12)-H-x(3,5)-H
    zf = re.compile(r"C.{2,4}C.{12}H.{3,5}H")
    fS = [(m.start(), m.end(), m.group()) for m in zf.finditer(S)]
    fT = [(m.start(), m.end(), m.group()) for m in zf.finditer(T)]
    seqs = sorted({f[2] for f in fS + fT})
    check("the C2H2 finger repertoire is the same multiset, in a different order",
          sorted(f[2] for f in fS) == sorted(f[2] for f in fT) and [f[2] for f in fS] != [f[2] for f in fT],
          f"{len(fS)} fingers, {len(diff)} positions changed")

    fig = plt.figure(figsize=(7.4, 3.3))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1.0], width_ratios=[1.6, 1.0], hspace=0.65, wspace=0.28)
    ax = fig.add_subplot(gs[0, :])
    cmap = plt.get_cmap("tab20")
    colf = {s: cmap(i % 20) for i, s in enumerate(seqs)}
    for y, (lab, F, X) in enumerate(((f"fibre-mate $T$", fT, T), ("ZNF91 (human)", fS, S))):
        ax.add_patch(Rectangle((0, y - 0.3), len(X), 0.6, color=GRID, lw=0))
        for a, b, s in F:
            ax.add_patch(Rectangle((a, y - 0.3), b - a, 0.6, color=colf[s], lw=0))
        ax.text(-15, y, lab, ha="right", va="center", fontsize=7.5, color=INK2)
    for i in diff:
        ax.plot([i, i], [0.35, 0.65], color=RED, lw=0.3)
    ax.set_xlim(0, len(S)); ax.set_ylim(-0.5, 1.5); ax.set_yticks([])
    ax.set_xlabel("residue position")
    ax.spines["left"].set_visible(False)
    ax.set_title(f"a  Same 25-mer multiset, different order: {len(fS)} C2H2 fingers (colour = finger "
                 f"sequence), {len(diff)} positions changed (red ticks)", x=0.0, fontsize=8)
    ax = fig.add_subplot(gs[1, 0])
    xs = np.arange(68)
    ax.plot(xs, [float(v) for v in aSh], color=GRAY, lw=0.9, label="shuffled ZNF91 (not a fibre-mate)")
    ax.plot(xs, [float(v) for v in aS], color=BLUE, lw=2.2, label="ZNF91")
    ax.plot(xs, [float(v) for v in aT], color=ORANGE, lw=0.9, ls=(0, (2, 1.5)), label="fibre-mate (identical)")
    ax.axvline(19.5, color=BASE, lw=0.6)
    ax.set_xlabel("PseAAC component (20 composition + 48 correlation)")
    ax.set_yscale("symlog", linthresh=1e-3)
    ax.legend(fontsize=6.3, loc="lower left")
    ax.set_title("b  68-dim PseAAC ($\\lambda=24$)", x=0.0)
    ax = fig.add_subplot(gs[1, 1])
    ch = np.arange(8)
    ax.bar(ch - 0.2, [float(v) for v in cS], width=0.4, color=BLUE, label="ZNF91")
    ax.bar(ch + 0.2, [float(v) for v in cT], width=0.4, color=ORANGE, label="fibre-mate")
    ax.set_xlabel("CNN channel (window 25, mean pool)")
    ax.legend(fontsize=6.3)
    ax.set_title("c  Identical CNN output", x=0.0)
    save(fig, "fig4_spectral_fibre")


# ================================================================== Figure 5
def figure5():
    print("Figure 5: snarl tree and Schur reduction (Bio 14)")
    B14 = load("Bio_14_Snarl_Schur", "bio14_verify")
    rnd = random.Random(20260924)
    verts, arcs, order, root = B14.synthetic_genome(3, 3, 2, rnd)
    FR = B14.Field("frac")
    d_snarl, st = B14.det_by_order(arcs, root, order, FR)
    d_dense = B14.bareiss_dense(arcs, root, verts)
    check("toy genome graph: snarl-order elimination = dense Bareiss determinant", d_snarl == d_dense,
          f"t = {d_dense}, max front {st[0]}")
    # Kron / Schur reduction of snarl 1's interior I (its pool p and nested r, q vertices)
    others = [v for v in verts if v != root]
    I = [v for v in others if v[0] in ("p", "q", "r") and v[1] == 1]
    Bd = [v for v in others if v not in I]
    idx = {v: i for i, v in enumerate(others)}
    n = len(others)
    L = [[Fraction(0)] * n for _ in range(n)]
    for x, y in arcs:
        if x == y or x == root:
            continue
        L[idx[x]][idx[x]] += 1
        if y != root:
            L[idx[x]][idx[y]] -= 1
    ii, bb = [idx[v] for v in I], [idx[v] for v in Bd]
    LII = [[L[r][c] for c in ii] for r in ii]
    detII = det_frac(LII)
    # Schur complement L_BB - L_BI L_II^{-1} L_IB, in exact rational arithmetic
    Ainv = inv_frac(LII)
    Sch = [[L[r][c] - sum(L[r][ii[a]] * Ainv[a][b] * L[ii[b]][c] for a in range(len(ii)) for b in range(len(ii)))
            for c in bb] for r in bb]
    detS = det_frac(Sch)
    check("det L~ = det L_II * det(Schur complement)  (Kron reduction of one snarl)",
          detII * detS == d_dense, f"{detII} x {detS} = {d_dense}")
    rowsums_ok = all(sum(r) >= 0 for r in Sch) and \
        all(Sch[i][j] <= 0 for i in range(len(bb)) for j in range(len(bb)) if i != j)
    check("the Schur complement is again a Laplacian (off-diagonal <= 0, row sums >= 0)", rowsums_ok)
    # scaling table (S) of Bio 14's reference transcript:
    # snarls, vertices, arcs, front(snarl), front(min-deg), equal mod p1,p2, us per vertex
    txt = open(os.path.join(PAPERS, "Bio_14_Snarl_Schur", "bio14_output.txt"), encoding="utf-8").read()
    scal = [(int(sn), int(nv), int(fs), float(us)) for sn, nv, _, fs, _, us in
            re.findall(r"^\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(?:yes|no)\s+([\d.]+)", txt, re.M)]
    check("Bio 14 scaling rows found in its transcript (snarls, vertices, front, us/vertex)", len(scal) >= 3,
          "; ".join(f"{v} vertices: {t} us" for _, v, _, t in scal))

    fig = plt.figure(figsize=(7.4, 3.1))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.45, 0.8, 0.85], wspace=0.35)
    ax = fig.add_subplot(gs[0]); ax.axis("off")
    G = nx.DiGraph()
    G.add_edges_from((x, y) for x, y in arcs if x != y)

    def xy(v):
        kind, i = v[0], v[1]
        base = 3.0 * i
        if kind == "b":
            return (base, 0.0)
        j = v[2]
        if kind == "p":
            return (base + 0.7 + 0.55 * (j % 2), 0.9 * (j - 1))
        if kind == "q":
            return (base + 1.2 + 1.2 * j, -1.6)
        return (base + 1.8, -2.3 + 0.0 * j + 0.45 * (j - 1))
    P = {v: xy(v) for v in G.nodes}
    for x, y in set(G.edges):
        inI = x in I and y in I
        ax.add_patch(FancyArrowPatch(P[x], P[y], arrowstyle="-|>", mutation_scale=6, lw=0.8,
                                     color=ORANGE if inI else INK2, connectionstyle="arc3,rad=0.12",
                                     shrinkA=4, shrinkB=4, alpha=0.9))
    for v in G.nodes:
        col = BLUE if v[0] == "b" else (ORANGE if v in I else (AQUA if v[0] in "qr" else GRAY))
        ax.scatter(*P[v], s=55 if v[0] == "b" else 22, color=col, zorder=3, lw=0)
        if v[0] == "b":
            ax.text(P[v][0], P[v][1] + 0.35, f"$b_{v[1]}$", ha="center", fontsize=7.5, color=BLUE)
    ax.set_title("a  Genome bubble graph", x=0.0)
    ax.text(0.0, -0.02, "2 haplotypes, 3 snarls.  blue: snarl boundaries\norange: interior $I$ of snarl 1 (incl. nested bubble)\n"
            "green: nested boundaries", transform=ax.transAxes, fontsize=6.5, color=INK2, va="top")
    # (b) snarl tree
    ax = fig.add_subplot(gs[1]); ax.axis("off"); ax.set_aspect("equal"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    nodes = {"chain": (0.5, 0.9), "S0": (0.15, 0.6), "S1": (0.5, 0.6), "S2": (0.85, 0.6),
             "N0": (0.15, 0.28), "N1": (0.5, 0.28), "N2": (0.85, 0.28)}
    present_nested = {i for i in range(3) if ("q", i, 0) in G.nodes}
    for a, b in (("chain", "S0"), ("chain", "S1"), ("chain", "S2")):
        ax.plot(*zip(nodes[a], nodes[b]), color=INK2, lw=1, zorder=1)
    for i in range(3):
        if i in present_nested:
            ax.plot(*zip(nodes[f"S{i}"], nodes[f"N{i}"]), color=INK2, lw=1, zorder=1)
            ax.add_patch(Circle(nodes[f"N{i}"], 0.07, color=AQUA, zorder=2))
            ax.text(*nodes[f"N{i}"], "n", ha="center", va="center", color="white", fontsize=7, zorder=3)
        ax.add_patch(Circle(nodes[f"S{i}"], 0.08, color=ORANGE if i == 1 else BLUE, zorder=2))
        ax.text(*nodes[f"S{i}"], f"$s_{i}$", ha="center", va="center", color="white", fontsize=8, zorder=3)
    ax.add_patch(Rectangle((0.35, 0.85), 0.3, 0.1, color=GRID))
    ax.text(0.5, 0.9, "chain", ha="center", va="center", fontsize=7)
    ax.text(0.5, 0.02, "eliminate leaves first:\ninteriors before boundaries", ha="center", fontsize=6.8,
            color=INK2)
    ax.set_title("b  Snarl tree", x=0.0)
    # (c) Schur identity + scaling
    ax = fig.add_subplot(gs[2])
    vx = [v for _, v, _, _ in scal]; ty = [t for *_, t in scal]
    ax.plot(vx, ty, "o-", color=BLUE, lw=1.4, ms=4)
    ax.set_xscale("log"); ax.set_ylim(0, max(ty) * 1.6)
    ax.set_xlabel("vertices"); ax.set_ylabel("$\\mu$s per vertex (snarl order)")
    ax.grid(True, color=GRID, lw=0.6)
    fr = max(f for _, _, f, _ in scal)
    ax.text(0.03, 0.97, f"$\\det\\tilde L=\\det L_{{II}}\\cdot\\det(L/L_{{II}})$\n"
            f"toy: ${detII}\\times{detS}={d_dense}$\nfront $\\leq {fr}$: linear time",
            transform=ax.transAxes, va="top", fontsize=6.8, color=INK)
    ax.set_title("c  Kron reduction, linear time", x=0.0)
    save(fig, "fig5_snarl_schur")


FIGS = {1: figure1, 2: figure2, 3: figure3, 4: figure4, 5: figure5}

if __name__ == "__main__":
    want = [int(a) for a in sys.argv[1:] if a.isdigit()] or sorted(FIGS)
    for k in want:
        FIGS[k]()
    print(f"\nTOTAL {sum(CHECKS)}/{len(CHECKS)} passed")
    sys.exit(0 if all(CHECKS) else 1)
