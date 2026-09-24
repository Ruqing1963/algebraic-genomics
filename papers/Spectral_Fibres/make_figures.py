# -*- coding: utf-8 -*-
"""Figures 2-4 of spectral_fibres.tex, from fibre_results.json and the fibre-mates of check (M).
Needs matplotlib (Anaconda).  Writes fig_resolution.pdf, fig_families.pdf, fig_mates.pdf."""
import json, os, random, collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import spectral_fibres as SF

HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, MUTED, GRID, BASE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, AQUA, GRAY = "#2a78d6", "#eb6834", "#1baf7a", "#b9b8b1"
plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 8.5, "axes.edgecolor": BASE, "axes.labelcolor": INK2,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.8, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.axisbelow": True, "legend.frameon": False, "pdf.fonttype": 42, "axes.titlesize": 9,
    "axes.titleweight": "bold", "axes.titlelocation": "left", "axes.titlecolor": INK})

res = json.load(open(os.path.join(HERE, "fibre_results.json")))
for r in res:
    r["logN"] = {int(k): v for k, v in r["logN"].items()}
n = len(res)


def group(r):
    f, g = r["fam"], (r["gene"] or "")
    if "C2H2-type zinc-finger" in f:
        return "KZFP"
    if any(x in f for x in ("NBPF", "NPIP", "GOLGA6", "POTE")):
        return "segdup"
    return "other"


# ---------------------------------------------------------------- Figure 2: resolution
fig, (a, b) = plt.subplots(1, 2, figsize=(6.6, 2.5), gridspec_kw=dict(wspace=0.32))
ks = SF.KREP
share = [100 * sum(1 for r in res if r["L"] >= k and (r["logN"].get(k) or 0) > 1e-9)
         / sum(1 for r in res if r["L"] >= k) for k in ks]
a.plot(ks, share, color=BLUE, lw=1.6, marker="o", ms=4, mec="white", mew=0.8, zorder=3)
a.set_xscale("log")
a.set_yscale("log")
a.set_xticks([1, 2, 3, 5, 10, 25, 50])
a.set_xticklabels(["1", "2", "3", "5", "10", "25", "50"])
a.set_yticks([0.1, 1, 10, 100])
a.set_yticklabels(["0.1%", "1%", "10%", "100%"])
a.set_xlabel("word length $k$")
a.set_ylabel("proteins with $N_k>1$")
a.set_title("a  Unresolved at word length $k$")
for k, s in zip(ks, share):
    if k in (3, 5, 10, 25):
        a.annotate(f"{s:.1f}%", (k, s), xytext=(5, 4), textcoords="offset points", color=INK2,
                   fontsize=7.5)
kap = sorted(r["kappa"] for r in res if r["kappa"] is not None)
ys = [100 * (i + 1) / len(kap) for i in range(len(kap))]
b.step(kap, ys, where="post", color=BLUE, lw=1.6)
b.set_xscale("log")
b.set_xlim(1.8, 1100)
b.set_ylim(0, 101)
b.set_xticks([2, 5, 10, 20, 50, 100, 1000])
b.set_xticklabels(["2", "5", "10", "20", "50", "100", "1000"])
b.set_xlabel(r"spectral resolution $\kappa$")
b.set_ylabel("cumulative share of proteins")
b.set_yticks([0, 25, 50, 75, 100])
b.set_yticklabels(["0%", "25%", "50%", "75%", "100%"])
b.set_title(r"b  Spectral resolution $\kappa$")
b.axvline(5, color=BASE, lw=0.8, ls=(0, (3, 3)), zorder=1)
b.annotate("median 5", (5, 53.5), xytext=(6, -14), textcoords="offset points", color=INK2, fontsize=7.5)
b.annotate("NBPF20\n$\\kappa=951$", (951, 100), xytext=(-40, -34), textcoords="offset points",
           color=INK2, fontsize=7.5, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))
fig.savefig(os.path.join(HERE, "fig_resolution.pdf"), bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- Figure 3: families
fam = collections.defaultdict(list)
for r in res:
    if r["fam"] and r["kappa"]:
        fam[r["fam"]].append(r["kappa"])
short = {"Krueppel C2H2-type zinc-finger protein family": "Krueppel C2H2 zinc finger",
         "Fibrillar collagen family": "fibrillar collagen",
         "Fibril-associated collagens with interrupted helices (FACIT) family": "FACIT collagen",
         "KRTAP type 4 family": "KRTAP type 4", "KRTAP type 5 family": "KRTAP type 5",
         "KRTAP type 10 family": "KRTAP type 10", "NBPF family": "NBPF", "GOLGA6 family": "GOLGA6",
         "NPIP family": "NPIP", "POTE family": "POTE", "KRTAP type 28 family": "KRTAP type 28", "KRTAP type 9 family": "KRTAP type 9"}
rows = [(sorted(v)[len(v) // 2], f) for f, v in fam.items() if len(v) >= 8]
rows = sorted(rows, reverse=True)[:10]
fig, (a, b) = plt.subplots(1, 2, figsize=(6.8, 3.0), gridspec_kw=dict(wspace=0.55, width_ratios=[1.15, 1]))
labels = []
rnd = random.Random(1)
for y, (med, f) in enumerate(reversed(rows)):
    v = fam[f]
    col = ORANGE if short.get(f, "") in ("NBPF", "GOLGA6", "NPIP", "POTE") else \
        (BLUE if "zinc" in f else AQUA)
    a.scatter(v, [y + rnd.uniform(-0.18, 0.18) for _ in v], s=9, color=col, alpha=0.55, lw=0, zorder=2)
    a.plot([med, med], [y - 0.32, y + 0.32], color=INK, lw=1.6, zorder=3)
    labels.append(f"{short.get(f, f)}  (n={len(v)})")
allmed = kap[len(kap) // 2]
a.axvline(allmed, color=BASE, lw=0.8, ls=(0, (3, 3)), zorder=1)
a.text(allmed * 1.05, -0.75, "all proteins: 5", color=INK2, fontsize=7)
a.set_yticks(range(len(rows)))
a.set_yticklabels(labels, fontsize=7.5, color=INK2)
a.set_xscale("log")
a.set_xticks([3, 10, 30, 100, 300, 1000])
a.set_xticklabels(["3", "10", "30", "100", "300", "1000"])
a.set_xlabel(r"spectral resolution $\kappa$ (bar: family median)")
a.set_title("a  Families with the largest median $\\kappa$")
a.grid(axis="y", visible=False)
# b: kappa against length, three groups
for gname, col, mk, lab in (("other", GRAY, "o", "other proteins"),
                            ("KZFP", BLUE, "s", "C2H2 zinc-finger proteins"),
                            ("segdup", ORANGE, "^", "NBPF, NPIP, GOLGA6, POTE")):
    pts = [(r["L"], r["kappa"]) for r in res if r["kappa"] and group(r) == gname]
    b.scatter([p[0] for p in pts], [p[1] for p in pts], s=6 if gname == "other" else 11, marker=mk,
              color=col, alpha=0.35 if gname == "other" else 0.8, lw=0,
              label=f"{lab} ({len(pts):,})", zorder=2 if gname == "other" else 3)
b.set_xscale("log")
b.set_yscale("log")
b.set_xlabel("protein length $L$")
b.set_ylabel(r"$\kappa$")
b.set_title("b  Resolution is not a length effect")
b.legend(loc="upper center", bbox_to_anchor=(0.45, -0.2), fontsize=7, handletextpad=0.3, markerscale=1.4)
b.set_yticks([2, 5, 10, 30, 100, 1000])
b.set_yticklabels(["2", "5", "10", "30", "100", "1000"])
fig.savefig(os.path.join(HERE, "fig_families.pdf"), bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- Figure 4: fibre-mates
items = SF.load()
byacc = {x[0]: x for x in items}
rr = random.Random(7)
Wc = [[{a_: rr.randrange(-3, 4) for a_ in "ACDEFGHIKLMNPQRSTVWY"} for _ in range(25)] for _ in range(8)]
bc = [rr.randrange(-5, 6) for _ in range(8)]
cand = sorted([r for r in res if (r["logN"].get(25) or 0) > 0 and r["L"] <= 2500], key=lambda r: r["L"])
mates = []
for r in cand[:: max(1, len(cand) // 12)][:12]:
    S = byacc[r["acc"]][3]
    T = None
    for _ in range(200):
        x = SF.random_trail(S, 25, rr)
        if x and x != S:
            T = x
            break
    if T is not None:
        mates.append((r["gene"] or r["acc"], S, T))
mates = sorted(mates, key=lambda m: -sum(1 for p, q in zip(m[1], m[2]) if p != q))
fig, ax = plt.subplots(figsize=(6.6, 3.1))
for y, (g, S, T) in enumerate(reversed(mates)):
    L = len(S)
    ax.plot([0, L], [y, y], color=GRID, lw=5.5, solid_capstyle="round", zorder=1)
    # contiguous blocks of changed positions
    i = 0
    while i < L:
        if S[i] != T[i]:
            j = i
            while j < L and S[j] != T[j]:
                j += 1
            ax.plot([i, max(j, i + 3)], [y, y], color=ORANGE, lw=5.5, solid_capstyle="butt", zorder=2)
            i = j
        else:
            i += 1
    d = sum(1 for p, q in zip(S, T) if p != q)
    ax.text(L + 25, y, f"{d} changed ({100*d/L:.1f}%)", va="center", fontsize=7, color=INK2)
ax.set_yticks(range(len(mates)))
ax.set_yticklabels([m[0] for m in reversed(mates)], fontsize=7.5, color=INK2)
ax.set_xlabel("residue position")
ax.set_xlim(0, 2500)
ax.grid(axis="y", visible=False)
ax.set_title("Positions at which a fibre-mate differs from the human protein (same 25-spectrum, "
             "identical PseAAC and CNN output)", fontsize=8.5)
fig.savefig(os.path.join(HERE, "fig_mates.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(HERE, "fig_mates.png"), dpi=130, bbox_inches="tight")
plt.close(fig)
print("figures written")
