# -*- coding: utf-8 -*-
"""Reproduce the 12 fibre-mates of check (M) of spectral_fibres.py (same seed, same order) and
compute the Smith-Waterman ratio for all of them, including L > 1500."""
import json, os, random, time
import spectral_fibres as SF

HERE = os.path.dirname(os.path.abspath(__file__))
items = SF.load()
res = json.load(open(os.path.join(HERE, "fibre_results.json")))
for r in res:
    r["logN"] = {int(k): v for k, v in r["logN"].items()}
byacc = {x[0]: x for x in items}
rr = random.Random(7)
# the same draws as in check (M): CNN weights and biases first
Wc = [[{a: rr.randrange(-3, 4) for a in "ACDEFGHIKLMNPQRSTVWY"} for _ in range(25)] for _ in range(8)]
bc = [rr.randrange(-5, 6) for _ in range(8)]
cand = sorted([r for r in res if (r["logN"].get(25) or 0) > 0 and r["L"] <= 2500], key=lambda r: r["L"])
for r in cand[:: max(1, len(cand) // 12)][:12]:
    S = byacc[r["acc"]][3]
    T = None
    for _ in range(200):
        x = SF.random_trail(S, 25, rr)
        if x and x != S:
            T = x
            break
    if T is None:
        continue
    t = time.time()
    diff = sum(1 for a, b in zip(S, T) if a != b)
    sw = SF.sw_score(S, T) / SF.sw_score(S, S)
    print(f"{r['gene'] or r['acc']:10s} L={len(S):5d} diff={diff:4d} SW={sw:.3f}  [{time.time()-t:.0f}s]",
          flush=True)
