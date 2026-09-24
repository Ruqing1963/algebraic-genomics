# Running the pan-genome computation on your own machine

`real_pangenome.py` is **standalone**: Python 3.8+ and `numpy`, nothing else. No `sympy`, and
no other file of the series — the surface pools of check (V) are built inside the script, and
the chromosomes are fetched from NCBI into `./pangenome/` on first use.

## What to copy

Just the one file. Copy `pangenome/` too if you would rather not re-download 34 MB from NCBI;
`NC_000913.3.fasta` (K-12 MG1655) lives in `Bio 5`, and the script will fetch it if it is not
beside the others.

## The commands

```bash
# everything: checks (A) and (V), then all eight chromosomes.  About three minutes.
python -u real_pangenome.py --out real_pangenome_output.txt

# one pool, no checks
python -u real_pangenome.py --strains 5 -k 31 --skip-checks --out real_pangenome_p5.txt

# a different word length
python -u real_pangenome.py --strains 8 -k 51 --skip-checks --out real_pangenome_k51.txt
```

| flag | what it does |
|---|---|
| `--strains N` | how many of the eight chromosomes to pool, in the order listed in `STRAINS` (default all 8) |
| `-k 31` | word length; comma-separated for several in one run (`-k 31,51`) |
| `--primes 2,3` | primes tested **only if** the unit-pivot certificate fails to complete; on the data here it always completes, so no prime is ever tested |
| `--cache` | store the reduced boundary matrices under `./cache/`; a re-run with further `--primes` then skips the reduction |
| `--out FILE` | tee the transcript to a file as well as the screen |
| `--skip-checks` | go straight to the chromosomes |
| `--quiet` | no progress lines during elimination |

Progress lines appear every 20,000 columns during a long elimination, with a rate, so you can
tell within a minute whether a run is going to finish.

## Times and memory

One core of a laptop, k = 31, cold (no cache). Every pool is certified free.

| strains | branch vtx | nerve `q=1 / q=2` | reduction | total | `Rec` |
|---|---|---|---|---|---|
| 2 | 1,509 | 86 / 0 | — | 26 s | `Z^85` |
| 3 | 63,245 | 58,454 / 28,820 | 0 s | 40 s | `Z^29632` |
| 4 | 137,375 | 166,011 / 134,401 | 4 s | 60 s | `Z^66447` |
| 5 | 180,256 | 260,369 / 277,947 | 9 s | 84 s | `Z^87858` |
| 6 | 214,033 | 409,007 / 588,230 | 8 s | 110 s | `Z^104495` |
| 7 | 237,525 | 575,831 / 1,053,725 | 16 s | 140 s | `Z^116056` |
| 8 | 241,592 | 744,598 / 1,683,865 | 28 s | 182 s | `Z^117848` |

Peak memory is about 6 GB at eight strains, nearly all of it the `q=2` boundary. Most of the
wall clock is now the branch-vertex pass over 40 Mbp, not the algebra.

Going beyond eight means supplying more accessions in `STRAINS`. The nerve's two-cell count
grows as `C(p,3)` times the fragmentation of a triple intersection, so it is the memory rather
than the time that will stop you: at eight strains `C_2` is already 1.7 million.

## Two fixes that made this size reachable

**Pivot ordering.** `unit_reduce` and `sparse_rank_mod` used to choose each pivot by
`sorted(alive, key=len)` — an `O(n log n)` scan of the whole live column set **per pivot**. At
278,000 columns that scan was essentially the entire cost of the script: the 5-strain reduction
took 14,080 s. Both now take the shortest column from a lazily-updated heap, so the total cost
is proportional to the fill-in. The same matrix now reduces in 9 s.

The better ordering also produces less fill-in, and that changed an answer. At 4 and 5 strains
the old ordering inflated entries past `±1` and stalled with a residual block, so the note
could only report "no 2- or 3-torsion". The new ordering eliminates every column, which
**certifies `R` free at every prime**. The rational ranks — 66,447 and 87,858, and the
`rank d_2` values 99,561 and 172,507 — are identical under both orderings, which is the
cross-check that the two runs compute the same map.

**Component representatives.** `boundary_cols` found each component's representative with
`next(v for v, cc in comp[s].items() if cc == c)`, a linear scan of the intersection for every
one of its components — quadratic, and at six strains it alone outweighed the whole reduction.
One pass now builds the representative table. The 3-strain nerve step went 95 s → 1 s.
