# Ancillary files for *The double-stranded lattice sum as a fermionic partition function*

Ruqing Chen, 22 September 2026. Bio 13 of the applied strand; numbered XXVI in the series.

## Contents

| file | what it is |
|---|---|
| `anccft26.tex` / `.pdf` | the note |
| `ds_fermion.py` | the Berezin identity, the fixed root, the halving, the fermionic sweep and the arborescence elimination; 7 checks. Standard library only |
| `ds_count.py` | the compiled enumeration of the k = 10 half-box; 4 checks and the number. Needs `numpy` and `numba`, imports `ds_fermion.py` |
| `ds_fermion_output.txt` | transcript of `python -u ds_fermion.py --quick` |
| `ds_count_output.txt` | transcript of `python -u ds_count.py` |
| `phix174_NC_001422.1.txt` | bacteriophage phiX174, RefSeq NC_001422.1, 5386 bp |

## Running it

```
python -u ds_fermion.py --quick     # about 10 minutes
python -u ds_count.py               # about 40 minutes on 16 threads; resumes from
                                    # ds_count_k10_checkpoint.json if interrupted
```

## The result

```
#DS(phiX174, 10) = 31,925,753,246,212   (about 3.19e13)
```

This fills the last blank of Bio 8's Table 2:

| k | rank Lambda^- | box points | connected | #DS |
|---|---|---|---|---|
| 12 | 4 | 10 | 8 | 10 |
| 11 | 12 | 528 | 218 | 864 |
| 10 | 48 | 2,481,687,900 | 909,103,210 | **31,925,753,246,212** |

It lies inside Bio 8's bounds 8,610,708,632 <= #DS <= 1.368e36. On average, a point with
connected support carries about 35,118 circular sequences.

## The checks

| key | file | claim |
|---|---|---|
| D | ds_fermion | the double cover reproduces Bio 8's Table 1 at k = 10, 11, 12 |
| F | ds_fermion | the top coefficient of prod_a (1 - c psibar_t (psi_t - psi_h)) prod_{v not in supp}(1 - psibar_v psi_v) is (-1)^{\|V\|-1} times the reduced determinant, on 93 points |
| R | ds_fermion | on the half-box f_{j0} > 0 the fixed root always lies in supp c |
| S | ds_fermion | Nseq(m-c) = Nseq(c), and the half-box carries half of #DS |
| V | ds_fermion | the fermionic sweep gives 10, 864 and the brute force on 12 random covers |
| A | ds_fermion | the arborescence elimination likewise, also as a sum over sub-boxes; forced-arc contraction agrees with the full determinant on 636 points |
| Q | ds_fermion | on 5 random k = 10 sub-boxes, the arborescence elimination equals the point-by-point sum |
| W | ds_count | the compiled enumeration gives 10, 864, the point and connected counts, and the brute force on 12 random covers |
| X | ds_count | on k = 10 sub-boxes it agrees with the arborescence elimination (an unrelated algorithm) |
| P | ds_count | the k = 10 half-box has exactly 1,240,843,950 points and 454,551,605 connected ones: half of Bio 12's counts |
| K | ds_count | #DS(phiX174, 10) = 31,925,753,246,212 |

All checks are exact. The enumeration works modulo the five largest primes below 2^31, and the
Chinese remainder theorem gives the exact integer because it lies below Bio 8's upper bound.

## What did not work, and why

* **The fermionic frontier.** Its state is linear: a vector in the exterior algebra of the open
  generators. At k = 10 the order has 35 open generators at its widest, and a single history
  already fills up to 3^s coordinates, where s is the number of open rows carrying an arc. The
  signed terms cancel only at the end. The sweep is exact at k = 11, 12, but it does not finish
  a sub-box at k = 10.
* **The arborescence frontier.** This is the positive form of the same weight, with one state
  per history. It is exact, but its states (candidate multisets and class pointers) are finer
  than Bio 12's connectivity partition, which already needed 9.7 million states. Sub-boxes of
  12-166 points already peak at 55,890-623,000 states.
* **What worked.** The least clever method, made cheap at each point: forced arcs are
  contracted first (Lemma 7), so each point needs only a determinant of small order.

## Traps recorded

* **The sign of the interleaved order** (Remark 2). Listing all psibar before all psi makes an
  implementation of the Berezin product return -#DS on some instances.
* **The capacity check has to cover every endpoint of the step,** not only the vertices whose
  imbalance changed. Their remaining capacity has shrunk even when their imbalance has not;
  skipping them costs pruning, not correctness.
* **A checkpoint is keyed by its signature.** The k = 10 run stores `k`, the root orbit and the
  primes, and it discards a checkpoint written with anything else.
