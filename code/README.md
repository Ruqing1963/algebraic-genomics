# code/

`generate_figures.py` draws the five overview figures of the monograph into `../figures/`. Every
number in a figure is computed at run time from the chapters' own modules and data, and checked; the
script ends with `TOTAL a/b passed`.

```bash
python code/generate_figures.py        # all five, ~1-2 min
python code/generate_figures.py 2 4    # only figures 2 and 4
```

## Chapter verification scripts

The verification scripts stay in their chapter folders, because they read data from beside themselves
and import each other through relative paths. Run them from there, or all at once with
`python build_all.py`.

| chapter | script(s) |
|---|---|
| Bio 1 | `papers/Bio_01_Critical_Groups/dbg_sandpile.py` |
| Bio 2 | `papers/Bio_02_Rotor_Routing/rotor_assembly.py` |
| Bio 3 | `papers/Bio_03_Repeat_Splitting/repeat_splitting.py` |
| Bio 4 | `papers/Bio_04_Multicopy_Repeats/multicopy.py` |
| Bio 5 | `papers/Bio_05_Ecoli_Decomposition/ecoli_sandpile.py` |
| Bio 6 | `papers/Bio_06_RC_Double_Cover/rc_double_cover.py` |
| Bio 7 | `papers/Bio_07_Signed_Quotient/bidirected_sandpile.py` |
| Bio 8 | `papers/Bio_08_DS_Lattice_Sum/ds_lattice.py`, `ds_flips.py` |
| Bio 9 | `papers/Bio_09_Polyploid_Phasing/polyploid_phasing.py` |
| Bio 10 | `papers/Bio_10_Pangenome_Sheaf/pangenome_sheaf.py`, `real_pangenome.py` |
| Bio 11 | `papers/Bio_11_Isoform_Polytope/isoform_polytope.py` |
| Bio 12 | `papers/Bio_12_Frontier_Elimination/ds_transfer_matrix.py` |
| Bio 13 | `papers/Bio_13_Fermionic_Partition/ds_fermion.py`, `ds_count.py` |
| Bio 14 | `papers/Bio_14_Snarl_Schur/bio14_verify.py` |
| Bio 17 | `papers/Bio_17_Inversions_Fibre/bio17_verify.py` |
| Bio 19 | `papers/Bio_19_SharpP_Hardness/bio19_verify.py` |
| Spectral Fibres | `papers/Spectral_Fibres/spectral_fibres.py`, `make_figures.py` |

Shared support module: `lib/pgl3_building.py`.
