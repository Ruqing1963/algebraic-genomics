# data/

| file | content |
|---|---|
| `phiX174.fasta` | reference genome NC_001422.1 (5386 bp, circular), the running example of Parts I–II |
| `tables/key_results.csv` | the headline numbers of the monograph, each with the script that verifies it |

**Where the other data live.** Each chapter's input data sits next to its verification script, because
the scripts read them from their own folder. The large files are:

| file | size | used by |
|---|---|---|
| `papers/Spectral_Fibres/human_reviewed.tsv.gz` | 7.0 MB | UniProtKB/Swiss-Prot human reviewed proteome (20 431 proteins) |
| `papers/Spectral_Fibres/fibre_results.json` | 7.2 MB | per-protein fibre sizes and κ (output of `spectral_fibres.py`) |
| `papers/Bio_05_Ecoli_Decomposition/NC_000913.3.fasta`, `.ft.txt` | 6.2 MB | *E. coli* K-12 MG1655 genome and features |
| `papers/Bio_10_Pangenome_Sheaf/pangenome/*.fasta` | 36 MB | seven *E. coli* strains for the real pan-genome run |
| `papers/Bio_10_Pangenome_Sheaf/phix174_isolates.fasta` | 0.2 MB | φX174 isolates |
| `papers/*/phix174_NC_001422.1.txt` | 5 kB each | φX174 reference, one copy per chapter that reads it |

Licensing of third-party data: see `../LICENSE-CC-BY-4.0.md`.
