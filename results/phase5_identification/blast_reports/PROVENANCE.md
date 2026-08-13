# Raw search output backing the Phase 5 report

Plain-text responses from the NCBI BLAST URL API, saved as returned. Reports
longer than 260 lines are truncated after the description table and the first
alignments, with a marker at the cut. Nothing else is edited.

Searches were run on 2026-08-13 against `core_nt` (128,547,243 sequences,
948,086,538,397 letters, posted 2026-07-18) and `nr`.

Whole test set, `blastn` vs `core_nt`, all 38 contigs in one batch:

- `all38_vnom_contigs_nt.txt` — 24 Peach latent mosaic viroid, 14 *Prunus*
  host, 0 unidentified

Phase 1 cluster centroids, run individually first, `blastn` vs `core_nt`:

- `NODE_36652_nt.txt` — Peach latent mosaic viroid, 98% over 333 nt, E=5e-155
- `NODE_43211_nt.txt` — *Prunus* genomic
- `NODE_52236_nt.txt` — *Prunus* genomic / 40S ribosomal protein mRNA
- `NODE_38543_nt.txt` — *Prunus* genomic

Phase 3 candidates:

- `k25_megablast_nt.txt` — megablast, no significant similarity
- `k25_blastn_nt.txt` — blastn, diverged hits to Alphaproteobacteria MAGs
- `k25_blastx_nr.txt` — FAD-dependent oxidoreductase, 97% aa
- `k141_3_blastn_nt.txt`, `k141_3_blastx_nr.txt` — MBL fold metallo-hydrolase, 73% aa
- `k141_10_blastn_nt.txt`, `k141_10_blastx_nr.txt` — DUF2007 protein, 88% aa

Library-type evidence for SRR13291825:

- `primers.txt` — the three dominant 5' 20-mers, blastn with word size 7
- `repreads.txt` — one full read per primer class, all 18S rRNA

Rfam Infernal cmscan results are in `../rfam/` (zero hits for every query).

Regenerate any of this with:

```
python scripts/phase5_identify.py --input <fasta> --outdir <dir>
```
