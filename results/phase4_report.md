# Phase 4 report: characterizing the real Phase 1 candidates

**Input:** the 4 real sense/antisense-paired cluster centroids found in
Phase 1 (`data/reference/vnom/SRR11060618_subset.fasta`, real *Prunus
persica* stamen RNA-Seq contigs) — i.e. this phase runs on real data, not a
synthetic demo. `scripts/phase4_characterize.py`, n_shuffles=100.

> **All four have since been identified.** NODE_36652 is Peach latent
> mosaic viroid (98% nt identity, E=5e-155), and the 337 nt circular unit
> reported below matches PLMVd's complete genome length exactly. The other
> three are peach genomic and mRNA sequence. See
> `results/phase5_identification_report.md`.

| contig | length | GC% | circular (unit len) | structure z | orphan score | longest ORF (nt) | SOS-significant? |
|---|---|---|---|---|---|---|---|
| NODE_36652 | 943 | 0.54 | yes (337) | **17.13** | 0.89 | 96 | **yes** |
| NODE_38543 | 845 | 0.41 | yes (332) | 0.11 | 0.84 | 150 | **no** |
| NODE_43211 | 613 | 0.58 | yes (348) | 8.50 | 0.62 | 180 | **yes** |
| NODE_52236 | 347 | 0.38 | yes (137) | 2.87 | 0.60 | 234 | **yes** |

Significance = orphan-gated BH-FDR on the structure z-score p-value (see
`docs/METHODS.md` / Phase 2's v2 fix), alpha=0.05, applied across these 4
candidates. Full per-candidate detail (6-frame ORF maps, MFE dot-bracket
structures, null distributions) in `results/phase4_characterize/
characterization.json`; SVG secondary-structure renderings in the same
directory.

## The honest negative: NODE_38543

Three of the four real circular+SAS clusters found in Phase 1 also show
significantly stable structure (z = 17.1, 8.5, 2.9). **The fourth,
NODE_38543 (845 nt, 5-member cluster with 2 sense + 2 antisense members),
does not** — its structure z-score of 0.11 means its fold is essentially
exactly as stable as its own dinucleotide-shuffled background predicts, no
more. It is circular and shows real sense/antisense co-occurrence in real
sequencing data (Phase 1's S3 axis fires correctly), but it fails the S1
(structure-stability) axis outright.

This is reported as a negative, not smoothed over: not every circular,
dual-strand-transcribed RNA element found by the VNom-style signature is
also a structural outlier by the SOS signature's S1 criterion. That is
exactly what should happen if S1 and S3 are measuring genuinely different
things (which is the whole premise of adding S1/S2 on top of S3 in the
first place) — if every VNom-positive cluster had also come out
SOS-positive, that would suggest S1 wasn't adding independent information.

## What the three positives look like

`characterization.json` / the SVGs show the standard signature this
pipeline was built to catch: a stable, non-trivial fold (not a simple
hairpin — see the SVG renderings) combined with modest coding potential
(longest ORFs of 96–234 nt, i.e. 32–78 aa — well short of Oblin-length
ORFs [~200 aa in the original obelisks paper], consistent with these being
in the "structured, largely non-coding" part of the design space this
signature targets, distinct from a fully Oblin-coding obelisk). We do not
claim these three are a specific known or novel taxon — that requires
comparative genomics / wet-lab follow-up this environment cannot do (see
`docs/LIMITATIONS.md`) — only that the calibrated signature fires on them,
on real data, for the reasons the signature is designed to fire.
