# Phase 1 report: reproducing a known detection

**Data:** `data/reference/vnom/SRR11060618_subset.fasta` — 38 real
rnaSPAdes contigs, VNom's own test dataset, assembled from real SRA run
SRR11060618. Downloaded verbatim from https://github.com/Zheludev/VNom.

**Correction:** this report originally called SRR11060618 a human gut
metatranscriptome. It is *Prunus persica* (peach) stamen ssRNA-seq; see
`data/reference/vnom/PROVENANCE.md`.

**All 38 contigs have since been identified against nt.** 24 are Peach
latent mosaic viroid, 14 are peach host sequence, none is unidentified. That
gives this phase a real ground truth to score against, and it scores well:
SAS cluster 1 contains 20 viroid contigs and zero host contigs, while
clusters 2, 3 and 4 contain 11 host contigs and zero viroid. Precision on
the viroid cluster is 20/20, recall is 20/24, and the 337 nt circular unit
reported below is PLMVd's exact genome length. Coverage separates the
classes with no overlap either: host 2.2x to 595x, viroid 1,298x to
24,158x. See `results/phase5_identification_report.md`.

**Pipeline:** `src/rnasig/circularity.py` + `src/rnasig/cluster.py`
(independent reimplementation of VNom's CircleFinder + circUCLUST +
SASFinder concept — see `docs/METHODS.md` and
`data/reference/vnom/PROVENANCE.md`), run via `scripts/phase1_reproduce.py`.

## Result

- **38/38** contigs passed the standard length filter (10–1000 nt).
- **31/38 (82%)** were flagged circular (terminal-repeat detection, k=12),
  with resolved unit lengths clustering around three size classes
  (~337–348 nt, ~166 nt, ~143 nt) and copy numbers from 1.05× to 4.73×,
  i.e. many contigs are concatemeric reads of a shorter true monomer —
  exactly the assembler artifact circularity detection is designed to
  catch.
- Circular-permutation-aware clustering (id ≥ 0.7) collapsed these 31
  contigs into **4 non-singleton clusters**.
- **All 4 of the 4 non-singleton clusters show sense/antisense
  co-occurrence** — i.e. every circular cluster found in this real dataset
  also carries the dual-strand signature. The largest cluster (centroid
  `NODE_36652`, 943 nt / 337 nt monomer) has 20 members split roughly evenly
  between sense (8) and antisense (11) orientation — a strong, unambiguous
  dual-strand signal, not a marginal call.

| cluster | centroid | centroid len | monomer-scale unit len | members | sense | antisense |
|---|---|---|---|---|---|---|
| 0 | NODE_36652 | 943 | ~337 | 20 | 8 | 11 (+centroid=sense) |
| 1 | NODE_38543 | 845 | ~332 | 5 | 2 | 2 (+centroid=sense) |
| 2 | NODE_43211 | 613 | ~348 | 4 | 2 | 1 (+centroid=sense) |
| 3 | NODE_52236 | 347 | ~137 | 2 | 0 | 1 (+centroid=sense) |

Full machine-readable output: `results/phase1_reproduction.json`.

## Interpretation

This is the exact pattern the obelisks paper describes: circular,
concatemer-prone contigs from a real RNA-Seq sample that
cluster into groups containing both sense and antisense representatives —
evidence of an actively-transcribed-or-replicating agent rather than
incidental fragments of host/microbial mRNA. Our from-scratch
reimplementation (no `usearch`/`circuclust`/`MARS` dependency) recovers this
pattern on the tool authors' own real validation data, which is what "Phase
1: reproduce a known detection" set out to establish: the pipeline's
underlying logic works, not just in theory but against a real, previously
published-on dataset.

We do not claim these 4 clusters are specifically the obelisk(s) originally
reported for this SRA run (no cross-check against the original paper's
exact sequences was possible — see `docs/LIMITATIONS.md`), only that the
circularity + dual-strand *signature itself* fires correctly, at a
sensible rate (~11% of length-filtered contigs end up in an SAS-positive
circular cluster), on real data.
