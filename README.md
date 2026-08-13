# rnasig

A pipeline for finding new RNA elements in public metatranscriptomes by hunting for a structural signature instead of relying on sequence homology.

## Why this exists

In 2024, Zheludev et al. found Obelisks: a completely new class of ~1kb circular RNA replicon, sitting in public human gut sequencing data that had been public for years. They didn't have better data. They just stopped searching for what things looked like (homology) and started searching for how things behaved (a structural signature: apparent circularity plus both sense and antisense strands of the same molecule showing up in one sample). Their tool is [VNom](https://github.com/Zheludev/VNom).

This repo does two things.

First, it reimplements VNom's core signature from scratch (no `usearch`/`circuclust`/`MARS`, just Python + Biopython) and confirms it actually works by rerunning it on VNom's own test data.

Second, it defines a second, different signature: stable RNA structure combined with low protein-coding potential and dual-strand co-occurrence. Then it calibrates that signature against synthetic controls before letting it loose on anything.

Obelisks themselves would score high on the "structured" and "dual-strand" parts of this second signature, but low on the "non-coding" part (Oblins are real proteins). So it aims at a nearby but distinct part of the same design space that nobody has systematically searched.

## What's in the box

Five scripts, one per phase.

`scripts/phase1_reproduce.py` runs the VNom-style pipeline against VNom's own real test data (SRA run SRR11060618, *Prunus persica* stamen ssRNA-seq). It finds 31 out of 38 contigs are circular, and those collapse into 4 sense/antisense-paired clusters. That's the sanity check that the reimplementation actually works.

Phase 5 has since identified all 38 of those contigs. 24 are **Peach latent mosaic viroid**, a real circular, non-coding, 337 nt RNA replicon; the other 14 are peach host sequence; none is unidentified. The top sense/antisense cluster holds 20 of the 24 viroid contigs and no host sequence at all, and none of the other three clusters holds any viroid. The circularity code resolved the repeat unit as 337 nt from sequence alone, before any database was consulted, which is PLMVd's exact genome length.

So the reimplementation has now been scored against known answers rather than against an expected pattern: 20/20 precision on the viroid cluster, 20/24 recall, genome length exact. That is the strongest evidence in this repo that the signature works.

`scripts/phase2_calibrate.py` builds a synthetic labeled corpus (real positives, coding decoys, structured-coding decoys, plain nulls), scores everything, and reports ROC-AUC, empirical FDR at BH-corrected alpha=0.05, and a power-vs-depth curve. It also documents a real specificity gap that the first version had, and how it got fixed.

`scripts/phase3_sweep.py` is the runner you point at real assembled contigs. It ships with both a `--demo` mode (synthetic corpus, for smoke-testing the pipeline) and a real end-to-end run against SRR13291825 pulled fresh from the AWS SRA Open Data mirror, assembled with MEGAHIT in-sandbox. That run exposed two methodology gaps that Phase 2's synthetic calibration missed: adapter-dimer contigs mimic the signature perfectly, and so do off-target PCR products. Both are now filtered and covered by tests. See `results/phase3_real_sweep_report.md`.

Its actual candidates did not survive Phase 5. **SRR13291825 turned out to be an 18S rRNA amplicon survey of soil DNA, so it could not have contained the kind of molecule this pipeline looks for.** Read `results/phase5_identification_report.md` before reusing anything from that sweep.

`scripts/phase4_characterize.py` is for anything that survives Phase 3: deeper look with more shuffles, a 6-frame ORF map, GC%, circularity resolution, and an SVG rendering of the predicted structure.

`scripts/phase5_identify.py` does the two things the earlier phases kept deferring. It preflights a run's library type before you sweep it, and it resolves surviving candidates against nt, nr and Rfam so "novel" has to survive contact with a database. Running it against the Phase 3 input fails that run in a few seconds:

```
python scripts/phase5_identify.py --run SRR13291825
```

Everything under `src/rnasig/` is the library. `tests/` has 50 pytest tests that run in a few seconds.

## The signature itself

Three axes. All three have to line up.

Structure stability. Fold with ViennaRNA MFE. Compare against dinucleotide-shuffled versions of the same sequence (Altschul-Erikson shuffle, so mono- and dinucleotide composition are preserved exactly). Take the z-score.

Orphan / non-coding. Operationalized reference-free as low ORF coverage plus low in-frame codon bias, so the screen can run without a database. Not a validated gene predictor, and I say so in the code. It is now known to be unreliable in the direction that matters: it rated all three Phase 3 candidates as borderline non-coding and all three encode proteins, one at 97% amino-acid identity. Use `phase5_identify.py` to check anything it likes.

Strand symmetry. Sense and antisense contigs of the same molecule co-occur in the same sample, and/or the contig has a circular terminal repeat.

Ranking uses all three. Significance testing gates on the orphan axis first, then applies BH-FDR on the structure z-score. Phase 2 v1 tested significance on structure alone and quietly let 5 of 8 structured-coding decoys through. Phase 2 v2 fixed that. Both runs are in `results/` so the fix is auditable, not just claimed.

## Setup

```
pip install -e .
pip install -r requirements.txt
pytest tests/ -q
```

Standard stack: biopython, numpy, scipy, pandas, scikit-learn, ViennaRNA, pytest, matplotlib.

## Running it

```
python scripts/phase5_identify.py --run <ACCESSION>          # do this first
python scripts/phase1_reproduce.py
python scripts/phase2_calibrate.py
python scripts/phase3_sweep.py --demo --outdir results/sweep_demo
python scripts/phase4_characterize.py --input candidate.fasta --outdir results/char_x
python scripts/phase5_identify.py --input candidate.fasta --outdir results/identify_x
```

Phase 5 comes first and last on purpose. Check the library type before you spend
an assembly on a run, and check the databases before you call anything novel.

Each phase drops its outputs in `results/`, one markdown report and one JSON.

## Layout

```
src/rnasig/     library
scripts/        the five phase runners
tests/          pytest
data/reference/ VNom's original code + test data, verbatim
data/motifs/    empty on purpose; motifs are generated in code
results/        actual outputs from running this session
docs/           METHODS.md and LIMITATIONS.md
```

## Before you trust the numbers

`docs/LIMITATIONS.md` says plainly what came from real data (Phase 1) and what came from synthetic corpora (Phase 2), and which conclusions have since been withdrawn. Nothing here claims a new organism was found.

The one thing this repo did find is worth stating plainly, because it is a negative result and negative results are easy to bury. Phase 3 ran on a real public run, produced a candidate that survived every filter, and the candidate was a bacterial FAD-dependent oxidoreductase gene fragment in a soil DNA amplicon library. The pipeline's non-coding score rated it, and the two other survivors, as borderline non-coding; all three encode proteins at 73-97% amino-acid identity. `results/phase5_identification_report.md` has the evidence and what changed as a result.
