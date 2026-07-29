# rnasig

A pipeline for finding new RNA elements in public metatranscriptomes by hunting for a structural signature instead of relying on sequence homology.

## Why this exists

In 2024, Zheludev et al. found Obelisks: a completely new class of ~1kb circular RNA replicon, sitting in public human gut sequencing data that had been public for years. They didn't have better data. They just stopped searching for what things looked like (homology) and started searching for how things behaved (a structural signature: apparent circularity plus both sense and antisense strands of the same molecule showing up in one sample). Their tool is [VNom](https://github.com/Zheludev/VNom).

This repo does two things.

First, it reimplements VNom's core signature from scratch (no `usearch`/`circuclust`/`MARS`, just Python + Biopython) and confirms it actually works by rerunning it on VNom's own test data.

Second, it defines a second, different signature: stable RNA structure combined with low protein-coding potential and dual-strand co-occurrence. Then it calibrates that signature against synthetic controls before letting it loose on anything.

Obelisks themselves would score high on the "structured" and "dual-strand" parts of this second signature, but low on the "non-coding" part (Oblins are real proteins). So it aims at a nearby but distinct part of the same design space that nobody has systematically searched.

## What's in the box

Four scripts, one per phase.

`scripts/phase1_reproduce.py` runs the VNom-style pipeline against VNom's own real gut metatranscriptome test data. It finds 31 out of 38 contigs are circular, and those collapse into 4 sense/antisense-paired clusters. That's the sanity check that the reimplementation actually works.

`scripts/phase2_calibrate.py` builds a synthetic labeled corpus (real positives, coding decoys, structured-coding decoys, plain nulls), scores everything, and reports ROC-AUC, empirical FDR at BH-corrected alpha=0.05, and a power-vs-depth curve. It also documents a real specificity gap that the first version had, and how it got fixed.

`scripts/phase3_sweep.py` is the runner you point at real assembled contigs. It ships with both a `--demo` mode (synthetic corpus, for smoke-testing the pipeline) and a real end-to-end run against SRR13291825 pulled fresh from the AWS SRA Open Data mirror, assembled with MEGAHIT in-sandbox. That real run found a legitimate signal (structure z=5, sense/antisense-paired, high-coverage RNA), and — importantly — also revealed a real methodology gap (adapter-dimer contigs mimic the signature perfectly) that Phase 2's synthetic calibration missed. That gap is now fixed (`src/rnasig/preprocess.py`, adapter/homopolymer filter, on by default) and covered by tests. See `results/phase3_real_sweep_report.md` for the full story.

`scripts/phase4_characterize.py` is for anything that survives Phase 3: deeper look with more shuffles, a 6-frame ORF map, GC%, circularity resolution, and an SVG rendering of the predicted structure.

Everything under `src/rnasig/` is the library. `tests/` has 20 pytest tests that run in a few seconds.

## The signature itself

Three axes. All three have to line up.

Structure stability. Fold with ViennaRNA MFE. Compare against dinucleotide-shuffled versions of the same sequence (Altschul-Erikson shuffle, so mono- and dinucleotide composition are preserved exactly). Take the z-score.

Orphan / non-coding. No BLAST in this sandbox, so "matches nothing known" gets operationalized reference-free as low ORF coverage plus low in-frame codon bias. Not a validated gene predictor, and I say so in the code.

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
python scripts/phase1_reproduce.py
python scripts/phase2_calibrate.py
python scripts/phase3_sweep.py --demo --outdir results/sweep_demo
python scripts/phase4_characterize.py --input candidate.fasta --outdir results/char_x
```

Each phase drops its outputs in `results/`, one markdown report and one JSON.

## Layout

```
src/rnasig/     library
scripts/        the four phase runners
tests/          pytest
data/reference/ VNom's original code + test data, verbatim
data/motifs/    empty on purpose; motifs are generated in code
results/        actual outputs from running this session
docs/           METHODS.md and LIMITATIONS.md
```

## Before you trust the numbers

`docs/LIMITATIONS.md` says plainly what came from real data (Phase 1) and what came from synthetic corpora (Phases 2 and 3), and why the network policy in this sandbox made a real Phase 3 sweep impossible. Nothing here claims a new organism was found. It's a candidate-generation pipeline, calibrated honestly, run on what data was actually reachable.
