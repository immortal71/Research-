# rnasig — a signature-based hunt for undescribed RNA biology

Sequence-homology search can only find things that resemble something
already in a database. Zheludev et al. (2024) found **Obelisks** — a whole
new class of ~1kb circular RNA replicon, hiding in public human gut
metatranscriptomes — by searching for a *structural signature* instead:
apparent circularity plus co-occurrence of sense and antisense strands of
the same molecule in one sample. Their tool is
[VNom](https://github.com/Zheludev/VNom).

This repo follows the same four-phase playbook on a **second, different**
signature, aimed at a part of the same design space obelisks don't occupy:
stable RNA structure with low protein-coding potential, still co-occurring
as sense/antisense pairs or circular contigs — "sequences whose predicted
structure is stable but whose sequence matches nothing."

**Read `docs/LIMITATIONS.md` first.** It states plainly which results here
come from real data and which are synthetic, and exactly why (this sandbox
has no route to NCBI/SRA/ENA/EBI/Zenodo/OSF/Dryad/figshare — confirmed by
direct testing, not assumed). Nothing in this repo claims a new organism
was found; it is a calibrated candidate-generation pipeline, run honestly
on what data was actually reachable.

## The four phases

| Phase | Script | What it does | Data |
|---|---|---|---|
| 1. Reproduce | `scripts/phase1_reproduce.py` | Reimplements VNom's circularity + circular-permutation clustering + sense/antisense detection; runs it on VNom's own real test data | Real (SRR11060618, human gut) |
| 2. Calibrate | `scripts/phase2_calibrate.py` | Defines the novel signature (see below) and measures its ROC-AUC, empirical FDR/power, and power-vs-depth curve against synthetic labeled positives + three flavors of matched decoys | Synthetic, labeled |
| 3. Sweep | `scripts/phase3_sweep.py` | Runs the calibrated pipeline against any FASTA of assembled contigs | `--demo` (synthetic; see limitations) or `--input <real contigs>` |
| 4. Characterize | `scripts/phase4_characterize.py` | Deep-dives a flagged candidate: refolded structure + z-score, 6-frame ORF map, circularity status, SVG structure plot | Whatever Phase 3 flags |

Full methodology: `docs/METHODS.md`. Results: `results/*.md` and
`results/*.json`.

## The novel signature: SOS (Structure-stable, Orphan-sequence, Strand-symmetric)

Three orthogonal axes, combined into one score and calibrated with
Benjamini-Hochberg FDR control:

- **S1 — structure stability**: MFE (ViennaRNA) z-score vs. dinucleotide-
  shuffled nulls of the same sequence (Altschul-Erikson shuffle).
- **S2 — orphan/coding-potential**: low ORF coverage + low in-frame codon
  bias (reference-free proxy — see caveat in `src/rnasig/orphan.py`).
- **S3 — strand-symmetry**: sense+antisense co-occurrence and/or circular
  terminal repeat.

Obelisks themselves would score high on S1 and S3, but low on S2 (Oblins
are real proteins) — this signature deliberately hunts the complementary,
non-coding-but-structured corner nobody has systematically screened for.

## Setup

```
pip install -e .
pip install -r requirements.txt   # or just: pip install -r requirements.txt
pytest tests/ -q
```

## Quick start

```
python scripts/phase1_reproduce.py
python scripts/phase2_calibrate.py
python scripts/phase3_sweep.py --demo --outdir results/sweep_demo
python scripts/phase4_characterize.py --input <candidate.fasta> --outdir results/char_x
```

## Repository layout

```
src/rnasig/        core library (seqio, circularity, cluster, structure, orphan, signature, simulate, nullmodel)
scripts/           the four phase-runner scripts
tests/             pytest unit tests (run: pytest tests/ -q)
data/reference/    VNom's own source + real test data (verbatim, see PROVENANCE.md)
data/motifs/       (no static files -- motifs are generated in code, see PROVENANCE.md)
results/           actual output of every phase run in this session
docs/METHODS.md          full methodology
docs/LIMITATIONS.md      what's real, what's synthetic, and why -- read this
```
