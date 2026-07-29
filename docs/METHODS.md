# Methods

## Background: why this is tractable at all

Zheludev et al. (2024) found Obelisks — a new class of ~1kb circular RNA
replicon, invisible to homology search — not by better sequencing but by
searching public gut metatranscriptomes for a *structural signature*
(apparent circularity + co-occurrence of sense and antisense strands of the
same molecule in one sample) instead of matching against known sequence
databases. Their tool, VNom, is on GitHub
(https://github.com/Zheludev/VNom, MIT license). This project reimplements
that conceptual signature independently (see
`data/reference/vnom/PROVENANCE.md` for exactly what was and wasn't reused),
validates the reimplementation against VNom's own real test data, then
defines and calibrates a second, deliberately *different* signature aimed at
a part of the same design space obelisks don't occupy.

## Phase 1 — Reproduce a known detection

`src/rnasig/circularity.py` and `src/rnasig/cluster.py` reimplement VNom's
three core ideas without VNom's external, non-free dependencies (`usearch`,
`circuclust`, `MARS`):

1. **Circularity** (`circularity.find_circularity`): an RNA circle that gets
   assembled linearly produces a contig whose 3' end duplicates an earlier
   stretch of its own 5' end (the assembler read through the origin before
   its overlap graph closed). We detect this by exact terminal k-mer
   matching (k=12 nt by default) and report the resolved "monomer" (one
   unit-length copy) and overhang length. `circularity.resolve_concatemer`
   collapses near-integer multiples (tandem duplicates) to their monomer.
2. **Circular-permutation-aware clustering** (`cluster.cluster_sequences`):
   two circular molecules can look unrelated as raw strings if the
   assembler started reading at different points around the circle. We
   align the query against the *target doubled with itself*
   (`target + target`) using Biopython's `PairwiseAligner` in local mode, in
   both orientations (forward and reverse-complement), and take the best
   identity — this tolerates arbitrary rotation and reports which strand
   orientation matched.
3. **Sense/antisense (SAS) detection** (`cluster.Cluster.has_sense_antisense`):
   within a cluster, if any member's best-matching orientation to the
   centroid is "antisense" while others are "sense", the cluster contains
   both polarities of the same molecule in one sample — the paper's
   proposed signature of active (not just incidental) transcription, or
   outright replication.

`scripts/phase1_reproduce.py` runs this on
`data/reference/vnom/SRR11060618_subset.fasta` — VNom's own real test
dataset (a subset of 38 rnaSPAdes contigs from a real SRA gut
metatranscriptome run). Results are in `results/phase1_reproduction.json`
and `results/phase1_report.md`.

## Phase 2 — Define and calibrate a novel signature

**The signature ("SOS": Structure-stable, Orphan-sequence,
Strand-symmetric)** has three orthogonal axes, deliberately chosen so that
*obelisks themselves would score high on two axes but not the third* — this
signature targets a different, still-unsearched region of the same design
space, per the project brief's suggestion of "sequences whose predicted
structure is stable but whose sequence matches nothing":

- **S1 — structure stability** (`structure.structure_zscore`): fold with
  ViennaRNA's `RNA.fold` (MFE), then compare the folding energy to the MFE
  distribution of the same sequence's dinucleotide-composition-matched
  shuffles (`nullmodel.dinuc_shuffle`, the Altschul-Erikson 1985
  algorithm). A large positive z-score means the fold is more stable than
  the sequence's own base/dinucleotide composition would predict by chance
  — the standard single-sequence structure-outlier statistic (cf. Workman &
  Krogh 1999; the per-sequence analogue of RNAz's z-score, needed here
  because there is no reference alignment for a genuinely novel element).
- **S2 — orphan/coding-potential score** (`orphan.orphan_score`): with no
  BLAST/nr access in this environment, "matches nothing" is operationalized
  reference-free as *low protein-coding potential* (short longest-ORF
  coverage + low in-frame codon-usage bias relative to the sequence's own
  other frames). This is a lightweight, unsupervised proxy, not a validated
  gene predictor — see the caveat in `orphan.py`'s docstring. It is what
  separates this signature from obelisks: Oblins are real proteins, so
  obelisks would score *low* on this axis; we are hunting the
  complementary, non-coding-but-structured space.
- **S3 — strand-symmetry / circularity** (reused from Phase 1): SAS
  co-occurrence and/or a circular terminal repeat.

Combined score = `structure_z * orphan_score`, +2 if SAS-positive, +1 if
circular (`signature.score_cluster`). A one-sided p-value is derived from
the structure z-score under the standard normal (the shuffles already
define the null), and Benjamini-Hochberg FDR control
(`signature.benjamini_hochberg`) turns p-values into a calibrated call set.

**Calibration** (`scripts/phase2_calibrate.py`,
`simulate.build_calibration_corpus`): synthetic, labeled contig pools built
from:

- true positive "elements": programmatically constructed rod-like
  multi-hairpin RNAs (`simulate.make_rod` — a guaranteed-structured, but not
  claimed-to-be-biological, stand-in for the "genome-spanning rod fold"
  reported for obelisks), each represented by 1–5 replicate sense/antisense
  contig pairs with mutation (`simulate.make_sas_pair`) to mimic real
  assembly noise, a fraction also circularized
  (`simulate.make_circular_contig`).
- **matched decoys**, specifically designed to stress-test each axis's
  necessity rather than just its sensitivity:
  - *coding decoys* (`simulate.make_orf_like`): strong single-ORF, codon-biased,
    unstructured — should fail S1/S2.
  - *structured-coding decoys*: an ORF fused to a real stable hairpin — tests
    that S2 correctly suppresses "just structured" sequences that are
    nonetheless ordinary coding transcripts (structure alone must not be
    sufficient; this is the specificity control for S1).
  - *plain nulls*: dinucleotide-shuffled Markov background, single copy, no
    partner strand.

Metrics reported in `results/phase2_calibration.json` /
`results/phase2_report.md`: ROC-AUC / PR-AUC of the combined score against
ground truth, empirical FDR and power under BH control at alpha=0.05, a
power-vs-depth curve (recovery rate as replicate/sequencing depth varies),
and a false-call breakdown by decoy type (does S2 actually suppress
structured-coding decoys that S1 alone would have flagged?).

## Phase 3 — Sweep an environment

`scripts/phase3_sweep.py` runs the identical, already-calibrated pipeline
(length filter → circularity → circular-permutation clustering → SAS →
S1/S2/S3 scoring → BH-FDR) against any FASTA of assembled contigs. See
`docs/LIMITATIONS.md` for why this session's actual "sweep" is a labeled
synthetic demonstration rather than a real environmental dataset, and for
exactly what a real run needs.

## Phase 4 — Characterize

`scripts/phase4_characterize.py` takes a candidate sequence and produces a
deeper report: refolded MFE structure (more shuffles, tighter z-score
estimate), a 6-frame ORF map, GC content, circularity/concatemer status,
and an SVG rendering of the predicted secondary structure via ViennaRNA's
`svg_rna_plot`.
