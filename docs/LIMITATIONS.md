# Limitations and honest scope

This is the single most important document in the repo. Read it before
trusting any headline number.

## Withdrawn: the Phase 3 real sweep and everything downstream of it

SRR13291825, the run Phase 3 was pointed at, is an **18S rRNA amplicon
survey of soil DNA** (`library_strategy=AMPLICON`, `library_selection=PCR`,
`library_source=METAGENOMIC`). It is a PCR metabarcoding library. It
contains no RNA, its coverage numbers reflect amplification rather than
abundance, and the concatemers PCR generates are read by the circularity
detector as terminal repeats.

The pipeline ran correctly on it and should never have been pointed at it.
All three surviving candidates are fragments of bacterial protein-coding
genes:

| contig | encodes | aa identity |
|---|---|---|
| k141_25 | FAD-dependent oxidoreductase (*Pseudolabrys* / *Bradyrhizobium*) | 97% |
| k141_10 | DUF2007 signal-transducing protein (Acidobacteriota) | 88% |
| k141_3 | MBL fold metallo-hydrolase (Tepidisphaeraceae) | 73% |

Rfam returns nothing for any of them. `results/phase5_identification_report.md`
carries the full evidence, including the raw-read work that establishes the
library type independently of the archive metadata.

What this costs the repo: k141_25 is withdrawn as a candidate novel
circular RNA, k141_10 and k141_3 are withdrawn entirely, and Phase 4's
characterization of them describes PCR products. What it does not cost:
Phase 1 (real RNA-Seq, unaffected) and Phase 2 (always explicitly
synthetic).

`scripts/phase5_identify.py --run <ACC>` now performs this check before a
sweep rather than after one.

## What was actually run, on what data

| Phase | Data | Real or synthetic? |
|---|---|---|
| 1. Reproduce | `SRR11060618_subset.fasta` | **Real** — VNom authors' own test set, 38 contigs from a real RNA-Seq run (*Prunus persica* stamen, not the human gut metatranscriptome earlier drafts claimed), downloaded verbatim from https://github.com/Zheludev/VNom |
| 2. Calibrate | synthetic labeled corpora | **Synthetic by design** — calibration requires dial-able ground truth (copy number, mutation rate, decoy type), which is exactly what a real sample cannot give you |
| 3. Sweep (demo) | labeled synthetic "demo" corpus | **Synthetic, explicitly labeled as such in the script's own output** — see below for the original network-blocked scoping |
| 3. Sweep (real) | SRR13291825 assembly | **Real data, wrong kind.** Pulled from the AWS SRA Open Data mirror and assembled with MEGAHIT, but the run is an 18S rRNA DNA amplicon library. Conclusions withdrawn, see above. |
| 4. Characterize | whatever Phase 3 flags | inherits Phase 3's provenance, so the SRR13291825 characterizations are withdrawn too |
| 5. Identify | ENA metadata, raw reads, nt / nr / Rfam | **Real** — live queries, reproducible with `scripts/phase5_identify.py` |

## Why Phase 3 did not run on real public metatranscriptome data

**This section describes the sandbox as it was when Phase 3 was written, and
it no longer holds.** From the environment Phase 5 ran in,
`blast.ncbi.nlm.nih.gov`, `eutils.ncbi.nlm.nih.gov`, `www.ebi.ac.uk`,
`ftp.sra.ebi.ac.uk` and `batch.rfam.org` are all reachable, which is what
made the identification work possible. The reachability notes below are kept
because they explain why SRR13291825 was chosen, which is the decision that
went wrong. Nothing about a restricted network justified skipping the
library-type check; that check needs one metadata request, and the archive
that served the data would have answered it.

This session's outbound network access goes through a policy-gated proxy.
We tested reachability directly (not guessing) before deciding how to scope
Phase 3:

**Blocked (HTTP CONNECT tunnel failed / 403, confirmed by direct test):**
`www.ncbi.nlm.nih.gov`, `pmc.ncbi.nlm.nih.gov`, `ftp.ebi.ac.uk`,
`ftp.sra.ebi.ac.uk`, `academic.oup.com` (hosts the marine-obelisks
supplementary FASTA we were specifically looking for), `zenodo.org`,
`api.datacite.org`, `portal.nersc.gov`, `figshare.com`, `osf.io`,
`datadryad.org`.

**Reachable (confirmed 2xx/expected response):**
`github.com`, `raw.githubusercontent.com`, PyPI/npm (used to install
biopython/numpy/scipy/scikit-learn/ViennaRNA/pytest/matplotlib), AND —
this turned out to matter more than I initially thought — the **AWS SRA
Open Data mirrors**: `sra-pub-run-odp.s3.amazonaws.com`,
`sra-pub-src-1.s3.amazonaws.com`, `sra-pub-metadata-us-east-1.s3.amazonaws.com`.
The AWS mirror serves any SRA accession directly (URL pattern:
`https://sra-pub-run-odp.s3.amazonaws.com/sra/<ACC>/<ACC>`), no NCBI
eutils needed. This is what enabled the real Phase 3 sweep of
SRR13291825; see `results/phase3_real_sweep_report.md`.

Net effect: every standard public sequence-archive route (SRA, ENA, GenBank
via eutils, PMC full text, journal supplementary data, Zenodo/OSF/Dryad/
figshare deposits) is unreachable from this sandbox. We could, and did, pull
one real dataset because it happened to be small enough to be checked into
a GitHub repo (`Zheludev/VNom`'s own test data) — that is what powers Phase
1. We found no comparable real, small, non-gut environmental
(soil/marine/extremophile) contig set mirrored on GitHub in the time
available; `microbiomedata/metaT_Assembly`'s soil test data, for instance,
is hosted at `portal.nersc.gov`, which is blocked.

**What this means concretely:** Phase 3's "sweep" in this repo is run with
`--demo`, against a synthetic, explicitly-labeled corpus
(`simulate.build_calibration_corpus`), specifically so that:

1. the pipeline itself is demonstrably exercised end-to-end exactly as it
   would be on real data,
2. no output anywhere in this repo could be mistaken for "we found a new
   organism" — there is no such claim here, and there shouldn't be, given
   what was actually possible to run.

## What a real Phase 3 run needs

`scripts/phase3_sweep.py --input <contigs.fasta> --outdir <dir>` is
data-agnostic and ready to run as soon as real assembled contigs are
available. Concretely, from an environment with SRA/ENA access:

0. **Preflight the run**: `python scripts/phase5_identify.py --run <ACC>`.
   The run must be shotgun RNA. An AMPLICON strategy, a PCR selection or a
   GENOMIC/METAGENOMIC source disqualifies it, and so does missing metadata.
   If you have the FASTQ, pass `--reads` as well, so a mislabelled run is
   still caught by its own primer distribution.
1. Pick an under-sampled environment/host category (soil, hydrothermal
   vent, insect, plant endophyte, extremophile culture, non-mammalian
   vertebrate gut — anything far from the human-gut-heavy sampling to date).
2. Pull stranded RNA-seq runs (`fasterq-dump`/ENA direct FASTQ), assemble
   with `rnaSPAdes` (matching VNom's own assumption, though any de Bruijn
   assembler that reports contig length in its ID, or any FASTA, works with
   minor adjustment).
3. `python scripts/phase3_sweep.py --input <assembly>.fasta --outdir results/sweep_<env>`
4. Anything BH-significant goes to `scripts/phase4_characterize.py` for a
   deeper look, then to `scripts/phase5_identify.py --input` so that any
   claim of novelty has to survive nt, nr and Rfam first, and only from
   there to wet-lab or comparative-genomics
   follow-up — this pipeline only ever produces *candidates*, never a
   discovery claim by itself, exactly as VNom itself is explicit about
   nominating candidates rather than confirming a new taxon.

## Other scope limits, stated plainly

- **The reimplementation is not VNom.** It captures the same conceptual
  steps (circularity, circular-permutation clustering, SAS detection) but
  uses different tools (Biopython's aligner instead of usearch/MARS/
  circuclust) and different thresholds. Phase 1's agreement with the
  expected obelisk-like pattern in VNom's own test data is evidence the
  *concept* transfers, not a claim of byte-for-byte equivalence to the
  published tool.
- **The orphan/coding-potential score (S2) is a lightweight, unsupervised
  proxy**, not a validated gene predictor. It is calibrated only against our
  own synthetic decoys (Phase 2), which is real evidence of what it does on
  those decoys and not a general-purpose coding-potential classifier.
  It has now met real sequence once, and it got all three cases wrong:
  k141_3, k141_10 and k141_25 scored 0.58, 0.55 and 0.52 (above 0.5 reads
  as non-coding) and all three encode proteins, the worst miss being the
  one at 97% amino-acid identity. Nucleotide identity for that contig is
  only 78%, because synonymous substitution hides conservation at the
  nucleotide level and not at the protein level. Treat a `dbcheck` blastx
  pass as required before believing an S2 score, not optional.
- **The motif library used for calibration positives is synthetic by
  construction** (`simulate.make_rod`, `make_stem_loop`,
  `make_cloverleaf_like`) — guaranteed to fold with real thermodynamic
  stability (verified computationally, not just asserted), but explicitly
  not a claim that these are natural sequences. This is deliberate: Phase 2
  calibration needs ground truth we can dial, which synthetic construction
  gives us and no real dataset could.
- **No wet-lab or comparative-genomic confirmation of anything** — this is
  a computational candidate-generation pipeline, exactly analogous to what
  VNom itself is (a "nominator," in its own authors' words), not a
  discovery-confirmation system.
