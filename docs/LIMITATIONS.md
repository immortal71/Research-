# Limitations and honest scope

This is the single most important document in the repo. Read it before
trusting any headline number.

## What was actually run, on what data

| Phase | Data | Real or synthetic? |
|---|---|---|
| 1. Reproduce | `SRR11060618_subset.fasta` | **Real** — VNom authors' own test set, 38 contigs assembled from a real human gut metatranscriptome SRA run, downloaded verbatim from https://github.com/Zheludev/VNom |
| 2. Calibrate | synthetic labeled corpora | **Synthetic by design** — calibration requires dial-able ground truth (copy number, mutation rate, decoy type), which is exactly what a real sample cannot give you |
| 3. Sweep | labeled synthetic "demo" corpus | **Synthetic, explicitly labeled as such in the script's own output** — see below for why |
| 4. Characterize | whatever Phase 3 flags | inherits Phase 3's provenance |

## Why Phase 3 did not run on real public metatranscriptome data

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
`github.com`, `raw.githubusercontent.com`, bare cloud-storage roots
(`s3.amazonaws.com`, `sra-pub-run-odp.s3.amazonaws.com`,
`storage.googleapis.com` — reachable as hosts, but we do not have exact
object keys for specific SRA runs without eutils/SRA metadata access,
which is itself blocked), and PyPI/npm (used to install
biopython/numpy/scipy/scikit-learn/ViennaRNA/pytest/matplotlib).

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

1. Pick an under-sampled environment/host category (soil, hydrothermal
   vent, insect, plant endophyte, extremophile culture, non-mammalian
   vertebrate gut — anything far from the human-gut-heavy sampling to date).
2. Pull stranded RNA-seq runs (`fasterq-dump`/ENA direct FASTQ), assemble
   with `rnaSPAdes` (matching VNom's own assumption, though any de Bruijn
   assembler that reports contig length in its ID, or any FASTA, works with
   minor adjustment).
3. `python scripts/phase3_sweep.py --input <assembly>.fasta --outdir results/sweep_<env>`
4. Anything BH-significant goes to `scripts/phase4_characterize.py` for a
   deeper look, and from there to actual wet-lab or comparative-genomics
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
  proxy**, not a validated gene predictor (no Pfam/Rfam/nr access to build
  or check one against). It is validated only against our own synthetic
  decoys (Phase 2), which is real evidence of what it does on those
  decoys, not a general-purpose coding-potential classifier.
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
