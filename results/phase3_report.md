# Phase 3 report: sweep

**This is a labeled synthetic demonstration, not a real environmental
sweep** — see `docs/LIMITATIONS.md` for the network-reachability testing
that led to this scoping decision (every standard sequence-archive route —
SRA, ENA, GenBank eutils, PMC, journal supplementary data, Zenodo/OSF/
Dryad/figshare — was unreachable from this sandbox; only GitHub-hosted
small files and bare cloud-storage roots were).

Run: `python scripts/phase3_sweep.py --demo --outdir results/sweep_demo
--n-shuffles 15`

## What the demo corpus contains

79 synthetic contigs: 3 true "novel elements" (1–2 replicate sense/
antisense contig pairs each, mimicking a rare, low-depth real detection —
deliberately the *hard* case, not the easy one used in Phase 2's main
run), 25 coding decoys, 8 structured-coding decoys, and 40 plain
background nulls. Ground truth labels are written to
`results/sweep_demo/demo_ground_truth_labels.json` for this sanity check
only — the script explicitly notes that a real sweep would have no such
file, because there is no ground truth for real data.

## Result

- 79 contigs → 70 clusters (9 non-singleton) after circular-permutation
  clustering.
- **4 candidates called significant at alpha=0.05.**
- **3 of the top 3 hits by combined score are the true synthetic
  elements** (`POS000`, `POS002`, `POS001` — all 3 planted elements
  recovered, ranked strictly above every decoy).
- **1 false call**: `NULL005`, a plain background contig, structure
  z=3.04 — a single-contig, non-SAS, non-circular false positive,
  consistent with the ~11% empirical FDR measured in Phase 2's
  calibration (this demo uses only 15 shuffles for speed, noisier than
  Phase 2's already-reduced 10–shuffle main run, so this is within the
  expected noise band, not a new failure mode).
- Every coding decoy and every structured-coding decoy in this run was
  correctly rejected (0 false calls from either class) — consistent with
  Phase 2's orphan-gate fix.

Full output: `results/sweep_demo/sweep_results.json`.

## Reading this result honestly

This demonstrates the pipeline does what Phase 1 and Phase 2 calibrated it
to do — including on a *harder* scenario (1–2 replicates instead of
Phase 2's 2–4) than the main calibration run, and it still recovered all
3 planted elements. It is not a discovery: the "elements" are synthetic by
construction (`docs/LIMITATIONS.md`). The honest claim this phase supports
is narrower and more useful: *given real assembled contigs from an
under-sampled environment, this exact command
(`scripts/phase3_sweep.py --input <contigs.fasta> --outdir <dir>`) will
run and report calibrated candidates, at roughly the sensitivity/
specificity operating point measured in Phase 2.*
