# Provenance

Files in this directory are copied verbatim (unmodified) from the official
VNom repository for reference and reproduction purposes:

- Source: https://github.com/Zheludev/VNom (MIT License)
- `VNom_original_source.py` <- `VNom.py` @ `main`
- `VNom_original_README.md` <- `README.md` @ `main`
- `SRR11060618_subset.fasta` <- `test_data/SRR11060618_subset.fasta` @ `main`

`SRR11060618_subset.fasta` is the tool authors' own test dataset: a subset of
38 rnaSPAdes-assembled contigs from real SRA run `SRR11060618` (human gut
metatranscriptome), used by the Zheludev et al. (2024, "Viroid-like colonists
of human microbiomes") team to validate VNom itself. We use it in Phase 1 to
check that our reimplementation recovers the same class of circular,
sense/antisense-paired ("obelisk-like") contigs that the original tool was
built to find.

We did not vendor or depend on VNom's external binaries (`usearch`,
`circuclust`, `MARS`) — those are non-free/unavailable in this environment.
`src/rnasig/circularity.py` and `src/rnasig/cluster.py` are independent
reimplementations of the same *conceptual* steps (terminal-repeat circularity
detection, circular-permutation-aware clustering, sense/antisense strand
assignment) using pure Python + Biopython's pairwise aligner. Thresholds and
exact algorithmic details differ from the original; see `docs/METHODS.md`.
