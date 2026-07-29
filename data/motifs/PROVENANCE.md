# Motif library provenance

This directory intentionally contains no static sequence files.

All calibration motifs (structured "rod" positives, hairpins, cloverleaf-like
multi-stem constructs, codon-biased ORF decoys) are generated
programmatically by `src/rnasig/simulate.py` at run time, from a seeded
`random.Random`, so every calibration run is exactly reproducible from code
alone. See `docs/METHODS.md` (Phase 2) and `docs/LIMITATIONS.md` for why
these are synthetic constructions rather than sequences taken from a
database, and what is/isn't claimed about their biological realism.
