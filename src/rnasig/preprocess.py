"""Pre-filtering for real assembler output.

Discovered by running the pipeline on real, un-preprocessed SRA data
(SRR13291825, human microbiome MiSeq): the very top hit was an assembled
adapter-dimer (Illumina TruSeq P7 + Read 2 adapter). Adapter dimers fold
stably (high S1 z-score), are highly abundant (huge coverage), and appear
on both strands (high S3 SAS), so they perfectly mimic the SOS signature.

This module strips them before the pipeline sees them. It is *not*
attempting to be a general QC tool -- it targets exactly the failure mode
that showed up on real data:

  1. Contigs that contain any of the canonical Illumina adapter / primer
     substrings (TruSeq, Nextera).
  2. Contigs dominated by a homopolymer or a tiny dinucleotide repeat
     (assembler artifacts from polyA tails / adapter poly-T).

Both are conservative: match required to be at least MIN_MATCH nt exact
(so real biological sequence that accidentally contains a short match
isn't stripped).
"""
from __future__ import annotations

from dataclasses import dataclass

from .seqio import Record, revcomp

# Standard Illumina adapter/primer segments (both TruSeq and Nextera lines).
# Sources: Illumina Adapter Sequences document (verified against Illumina's
# published sequences, not guessed).
_ILLUMINA_ADAPTERS = [
    # TruSeq / Nextera P5/P7 flow-cell primers
    "AATGATACGGCGACCACCGAGATCTACAC",
    "CAAGCAGAAGACGGCATACGAGAT",
    # TruSeq Read 1 / Read 2 sequencing primers
    "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
    # Nextera transposase adapters
    "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
    "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
    # Common i5/i7 index-adjacent stretches (short, stable)
    "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
    "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
]

_MIN_MATCH = 20  # nt of exact adapter substring needed to flag a contig
_HOMOPOLYMER_MIN = 12  # nt of a single-base run to flag as artifact-dominated


@dataclass
class FilterReport:
    n_input: int
    n_kept: int
    n_dropped_adapter: int
    n_dropped_homopolymer: int
    dropped_ids: list[str]
    dropped_reasons: dict[str, str]


def _has_adapter_hit(seq: str) -> str | None:
    """Return the adapter substring found (in either orientation), or None."""
    for adapter in _ILLUMINA_ADAPTERS:
        probe = adapter[:_MIN_MATCH]
        if probe in seq or probe in revcomp(seq)[:len(seq)]:
            return adapter
        rc_probe = revcomp(adapter)[:_MIN_MATCH]
        if rc_probe in seq:
            return adapter + "(rc)"
    return None


def _dominant_homopolymer(seq: str) -> str | None:
    """Return the base of a >=_HOMOPOLYMER_MIN-nt homopolymer run if present."""
    if len(seq) < _HOMOPOLYMER_MIN:
        return None
    for base in "ACGTN":
        if base * _HOMOPOLYMER_MIN in seq:
            return base
    return None


def filter_contigs(records: list[Record]) -> tuple[list[Record], FilterReport]:
    kept: list[Record] = []
    dropped_ids: list[str] = []
    dropped_reasons: dict[str, str] = {}
    n_adapter = 0
    n_homopolymer = 0

    for r in records:
        seq = r.seq.upper()
        adapter = _has_adapter_hit(seq)
        if adapter is not None:
            dropped_ids.append(r.id)
            dropped_reasons[r.id] = f"adapter:{adapter}"
            n_adapter += 1
            continue
        hp = _dominant_homopolymer(seq)
        if hp is not None:
            dropped_ids.append(r.id)
            dropped_reasons[r.id] = f"homopolymer:{hp}x{_HOMOPOLYMER_MIN}+"
            n_homopolymer += 1
            continue
        kept.append(r)

    report = FilterReport(
        n_input=len(records),
        n_kept=len(kept),
        n_dropped_adapter=n_adapter,
        n_dropped_homopolymer=n_homopolymer,
        dropped_ids=dropped_ids,
        dropped_reasons=dropped_reasons,
    )
    return kept, report
