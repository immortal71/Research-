"""Axis S2 of the novel signature: "sequence matches nothing" (recognizable).

We have no BLAST/nr access in this environment, so "matches nothing" cannot
mean "no BLAST hit" here -- it is operationalized reference-free, as low
protein-coding potential: viroid-like/ribozyme-like elements are
non-coding, so they lack the long, statistically biased open reading
frames that mark real protein-coding sequence. This is a lightweight,
unsupervised proxy (ORF-coverage + in-frame codon-bias), *not* a validated
gene predictor -- it is meant only to separate "looks like it encodes a
protein" from "looks like bare structured/regulatory RNA" within our own
calibration, and to flag candidates for closer characterization (Phase 4),
not to make a definitive coding/non-coding call.
"""
from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass

from .seqio import dna, revcomp

_STOPS = {"TAA", "TAG", "TGA"}
_START = "ATG"


def _codons(seq: str, frame: int) -> list[str]:
    seq = seq[frame:]
    return [seq[i : i + 3] for i in range(0, len(seq) - 2, 3)]


def longest_orf_nt(seq: str) -> int:
    """Longest open reading frame (ATG...stop, in nt including stop) across
    all 6 frames (3 forward, 3 reverse-complement)."""
    seq = dna(seq)
    best = 0
    for strand_seq in (seq, revcomp(seq)):
        for frame in range(3):
            codons = _codons(strand_seq, frame)
            start_idx = None
            for i, c in enumerate(codons):
                if start_idx is None:
                    if c == _START:
                        start_idx = i
                elif c in _STOPS:
                    orf_len = (i - start_idx + 1) * 3
                    best = max(best, orf_len)
                    start_idx = None
            # unterminated ORF running to the end of the frame: don't count
            # (assemblies are often fragments, but we want a conservative
            # coding-potential estimate, not to reward truncation)
    return best


def codon_bias_chisq(seq: str) -> float:
    """Chi-square statistic comparing the best-ORF-frame's codon usage to the
    average codon usage pooled across all 3 forward frames of the same
    sequence (its own background). Real coding sequence in the correct frame
    typically shows a markedly more skewed codon distribution than the
    sequence's other, non-coding frames."""
    seq = dna(seq)
    if len(seq) < 30:
        return 0.0

    frame_codon_counts = [Counter(_codons(seq, f)) for f in range(3)]
    frame_totals = [sum(c.values()) for c in frame_codon_counts]

    best_frame = max(range(3), key=lambda f: frame_totals[f] and _skew(frame_codon_counts[f]))
    pooled = Counter()
    for f in range(3):
        pooled.update(frame_codon_counts[f])
    pooled_total = sum(pooled.values())
    if pooled_total == 0:
        return 0.0

    chisq = 0.0
    observed = frame_codon_counts[best_frame]
    obs_total = frame_totals[best_frame] or 1
    for codon, obs_count in observed.items():
        expected_freq = pooled[codon] / pooled_total
        expected_count = expected_freq * obs_total
        if expected_count > 0:
            chisq += (obs_count - expected_count) ** 2 / expected_count
    return chisq


def _skew(counter: Counter) -> float:
    total = sum(counter.values())
    if total == 0:
        return 0.0
    probs = [c / total for c in counter.values()]
    return -sum(p * math.log2(p) for p in probs if p > 0)  # entropy; lower = more skewed


@dataclass
class OrphanScore:
    longest_orf_nt: int
    orf_coverage: float       # longest_orf / len(seq), capped at 1
    codon_bias_chisq: float
    orphan_score: float       # higher = more "orphan"/non-coding-like


def orphan_score(seq: str) -> OrphanScore:
    n = len(seq)
    orf_len = longest_orf_nt(seq)
    coverage = min(1.0, orf_len / n) if n else 0.0
    chisq = codon_bias_chisq(seq)
    # Combine: low ORF coverage AND low codon bias -> high orphan score.
    # chisq is unbounded; squash with a simple saturating transform so it
    # contributes on a comparable [0,1] scale.
    chisq_component = 1.0 - (chisq / (chisq + 50.0))
    score = 0.5 * (1.0 - coverage) + 0.5 * chisq_component
    return OrphanScore(
        longest_orf_nt=orf_len,
        orf_coverage=coverage,
        codon_bias_chisq=chisq,
        orphan_score=score,
    )
