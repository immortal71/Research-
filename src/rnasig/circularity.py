"""Circularity detection: reimplementation of the core idea behind VNom's
CircleFinder (Zheludev et al. 2024) using pure Python.

Rationale (from VNom's own README/help text): a circular RNA molecule that
gets sequenced and assembled through its origin produces a *linear* contig
whose 3' end duplicates a stretch of its 5' end (the assembler read through
the junction once or more before its overlap graph closed the loop). So:
apparent circularity <=> the contig's tail is a repeat of some earlier
stretch of the same contig.

We detect this via exact terminal k-mer matching (a fast, alignment-free
proxy for VNom's terminal-repeat scan), then resolve the "monomer" (the
de-duplicated, one-copy-per-unit sequence) and flag concatemers (>=2 tandem
copies).

This is an independent reimplementation, not a port of VNom.py -- see
data/reference/vnom/PROVENANCE.md.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class CircularityResult:
    is_circular: bool
    overhang: int = 0          # length of the duplicated (repeated) tail
    unit_length: int | None = None   # length of one monomer copy
    n_copies: float | None = None    # len(seq) / unit_length
    monomer: str | None = None       # de-duplicated single-copy sequence


def find_circularity(seq: str, k: int = 12, max_mismatch: int = 0) -> CircularityResult:
    """Look for a terminal repeat: does seq end with a k-mer that also occurs
    earlier in the sequence? If several candidate offsets exist, prefer the
    one closest to the 3' end (smallest overhang), matching the typical
    assembler behaviour of stopping as soon as the graph closes.

    max_mismatch allows the terminal k-mer to differ from its earlier copy.
    Exact matching assumes the two passes through the origin are identical,
    which holds for a stable template and fails for a mutable one. On the
    VNom test set, NODE_43354 and NODE_43485 are both Peach latent mosaic
    viroid and both carry a terminal repeat with exactly one substitution,
    implying a 338 nt unit; exact matching misses both. PLMVd replicates as
    a quasispecies, so within-contig variation is the expected case rather
    than an error. See results/phase5_identification_report.md.

    The cost is a higher false-positive rate, since a mismatch-tolerant
    k-mer is a weaker constraint. Measure it with
    circularity_false_positive_rate before raising this on a new dataset.
    """
    n = len(seq)
    if n < 2 * k:
        return CircularityResult(is_circular=False)

    tail = seq[-k:]
    best_i = None
    if max_mismatch <= 0:
        # search for earlier occurrences of the tail k-mer, excluding the
        # trivial position at the very end itself
        start = 0
        while True:
            i = seq.find(tail, start, n - 1)
            if i == -1:
                break
            best_i = i  # keep overwriting -> last (rightmost, closest to 3') match wins
            start = i + 1
    else:
        # Rightmost position whose k-mer is within max_mismatch of the tail.
        # Scanned from the 3' end so the first acceptable hit is the closest
        # one, matching the exact-match branch's preference.
        for i in range(n - k - 1, -1, -1):
            mismatches = 0
            for a, b in zip(tail, seq[i : i + k]):
                if a != b:
                    mismatches += 1
                    if mismatches > max_mismatch:
                        break
            else:
                best_i = i
                break

    if best_i is None:
        return CircularityResult(is_circular=False)

    overhang = best_i + k
    unit_length = n - overhang
    if unit_length < k:
        return CircularityResult(is_circular=False)

    monomer = seq[:unit_length]
    n_copies = n / unit_length

    return CircularityResult(
        is_circular=True,
        overhang=overhang,
        unit_length=unit_length,
        n_copies=n_copies,
        monomer=monomer,
    )


def resolve_concatemer(result: CircularityResult, tandem_thr: float = 0.1) -> CircularityResult:
    """If a circular contig is actually ~N tandem copies of a shorter unit
    (a "concatemer"), collapse it to the single monomer. VNom calls this
    tandem resolution; here we just check whether n_copies is close to an
    integer >= 2 within tandem_thr fractional tolerance."""
    if not result.is_circular or result.n_copies is None:
        return result
    nearest_int = round(result.n_copies)
    if nearest_int >= 2 and abs(result.n_copies - nearest_int) / nearest_int <= tandem_thr:
        true_unit_len = len(result.monomer) if result.monomer else None
        # unit_length already reflects one copy; nothing further to do,
        # but flag copy number as the rounded integer for downstream reporting
        return CircularityResult(
            is_circular=True,
            overhang=result.overhang,
            unit_length=result.unit_length,
            n_copies=float(nearest_int),
            monomer=result.monomer,
        )
    return result


def circularity_false_positive_rate(null_seqs: list[str], k: int = 12, max_mismatch: int = 0) -> float:
    """Fraction of a null (non-circular, e.g. shuffled) sequence set flagged
    as circular -- used to calibrate k and max_mismatch and report the
    empirical FPR of the detector itself."""
    if not null_seqs:
        return 0.0
    hits = sum(1 for s in null_seqs if find_circularity(s, k=k, max_mismatch=max_mismatch).is_circular)
    return hits / len(null_seqs)
