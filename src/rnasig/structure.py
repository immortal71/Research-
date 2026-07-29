"""Axis S1 of the novel signature: is this sequence's predicted secondary
structure implausibly stable for its base composition?

We fold with ViennaRNA's RNAfold (MFE), then compare against the MFE
distribution of dinucleotide-composition-matched shuffles of the same
sequence. This is the standard single-sequence "thermodynamic z-score"
approach for flagging structured non-coding RNA candidates without needing
a comparative alignment (Workman & Krogh 1999; Rivas & Eddy 2000; the same
logic underlies RNAz's z-score component, but computed per-sequence rather
than per-alignment column, since we have no reference genome to align
against for a genuinely novel element).
"""
from __future__ import annotations

import random
import statistics
from dataclasses import dataclass

import RNA

from .nullmodel import dinuc_shuffle
from .seqio import rna


def mfe(seq: str) -> tuple[str, float]:
    """Minimum free energy structure and energy (kcal/mol) for `seq`."""
    structure, energy = RNA.fold(rna(seq))
    return structure, energy


@dataclass
class StructureStability:
    energy: float
    structure: str
    null_mean: float
    null_std: float
    z_score: float
    n_shuffles: int


def structure_zscore(
    seq: str,
    n_shuffles: int = 100,
    rng: random.Random | None = None,
) -> StructureStability:
    """z-score of structural stability relative to dinucleotide-shuffled
    nulls. Positive z = more stable (more negative MFE) than expected by
    chance given the sequence's own composition -- i.e. "the fold isn't
    just an artifact of base content."""
    rng = rng or random.Random()
    structure, energy = mfe(seq)

    null_energies = []
    for _ in range(n_shuffles):
        shuffled = dinuc_shuffle(seq, rng=rng)
        _, e = mfe(shuffled)
        null_energies.append(e)

    null_mean = statistics.fmean(null_energies)
    null_std = statistics.pstdev(null_energies) or 1e-9
    z = (null_mean - energy) / null_std  # more negative real energy -> larger positive z

    return StructureStability(
        energy=energy,
        structure=structure,
        null_mean=null_mean,
        null_std=null_std,
        z_score=z,
        n_shuffles=n_shuffles,
    )
