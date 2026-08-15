"""Rod-likeness: does a fold look like a viroid, or like one hairpin?

The structure z-score asks whether a sequence folds better than its own
shuffle. That turns out to be necessary but not sufficient, and the gap
showed up the first time the Phase 6 sweep produced candidates.

SRR5949183_contig_204 scored z=5.05, comfortably over the shortlist bar,
on the strength of a single 17 bp hairpin with only 44% of its bases paired
and the rest dangling. That is what a transcription terminator or a
transposon's terminal inverted repeat looks like, and there are a great many
of both in any metatranscriptome. Peach latent mosaic viroid, recovered from
real reads by this same pipeline, is 68% paired with an MFE of -0.473
kcal/mol per nucleotide. Its 943 nt concatemer gives the same two numbers,
66% and -0.473, so they describe the molecule rather than the contig.

Viroids and obelisks are rods: base-paired along essentially their whole
length rather than locally. These two measures capture that and a z-score
does not, because a shuffle of a sequence containing one strong inverted
repeat also folds poorly, which makes the z look impressive.

Thresholds are set from the PLMVd reference with margin, not tuned to make
anything in particular pass.
"""
from __future__ import annotations

from dataclasses import dataclass

from .seqio import rna

# Calibrated against nine described viroids (data/reference/viroids), not
# against PLMVd alone. The first version of this module set the cutoffs from
# PLMVd only and was wrong in a way worth recording: Avocado sunblotch
# viroid is AT-rich (GC 0.38) and folds to -0.294 kcal/mol/nt, so a -0.30
# cutoff rejected a real viroid.
#
# Measured across the nine: paired 65-74%, MFE/nt -0.294 to -0.482. Both
# cutoffs sit outside that range with margin.
MIN_PAIRED_FRACTION = 0.60
MAX_MFE_PER_NT = -0.25

# Distinct 6-mers as a fraction of positions. This replaces the structure
# z-score as the guard against low-complexity sequence, and it does that job
# far better.
#
# z was doing two things at once: rewarding real structure and rejecting
# repeats. It is unreliable at the first, because viroid base composition is
# itself self-complementary, so a dinucleotide shuffle folds about as well as
# the molecule. Pear blister canker viroid scores z = -0.02: a genuine viroid
# that folds no better than its own shuffle. Apple scar skin viroid scores
# 3.06. A z >= 3 gate rejected three of seven real viroids.
#
# Complexity separates the two cases cleanly and without ambiguity. The nine
# viroids run 0.90-0.98; poly-AT, poly-GC, (ACGT)n and (AAGCT)n all sit at
# 0.01. Random sequence scores 0.95 on complexity but is caught by paired
# fraction at 58%, so the two measures cover each other.
MIN_COMPLEXITY = 0.50
COMPLEXITY_K = 6


@dataclass
class RodProfile:
    length: int
    mfe: float
    mfe_per_nt: float
    paired_fraction: float
    longest_helix: int
    complexity: float = 1.0
    structure: str = ""

    @property
    def is_rodlike(self) -> bool:
        """Rod-shaped, thermodynamically stable, and not a repeat."""
        return (
            self.paired_fraction >= MIN_PAIRED_FRACTION
            and self.mfe_per_nt <= MAX_MFE_PER_NT
            and self.complexity >= MIN_COMPLEXITY
        )

    def describe(self) -> str:
        verdict = "rod-like" if self.is_rodlike else "not rod-like"
        return (
            f"len={self.length} MFE={self.mfe:.1f} ({self.mfe_per_nt:.3f}/nt) "
            f"paired={self.paired_fraction:.0%} complexity={self.complexity:.2f} "
            f"longest_helix={self.longest_helix} -> {verdict}"
        )


def sequence_complexity(seq: str, k: int = COMPLEXITY_K) -> float:
    """Distinct k-mers as a fraction of positions. 1.0 is fully non-repetitive."""
    if len(seq) < k:
        return 0.0
    positions = len(seq) - k + 1
    return len({seq[i : i + k] for i in range(positions)}) / positions


def rod_profile(seq: str) -> RodProfile:
    """Fold a sequence and measure how rod-like the result is."""
    import RNA  # imported lazily so the module can be inspected without ViennaRNA

    structure, mfe = RNA.fold(rna(seq))
    paired = structure.count("(") * 2
    longest = 0
    run = 0
    for char in structure:
        if char == "(":
            run += 1
            longest = max(longest, run)
        else:
            run = 0
    n = max(len(seq), 1)
    return RodProfile(
        length=len(seq),
        mfe=mfe,
        mfe_per_nt=mfe / n,
        paired_fraction=paired / n,
        longest_helix=longest,
        complexity=sequence_complexity(seq),
        structure=structure,
    )
