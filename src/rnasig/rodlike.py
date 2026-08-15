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

# From PLMVd measured twice (337 nt monomer and 943 nt concatemer): 66-68%
# paired, -0.473 kcal/mol/nt. Both cutoffs sit well below that so a less
# extreme replicon still passes, while a lone hairpin does not.
MIN_PAIRED_FRACTION = 0.55
MAX_MFE_PER_NT = -0.30


@dataclass
class RodProfile:
    length: int
    mfe: float
    mfe_per_nt: float
    paired_fraction: float
    longest_helix: int
    structure: str = ""

    @property
    def is_rodlike(self) -> bool:
        return (
            self.paired_fraction >= MIN_PAIRED_FRACTION
            and self.mfe_per_nt <= MAX_MFE_PER_NT
        )

    def describe(self) -> str:
        verdict = "rod-like" if self.is_rodlike else "not rod-like"
        return (
            f"len={self.length} MFE={self.mfe:.1f} ({self.mfe_per_nt:.3f}/nt) "
            f"paired={self.paired_fraction:.0%} longest_helix={self.longest_helix} -> {verdict}"
        )


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
        structure=structure,
    )
