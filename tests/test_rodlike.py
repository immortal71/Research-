"""Tests for rod-likeness, calibrated against a real viroid.

The PLMVd numbers here are measured, not invented: they come from the
337 nt monomer this pipeline recovered from SRR11060618 and from VNom's
943 nt NODE_36652 contig of the same molecule.
"""
import pytest

from rnasig.rodlike import MAX_MFE_PER_NT, MIN_PAIRED_FRACTION, rod_profile
from rnasig.seqio import read_fasta

RNA = pytest.importorskip("RNA", reason="ViennaRNA not installed")

# The 337 nt PLMVd monomer recovered by scripts/phase6_hunt.py.
PLMVD = (
    "AGCGGTCGAACCCAGGGGGAGTGTGACCCAGGTACCGCCGTAGAAACTGGGTTACGACGCCTACCCGGGA"
    "TCGACGTCGATACCATGACGATGCGCTCGTCGGTCCCCGGCCTTTACGTCGCCGGCGGTGTTGGCGGGCA"
    "CAGCAACGGGCTGATTGGGCTTGCGACCTACGGCAAGACATTTCGCCAATTGGGCCTGCAATTCCCTCAG"
    "CAGTTGGTCGAGACCGCGGTCACTTCACACTATCGGCAGGGACAGACCTCGCAGCTCACCACAGCGGTAC"
    "CTGGGTCACACTCCCCCTGGGTTCGACCGCT"
)


def _hairpin_with_tail():
    """One strong hairpin, then unpaired sequence: the contig_204 shape."""
    stem = "GGACGGATACAAATTC"
    loop = "TTTACGAA"
    tail = "TTTATCATATTTTATTTGCTTTATTATCATCCAATCCAAATCTATCATACAACTTTATTGATTTGAATAT"
    from rnasig.seqio import revcomp

    return stem + loop + revcomp(stem) + tail * 3


def test_plmvd_is_rodlike():
    profile = rod_profile(PLMVD)
    assert profile.is_rodlike
    assert profile.paired_fraction > 0.6
    assert profile.mfe_per_nt < -0.4
    assert "rod-like" in profile.describe()


def test_single_hairpin_with_tail_is_not_rodlike():
    profile = rod_profile(_hairpin_with_tail())
    assert not profile.is_rodlike
    assert profile.paired_fraction < MIN_PAIRED_FRACTION


def test_thresholds_sit_below_the_plmvd_reference():
    """Cutoffs must have margin under the real molecule, not hug it."""
    profile = rod_profile(PLMVD)
    assert MIN_PAIRED_FRACTION < profile.paired_fraction
    assert MAX_MFE_PER_NT > profile.mfe_per_nt


def test_profile_reports_the_longest_helix():
    profile = rod_profile(PLMVD)
    assert profile.longest_helix >= 15
    assert profile.structure and len(profile.structure) == len(PLMVD)


def test_low_complexity_repeat_is_rodlike_and_needs_the_z_score_too():
    """Rod-likeness alone is not enough, and this is why.

    A poly-AT repeat folds into one perfect long helix: 98% paired at
    -0.54 kcal/mol/nt, which passes both rod thresholds. It is caught by
    the structure z-score instead, because a shuffle of it folds just as
    well, giving z=0. The two measures are complementary and the pipeline
    requires both.
    """
    from rnasig.structure import structure_zscore

    profile = rod_profile("AT" * 80)
    assert profile.is_rodlike
    assert structure_zscore("AT" * 80, n_shuffles=50).z_score < 1.0


def test_recovered_plmvd_matches_the_reference_profile(tmp_path):
    """The pipeline's own output should reproduce the reference numbers."""
    fasta = tmp_path / "p.fasta"
    fasta.write_text(f">plmvd\n{PLMVD}\n")
    profile = rod_profile(read_fasta(str(fasta))[0].seq)
    assert 0.60 <= profile.paired_fraction <= 0.75
    assert -0.55 <= profile.mfe_per_nt <= -0.40
