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


def test_low_complexity_repeats_are_rejected_by_complexity():
    """A poly-AT repeat folds into one perfect helix and must still fail.

    98% paired at -0.54 kcal/mol/nt passes both shape thresholds. The
    structure z-score used to be what caught it, but z turned out to reject
    real viroids too, so complexity does this job now: repeats sit at 0.01
    while the nine described viroids run 0.90-0.98.
    """
    for repeat in ("AT" * 80, "GC" * 80, "ACGT" * 40, "AAGCT" * 32):
        profile = rod_profile(repeat)
        assert profile.complexity < 0.10
        assert not profile.is_rodlike, f"{repeat[:10]}... slipped through"


def test_every_described_viroid_passes():
    """The calibration set. A filter that rejects real viroids is broken.

    z >= 3 rejected three of these, which is why it is no longer a gate.
    """
    import os

    path = os.path.join("data", "reference", "viroids", "described_viroids.fasta")
    if not os.path.exists(path):
        pytest.skip("viroid reference not present")
    viroids = [r for r in read_fasta(path) if "viroid" in r.id.lower()]
    assert len(viroids) >= 8
    failed = [r.id.split()[0] for r in viroids if not rod_profile(r.seq).is_rodlike]
    assert not failed, f"filters reject real viroids: {failed}"


def test_avocado_sunblotch_the_at_rich_case_passes():
    """The viroid that broke the PLMVd-only calibration.

    ASBVd is GC 0.38 and folds to -0.294 kcal/mol/nt, so the original
    -0.30 cutoff rejected it.
    """
    import os

    path = os.path.join("data", "reference", "viroids", "described_viroids.fasta")
    if not os.path.exists(path):
        pytest.skip("viroid reference not present")
    asbvd = [r for r in read_fasta(path) if "sunblotch" in r.id.lower()]
    if not asbvd:
        pytest.skip("ASBVd not in reference")
    profile = rod_profile(asbvd[0].seq)
    assert profile.is_rodlike
    assert profile.mfe_per_nt > MAX_MFE_PER_NT - 0.06


def test_recovered_plmvd_matches_the_reference_profile(tmp_path):
    """The pipeline's own output should reproduce the reference numbers."""
    fasta = tmp_path / "p.fasta"
    fasta.write_text(f">plmvd\n{PLMVD}\n")
    profile = rod_profile(read_fasta(str(fasta))[0].seq)
    assert 0.60 <= profile.paired_fraction <= 0.75
    assert -0.55 <= profile.mfe_per_nt <= -0.40
