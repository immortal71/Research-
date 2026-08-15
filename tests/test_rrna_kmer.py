"""Tests for reference-based rRNA exclusion.

The case that motivated this module: a 342 nt contig from oral
metatranscriptome SRR19432462, assembled at 6,752x coverage, which cleared
both the structure z bar and the rod-likeness bar and turned out on blastn
to be 28S rRNA at 100% identity. Ribosomal RNA is structured, abundant and
assembles with terminal repeats, so it clears every axis of this signature
at once.

No network is used here: the reference is a short synthetic stand-in, and
what is being tested is the screening logic rather than NCBI.
"""
from rnasig.rrna_kmer import (
    DEFAULT_K,
    build_reference_kmers,
    rrna_kmer_fraction,
    screen_records,
)
from rnasig.seqio import Record, revcomp

_REF_SEQ = (
    "GGTTAAGCGACTAAGCGTACACGGTGGATGCCCTGGCAGTCAGAGGCGATGAAGGACGTGCTAATCTGCG"
    "ATAAGCGTCGGTAAGGTGATATGAACCGTTATAACCGGCGATTTCCGAATGGGGAAACCCAGTGTGTTTC"
    "GACACACTATCATTAACTGAATCCATAGGTTAATGAGGCGAACCGGGGGAACTGAAACATCTAAGTACCC"
)
_REFERENCE = f">ref_rrna\n{_REF_SEQ}\n"


def test_reference_kmers_cover_both_strands():
    kmers = build_reference_kmers(_REFERENCE, k=DEFAULT_K)
    assert _REF_SEQ[:DEFAULT_K] in kmers
    assert revcomp(_REF_SEQ[:DEFAULT_K]) in kmers


def test_a_reference_fragment_scores_one():
    ref = build_reference_kmers(_REFERENCE)
    assert rrna_kmer_fraction(_REF_SEQ[30:180], ref) == 1.0


def test_reverse_complement_fragment_is_also_caught():
    """An rRNA contig assembled on the opposite strand is still rRNA."""
    ref = build_reference_kmers(_REFERENCE)
    assert rrna_kmer_fraction(revcomp(_REF_SEQ[30:180]), ref) == 1.0


def test_unrelated_sequence_scores_zero():
    ref = build_reference_kmers(_REFERENCE)
    unrelated = "ACGTTGCAATCGGATCCTAGGCATTACGGCATTACGGCATTACGGCATTACGGCATTACGGCATTACGGC"
    assert rrna_kmer_fraction(unrelated, ref) == 0.0


def test_screen_drops_rrna_and_keeps_the_rest():
    ref = build_reference_kmers(_REFERENCE)
    records = [
        Record("rrna_contig", _REF_SEQ[20:200]),
        Record("real_candidate", "ACGTTGCAATCGGATCCTAGG" * 8),
    ]
    kept, report = screen_records(records, ref)
    assert [r.id for r in kept] == ["real_candidate"]
    assert report.n_dropped == 1
    assert report.dropped["rrna_contig"] == 1.0
    assert "1 dropped" in report.describe()


def test_partial_rrna_still_exceeds_the_threshold():
    """A contig that is only part rRNA must still be dropped."""
    ref = build_reference_kmers(_REFERENCE)
    chimeric = "ACGTTGCAATCGGATCCTAGGCATTACGGCATTACGG" * 3 + _REF_SEQ[:120]
    fraction = rrna_kmer_fraction(chimeric, ref)
    assert fraction >= 0.10


def test_missing_reference_fails_open():
    """With no reference the screen must pass everything, not drop everything."""
    records = [Record("a", "ACGT" * 40), Record("b", "TTGC" * 40)]
    kept, report = screen_records(records, set())
    assert len(kept) == 2
    assert report.n_dropped == 0
    assert not report.reference_loaded
    assert "skipped" in report.describe()


def test_sequences_shorter_than_k_score_zero():
    ref = build_reference_kmers(_REFERENCE)
    assert rrna_kmer_fraction("ACGT", ref) == 0.0
