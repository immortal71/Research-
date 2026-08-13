from rnasig.preprocess import filter_contigs
from rnasig.seqio import Record


def test_flags_illumina_p7_adapter_contig():
    # This is the exact top hit that was flagged from real SRA run SRR13291825
    # -- the assembled Illumina adapter dimer that motivated this filter.
    adapter_contig = (
        "TTTTTTTTTTTTCAAGCAGAAGACGGCATACGAGATTGAGCTAGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTA"
        "GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
    )
    good = "ACGTACGTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
    records = [Record("adapter", adapter_contig), Record("real", good)]
    kept, report = filter_contigs(records)
    assert [r.id for r in kept] == ["real"]
    assert report.n_dropped_adapter == 1
    assert "adapter" in report.dropped_reasons["adapter"]


def test_flags_homopolymer_contig():
    records = [
        Record("polyA", "A" * 30 + "CGTACGTACGTACGTACGT"),
        Record("clean", "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"),
    ]
    kept, report = filter_contigs(records)
    assert [r.id for r in kept] == ["clean"]
    assert report.n_dropped_homopolymer == 1


def test_flags_truncated_adapter_readthrough_at_contig_end():
    """k141_3 from the SRR13291825 sweep, verbatim.

    It ends in AGATCGGAAG -- 10 nt of the universal Illumina adapter, cut
    short where the assembler ran out of coverage. The 20-nt interior rule
    cannot see that, so this contig was ranked as a novel candidate twice
    before the terminal check existed. Its 3' junction is present in 894
    raw reads, so it is read-through, not a coincidence.
    """
    k141_3 = (
        "CACCCGCCGCGCTTCCGACTGTATAACCTTCTTTCGTTCTTGATCACATTGAATCGGCGGACGTCGTCGG"
        "GGCCCATCGGGTGGTCGGCGGGCGCGTGCGTCATCACCACGGCGTCGATCATGTCGACGCCGTTGACCAC"
        "GCACTGCAGGCGAAGCTCGGGCGTGGTGTCAACCAGGACACGGGTTTCGCCCTGGCTGATGAGCACCGAC"
        "GGGCGGGTGCGTTTGTCATGCGAGTCGGTGCTGCGGCACACCTCGCAGTGGCAGCCGATCATGGGAATTA"
        "CCGCGGCTGCTGGAGTAGAGAAGATCGGAAG"
    )
    kept, report = filter_contigs([Record("k141_3", k141_3)])
    assert kept == []
    assert report.n_dropped_adapter == 1
    assert "adapter_terminal" in report.dropped_reasons["k141_3"]
    assert "AGATCGGAAG" in report.dropped_reasons["k141_3"]


def test_flags_reverse_complement_adapter_at_contig_start():
    rc_tail = "GCTCTTCCGATCT"  # 3' end of the TruSeq Read 1 primer
    contig = rc_tail + "ACGTGGCATTACGATCGGCATTACGGCATTACGGCATTACGGCAT"
    kept, report = filter_contigs([Record("rc_lead", contig)])
    assert kept == []
    assert report.n_dropped_adapter == 1


def test_terminal_check_does_not_strip_short_incidental_matches():
    """Below the 10-nt terminal floor, sequence is left alone."""
    contig = "GCATTACGGCATTACGGCATTACGGCATTACGGCATTACGGCATTACGG" + "AGATCGG"
    kept, report = filter_contigs([Record("incidental", contig)])
    assert [r.id for r in kept] == ["incidental"]
    assert report.n_dropped_adapter == 0


def test_keeps_legitimate_sequence():
    records = [
        Record("real1", "CGATCGATCGATCGATCGATCGATCGATCGAT" * 5),
        Record("real2", "GAGTAGTAGTCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"),
    ]
    kept, report = filter_contigs(records)
    assert len(kept) == 2
    assert report.n_dropped_adapter == 0
    assert report.n_dropped_homopolymer == 0
