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


def test_keeps_legitimate_sequence():
    records = [
        Record("real1", "CGATCGATCGATCGATCGATCGATCGATCGAT" * 5),
        Record("real2", "GAGTAGTAGTCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"),
    ]
    kept, report = filter_contigs(records)
    assert len(kept) == 2
    assert report.n_dropped_adapter == 0
    assert report.n_dropped_homopolymer == 0
