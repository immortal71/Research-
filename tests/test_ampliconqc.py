"""Tests for reference-free amplicon detection.

The point of this module is to catch a targeted library even when the
archive metadata is missing or wrong, so the tests contrast a synthetic
shotgun read set against a synthetic primed one, and pin the real
SRR13291825 numbers that motivated the threshold.
"""
import gzip
import random

from rnasig.ampliconqc import profile_fastq, profile_prefixes


def _shotgun_reads(n=2000, length=150, seed=0):
    """Reads starting at arbitrary positions -- 5' prefixes should be flat."""
    rng = random.Random(seed)
    return ["".join(rng.choice("ACGT") for _ in range(length)) for _ in range(n)]


def _amplicon_reads(n=2000, length=150, seed=1):
    """Reads built from three fixed primers, as a PCR library would be."""
    rng = random.Random(seed)
    primers = ["AGCTCCAATAGCGTATATTA", "AGAAGACATCCTTGGTGAAT", "TGAAAACATCCTTGGCAAAT"]
    reads = []
    for i in range(n):
        primer = primers[i % len(primers)]
        body = "".join(rng.choice("ACGT") for _ in range(length - len(primer)))
        reads.append(primer + body)
    return reads


def test_shotgun_library_is_not_flagged():
    profile = profile_prefixes(_shotgun_reads())
    assert not profile.is_amplicon
    assert profile.top_fraction < 0.01
    assert profile.n_reads == 2000


def test_primer_dominated_library_is_flagged():
    profile = profile_prefixes(_amplicon_reads())
    assert profile.is_amplicon
    assert profile.top_fraction > 0.99
    assert len(profile.top_prefixes) == 3


def test_srr13291825_observed_fractions_cross_the_threshold():
    """The three real 5'-20mers from SRR13291825 cover 68.5% of R1.

    Rebuilt here at the observed proportions rather than shipping the
    6.6 MB FASTQ. Counts are from the ENA copy of the run.
    """
    observed = {
        "AGCTCCAATAGCGTATATTA": 16324,
        "AGAAGACATCCTTGGTGAAT": 11542,
        "TGAAAACATCCTTGGCAAAT": 11104,
    }
    total = 56865
    reads = []
    for prefix, count in observed.items():
        reads.extend([prefix + "A" * 130] * count)
    rng = random.Random(7)
    for _ in range(total - sum(observed.values())):
        reads.append("".join(rng.choice("ACGT") for _ in range(150)))

    profile = profile_prefixes(reads)
    assert profile.is_amplicon
    assert 0.68 < profile.top_fraction < 0.69


def test_mixed_library_below_threshold_is_not_flagged():
    reads = _shotgun_reads(n=1900) + _amplicon_reads(n=100)
    profile = profile_prefixes(reads)
    assert not profile.is_amplicon


def test_reads_shorter_than_prefix_are_skipped():
    profile = profile_prefixes(["ACGT", "AC"], prefix_len=20)
    assert profile.n_reads == 0
    assert not profile.is_amplicon
    assert "no reads" in profile.describe()


def test_profile_fastq_reads_plain_and_gzipped(tmp_path):
    reads = _amplicon_reads(n=60)
    lines = []
    for i, seq in enumerate(reads):
        lines += [f"@r{i}", seq, "+", "I" * len(seq)]
    payload = "\n".join(lines) + "\n"

    plain = tmp_path / "reads.fastq"
    plain.write_text(payload)
    assert profile_fastq(str(plain)).is_amplicon

    gzipped = tmp_path / "reads.fastq.gz"
    with gzip.open(gzipped, "wt") as fh:
        fh.write(payload)
    gz_profile = profile_fastq(str(gzipped))
    assert gz_profile.is_amplicon
    assert gz_profile.n_reads == 60


def test_describe_mentions_verdict_and_prefixes():
    text = profile_prefixes(_amplicon_reads(n=300)).describe()
    assert "amplicon-like" in text
    assert "AGCTCCAATAGCGTATATTA" in text
