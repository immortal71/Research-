"""Tests for sense/antisense co-occurrence detection.

This axis is half the signature and was missing from phase6 entirely. It
matters most where circularity fails: on Streptococcus sanguinis SK36 the
assembler never closed the obelisk circle, so the circularity test found
nothing, while SAS surfaced a 1011 nt rod-like contig with no nt or
protein hit whose antisense strand sits at 29,165x.
"""
import random

from rnasig.sasfind import find_sas_pairs
from rnasig.seqio import Record, revcomp


def _rand(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _background(n=150, length=600, seed=99):
    rng = random.Random(seed)
    return [
        Record(f"bg{i} multi=8.0", "".join(rng.choice("ACGT") for _ in range(length)))
        for i in range(n)
    ]


def test_finds_a_molecule_present_on_both_strands():
    molecule = _rand(900, 5)
    records = [
        Record("sense multi=120.0", molecule),
        Record("antisense multi=95.0", revcomp(molecule)),
    ] + _background()
    report = find_sas_pairs(records)
    assert len(report.pairs) == 1
    pair = report.pairs[0]
    assert pair.shared_fraction > 0.9
    assert {pair.sense_id, pair.antisense_id} == {"sense", "antisense"}


def test_single_stranded_molecule_does_not_pair():
    records = [Record("solo multi=100.0", _rand(900, 5))] + _background()
    assert find_sas_pairs(records).pairs == []


def test_unrelated_contigs_do_not_pair():
    assert find_sas_pairs(_background(n=200)).pairs == []


def test_partial_antisense_still_pairs():
    """A quasispecies will not give an identical reverse complement."""
    molecule = _rand(900, 7)
    rc = list(revcomp(molecule))
    rng = random.Random(8)
    for _ in range(25):
        pos = rng.randrange(len(rc))
        rc[pos] = rng.choice("ACGT")
    records = [
        Record("sense multi=50.0", molecule),
        Record("antisense multi=40.0", "".join(rc)),
    ] + _background()
    assert len(find_sas_pairs(records).pairs) == 1


def test_short_contigs_are_ignored():
    molecule = _rand(100, 9)
    records = [Record("a", molecule), Record("b", revcomp(molecule))]
    assert find_sas_pairs(records, min_length=150).pairs == []


def test_coverage_is_parsed_from_either_header_style():
    molecule = _rand(600, 11)
    for header_a, header_b in ((("a multi=33.0"), ("b multi=12.0")),
                               (("a cov=33.0"), ("b cov=12.0"))):
        records = [Record(header_a, molecule), Record(header_b, revcomp(molecule))]
        pair = find_sas_pairs(records).pairs[0]
        assert max(pair.sense_coverage, pair.antisense_coverage) == 33.0


def test_pairs_are_ranked_by_coverage():
    lo, hi = _rand(700, 13), _rand(700, 14)
    records = [
        Record("lo_s multi=5.0", lo), Record("lo_a multi=4.0", revcomp(lo)),
        Record("hi_s multi=900.0", hi), Record("hi_a multi=800.0", revcomp(hi)),
    ]
    pairs = find_sas_pairs(records).pairs
    assert len(pairs) == 2
    assert max(pairs[0].sense_coverage, pairs[0].antisense_coverage) == 900.0
