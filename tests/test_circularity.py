import random

from rnasig.circularity import find_circularity, circularity_false_positive_rate
from rnasig.nullmodel import markov_background


def test_detects_planted_circularity():
    rng = random.Random(1)
    monomer = "".join(rng.choice("ACGT") for _ in range(500))
    overhang = 40
    contig = monomer + monomer[:overhang]
    result = find_circularity(contig, k=12)
    assert result.is_circular
    assert result.unit_length == len(monomer)
    assert result.monomer == monomer


def test_rejects_linear_sequence():
    rng = random.Random(2)
    seq = "".join(rng.choice("ACGT") for _ in range(600))
    result = find_circularity(seq, k=12)
    assert not result.is_circular


def test_false_positive_rate_is_low_on_random_backgrounds():
    rng = random.Random(3)
    nulls = [markov_background(700, gc_content=0.5, rng=rng) for _ in range(200)]
    fpr = circularity_false_positive_rate(nulls, k=12)
    assert fpr < 0.05


def test_short_sequences_never_flagged():
    result = find_circularity("ACGTACGTACGT", k=12)
    assert not result.is_circular
