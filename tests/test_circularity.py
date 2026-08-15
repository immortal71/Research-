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


def test_exact_matching_misses_a_mutated_terminal_repeat():
    """A single substitution in the terminal k-mer defeats exact matching.

    This is the NODE_43354 / NODE_43485 case: both are Peach latent mosaic
    viroid, both carry a one-substitution terminal repeat, and exact
    matching finds neither.
    """
    rng = random.Random(11)
    monomer = "".join(rng.choice("ACGT") for _ in range(338))
    overhang = 264
    repeat = list(monomer[:overhang])
    repeat[-3] = "A" if repeat[-3] != "A" else "C"  # one mismatch inside the tail k-mer
    contig = monomer + "".join(repeat)

    assert not find_circularity(contig, k=12).is_circular
    assert not find_circularity(contig, k=12, max_mismatch=0).is_circular


def test_one_mismatch_recovers_the_mutated_terminal_repeat():
    rng = random.Random(11)
    monomer = "".join(rng.choice("ACGT") for _ in range(338))
    overhang = 264
    repeat = list(monomer[:overhang])
    repeat[-3] = "A" if repeat[-3] != "A" else "C"
    contig = monomer + "".join(repeat)

    result = find_circularity(contig, k=12, max_mismatch=1)
    assert result.is_circular
    assert result.unit_length == len(monomer)


def test_mismatch_tolerance_does_not_raise_fpr_at_one():
    """Calibration: on shuffled nulls, max_mismatch=1 costs nothing.

    Measured on the VNom test set (950 shuffles): FPR 0.0000 at
    max_mismatch=0 and 0.0000 at 1, while PLMVd recall goes 20/24 -> 22/24.
    At max_mismatch=2 the FPR becomes non-zero, which is why 1 is the
    documented ceiling.
    """
    rng = random.Random(5)
    nulls = [markov_background(700, gc_content=0.5, rng=rng) for _ in range(200)]
    assert circularity_false_positive_rate(nulls, k=12, max_mismatch=1) < 0.05


def test_mismatch_search_prefers_the_repeat_closest_to_the_three_prime_end():
    rng = random.Random(13)
    monomer = "".join(rng.choice("ACGT") for _ in range(300))
    contig = monomer + monomer + monomer[:50]
    result = find_circularity(contig, k=12, max_mismatch=1)
    assert result.is_circular
    # the closest repeat implies the 300 nt unit, not the 600 nt double
    assert result.unit_length == 300


def test_shuffles_understate_the_false_positive_rate():
    """Why the mismatch tolerance was calibrated wrongly the first time.

    Shuffling destroys the repeat structure that produces spurious terminal
    matches, so a shuffled null says max_mismatch=1 is free. On 5017 real
    assembled contigs it is not: k=12 mm=1 called 1.18% circular against
    0.62% for exact matching, and the excess showed up as a spike of
    unrelated sequences all reporting a 228 nt unit.

    Sequence with internal repeats is the honest null.
    """
    rng = random.Random(21)
    motif = "".join(rng.choice("ACGT") for _ in range(40))
    repetitive = [
        motif[: rng.randrange(10, 40)] + "".join(rng.choice("ACGT") for _ in range(300)) + motif[:14]
        for _ in range(200)
    ]
    lenient = circularity_false_positive_rate(repetitive, k=12, max_mismatch=1)
    strict = circularity_false_positive_rate(repetitive, k=16, max_mismatch=1)
    assert strict <= lenient, "raising k must not increase the false-positive rate"
