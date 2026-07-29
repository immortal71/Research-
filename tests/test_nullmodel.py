import random
from collections import Counter

from rnasig.nullmodel import dinuc_shuffle, markov_background


def _dinuc_counts(seq: str) -> Counter:
    return Counter(seq[i : i + 2] for i in range(len(seq) - 1))


def test_dinuc_shuffle_preserves_dinucleotide_composition_exactly():
    rng = random.Random(42)
    for trial in range(20):
        length = rng.randint(20, 500)
        seq = "".join(rng.choice("ACGT") for _ in range(length))
        shuffled = dinuc_shuffle(seq, rng=rng)
        assert len(shuffled) == len(seq)
        assert Counter(seq) == Counter(shuffled)
        assert _dinuc_counts(seq) == _dinuc_counts(shuffled)


def test_dinuc_shuffle_handles_highly_repetitive_sequences():
    # regression: naive bucket shuffling can strand the reconstruction walk
    # on low-diversity/repetitive sequences unless the last-edge-per-node
    # is preserved (Altschul-Erikson).
    rng = random.Random(7)
    seq = "AT" * 100 + "GC" * 50
    for _ in range(50):
        shuffled = dinuc_shuffle(seq, rng=rng)
        assert Counter(seq) == Counter(shuffled)


def test_dinuc_shuffle_actually_permutes():
    rng = random.Random(3)
    seq = "".join(rng.choice("ACGT") for _ in range(300))
    shuffled = dinuc_shuffle(seq, rng=rng)
    assert shuffled != seq


def test_markov_background_matches_target_gc():
    rng = random.Random(5)
    seq = markov_background(5000, gc_content=0.65, rng=rng)
    gc = (seq.count("G") + seq.count("C")) / len(seq)
    assert abs(gc - 0.65) < 0.03
