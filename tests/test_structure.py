import random

from rnasig.structure import structure_zscore
from rnasig.simulate import make_rod
from rnasig.nullmodel import markov_background


def test_designed_rod_scores_higher_than_random_background():
    rng = random.Random(9)
    rod = make_rod(rng, n_hairpins=6, stem_len=14, loop_len=5)
    background = markov_background(len(rod), gc_content=0.5, rng=rng)

    rod_stability = structure_zscore(rod, n_shuffles=30, rng=rng)
    bg_stability = structure_zscore(background, n_shuffles=30, rng=rng)

    assert rod_stability.z_score > bg_stability.z_score


def test_zscore_is_finite_and_reasonable():
    rng = random.Random(10)
    seq = make_rod(rng, n_hairpins=4)
    result = structure_zscore(seq, n_shuffles=20, rng=rng)
    assert -20 < result.z_score < 20
    assert result.n_shuffles == 20
