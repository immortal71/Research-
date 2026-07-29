import random

from rnasig.orphan import orphan_score
from rnasig.simulate import make_orf_like, make_rod


def test_orf_like_sequence_scores_low_orphan():
    rng = random.Random(11)
    coding = make_orf_like(rng, n_codons=200)
    result = orphan_score(coding)
    assert result.orf_coverage > 0.8
    assert result.orphan_score < 0.4


def test_structured_noncoding_sequence_scores_higher_orphan():
    rng = random.Random(12)
    noncoding = make_rod(rng, n_hairpins=6)
    result = orphan_score(noncoding)
    assert result.orphan_score > 0.4
