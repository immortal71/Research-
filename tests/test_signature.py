from rnasig.signature import benjamini_hochberg


def test_bh_fdr_all_significant_when_all_pvalues_tiny():
    pvals = [1e-6, 1e-5, 1e-4, 1e-6]
    reject, q = benjamini_hochberg(pvals, alpha=0.05)
    assert all(reject)


def test_bh_fdr_none_significant_when_all_pvalues_large():
    pvals = [0.9, 0.8, 0.95, 0.99]
    reject, q = benjamini_hochberg(pvals, alpha=0.05)
    assert not any(reject)


def test_bh_fdr_mixed():
    pvals = [0.001, 0.02, 0.5, 0.7, 0.9]
    reject, q = benjamini_hochberg(pvals, alpha=0.1)
    assert reject[0]  # smallest p-value should survive
    assert not reject[-1]  # largest should not
