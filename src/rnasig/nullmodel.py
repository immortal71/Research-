"""Composition-matched null sequence generators used throughout calibration."""
from __future__ import annotations

import random


def dinuc_shuffle(seq: str, rng: random.Random | None = None) -> str:
    """Shuffle `seq` while exactly preserving its dinucleotide (and therefore
    mononucleotide) composition.

    Method: the Altschul-Erikson (1985) algorithm. Bucket each character's
    outgoing "next character" by source character (i.e. group the
    sequence's dinucleotide edges by their first base), shuffle each
    bucket *except* keep each character's last-occurring outgoing edge
    fixed in the last slot, then walk from the original first character
    popping from the shuffled buckets. Naively shuffling every edge
    (including the last one) can strand the walk part-way through
    reconstruction on repetitive/structured sequences -- this is exactly
    the failure mode Altschul-Erikson's construction avoids: preserving
    the true last exit edge of every node guarantees the walk always
    completes, while the dinucleotide count matrix is still preserved
    exactly. This is the standard "doublet-preserving" shuffle used as a
    null model for RNA structure statistics (e.g. Workman & Krogh 1999,
    Rivas & Eddy 2000).
    """
    rng = rng or random
    n = len(seq)
    if n < 2:
        return seq

    buckets: dict[str, list[str]] = {}
    for a, b in zip(seq[:-1], seq[1:]):
        buckets.setdefault(a, []).append(b)
    for c, edges in buckets.items():
        if len(edges) > 1:
            last = edges[-1]
            rest = edges[:-1]
            rng.shuffle(rest)
            buckets[c] = rest + [last]

    cursors = {k: 0 for k in buckets}
    out = [seq[0]]
    current = seq[0]
    for _ in range(n - 1):
        bucket = buckets[current]
        idx = cursors[current]
        nxt = bucket[idx]
        cursors[current] += 1
        out.append(nxt)
        current = nxt
    return "".join(out)


def mononuc_shuffle(seq: str, rng: random.Random | None = None) -> str:
    """Simple order-0 shuffle: preserves base composition only, destroys all
    positional/dinucleotide structure. Useful as a stricter/weaker null for
    sensitivity comparisons."""
    rng = rng or random
    chars = list(seq)
    rng.shuffle(chars)
    return "".join(chars)


_BASES = "ACGT"


def markov_background(
    length: int,
    gc_content: float = 0.5,
    rng: random.Random | None = None,
) -> str:
    """Generate an i.i.d. background sequence at a given GC content. Used as
    generic "surrounding" sequence context for injection experiments -- NOT a
    claim about any real organism's transcriptome, just a neutral canvas."""
    rng = rng or random
    p_gc = gc_content / 2
    p_at = (1 - gc_content) / 2
    weights = [p_at, p_gc, p_gc, p_at]  # A C G T
    return "".join(rng.choices(_BASES, weights=weights, k=length))
