"""The novel signature ("SOS": Structure-stable, Orphan-sequence,
Strand-symmetric) and calibration statistics.

Three orthogonal axes, deliberately chosen to be complementary to the
obelisk/VNom signature rather than a re-detection of it:

  S1 -- Structure stability z-score (structure.py): predicted fold is a
        statistical outlier of stability relative to composition-matched
        shuffles ("stable structure").
  S2 -- Orphan/coding-potential score (orphan.py): sequence shows little
        evidence of protein-coding potential ("matches nothing" in the
        reference-free sense available offline -- see orphan.py docstring).
  S3 -- Strand-symmetry / circularity evidence (cluster.py, circularity.py):
        both sense and antisense contigs co-occur in a cluster, and/or the
        contig shows a circular terminal repeat -- evidence of active,
        possibly self-replicating, expression rather than incidental
        transcription of one strand of host DNA.

Obelisks themselves score high on S1 and S3 but are *coding* (Oblins), so
score low on S2 -- this signature is intentionally aimed at a different,
still-unsearched corner of the same design space: stable, strand-symmetric,
but non-coding elements.
"""
from __future__ import annotations

import random
import statistics
from dataclasses import dataclass

from scipy.stats import norm

from .cluster import Cluster
from .circularity import find_circularity
from .orphan import orphan_score as compute_orphan_score
from .structure import structure_zscore

# Minimum orphan/coding-potential score (S2) required before a candidate's
# structure stability (S1) is even eligible for significance testing --
# see the rationale in score_cluster() below.
ORPHAN_GATE = 0.5


@dataclass
class Candidate:
    cluster_id: str
    representative_id: str
    length: int
    structure_z: float
    orphan_score: float
    has_sas: bool
    is_circular: bool
    combined_score: float
    p_value: float


def score_cluster(
    cluster: Cluster,
    cluster_id: str,
    n_shuffles: int = 100,
    rng: random.Random | None = None,
) -> Candidate:
    rng = rng or random.Random()
    rep = cluster.centroid
    circ = find_circularity(rep.seq)
    struct = structure_zscore(rep.seq, n_shuffles=n_shuffles, rng=rng)
    orphan = compute_orphan_score(rep.seq)

    combined = struct.z_score * orphan.orphan_score
    if cluster.has_sense_antisense():
        combined += 2.0
    if circ.is_circular:
        combined += 1.0

    # Significance testing: gate on the orphan/coding-potential axis (S2)
    # before testing structure stability (S1).
    #
    # v1 of this pipeline derived the BH-FDR p-value from the structure
    # z-score alone. Calibration (results/phase2_calibration_v1_*.json)
    # caught a real specificity failure this produces: 5/8
    # structured-coding decoys (real ORF + real stable hairpin -- i.e.
    # ordinary structured mRNA, not a novel element) were called
    # significant at alpha=0.05, because a stable fold is a stable fold
    # regardless of whether the sequence also encodes a protein. Since this
    # signature is explicitly defined as *non-coding* structure (S1 AND S2,
    # not S1 alone -- see docs/METHODS.md), the significance test needs to
    # reflect that conjunction, not just rank by it. We therefore require
    # orphan_score >= ORPHAN_GATE (a non-coding-like majority: less than
    # half the sequence explained by coding-potential evidence) before a
    # candidate is even eligible; candidates failing the gate get p=1.0
    # (never significant under BH-FDR) regardless of how stable their fold
    # is. This is a hard biological requirement of the signature's own
    # definition, not a post-hoc statistical fudge.
    if orphan.orphan_score >= ORPHAN_GATE:
        p = float(norm.sf(struct.z_score))
    else:
        p = 1.0

    return Candidate(
        cluster_id=cluster_id,
        representative_id=rep.id,
        length=len(rep.seq),
        structure_z=struct.z_score,
        orphan_score=orphan.orphan_score,
        has_sas=cluster.has_sense_antisense(),
        is_circular=circ.is_circular,
        combined_score=combined,
        p_value=p,
    )


def benjamini_hochberg(p_values: list[float], alpha: float = 0.05) -> tuple[list[bool], list[float]]:
    """Standard BH step-up FDR control. Returns (reject, q_values) in the
    original input order."""
    n = len(p_values)
    if n == 0:
        return [], []
    order = sorted(range(n), key=lambda i: p_values[i])
    ranked_p = [p_values[i] for i in order]

    q_sorted = [0.0] * n
    min_q = 1.0
    for rank in range(n, 0, -1):
        i = rank - 1
        q = ranked_p[i] * n / rank
        min_q = min(min_q, q)
        q_sorted[i] = min_q

    q_values = [0.0] * n
    for rank_pos, orig_i in enumerate(order):
        q_values[orig_i] = q_sorted[rank_pos]

    reject = [q <= alpha for q in q_values]
    return reject, q_values
