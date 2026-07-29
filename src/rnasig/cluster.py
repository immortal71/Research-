"""Circular-permutation-aware clustering + sense/antisense (SAS) strand
assignment -- reimplements the conceptual role of VNom's circUCLUST +
SASFinder steps using Biopython's PairwiseAligner instead of usearch
(which is proprietary and unavailable here). See
data/reference/vnom/PROVENANCE.md.

Two circular sequences can look completely different as raw strings even
when they represent "the same" molecule, because the assembler can start
the contig at any point around the circle. The standard trick -- and the
one used here -- is to align the *query* against the *target doubled with
itself* (target+target), so any rotation of the query can still find a
contiguous local match.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from Bio.Align import PairwiseAligner, substitution_matrices

from .seqio import Record, revcomp


def _aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 1.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -2.0
    aligner.extend_gap_score = -0.5
    return aligner


_ALIGNER = _aligner()


def best_circular_identity(query: str, target: str) -> tuple[float, str]:
    """Return (identity, orientation) of the best local alignment of `query`
    (in either orientation) against `target` doubled to tolerate circular
    permutation. identity is matches / min(len(query), len(target))."""
    doubled = target + target
    denom = min(len(query), len(target))
    if denom == 0:
        return 0.0, "sense"

    best_score = -1.0
    best_orientation = "sense"
    for orientation, q in (("sense", query), ("antisense", revcomp(query))):
        alignment = _ALIGNER.align(q, doubled)[0]
        matches = _count_matches(alignment)
        identity = matches / denom
        if identity > best_score:
            best_score = identity
            best_orientation = orientation
    return best_score, best_orientation


def _count_matches(alignment) -> int:
    a_str, b_str = str(alignment[0]), str(alignment[1])
    return sum(1 for x, y in zip(a_str, b_str) if x == y and x != "-")


@dataclass
class Cluster:
    centroid: Record
    members: list[Record] = field(default_factory=list)
    orientations: list[str] = field(default_factory=list)  # per member, relative to centroid

    def has_sense_antisense(self) -> bool:
        """A member other than the centroid itself assigned 'antisense', in
        the presence of at least one 'sense' assignment (the centroid, by
        convention, defines 'sense'), is evidence of both-strand
        co-occurrence -- the hallmark obelisks/VNom signature of an
        actively-replicating (rather than merely transcribed) element."""
        return "antisense" in self.orientations and "sense" in (["sense"] + self.orientations)


def cluster_sequences(records: list[Record], id_threshold: float = 0.7) -> list[Cluster]:
    """Greedy circular-permutation-aware clustering. Longest sequences become
    centroids first (mirrors circUCLUST's length-sorted greedy behaviour)."""
    ordered = sorted(records, key=lambda r: len(r), reverse=True)
    clusters: list[Cluster] = []

    for rec in ordered:
        best_cluster = None
        best_identity = 0.0
        best_orientation = "sense"
        for cl in clusters:
            identity, orientation = best_circular_identity(rec.seq, cl.centroid.seq)
            if identity > best_identity:
                best_identity = identity
                best_cluster = cl
                best_orientation = orientation
        if best_cluster is not None and best_identity >= id_threshold:
            best_cluster.members.append(rec)
            best_cluster.orientations.append(best_orientation)
        else:
            clusters.append(Cluster(centroid=rec))
    return clusters


def non_singleton_clusters(clusters: list[Cluster]) -> list[Cluster]:
    return [c for c in clusters if c.members]


def sas_clusters(clusters: list[Cluster]) -> list[Cluster]:
    return [c for c in non_singleton_clusters(clusters) if c.has_sense_antisense()]
