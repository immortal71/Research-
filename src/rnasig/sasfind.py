"""Sense/antisense co-occurrence, found by k-mer sharing rather than alignment.

`cluster.py` already does this properly, with circular-permutation-aware
pairwise alignment, and it is the right tool on the tens of contigs Phase 1
deals with. It is quadratic in pairwise alignments, so it does not survive
contact with the 10,000-25,000 contigs a real sweep produces, and
`phase6_hunt` consequently shipped without any SAS test at all.

That was a real gap rather than an optimisation. The signature this
pipeline hunts is stable structure plus low coding potential plus *either*
circularity *or* both strands of the same molecule turning up in one
sample. Phase 6 tested only circularity, which is the harder of the two to
observe: it needs the assembler to read through the junction and close the
loop, and on Obelisk-S.s that is exactly what failed. Sense/antisense
co-occurrence needs no such thing. It only needs the molecule to be
transcribed both ways, which is what a replicating RNA does.

`assemble` is deliberately strand-specific and does not canonicalise
k-mers, so both strands of a molecule arrive as separate contigs and the
evidence is already there to read. This module reads it in one pass over
the k-mers instead of a pass over every pair.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field

from .seqio import Record, revcomp

DEFAULT_K = 21
# Fraction of one contig's k-mers that must appear, reverse-complemented, in
# its partner. Set high enough that incidental sharing between unrelated
# transcripts does not register, low enough to tolerate the variation a
# quasispecies carries between its two strands.
DEFAULT_MIN_SHARED = 0.30


@dataclass
class SASPair:
    sense_id: str
    antisense_id: str
    sense_len: int
    antisense_len: int
    shared_fraction: float
    sense_coverage: float = 0.0
    antisense_coverage: float = 0.0

    def describe(self) -> str:
        return (
            f"{self.sense_id} ({self.sense_len} nt, {self.sense_coverage:.0f}x) "
            f"<-> {self.antisense_id} ({self.antisense_len} nt, "
            f"{self.antisense_coverage:.0f}x) sharing {self.shared_fraction:.0%}"
        )


@dataclass
class SASReport:
    n_contigs: int = 0
    pairs: list[SASPair] = field(default_factory=list)

    def describe(self) -> str:
        return f"SAS: {len(self.pairs)} sense/antisense pairs among {self.n_contigs} contigs"


def _coverage_of(record: Record) -> float:
    if "multi=" in record.id:
        try:
            return float(record.id.split("multi=")[1].split()[0])
        except (IndexError, ValueError):
            return 0.0
    if "cov=" in record.id:
        try:
            return float(record.id.split("cov=")[1].split()[0])
        except (IndexError, ValueError):
            return 0.0
    return 0.0


def find_sas_pairs(
    records: list[Record],
    k: int = DEFAULT_K,
    min_shared: float = DEFAULT_MIN_SHARED,
    min_length: int = 150,
) -> SASReport:
    """Find contigs that are the reverse complement of another contig.

    One pass builds a k-mer index; a second looks each contig's reverse
    complement up in it. Cost is linear in total sequence length rather
    than quadratic in contig count, which is what makes this usable on a
    whole assembly.
    """
    usable = [r for r in records if len(r.seq) >= min_length]
    report = SASReport(n_contigs=len(usable))
    if not usable:
        return report

    index: dict[str, set[int]] = defaultdict(set)
    kmer_sets: list[set[str]] = []
    for i, record in enumerate(usable):
        seq = record.seq.upper()
        kmers = {seq[j : j + k] for j in range(len(seq) - k + 1)}
        kmer_sets.append(kmers)
        for kmer in kmers:
            index[kmer].add(i)

    seen: set[tuple[int, int]] = set()
    for i, record in enumerate(usable):
        rc_kmers = {revcomp(km) for km in kmer_sets[i]}
        hits: dict[int, int] = defaultdict(int)
        for kmer in rc_kmers:
            for j in index.get(kmer, ()):  # partners carrying this k-mer
                if j != i:
                    hits[j] += 1
        for j, shared in hits.items():
            pair = (min(i, j), max(i, j))
            if pair in seen:
                continue
            fraction = shared / min(len(kmer_sets[i]), len(kmer_sets[j]))
            if fraction < min_shared:
                continue
            seen.add(pair)
            report.pairs.append(
                SASPair(
                    sense_id=usable[i].id.split()[0],
                    antisense_id=usable[j].id.split()[0],
                    sense_len=len(usable[i].seq),
                    antisense_len=len(usable[j].seq),
                    shared_fraction=round(fraction, 3),
                    sense_coverage=_coverage_of(usable[i]),
                    antisense_coverage=_coverage_of(usable[j]),
                )
            )

    report.pairs.sort(key=lambda p: -max(p.sense_coverage, p.antisense_coverage))
    return report
