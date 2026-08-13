"""Reference-free detection of amplicon (targeted PCR) read sets.

`srameta.assess_library` catches targeted libraries when the archive
metadata is present and honest. This module catches them when it is not,
which matters because metadata is frequently missing, and because a
pipeline that only trusts a submitter's annotation has no way to notice
when the annotation is wrong.

The observation it keys on is mechanical. In a shotgun library, reads
start wherever a fragment happened to break, so the distribution of 5'
k-mers is close to flat and the most common prefix accounts for a
vanishing fraction of reads. In an amplicon library, every molecule was
copied from the same handful of primers, so a very small number of 5'
k-mers accounts for most of the run.

On SRR13291825 -- the run this repo swept in Phase 3 -- three 20-mers
cover 68.5% of all 56,865 reads, in both mates, with degenerate variants
at the expected positions. All three match 18S rRNA. For comparison, a
shotgun library's top prefix is normally well under 1%.

No reference data and no network are needed: this is a property of the
reads themselves.
"""
from __future__ import annotations

import gzip
import io
from collections import Counter
from dataclasses import dataclass, field

# A run whose top few 5'-mers cover this much of the library was primed,
# not fragmented. Chosen with a wide margin: shotgun runs sit near zero and
# amplicon runs near one, so the threshold is not a delicate call.
DEFAULT_TOP_N = 3
DEFAULT_PREFIX_LEN = 20
DEFAULT_THRESHOLD = 0.20


@dataclass
class PrimerProfile:
    n_reads: int
    prefix_len: int
    top_n: int
    top_prefixes: list[tuple[str, int]] = field(default_factory=list)
    top_fraction: float = 0.0
    is_amplicon: bool = False
    threshold: float = DEFAULT_THRESHOLD

    def describe(self) -> str:
        if not self.n_reads:
            return "no reads examined"
        verdict = "amplicon-like (primer-dominated)" if self.is_amplicon else "shotgun-like"
        lines = [
            f"{self.n_reads} reads, top {self.top_n} 5'-{self.prefix_len}mers "
            f"cover {self.top_fraction:.1%} (threshold {self.threshold:.0%}) -> {verdict}"
        ]
        for seq, count in self.top_prefixes:
            lines.append(f"  {seq}  {count}  ({count / self.n_reads:.1%})")
        return "\n".join(lines)


def profile_prefixes(
    sequences,
    prefix_len: int = DEFAULT_PREFIX_LEN,
    top_n: int = DEFAULT_TOP_N,
    threshold: float = DEFAULT_THRESHOLD,
) -> PrimerProfile:
    """Score an iterable of read sequences for primer dominance."""
    counts: Counter[str] = Counter()
    n = 0
    for seq in sequences:
        seq = seq.strip().upper()
        if len(seq) < prefix_len:
            continue
        counts[seq[:prefix_len]] += 1
        n += 1

    if not n:
        return PrimerProfile(n_reads=0, prefix_len=prefix_len, top_n=top_n, threshold=threshold)

    top = counts.most_common(top_n)
    fraction = sum(c for _, c in top) / n
    return PrimerProfile(
        n_reads=n,
        prefix_len=prefix_len,
        top_n=top_n,
        top_prefixes=top,
        top_fraction=fraction,
        is_amplicon=fraction >= threshold,
        threshold=threshold,
    )


def _open_maybe_gzip(path: str):
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def iter_fastq_sequences(path: str, limit: int | None = None):
    """Yield sequence lines from a FASTQ file, transparently handling .gz."""
    with _open_maybe_gzip(path) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                yield line.strip()
                if limit is not None and (i // 4) + 1 >= limit:
                    return


def profile_fastq(
    path: str,
    prefix_len: int = DEFAULT_PREFIX_LEN,
    top_n: int = DEFAULT_TOP_N,
    threshold: float = DEFAULT_THRESHOLD,
    limit: int | None = None,
) -> PrimerProfile:
    """Profile a FASTQ (optionally gzipped) for primer dominance."""
    return profile_prefixes(
        iter_fastq_sequences(path, limit=limit),
        prefix_len=prefix_len,
        top_n=top_n,
        threshold=threshold,
    )
