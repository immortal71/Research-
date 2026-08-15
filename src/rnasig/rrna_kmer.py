"""Reference-based rRNA exclusion that does not need barrnap.

`rrna_filter` wraps barrnap and tRNAscan-SE, which is the right tool when
they are installed. They are apt packages with no Windows build, so on the
machine this pipeline was developed on the rRNA stage silently did nothing,
and `phase6_hunt` ran without it entirely.

That is not a cosmetic gap. Ribosomal RNA is the dominant false positive for
this signature and hits all three axes at once: it is extremely structured,
it is the most abundant thing in any metatranscriptome, and an rRNA operon
assembles with terminal repeats that read as circular. The first oral sweep
produced a 342 nt contig at 6,752x coverage that cleared both the structure
and rod-likeness bars, and blastn put it at 100% identity to 28S rRNA.

This module screens against reference rRNA by exact k-mer sharing instead.
A contig that shares a meaningful fraction of its k-mers with a reference is
rRNA regardless of what any structure score says. Exact matching is crude
next to a covariance model and is deliberately so: it needs no binary, no
alignment, and no network once the reference is cached, and rRNA is
conserved enough that k-mer sharing is a strong signal.

References are fetched once from NCBI and cached on disk. With no cache and
no network the filter fails open, reporting that it did nothing, in the same
spirit as rrna_filter's soft failure when barrnap is absent.
"""
from __future__ import annotations

import os
import urllib.error
import urllib.request
from dataclasses import dataclass, field

from .seqio import Record, revcomp

# SSU and LSU across bacteria, archaea and eukaryotes, plus 5S/5.8S. These
# are nuccore accessions, fetched as FASTA.
_REFERENCE_ACCESSIONS = [
    "NR_046235.3",   # Homo sapiens 45S pre-rRNA (18S, 5.8S, 28S)
    "J01695.2",      # Escherichia coli rrnB operon (16S, 23S, 5S)
    "NR_044838.1",   # Methanocaldococcus jannaschii 16S
    "NR_029158.1",   # Homo sapiens 5S
    "CP000411.2",    # Oenococcus oeni: another bacterial rRNA source
]

_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

DEFAULT_K = 21
# A contig sharing this fraction of its k-mers with reference rRNA is rRNA.
# Conserved regions alone will push a genuine rRNA fragment well past this,
# while an unrelated sequence shares essentially nothing at k=21.
DEFAULT_THRESHOLD = 0.10


@dataclass
class RRNAScreen:
    n_input: int = 0
    n_kept: int = 0
    n_dropped: int = 0
    reference_loaded: bool = False
    dropped: dict[str, float] = field(default_factory=dict)

    def describe(self) -> str:
        if not self.reference_loaded:
            return "rRNA screen skipped: no reference available"
        return f"rRNA screen: {self.n_input} -> {self.n_kept} kept, {self.n_dropped} dropped"


def fetch_reference(cache_path: str, timeout: float = 60.0) -> str | None:
    """Return reference rRNA FASTA, downloading and caching it if needed."""
    if os.path.exists(cache_path) and os.path.getsize(cache_path) > 1000:
        with open(cache_path, encoding="utf-8", errors="replace") as fh:
            return fh.read()

    ids = ",".join(_REFERENCE_ACCESSIONS)
    url = f"{_EUTILS}?db=nuccore&id={ids}&rettype=fasta&retmode=text"
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            text = resp.read().decode("utf-8", "replace")
    except (urllib.error.URLError, OSError):
        return None
    if ">" not in text:
        return None

    os.makedirs(os.path.dirname(cache_path) or ".", exist_ok=True)
    with open(cache_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)
    return text


def build_reference_kmers(fasta_text: str, k: int = DEFAULT_K) -> set[str]:
    """Collect k-mers from reference rRNA, both strands."""
    kmers: set[str] = set()
    seq_parts: list[str] = []

    def absorb(seq: str) -> None:
        seq = seq.upper()
        for strand in (seq, revcomp(seq)):
            for i in range(len(strand) - k + 1):
                kmer = strand[i : i + k]
                if "N" not in kmer:
                    kmers.add(kmer)

    for line in fasta_text.splitlines():
        if line.startswith(">"):
            if seq_parts:
                absorb("".join(seq_parts))
                seq_parts = []
        else:
            seq_parts.append(line.strip())
    if seq_parts:
        absorb("".join(seq_parts))
    return kmers


def rrna_kmer_fraction(seq: str, reference: set[str], k: int = DEFAULT_K) -> float:
    """Fraction of a sequence's k-mers that appear in reference rRNA."""
    seq = seq.upper()
    if len(seq) < k or not reference:
        return 0.0
    total = len(seq) - k + 1
    hits = sum(1 for i in range(total) if seq[i : i + k] in reference)
    return hits / total


def screen_records(
    records: list[Record],
    reference: set[str],
    k: int = DEFAULT_K,
    threshold: float = DEFAULT_THRESHOLD,
) -> tuple[list[Record], RRNAScreen]:
    """Drop records that share too many k-mers with reference rRNA."""
    report = RRNAScreen(n_input=len(records), reference_loaded=bool(reference))
    if not reference:
        report.n_kept = len(records)
        return records, report

    kept: list[Record] = []
    for record in records:
        fraction = rrna_kmer_fraction(record.seq, reference, k=k)
        if fraction >= threshold:
            report.dropped[record.id.split()[0]] = round(fraction, 4)
            report.n_dropped += 1
        else:
            kept.append(record)
    report.n_kept = len(kept)
    return kept, report
