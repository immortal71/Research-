"""Library-type preflight for public sequencing runs.

Why this module exists: the Phase 3 "real" sweep in this repo was run on
SRR13291825 on the strength of it being a public environmental run that
the sandbox could actually reach. Nobody checked what kind of library it
was. It is an AMPLICON / PCR / METAGENOMIC run -- an 18S rRNA
metabarcoding survey of soil DNA. There is no RNA in it at all, so every
RNA-level statement downstream of that sweep was void at the source. See
results/phase5_identification_report.md.

The signature this pipeline hunts (stable structure + non-coding + dual
strand / circularity) is only meaningful on shotgun RNA. Targeted PCR
libraries break it in three separate ways at once:

  - they are DNA, so "RNA structure" is not a property of the molecule;
  - one locus is amplified to saturation, so coverage carries no
    abundance information and "high coverage" stops being evidence;
  - PCR generates chimeras and concatemers, which the circularity
    detector reads as terminal repeats.

So this is a hard gate, checked before a sweep starts rather than argued
about afterwards. `assess_library` is a pure function over metadata and
is what the tests exercise; `fetch_run_metadata` is the network half and
fails soft so an offline run degrades to "unknown" instead of crashing.
"""
from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"

_FIELDS = (
    "run_accession,library_strategy,library_selection,library_source,"
    "instrument_platform,scientific_name,read_count,fastq_ftp,study_accession"
)

# Library strategies whose reads can carry the SOS signature. Anything that
# is not shotgun RNA is rejected: the signature assumes untargeted sampling
# of an RNA population.
_ALLOWED_STRATEGY = {"RNA-SEQ", "SSRNA-SEQ", "MRNA-SEQ", "TOTAL-RNA-SEQ", "NCRNA-SEQ", "MIRNA-SEQ"}
_ALLOWED_SOURCE = {"TRANSCRIPTOMIC", "METATRANSCRIPTOMIC", "TRANSCRIPTOMIC SINGLE CELL"}

# Selection methods that mean "one locus was amplified on purpose". Coverage
# and circularity are both uninterpretable under these.
_TARGETED_SELECTION = {"PCR", "RT-PCR", "REPEAT FRACTIONATION", "5-METHYLCYTIDINE ANTIBODY",
                       "MF", "MSLL", "HMPR", "CF-S", "CF-M", "CF-H", "CF-T"}


@dataclass
class RunMetadata:
    accession: str
    library_strategy: str = ""
    library_selection: str = ""
    library_source: str = ""
    instrument_platform: str = ""
    scientific_name: str = ""
    read_count: int | None = None
    study_accession: str = ""
    fetched: bool = False

    @property
    def summary(self) -> str:
        return (f"{self.accession}: strategy={self.library_strategy or '?'} "
                f"selection={self.library_selection or '?'} "
                f"source={self.library_source or '?'} "
                f"({self.scientific_name or 'unknown sample'})")


@dataclass
class LibraryVerdict:
    accession: str
    compatible: bool
    unknown: bool = False
    reasons: list[str] = field(default_factory=list)

    def __bool__(self) -> bool:
        return self.compatible


def assess_library(meta: RunMetadata) -> LibraryVerdict:
    """Decide whether a run's library type can support the SOS signature.

    Returns compatible=False with an explicit reason list when it cannot.
    A run with no usable metadata comes back unknown=True and
    compatible=False -- absent metadata is not treated as permission,
    because that is exactly how SRR13291825 got swept in the first place.
    """
    strategy = (meta.library_strategy or "").strip().upper()
    selection = (meta.library_selection or "").strip().upper()
    source = (meta.library_source or "").strip().upper()

    if not any((strategy, selection, source)):
        return LibraryVerdict(
            accession=meta.accession, compatible=False, unknown=True,
            reasons=["no library metadata available; cannot confirm this is shotgun RNA"],
        )

    reasons: list[str] = []
    if strategy == "AMPLICON":
        reasons.append(
            "library_strategy=AMPLICON: a targeted marker-gene survey. One locus is "
            "amplified to saturation, so coverage is not abundance and PCR concatemers "
            "mimic circularity."
        )
    elif strategy and strategy not in _ALLOWED_STRATEGY:
        reasons.append(f"library_strategy={meta.library_strategy}: not a shotgun RNA strategy")

    if selection in _TARGETED_SELECTION:
        reasons.append(
            f"library_selection={meta.library_selection}: template was amplified from "
            "specific primers, so contigs are PCR products rather than sampled molecules"
        )

    if source in {"GENOMIC", "METAGENOMIC", "GENOMIC SINGLE CELL"}:
        reasons.append(
            f"library_source={meta.library_source}: the molecules sequenced are DNA, so "
            "RNA secondary-structure scores describe a molecule that was never present"
        )
    elif source and source not in _ALLOWED_SOURCE:
        reasons.append(f"library_source={meta.library_source}: not an RNA source")

    return LibraryVerdict(accession=meta.accession, compatible=not reasons, reasons=reasons)


def fetch_run_metadata(accession: str, timeout: float = 30.0, attempts: int = 3) -> RunMetadata:
    """Pull library metadata for an SRA/ENA run accession from the ENA portal.

    Retries before giving up. A single dropped request is indistinguishable
    from a run with no metadata once it reaches assess_library, and both are
    refused, so a transient failure silently costs a run that would have
    passed. That happened three times over one sweep of two dozen
    accessions, which is often enough to be worth a retry.

    Fails soft after the last attempt: any network or parse problem returns
    an unfetched RunMetadata, which assess_library reports as unknown and
    refuses. Absent metadata is never treated as permission.
    """
    query = urllib.parse.urlencode(
        {"accession": accession, "result": "read_run", "fields": _FIELDS, "format": "tsv"}
    )
    text = None
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(f"{ENA_PORTAL}?{query}", timeout=timeout) as resp:
                text = resp.read().decode("utf-8", "replace")
            break
        except (urllib.error.URLError, OSError, ValueError):
            if attempt == attempts - 1:
                return RunMetadata(accession=accession)
            time.sleep(2.0 * (attempt + 1))
    if text is None:
        return RunMetadata(accession=accession)

    lines = [ln for ln in text.splitlines() if ln.strip()]
    if len(lines) < 2:
        return RunMetadata(accession=accession)

    header = lines[0].split("\t")
    values = lines[1].split("\t")
    row = dict(zip(header, values))

    try:
        read_count = int(row.get("read_count", "") or 0) or None
    except ValueError:
        read_count = None

    return RunMetadata(
        accession=row.get("run_accession", accession),
        library_strategy=row.get("library_strategy", ""),
        library_selection=row.get("library_selection", ""),
        library_source=row.get("library_source", ""),
        instrument_platform=row.get("instrument_platform", ""),
        scientific_name=row.get("scientific_name", ""),
        read_count=read_count,
        study_accession=row.get("study_accession", ""),
        fetched=True,
    )


def preflight(accession: str, timeout: float = 30.0) -> tuple[RunMetadata, LibraryVerdict]:
    """Convenience wrapper: fetch metadata for a run and judge it."""
    meta = fetch_run_metadata(accession, timeout=timeout)
    return meta, assess_library(meta)


def to_json(meta: RunMetadata, verdict: LibraryVerdict) -> str:
    return json.dumps(
        {
            "accession": meta.accession,
            "library_strategy": meta.library_strategy,
            "library_selection": meta.library_selection,
            "library_source": meta.library_source,
            "scientific_name": meta.scientific_name,
            "study_accession": meta.study_accession,
            "read_count": meta.read_count,
            "metadata_fetched": meta.fetched,
            "compatible": verdict.compatible,
            "unknown": verdict.unknown,
            "reasons": verdict.reasons,
        },
        indent=2,
    )
