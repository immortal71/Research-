"""Known-ncRNA exclusion filter, using barrnap (rRNAs) and tRNAscan-SE
(tRNAs) if either is available on the system PATH.

Motivation: on real environmental data the SOS signature fires most
strongly on rRNAs -- they are structured, dual-strand-abundant, and
high-coverage. That was confirmed empirically on SRR13291825 (see
results/phase3_real_sweep_report.md): the top surviving hit after the
adapter filter was 18S rRNA. To reach genuinely novel candidates, known-
ncRNA contigs need to be excluded before the sweep flags anything.

Both tools use HMM libraries shipped with their apt packages -- no online
lookup is done. If a tool is not installed, this filter silently skips
that class (with a note in the report), rather than failing.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field

from .seqio import Record, write_fasta


@dataclass
class NcRNAFilterReport:
    n_input: int
    n_kept: int
    n_dropped: int
    barrnap_available: bool
    trnascan_available: bool
    known_by_id: dict[str, list[str]] = field(default_factory=dict)


def _run_barrnap(fasta_path: str, kingdom: str, reject: float = 0.05, evalue: str = "1e-3") -> list[str]:
    if not shutil.which("barrnap"):
        return []
    r = subprocess.run(
        ["barrnap", "--kingdom", kingdom, "--reject", str(reject),
         "--evalue", evalue, "--quiet", fasta_path],
        capture_output=True, text=True, check=False,
    )
    hits = []
    for line in r.stdout.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) >= 1:
            hits.append(parts[0])
    return hits


def _run_trnascan(fasta_path: str) -> list[str]:
    if not shutil.which("tRNAscan-SE"):
        return []
    with tempfile.NamedTemporaryFile(mode="r", suffix=".trna", delete=False) as fh:
        out_path = fh.name
    try:
        subprocess.run(
            ["tRNAscan-SE", "-G", "-q", "-o", out_path, fasta_path],
            capture_output=True, text=True, check=False,
        )
        hits = []
        for line in open(out_path):
            line = line.strip()
            # tRNAscan-SE header lines start with "Sequence", "Name", "-" etc.
            if not line or line.startswith(("Sequence", "Name", "---")):
                continue
            parts = line.split()
            if parts:
                hits.append(parts[0])
        return hits
    finally:
        try:
            os.unlink(out_path)
        except OSError:
            pass


def filter_known_ncrna(records: list[Record], include_trna: bool = False) -> tuple[list[Record], NcRNAFilterReport]:
    """Drop contigs identified as known ncRNAs.

    include_trna: opt-in only. tRNAscan-SE loads a full covariance-model
    database on every invocation and is slow (tens of seconds even on tiny
    inputs). On our real SRR13291825 sweep it found zero tRNAs, so it is
    off by default; enable when you specifically expect tRNA-length noise.
    """
    barrnap_ok = shutil.which("barrnap") is not None
    trnascan_ok = include_trna and (shutil.which("tRNAscan-SE") is not None)

    if not records or (not barrnap_ok and not trnascan_ok):
        return records, NcRNAFilterReport(
            n_input=len(records), n_kept=len(records), n_dropped=0,
            barrnap_available=barrnap_ok, trnascan_available=trnascan_ok,
        )

    known: dict[str, list[str]] = {}
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as fh:
        fasta_path = fh.name
    try:
        write_fasta(fasta_path, records)
        if barrnap_ok:
            for kingdom in ("bac", "arc", "euk", "mito"):
                for cid in _run_barrnap(fasta_path, kingdom):
                    known.setdefault(cid, []).append(f"rRNA({kingdom})")
        if trnascan_ok:
            for cid in _run_trnascan(fasta_path):
                known.setdefault(cid, []).append("tRNA")
    finally:
        try:
            os.unlink(fasta_path)
        except OSError:
            pass

    # FASTA convention (and barrnap/tRNAscan-SE follow it): the sequence ID
    # is the first whitespace-delimited token of the header. Assembler outputs
    # (MEGAHIT, SPAdes) put extra metadata after a space, so records loaded by
    # rnasig.seqio carry the full header string. Compare on the leading token.
    def leading_token(s: str) -> str:
        return s.split()[0] if s else s

    kept = [r for r in records if leading_token(r.id) not in known]
    return kept, NcRNAFilterReport(
        n_input=len(records),
        n_kept=len(kept),
        n_dropped=len(known),
        barrnap_available=barrnap_ok,
        trnascan_available=trnascan_ok,
        known_by_id={k: sorted(set(v)) for k, v in known.items()},
    )
