"""Remote homology resolution for surviving candidates.

`docs/LIMITATIONS.md` used to end every candidate discussion the same
way: whether a contig is novel or merely unfamiliar cannot be settled
without nt / nr / Rfam, and those were unreachable. That is the step this
module performs, so it stops being a promise and becomes a pipeline
stage that anyone can re-run.

Three independent questions, three services:

  - Is it a known *sequence*?   blastn against core_nt (NCBI URL API).
  - Does it encode a known *protein*?  blastx against nr. This is the one
    that mattered: k141_25, which the pipeline scored as "borderline
    non-coding", is 97% identical at the amino-acid level to bacterial
    FAD-dependent oxidoreductase. A nucleotide search alone would have
    left that ambiguous, because the nucleotide identity to any single
    genome is only 78%.
  - Is it a known structured *RNA family*?  Infernal cmscan against Rfam
    (batch.rfam.org).

Everything here is polite to the public services: one submission at a
time, polling on the interval NCBI asks for, and an identifying email on
each request. Every call fails soft, so an offline run reports
"unresolved" rather than raising -- the same convention rrna_filter uses
for missing binaries.
"""
from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field

from .seqio import Record

NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
RFAM_SUBMIT_URL = "https://batch.rfam.org/submit-job"
RFAM_RESULT_URL = "https://batch.rfam.org/result"

# NCBI asks that automated clients identify themselves and poll no faster
# than every 60s. We honour both.
CONTACT_EMAIL = "newaashish190@gmail.com"
POLL_INTERVAL = 20.0
MAX_WAIT = 900.0


@dataclass
class Hit:
    accession: str
    description: str
    score: float
    evalue: float
    identity: float | None = None

    def describe(self) -> str:
        ident = f", {self.identity:.0f}% id" if self.identity is not None else ""
        return f"{self.accession} {self.description} (E={self.evalue:g}{ident})"


@dataclass
class Identification:
    query_id: str
    length: int
    nt_hits: list[Hit] = field(default_factory=list)
    protein_hits: list[Hit] = field(default_factory=list)
    rfam_hits: list[dict] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)

    @property
    def resolved(self) -> bool:
        return bool(self.nt_hits or self.protein_hits or self.rfam_hits)

    def verdict(self) -> str:
        """One line describing what this sequence turned out to be."""
        if self.rfam_hits:
            names = ", ".join(sorted({h.get("id", "?") for h in self.rfam_hits}))
            return f"known structured RNA family (Rfam: {names})"
        if self.protein_hits:
            best = self.protein_hits[0]
            ident = f" at {best.identity:.0f}% aa identity" if best.identity is not None else ""
            return f"protein-coding: {best.description}{ident}"
        if self.nt_hits:
            best = self.nt_hits[0]
            return f"matches known sequence: {best.description} (E={best.evalue:g})"
        if self.errors:
            return "unresolved (lookup failed)"
        return "no significant database match"

    def to_dict(self) -> dict:
        return {
            "query_id": self.query_id,
            "length": self.length,
            "verdict": self.verdict(),
            "nt_hits": [vars(h) for h in self.nt_hits],
            "protein_hits": [vars(h) for h in self.protein_hits],
            "rfam_hits": self.rfam_hits,
            "errors": self.errors,
        }


def _post(url: str, data: dict, timeout: float = 120.0) -> str:
    body = urllib.parse.urlencode(data).encode()
    req = urllib.request.Request(url, data=body, headers={"User-Agent": f"rnasig ({CONTACT_EMAIL})"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8", "replace")


def _get(url: str, timeout: float = 120.0) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": f"rnasig ({CONTACT_EMAIL})"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8", "replace")


def _parse_rid(text: str) -> str | None:
    for line in text.splitlines():
        if "RID = " in line:
            return line.split("RID = ", 1)[1].strip()
    return None


def parse_blast_descriptions(report: str, max_hits: int = 10) -> list[Hit]:
    """Pull the description table out of a plain-text BLAST report.

    The table sits between the 'Sequences producing significant alignments'
    header and the 'ALIGNMENTS' marker; each row ends with score, E-value
    and percent identity.
    """
    hits: list[Hit] = []
    in_table = False
    for line in report.splitlines():
        stripped = line.strip()
        if stripped.startswith("Sequences producing significant alignments"):
            in_table = True
            continue
        if not in_table:
            continue
        if stripped.startswith("ALIGNMENTS") or stripped.startswith(">"):
            break
        if not stripped or stripped.startswith("("):
            continue

        parts = stripped.split()
        if len(parts) < 3:
            continue
        identity = None
        tail = parts[-1]
        if tail.endswith("%"):
            try:
                identity = float(tail.rstrip("%"))
            except ValueError:
                identity = None
            parts = parts[:-1]
        if len(parts) < 3:
            continue
        try:
            evalue = float(parts[-1].replace("e", "E"))
            score = float(parts[-2])
        except ValueError:
            continue
        accession = parts[0]
        description = " ".join(parts[1:-2])
        hits.append(Hit(accession=accession, description=description,
                        score=score, evalue=evalue, identity=identity))
        if len(hits) >= max_hits:
            break
    return hits


def blast_remote(
    record: Record,
    program: str = "blastn",
    database: str = "core_nt",
    megablast: bool = False,
    max_hits: int = 10,
    poll_interval: float = POLL_INTERVAL,
    max_wait: float = MAX_WAIT,
) -> tuple[list[Hit], str | None]:
    """Run one query through the NCBI BLAST URL API. Returns (hits, error)."""
    fasta = f">{record.id.split()[0]}\n{record.seq}\n"
    submit = {
        "CMD": "Put",
        "PROGRAM": program,
        "DATABASE": database,
        "QUERY": fasta,
        "HITLIST_SIZE": str(max(max_hits, 20)),
        "EMAIL": CONTACT_EMAIL,
    }
    if megablast:
        submit["MEGABLAST"] = "on"

    try:
        rid = _parse_rid(_post(NCBI_BLAST_URL, submit))
    except (urllib.error.URLError, OSError) as exc:
        return [], f"{program}: submission failed ({exc})"
    if not rid:
        return [], f"{program}: no RID returned"

    waited = 0.0
    status_url = f"{NCBI_BLAST_URL}?CMD=Get&FORMAT_OBJECT=SearchInfo&RID={rid}"
    while waited < max_wait:
        time.sleep(poll_interval)
        waited += poll_interval
        try:
            info = _get(status_url)
        except (urllib.error.URLError, OSError) as exc:
            return [], f"{program}: polling failed ({exc})"
        if "Status=FAILED" in info or "Status=UNKNOWN" in info:
            return [], f"{program}: search failed server-side (RID {rid})"
        if "Status=READY" in info:
            break
    else:
        return [], f"{program}: timed out after {max_wait:.0f}s (RID {rid})"

    try:
        report = _get(
            f"{NCBI_BLAST_URL}?CMD=Get&RID={rid}&FORMAT_TYPE=Text"
            f"&DESCRIPTIONS={max(max_hits, 20)}&ALIGNMENTS=5"
        )
    except (urllib.error.URLError, OSError) as exc:
        return [], f"{program}: result fetch failed ({exc})"

    if "No significant similarity found" in report:
        return [], None
    return parse_blast_descriptions(report, max_hits=max_hits), None


def rfam_scan(record: Record, timeout: float = 60.0, max_wait: float = 180.0) -> tuple[list[dict], str | None]:
    """Run Infernal cmscan against Rfam via the public batch service."""
    fasta = f">{record.id.split()[0]}\n{record.seq}\n"
    boundary = "----rnasig-rfam-boundary"
    body = (
        f"--{boundary}\r\n"
        f'Content-Disposition: form-data; name="sequence_file"; filename="query.fasta"\r\n'
        f"Content-Type: text/plain\r\n\r\n{fasta}\r\n--{boundary}--\r\n"
    ).encode()
    req = urllib.request.Request(
        RFAM_SUBMIT_URL,
        data=body,
        headers={
            "Content-Type": f"multipart/form-data; boundary={boundary}",
            "Accept": "application/json",
            "User-Agent": f"rnasig ({CONTACT_EMAIL})",
        },
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            submitted = json.loads(resp.read().decode("utf-8", "replace"))
    except (urllib.error.URLError, OSError, ValueError) as exc:
        return [], f"rfam: submission failed ({exc})"

    job_id = submitted.get("jobId")
    if not job_id:
        return [], "rfam: no jobId returned"

    waited = 0.0
    while waited < max_wait:
        time.sleep(5.0)
        waited += 5.0
        try:
            payload = json.loads(_get(f"{RFAM_RESULT_URL}/{job_id}", timeout=timeout))
        except (urllib.error.URLError, OSError, ValueError):
            continue
        if "numHits" not in payload:
            continue
        hits = payload.get("hits") or {}
        flat: list[dict] = []
        if isinstance(hits, dict):
            for family, entries in hits.items():
                for entry in entries if isinstance(entries, list) else [entries]:
                    record_hit = dict(entry) if isinstance(entry, dict) else {"raw": entry}
                    record_hit.setdefault("id", family)
                    flat.append(record_hit)
        elif isinstance(hits, list):
            flat = [h for h in hits if isinstance(h, dict)]
        return flat, None

    return [], "rfam: timed out"


def identify(
    record: Record,
    do_nt: bool = True,
    do_protein: bool = True,
    do_rfam: bool = True,
    max_hits: int = 10,
) -> Identification:
    """Resolve one candidate against nt, nr and Rfam."""
    result = Identification(query_id=record.id.split()[0], length=len(record.seq))

    if do_rfam:
        hits, err = rfam_scan(record)
        result.rfam_hits = hits
        if err:
            result.errors.append(err)

    if do_nt:
        hits, err = blast_remote(record, program="blastn", database="core_nt", max_hits=max_hits)
        result.nt_hits = hits
        if err:
            result.errors.append(err)

    if do_protein:
        hits, err = blast_remote(record, program="blastx", database="nr", max_hits=max_hits)
        result.protein_hits = hits
        if err:
            result.errors.append(err)

    return result
