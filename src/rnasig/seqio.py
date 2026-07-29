"""FASTA I/O and basic sequence utilities."""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Iterator

_COMPLEMENT = str.maketrans("ACGTUacgtuNn", "TGCAAtgcaaNn")


@dataclass
class Record:
    id: str
    seq: str

    def __len__(self) -> int:
        return len(self.seq)


def read_fasta(path: str) -> list[Record]:
    records: list[Record] = []
    header = None
    chunks: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(Record(header, "".join(chunks).upper()))
                header = line[1:]
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            records.append(Record(header, "".join(chunks).upper()))
    return records


def write_fasta(path: str, records: list[Record], wrap: int = 70) -> None:
    with open(path, "w") as fh:
        for rec in records:
            fh.write(f">{rec.id}\n")
            seq = rec.seq
            for i in range(0, len(seq), wrap):
                fh.write(seq[i : i + wrap] + "\n")


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def dna(seq: str) -> str:
    return seq.upper().replace("U", "T")


def gc_content(seq: str) -> float:
    seq = seq.upper()
    if not seq:
        return 0.0
    gc = sum(seq.count(b) for b in "GC")
    return gc / len(seq)


_LEN_RE = re.compile(r"_length_(\d+)_")


def spades_length(seq_id: str) -> int | None:
    """Extract the assembler-reported length from an rnaSPAdes-style contig id
    (e.g. NODE_1_length_943_cov_17.4_g0_i0). Falls back to None if absent."""
    m = _LEN_RE.search(seq_id)
    return int(m.group(1)) if m else None
