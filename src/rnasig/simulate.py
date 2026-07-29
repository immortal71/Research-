"""Synthetic sequence construction for calibration: matched nulls and
signal injections.

Everything here is explicitly synthetic/constructed (not claimed to be a
real organism's sequence) -- see docs/METHODS.md and
data/motifs/PROVENANCE.md. This is deliberate: Phase 2 calibration needs
ground truth we can dial (copy number, mutation rate, GC content, decoy
type) to measure power and false-discovery rate honestly, which requires
synthetic control over the "true" label -- exactly what synthetic
injection/matched-null calibration means.
"""
from __future__ import annotations

import random
from dataclasses import dataclass, field

from .nullmodel import dinuc_shuffle, markov_background
from .seqio import revcomp

_BASES = "ACGT"
_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}

# All 61 sense codons (no stops), used to build ORF-like decoys with
# genuine, dial-able codon-usage bias.
_STOPS = {"TAA", "TAG", "TGA"}
_SENSE_CODONS = [
    a + b + c for a in _BASES for b in _BASES for c in _BASES if a + b + c not in _STOPS
]


def _complement_strand(arm: str, rng: random.Random, mismatch_rate: float = 0.0) -> str:
    """Reverse-complement an arm, optionally introducing mismatches (to make
    imperfect, more realistic stems)."""
    comp = []
    for base in reversed(arm):
        b = _COMPLEMENT[base]
        if mismatch_rate and rng.random() < mismatch_rate:
            b = rng.choice([x for x in _BASES if x != b])
        comp.append(b)
    return "".join(comp)


def make_stem_loop(
    rng: random.Random,
    stem_len: int = 15,
    loop_len: int = 6,
    mismatch_rate: float = 0.05,
) -> str:
    """A single hairpin: 5' stem arm + loop + complementary 3' stem arm."""
    arm5 = "".join(rng.choice(_BASES) for _ in range(stem_len))
    loop = "".join(rng.choice(_BASES) for _ in range(loop_len))
    arm3 = _complement_strand(arm5, rng, mismatch_rate)
    return arm5 + loop + arm3


def make_rod(
    rng: random.Random,
    n_hairpins: int = 6,
    stem_len: int = 14,
    loop_len: int = 5,
    linker_len: int = 3,
) -> str:
    """Concatenate several stem-loops with short linkers -- a simplified
    stand-in for the "rod-like, genome-spanning" secondary structure
    reported for obelisks: not a claim of structural equivalence, just a
    sequence guaranteed to fold with multiple stacked stable helices rather
    than being unstructured."""
    parts = []
    for i in range(n_hairpins):
        parts.append(make_stem_loop(rng, stem_len, loop_len))
        if i < n_hairpins - 1 and linker_len:
            parts.append("".join(rng.choice(_BASES) for _ in range(linker_len)))
    return "".join(parts)


def make_cloverleaf_like(rng: random.Random, arm_len: int = 10, loop_len: int = 5) -> str:
    """Four stem-loops back to back -- a coarse stand-in for tRNA-like
    multi-branch (cloverleaf) topology; used as a *different* structure
    class from the rod, for negative/diversity controls."""
    return make_rod(rng, n_hairpins=4, stem_len=arm_len, loop_len=loop_len, linker_len=2)


def make_orf_like(rng: random.Random, n_codons: int = 150, bias_strength: float = 3.0) -> str:
    """A sequence dominated by one long, codon-biased ORF -- a genuine
    "looks coding" decoy so the orphan-score axis (S2) is tested for
    specificity, not just sensitivity."""
    # Zipf-like weights over the 61 sense codons for a reproducible, strong bias.
    weights = [1.0 / (rank + 1) ** (bias_strength / 10.0) for rank in range(len(_SENSE_CODONS))]
    shuffled_codons = _SENSE_CODONS[:]
    rng.shuffle(shuffled_codons)
    codons = rng.choices(shuffled_codons, weights=weights, k=n_codons)
    return "ATG" + "".join(codons) + rng.choice(list(_STOPS))


def make_circular_contig(monomer: str, overhang_frac: float, rng: random.Random) -> str:
    """Simulate an assembler's read-through-the-origin artifact: append a
    short duplicate of the 5' end onto the 3' end, exactly the terminal
    repeat that circularity.find_circularity looks for."""
    overhang = max(1, int(len(monomer) * overhang_frac))
    return monomer + monomer[:overhang]


def make_sas_pair(
    monomer: str, rng: random.Random, mutation_rate: float = 0.03
) -> tuple[str, str]:
    """Simulate independently-assembled sense and antisense contigs of the
    same replicating element (as stranded RNA-seq would yield for an
    actively bidirectionally-transcribed/replicating agent), with a small
    amount of assembly/sequencing-error-like divergence between them."""

    def mutate(s: str) -> str:
        chars = list(s)
        for i in range(len(chars)):
            if rng.random() < mutation_rate:
                chars[i] = rng.choice([b for b in _BASES if b != chars[i]])
        return "".join(chars)

    sense = mutate(monomer)
    antisense = revcomp(mutate(monomer))
    return sense, antisense


@dataclass
class CorpusEntry:
    id: str
    seq: str
    is_true_positive: bool
    element_id: str
    kind: str  # e.g. "rod_sas", "coding_decoy", "structured_coding_decoy", "plain_null"


@dataclass
class Corpus:
    entries: list[CorpusEntry] = field(default_factory=list)

    def fasta_records(self):
        from .seqio import Record

        return [Record(e.id, e.seq) for e in self.entries]

    def labels(self) -> dict[str, bool]:
        return {e.id: e.is_true_positive for e in self.entries}


def build_calibration_corpus(
    rng: random.Random,
    n_positive_elements: int = 12,
    n_coding_decoys: int = 60,
    n_structured_coding_decoys: int = 20,
    n_plain_nulls: int = 60,
    replicate_range: tuple[int, int] = (2, 4),
    circularize_frac: float = 0.5,
    background_gc: float = 0.5,
) -> Corpus:
    """Build a labeled synthetic contig pool mimicking one metatranscriptome
    assembly's contig set, containing:

    - true positives: rod-like structured elements, each represented by
      2-4 sense/antisense contigs (simulating real sequencing depth of an
      actively-transcribed/replicating element), a fraction of which are
      also circularized.
    - coding decoys: strong single-ORF sequences, unstructured otherwise
      (should score low on S1 and low orphan/S2 -> correctly rejected).
    - structured-coding decoys: an ORF fused to a stable hairpin (tests
      that S2 correctly suppresses real coding+structured sequences, i.e.
      structure alone must not be sufficient).
    - plain nulls: dinucleotide-shuffled Markov background, single copy,
      no strand partner, no structure (baseline).
    """
    corpus = Corpus()

    for i in range(n_positive_elements):
        eid = f"POS{i:03d}"
        n_hairpins = rng.randint(3, 5)
        monomer = make_rod(rng, n_hairpins=n_hairpins, stem_len=10, loop_len=4, linker_len=2)
        n_replicates = rng.randint(*replicate_range)
        circularize = rng.random() < circularize_frac
        for r in range(n_replicates):
            sense, antisense = make_sas_pair(monomer, rng)
            for strand_seq, pol in ((sense, "sense"), (antisense, "antisense")):
                contig = (
                    make_circular_contig(strand_seq, 0.08, rng) if circularize else strand_seq
                )
                cid = f"{eid}_rep{r}_{pol}"
                corpus.entries.append(
                    CorpusEntry(cid, contig, True, eid, "rod_sas_circular" if circularize else "rod_sas")
                )

    for i in range(n_coding_decoys):
        cid = f"CODE{i:03d}"
        length_codons = rng.randint(30, 90)
        seq = make_orf_like(rng, n_codons=length_codons)
        corpus.entries.append(CorpusEntry(cid, seq, False, cid, "coding_decoy"))

    for i in range(n_structured_coding_decoys):
        cid = f"STRCODE{i:03d}"
        orf = make_orf_like(rng, n_codons=rng.randint(30, 70))
        hairpin = make_stem_loop(rng, stem_len=12, loop_len=5)
        seq = orf + hairpin if rng.random() < 0.5 else hairpin + orf
        corpus.entries.append(CorpusEntry(cid, seq, False, cid, "structured_coding_decoy"))

    for i in range(n_plain_nulls):
        cid = f"NULL{i:03d}"
        length = rng.randint(120, 400)
        seq = markov_background(length, gc_content=background_gc, rng=rng)
        corpus.entries.append(CorpusEntry(cid, seq, False, cid, "plain_null"))

    rng.shuffle(corpus.entries)
    return corpus
