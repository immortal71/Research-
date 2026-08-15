"""Tests for the abundance-targeted de Bruijn assembler.

The synthetic cases here mirror the real one that drove the design: a small
circular molecule sampled at high coverage, carrying substitutions across a
large fraction of the population. That is what Peach latent mosaic viroid
looks like in SRR11060618, and it is the case strict unitig extension cannot
handle.
"""
import random

import pytest

from rnasig.assemble import (
    assemble,
    build_unitigs,
    count_kmers,
    greedy_contigs,
    remove_tips,
)
from rnasig.circularity import find_circularity


def _circular_reads(monomer, n=4000, read_len=120, variant_rate=0.0, seed=3):
    """Sample reads uniformly around a circular template."""
    rng = random.Random(seed)
    doubled = monomer + monomer
    reads = []
    for _ in range(n):
        start = rng.randrange(len(monomer))
        read = list(doubled[start : start + read_len])
        if variant_rate and rng.random() < variant_rate:
            for _ in range(2):
                pos = rng.randrange(len(read))
                read[pos] = rng.choice("ACGT")
        reads.append("".join(read))
    return reads


def _random_seq(n, seed=0):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_min_count_discards_the_low_abundance_background():
    monomer = _random_seq(337, seed=1)
    rng = random.Random(2)
    reads = _circular_reads(monomer, n=3000)
    reads += ["".join(rng.choice("ACGT") for _ in range(120)) for _ in range(3000)]

    graph, n_reads, total = count_kmers(reads, k=31, min_count=5)
    assert n_reads == 6000
    # only the circular template survives; the background is all singletons
    assert len(graph) == 337
    assert total > 100_000


def test_assembles_a_clean_circular_molecule_and_recovers_its_unit():
    monomer = _random_seq(337, seed=4)
    result = assemble(_circular_reads(monomer, n=3000), k=31, min_count=5, min_contig_len=200)

    assert result.stats.n_contigs == 1
    circ = find_circularity(result.contigs[0].seq, k=12, max_mismatch=1)
    assert circ.is_circular
    assert circ.unit_length == 337
    assert circ.monomer in monomer + monomer


def test_greedy_extension_recovers_a_quasispecies_that_unitigs_cannot():
    """The PLMVd case, in miniature.

    With variants across half the population the graph becomes a nested
    bubble field. Strict unitig extension shreds it; greedy extension follows
    the dominant path straight through.
    """
    monomer = _random_seq(337, seed=5)
    reads = _circular_reads(monomer, n=6000, variant_rate=0.5, seed=6)

    graph, _, _ = count_kmers(reads, k=25, min_count=5)
    remove_tips(graph, 25)
    unitigs = [u for u in build_unitigs(dict(graph), 25) if len(u) >= 200]
    assert unitigs == []  # this is why greedy exists

    greedy = [c for c in greedy_contigs(dict(graph), 25) if len(c) >= 200]
    assert greedy
    circ = find_circularity(max(greedy, key=len), k=12, max_mismatch=1)
    assert circ.is_circular
    assert circ.unit_length == 337


def test_greedy_seeds_from_the_deepest_kmer_first():
    """An abundant molecule must not be truncated by a shallower one."""
    deep = _random_seq(400, seed=7)
    shallow = _random_seq(400, seed=8)
    reads = _circular_reads(deep, n=4000, seed=9) + _circular_reads(shallow, n=200, seed=10)

    result = assemble(reads, k=31, min_count=5, min_contig_len=200)
    lengths = {len(c.seq) for c in result.contigs}
    assert any(l >= 400 for l in lengths)
    circs = [find_circularity(c.seq, k=12, max_mismatch=1) for c in result.contigs]
    assert any(c.is_circular and c.unit_length == 400 for c in circs)


def test_two_molecules_are_not_fused_into_one_contig():
    a = _random_seq(500, seed=11)
    b = _random_seq(500, seed=12)
    reads = _circular_reads(a, n=2000, seed=13) + _circular_reads(b, n=2000, seed=14)

    result = assemble(reads, k=31, min_count=5, min_contig_len=200)
    # each molecule should come back at roughly its own length, not glued
    assert all(len(c.seq) < 900 for c in result.contigs)
    assert len(result.contigs) >= 2


def test_max_table_prunes_without_losing_the_abundant_molecule():
    monomer = _random_seq(337, seed=15)
    rng = random.Random(16)
    reads = _circular_reads(monomer, n=3000, seed=17)
    reads += ["".join(rng.choice("ACGT") for _ in range(120)) for _ in range(5000)]

    result = assemble(reads, k=31, min_count=5, min_contig_len=200, max_table=50_000)
    circs = [find_circularity(c.seq, k=12, max_mismatch=1) for c in result.contigs]
    assert any(c.is_circular and c.unit_length == 337 for c in circs)


def test_reads_with_ambiguous_bases_do_not_create_kmers():
    graph, _, _ = count_kmers(["ACGT" * 10, "ACGN" * 10], k=31, min_count=1)
    assert all("N" not in kmer for kmer in graph)


def test_empty_input_returns_an_empty_assembly():
    result = assemble([], k=31, min_count=5)
    assert result.contigs == []
    assert result.stats.n_contigs == 0
    assert "contigs=0" in result.stats.summary()


def test_contig_headers_carry_coverage_the_sweep_can_read():
    """phase3_sweep reads multi= out of assembler headers."""
    monomer = _random_seq(337, seed=18)
    result = assemble(_circular_reads(monomer, n=3000, seed=19), k=31, min_count=5)
    header = result.contigs[0].id
    assert "multi=" in header and "len=" in header
    assert float(header.split("multi=")[1].split()[0]) > 100


@pytest.mark.xfail(
    reason="known limitation: greedy extension cannot resolve shared sequence. "
           "At a node two molecules share, the deeper one wins and the shallower "
           "path is hijacked into it. Fixing this needs paired-end links or "
           "coverage-consistency checks, neither of which this assembler has.",
    strict=True,
)
def test_two_molecules_sharing_a_repeat_both_assemble():
    """Documents what is still broken, not what was fixed.

    Removing the global `visited` set from extension stopped truncation at
    already-claimed nodes, which is a real bug and is fixed. It does not
    solve the harder case below: two circles sharing a stretch, where greedy
    follows whichever branch is deeper and the shallow molecule is lost.

    Kept as a strict xfail so that if a future change does resolve shared
    sequence, this fails loudly and gets promoted to a passing test.
    """
    shared = _random_seq(60, seed=41)
    a = _random_seq(420, seed=42) + shared
    b = _random_seq(400, seed=43) + shared

    reads = _circular_reads(a, n=6000, seed=44) + _circular_reads(b, n=1200, seed=45)
    result = assemble(reads, k=31, min_count=5, min_contig_len=200)

    circs = [find_circularity(c.seq, k=12, max_mismatch=1) for c in result.contigs]
    units = {c.unit_length for c in circs if c.is_circular}
    assert len(a) in units, f"deep molecule lost; units={units}"
    assert len(b) in units, f"shallow molecule truncated by the deep one; units={units}"


def test_contained_duplicates_are_collapsed():
    """Unblocking extension lets seeds re-trace a molecule; dedup handles it."""
    monomer = _random_seq(500, seed=46)
    result = assemble(_circular_reads(monomer, n=4000, seed=47), k=31, min_count=5,
                      min_contig_len=200)
    seqs = [c.seq for c in result.contigs]
    for i, s in enumerate(seqs):
        for j, t in enumerate(seqs):
            if i != j:
                assert s not in t, "a contig is contained in another; dedup failed"
