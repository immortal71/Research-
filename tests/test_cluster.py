import random

from rnasig.cluster import cluster_sequences, sas_clusters, best_circular_identity
from rnasig.seqio import Record, revcomp
from rnasig.simulate import make_rod, make_sas_pair


def test_identical_sequence_gets_high_identity():
    rng = random.Random(4)
    seq = make_rod(rng, n_hairpins=5)
    identity, orientation = best_circular_identity(seq, seq)
    assert identity > 0.95
    assert orientation == "sense"


def test_antisense_detected_via_revcomp():
    rng = random.Random(5)
    seq = make_rod(rng, n_hairpins=5)
    identity, orientation = best_circular_identity(revcomp(seq), seq)
    assert identity > 0.95
    assert orientation == "antisense"


def test_circular_permutation_still_clusters():
    rng = random.Random(6)
    seq = make_rod(rng, n_hairpins=6)
    rotated = seq[100:] + seq[:100]
    identity, _ = best_circular_identity(rotated, seq)
    assert identity > 0.9


def test_sas_cluster_detected_for_replicating_element():
    rng = random.Random(7)
    monomer = make_rod(rng, n_hairpins=6)
    records = []
    for r in range(3):
        sense, antisense = make_sas_pair(monomer, rng, mutation_rate=0.02)
        records.append(Record(f"rep{r}_sense", sense))
        records.append(Record(f"rep{r}_antisense", antisense))

    clusters = cluster_sequences(records, id_threshold=0.7)
    sas = sas_clusters(clusters)
    assert len(sas) >= 1
    assert len(sas[0].members) >= 1


def test_unrelated_sequences_do_not_cluster():
    rng = random.Random(8)
    a = make_rod(rng, n_hairpins=5)
    b = make_rod(rng, n_hairpins=5)
    records = [Record("a", a), Record("b", b)]
    clusters = cluster_sequences(records, id_threshold=0.7)
    assert len(clusters) == 2
