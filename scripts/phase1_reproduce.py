#!/usr/bin/env python3
"""Phase 1: reproduce a known detection.

Runs our reimplementation of VNom's core logic (circularity ->
circular-permutation-aware clustering -> sense/antisense co-occurrence)
against the VNom authors' own real test dataset (SRR11060618_subset.fasta,
38 rnaSPAdes contigs assembled from a real human gut metatranscriptome run)
and reports which contigs are flagged circular and which clusters show
sense+antisense co-occurrence -- the exact "obelisk-like" signature the
original tool was built to find.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rnasig.seqio import read_fasta, spades_length
from rnasig.circularity import find_circularity
from rnasig.cluster import cluster_sequences, non_singleton_clusters, sas_clusters


def main():
    fasta_path = ROOT / "data/reference/vnom/SRR11060618_subset.fasta"
    records = read_fasta(str(fasta_path))
    print(f"loaded {len(records)} contigs from {fasta_path.name}")

    max_len, min_len = 1000, 10
    filtered = [r for r in records if min_len <= len(r) <= max_len]
    print(f"length-filtered ({min_len}-{max_len}nt): {len(filtered)} contigs remain")

    circ_hits = []
    for r in filtered:
        result = find_circularity(r.seq, k=12)
        if result.is_circular:
            circ_hits.append((r, result))

    print(f"\ncircularity: {len(circ_hits)}/{len(filtered)} contigs flagged circular (k=12)")
    for r, result in circ_hits:
        print(
            f"  {r.id}: len={len(r)} unit_length={result.unit_length} "
            f"overhang={result.overhang} n_copies={result.n_copies:.2f}"
        )

    circular_records = [r for r, _ in circ_hits]
    clusters = cluster_sequences(circular_records, id_threshold=0.7)
    nontrivial = non_singleton_clusters(clusters)
    print(f"\nclustering (circular-permutation-aware, id>=0.7): "
          f"{len(clusters)} total clusters, {len(nontrivial)} with >1 member")
    for i, cl in enumerate(nontrivial):
        members = [cl.centroid.id] + [m.id for m in cl.members]
        print(f"  cluster {i}: centroid={cl.centroid.id} (len={len(cl.centroid)}), "
              f"members={len(cl.members)}, orientations={cl.orientations}")

    sas = sas_clusters(clusters)
    print(f"\nsense/antisense (SAS) clusters: {len(sas)}")
    result_summary = {
        "input_contigs": len(records),
        "length_filtered": len(filtered),
        "circular_hits": len(circ_hits),
        "circular_contig_ids": [r.id for r, _ in circ_hits],
        "nonsingleton_clusters": len(nontrivial),
        "sas_clusters": len(sas),
        "sas_cluster_details": [
            {
                "centroid": cl.centroid.id,
                "centroid_length": len(cl.centroid),
                "members": [m.id for m in cl.members],
                "orientations": cl.orientations,
            }
            for cl in sas
        ],
    }

    out_path = ROOT / "results/phase1_reproduction.json"
    out_path.parent.mkdir(exist_ok=True)
    out_path.write_text(json.dumps(result_summary, indent=2))
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()
