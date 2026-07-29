#!/usr/bin/env python3
"""Phase 3: sweep an environment for the novel signature.

Usage (real data, once you have assembled contigs from an environment of
interest -- e.g. rnaSPAdes output from soil, hydrothermal vent, insect gut,
or any other under-sampled metatranscriptome):

    python scripts/phase3_sweep.py --input contigs.fasta --outdir results/sweep_soil

Usage (demo mode -- see docs/LIMITATIONS.md for why this is a *labeled
synthetic* stand-in rather than a real environmental sweep in this
sandboxed session: NCBI/SRA/ENA/EBI/OUP/PMC/NERSC/Zenodo/OSF/Dryad/figshare
were all unreachable; only github.com/raw.githubusercontent.com and bare
cloud-storage roots were):

    python scripts/phase3_sweep.py --demo --outdir results/sweep_demo

This script does NOT invent findings -- it runs the identical pipeline
exercised and calibrated in Phase 1/2 on whatever FASTA it is pointed at,
and reports whatever it finds (including "nothing survived FDR control",
which is itself the honest, expected answer for a demo/null-background
run).
"""
from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rnasig.seqio import read_fasta
from rnasig.cluster import cluster_sequences, non_singleton_clusters
from rnasig.signature import score_cluster, benjamini_hochberg
from rnasig.simulate import build_calibration_corpus
from rnasig.preprocess import filter_contigs


def run_sweep(records, outdir: Path, n_shuffles=100, id_threshold=0.7, alpha=0.05, max_len=1000, min_len=10, skip_preprocess=False):
    rng = random.Random(0)
    filtered = [r for r in records if min_len <= len(r) <= max_len]
    if skip_preprocess:
        pre_report = None
    else:
        filtered, pre_report = filter_contigs(filtered)
        print(f"preprocess: {pre_report.n_input} -> {pre_report.n_kept} contigs "
              f"({pre_report.n_dropped_adapter} adapter, "
              f"{pre_report.n_dropped_homopolymer} homopolymer)")
    clusters = cluster_sequences(filtered, id_threshold=id_threshold)

    candidates = []
    for i, cl in enumerate(clusters):
        cand = score_cluster(cl, cluster_id=f"cluster_{i}", n_shuffles=n_shuffles, rng=rng)
        candidates.append((cl, cand))

    pvals = [c.p_value for _, c in candidates]
    reject, qvals = benjamini_hochberg(pvals, alpha=alpha)

    hits = []
    for (cl, cand), r, q in zip(candidates, reject, qvals):
        member_ids = [cl.centroid.id] + [m.id for m in cl.members]
        hits.append({
            "cluster_id": cand.cluster_id,
            "representative_id": cand.representative_id,
            "cluster_size": len(member_ids),
            "member_ids": member_ids,
            "length": cand.length,
            "structure_z": cand.structure_z,
            "orphan_score": cand.orphan_score,
            "has_sas": cand.has_sas,
            "is_circular": cand.is_circular,
            "combined_score": cand.combined_score,
            "p_value": cand.p_value,
            "q_value": q,
            "significant_at_alpha": bool(r),
        })

    hits.sort(key=lambda h: h["combined_score"], reverse=True)

    outdir.mkdir(parents=True, exist_ok=True)
    summary = {
        "n_input_contigs": len(records),
        "n_length_filtered": len(filtered),
        "n_clusters": len(clusters),
        "n_nonsingleton_clusters": len(non_singleton_clusters(clusters)),
        "alpha": alpha,
        "n_significant": sum(reject),
        "hits": hits,
        "preprocess": {
            "enabled": pre_report is not None,
            "n_dropped_adapter": pre_report.n_dropped_adapter if pre_report else 0,
            "n_dropped_homopolymer": pre_report.n_dropped_homopolymer if pre_report else 0,
            "dropped_reasons": pre_report.dropped_reasons if pre_report else {},
        },
    }
    (outdir / "sweep_results.json").write_text(json.dumps(summary, indent=2))
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=str, help="FASTA file of assembled contigs to sweep")
    ap.add_argument("--demo", action="store_true", help="run on a labeled synthetic demo corpus instead")
    ap.add_argument("--outdir", type=str, required=True)
    ap.add_argument("--n-shuffles", type=int, default=100)
    ap.add_argument("--id-threshold", type=float, default=0.7)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--skip-preprocess", action="store_true",
                    help="skip adapter/homopolymer contig filter (default: filter on -- see src/rnasig/preprocess.py)")
    args = ap.parse_args()

    outdir = Path(args.outdir)

    if args.demo:
        print("DEMO MODE: input is a labeled SYNTHETIC corpus, not real environmental "
              "sequence data -- see docs/LIMITATIONS.md for why. Ground truth labels are "
              "written alongside results for a sanity check, but a real sweep's output "
              "would NOT include a labels file (there is no ground truth for real data).")
        rng = random.Random(42)
        corpus = build_calibration_corpus(
            rng,
            n_positive_elements=3,
            n_coding_decoys=25,
            n_structured_coding_decoys=8,
            n_plain_nulls=40,
            replicate_range=(1, 2),  # sparse, low-depth: mimics a real rare-element scenario
            circularize_frac=0.6,
        )
        records = corpus.fasta_records()
        outdir.mkdir(parents=True, exist_ok=True)
        (outdir / "demo_ground_truth_labels.json").write_text(
            json.dumps(corpus.labels(), indent=2)
        )
    elif args.input:
        records = read_fasta(args.input)
    else:
        print("must pass --input <fasta> or --demo", file=sys.stderr)
        sys.exit(1)

    summary = run_sweep(
        records, outdir, n_shuffles=args.n_shuffles, id_threshold=args.id_threshold,
        alpha=args.alpha, skip_preprocess=args.skip_preprocess,
    )
    print(f"input contigs: {summary['n_input_contigs']} "
          f"(length-filtered: {summary['n_length_filtered']})")
    print(f"clusters: {summary['n_clusters']} ({summary['n_nonsingleton_clusters']} non-singleton)")
    print(f"significant at alpha={args.alpha}: {summary['n_significant']}")
    for h in summary["hits"][:10]:
        print(f"  {h['representative_id']}: score={h['combined_score']:.2f} "
              f"z={h['structure_z']:.2f} orphan={h['orphan_score']:.2f} "
              f"sas={h['has_sas']} circular={h['is_circular']} "
              f"sig={h['significant_at_alpha']}")
    print(f"\nwrote {outdir / 'sweep_results.json'}")


if __name__ == "__main__":
    main()
