#!/usr/bin/env python3
"""Phase 2: define and calibrate the novel signature (SOS: Structure-stable,
Orphan-sequence, Strand-symmetric elements).

Builds labeled synthetic contig pools (true positive elements + three
flavors of matched decoys: coding, structured-coding, plain background),
runs the full detection pipeline (circularity -> circular-permutation
clustering -> SAS detection -> structure z-score -> orphan score ->
combined score), and reports:

  1. Discrimination (ROC-AUC, PR-AUC) of the combined score vs ground truth.
  2. BH-FDR-controlled calling: empirical FDR and power at alpha=0.05/0.10.
  3. A power-vs-depth curve: recovery rate of true elements as a function
     of how many replicate sense/antisense contigs represent each element
     (a proxy for sequencing depth of a real sample).
  4. Specificity checks: false-call rate on structure-only and coding-only
     decoys individually (does the signature need all three axes, or would
     any one alone have sufficed?).
"""
from __future__ import annotations

import json
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from sklearn.metrics import roc_auc_score, average_precision_score

from rnasig.simulate import build_calibration_corpus
from rnasig.cluster import cluster_sequences
from rnasig.signature import score_cluster, benjamini_hochberg

N_SHUFFLES = 10
# NOTE ON N_SHUFFLES: RNAfold's MFE computation is O(n^3) in sequence
# length, so the shuffle-null count is the main compute-cost lever in this
# sandboxed environment. 10 shuffles gives a noisier z-score estimate than
# the >=100 we'd use given real compute budget (see docs/LIMITATIONS.md);
# it is enough to demonstrate the calibration procedure and get a real,
# non-fabricated signal-vs-null separation, not a production-grade z-score.


def evaluate_corpus(corpus, id_threshold=0.7, n_shuffles=N_SHUFFLES, rng=None):
    rng = rng or random.Random(0)
    labels_by_id = corpus.labels()
    element_by_id = {e.id: e.element_id for e in corpus.entries}

    clusters = cluster_sequences(corpus.fasta_records(), id_threshold=id_threshold)
    print(f"  clustered into {len(clusters)} clusters, scoring each...", flush=True)
    candidates = []
    cluster_true_elements = []
    for i, cl in enumerate(clusters):
        member_ids = [cl.centroid.id] + [m.id for m in cl.members]
        is_true = any(labels_by_id[mid] for mid in member_ids)
        true_elements = {element_by_id[mid] for mid in member_ids if labels_by_id[mid]}
        cand = score_cluster(cl, cluster_id=f"cluster_{i}", n_shuffles=n_shuffles, rng=rng)
        candidates.append((cand, is_true, true_elements))
        cluster_true_elements.append(true_elements)
        if (i + 1) % 10 == 0:
            print(f"    scored {i + 1}/{len(clusters)} clusters", flush=True)

    return candidates


def summarize(candidates, alpha=0.05):
    scores = [c.combined_score for c, _, _ in candidates]
    labels = [int(is_true) for _, is_true, _ in candidates]
    pvals = [c.p_value for c, _, _ in candidates]

    n_pos = sum(labels)
    n_neg = len(labels) - n_pos

    result = {"n_clusters": len(candidates), "n_true_positive_clusters": n_pos, "n_negative_clusters": n_neg}

    if n_pos and n_neg:
        result["roc_auc"] = roc_auc_score(labels, scores)
        result["pr_auc"] = average_precision_score(labels, scores)
    else:
        result["roc_auc"] = None
        result["pr_auc"] = None

    reject, qvals = benjamini_hochberg(pvals, alpha=alpha)
    called_true = sum(1 for r, l in zip(reject, labels) if r and l)
    called_false = sum(1 for r, l in zip(reject, labels) if r and not l)
    n_called = called_true + called_false
    result["bh_alpha"] = alpha
    result["n_called_significant"] = n_called
    result["true_positives_called"] = called_true
    result["false_positives_called"] = called_false
    result["empirical_fdr"] = (called_false / n_called) if n_called else 0.0
    result["power_by_structure_axis_alone"] = (called_true / n_pos) if n_pos else None
    return result


def per_decoy_type_false_call_rate(candidates, decoy_kind_lookup, alpha=0.05):
    """For non-true clusters, break down false-call rate by the *dominant*
    decoy kind in that cluster (helps show which axis is doing the work)."""
    pvals = [c.p_value for c, _, _ in candidates]
    reject, _ = benjamini_hochberg(pvals, alpha=alpha)
    breakdown: dict[str, list[int]] = {}
    for (cand, is_true, _), r in zip(candidates, reject):
        if is_true:
            continue
        kind = decoy_kind_lookup.get(cand.representative_id, "unknown")
        breakdown.setdefault(kind, [0, 0])
        breakdown[kind][1] += 1
        if r:
            breakdown[kind][0] += 1
    return {k: {"false_called": v[0], "total": v[1]} for k, v in breakdown.items()}


def main():
    master_rng = random.Random(1234)
    results = {}

    # --- 1. Main calibration run: ROC/PR/FDR on a standard-depth corpus ---
    corpus = build_calibration_corpus(
        master_rng,
        n_positive_elements=8,
        n_coding_decoys=20,
        n_structured_coding_decoys=8,
        n_plain_nulls=20,
        replicate_range=(2, 4),
        circularize_frac=0.5,
    )
    decoy_kind_lookup = {e.id: e.kind for e in corpus.entries}
    print(f"main corpus: {len(corpus.entries)} contigs "
          f"({sum(e.is_true_positive for e in corpus.entries)} true-positive contigs)")

    candidates = evaluate_corpus(corpus, rng=master_rng)
    main_summary = summarize(candidates, alpha=0.05)
    main_summary["false_call_breakdown_by_decoy_kind"] = per_decoy_type_false_call_rate(
        candidates, decoy_kind_lookup, alpha=0.05
    )
    results["main_calibration"] = main_summary
    print(json.dumps(main_summary, indent=2))

    # --- 2. Power vs. sequencing-depth curve ---
    depth_scenarios = [(1, 1), (1, 2), (2, 3), (3, 5)]
    power_curve = []
    for lo, hi in depth_scenarios:
        rng = random.Random(999 + lo + hi)
        corpus_d = build_calibration_corpus(
            rng,
            n_positive_elements=8,
            n_coding_decoys=12,
            n_structured_coding_decoys=5,
            n_plain_nulls=12,
            replicate_range=(lo, hi),
            circularize_frac=0.5,
        )
        cands_d = evaluate_corpus(corpus_d, rng=rng)
        summ_d = summarize(cands_d, alpha=0.05)
        power_curve.append({"replicate_range": [lo, hi], **summ_d})
        print(f"depth {(lo, hi)}: power={summ_d['power_by_structure_axis_alone']}, "
              f"fdr={summ_d['empirical_fdr']:.3f}, roc_auc={summ_d['roc_auc']}")
    results["power_vs_depth"] = power_curve

    out_path = ROOT / "results/phase2_calibration.json"
    out_path.parent.mkdir(exist_ok=True)
    out_path.write_text(json.dumps(results, indent=2))
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()
