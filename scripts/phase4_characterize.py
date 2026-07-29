#!/usr/bin/env python3
"""Phase 4: characterize a flagged candidate in depth.

Given a FASTA of one representative sequence (typically the centroid of a
cluster flagged significant in Phase 3), produce:
  - length, GC%, MFE structure + dot-bracket, structure z-score (with more
    shuffles than the screening pass, for a tighter estimate)
  - ORF map across all 6 frames
  - orphan-score breakdown
  - circularity / concatemer status
  - an SVG rendering of the predicted secondary structure

Usage:
    python scripts/phase4_characterize.py --input candidate.fasta --outdir results/characterize_X
"""
from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import RNA

from rnasig.seqio import read_fasta, gc_content, rna
from rnasig.circularity import find_circularity, resolve_concatemer
from rnasig.structure import structure_zscore
from rnasig.orphan import orphan_score, longest_orf_nt, codon_bias_chisq


def _orf_map(seq: str) -> list[dict]:
    from rnasig.seqio import dna, revcomp

    seq = dna(seq)
    orfs = []
    for strand, strand_seq in (("+", seq), ("-", revcomp(seq))):
        for frame in range(3):
            codons_seq = strand_seq[frame:]
            codons = [codons_seq[i : i + 3] for i in range(0, len(codons_seq) - 2, 3)]
            start_idx = None
            for i, c in enumerate(codons):
                if start_idx is None and c == "ATG":
                    start_idx = i
                elif start_idx is not None and c in {"TAA", "TAG", "TGA"}:
                    orfs.append({
                        "strand": strand,
                        "frame": frame,
                        "nt_start": frame + start_idx * 3,
                        "nt_end": frame + (i + 1) * 3,
                        "length_aa": i - start_idx,
                    })
                    start_idx = None
    orfs.sort(key=lambda o: o["length_aa"], reverse=True)
    return orfs


def characterize(rec_id: str, seq: str, n_shuffles: int = 200, rng=None) -> dict:
    rng = rng or random.Random(0)
    circ = find_circularity(seq)
    circ = resolve_concatemer(circ)
    struct = structure_zscore(seq, n_shuffles=n_shuffles, rng=rng)
    orphan = orphan_score(seq)
    orfs = _orf_map(seq)

    return {
        "id": rec_id,
        "length": len(seq),
        "gc_content": gc_content(seq),
        "circular": circ.is_circular,
        "unit_length": circ.unit_length,
        "n_copies": circ.n_copies,
        "mfe_energy": struct.energy,
        "mfe_structure": struct.structure,
        "structure_z_score": struct.z_score,
        "structure_null_mean": struct.null_mean,
        "structure_null_std": struct.null_std,
        "n_shuffles": struct.n_shuffles,
        "orphan_score": orphan.orphan_score,
        "orf_coverage": orphan.orf_coverage,
        "longest_orf_nt": orphan.longest_orf_nt,
        "codon_bias_chisq": orphan.codon_bias_chisq,
        "top_orfs": orfs[:5],
        "n_orfs_total": len(orfs),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=str, required=True)
    ap.add_argument("--outdir", type=str, required=True)
    ap.add_argument("--n-shuffles", type=int, default=200)
    args = ap.parse_args()

    records = read_fasta(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    reports = []
    for rec in records:
        report = characterize(rec.id, rec.seq, n_shuffles=args.n_shuffles)
        reports.append(report)

        safe_id = "".join(c if c.isalnum() else "_" for c in rec.id)[:80]
        svg_path = outdir / f"{safe_id}_structure.svg"
        try:
            RNA.svg_rna_plot(rna(rec.seq), report["mfe_structure"], str(svg_path))
        except Exception as e:
            print(f"  (structure plot failed for {rec.id}: {e})")

        print(f"{rec.id}: len={report['length']} gc={report['gc_content']:.2f} "
              f"z={report['structure_z_score']:.2f} orphan={report['orphan_score']:.2f} "
              f"circular={report['circular']} longest_orf_nt={report['longest_orf_nt']}")

    (outdir / "characterization.json").write_text(json.dumps(reports, indent=2))
    print(f"\nwrote {outdir / 'characterization.json'}")


if __name__ == "__main__":
    main()
