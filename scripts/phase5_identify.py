#!/usr/bin/env python3
"""Phase 5: find out what a candidate actually is.

Phases 1-4 nominate. Nothing in them can tell you whether a contig is
undescribed or merely unfamiliar, because that question is answered by
databases and the earlier work had no route to one. This script is that
route. It does two separable jobs:

  --run <ACCESSION>
      Preflight a sequencing run before you sweep it. Checks the archive's
      library metadata (srameta) and, if raw reads are supplied, checks the
      reads themselves for primer dominance (ampliconqc). Either check can
      fail the run on its own.

  --input <FASTA>
      Resolve candidates against core_nt, nr and Rfam (dbcheck), and print
      a one-line verdict per sequence.

Running the preflight against SRR13291825 is what this repo should have
done before Phase 3, and is why the module exists:

    python scripts/phase5_identify.py --run SRR13291825

Usage:
    python scripts/phase5_identify.py --run SRR13291825
    python scripts/phase5_identify.py --run SRR13291825 --reads r1.fastq.gz
    python scripts/phase5_identify.py --input candidates.fasta --outdir results/identify
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rnasig.ampliconqc import profile_fastq
from rnasig.dbcheck import identify
from rnasig.seqio import read_fasta
from rnasig.srameta import assess_library, fetch_run_metadata, to_json


def run_preflight(accession: str, reads: list[str] | None, outdir: Path | None) -> int:
    print(f"== library preflight: {accession} ==")
    meta = fetch_run_metadata(accession)
    verdict = assess_library(meta)

    if not meta.fetched:
        print("  metadata: could not be retrieved (offline or unknown accession)")
    else:
        print(f"  {meta.summary}")
        if meta.read_count:
            print(f"  reads: {meta.read_count}")

    if verdict.compatible:
        print("  VERDICT: compatible with the SOS signature")
    else:
        label = "UNKNOWN" if verdict.unknown else "INCOMPATIBLE"
        print(f"  VERDICT: {label} -- do not sweep this run")
        for reason in verdict.reasons:
            print(f"    - {reason}")

    payload = json.loads(to_json(meta, verdict))

    if reads:
        payload["read_profiles"] = []
        for path in reads:
            print(f"\n== read profile: {path} ==")
            profile = profile_fastq(path)
            print("  " + profile.describe().replace("\n", "\n  "))
            payload["read_profiles"].append(
                {
                    "path": path,
                    "n_reads": profile.n_reads,
                    "top_fraction": profile.top_fraction,
                    "is_amplicon": profile.is_amplicon,
                    "top_prefixes": profile.top_prefixes,
                }
            )
            if profile.is_amplicon:
                print("  VERDICT: primer-dominated -- targeted library, not shotgun RNA")

    if outdir:
        outdir.mkdir(parents=True, exist_ok=True)
        path = outdir / f"preflight_{accession}.json"
        path.write_text(json.dumps(payload, indent=2))
        print(f"\nwrote {path}")

    incompatible = not verdict.compatible or any(
        p["is_amplicon"] for p in payload.get("read_profiles", [])
    )
    return 1 if incompatible else 0


def run_identify(fasta: str, outdir: Path | None, skip: set[str]) -> int:
    records = read_fasta(fasta)
    print(f"== identifying {len(records)} sequence(s) from {fasta} ==")
    print("   (each query is a live search against NCBI/Rfam; expect a few minutes each)\n")

    results = []
    for record in records:
        name = record.id.split()[0]
        print(f"-- {name} ({len(record.seq)} nt)")
        ident = identify(
            record,
            do_nt="nt" not in skip,
            do_protein="protein" not in skip,
            do_rfam="rfam" not in skip,
        )
        print(f"   {ident.verdict()}")
        for hit in ident.protein_hits[:3]:
            print(f"     protein: {hit.describe()}")
        for hit in ident.nt_hits[:3]:
            print(f"     nt:      {hit.describe()}")
        for err in ident.errors:
            print(f"     ! {err}")
        results.append(ident.to_dict())
        print()

    if outdir:
        outdir.mkdir(parents=True, exist_ok=True)
        path = outdir / "identification.json"
        path.write_text(json.dumps({"source": fasta, "results": results}, indent=2))
        print(f"wrote {path}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", help="SRA/ENA run accession to preflight before sweeping")
    ap.add_argument("--reads", nargs="*", help="FASTQ(.gz) files for reference-free amplicon detection")
    ap.add_argument("--input", help="FASTA of candidates to resolve against nt/nr/Rfam")
    ap.add_argument("--outdir", help="directory for JSON output")
    ap.add_argument("--skip", nargs="*", default=[], choices=["nt", "protein", "rfam"],
                    help="database searches to skip")
    args = ap.parse_args()

    if not args.run and not args.input:
        ap.error("give --run (preflight a run) or --input (resolve candidates)")

    outdir = Path(args.outdir) if args.outdir else None
    status = 0
    if args.run:
        status |= run_preflight(args.run, args.reads, outdir)
    if args.input:
        status |= run_identify(args.input, outdir, set(args.skip))
    return status


if __name__ == "__main__":
    raise SystemExit(main())
