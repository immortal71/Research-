#!/usr/bin/env python3
"""Phase 6: sweep public metatranscriptomes for circular RNA candidates.

This is the run the repo never actually did. Phases 1 to 5 built and
validated the parts; this points them at data chosen on purpose rather than
on availability, and it is the only script here that starts from an
accession instead of a file someone already downloaded.

For each accession it:

  1. preflights the library type (srameta) and skips anything that is not
     shotgun RNA,
  2. streams a read subset straight from ENA without saving the full run,
     which matters because these files are 5-15 GB each,
  3. assembles the abundant fraction (assemble),
  4. keeps contigs that look circular in the size range where viroids and
     obelisks live, then scores each one's structure and coding potential
     and shortlists only those folding better than shuffled sequence.

Survivors are written per-run for `phase5_identify.py --input` to resolve.
Nothing here decides that anything is novel; that takes a database, and
Phase 5 is where it happens.

Usage:
    python scripts/phase6_hunt.py --accessions accs.txt --outdir results/hunt
    python scripts/phase6_hunt.py --accessions SRR4436415 --reads 400000
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rnasig.assemble import assemble
from rnasig.circularity import find_circularity
from rnasig.orphan import orphan_score
from rnasig.rodlike import rod_profile
from rnasig.rrna_kmer import build_reference_kmers, fetch_reference, rrna_kmer_fraction
from rnasig.seqio import Record, gc_content, write_fasta
from rnasig.srameta import assess_library, fetch_run_metadata
from rnasig.structure import structure_zscore

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"

# Viroids run ~250-400 nt, obelisks ~1 kb. Outside this band a circular
# signal is far more likely to be a repeat or an assembly artifact.
MIN_UNIT = 150
MAX_UNIT = 2500

# Circularity on its own is close to worthless as evidence, which the first
# version of this script got wrong. A de Bruijn graph cycles wherever a
# low-complexity repeat does, so an AT-rich tandem repeat looks exactly like
# a circular molecule to the terminal-repeat test. The first hot-springs run
# swept here returned six "circular" contigs scoring structure z of -0.55 to
# 1.51, i.e. no more structured than shuffled sequence.
#
# Structure is what separates them. PLMVd, the one known replicon this
# pipeline has recovered from real data, scores z=17.13. A candidate has to
# clear this bar before it is worth a database lookup.
MIN_STRUCT_Z = 3.0

# Structure alone was still not enough. Two further filters, both added
# after real candidates exposed the gap:
#
#   Rod-likeness. SRR5949183_contig_204 scored z=5.05 on one 17 bp hairpin
#   with 44% of its bases paired. Viroids are rods; PLMVd is 68% paired at
#   -0.473 kcal/mol/nt. See rnasig.rodlike.
#
#   rRNA. SRR19432462 yielded a 342 nt contig at 6752x coverage that passed
#   structure and rod-likeness and is 28S rRNA at 100% identity. Ribosomal
#   RNA is structured, abundant, and assembles with terminal repeats, so it
#   clears every axis of this signature at once. See rnasig.rrna_kmer.
MAX_RRNA_FRACTION = 0.10
RRNA_CACHE = str(ROOT / "data" / "reference" / "rrna" / "rrna_reference.fasta")


def fastq_urls(accession: str, timeout: float = 30.0) -> list[str]:
    query = f"?accession={accession}&result=read_run&fields=fastq_ftp&format=tsv"
    try:
        with urllib.request.urlopen(ENA_PORTAL + query, timeout=timeout) as resp:
            text = resp.read().decode("utf-8", "replace")
    except (urllib.error.URLError, OSError):
        return []
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if len(lines) < 2:
        return []
    paths = lines[1].split("\t")[-1].split(";")
    return [f"https://{p}" for p in paths if p.strip()]


def stream_reads(url: str, limit: int, timeout: float = 300.0, min_useful: int = 100_000):
    """Yield sequences from a remote gzipped FASTQ, stopping at limit.

    The connection is closed as soon as enough reads have been seen, so a
    12 GB run costs only the first few hundred MB of transfer.

    A deep subset means a transfer of hundreds of MB, and those drop: the
    server closes early, gzip hits the end of its data without a trailer,
    and raises. Treating that as a failed run throws away everything already
    received, which on the first deep run was most of a 280 MB download. So
    a truncated stream ends the iteration instead, and the caller assembles
    whatever arrived. Anything shorter than min_useful is re-raised, since
    that is a genuinely failed fetch rather than a short one.
    """
    req = urllib.request.Request(url, headers={"User-Agent": "rnasig phase6"})
    seen = 0
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        with gzip.GzipFile(fileobj=resp) as gz:
            text = io.TextIOWrapper(gz, encoding="utf-8", errors="replace")
            try:
                for i, line in enumerate(text):
                    if i % 4 == 1:
                        seen += 1
                        yield line.strip()
                        if seen >= limit:
                            return
            except (EOFError, OSError) as exc:
                if seen < min_useful:
                    raise
                print(f"      (stream truncated after {seen} reads: {exc}; assembling those)")
                return


def hunt_one(
    accession: str,
    n_reads: int,
    k: int,
    min_count: int,
    outdir: Path,
    max_table: int,
    n_shuffles: int = 200,
    rrna_ref: set[str] | None = None,
) -> dict:
    record: dict = {"accession": accession}
    meta = fetch_run_metadata(accession)
    verdict = assess_library(meta)
    record["library"] = meta.summary
    record["compatible"] = verdict.compatible
    if not verdict.compatible:
        record["skipped"] = verdict.reasons
        print(f"  SKIP {accession}: {verdict.reasons[0][:80]}")
        return record
    record["sample"] = meta.scientific_name

    urls = fastq_urls(accession)
    if not urls:
        record["error"] = "no fastq url"
        print(f"  SKIP {accession}: no FASTQ url")
        return record

    # Stream straight into the counter rather than materialising the reads.
    # A 3M-read subset held as a list is several hundred MB before the k-mer
    # table is even allocated, and deep subsets are the whole point: a
    # 300k-read slice of a 100M-read run only reaches elements above ~2000x,
    # which is far above where a replicon is likely to sit.
    t0 = time.time()
    try:
        asm = assemble(
            stream_reads(urls[0], n_reads),
            k=k, min_count=min_count, min_contig_len=MIN_UNIT, max_table=max_table,
        )
    except (urllib.error.URLError, OSError, EOFError) as exc:
        record["error"] = f"stream failed: {exc}"
        print(f"  FAIL {accession}: stream failed ({exc})")
        return record
    t_asm = time.time() - t0
    record["stats"] = asm.stats.summary()

    hits = []
    for contig in asm.contigs:
        circ = find_circularity(contig.seq, k=12, max_mismatch=1)
        if not circ.is_circular or circ.unit_length is None:
            continue
        if not (MIN_UNIT <= circ.unit_length <= MAX_UNIT):
            continue
        monomer = circ.monomer or contig.seq
        stability = structure_zscore(monomer, n_shuffles=n_shuffles)
        orphan = orphan_score(monomer)
        rod = rod_profile(monomer)
        rrna = rrna_kmer_fraction(monomer, rrna_ref) if rrna_ref else 0.0
        cov = float(contig.id.split("multi=")[1].split()[0])
        hits.append(
            {
                "contig": contig.id.split()[0],
                "length": len(contig.seq),
                "unit_length": circ.unit_length,
                "n_copies": round(circ.n_copies or 0, 2),
                "coverage": round(cov, 1),
                "gc": round(gc_content(monomer), 3),
                "struct_z": round(stability.z_score, 2),
                "orphan": round(orphan.orphan_score, 2),
                "paired_fraction": round(rod.paired_fraction, 3),
                "mfe_per_nt": round(rod.mfe_per_nt, 3),
                "rodlike": rod.is_rodlike,
                "rrna_fraction": round(rrna, 3),
                "longest_orf_nt": orphan.longest_orf_nt,
                "monomer": monomer,
            }
        )

    hits.sort(key=lambda h: -h["struct_z"])
    shortlist = [h for h in hits if h["struct_z"] >= MIN_STRUCT_Z
                 and h["rodlike"] and h["rrna_fraction"] < MAX_RRNA_FRACTION]
    record["n_circular"] = len(hits)
    record["n_structured"] = len(shortlist)
    record["circular"] = [{kk: v for kk, v in h.items() if kk != "monomer"} for h in hits[:25]]
    record["timing"] = {"stream_and_assemble_s": round(t_asm)}

    print(
        f"  {accession} [{meta.scientific_name}] {asm.stats.n_reads} reads "
        f"-> {asm.stats.n_contigs} contigs, {len(hits)} circular, "
        f"{len(shortlist)} with z>={MIN_STRUCT_Z} "
        f"({t_asm:.0f}s stream+asm)"
    )
    for h in hits[:5]:
        if h["rrna_fraction"] >= MAX_RRNA_FRACTION:
            mark = "  (rRNA)"
        elif not h["rodlike"]:
            mark = "  (not rod-like)"
        elif h["struct_z"] < MIN_STRUCT_Z:
            mark = ""
        else:
            mark = "  <== SHORTLIST"
        print(
            f"      unit={h['unit_length']:5d} cov={h['coverage']:8.1f} GC={h['gc']:.2f} "
            f"z={h['struct_z']:6.2f} paired={h['paired_fraction']:.0%} "
            f"rRNA={h['rrna_fraction']:.2f}{mark}"
        )

    if shortlist:
        outdir.mkdir(parents=True, exist_ok=True)
        recs = [
            Record(f"{accession}_{h['contig']} unit={h['unit_length']} cov={h['coverage']}", h["monomer"])
            for h in shortlist
            if h["monomer"]
        ]
        write_fasta(str(outdir / f"{accession}_circular.fasta"), recs)
    return record


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accessions", required=True,
                    help="file of run accessions (one per line), or a single accession")
    ap.add_argument("--reads", type=int, default=300_000, help="reads to stream per run")
    ap.add_argument("--k", type=int, default=25)
    ap.add_argument("--min-count", type=int, default=5)
    ap.add_argument("--max-table", type=int, default=12_000_000)
    ap.add_argument("--shuffles", type=int, default=200,
                    help="shuffles per structure z-score")
    ap.add_argument("--outdir", default="results/hunt")
    args = ap.parse_args()

    path = Path(args.accessions)
    accessions = ([ln.strip() for ln in path.read_text().splitlines() if ln.strip()]
                  if path.exists() else [args.accessions])

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"hunting {len(accessions)} run(s), {args.reads} reads each, k={args.k} min_count={args.min_count}\n")

    rrna_ref = build_reference_kmers(fetch_reference(RRNA_CACHE) or "")
    print(f"rRNA reference: {len(rrna_ref)} k-mers"
          if rrna_ref else "rRNA reference unavailable; that filter is off")

    results = []
    for i, acc in enumerate(accessions, 1):
        print(f"[{i}/{len(accessions)}]")
        try:
            results.append(hunt_one(acc, args.reads, args.k, args.min_count, outdir,
                                    args.max_table, args.shuffles, rrna_ref))
        except KeyboardInterrupt:
            raise
        except Exception as exc:  # a bad run must not end the sweep
            print(f"  ERROR {acc}: {type(exc).__name__}: {exc}")
            results.append({"accession": acc, "error": f"{type(exc).__name__}: {exc}"})
        (outdir / "hunt_results.json").write_text(json.dumps(results, indent=2))

    total = sum(r.get("n_circular", 0) for r in results)
    print(f"\ndone: {total} circular candidates across {len(results)} runs")
    print(f"wrote {outdir / 'hunt_results.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
