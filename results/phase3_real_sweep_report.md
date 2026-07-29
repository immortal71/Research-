# Phase 3 report: real environmental sweep

**Update to the earlier synthetic-only Phase 3.** After the initial writeup
noted that "every standard sequence-archive route was blocked by the
sandbox network policy," I retested more carefully and found the **AWS
Open Data mirror of SRA (`sra-pub-run-odp.s3.amazonaws.com`) is
reachable** from this sandbox — direct file downloads by accession, no
NCBI eutils needed. This report is what happened when I actually ran the
pipeline end-to-end on real, from-scratch environmental data.

## What ran

1. **Pull** raw SRA file from AWS Open Data mirror:
   `curl https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR13291825/SRR13291825`
   → 18 MB, real Illumina MiSeq paired-end run (57k spots, 2×275 nt reads,
   verified from the .sra metadata, multi-barcode multiplexed).
2. **Convert** `.sra` → paired FASTQ with `fasterq-dump` (from apt-installed
   `sra-toolkit`).
3. **Assemble** with `megahit` (apt-installed, `--min-contig-len 200`):
   26 contigs, N50 383 bp, max 506 bp. Small, but real.
4. **Sweep** with `python scripts/phase3_sweep.py --input <contigs.fa>`:
   the exact command Phase 2 was calibrated for.

## What we found (first pass, no preprocess)

The pipeline flagged 9/14 clusters significant at α=0.05 BH-FDR. Top hit:

| rank | contig | length | coverage | struct z | orphan | SAS | circular | sig |
|---|---|---|---|---|---|---|---|---|
| 1 | k141_0 | 379 | 137× | 21.03 | 0.57 | ✅ | no | ✅ |
| 2 | k141_11 | 454 | 102× | 3.26 | 0.52 | ✅ | no | ✅ |
| 3 | k141_10 | 422 | 5× | 3.69 | 0.55 | no | no | ✅ |

Full output: `results/real_sweep_SRR13291825/sweep_results.json`.

**But** inspecting the top hit `k141_0`'s sequence:

```
TTTTTTTTTTTTCAAGCAGAAGACGGCATACGAGATTGAGCTAGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT...
             └── Illumina P7 primer ──┘             └── TruSeq Read 2 adapter ─┘
```

It's an **assembled Illumina adapter dimer**, not real biology. Adapter
sequences fold stably (large z), are highly abundant (137× coverage),
and appear on both strands (R1 and R2 adapters are reverse-complementary),
so they perfectly mimic the SOS signature. Phase 2's synthetic decoys
didn't include this failure mode, because I hadn't run on real data yet.

## Fix, calibrated on the same real data

Added `src/rnasig/preprocess.py`: strips any contig containing ≥20 nt of
any canonical Illumina TruSeq/Nextera adapter (either orientation), or
any ≥12 nt homopolymer run. Wired into `scripts/phase3_sweep.py` as
default-on; tests in `tests/test_preprocess.py`. This is a real
methodology change discovered *by running on real data*, not an
after-the-fact rationalization.

**Re-run of exactly the same sweep with filter on:**

- Preprocess: 26 → 20 contigs (**6 adapter contigs dropped**,
  including k141_0 the false top hit).
- Clusters: 10, 1 non-singleton.
- Significant at α=0.05: **4** (down from 9).
- New top hit: `k141_11` (length 454, z=5.04, orphan=0.52, SAS=✅,
  coverage 102×) — the same contig that ranked #2 unfiltered, now
  correctly promoted after the adapter noise was removed.

Full output: `results/real_sweep_SRR13291825_filtered/sweep_results.json`.

## Characterization of the new top hit

`k141_11` at 200 shuffles: length=454, GC=45%, structure z=**4.58**
(highly stable), orphan=0.52, longest ORF=258 nt (86 aa, non-trivial),
predicted MFE structure has multiple stacked helices (SVG in
`results/real_char_k141_11/`).

**What this contig probably is — honest reading:** the sequence's first
~60 nt are `CAGCTCTAATAGCGTATACTAAAGTTGTTGCAGTTAAAAAGCTCGT...`, which
matches the well-known conserved-region pattern of eukaryotic 18S small-
subunit ribosomal RNA (specifically around the V4 variable region). If
so, this hit is:
- **Real** ✓ (it's a real, structured, dual-strand RNA in the sample)
- **Not novel** ✗ (ribosomal RNAs are the most abundant structured RNAs in
  any transcriptome and are what the SOS signature will always fire on
  first, in the absence of an rRNA-exclusion step)

I cannot BLAST/SILVA-confirm this from inside the sandbox (no route to
BLAST or SILVA), so I'm calling it a *likely* rRNA fragment based on the
sequence match at the terminus, not a confirmed identification.

## What this run actually established

This is not a novel-organism discovery, and I'm not going to pretend it is.
It's four concrete, honest things:

1. **The pipeline runs end-to-end on real public sequencing data**, pulled
   fresh from AWS, assembled from scratch, in this sandbox — not just on
   pre-packaged toy data.
2. **The signature fires on real data** in exactly the way the definition
   predicts: structured, high-coverage, sense/antisense-paired contigs.
3. **The first honest run revealed a real methodology gap** (adapter-dimer
   false positive) that Phase 2's synthetic calibration missed, and that
   gap is now fixed and tested (`src/rnasig/preprocess.py`,
   `tests/test_preprocess.py`).
4. **The top surviving hit is likely rRNA**, which is the correct thing to
   surface first in a signature-only pipeline. To find genuinely novel
   elements the next fix would be an rRNA/tRNA-exclusion step against
   SILVA/Rfam — impossible from inside this sandbox, but a one-line
   `easel-esl-sfetch` + `cmscan` add-on for anyone running this with
   real database access.

## What a follow-on run would need

Same pipeline, one additional filter step:

```
# after phase3_sweep, before phase4:
cmscan --rfam --tblout hits.tbl candidates.fa
# drop any candidate with an Rfam hit E<1e-5 in the tRNA/rRNA/snoRNA families
```

Then whatever survives goes to `phase4_characterize.py` for the deep dive,
and from there to actual comparative-genomics or wet-lab confirmation.
This last mile is what a genuine "found a new class of biology" claim
needs, and it does not happen inside this sandbox.
