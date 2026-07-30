# Phase 3 full-filter report: after the ncRNA exclusion step

Follow-up to `results/phase3_real_sweep_report.md`. The earlier report
noted that "the top surviving hit is likely rRNA, which is the correct
thing to surface first in a signature-only pipeline. To find genuinely
novel elements the next fix would be an rRNA/tRNA-exclusion step against
SILVA/Rfam." That was the last step. This report is that step, done.

## What was added

`src/rnasig/rrna_filter.py`: wraps `barrnap` (rRNAs, apt package, ships
its own bacterial/archaeal/eukaryotic/mitochondrial rRNA HMMs) and
optionally `tRNAscan-SE` (opt-in — its CM database load is very slow,
zero hits on this sample, so off by default). Called after the adapter
filter, before clustering. Fails soft if either tool is not installed.
Tests in `tests/test_rrna_filter.py`; wired into `phase3_sweep.py` as
default-on with `--skip-ncrna` to disable.

Regression fix: FASTA headers from MEGAHIT/SPAdes carry metadata after a
space (`k141_11 flag=1 multi=101 len=454`), but barrnap reports the
leading token only (`k141_11`). The first version of this filter
compared full header vs. token and silently ignored every hit; caught
during the first run because the "12 dropped" count didn't match the
"20 → 8 kept" line and top-hit list still included known rRNAs. Fixed
by comparing on the leading token. `test_matches_barrnap_short_ids_against_fasta_full_headers`
in `tests/test_rrna_filter.py` locks the fix.

## What survives on the real SRR13291825 sweep, all filters on

```
preprocess: 26 -> 20 contigs (6 adapter, 0 homopolymer)
ncRNA filter: 20 -> 8 contigs (12 known-ncRNA dropped)
input contigs: 26 (length-filtered: 8)
clusters: 8 (0 non-singleton)
significant at alpha=0.05: 2
```

The 12 barrnap-flagged rRNAs (including k141_11 and k141_7, the earlier
top hits) are all dropped. Full JSON:
`results/real_sweep_SRR13291825_full_filter/sweep_results.json`.

Ranked surviving candidates:

| rank | contig | length | coverage | struct z (n=200) | orphan | SAS | circular | BH-sig | notes |
|---|---|---|---|---|---|---|---|---|---|
| 1 | k141_3 | 311 | 2× | 2.88 | 0.58 | no | no | **yes** | weak 18S rRNA at loose nhmmer E≤1 — probably still rRNA |
| 2 | k141_10 | 422 | 5× | 3.01 | 0.55 | no | no | **yes** | no rRNA match at any threshold — genuine non-rRNA candidate |
| 3 | k141_25 | 327 | 25× | 1.62 | 0.52 | no | **yes** | no | no rRNA match at any threshold — **circular** + moderately structured, sub-BH-threshold but the circularity flag is the most obelisk-like feature we found on real data |

`nhmmer` at lenient threshold (E ≤ 1) against barrnap's shipped rRNA HMMs
found a real 18S conserved-motif match in k141_3 at position 293-274:

```
  18S_rRNA 567 ccAGcAgccgcGGuAAuucc 586
                ccAGcAgccgcGGuAAuucc
    k141_3 293 CCAGCAGCCGCGGUAAUUCC 274
```

`CCAGCAGCCGCGGTAATTCC` is one of the most conserved motifs in the 18S
rRNA (present across essentially all eukaryotes). barrnap's default
`--reject 0.05 --evalue 1e-3` missed it because the aligned fraction is
tiny; nhmmer at `E ≤ 1` catches it. So k141_3 is *likely* another 18S
rRNA fragment escaping the standard filter, not a novel element.

k141_10 and k141_25 show no such match at any threshold. They are the
genuine "signature-positive AND known-ncRNA-negative" candidates from
this sweep.

## Honest read of what we have

**k141_10**: 422 nt, structure z=3.01 (real but modest — one-sided p ≈
0.001), orphan=0.55 (borderline non-coding), coverage 5× (low —
essentially 3 read pairs, so any structure claim on it is on thin
statistical ice from the assembly side). No sense/antisense partner. Not
circular. This is a "maybe" — the pipeline flags it, but 5× coverage on
one contig from 26 total is not remotely enough to claim a novel
element without an independent replicate or wet-lab validation.

**k141_25**: 327 nt, structure z=1.62 (marginal), orphan=0.52, coverage
25× (adequate), **circular** (terminal-repeat detected — 3-copy tandem
resolution), **not** SAS-positive. Its structural stability is only
marginally above shuffle background (well below BH-significance), but
its circular signature is the strongest obelisk-like feature this whole
sweep produced, and it has enough coverage to be more than assembly
noise. Even so: one circular contig with a marginal z-score in one
smallish sample is a candidate for follow-up, not a discovery.

## What this actually establishes

Not: a novel organism.

Actually:

1. The full three-stage filter pipeline (adapter → rRNA/tRNA → SOS
   scoring → BH-FDR) runs end-to-end on real public sequencing data
   pulled fresh from a public archive, and drops the right things at
   the right stages.
2. Each filter stage was validated *by finding a real failure it fixes*
   — adapter dimers (v1), then rRNAs (v2). Both fixes are tested and
   committed.
3. The one surviving candidate that shows an obelisk-like signature
   (**k141_25**, circular, moderately structured, 25× coverage, no
   rRNA match) is a real hit from a real environmental sample. It is
   not "we found a new class of biology". It is exactly what a
   candidate-generation pipeline is supposed to produce: a
   short-list one item long, from a small sweep, for follow-up.

## What the honest next step would be

- Sweep more samples (this was 26 contigs from a 18 MB run; a real
  discovery sweep would cover many orders of magnitude more data).
- BLAST/CMSCAN k141_25 against nt / Rfam / RefSeq to rule out known
  viroids, small RNA viruses, or genomic origins. Rfam full CM database
  wasn't reachable from this sandbox at push time.
- If it survives that, pull the raw reads for k141_25's monomer, verify
  the circularity is real (paired-end reads should span the junction if
  it's a real circle), and confirm coverage on both strands.
- Only then does "possible novel structured RNA element" become a
  claim worth making, and even then only in a paper with the wet-lab
  or comparative-genomics confirmation this sandbox can't do.
