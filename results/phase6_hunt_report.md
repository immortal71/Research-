# Phase 6: a validated hunt, and what it did not find

Phases 1 to 5 built the parts and established what the earlier candidates
really were. What the repo still could not do was start from an accession
and end at a candidate, because getting contigs meant MEGAHIT or rnaSPAdes
and neither runs on Windows. That gap is why Phase 3 swept whichever dataset
happened to be reachable.

This phase closes it and then uses it. The headline is that the pipeline now
demonstrably works end to end, and that the environments swept with it came
back empty.

## The pipeline recovers a real viroid from an accession alone

Given only `SRR11060618` and no other input, `phase6_hunt`:

1. preflighted the library (RNA-Seq / Oligo-dT / TRANSCRIPTOMIC, accepted),
2. streamed 300,000 reads from ENA without downloading the 6 GB run,
3. assembled them into 685 contigs with `rnasig.assemble`,
4. found 7 circular,
5. shortlisted exactly **one** at structure z >= 3.

That one is Peach latent mosaic viroid: unit length **337 nt**, which is
PLMVd's exact genome length, at 161x coverage and z=9.78. Confirmed
afterwards by blastn against core_nt at **96% identity over 282 nt,
E=5e-125** (`blast_plmvd.txt`).

The six rejected contigs matter as much as the accepted one. They scored
0.94 to 1.69, i.e. no better than their own shuffles, and every one of them
was circular. On this sample the structure filter is 1/1 on sensitivity and
6/6 on specificity.

This is the first time anything in this repo has gone from an accession to a
confirmed molecule without a human choosing the dataset, the assembly, or
the candidate.

## What it did not find

Two sweeps, both negative.

**Breadth**, 300k reads per run across eight environments:

| environment | runs analysed | circular contigs | cleared z >= 3 |
|---|---|---|---|
| freshwater | 3 | 2 | 0 |
| permafrost | 3 | 2 | 0 |
| peat | 2 | 2 | 0 |
| hot springs | 2 | 1 | 0 |
| plant | 2 | 1 | 0 |
| **total** | **12** | **8** | **0** |

Eleven further runs were skipped, nearly all because ENA metadata requests
failed while three sweeps were hammering the API at once. Those are lost
coverage, not evidence of anything.

**Depth**, 2M reads per run on the smallest runs, where a fixed subset is a
larger fraction of the data:

| run | environment | contigs | circular | z >= 3 |
|---|---|---|---|---|
| SRR23637081 | permafrost | 2223 | 0 | 0 |
| SRR23637102 | permafrost | 2517 | 2 | 0 |

Both circular hits scored below zero.

## What the negative is worth

This is the part that needed the positive control. A sweep that finds
nothing is uninterpretable unless you know what it could have found.

The breadth sweep is weak evidence. A 300k-read subset of a 100M-read run
only reaches elements present above roughly 2000x in the full run. PLMVd sat
at 17,349x and is unusual; most replicons would be invisible at that depth.
Read the breadth table as "nothing extremely abundant", not "nothing".

The depth sweep is worth more. At 2M reads from a 16M-read run the subset is
around 12%, which puts the detection floor near 30-60x. Combined with a
pipeline demonstrated to recover a viroid at 161x and to reject six circular
decoys, "no structured circular RNA in these permafrost runs above roughly
50x" is a real statement.

It is also a small one. Two runs from one biome is not a survey, and
permafrost and hot springs are low-biomass environments where a replicating
RNA element may simply not be present at all.

## Circularity on its own is nearly worthless

The first version of this script shortlisted any circular contig, which was
wrong and worth recording. A de Bruijn graph cycles wherever a
low-complexity repeat does, so an AT-rich tandem repeat is indistinguishable
from a circular molecule by terminal repeat alone.

The first hot-springs run swept here produced six such contigs, GC 0.27 to
0.36, structure z from -0.55 to 1.51. Every circular contig found in either
sweep behaved the same way. Structure is the discriminator, and the gap
between the decoys and PLMVd is not subtle: 1.51 against 9.78 in the same
size range, from the same code path.

## MGnify, and the same trap a second time

Assembling each run in Python is slow, so pre-assembled contigs looked like
a much better source. MGnify serves them: one analysis of study
MGYS00006862 yielded 75,068 contigs and 602 circular ones in the
viroid/obelisk size band, at no assembly cost.

The study was reached through MGnify's own `experiment_type=metatranscriptomic`
filter. Its underlying runs are `library_strategy=WGA`,
`library_selection=RANDOM`, `library_source=METAGENOMIC`: whole-genome
amplified DNA. `srameta.assess_library` rejects it on two independent
counts, and the scoring run was killed before it finished.

Had those 602 been scored and reported, this repo would have published a
second set of "circular RNA candidates" containing no RNA, from a dataset
labelled as transcriptomic by the archive serving it. The gate written for
the Phase 5 post-mortem is the only reason that did not happen, which is a
better argument for preflighting than anything in the Phase 5 report.

The lesson generalises past this repo: an archive's own experiment-type
label is not a substitute for the library metadata of the underlying runs.

## Honest bottom line

No novel RNA element was found. What exists now is a search whose
sensitivity can be stated rather than guessed: it recovers a known viroid
from raw reads, rejects circular decoys, and reports nothing in five
under-sampled environments at a depth where only an unusually abundant
element would have shown up.

The obvious next step is more runs at the depth of the second sweep rather
than more environments at the depth of the first, since detection floor, not
breadth, is what is limiting. Obelisks were reported in a minority of gut
metatranscriptomes even when searched for properly, so a dozen runs was
always a long shot, and it should be described as one.
