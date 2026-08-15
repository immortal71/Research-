# The obelisk control failed, and that bounds everything else

`phase6_oral_gut_report.md` ended by flagging that the pipeline had never
been tested on the actual target class. It was validated on Peach latent
mosaic viroid, a 337 nt viroid, while obelisks are around 1 kb,
protein-coding, and structured differently. This is that test, and it did
not pass.

## The setup

Zheludev et al. reported Obelisk-S.s, 1137 nt, as highly abundant in
*Streptococcus sanguinis* SK36. The follow-up paper (PMC12198308) names the
exact runs it analysed, which removes the main ambiguity of the earlier
attempt: the previous *S. sanguinis* control used runs that were not
labelled SK36 at all, so its failure said nothing.

Four SK36 RNA-Seq runs were swept, all passing the library gate, 2M reads
each, which is 24% of the smallest run.

## The result

`SRR20627698`, 2M reads: 1489 contigs, 4 circular, **none in the obelisk
size range** and none clearing z >= 3. The circular contigs were 327, 297,
252 and 486 nt.

Since a null result is uninformative without knowing whether the molecule
could have been assembled at all, the assembly was inspected directly for
anything near 1137 nt. Two contigs fell in the 900-1600 nt window:

| length | coverage | circular | identification |
|---|---|---|---|
| 1067 nt | 157x | no | *S. sanguinis* SK36 chromosome, 100% nt |
| 909 nt | 204x | no | *S. sanguinis* SK36 chromosome, 100% nt |

Both are host transcripts. Obelisk-S.s is not in the assembly.

## Why: the assembler, not the filters

Three parameter settings were tried on the same reads:

| k | min_count | contigs | N50 | max | circular in range |
|---|---|---|---|---|---|
| 21 | 3 | 2339 | 244 | 1318 | 0 |
| 25 | 5 | 1489 | 222 | 1067 | 0 |
| 25 | 10 | 493 | 210 | 1039 | 0 |
| 31 | 5 | 594 | 190 | 708 | 0 |

N50 sits near 200 nt regardless. For a bacterial monoculture transcriptome
at this depth that is poor; a production assembler would reach kilobase
contigs. The filters were never the constraint, because the molecule never
reached them.

One plausible cause was ruled out rather than assumed. `assemble` does not
canonicalise k-mers, deliberately, so that sense and antisense stay
separable for the SAS axis. On an unstranded library that would split every
molecule's coverage across two graphs and fragment the assembly. It is not
the explanation here: only 1% of k-mers in this run have their reverse
complement also present, so the library is stranded and strand-specific
assembly is the right choice. Collapsing to canonical k-mers produced zero
contigs, which is worse, not better.

The fragmentation is intrinsic to greedy extension over a transcriptome
graph, where overlapping transcripts and residual error k-mers branch the
path constantly.

## What this does to the earlier results

The PLMVd recovery still stands, but it is now clear what it demonstrated
and what it did not. PLMVd is 337 nt and sat at 17,349x in its run. It is
short enough for a greedy walk to close and abundant enough that its k-mers
dominate every branch point. That is close to the easiest possible case.

So the pipeline is validated for **short, extremely abundant circular
RNAs** and demonstrably fails on a **1 kb element at moderate abundance**,
which is the size class obelisks occupy.

Every negative in `phase6_hunt_report.md` and
`phase6_oral_gut_report.md` has to be read against that. The detection
floors quoted there, ~8-30x, are floors on *read coverage*, and they assumed
the assembly step was not itself a bottleneck. For anything much longer than
a few hundred nucleotides it is. "No obelisk-like element in 26 oral and gut
metatranscriptomes" is not a claim this work can support; "no short,
highly abundant circular RNA" is closer to what was actually tested.

The candidate identifications are unaffected. Every contig that was found
and identified really is what blastn and blastx say it is.

## What would fix it

The assembler is the bottleneck and it is the newest, least tested
component. The honest options, in order of cost:

1. Run the sweeps on Linux with rnaSPAdes or MEGAHIT and use `assemble` only
   where those cannot run. The pipeline already takes a FASTA at
   `phase3_sweep.py --input`, so nothing else has to change.
2. Improve extension: paired-end information, a proper bubble-popping pass
   over nested bubbles, and error correction before graph construction.
   These are what a real assembler does and this one does not.
3. Re-run this exact control after either. It is a clean pass/fail on a
   published, documented case, which makes it the right regression test for
   any future change.

Until one of those happens, obelisk-scale absence should not be reported as
a finding from this pipeline.
