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

All four runs, 2M reads each: **57 circular contigs, none clearing z >= 3,
and no obelisk**. `SRR20627698` alone gave 1489 contigs and 4 circular ones
at 327, 297, 252 and 486 nt.

Obelisk-S.s cannot have been missed for lack of depth. The follow-up paper
measures it at **1-13.3% of total mapped reads**, present in all 17 datasets
it examined and more abundant than any mRNA in 11 of them. At 2M reads that
is tens of thousands of reads for a 1137 nt molecule. Whatever went wrong
was not sampling.

Since a null result is uninformative without knowing whether the molecule
could have been assembled at all, the assembly was inspected directly for
anything near 1137 nt. Two contigs fell in the 900-1600 nt window:

| length | coverage | circular | identification |
|---|---|---|---|
| 1067 nt | 157x | no | *S. sanguinis* SK36 chromosome, 100% nt |
| 909 nt | 204x | no | *S. sanguinis* SK36 chromosome, 100% nt |

Both are host transcripts. Obelisk-S.s is not in the assembly.

## It was not missed, it was shattered

Since the molecule is far too abundant to have been sampled away, the
assembly was searched for it directly. Three of the highest-coverage
contigs in the run have no significant nucleotide hit at all:

| length | coverage | best nt hit |
|---|---|---|
| 194 nt | 38,911x | none, E=5.6 |
| 566 nt | 21,390x | none, E=0.45 |
| 282 nt | 15,228x | none, E=0.72 |

Unidentified sequence at that abundance in a monoculture is what an
obelisk looks like, since obelisk sequences are not deposited in GenBank
under any searchable name. The molecule was in the assembly as fragments
and never as a circle, so the circularity test could not fire and nothing
downstream ever saw it.

The mechanism was then read straight off the graph. Every one of those
fragments still had successors available, meaning extension had stopped
with somewhere left to go rather than running out of sequence.
`greedy_contigs` filtered its extension options through the global
`visited` set, so a molecule whose path crossed nodes an earlier,
higher-coverage contig had already claimed was cut dead at the crossing.
Abundant molecules cross claimed nodes constantly, and a quasispecies does
it at every shared variant.

Two other explanations were tested first and rejected rather than assumed:

- **Strandedness.** `assemble` deliberately does not canonicalise k-mers,
  which would fragment an unstranded library. Only 1% of k-mers in this run
  have their reverse complement present, so the library is stranded and the
  choice is right. Collapsing to canonical k-mers gives zero contigs.
- **Error clouds.** At 10,000x-plus coverage, sequencing errors sit far above
  a `min_count` of 5. Raising it to 50, 200, 500 and 1000 made the assembly
  strictly worse at every step, so the extra k-mers were real, not noise.

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

## The bug was real. Fixing it did not recover the obelisk.

The truncation described above is a genuine defect and it is fixed.
`greedy_contigs` no longer filters extension through the global `visited`
set; a path may cross claimed nodes for a bounded run of 2k before it
stops, which keeps the walk linear. Unblocking crossings outright also
works and is unusable, because every seed then re-traverses the graph and
a 1.5M-read assembly does not finish. The dedup this required was
quadratic at first and is now a single k-mer-coverage pass.

On synthetic data the fix does exactly what it should. A 1137 nt circle,
the obelisk's own length, carrying variants in 30% of reads and sitting in
300 background transcripts, assembles in 3.9 s and returns unit=1137
exactly. Before the fix that case fragmented.

On the real SK36 run it changes nothing. The same three circular contigs
come back, 1033, 1212 and 1180 nt, all previously identified as
*S. sanguinis* chromosome. Obelisk-S.s is still absent.

The PLMVd control was re-run to confirm nothing was broken in the process:
733 contigs, 8 circular, exactly one shortlisted at unit=337, z=8.13, 68%
paired, with the other seven correctly rejected. Unchanged.

One observation argues against the reading offered above, and it belongs
here rather than in a footnote. The three unidentified high-coverage
fragments sit at 38,911x, 21,390x and 15,228x, a 2.5-fold spread.
Fragments of one circular molecule should have broadly similar coverage.
That is more consistent with several distinct abundant RNAs than with one
shattered obelisk. Without the reference sequence, which is not deposited
in GenBank under any searchable name, this cannot be settled either way,
and the favourable reading should not be assumed.

So: a real bug was found and fixed, and it was not the reason the control
fails. Three hypotheses have now been tested and rejected (error clouds,
strandedness, visited-set truncation), and the molecule still does not
come out.

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
