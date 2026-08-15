# What this search found, and what it found out about itself

A consolidated record. The phase reports have the detail; this is the
summary and the honest accounting.

## Headline

No novel RNA element was found. Every candidate that survived filtering and
reached a database was identified as host sequence, bacterial sequence,
ribosomal RNA, or the deliberately planted positive control.

The substantive results are about the method. Four separate defects were
found in the search itself, three of them by validating against real
reference molecules rather than against synthetic nulls, and each one
would have made a "nothing found" result meaningless.

## What was searched

| sweep | target | runs | outcome |
|---|---|---|---|
| Phase 3 (inherited) | soil "metatranscriptome" | 1 | void: the run is a DNA amplicon library |
| Phase 6 breadth | freshwater, permafrost, peat, hot springs, plant | 12 | 8 circular, 0 structured |
| Phase 6 depth | permafrost | 2 | 2 circular, both unstructured |
| Phase 6b | human oral + gut | 23 | 27 circular, 11 structured, all identified |
| Obelisk control | *S. sanguinis* SK36 | 4 | 57 circular, obelisk not recovered |
| Plants | 28 species, viroid-scale | 28 | in progress at time of writing |

Every run passed a library-type gate before being swept.

## Everything that was identified

Across all sweeps, every candidate reaching blastn or blastx:

| candidate | identification | identity |
|---|---|---|
| SRR11060618 337 nt | **Peach latent mosaic viroid** (control) | 96% nt |
| k141_25 | FAD-dependent oxidoreductase, *Pseudolabrys* | 97% aa |
| k141_10 | DUF2007 protein, Acidobacteriota | 88% aa |
| k141_3 | MBL fold metallo-hydrolase, Tepidisphaeraceae | 73% aa |
| SRR19432446 264 nt | *Streptococcus oralis* hypothetical protein | 89% aa |
| SRR8204361 1481 nt | *Cutibacterium* sp. genome | 100% nt |
| SRR8204361 929 nt | *Escherichia coli* genome | 100% nt |
| SRR19432464 810 nt | *Streptococcus cristatus* | 98% nt |
| SRR11790910 549 nt | *Streptococcus sanguinis* genome | 100% nt |
| SRR19432494 425 nt | *Prevotella melaninogenica* | 100% nt |
| SRR27841025 462 nt | *Malus domestica*, metallothionein mRNA | 97% nt |
| SRR19432462 342 nt | 28S rRNA | 100% nt |
| SRR19432484 336 nt | *Abiotrophia defectiva* | 92% nt |
| SRR19432522 294 nt | *Capnocytophaga leadbetteri* | 96% nt |
| SRR19432483 258 nt | *Streptococcus* sp. | 93% nt |
| SRR19432522 252 nt | *Corynebacterium matruchotii* | 83% nt |
| SRR19432553 228 nt | *Neisseria perflava* | 97% nt |
| SRR5105111 165 nt | *Syzygium* mRNA (both Myrtaceae) | 88% nt |
| SRR19432483 159 nt | *Neisseria mucosa* | 98% nt |
| DRR047292 159 nt | *Eragrostis tef* chromosome (both grasses) | 87% nt |

Nothing unexplained. The single molecule that is what the pipeline was
built to find is the one that was put there on purpose.

## The four defects

**1. The input was the wrong kind of data.** Phase 3's "real sweep" ran on
SRR13291825, which is an 18S rRNA amplicon survey of soil DNA. No RNA, so
every RNA-level statement downstream was void at the source. Established
three independent ways: archive metadata, three 20-mers covering 68.5% of
reads, and full reads hitting 18S at 100%. This is what `srameta` now
prevents, and it caught a repeat of the same mistake on MGnify: a study
returned by that archive's own `experiment_type=metatranscriptomic` filter
whose underlying runs are WGA/METAGENOMIC. 602 circular contigs from it
were about to be scored as RNA.

**2. The orphan score is 0 for 4 on real coding sequence.** k141_3,
k141_10, k141_25 and SRR19432446 scored 0.58, 0.55, 0.52 and 0.65, all
above the 0.5 that reads as non-coding, and all four encode proteins at
73-97% amino-acid identity. The failure is structural: S2 wants a clean
ATG-to-stop ORF and an internal gene fragment has none. Two of those four
were also non-significant by blastn and only resolved by blastx, because
synonymous substitution hides conservation at the nucleotide level.

**3. The structure z-score rejects real viroids.** z >= 3 was calibrated on
PLMVd alone, which scores 9.4. Against a panel fetched from NCBI it
rejected Pear blister canker viroid (z = -0.02), Apple scar skin viroid
(3.06) and Avocado sunblotch viroid, three of seven. Viroid base
composition is itself self-complementary, so a dinucleotide shuffle folds
about as well as the molecule and z collapses. Replaced by sequence
complexity, which does the job z was actually needed for: nine described
viroids score 0.90-0.98, repeat decoys score 0.01. New filter 9/9 viroids
and 0/5 decoys, against 4/7 before.

**4. Circularity was calibrated against the wrong null.** Allowing one
mismatch in the terminal k-mer, added to recover two real PLMVd contigs,
measured a false-positive rate of 0.00 against shuffled sequence and 1.18%
against 5017 real assembled contigs. Shuffling destroys the repeat
structure that produces spurious terminal matches. The excess showed up as
a spike of "circular" 228 nt units across *Sideroxylon*, *Brassica*,
*Melaleuca*, *Citrus*, human oral samples and *Streptococcus sanguinis*;
the sequences behind it share zero 20-mers with each other. Same length,
unrelated sequence. Raising k to 16 keeps both real PLMVd contigs and
halves the rate.

Defects 3 and 4 were found within an hour of each other and share a cause:
both parameters had been calibrated against synthetic nulls instead of
against real molecules of the class being hunted.

## What the pipeline can and cannot do

Can: recover Peach latent mosaic viroid from an accession alone, with no
human choosing the dataset, the assembly or the candidate. 685 contigs, 7
circular, exactly one shortlisted, confirmed at 96% identity. The six
rejected circular contigs are the other half of that result.

Cannot: recover Obelisk-S.s from *S. sanguinis* SK36, where it is
documented at 1-13.3% of mapped reads. Three hypotheses were tested and
rejected (error clouds, strandedness, and a real truncation bug in greedy
extension that was found and fixed). The molecule still does not come out.

So obelisk-scale absence is not a claim this work supports. Viroid-scale
absence is supported only for the runs actually swept, at the depths
recorded, and only since the recalibration.

## The honest summary

The instrument had more wrong with it than the data did. That is worth
saying plainly, because the temptation in a search like this is to report
the candidate count and let it imply something. Every count this pipeline
produced before the recalibration was dominated by artifacts of its own
filters, and the one time a striking pattern appeared, a 228 nt element
recurring across seven plant genera, it was a detector artifact that took
sequence-level comparison to rule out.

The nulls are now bounded by a stated detection floor, a working positive
control on the viroid class, a documented failure on the obelisk class, and
filters calibrated against nine real viroids rather than one.
