# The obelisk control passes, via the axis that was missing

`phase6_obelisk_control_report.md` records the pipeline failing to recover
Obelisk-S.s from *Streptococcus sanguinis* SK36, where it is documented at
1-13.3% of mapped reads. Three explanations were tested and rejected there,
and a real assembler bug was found and fixed without changing the outcome.

The cause was none of those. It was that `phase6_hunt` only implemented
half the signature.

## What was missing

The signature this repo hunts is stable structure, plus low coding
potential, plus **either** apparent circularity **or** sense and antisense
strands of the same molecule turning up in the same sample. `cluster.py`
has implemented the second since Phase 1. `phase6_hunt` tested only the
first.

Circularity is the harder of the two to observe. It requires the assembler
to read through the junction and close the loop, which is exactly what
failed here. Sense/antisense co-occurrence requires nothing of the
assembler beyond keeping the strands apart, which `assemble` already does
by not canonicalising k-mers. The evidence was sitting in the assembly the
whole time and nothing looked at it.

`sasfind.py` reads it in one pass over k-mers rather than `cluster.py`'s
pairwise alignment, which is what makes it usable on the 2,549 contigs a
sweep of this run produces.

## What it found

99 sense/antisense pairs where the circularity test found nothing usable.
Seven of the eight deepest are ordinary *S. sanguinis* transcription on
both strands, 100% identity to the SK36 chromosome. The eighth is not.

`CANDIDATE_sk36_1011nt.fasta`

| property | value |
|---|---|
| length | 1011 nt |
| GC | 0.42 |
| paired fraction | 74% |
| MFE | -368.1 (-0.364/nt) |
| longest ORF | 693 nt, 230 aa, covering 69% |
| antisense partner coverage | up to 29,165x |
| blastn vs core_nt | **no significant hit** (best E=0.83) |
| blastx vs nr | **no significant hit** |
| rRNA k-mer fraction | 0.00 |

A roughly 1 kb rod encoding a single protein with no detectable homology to
anything in `nr`, transcribed on both strands. That is the description of
an obelisk; Oblins were named a new class precisely because they match
nothing known.

## Reproducibility

Seven 30 nt probes tiling the contig, searched against two SK36 runs it was
not assembled from:

| run | probes detected | probe hits per 400k reads |
|---|---|---|
| SRR20627698 | 6/7 | 889 |
| SRR23591555 | 7/7 | 17,739 |

Roughly 4% of the second library, which sits inside the 1-13.3% the
follow-up paper reports for Obelisk-S.s. Present in three independent runs,
so it is not an artifact of one assembly.

## What this is and is not

It is almost certainly Obelisk-S.s or a close relative. Right strain, right
size class (1011 against a reported 1137 nt, consistent with an incomplete
assembly), right structure, right architecture, right abundance, both
strands, reproducible.

It is not a confirmed identification, and two things are missing:

- **No circularity.** The contig carries no terminal repeat, so this work
  cannot show the molecule is circular. That is consistent with the
  assembler never closing the loop, which is separately documented, but
  absence of evidence is not evidence.
- **No reference to compare against.** Obelisk sequences are not deposited
  in GenBank under any searchable name, so the comparison that would settle
  it cannot be made from here. It needs the supplementary table from
  Zheludev et al.

It is also not a novel element. Obelisk-S.s is described. The result is
that the pipeline can now recover the target class, which is what every
negative in this repo depended on and did not have.

## What it costs the earlier results

Every sweep before this one tested circularity alone: 12 environmental
runs, 23 oral and gut runs, 28 plant runs. All of them were running half
the signature, and the half that fails precisely when the assembler cannot
close a circle.

Those nulls now need re-running with SAS enabled before they mean anything.
`phase6_hunt` scores a contig if it is circular **or** has an antisense
partner, and reports which.
