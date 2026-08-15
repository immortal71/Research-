# Phase 6b: hunting where the elements actually are

The first Phase 6 sweeps went to permafrost, hot springs, peat and
freshwater because `docs/LIMITATIONS.md` recommended under-sampled
environments. Under-sampled is not the same as likely. Obelisks were
reported from host-associated microbiomes, and Obelisk-S.s specifically
from *Streptococcus sanguinis*, an oral bacterium. This sweep went there
instead.

It found nothing novel, and this time that statement is worth something,
because every candidate it produced was run to ground rather than left as
a count.

## What was searched

26 runs, all passing the library-type gate as shotgun metatranscriptomes or
transcriptomes:

| set | source | runs | subset | detection floor |
|---|---|---|---|---|
| batch 1 | human gut + oral | 13 | 2M reads | ~10x |
| batch 2 | human oral | 10 | 2M reads | ~8x |
| control | *S. sanguinis* | 3 | 2M reads | whole run |

These runs are 3-4M reads each, so a 2M-read subset is half to two thirds
of the data. That is the difference from the first sweep, where a 300k
subset of a 100M run reached only elements above ~2000x.

## What came out

27 circular contigs in the viroid/obelisk size band. 11 cleared both
structure z >= 3 and rod-likeness. Every one of the 11 was then identified:

| candidate | length | cov | z | identification | identity |
|---|---|---|---|---|---|
| SRR11060618 | 337 | 161 | 9.04 | **Peach latent mosaic viroid** (control) | 96% nt |
| SRR8204361 | 1481 | 37.5 | 3.76 | *Cutibacterium* sp. genome | 100% nt |
| SRR8204361 | 929 | 22.8 | 4.90 | *Escherichia coli* genome | 100% nt |
| SRR11790910 | 549 | 63.5 | 3.73 | *Streptococcus sanguinis* genome | 100% nt |
| SRR19432462 | 342 | 6752 | 3.43 | 28S rRNA | 100% nt |
| SRR19432522 | 294 | 10.6 | 3.65 | *Capnocytophaga leadbetteri* | 96% nt |
| SRR19432483 | 258 | 36.9 | 3.22 | *Streptococcus* sp. | 93% nt |
| SRR19432522 | 252 | 11.6 | 4.55 | *Corynebacterium matruchotii* | 83% nt |
| SRR19432553 | 228 | 14.1 | 3.51 | *Neisseria perflava* | 97% nt |
| SRR19432483 | 159 | 8.5 | 6.18 | *Neisseria mucosa* | 98% nt |
| SRR19432446 | 264 | 70.6 | 3.19 | *Streptococcus oralis* protein | **89% aa** |

Ten bacterial or ribosomal sequences, all of them ordinary members of the
oral and gut communities being sampled, and one viroid that was put there
deliberately as a control. No novel element.

The positive control behaving correctly is what makes the rest
interpretable. PLMVd scored the highest z of anything in the sweep and was
confirmed at 96% identity, from an accession, with nobody choosing it.

## The last candidate, and why nucleotide search was not enough

SRR19432446_264nt survived every filter and every nucleotide search. Its
best blastn hit was E=0.001 at 75% identity, and a sensitive re-run at word
size 7 with E=100 did no better: marginal hits to *Staphylococcus* and
*Streptococcus*, nothing convincing. On nucleotide evidence alone it looked
unidentified.

blastx settles it in one line: a hypothetical protein from *Streptococcus
oralis*, **89% amino-acid identity**, 95% positives, E=1e-20, one frame,
covering the contig end to end.

This is the k141_25 lesson exactly repeated. Synonymous substitution hides
conservation at the nucleotide level and not at the protein level, so a
sequence can be 75% and non-significant by blastn while being 89% and
unambiguous by blastx. Any claim of novelty that rests on a nucleotide
search alone is not safe.

The protein also explains everything else about the contig. It is a
low-complexity acidic repeat, `...DLDVLSLALTDSEADLLAEALSEADLDKLSDSELDLLVDVDLE...`,
so the underlying DNA is repetitive, which makes it partly self-
complementary, which is why it folded into a 66%-paired rod and cleared the
rod-likeness bar. A repeat protein looks like a structured RNA to a folding
program.

## The orphan score is now 0 for 4

Every real protein-coding sequence this pipeline has met, it has called
non-coding:

| sequence | orphan score | truth |
|---|---|---|
| k141_3 | 0.58 | MBL fold metallo-hydrolase, 73% aa |
| k141_10 | 0.55 | DUF2007 protein, 88% aa |
| k141_25 | 0.52 | FAD-dependent oxidoreductase, 97% aa |
| SRR19432446 | **0.65** | *S. oralis* hypothetical protein, 89% aa |

Above 0.5 reads as non-coding. The newest miss is the worst score of the
four. S2 looks for a clean ATG-to-stop ORF and in-frame codon bias, and a
contig that is an internal fragment of a gene has neither, so a genuine
coding fragment scores as orphan. It should be treated as a weak prior at
best, never as evidence, and always paired with a blastx check.

## The obelisk control did not work, and that is a real limitation

Three *S. sanguinis* transcriptomes were swept specifically to see whether
this pipeline can recover an obelisk, since Obelisk-S.s was reported from
that species and is the actual target class rather than a related one.

It did not find one. The only candidate that cleared the filters,
SRR11790910_549nt, is 100% identical to the *S. sanguinis* chromosome.

Two explanations, and this sweep cannot separate them. The runs chosen may
simply not be strain SK36 or may not carry the element, in which case the
absence says nothing about the method. Or the pipeline may not detect
obelisks, in which case every negative in this report is weaker than it
looks. The PLMVd control shows the pipeline finds a viroid, but obelisks
are around 1 kb, protein-coding, and structured differently, so viroid
recovery does not transfer automatically.

Until an obelisk is recovered from a dataset known to contain one, the
honest statement is that this pipeline is validated on viroid-like
molecules and untested on obelisk-like ones.

## Bottom line

26 runs in the environments where these elements are actually reported, at
a depth that reaches anything above roughly 10x, produced 27 circular
contigs, 11 structured rods, and zero novel elements. Everything found was
a bacterium, a ribosome, or the control.

That is a real negative, bounded by a stated detection floor and a working
positive control, on the highest-prior environments available. It is also a
small sample, and the target-class control is still missing.
