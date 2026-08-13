# Phase 5: what the candidates actually were

Every earlier report in this repo ends at the same wall. `phase3_full_filter_report.md`
lists it as the next step to take, "BLAST/CMSCAN k141_25 against nt / Rfam / RefSeq
to rule out known viroids, small RNA viruses, or genomic origins", and
`real_k141_25_confirmation/README.md` files the same question under *not
established*. Both note that the databases were unreachable at push time.

They are reachable now. Running every candidate this repo ever produced
against nt, nr and Rfam answers two questions at once, and the answers point
in opposite directions.

## Short version

**Phase 1 found a real viroid, and the clustering was exactly right.** 24 of
the 38 test contigs are Peach latent mosaic viroid, a genuine circular,
non-coding, self-replicating RNA of 337 nt. The pipeline's top
sense/antisense cluster holds 20 of them and no host sequence at all, and
none of its other three clusters holds any viroid. The circularity code
independently resolved the repeat unit as 337 nt from sequence alone, with
no database involved, which is PLMVd's exact genome length. This is the
pipeline recovering exactly the class of molecule it was designed to find,
ranked first, on real data.

**Phase 3 was pointed at the wrong kind of data.** SRR13291825 is an 18S
rRNA amplicon survey of soil DNA. The molecules sequenced were DNA, and they
were chosen by a primer rather than sampled at random, so neither the
structure scores nor the coverage numbers describe anything that existed.
All three surviving candidates are fragments of bacterial protein-coding
genes. The arithmetic was fine. It was arithmetic about molecules that were
never in the tube.

The same preflight passes the first run and fails the second, which is the
distinction the pipeline previously had no way to draw.

## Phase 1: the true positive

SRR11060618 is real RNA-Seq (`library_strategy=RNA-Seq`,
`library_selection=Oligo-dT`, `library_source=TRANSCRIPTOMIC`). All 38
contigs in the test set were searched, not just the four cluster centroids,
which turns this into a benchmark with real ground truth rather than a spot
check.

**24 of 38 contigs are Peach latent mosaic viroid. The other 14 are
*Prunus* host sequence. Nothing is unidentified.**

That makes it possible to score the clustering exactly:

| SAS cluster | centroid | n | PLMVd | host |
|---|---|---|---|---|
| 1 | NODE_36652 | 20 | **20** | 0 |
| 2 | NODE_38543 | 5 | 0 | 5 |
| 3 | NODE_43211 | 4 | 0 | 4 |
| 4 | NODE_52236 | 2 | 0 | 2 |

The circular-permutation clustering put 20 viroid contigs in one cluster and
no host sequence in it, and put no viroid in any of the other three.
Precision on the viroid cluster is 20/20. Recall across all PLMVd contigs is
20/24, the four misses being contigs that did not land in a
sense/antisense-paired cluster at all. On a real dataset with a real answer,
the method separates the replicating agent from its host cleanly.

Coverage separates the two classes just as cleanly, with no overlap: host
contigs run 2.2x to 595x, viroid contigs 1,298x to 24,158x. The gap is a
factor of 2.2 with nothing inside it. A viroid replicating in host tissue
should be far more abundant than any single host transcript, and it is.

NODE_36652 specifically scored structure z=17.13, orphan 0.89, circular, at
17,349x coverage. Those are the numbers a viroid should produce: viroids are
among the most thermodynamically structured RNAs known, they encode no
protein, and they replicate as covalently closed circles.

The circular unit length is the part worth dwelling on. `phase4_report.md`
recorded a 337 nt unit for NODE_36652, derived by `circularity.find_circularity`
from terminal-repeat structure in a 943 nt concatemeric contig. PLMVd's
complete genome, MZ220895.1, is 337 nt. The reimplementation reconstructed a
viroid genome to the base, months before anyone checked what it was.

## Phase 1's data is not what the repo said it was

While confirming the above: SRR11060618 is *Prunus persica* stamen
ssRNA-seq (study PRJNA493230), not the human gut metatranscriptome that
`PROVENANCE.md`, `phase1_report.md`, `phase4_report.md`, `LIMITATIONS.md`
and the README all described. The Obelisks paper is about human microbiomes,
and that description was carried across to the test file without being
checked against the archive.

It makes the Phase 1 result better rather than worse. A peach transcriptome
is the natural place to find a peach viroid, the contigs are peach sequence
throughout, and the run is still real RNA-Seq. But it is the second time in
this repo that an unchecked assumption about a dataset propagated into five
documents, which is the same failure that produced the Phase 3 problem
below. Those files are now corrected.

## Phase 3: the false positives

| contig | encodes | closest taxon | aa identity | E |
|---|---|---|---|---|
| k141_25 | FAD-dependent oxidoreductase | *Pseudolabrys* / *Bradyrhizobium* (Alphaproteobacteria) | 97% | 8e-13 |
| k141_10 | DUF2007 signal-transducing protein | Acidobacteriota / Bryobacteraceae | 88% | 2e-56 |
| k141_3 | MBL fold metallo-hydrolase | Tepidisphaeraceae (Planctomycetota) | 73% | 5e-38 |

Rfam cmscan returns zero hits for all of them. None is a structured RNA of
any known family, because none of them is RNA.

### How the library type was established

Three independent lines, none depending on the other two.

**1. Archive metadata.** ENA and the SRA runinfo service both report
SRR13291825 as `library_strategy=AMPLICON`, `library_selection=PCR`,
`library_source=METAGENOMIC`, `scientific_name=soil metagenome`, library
name `eGT55/10-20 cm_a` (a soil depth horizon), 56,865 read pairs on a
MiSeq. Study PRJNA685954, GFZ Potsdam.

**2. The reads are primer-led.** Three 20-mers account for 68.5% of R1 and
68.7% of R2:

```
AGCTCCAATAGCGTATATTA   16324   28.7%
AGAAGACATCCTTGGTGAAT   11542   20.3%
TGAAAACATCCTTGGCAAAT   11104   19.5%
```

Each has lower-frequency variants differing at one position
(`...TATTT`, `AGCTCTAAT...`, `TGAAAACGTCC...`), which is what degenerate
primer positions look like in read data. A shotgun library's most common 5'
20-mer sits far below 1%, because fragments break wherever they break.

**3. The reads are 18S.** Full-length representative reads for each prefix
hit small-subunit rRNA:

| prefix | best nt hit | identity | E |
|---|---|---|---|
| primer C | uncultured *Trebouxia* / Trebouxiophyceae 18S | 100% | 7e-138 |
| primer B | uncultured Cercozoa / *Polymyxa* 18S | 100% | 1e-139 |
| primer A | Nanorchestidae sp. 18S (soil mites) | 80-90% | 4e-40 |

Soil algae, soil protists and soil microarthropods, which is the expected
content of a eukaryotic V4 metabarcoding survey. The primer A read begins
`AGCTCCAATAGCGTATATTAAAGTTGTTGCGGGTAAAAAGCTCGTAGTTGGATCTCT`, the canonical
eukaryotic 18S V4 opening.

### Why the pipeline produced a candidate anyway

The filters worked. The 12 contigs barrnap dropped as known rRNA were the
*on-target* amplicons, which is what an 18S survey is supposed to produce.
What survived to be ranked was the off-target residue: non-specific PCR
products and genomic carryover, which is bacterial DNA and therefore neither
rRNA nor a novel RNA.

k141_3 shows the mechanism in one contig. Its 3' end reads

```
[ MBL hydrolase gene fragment ] GGAATTACCGCGGCTGCTGG AGATCGGAAG
                                 ^ universal SSU primer  ^ Illumina adapter
```

Insert, then primer, then adapter, which is how an amplicon read is built
and not how a transcript is. That exact junction appears in 894 raw reads,
so it is real chemistry rather than an assembly glitch. The earlier report
read the `CCAGCAGCCGCGGTAATTCC` motif in this contig as evidence of a stray
18S fragment. The motif is genuinely there and is genuinely SSU, but it is
there as a primer site fused to bacterial DNA. The conclusion "probably
still rRNA" was right that the contig was not novel and wrong about what it
was.

k141_25's circularity has the same explanation. The 42 junction-spanning
reads and the 1.76-copy tandem are not fabricated. PCR concatemers are a
routine amplicon artifact, and a concatemer is indistinguishable from a
circle by read topology alone. `real_k141_25_confirmation/README.md` already
listed "a tandem-repeat DNA locus" as the alternative it could not exclude.
The library type excludes the RNA branch outright.

## The orphan score failed on the only real test it has had

`docs/LIMITATIONS.md` flagged S2 as "a lightweight, unsupervised proxy...
validated only against our own synthetic decoys". This is the first time it
met real sequence, and it scored all three protein-coding contigs as
borderline non-coding:

| contig | orphan score (>0.5 reads as non-coding) | truth |
|---|---|---|
| k141_3 | 0.58 | protein-coding, 73% aa identity |
| k141_10 | 0.55 | protein-coding, 88% aa identity |
| k141_25 | 0.52 | protein-coding, 97% aa identity |

Three out of three wrong, and the worst miss is on the strongest protein.
k141_25's nucleotide identity to any single genome is only 78%, so a
nucleotide search alone leaves it ambiguous. The protein search settles it,
because synonymous substitution hides conservation at the nucleotide level
and not at the amino-acid level. Any future use of S2 should treat a blastx
check as mandatory rather than optional.

Worth noting against this: S2 gave NODE_36652 a 0.89, its highest score on
any real contig, and that one really is non-coding. The proxy is not
useless. It is unreliable in the direction that matters most, which is
calling something non-coding when it is not.

## What changed in the code

`srameta.py` preflights a run against the ENA portal and refuses AMPLICON
strategies, PCR-style selections and DNA sources before a sweep starts. It
also refuses runs whose metadata is missing, because treating an absent
annotation as permission is how this run got swept.

`ampliconqc.py` reaches the same verdict without metadata, by measuring how
concentrated the 5' prefixes of the reads are. A submitter's annotation can
be wrong or absent, and this test depends on nothing except the reads
themselves.

`dbcheck.py` resolves sequences against nt, nr and Rfam through the NCBI URL
API and Rfam's batch service, which turns "is this known?" from a promise in
a limitations file into a stage you can re-run.

`preprocess.py` gains a terminal adapter check. Read-through truncated at a
contig end leaves fewer than 20 nt of adapter, which the existing interior
rule could not see, and that gap is what carried k141_3 through two rounds
of filtering.

`scripts/phase5_identify.py` runs both halves.

The preflight separates the two runs in this repo correctly, in a few
seconds each:

```
$ python scripts/phase5_identify.py --run SRR11060618
  SRR11060618: strategy=RNA-Seq selection=Oligo-dT source=TRANSCRIPTOMIC (Prunus persica)
  VERDICT: compatible with the SOS signature

$ python scripts/phase5_identify.py --run SRR13291825
  SRR13291825: strategy=AMPLICON selection=PCR source=METAGENOMIC (soil metagenome)
  VERDICT: INCOMPATIBLE -- do not sweep this run
    - library_strategy=AMPLICON: a targeted marker-gene survey...
    - library_selection=PCR: template was amplified from specific primers...
    - library_source=METAGENOMIC: the molecules sequenced are DNA...
```

## What stands and what does not

Withdrawn:

- k141_25 as a candidate novel circular RNA element. It is a bacterial
  FAD-dependent oxidoreductase gene fragment.
- k141_10 and k141_3 as candidates of any kind.
- The Phase 3 real sweep and Phase 4's characterization of its output as
  evidence about RNA. They stand as a demonstration that the code runs end
  to end on real archive data.
- The description of SRR11060618 as a human gut metatranscriptome.
- The claim in `real_k141_25_confirmation/README.md` that k141_25 is "a
  real, unclassified ~186 nt sequence element". It is classified.

Strengthened:

- Phase 1, considerably. It is no longer "the reimplementation agrees with
  VNom's expected pattern" but a scored benchmark against known answers:
  20/20 precision on the viroid cluster, 20/24 recall, no host sequence
  misplaced into it, and the viroid's genome length reconstructed exactly.
- The circularity module, which produced that 337 nt unit unaided.
- The adapter and rRNA filters, which did their jobs.

## The honest lesson

The pipeline never had a way to ask whether its input could contain the
thing it searches for. Three filters were built to catch contigs that mimic
the signature, and each of them caught something real, but all three run
after assembly, on sequences that should not have been assembled for this
purpose at all. Structure z-scores, coverage, dual-strand co-occurrence and
circularity were all computed correctly, and on SRR13291825 all of them
described a PCR product.

The reason that run got swept was that it was the one reachable from a
network-restricted sandbox. Reachability is not a sampling criterion, and
nothing in the code was in a position to say so. Now something is, and it
takes a few seconds to run.

The Phase 1 result is the other half of the same lesson. Given real RNA from
a host that carries a circular RNA replicon, this signature finds it, ranks
it first, and gets its genome length right without a reference. The method
was never the problem.
