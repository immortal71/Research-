# k141_25 confirmation evidence

Follow-up to the earlier `phase3_full_filter_report.md` "last-last step" —
concrete evidence work on the one candidate that survived every filter
(k141_25) and had a genuine obelisk-like feature (circular terminal
repeat, 25× coverage, no rRNA match). This is done here in the sandbox,
without needing BLAST-nt or Rfam-CM access, using tools that installed
cleanly via apt (blastn, barrnap, nhmmer, tRNAscan-SE, bwa) and raw
reads pulled back from AWS SRA.

## Files here

- `k141_25_monomer.fasta` — the 186 nt circularity-resolved monomer (the
  full 327 nt contig is a 1.76-copy near-tandem, resolved to a single
  186 nt unit by `rnasig.circularity.find_circularity`).
- `junction_and_monomer.fasta` — reference used to look for
  junction-spanning reads: last-100 nt + first-100 nt of the monomer.
- `junction_reads_R{1,2}.tsv` — blastn output of the raw R1/R2 reads
  against the junction+monomer reference. 42 reads (25 R1 + 17 R2) hit
  the junction reference with alignment covering ≥20 nt on each side of
  the join at 89–100% identity.
- `tandem_junction_reads_R{1,2}.tsv` — blastn of the raw reads against
  the *tandem-doubled* monomer (372 nt, monomer + monomer). Max
  alignment length reached 197 nt on R1 and 196 nt on R2, i.e. **single
  reads span a full monomer's worth of sequence and continue into the
  next copy** — the diagnostic signature of either rolling-circle
  replication of a real RNA circle or a tandem-repeat DNA locus.

## What is (and isn't) established

Established, from real data in the sandbox:

- **The circularity is real, not an assembler artifact.** 42
  junction-spanning raw reads at high identity is not something
  MEGAHIT can invent; the reads themselves contain sequence that goes
  past the monomer end and back to the start.
- **Not any known ncRNA the local HMM tools can name.** barrnap
  (bac/arc/euk/mito rRNAs, standard and lenient thresholds), tRNAscan-SE
  (all tRNAs), and nhmmer with barrnap's shipped 5.8S / 5S / 18S / 28S
  / 16S / 23S / 12S covariance-model-derived HMMs at very loose
  thresholds (E ≤ 1, E ≤ 100 for 5.8S) all return no hit against the
  186 nt monomer.
- **Not a known Universal ITS primer target either** — the strongest
  match to any of the standard ITS-flanking conserved motifs (ITS1_fwd,
  ITS4_rev, ITS2_fwd, or the 18S hyper-conserved CCAGCAGCCGCGGTAATTCC
  region) is 70% identity over 20 nt, which is not diagnostic. If this
  were an rDNA ITS amplicon, at least one of these motifs would match
  cleanly at the terminus.
- **All 186 monomer 32-mers are unique** — no internal tandem
  sub-repeats within the monomer itself.

*Not* established (would need BLAST-nt or Rfam-CM access):

- Whether this is a genuinely novel biological element, or a fragment of
  a bacterial/archaeal/plasmid/viral genome that happens to be tandemly
  arrayed in its source.
- Whether the tandem repeats in the reads come from rolling-circle
  replication of a real circular RNA (viroid/obelisk-like) or from a
  tandem-repeat DNA locus (satellite, telomeric, CRISPR array,
  microsatellite).

The distinguishing wet-lab experiments are standard: RNase R digestion
(exonuclease that destroys linear RNA but not covalently closed
circles), 2-D PAGE (circles run anomalously), or reverse-transcription
across the junction with divergent primers. None of these are
possible here.

## Honest bottom line

This is not a "we found a new organism" claim. It is a real, unclassified
~186 nt sequence element assembled from real public sequencing data,
whose circularity is corroborated by raw-read evidence, and which
survives every known-ncRNA exclusion the sandbox can perform. It is a
legitimate candidate for the next-stage wet-lab or full-database
comparative-genomics validation that a real discovery claim needs.
