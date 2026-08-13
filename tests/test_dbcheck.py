"""Tests for the remote-homology module.

No network is touched: the report fixtures below are trimmed verbatim
from the real NCBI responses for k141_25 and k141_10, so the parser is
tested against the format the service actually returns rather than one
invented to match the parser.
"""
from rnasig.dbcheck import Hit, Identification, parse_blast_descriptions

# Real blastx-vs-nr response for the 186 nt k141_25 monomer.
K25_BLASTX_REPORT = """
Query= k141_25_monomer_186nt

Length=186


                                                                   Score     E     Max
Sequences producing significant alignments:                       (Bits)  Value  Ident

HZP76702.1 FAD-binding protein [Pseudolabrys sp.]                  73.6    8e-13  97%
MHH8797386.1 hypothetical protein [Xanthobacteraceae bacterium]    73.2    9e-13  94%
HEV8464776.1 FAD-binding protein [Pseudolabrys sp.]                73.2    1e-12  97%

ALIGNMENTS
>HZP76702.1 FAD-binding protein [Pseudolabrys sp.]
Length=569
"""

# Real blastn-vs-core_nt response for k141_10: nothing reaches significance.
K10_BLASTN_REPORT = """
Query= k141_10

Length=422


Sequences producing significant alignments:                       (Bits)  Value  Ident

CP101156.1 Natronosalvus amylolyticus strain WLHSJ1 chromosome...  49.1    1.1    80%
CP191906.1 Kutzneria sp. CA-250172 chromosome, complete genome     49.1    1.1    87%

ALIGNMENTS
"""

EMPTY_REPORT = """
Query= something

Length=100

***** No hits found *****
"""


def test_parses_accession_description_and_scores():
    hits = parse_blast_descriptions(K25_BLASTX_REPORT)
    assert len(hits) == 3
    top = hits[0]
    assert top.accession == "HZP76702.1"
    assert "FAD-binding protein" in top.description
    assert top.evalue == 8e-13
    assert top.score == 73.6
    assert top.identity == 97.0


def test_stops_at_alignments_marker():
    hits = parse_blast_descriptions(K25_BLASTX_REPORT)
    assert all(not h.accession.startswith(">") for h in hits)
    assert "Length=569" not in " ".join(h.description for h in hits)


def test_respects_max_hits():
    assert len(parse_blast_descriptions(K25_BLASTX_REPORT, max_hits=2)) == 2


def test_report_with_no_table_yields_nothing():
    assert parse_blast_descriptions(EMPTY_REPORT) == []


def test_verdict_prefers_protein_evidence_over_nucleotide():
    """The k141_25 case: 78% nt identity is ambiguous, 97% aa identity is not."""
    ident = Identification(
        query_id="k141_25",
        length=186,
        nt_hits=[Hit("OZ375105.1", "MAG: uncultured Alphaproteobacteria", 94.2, 1e-14, 78.0)],
        protein_hits=parse_blast_descriptions(K25_BLASTX_REPORT),
    )
    verdict = ident.verdict()
    assert verdict.startswith("protein-coding")
    assert "FAD-binding protein" in verdict
    assert "97% aa identity" in verdict
    assert ident.resolved


def test_verdict_reports_rfam_family_first_when_present():
    ident = Identification(
        query_id="x", length=120,
        nt_hits=[Hit("A1", "some genome", 50.0, 1e-5, 90.0)],
        rfam_hits=[{"id": "RF00001"}],
    )
    assert "Rfam: RF00001" in ident.verdict()


def test_verdict_falls_back_to_nucleotide_match():
    ident = Identification(
        query_id="x", length=120,
        nt_hits=[Hit("A1", "Humisphaera borealis chromosome", 154.0, 8e-33, 74.0)],
    )
    assert ident.verdict().startswith("matches known sequence")


def test_no_match_is_distinguished_from_lookup_failure():
    clean = Identification(query_id="k141_10", length=422)
    assert clean.verdict() == "no significant database match"
    assert not clean.resolved

    broken = Identification(query_id="k141_10", length=422, errors=["blastn: timed out"])
    assert broken.verdict() == "unresolved (lookup failed)"


def test_to_dict_round_trips_the_verdict():
    ident = Identification(
        query_id="k141_25", length=186,
        protein_hits=parse_blast_descriptions(K25_BLASTX_REPORT),
    )
    payload = ident.to_dict()
    assert payload["query_id"] == "k141_25"
    assert payload["length"] == 186
    assert payload["verdict"] == ident.verdict()
    assert payload["protein_hits"][0]["accession"] == "HZP76702.1"


def test_marginal_nucleotide_hits_still_parse():
    hits = parse_blast_descriptions(K10_BLASTN_REPORT)
    assert len(hits) == 2
    assert hits[0].evalue == 1.1
