"""Tests for the library-type preflight.

These exercise assess_library, which is pure -- no network is touched.
The SRR13291825 case is pinned explicitly because it is the run whose
library type this whole guard exists to have caught.
"""
from rnasig.srameta import RunMetadata, assess_library


def _rnaseq(**kw) -> RunMetadata:
    base = dict(
        accession="SRRTEST",
        library_strategy="RNA-Seq",
        library_selection="RANDOM",
        library_source="METATRANSCRIPTOMIC",
    )
    base.update(kw)
    return RunMetadata(**base)


def test_shotgun_metatranscriptome_is_accepted():
    verdict = assess_library(_rnaseq())
    assert verdict.compatible
    assert verdict.reasons == []
    assert bool(verdict) is True


def test_transcriptomic_source_is_accepted():
    assert assess_library(_rnaseq(library_source="TRANSCRIPTOMIC")).compatible


def test_srr13291825_amplicon_run_is_rejected():
    """The actual metadata for the run swept in Phase 3, from ENA."""
    meta = RunMetadata(
        accession="SRR13291825",
        library_strategy="AMPLICON",
        library_selection="PCR",
        library_source="METAGENOMIC",
        scientific_name="soil metagenome",
    )
    verdict = assess_library(meta)
    assert not verdict.compatible
    assert not verdict.unknown
    # All three independent problems should be reported, not just the first.
    joined = " ".join(verdict.reasons)
    assert "AMPLICON" in joined
    assert "PCR" in joined
    assert "METAGENOMIC" in joined
    assert len(verdict.reasons) == 3


def test_amplicon_strategy_alone_is_rejected():
    verdict = assess_library(_rnaseq(library_strategy="AMPLICON"))
    assert not verdict.compatible
    assert any("AMPLICON" in r for r in verdict.reasons)


def test_dna_source_alone_is_rejected():
    verdict = assess_library(_rnaseq(library_source="GENOMIC"))
    assert not verdict.compatible
    assert any("DNA" in r for r in verdict.reasons)


def test_pcr_selection_alone_is_rejected():
    verdict = assess_library(_rnaseq(library_selection="PCR"))
    assert not verdict.compatible
    assert any("primers" in r for r in verdict.reasons)


def test_wgs_strategy_is_rejected():
    assert not assess_library(_rnaseq(library_strategy="WGS")).compatible


def test_missing_metadata_is_not_treated_as_permission():
    """Absent metadata must not read as 'fine to sweep'."""
    verdict = assess_library(RunMetadata(accession="SRRUNKNOWN"))
    assert not verdict.compatible
    assert verdict.unknown
    assert verdict.reasons


def test_case_and_whitespace_are_normalised():
    verdict = assess_library(
        RunMetadata(
            accession="SRRTEST",
            library_strategy="  amplicon ",
            library_selection=" pcr",
            library_source="metagenomic ",
        )
    )
    assert not verdict.compatible
    assert len(verdict.reasons) == 3
