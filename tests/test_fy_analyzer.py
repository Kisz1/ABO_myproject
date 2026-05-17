"""Unit + smoke tests for FYAnalyzer against the NG_011626.3 genomic reference.

Coverage:
  - Reference loads as genomic, ACKR1 features (not CADM3) are picked up.
  - Each diagnostic SNP's reference base matches the audit comment (regression
    sentinel against NCBI revisions or paralog mis-routing).
  - Pure FY*A (reference) synthesis -> phenotype Fy(a+b-), single haplotype.
  - Pure FY*B synthesis (c.125 G->A) -> phenotype Fy(a-b+).
  - IUPAC het at c.125 (R = G/A) -> phenotype Fy(a+b+).
  - GATA hom-alt (-67 T->C) -> phenotype Fy(a-b-) regardless of A/B.
  - GATA het + A/B het -> two phase-ambiguous allele options listed.
  - FY*X (c.265 het) on top of FY*B -> "weak" tag in serology.
  - Reverse-complement input -> same phenotype as forward (probe-based caller).
  - Multi-read concordance -> HIGH confidence at the primary SNP.
  - cDNA-only reference -> GATA SNP is uncoverable with a clear reason.
  - Low-identity random read -> uncallable, no SNPs trusted.
  - Real patient FY12 amplicons -> non-Indeterminate phenotype (skipped when
    the patient data directory is absent so the file is portable).
"""

import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from Bio import SeqIO  # noqa: E402
from Bio.Seq import Seq  # noqa: E402

from utils.fy_analyzer import (  # noqa: E402
    FYAnalyzer,
    FYReferenceMissingError,
    FY_DIAGNOSTIC_SNPS,
)


# --------------------------------------------------------------------------- #
# Fixtures                                                                    #
# --------------------------------------------------------------------------- #

REPO_ROOT = Path(__file__).resolve().parent.parent
GENBANK_PATH = REPO_ROOT / "utils" / "data" / "fy_referance.gb"

PATIENT_ROOT = Path(r"C:\Users\ExPertComputer\Desktop\blood_group")
PATIENT_SAMPLES = ["Nampeung", "WHOreference 1", "WHOreference 5"]


@pytest.fixture(scope="module")
def analyzer() -> FYAnalyzer:
    a = FYAnalyzer()
    assert a.is_loaded(), "NG_011626.3 GenBank missing — cannot run FY tests"
    assert a.reference_kind == "genomic"
    return a


@pytest.fixture(scope="module")
def reference_seq(analyzer: FYAnalyzer) -> str:
    return analyzer.reference_seq


def _patch_at_genomic_index(seq: str, genomic_index: int, new_base: str) -> str:
    """Replace one base in the genomic reference at a 0-based index."""
    return seq[:genomic_index] + new_base + seq[genomic_index + 1:]


def _patch_snp(analyzer: FYAnalyzer, seq: str, snp_name: str, new_base: str) -> str:
    """Patch the analyzer-located position for a diagnostic SNP."""
    info = FY_DIAGNOSTIC_SNPS[snp_name]
    idx = analyzer._position_to_ref_index(info)
    assert idx is not None, f"could not locate {snp_name} in reference"
    return _patch_at_genomic_index(seq, idx, new_base)


# --------------------------------------------------------------------------- #
# Reference sanity                                                            #
# --------------------------------------------------------------------------- #

def test_reference_loads_as_ng011626_genomic(analyzer: FYAnalyzer):
    assert analyzer.reference_kind == "genomic"
    # NG_011626.3 has 2 ACKR1 exons; the CDS start must be ACKR1's, not the
    # neighbouring CADM3's (which would put cds_start at ~738 instead).
    assert len(analyzer.exon_map) == 2, analyzer.exon_map
    assert analyzer.cds_start_genomic == 5947, (
        f"cds_start_genomic={analyzer.cds_start_genomic} — looks like the "
        "loader picked up the wrong gene (CADM3?) from the RefSeqGene record."
    )


def test_file_base_at_c125_matches_audit_note(analyzer: FYAnalyzer, reference_seq: str):
    """NG_011626.3 carries the FY*A base G at c.125 (no parity inversion)."""
    idx = analyzer._position_to_ref_index(FY_DIAGNOSTIC_SNPS["c.125G>A"])
    assert reference_seq[idx] == "G", (
        f"c.125 base in NG_011626.3 is {reference_seq[idx]!r}, expected 'G' — "
        "see audit note in utils/data/FY_REFERENCE_README.md"
    )


def test_file_base_at_gata_matches_audit_note(analyzer: FYAnalyzer, reference_seq: str):
    """NG_011626.3 carries the wild-type GATA-1 site base T at -67."""
    idx = analyzer._position_to_ref_index(FY_DIAGNOSTIC_SNPS["-67T>C (GATA)"])
    assert reference_seq[idx] == "T", (
        f"GATA -67 base is {reference_seq[idx]!r}, expected 'T'"
    )


def test_file_base_at_c265_matches_audit_note(analyzer: FYAnalyzer, reference_seq: str):
    """NG_011626.3 carries the FY*X reference (non-weak) base C at c.265."""
    idx = analyzer._position_to_ref_index(FY_DIAGNOSTIC_SNPS["c.265C>T"])
    assert reference_seq[idx] == "C", (
        f"c.265 base is {reference_seq[idx]!r}, expected 'C'"
    )


# --------------------------------------------------------------------------- #
# Pure homozygous calls                                                       #
# --------------------------------------------------------------------------- #

def test_pure_fya_reference_gives_fy_a_pos_b_neg(analyzer: FYAnalyzer, reference_seq: str):
    # Reference is already FY*A; no patching needed.
    result = analyzer.analyze([("pure_fya", reference_seq)])

    assert result["phenotype"] == "Fy(a+b-)"
    assert result["a_b_call"]["genotype"] == "AA"
    assert result["allele_options"] == [
        {"haplotypes": "FY*A/FY*A", "isbt": "FY*01/FY*01", "serology": "Fy(a+b-)"},
    ]
    primary = result["snp_consensus"]["c.125G>A"]
    assert primary["consensus"] == "hom_ref"
    assert primary["call"] == "A"


def test_pure_fyb_homozygous_gives_fy_a_neg_b_pos(analyzer: FYAnalyzer, reference_seq: str):
    seq = _patch_snp(analyzer, reference_seq, "c.125G>A", "A")
    result = analyzer.analyze([("pure_fyb", seq)])

    assert result["phenotype"] == "Fy(a-b+)"
    assert result["a_b_call"]["genotype"] == "BB"
    assert result["allele_options"][0]["isbt"] == "FY*02/FY*02"
    primary = result["snp_consensus"]["c.125G>A"]
    assert primary["consensus"] == "hom_alt"
    assert primary["call"] == "B"


# --------------------------------------------------------------------------- #
# Heterozygous (IUPAC) decoding                                               #
# --------------------------------------------------------------------------- #

def test_iupac_het_at_c125_gives_fy_a_pos_b_pos(analyzer: FYAnalyzer, reference_seq: str):
    # R = G/A het at c.125.
    seq = _patch_snp(analyzer, reference_seq, "c.125G>A", "R")
    result = analyzer.analyze([("het_c125", seq)])

    assert result["phenotype"] == "Fy(a+b+)"
    assert result["a_b_call"]["genotype"] == "AB"
    assert result["allele_options"][0]["isbt"] == "FY*01/FY*02"


# --------------------------------------------------------------------------- #
# GATA modifier — FY*Null                                                     #
# --------------------------------------------------------------------------- #

def test_gata_hom_alt_overrides_to_fy_null(analyzer: FYAnalyzer, reference_seq: str):
    """Both haplotypes carry the GATA mutation -> Fy(a-b-) regardless of A/B."""
    # Leave c.125 = FY*A (so we'd otherwise call Fy(a+b-)); flip GATA T->C.
    seq = _patch_snp(analyzer, reference_seq, "-67T>C (GATA)", "C")
    result = analyzer.analyze([("gata_hom", seq)])

    assert result["phenotype"] == "Fy(a-b-)"
    assert result["allele_options"][0]["haplotypes"] == "FY*0/FY*0"
    assert "silenced" in result["allele_options"][0]["serology"].lower()


def test_gata_het_with_ab_het_yields_phase_ambiguous_options(
    analyzer: FYAnalyzer, reference_seq: str
):
    """GATA het + A/B het cannot be phased by Sanger — both options listed."""
    seq = _patch_snp(analyzer, reference_seq, "c.125G>A", "R")          # A/B het
    seq = _patch_snp(analyzer, seq, "-67T>C (GATA)", "Y")               # GATA het
    result = analyzer.analyze([("gata_phase_ambig", seq)])

    assert len(result["allele_options"]) == 2, result["allele_options"]
    isbts = {opt["isbt"] for opt in result["allele_options"]}
    assert isbts == {"FY*01N.01 / FY*02", "FY*01 / FY*02N.01"}
    assert "phase ambiguous" in result["phenotype"].lower()


# --------------------------------------------------------------------------- #
# FY*X weak-B modifier                                                        #
# --------------------------------------------------------------------------- #

def test_fy_x_present_on_fyb_carrier_adds_weak_tag(
    analyzer: FYAnalyzer, reference_seq: str
):
    """A/B = het (Fy(a+b+)) + c.265 het (FY*X carrier) -> weak-B tag in serology."""
    seq = _patch_snp(analyzer, reference_seq, "c.125G>A", "R")          # A/B het
    seq = _patch_snp(analyzer, seq, "c.265C>T", "Y")                    # FY*X het
    result = analyzer.analyze([("fy_x_het", seq)])

    assert result["a_b_call"]["genotype"] == "AB"
    serology = result["allele_options"][0]["serology"].lower()
    assert "weak" in serology and "fy*x" in serology


# --------------------------------------------------------------------------- #
# Robustness                                                                  #
# --------------------------------------------------------------------------- #

def test_reverse_complement_input_matches_forward(
    analyzer: FYAnalyzer, reference_seq: str
):
    fwd = reference_seq
    rev = str(Seq(fwd).reverse_complement())

    result_fwd = analyzer.analyze([("fwd", fwd)])
    result_rev = analyzer.analyze([("rev", rev)])

    assert result_fwd["phenotype"] == result_rev["phenotype"] == "Fy(a+b-)"


def test_multi_read_concordance_yields_high_confidence(
    analyzer: FYAnalyzer, reference_seq: str
):
    result = analyzer.analyze([
        ("read_a", reference_seq),
        ("read_b", reference_seq),
        ("read_c", reference_seq),
    ])
    assert result["snp_consensus"]["c.125G>A"]["confidence"] == "HIGH"
    assert result["reads_callable"] == 3
    assert result["overall_confidence"] == "HIGH"


# --------------------------------------------------------------------------- #
# Phred Q-score gate at SNP column                                            #
# --------------------------------------------------------------------------- #

def test_phred_gate_blocks_call_when_snp_base_below_q30(
    analyzer: FYAnalyzer, reference_seq: str
):
    info = FY_DIAGNOSTIC_SNPS["c.125G>A"]
    snp_idx = analyzer._position_to_ref_index(info)
    quality = [40] * len(reference_seq)
    quality[snp_idx] = 15
    result = analyzer.analyze([("low_q_at_snp", reference_seq, quality)])
    primary = result["snp_consensus"]["c.125G>A"]
    assert primary["consensus"] == "no_call"
    snp_call = result["per_read_details"][0]["snp_calls"]["c.125G>A"]
    assert snp_call.get("phred_at_snp") == 15
    assert "Phred" in snp_call.get("reason", "")
    assert result["per_read_details"][0].get("low_phred_at_snp") == 1


def test_no_quality_supplied_skips_phred_gate_backward_compat(
    analyzer: FYAnalyzer, reference_seq: str
):
    result = analyzer.analyze([("legacy_2tuple", reference_seq)])
    assert result["a_b_call"]["genotype"] == "AA"
    snp_call = result["per_read_details"][0]["snp_calls"]["c.125G>A"]
    assert snp_call.get("phred_at_snp") is None


def test_random_low_identity_read_is_uncallable(analyzer: FYAnalyzer):
    rng = random.Random(2026)
    junk = "".join(rng.choices("ACGT", k=1000))
    result = analyzer.analyze_single_read(junk, read_id="junk")
    assert result["callable"] is False
    for call in result["snp_calls"].values():
        if call.get("covered"):
            assert call.get("call") == "no_call", call
        else:
            assert call.get("covered") is False
    assert result["trusted_snps"] == 0


# --------------------------------------------------------------------------- #
# cDNA-only reference: GATA SNP must surface as uncoverable                   #
# --------------------------------------------------------------------------- #

def test_gata_snp_is_uncoverable_against_cdna_reference():
    """Mimic a cDNA-only reference (no promoter coverage) — GATA must return
    covered=False with a clear reason, not a silent miscalled hom_ref."""
    # Synthetic cDNA: a chunk of upstream-of-ATG context + ATG + 300 bp of CDS.
    # The promoter SNP at -67 from ATG is intentionally NOT in this sequence.
    cdna = (
        "ATGGGCAACTGCCTGCACCGCGCCGAGCTGAGCCCGAGCACG"
        + "ACGTGCCCAGGCTGCGCACGCAGGAGCACG" * 10
    )
    a = FYAnalyzer(reference_seq=cdna)
    # Sanity: reference_kind reset to cdna; cds_start_genomic stays None.
    assert a.reference_kind == "cdna"
    assert a.cds_start_genomic is None

    info = FY_DIAGNOSTIC_SNPS["-67T>C (GATA)"]
    assert a._position_to_ref_index(info) is None, (
        "GATA SNP should be unreachable against a cDNA-only reference"
    )


# --------------------------------------------------------------------------- #
# Missing reference raises a clear, instruction-bearing error                 #
# --------------------------------------------------------------------------- #

def test_missing_reference_raises_with_download_hint(tmp_path: Path):
    a = FYAnalyzer(gb_path=str(tmp_path / "absent.gb"),
                   fasta_path=str(tmp_path / "absent.fasta"))
    assert not a.is_loaded()
    with pytest.raises(FYReferenceMissingError) as exc_info:
        a.analyze([("x", "ACGT" * 50)])
    msg = str(exc_info.value)
    assert "FY_REFERENCE_README" in msg
    assert "GenBank" in msg or "FASTA" in msg


# --------------------------------------------------------------------------- #
# Real patient FY12 amplicons (smoke test, skipped if data missing)           #
# --------------------------------------------------------------------------- #

def _patient_fy_fastas(sample: str) -> list:
    """Return list of (read_id, sequence) tuples for FY12 FASTAs."""
    out = []
    d = PATIENT_ROOT / sample / "FY12"
    if not d.is_dir():
        return out
    for fa in sorted(d.glob("*.fasta")):
        try:
            rec = next(SeqIO.parse(str(fa), "fasta"))
        except StopIteration:
            continue
        out.append((f"{sample}/FY12/{fa.stem}", str(rec.seq)))
    return out


@pytest.mark.skipif(
    not PATIENT_ROOT.exists(),
    reason=f"Real patient data not present at {PATIENT_ROOT}",
)
@pytest.mark.parametrize("sample", PATIENT_SAMPLES)
def test_real_patient_analyze_smoke(analyzer: FYAnalyzer, sample: str):
    reads = _patient_fy_fastas(sample)
    if not reads:
        pytest.skip(f"No FY12 FASTAs found for sample {sample}")

    result = analyzer.analyze(reads)

    assert set(result) >= {
        "phenotype",
        "allele_options",
        "a_b_call",
        "gata_call",
        "fy_x_call",
        "snp_consensus",
        "overall_confidence",
        "reads_total",
        "reads_callable",
        "reads_with_trusted_snps",
        "per_read_details",
    }
    assert result["reads_total"] == len(reads)
    assert isinstance(result["phenotype"], str) and result["phenotype"]

    assert result["reads_with_trusted_snps"] >= 1, (
        f"No SNP cleared the probe-identity gate for {sample}: "
        f"{[(r['read_id'], r['identity'], r.get('trusted_snps')) for r in result['per_read_details']]}"
    )
    # FY12 amplicons cover c.125 (primary), so the A/B axis must resolve.
    assert result["phenotype"] != "Indeterminate", (
        f"{sample} indeterminate: reason={result['reason']}; "
        f"primary={result['snp_consensus']['c.125G>A']}"
    )
    assert result["overall_confidence"] != "NONE", (
        f"{sample} overall_confidence still NONE: {result['reason']}"
    )
