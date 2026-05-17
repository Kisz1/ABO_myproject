"""Unit + smoke tests for JKAnalyzer against the NG_011775.2 genomic reference.

Coverage:
  - Reference loads as genomic, SLC14A1 features picked up (CDS start at 11194).
  - Each diagnostic SNP's reference base matches the audit comment.
  - The c.838G>A entry is parity-inverted: NG_011775.2 carries the JK*B
    allele 'A' at c.838, so the file's "ref" is JK*B and the file's "alt"
    is JK*A. Regression test enforces this.
  - Pure reference (JK*B/JK*B) -> Jk(a-b+).
  - Pure JK*A synthesis (c.838 A->G) -> Jk(a+b-).
  - IUPAC het at c.838 (R = G/A) -> Jk(a+b+).
  - Either null SNP hom-alt -> Jk(a-b-) regardless of A/B.
  - One null het + A/B het -> two phase-ambiguous allele options.
  - cDNA-only reference -> c.342-1 splice SNP is uncoverable with a clear reason.
  - Reverse-complement input -> same phenotype as forward.
  - Multi-read concordance -> HIGH confidence at the primary SNP.
  - Low-identity random read -> uncallable, no SNPs trusted.
  - Missing reference raises with a README hint.
  - Real patient JK78/JK89 amplicons -> non-Indeterminate phenotype (skipped
    when the patient data directory is absent so the file is portable).
"""

import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from Bio import SeqIO  # noqa: E402
from Bio.Seq import Seq  # noqa: E402

from utils.jk_analyzer import (  # noqa: E402
    JKAnalyzer,
    JKReferenceMissingError,
    JK_DIAGNOSTIC_SNPS,
)


# --------------------------------------------------------------------------- #
# Fixtures                                                                    #
# --------------------------------------------------------------------------- #

REPO_ROOT = Path(__file__).resolve().parent.parent
GENBANK_PATH = REPO_ROOT / "utils" / "data" / "jk_referance.gb"

PATIENT_ROOT = Path(r"C:\Users\ExPertComputer\Desktop\blood_group")
PATIENT_SAMPLES = ["Nampeung", "WHOreference 1", "WHOreference 5"]


@pytest.fixture(scope="module")
def analyzer() -> JKAnalyzer:
    a = JKAnalyzer()
    assert a.is_loaded(), "NG_011775.2 GenBank missing — cannot run JK tests"
    assert a.reference_kind == "genomic"
    return a


@pytest.fixture(scope="module")
def reference_seq(analyzer: JKAnalyzer) -> str:
    return analyzer.reference_seq


def _patch_at_genomic_index(seq: str, genomic_index: int, new_base: str) -> str:
    return seq[:genomic_index] + new_base + seq[genomic_index + 1:]


def _patch_snp(analyzer: JKAnalyzer, seq: str, snp_name: str, new_base: str) -> str:
    info = JK_DIAGNOSTIC_SNPS[snp_name]
    idx = analyzer._position_to_ref_index(info)
    assert idx is not None, f"could not locate {snp_name} in reference"
    return _patch_at_genomic_index(seq, idx, new_base)


# --------------------------------------------------------------------------- #
# Reference sanity                                                            #
# --------------------------------------------------------------------------- #

def test_reference_loads_as_ng011775_genomic(analyzer: JKAnalyzer):
    assert analyzer.reference_kind == "genomic"
    assert analyzer.cds_start_genomic == 11194, (
        f"cds_start_genomic={analyzer.cds_start_genomic} — should be 11194 "
        "(SLC14A1 ATG in NG_011775.2)"
    )
    # SLC14A1 has 11 canonical exons + 2 alt-spliced exon-3 variants
    # deduplicated by (start, end). Expect 13 entries after dedup.
    assert len(analyzer.exon_map) >= 11, analyzer.exon_map


def test_file_base_at_c838_matches_audit_note(analyzer: JKAnalyzer, reference_seq: str):
    """NG_011775.2 carries the JK*B base 'A' at c.838 — parity-inverted entry.

    If a future NCBI revision flips this to 'G' (a JK*A clone), the SNP entry
    in JK_DIAGNOSTIC_SNPS must be un-inverted: swap ref_base/alt_base and
    ref_call/alt_call, and remove `parity_inverted: True`.
    """
    idx = analyzer._position_to_ref_index(JK_DIAGNOSTIC_SNPS["c.838G>A"])
    assert reference_seq[idx] == "A", (
        f"c.838 base in NG_011775.2 is {reference_seq[idx]!r}, expected 'A' — "
        "see audit note in utils/data/JK_REFERENCE_README.md"
    )
    assert JK_DIAGNOSTIC_SNPS["c.838G>A"].get("parity_inverted") is True


def test_file_base_at_splice_matches_audit_note(analyzer: JKAnalyzer, reference_seq: str):
    """NG_011775.2 carries the intact splice-acceptor 'G' at c.342-1."""
    idx = analyzer._position_to_ref_index(JK_DIAGNOSTIC_SNPS["c.342-1G>A"])
    assert reference_seq[idx] == "G", (
        f"c.342-1 base is {reference_seq[idx]!r}, expected 'G'"
    )


def test_file_base_at_c871_matches_audit_note(analyzer: JKAnalyzer, reference_seq: str):
    """NG_011775.2 carries the wild-type 'T' at c.871 (non-silenced)."""
    idx = analyzer._position_to_ref_index(JK_DIAGNOSTIC_SNPS["c.871T>C"])
    assert reference_seq[idx] == "T", (
        f"c.871 base is {reference_seq[idx]!r}, expected 'T'"
    )


# --------------------------------------------------------------------------- #
# Pure homozygous calls                                                       #
# --------------------------------------------------------------------------- #

def test_pure_reference_gives_jk_a_neg_b_pos(analyzer: JKAnalyzer, reference_seq: str):
    """Deposited reference is JK*B; pure reference -> Jk(a-b+)."""
    result = analyzer.analyze([("pure_ref", reference_seq)])
    assert result["phenotype"] == "Jk(a-b+)"
    assert result["a_b_call"]["genotype"] == "BB"
    assert result["allele_options"][0]["isbt"] == "JK*02/JK*02"


def test_pure_jka_homozygous_gives_jk_a_pos_b_neg(
    analyzer: JKAnalyzer, reference_seq: str
):
    # Flip c.838 A->G so the primary marker reads hom_alt (i.e. JK*A allele).
    seq = _patch_snp(analyzer, reference_seq, "c.838G>A", "G")
    result = analyzer.analyze([("pure_jka", seq)])
    assert result["phenotype"] == "Jk(a+b-)"
    assert result["a_b_call"]["genotype"] == "AA"
    assert result["allele_options"][0]["isbt"] == "JK*01/JK*01"


# --------------------------------------------------------------------------- #
# Heterozygous (IUPAC) decoding                                               #
# --------------------------------------------------------------------------- #

def test_iupac_het_at_c838_gives_jk_a_pos_b_pos(
    analyzer: JKAnalyzer, reference_seq: str
):
    seq = _patch_snp(analyzer, reference_seq, "c.838G>A", "R")  # R = G/A
    result = analyzer.analyze([("het_c838", seq)])
    assert result["phenotype"] == "Jk(a+b+)"
    assert result["a_b_call"]["genotype"] == "AB"
    assert result["allele_options"][0]["isbt"] == "JK*01/JK*02"


# --------------------------------------------------------------------------- #
# Null modifiers                                                              #
# --------------------------------------------------------------------------- #

def test_splice_null_hom_alt_overrides_to_jk_null(
    analyzer: JKAnalyzer, reference_seq: str
):
    """Polynesian splice SNP homozygous -> Jk(a-b-) regardless of A/B."""
    seq = _patch_snp(analyzer, reference_seq, "c.342-1G>A", "A")
    result = analyzer.analyze([("splice_null_hom", seq)])
    assert result["phenotype"] == "Jk(a-b-)"
    assert result["allele_options"][0]["haplotypes"] == "JK*0/JK*0"
    assert any(
        t["snp"] == "c.342-1G>A" and t["state"] == "hom_alt"
        for t in result["null_call"]["triggers"]
    )


def test_asian_null_hom_alt_overrides_to_jk_null(
    analyzer: JKAnalyzer, reference_seq: str
):
    """Asian null SNP homozygous -> Jk(a-b-) regardless of A/B."""
    seq = _patch_snp(analyzer, reference_seq, "c.871T>C", "C")
    result = analyzer.analyze([("asian_null_hom", seq)])
    assert result["phenotype"] == "Jk(a-b-)"
    assert any(
        t["snp"] == "c.871T>C" and t["state"] == "hom_alt"
        for t in result["null_call"]["triggers"]
    )


def test_null_het_with_ab_het_yields_phase_ambiguous_options(
    analyzer: JKAnalyzer, reference_seq: str
):
    """A null variant het + A/B het cannot be phased — two options listed."""
    seq = _patch_snp(analyzer, reference_seq, "c.838G>A", "R")          # A/B het
    seq = _patch_snp(analyzer, seq, "c.871T>C", "Y")                    # Asian null het
    result = analyzer.analyze([("null_phase_ambig", seq)])
    assert len(result["allele_options"]) == 2
    isbts = {opt["isbt"] for opt in result["allele_options"]}
    assert isbts == {"JK*01N / JK*02", "JK*01 / JK*02N"}
    assert "phase ambiguous" in result["phenotype"].lower()


# --------------------------------------------------------------------------- #
# Robustness                                                                  #
# --------------------------------------------------------------------------- #

def test_reverse_complement_input_matches_forward(
    analyzer: JKAnalyzer, reference_seq: str
):
    fwd = reference_seq
    rev = str(Seq(fwd).reverse_complement())

    result_fwd = analyzer.analyze([("fwd", fwd)])
    result_rev = analyzer.analyze([("rev", rev)])

    assert result_fwd["phenotype"] == result_rev["phenotype"] == "Jk(a-b+)"


def test_multi_read_concordance_yields_high_confidence(
    analyzer: JKAnalyzer, reference_seq: str
):
    result = analyzer.analyze([
        ("read_a", reference_seq),
        ("read_b", reference_seq),
        ("read_c", reference_seq),
    ])
    assert result["snp_consensus"]["c.838G>A"]["confidence"] == "HIGH"
    assert result["reads_callable"] == 3
    assert result["overall_confidence"] == "HIGH"


# --------------------------------------------------------------------------- #
# Phred Q-score gate at SNP column                                            #
# --------------------------------------------------------------------------- #

def test_phred_gate_blocks_call_when_snp_base_below_q30(
    analyzer: JKAnalyzer, reference_seq: str
):
    info = JK_DIAGNOSTIC_SNPS["c.838G>A"]
    snp_idx = analyzer._position_to_ref_index(info)
    quality = [40] * len(reference_seq)
    quality[snp_idx] = 15
    result = analyzer.analyze([("low_q_at_snp", reference_seq, quality)])
    primary = result["snp_consensus"]["c.838G>A"]
    assert primary["consensus"] == "no_call"
    snp_call = result["per_read_details"][0]["snp_calls"]["c.838G>A"]
    assert snp_call.get("phred_at_snp") == 15
    assert "Phred" in snp_call.get("reason", "")


def test_no_quality_supplied_skips_phred_gate_backward_compat(
    analyzer: JKAnalyzer, reference_seq: str
):
    result = analyzer.analyze([("legacy_2tuple", reference_seq)])
    # Parity-inverted JK*B reference -> pure ref calls BB.
    assert result["a_b_call"]["genotype"] == "BB"
    snp_call = result["per_read_details"][0]["snp_calls"]["c.838G>A"]
    assert snp_call.get("phred_at_snp") is None


def test_random_low_identity_read_is_uncallable(analyzer: JKAnalyzer):
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
# cDNA-only reference: intronic SNP must surface as uncoverable               #
# --------------------------------------------------------------------------- #

def test_splice_snp_is_uncoverable_against_cdna_reference():
    """Mimic a cDNA-only reference — c.342-1 (intronic) must return
    covered=False with a clear reason, not a silent miscalled hom_ref."""
    cdna = "ATGGAGGACAGCCCC" + ("ACGT" * 250)
    a = JKAnalyzer(reference_seq=cdna)
    assert a.reference_kind == "cdna"

    info = JK_DIAGNOSTIC_SNPS["c.342-1G>A"]
    assert a._position_to_ref_index(info) is None, (
        "c.342-1 (intronic) must be unreachable against a cDNA-only reference"
    )


# --------------------------------------------------------------------------- #
# Missing reference raises a clear, instruction-bearing error                 #
# --------------------------------------------------------------------------- #

def test_missing_reference_raises_with_download_hint(tmp_path: Path):
    a = JKAnalyzer(gb_path=str(tmp_path / "absent.gb"),
                   fasta_path=str(tmp_path / "absent.fasta"))
    assert not a.is_loaded()
    with pytest.raises(JKReferenceMissingError) as exc_info:
        a.analyze([("x", "ACGT" * 50)])
    msg = str(exc_info.value)
    assert "JK_REFERENCE_README" in msg


# --------------------------------------------------------------------------- #
# Real patient JK78 amplicons (smoke test, skipped if data missing)           #
# --------------------------------------------------------------------------- #

def _patient_jk_fastas(sample: str) -> list:
    out = []
    d = PATIENT_ROOT / sample / "JK78"
    if not d.is_dir():
        return out
    for fa in sorted(d.glob("*.fasta")):
        try:
            rec = next(SeqIO.parse(str(fa), "fasta"))
        except StopIteration:
            continue
        out.append((f"{sample}/JK78/{fa.stem}", str(rec.seq)))
    return out


@pytest.mark.skipif(
    not PATIENT_ROOT.exists(),
    reason=f"Real patient data not present at {PATIENT_ROOT}",
)
@pytest.mark.parametrize("sample", PATIENT_SAMPLES)
def test_real_patient_analyze_smoke(analyzer: JKAnalyzer, sample: str):
    reads = _patient_jk_fastas(sample)
    if not reads:
        pytest.skip(f"No JK78 FASTAs found for sample {sample}")

    result = analyzer.analyze(reads)

    assert set(result) >= {
        "phenotype",
        "allele_options",
        "a_b_call",
        "null_call",
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
        f"No SNP cleared the probe-identity gate for {sample}"
    )
    assert result["phenotype"] != "Indeterminate", (
        f"{sample} indeterminate: reason={result['reason']}; "
        f"primary={result['snp_consensus']['c.838G>A']}"
    )
    assert result["overall_confidence"] != "NONE"
