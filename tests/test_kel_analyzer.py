"""Unit + smoke tests for KELAnalyzer against the NM_000420.3 cDNA reference.

Coverage:
  - Reference loads and the NM_000420.3 ATG offset is detected.
  - The file's base at c.578 matches the audit comment (regression sentinel).
  - Pure k-allele synthesis -> phenotype 'kk', single ISBT haplotype option.
  - Pure K-allele synthesis (c.578 C->T) -> phenotype 'KK'.
  - IUPAC heterozygote (c.578 = Y) -> phenotype 'Kk' with both alleles.
  - Low-identity random read -> uncallable, no SNPs trusted.
  - Reverse-complement input -> same phenotype as forward.
  - Multi-read concordance -> HIGH confidence at the primary SNP.
  - Real patient KEL56 amplicons -> non-Indeterminate phenotype (skipped when
    the patient data directory is absent so the file is portable).
"""

import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from Bio import SeqIO  # noqa: E402
from Bio.Seq import Seq  # noqa: E402

from utils.kel_analyzer import (  # noqa: E402
    KELAnalyzer,
    KEL_DIAGNOSTIC_SNPS,
)


# --------------------------------------------------------------------------- #
# Fixtures                                                                    #
# --------------------------------------------------------------------------- #

REPO_ROOT = Path(__file__).resolve().parent.parent
FASTA_PATH = REPO_ROOT / "utils" / "data" / "kel_cdna_reference.fasta"

PATIENT_ROOT = Path(r"C:\Users\ExPertComputer\Desktop\blood_group")
PATIENT_SAMPLES = ["Nampeung", "WHOreference 1", "WHOreference 5"]


@pytest.fixture(scope="module")
def analyzer() -> KELAnalyzer:
    a = KELAnalyzer()
    assert a.is_loaded(), "NM_000420.3 FASTA missing — cannot run KEL tests"
    assert a.reference_kind == "cdna"
    return a


@pytest.fixture(scope="module")
def reference_seq(analyzer: KELAnalyzer) -> str:
    return analyzer.reference_seq


def _patch(seq: str, cdna_pos: int, new_base: str, offset: int) -> str:
    idx = (cdna_pos - 1) + offset
    return seq[:idx] + new_base + seq[idx + 1 :]


# --------------------------------------------------------------------------- #
# Reference sanity                                                            #
# --------------------------------------------------------------------------- #

def test_reference_loads_as_nm000420_cdna(analyzer: KELAnalyzer):
    assert analyzer.reference_kind == "cdna"
    ctx = analyzer.reference_seq[analyzer.cdna_offset : analyzer.cdna_offset + 15]
    assert ctx == "ATGGAAGGTGGGGAC", f"unexpected ATG context: {ctx}"


def test_file_base_at_c578_matches_audit_note(analyzer: KELAnalyzer, reference_seq: str):
    """NM_000420.3 carries the k-allele base C at c.578 (no parity inversion needed).

    If NCBI ever revises the record and this changes, swap ref_base/alt_base
    + ref_call/alt_call in KEL_DIAGNOSTIC_SNPS and add parity_inverted: True
    — mirroring the RHCE c.178 treatment.
    """
    info = KEL_DIAGNOSTIC_SNPS["c.578C>T"]
    idx = (info["cDNA_position"] - 1) + analyzer.cdna_offset
    assert reference_seq[idx] == "C", (
        f"c.578 file base is {reference_seq[idx]}, expected C — see audit "
        f"note in utils/data/KEL_REFERENCE_README.md"
    )


# --------------------------------------------------------------------------- #
# Pure homozygous calls                                                       #
# --------------------------------------------------------------------------- #

def test_pure_k_allele_homozygous_gives_kk(analyzer: KELAnalyzer, reference_seq: str):
    # Reference is already a k-allele; no patching needed.
    result = analyzer.analyze([("pure_k", reference_seq)])

    assert result["phenotype"] == "kk"
    assert result["k_axis_call"]["genotype"] == "kk"
    assert result["allele_options"] == [
        {"haplotypes": "k/k", "isbt": "KEL*02/KEL*02", "serology": "K-k+"},
    ]
    assert result["k_axis_call"]["discordant_snps"] == []
    primary = result["snp_consensus"]["c.578C>T"]
    assert primary["consensus"] == "hom_ref"
    assert primary["call"] == "k"


def test_pure_capital_K_homozygous_gives_KK(analyzer: KELAnalyzer, reference_seq: str):
    # Patch c.578 C->T so the primary marker reads hom_alt.
    seq = _patch(reference_seq, 578, "T", analyzer.cdna_offset)
    result = analyzer.analyze([("pure_K", seq)])

    assert result["phenotype"] == "KK"
    assert result["k_axis_call"]["genotype"] == "KK"
    assert result["allele_options"] == [
        {"haplotypes": "K/K", "isbt": "KEL*01/KEL*01", "serology": "K+k-"},
    ]
    primary = result["snp_consensus"]["c.578C>T"]
    assert primary["consensus"] == "hom_alt"
    assert primary["call"] == "K"


# --------------------------------------------------------------------------- #
# Heterozygous (IUPAC) decoding                                               #
# --------------------------------------------------------------------------- #

def test_iupac_het_at_c578_gives_Kk(analyzer: KELAnalyzer, reference_seq: str):
    # Y = C/T het at c.578.
    seq = _patch(reference_seq, 578, "Y", analyzer.cdna_offset)
    result = analyzer.analyze([("het_c578", seq)])

    assert result["phenotype"] == "Kk"
    assert result["k_axis_call"]["genotype"] == "Kk"
    assert result["allele_options"] == [
        {"haplotypes": "K/k", "isbt": "KEL*01/KEL*02", "serology": "K+k+"},
    ]
    primary = result["snp_consensus"]["c.578C>T"]
    assert primary["consensus"] == "het"


# --------------------------------------------------------------------------- #
# Robustness                                                                  #
# --------------------------------------------------------------------------- #

def test_random_low_identity_read_is_uncallable(analyzer: KELAnalyzer):
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


def test_reverse_complement_input_matches_forward(
    analyzer: KELAnalyzer, reference_seq: str
):
    fwd = reference_seq
    rev = str(Seq(fwd).reverse_complement())

    result_fwd = analyzer.analyze([("fwd", fwd)])
    result_rev = analyzer.analyze([("rev", rev)])

    assert result_fwd["phenotype"] == result_rev["phenotype"] == "kk"
    # The probe caller is what makes RC work — record the strand it found.
    rev_call = result_rev["snp_consensus"]["c.578C>T"]
    assert rev_call["consensus"] == "hom_ref"


def test_multi_read_concordance_yields_high_confidence(
    analyzer: KELAnalyzer, reference_seq: str
):
    result = analyzer.analyze([
        ("read_a", reference_seq),
        ("read_b", reference_seq),
        ("read_c", reference_seq),
    ])
    assert result["snp_consensus"]["c.578C>T"]["confidence"] == "HIGH"
    assert result["reads_callable"] == 3
    assert result["overall_confidence"] == "HIGH"


# --------------------------------------------------------------------------- #
# Phred Q-score gate at SNP column                                            #
# --------------------------------------------------------------------------- #

def test_phred_gate_blocks_call_when_snp_base_below_q30(
    analyzer: KELAnalyzer, reference_seq: str
):
    """Per-base Phred at the SNP column < Q30 must produce no_call with a
    'low Phred at SNP position' reason, even if the base itself is correct.

    This is the lab-standard quality gate: a Q20 base is reliable enough
    for alignment context but not for variant calling at the called site.
    """
    # All bases except c.578 set to Q40 (high). The c.578 base set to Q15
    # (low). The base call itself is still 'C' (matches reference), so
    # without Phred gating this would call 'kk' hom_ref normally.
    info = KEL_DIAGNOSTIC_SNPS["c.578C>T"]
    snp_idx = (info["cDNA_position"] - 1) + analyzer.cdna_offset
    quality = [40] * len(reference_seq)
    quality[snp_idx] = 15

    result = analyzer.analyze([("low_q_at_snp", reference_seq, quality)])
    primary = result["snp_consensus"]["c.578C>T"]
    assert primary["consensus"] == "no_call", (
        f"Expected Phred gate to produce no_call, got {primary['consensus']}"
    )
    # Per-read detail should record the Phred Q-score that triggered the gate.
    per_read = result["per_read_details"][0]
    snp_call = per_read["snp_calls"]["c.578C>T"]
    assert snp_call.get("phred_at_snp") == 15
    assert "Phred" in snp_call.get("reason", ""), snp_call.get("reason")
    assert per_read.get("low_phred_at_snp") == 1


def test_phred_gate_passes_call_when_snp_base_at_q30_or_above(
    analyzer: KELAnalyzer, reference_seq: str
):
    """Q30 at the SNP column should pass the gate and call normally."""
    info = KEL_DIAGNOSTIC_SNPS["c.578C>T"]
    snp_idx = (info["cDNA_position"] - 1) + analyzer.cdna_offset
    quality = [40] * len(reference_seq)
    quality[snp_idx] = 30  # exactly at threshold — should pass

    result = analyzer.analyze([("q30_at_snp", reference_seq, quality)])
    primary = result["snp_consensus"]["c.578C>T"]
    assert primary["consensus"] == "hom_ref"
    assert result["phenotype"] == "kk"


def test_no_quality_supplied_skips_phred_gate_backward_compat(
    analyzer: KELAnalyzer, reference_seq: str
):
    """2-tuple (rid, seq) input — no quality — must still work as before."""
    result = analyzer.analyze([("legacy_2tuple", reference_seq)])
    assert result["phenotype"] == "kk"
    # phred_at_snp should be absent (not added) when no quality was supplied.
    snp_call = result["per_read_details"][0]["snp_calls"]["c.578C>T"]
    assert snp_call.get("phred_at_snp") is None
    assert result["per_read_details"][0].get("low_phred_at_snp") == 0


# --------------------------------------------------------------------------- #
# Real patient KEL56 amplicons (smoke test, skipped if data missing)          #
# --------------------------------------------------------------------------- #

def _patient_kel_fastas(sample: str) -> list:
    """Return list of (read_id, sequence) tuples for KEL56 FASTAs."""
    out = []
    d = PATIENT_ROOT / sample / "KEL56"
    if not d.is_dir():
        return out
    for fa in sorted(d.glob("*.fasta")):
        try:
            rec = next(SeqIO.parse(str(fa), "fasta"))
        except StopIteration:
            continue
        out.append((f"{sample}/KEL56/{fa.stem}", str(rec.seq)))
    return out


@pytest.mark.skipif(
    not PATIENT_ROOT.exists(),
    reason=f"Real patient data not present at {PATIENT_ROOT}",
)
@pytest.mark.parametrize("sample", PATIENT_SAMPLES)
def test_real_patient_analyze_smoke(analyzer: KELAnalyzer, sample: str):
    reads = _patient_kel_fastas(sample)
    if not reads:
        pytest.skip(f"No KEL56 FASTAs found for sample {sample}")

    result = analyzer.analyze(reads)

    assert set(result) >= {
        "phenotype",
        "allele_options",
        "k_axis_call",
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
    # KEL56 amplicons cover c.578, so the primary marker must resolve.
    assert result["phenotype"] != "Indeterminate", (
        f"{sample} indeterminate: reason={result['reason']}; "
        f"primary={result['snp_consensus']['c.578C>T']}"
    )
    assert result["overall_confidence"] != "NONE", (
        f"{sample} overall_confidence still NONE: {result['reason']}"
    )
