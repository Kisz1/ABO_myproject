"""Unit + smoke tests for MNSAnalyzer against NG_007470.3 (GYPA) and NG_007483.3 (GYPB).

Coverage:
  - Both references load; gene_refs has GYPA + GYPB entries with correct CDS starts.
  - Each diagnostic SNP's reference base matches the audit comment.
  - GYPA c.59 audit: ref = C (M allele).
  - GYPB c.143 audit: ref = C (s allele) — parity-inverted entry.
  - Pure GYPA reference -> M+N- on the M/N axis (S/s indeterminate, that's fine).
  - Patched c.59 C->T -> N/N (parity-safe genotype mapping).
  - Patched c.59 C->Y (IUPAC het) -> MN.
  - GYPB pure reference -> ss on S/s axis.
  - Patched c.143 C->T (against GYPB ref) -> SS.
  - Reverse-complement input -> same axis call.
  - Multi-read concordance -> HIGH confidence at the contributing SNP.
  - Random low-identity read -> uncallable.
  - Missing reference (one or both) -> raises with README hint.
  - Real patient MIA234 amplicons: smoke run only — does not assert on
    phenotype since MIA234 may not cover the M/N or S/s SNP positions.
"""

import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from Bio import SeqIO  # noqa: E402
from Bio.Seq import Seq  # noqa: E402

from utils.mns_analyzer import (  # noqa: E402
    MNSAnalyzer,
    MNSReferenceMissingError,
    MNS_DIAGNOSTIC_SNPS,
)


REPO_ROOT = Path(__file__).resolve().parent.parent
GYPA_PATH = REPO_ROOT / "utils" / "data" / "mns_gypa_referance.gb"
GYPB_PATH = REPO_ROOT / "utils" / "data" / "mns_gypb_referance.gb"

PATIENT_ROOT = Path(r"C:\Users\ExPertComputer\Desktop\blood_group")
PATIENT_SAMPLES = ["Nampeung", "WHOreference 1", "WHOreference 5"]


# --------------------------------------------------------------------------- #
# Fixtures                                                                    #
# --------------------------------------------------------------------------- #

@pytest.fixture(scope="module")
def analyzer() -> MNSAnalyzer:
    a = MNSAnalyzer()
    assert a.is_loaded(), "GYPA / GYPB GenBank files missing — cannot run MNS tests"
    return a


def _patch_at_genomic_index(seq: str, idx: int, new_base: str) -> str:
    return seq[:idx] + new_base + seq[idx + 1:]


def _patch_snp(analyzer: MNSAnalyzer, snp_name: str, new_base: str) -> str:
    """Patch a SNP in the relevant gene's reference and return the new
    full-gene reference sequence."""
    info = MNS_DIAGNOSTIC_SNPS[snp_name]
    gene_seq = analyzer.gene_refs[info["gene"]]["seq"]
    idx = analyzer._position_to_ref_index(info)
    assert idx is not None, f"could not locate {snp_name} in reference"
    return _patch_at_genomic_index(gene_seq, idx, new_base)


# --------------------------------------------------------------------------- #
# Reference sanity                                                            #
# --------------------------------------------------------------------------- #

def test_reference_loads_both_genes(analyzer: MNSAnalyzer):
    assert "GYPA" in analyzer.gene_refs
    assert "GYPB" in analyzer.gene_refs
    assert analyzer.gene_refs["GYPA"]["cds_start_genomic"] == 5116
    assert analyzer.gene_refs["GYPB"]["cds_start_genomic"] == 7574


def test_file_base_at_c59_matches_audit_note(analyzer: MNSAnalyzer):
    """NG_007470.3 carries the M allele 'C' at GYPA c.59."""
    info = MNS_DIAGNOSTIC_SNPS["c.59C>T (GYPA)"]
    idx = analyzer._position_to_ref_index(info)
    assert analyzer.gene_refs["GYPA"]["seq"][idx] == "C"


def test_file_base_at_c143_matches_audit_note(analyzer: MNSAnalyzer):
    """NG_007483.3 carries the s allele 'C' at GYPB c.143 — parity-inverted entry."""
    info = MNS_DIAGNOSTIC_SNPS["c.143C>T (GYPB)"]
    idx = analyzer._position_to_ref_index(info)
    assert analyzer.gene_refs["GYPB"]["seq"][idx] == "C"
    assert info.get("parity_inverted") is True


# --------------------------------------------------------------------------- #
# M/N axis (GYPA)                                                             #
# --------------------------------------------------------------------------- #

def test_pure_gypa_reference_calls_M_axis_only(analyzer: MNSAnalyzer):
    """GYPA reference alone -> MM on the M/N axis (S/s remains indeterminate
    because no GYPB sequence is in the read)."""
    gypa_ref = analyzer.gene_refs["GYPA"]["seq"]
    result = analyzer.analyze([("pure_gypa", gypa_ref)])
    assert result["m_n_call"]["genotype"] == "MM"
    # S/s axis can't resolve from a GYPA-only read.
    assert result["s_s_call"]["genotype"] is None


def test_patched_jc59_to_T_yields_NN(analyzer: MNSAnalyzer):
    seq = _patch_snp(analyzer, "c.59C>T (GYPA)", "T")
    result = analyzer.analyze([("nn", seq)])
    assert result["m_n_call"]["genotype"] == "NN"


def test_patched_c59_to_Y_yields_MN(analyzer: MNSAnalyzer):
    """IUPAC Y = C/T het at c.59 → MN (sorted alphabetically)."""
    seq = _patch_snp(analyzer, "c.59C>T (GYPA)", "Y")
    result = analyzer.analyze([("mn", seq)])
    assert result["m_n_call"]["genotype"] == "MN"


# --------------------------------------------------------------------------- #
# S/s axis (GYPB) — parity-inverted                                           #
# --------------------------------------------------------------------------- #

def test_pure_gypb_reference_calls_S_axis_only(analyzer: MNSAnalyzer):
    """GYPB reference alone -> ss (parity-inverted: ref=C encodes s)."""
    gypb_ref = analyzer.gene_refs["GYPB"]["seq"]
    result = analyzer.analyze([("pure_gypb", gypb_ref)])
    assert result["s_s_call"]["genotype"] == "ss"


def test_patched_c143_to_T_yields_SS(analyzer: MNSAnalyzer):
    """Patched c.143 C->T (alt) -> SS (parity-inverted: alt=T encodes S)."""
    seq = _patch_snp(analyzer, "c.143C>T (GYPB)", "T")
    result = analyzer.analyze([("ss_pos", seq)])
    assert result["s_s_call"]["genotype"] == "SS"


def test_patched_c143_to_Y_yields_Ss(analyzer: MNSAnalyzer):
    seq = _patch_snp(analyzer, "c.143C>T (GYPB)", "Y")
    result = analyzer.analyze([("ss_het", seq)])
    assert result["s_s_call"]["genotype"] == "Ss"


# --------------------------------------------------------------------------- #
# Combined phenotype + phase ambiguity                                        #
# --------------------------------------------------------------------------- #

def test_double_het_yields_two_phase_options(analyzer: MNSAnalyzer):
    """MN (c.59 = Y) + Ss (c.143 = Y) → both haplotype phases listed."""
    gypa_seq = _patch_snp(analyzer, "c.59C>T (GYPA)", "Y")
    gypb_seq = _patch_snp(analyzer, "c.143C>T (GYPB)", "Y")
    result = analyzer.analyze([("gypa_mn", gypa_seq), ("gypb_ss", gypb_seq)])
    assert result["m_n_call"]["genotype"] == "MN"
    assert result["s_s_call"]["genotype"] == "Ss"
    assert len(result["allele_options"]) == 2
    haps = {opt["haplotypes"] for opt in result["allele_options"]}
    assert haps == {"MS/Ns", "Ms/NS"}


# --------------------------------------------------------------------------- #
# Robustness                                                                  #
# --------------------------------------------------------------------------- #

def test_reverse_complement_input_matches_forward(analyzer: MNSAnalyzer):
    gypa_ref = analyzer.gene_refs["GYPA"]["seq"]
    rev = str(Seq(gypa_ref).reverse_complement())
    fwd_call = analyzer.analyze([("fwd", gypa_ref)])["m_n_call"]["genotype"]
    rev_call = analyzer.analyze([("rev", rev)])["m_n_call"]["genotype"]
    assert fwd_call == rev_call == "MM"


def test_multi_read_concordance_yields_high_confidence(analyzer: MNSAnalyzer):
    gypa_ref = analyzer.gene_refs["GYPA"]["seq"]
    result = analyzer.analyze([
        ("a", gypa_ref), ("b", gypa_ref), ("c", gypa_ref),
    ])
    assert result["snp_consensus"]["c.59C>T (GYPA)"]["confidence"] == "HIGH"


# --------------------------------------------------------------------------- #
# Phred Q-score gate at SNP column                                            #
# --------------------------------------------------------------------------- #

def test_phred_gate_blocks_call_when_snp_base_below_q30(analyzer: MNSAnalyzer):
    """Q-score gate on the M/N axis (GYPA)."""
    info = MNS_DIAGNOSTIC_SNPS["c.59C>T (GYPA)"]
    gypa_ref = analyzer.gene_refs["GYPA"]["seq"]
    snp_idx = analyzer._position_to_ref_index(info)
    quality = [40] * len(gypa_ref)
    quality[snp_idx] = 15
    result = analyzer.analyze([("low_q_at_snp", gypa_ref, quality)])
    primary = result["snp_consensus"]["c.59C>T (GYPA)"]
    assert primary["consensus"] == "no_call"
    snp_call = result["per_read_details"][0]["snp_calls"]["c.59C>T (GYPA)"]
    assert snp_call.get("phred_at_snp") == 15
    assert "Phred" in snp_call.get("reason", "")


def test_no_quality_supplied_skips_phred_gate_backward_compat(analyzer: MNSAnalyzer):
    gypa_ref = analyzer.gene_refs["GYPA"]["seq"]
    result = analyzer.analyze([("legacy_2tuple", gypa_ref)])
    assert result["m_n_call"]["genotype"] == "MM"
    snp_call = result["per_read_details"][0]["snp_calls"]["c.59C>T (GYPA)"]
    assert snp_call.get("phred_at_snp") is None


def test_random_low_identity_read_is_uncallable(analyzer: MNSAnalyzer):
    rng = random.Random(2026)
    junk = "".join(rng.choices("ACGT", k=1000))
    result = analyzer.analyze_single_read(junk, read_id="junk")
    assert result["callable"] is False
    assert result["trusted_snps"] == 0


# --------------------------------------------------------------------------- #
# Missing references                                                          #
# --------------------------------------------------------------------------- #

def test_missing_both_references_raises_with_hint(tmp_path: Path):
    a = MNSAnalyzer(gypa_gb_path=str(tmp_path / "absent_a.gb"),
                    gypb_gb_path=str(tmp_path / "absent_b.gb"))
    assert not a.is_loaded()
    with pytest.raises(MNSReferenceMissingError) as exc_info:
        a.analyze([("x", "ACGT" * 50)])
    msg = str(exc_info.value)
    assert "MNS_REFERENCE_README" in msg
    assert "GYPA" in msg
    assert "GYPB" in msg


def test_missing_only_gypb_still_raises(tmp_path: Path):
    a = MNSAnalyzer(gypa_gb_path=str(GYPA_PATH),
                    gypb_gb_path=str(tmp_path / "absent_b.gb"))
    assert "GYPA" in a.gene_refs and "GYPB" not in a.gene_refs
    assert not a.is_loaded()
    with pytest.raises(MNSReferenceMissingError) as exc_info:
        a.analyze([("x", "ACGT" * 50)])
    assert "GYPB" in str(exc_info.value)


# --------------------------------------------------------------------------- #
# Real patient MIA234 amplicons (smoke only — region coverage unknown)        #
# --------------------------------------------------------------------------- #

def _patient_mns_fastas(sample: str) -> list:
    out = []
    d = PATIENT_ROOT / sample / "MIA234"
    if not d.is_dir():
        return out
    for fa in sorted(d.glob("*.fasta")):
        try:
            rec = next(SeqIO.parse(str(fa), "fasta"))
        except StopIteration:
            continue
        out.append((f"{sample}/MIA234/{fa.stem}", str(rec.seq)))
    return out


@pytest.mark.skipif(
    not PATIENT_ROOT.exists(),
    reason=f"Real patient data not present at {PATIENT_ROOT}",
)
@pytest.mark.parametrize("sample", PATIENT_SAMPLES)
def test_real_patient_analyze_smoke(analyzer: MNSAnalyzer, sample: str):
    """MIA234 is a Miltenberger-region amplicon; it may not cover c.59 or
    c.143. The smoke test just exercises the pipeline shape — it doesn't
    assert phenotype, only that the result dict is well-formed and
    `reads_total` matches input."""
    reads = _patient_mns_fastas(sample)
    if not reads:
        pytest.skip(f"No MIA234 FASTAs found for sample {sample}")

    result = analyzer.analyze(reads)

    assert set(result) >= {
        "phenotype",
        "allele_options",
        "m_n_call",
        "s_s_call",
        "snp_consensus",
        "overall_confidence",
        "reads_total",
        "reads_callable",
        "reads_with_trusted_snps",
        "per_read_details",
    }
    assert result["reads_total"] == len(reads)
    assert isinstance(result["phenotype"], str) and result["phenotype"]
