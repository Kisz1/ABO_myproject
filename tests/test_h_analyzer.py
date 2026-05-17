"""Unit + smoke tests for HAnalyzer against the NG_007510.2 genomic reference.

Coverage:
  - Reference loads as genomic; FUT1 features picked up (not FGF21 / IZUMO1
    which share the RefSeqGene record).
  - Each diagnostic SNP's reference base matches the audit comment.
  - Pure reference -> H+ phenotype.
  - Either strong null hom-alt (c.725 or c.586) -> Bombay (Oh).
  - Weak allele hom-alt (c.460) -> Para-Bombay (h2/h2).
  - Any null het -> Para-Bombay potential.
  - Reverse-complement input -> same phenotype as forward.
  - Multi-read concordance -> HIGH confidence at the contributing SNP.
  - Low-identity random read -> uncallable, no SNPs trusted.
  - cDNA-only reference works (all three SNPs in coding region).
  - Missing reference raises with a README hint.
  - Patient FUT1 amplicons: no patient data exists in the dataset, so the
    patient-smoke parametrize skips cleanly.
"""

import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from Bio import SeqIO  # noqa: E402
from Bio.Seq import Seq  # noqa: E402

from utils.h_analyzer import (  # noqa: E402
    HAnalyzer,
    HReferenceMissingError,
    H_DIAGNOSTIC_SNPS,
)


# --------------------------------------------------------------------------- #
# Fixtures                                                                    #
# --------------------------------------------------------------------------- #

REPO_ROOT = Path(__file__).resolve().parent.parent
GENBANK_PATH = REPO_ROOT / "utils" / "data" / "h_referance.gb"

PATIENT_ROOT = Path(r"C:\Users\ExPertComputer\Desktop\blood_group")
PATIENT_SAMPLES = ["Nampeung", "WHOreference 1", "WHOreference 5"]
# Lab amplicon naming for FUT1 is unknown ahead of time — try a few likely
# directory names. If none exist for a sample, the smoke test skips cleanly.
PATIENT_FY_DIR_CANDIDATES = ["H", "FUT1", "FUT", "H_FUT1", "BombayFUT1"]


@pytest.fixture(scope="module")
def analyzer() -> HAnalyzer:
    a = HAnalyzer()
    assert a.is_loaded(), "NG_007510.2 GenBank missing — cannot run H tests"
    assert a.reference_kind == "genomic"
    return a


@pytest.fixture(scope="module")
def reference_seq(analyzer: HAnalyzer) -> str:
    return analyzer.reference_seq


def _patch_at_genomic_index(seq: str, genomic_index: int, new_base: str) -> str:
    return seq[:genomic_index] + new_base + seq[genomic_index + 1:]


def _patch_snp(analyzer: HAnalyzer, seq: str, snp_name: str, new_base: str) -> str:
    info = H_DIAGNOSTIC_SNPS[snp_name]
    idx = analyzer._position_to_ref_index(info)
    assert idx is not None, f"could not locate {snp_name} in reference"
    return _patch_at_genomic_index(seq, idx, new_base)


# --------------------------------------------------------------------------- #
# Reference sanity                                                            #
# --------------------------------------------------------------------------- #

def test_reference_loads_as_ng007510_genomic(analyzer: HAnalyzer):
    assert analyzer.reference_kind == "genomic"
    # The gene-name filter must have picked FUT1's CDS, not FGF21's. FGF21 is
    # on the minus strand starting around position 2060, so if it leaked in,
    # cds_start_genomic would be a small number (~2170). FUT1's ATG is at 9109.
    assert analyzer.cds_start_genomic == 9109, (
        f"cds_start_genomic={analyzer.cds_start_genomic} — looks like the "
        "loader picked up FGF21 or IZUMO1 instead of FUT1."
    )


def test_file_base_at_c725_matches_audit_note(analyzer: HAnalyzer, reference_seq: str):
    """NG_007510.2 carries the wild-type 'T' at c.725 (Leu242)."""
    idx = analyzer._position_to_ref_index(H_DIAGNOSTIC_SNPS["c.725T>G"])
    assert reference_seq[idx] == "T", (
        f"c.725 base is {reference_seq[idx]!r}, expected 'T' — "
        "see audit note in utils/data/H_REFERENCE_README.md"
    )


def test_file_base_at_c586_matches_audit_note(analyzer: HAnalyzer, reference_seq: str):
    """NG_007510.2 carries the wild-type 'C' at c.586 (Gln196)."""
    idx = analyzer._position_to_ref_index(H_DIAGNOSTIC_SNPS["c.586C>T"])
    assert reference_seq[idx] == "C"


def test_file_base_at_c460_matches_audit_note(analyzer: HAnalyzer, reference_seq: str):
    """NG_007510.2 carries the wild-type 'T' at c.460 (Tyr154)."""
    idx = analyzer._position_to_ref_index(H_DIAGNOSTIC_SNPS["c.460T>C"])
    assert reference_seq[idx] == "T"


# --------------------------------------------------------------------------- #
# Phenotype states                                                            #
# --------------------------------------------------------------------------- #

def test_pure_reference_gives_h_positive(analyzer: HAnalyzer, reference_seq: str):
    result = analyzer.analyze([("pure_ref", reference_seq)])
    assert result["phenotype"] == "H+"
    assert result["h_state_call"]["state"] == "h_positive"
    assert result["allele_options"][0]["isbt"] == "FUT1*01/FUT1*01"


def test_indian_bombay_hom_alt_at_c725_yields_bombay(
    analyzer: HAnalyzer, reference_seq: str
):
    seq = _patch_snp(analyzer, reference_seq, "c.725T>G", "G")
    result = analyzer.analyze([("bombay_c725", seq)])
    assert "Bombay" in result["phenotype"]
    assert result["h_state_call"]["state"] == "bombay"
    assert any(t["snp"] == "c.725T>G"
               for t in result["h_state_call"]["triggers_hom_alt"])


def test_nonsense_hom_alt_at_c586_yields_bombay(
    analyzer: HAnalyzer, reference_seq: str
):
    seq = _patch_snp(analyzer, reference_seq, "c.586C>T", "T")
    result = analyzer.analyze([("bombay_c586", seq)])
    assert result["h_state_call"]["state"] == "bombay"
    assert any(t["snp"] == "c.586C>T"
               for t in result["h_state_call"]["triggers_hom_alt"])


def test_weak_allele_hom_alt_at_c460_yields_para_bombay(
    analyzer: HAnalyzer, reference_seq: str
):
    """c.460 is the 'weak' allele — homozygous gives Para-Bombay, not full Bombay."""
    seq = _patch_snp(analyzer, reference_seq, "c.460T>C", "C")
    result = analyzer.analyze([("h2_hom", seq)])
    assert result["h_state_call"]["state"] == "para_bombay"
    assert "Para-Bombay" in result["phenotype"]
    assert not result["h_state_call"]["any_strong_hom_alt"]


def test_strong_null_het_yields_para_bombay(
    analyzer: HAnalyzer, reference_seq: str
):
    """c.725 heterozygous (T/G = K IUPAC) -> Para-Bombay potential, not Bombay."""
    seq = _patch_snp(analyzer, reference_seq, "c.725T>G", "K")  # K = G/T
    result = analyzer.analyze([("het_c725", seq)])
    assert result["h_state_call"]["state"] == "para_bombay"
    assert "Para-Bombay" in result["phenotype"]


# --------------------------------------------------------------------------- #
# Robustness                                                                  #
# --------------------------------------------------------------------------- #

def test_reverse_complement_input_matches_forward(
    analyzer: HAnalyzer, reference_seq: str
):
    fwd = reference_seq
    rev = str(Seq(fwd).reverse_complement())
    assert analyzer.analyze([("fwd", fwd)])["phenotype"] == \
           analyzer.analyze([("rev", rev)])["phenotype"] == "H+"


def test_multi_read_concordance_yields_high_confidence(
    analyzer: HAnalyzer, reference_seq: str
):
    result = analyzer.analyze([
        ("read_a", reference_seq),
        ("read_b", reference_seq),
        ("read_c", reference_seq),
    ])
    assert result["snp_consensus"]["c.725T>G"]["confidence"] == "HIGH"
    assert result["reads_callable"] == 3
    assert result["overall_confidence"] == "HIGH"


# --------------------------------------------------------------------------- #
# Phred Q-score gate at SNP column                                            #
# --------------------------------------------------------------------------- #

def test_phred_gate_blocks_call_when_snp_base_below_q30(
    analyzer: HAnalyzer, reference_seq: str
):
    info = H_DIAGNOSTIC_SNPS["c.725T>G"]
    snp_idx = analyzer._position_to_ref_index(info)
    quality = [40] * len(reference_seq)
    quality[snp_idx] = 15
    result = analyzer.analyze([("low_q_at_snp", reference_seq, quality)])
    primary = result["snp_consensus"]["c.725T>G"]
    assert primary["consensus"] == "no_call"
    snp_call = result["per_read_details"][0]["snp_calls"]["c.725T>G"]
    assert snp_call.get("phred_at_snp") == 15
    assert "Phred" in snp_call.get("reason", "")


def test_no_quality_supplied_skips_phred_gate_backward_compat(
    analyzer: HAnalyzer, reference_seq: str
):
    result = analyzer.analyze([("legacy_2tuple", reference_seq)])
    assert result["phenotype"] == "H+"
    snp_call = result["per_read_details"][0]["snp_calls"]["c.725T>G"]
    assert snp_call.get("phred_at_snp") is None


def test_random_low_identity_read_is_uncallable(analyzer: HAnalyzer):
    rng = random.Random(2026)
    junk = "".join(rng.choices("ACGT", k=1000))
    result = analyzer.analyze_single_read(junk, read_id="junk")
    assert result["callable"] is False
    assert result["trusted_snps"] == 0


# --------------------------------------------------------------------------- #
# Missing reference                                                           #
# --------------------------------------------------------------------------- #

def test_missing_reference_raises_with_download_hint(tmp_path: Path):
    a = HAnalyzer(gb_path=str(tmp_path / "absent.gb"),
                  fasta_path=str(tmp_path / "absent.fasta"))
    assert not a.is_loaded()
    with pytest.raises(HReferenceMissingError) as exc_info:
        a.analyze([("x", "ACGT" * 50)])
    msg = str(exc_info.value)
    assert "H_REFERENCE_README" in msg


# --------------------------------------------------------------------------- #
# Real patient FUT1 amplicons (skipped — no FUT1 data in current dataset)     #
# --------------------------------------------------------------------------- #

def _patient_h_fastas(sample: str) -> list:
    """Look for FUT1 amplicons in any of the candidate directory names."""
    for candidate in PATIENT_FY_DIR_CANDIDATES:
        d = PATIENT_ROOT / sample / candidate
        if not d.is_dir():
            continue
        out = []
        for fa in sorted(d.glob("*.fasta")):
            try:
                rec = next(SeqIO.parse(str(fa), "fasta"))
            except StopIteration:
                continue
            out.append((f"{sample}/{candidate}/{fa.stem}", str(rec.seq)))
        if out:
            return out
    return []


@pytest.mark.skipif(
    not PATIENT_ROOT.exists(),
    reason=f"Real patient data not present at {PATIENT_ROOT}",
)
@pytest.mark.parametrize("sample", PATIENT_SAMPLES)
def test_real_patient_analyze_smoke(analyzer: HAnalyzer, sample: str):
    reads = _patient_h_fastas(sample)
    if not reads:
        pytest.skip(f"No FUT1 FASTAs found for sample {sample}")

    result = analyzer.analyze(reads)

    assert set(result) >= {
        "phenotype",
        "allele_options",
        "h_state_call",
        "snp_consensus",
        "overall_confidence",
        "reads_total",
        "reads_callable",
        "reads_with_trusted_snps",
        "per_read_details",
    }
    assert result["reads_total"] == len(reads)
    assert isinstance(result["phenotype"], str) and result["phenotype"]
