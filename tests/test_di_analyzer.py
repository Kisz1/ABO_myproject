"""Unit + smoke tests for DIAnalyzer against the NG_007498.1 genomic reference.

Coverage:
  - Reference loads as genomic; SLC4A1 features picked up.
  - c.2561 audit: ref = C (Di-b allele) — parity-inverted entry.
  - Pure reference -> Di(a-b+) (parity-inverted: hom_ref = bb).
  - Patched c.2561 C->T -> Di(a+b-) (hom_alt = aa).
  - IUPAC het Y -> Di(a+b+) (genotype 'ab', sorted alphabetically).
  - Reverse-complement input -> same phenotype as forward.
  - Multi-read concordance -> HIGH confidence.
  - Random low-identity read -> uncallable.
  - Missing reference raises with README hint.
  - Real patient DI1819 amplicons -> non-Indeterminate phenotype.
"""

import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from Bio import SeqIO  # noqa: E402
from Bio.Seq import Seq  # noqa: E402

from utils.di_analyzer import (  # noqa: E402
    DIAnalyzer,
    DIReferenceMissingError,
    DI_DIAGNOSTIC_SNPS,
)


REPO_ROOT = Path(__file__).resolve().parent.parent
GENBANK_PATH = REPO_ROOT / "utils" / "data" / "di_referance.gb"

PATIENT_ROOT = Path(r"C:\Users\ExPertComputer\Desktop\blood_group")
PATIENT_SAMPLES = ["Nampeung", "WHOreference 1", "WHOreference 5"]


@pytest.fixture(scope="module")
def analyzer() -> DIAnalyzer:
    a = DIAnalyzer()
    assert a.is_loaded(), "NG_007498.1 GenBank missing — cannot run DI tests"
    assert a.reference_kind == "genomic"
    return a


@pytest.fixture(scope="module")
def reference_seq(analyzer: DIAnalyzer) -> str:
    return analyzer.reference_seq


def _patch_at_genomic_index(seq: str, idx: int, new_base: str) -> str:
    return seq[:idx] + new_base + seq[idx + 1:]


def _patch_snp(analyzer: DIAnalyzer, seq: str, snp_name: str, new_base: str) -> str:
    info = DI_DIAGNOSTIC_SNPS[snp_name]
    idx = analyzer._position_to_ref_index(info)
    assert idx is not None, f"could not locate {snp_name} in reference"
    return _patch_at_genomic_index(seq, idx, new_base)


# --------------------------------------------------------------------------- #
# Reference sanity                                                            #
# --------------------------------------------------------------------------- #

def test_reference_loads_as_ng007498_genomic(analyzer: DIAnalyzer):
    assert analyzer.reference_kind == "genomic"
    assert analyzer.cds_start_genomic == 10268
    assert len(analyzer.exon_map) >= 19, analyzer.exon_map


def test_file_base_at_c2561_matches_audit_note(analyzer: DIAnalyzer, reference_seq: str):
    """NG_007498.1 carries the Di-b base 'C' at c.2561 (parity-inverted entry).

    If a future NCBI revision flips this to 'T' (a Di-a clone), the SNP
    entry must be un-inverted: swap ref_base/alt_base + ref_call/alt_call,
    and remove `parity_inverted: True`.
    """
    info = DI_DIAGNOSTIC_SNPS["c.2561T>C"]
    idx = analyzer._position_to_ref_index(info)
    assert reference_seq[idx] == "C", (
        f"c.2561 base in NG_007498.1 is {reference_seq[idx]!r}, expected 'C' — "
        "see audit note in utils/data/DI_REFERENCE_README.md"
    )
    assert info.get("parity_inverted") is True


# --------------------------------------------------------------------------- #
# Phenotype calls (parity-inverted: ref = Di-b)                               #
# --------------------------------------------------------------------------- #

def test_pure_reference_gives_di_a_neg_b_pos(analyzer: DIAnalyzer, reference_seq: str):
    """Pure ref (Di-b/Di-b) -> Di(a-b+)."""
    result = analyzer.analyze([("pure_ref", reference_seq)])
    assert result["phenotype"] == "Di(a-b+)"
    assert result["di_axis_call"]["genotype"] == "bb"


def test_patched_c2561_to_T_yields_di_a_pos_b_neg(
    analyzer: DIAnalyzer, reference_seq: str
):
    """Flip c.2561 C->T (alt) -> Di(a+b-) (hom_alt aa)."""
    seq = _patch_snp(analyzer, reference_seq, "c.2561T>C", "T")
    result = analyzer.analyze([("pure_dia", seq)])
    assert result["phenotype"] == "Di(a+b-)"
    assert result["di_axis_call"]["genotype"] == "aa"


def test_iupac_het_at_c2561_yields_di_a_pos_b_pos(
    analyzer: DIAnalyzer, reference_seq: str
):
    """IUPAC Y = C/T het -> 'ab' (alphabetical) -> Di(a+b+)."""
    seq = _patch_snp(analyzer, reference_seq, "c.2561T>C", "Y")
    result = analyzer.analyze([("het", seq)])
    assert result["phenotype"] == "Di(a+b+)"
    assert result["di_axis_call"]["genotype"] == "ab"


# --------------------------------------------------------------------------- #
# Robustness                                                                  #
# --------------------------------------------------------------------------- #

def test_reverse_complement_input_matches_forward(
    analyzer: DIAnalyzer, reference_seq: str
):
    fwd = reference_seq
    rev = str(Seq(fwd).reverse_complement())
    assert (analyzer.analyze([("fwd", fwd)])["phenotype"]
            == analyzer.analyze([("rev", rev)])["phenotype"]
            == "Di(a-b+)")


def test_multi_read_concordance_yields_high_confidence(
    analyzer: DIAnalyzer, reference_seq: str
):
    result = analyzer.analyze([
        ("a", reference_seq), ("b", reference_seq), ("c", reference_seq),
    ])
    assert result["snp_consensus"]["c.2561T>C"]["confidence"] == "HIGH"
    assert result["overall_confidence"] == "HIGH"


# --------------------------------------------------------------------------- #
# Phred Q-score gate at SNP column                                            #
# --------------------------------------------------------------------------- #

def test_phred_gate_blocks_call_when_snp_base_below_q30(
    analyzer: DIAnalyzer, reference_seq: str
):
    info = DI_DIAGNOSTIC_SNPS["c.2561T>C"]
    snp_idx = analyzer._position_to_ref_index(info)
    quality = [40] * len(reference_seq)
    quality[snp_idx] = 15
    result = analyzer.analyze([("low_q_at_snp", reference_seq, quality)])
    primary = result["snp_consensus"]["c.2561T>C"]
    assert primary["consensus"] == "no_call"
    snp_call = result["per_read_details"][0]["snp_calls"]["c.2561T>C"]
    assert snp_call.get("phred_at_snp") == 15
    assert "Phred" in snp_call.get("reason", "")


def test_no_quality_supplied_skips_phred_gate_backward_compat(
    analyzer: DIAnalyzer, reference_seq: str
):
    result = analyzer.analyze([("legacy_2tuple", reference_seq)])
    # Parity-inverted Di-b reference -> pure ref calls bb (Di(a-b+)).
    assert result["di_axis_call"]["genotype"] == "bb"
    snp_call = result["per_read_details"][0]["snp_calls"]["c.2561T>C"]
    assert snp_call.get("phred_at_snp") is None


def test_random_low_identity_read_is_uncallable(analyzer: DIAnalyzer):
    rng = random.Random(2026)
    junk = "".join(rng.choices("ACGT", k=1000))
    result = analyzer.analyze_single_read(junk, read_id="junk")
    assert result["callable"] is False
    assert result["trusted_snps"] == 0


# --------------------------------------------------------------------------- #
# Missing reference                                                           #
# --------------------------------------------------------------------------- #

def test_missing_reference_raises_with_download_hint(tmp_path: Path):
    a = DIAnalyzer(gb_path=str(tmp_path / "absent.gb"),
                   fasta_path=str(tmp_path / "absent.fasta"))
    assert not a.is_loaded()
    with pytest.raises(DIReferenceMissingError) as exc_info:
        a.analyze([("x", "ACGT" * 50)])
    assert "DI_REFERENCE_README" in str(exc_info.value)


# --------------------------------------------------------------------------- #
# Real patient DI1819 amplicons                                               #
# --------------------------------------------------------------------------- #

def _patient_di_fastas(sample: str) -> list:
    out = []
    d = PATIENT_ROOT / sample / "DI1819"
    if not d.is_dir():
        return out
    for fa in sorted(d.glob("*.fasta")):
        try:
            rec = next(SeqIO.parse(str(fa), "fasta"))
        except StopIteration:
            continue
        out.append((f"{sample}/DI1819/{fa.stem}", str(rec.seq)))
    return out


@pytest.mark.skipif(
    not PATIENT_ROOT.exists(),
    reason=f"Real patient data not present at {PATIENT_ROOT}",
)
@pytest.mark.parametrize("sample", PATIENT_SAMPLES)
def test_real_patient_analyze_smoke(analyzer: DIAnalyzer, sample: str):
    reads = _patient_di_fastas(sample)
    if not reads:
        pytest.skip(f"No DI1819 FASTAs found for sample {sample}")

    result = analyzer.analyze(reads)

    assert set(result) >= {
        "phenotype",
        "allele_options",
        "di_axis_call",
        "snp_consensus",
        "overall_confidence",
        "reads_total",
        "reads_callable",
        "reads_with_trusted_snps",
        "per_read_details",
    }
    assert result["reads_total"] == len(reads)
    assert isinstance(result["phenotype"], str) and result["phenotype"]
    # DI1819 covers c.2561, so the primary axis should resolve.
    assert result["phenotype"] != "Indeterminate", (
        f"{sample} indeterminate: reason={result['reason']}; "
        f"primary={result['snp_consensus']['c.2561T>C']}"
    )
