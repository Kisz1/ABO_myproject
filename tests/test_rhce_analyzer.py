"""Unit + smoke tests for RHCEAnalyzer against the NM_020485.7 cDNA reference.

Coverage:
  - Reference loads and the NM_020485.7 ATG offset is detected.
  - The file's actual base at each diagnostic SNP matches the audit comment in
    rhce_analyzer.py (this is a regression sentinel — if NCBI ever revises
    NM_020485.7, the parity-inverted SNP table needs reviewing).
  - Pure c/c, e/e (clean c-allele synthesis) -> phenotype 'ccee' with all
    confirmatory markers supporting.
  - Pure C/C, E/E (clean C-allele synthesis) -> phenotype 'CCEE' with all
    confirmatories supporting and no partial-E flag.
  - Asian partial-E Category EII detected via c.1025 inverted-parity logic.
  - Asian partial-E Category EIII detected via c.1226 (non-inverted) logic.
  - Double-heterozygote CcEe via IUPAC codes -> two haplotype phase options.
  - IUPAC het at primary c.48 -> Cc genotype.
  - Low-identity (random) read -> marked uncallable, no SNP calls trusted.
  - Multi-read concordance -> HIGH confidence at the primary SNPs.
  - Reverse-complement input -> same phenotype as forward input.
  - Real patient FASTAs (if present on this machine) -> analyzer returns a
    well-formed result dict without crashing. Skipped if the data path is
    absent so the test file is portable.
"""

import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from Bio import SeqIO  # noqa: E402
from Bio.Seq import Seq  # noqa: E402

from utils.rhce_analyzer import (  # noqa: E402
    RHCEAnalyzer,
    RHCE_DIAGNOSTIC_SNPS,
)


# --------------------------------------------------------------------------- #
# Fixtures                                                                    #
# --------------------------------------------------------------------------- #

REPO_ROOT = Path(__file__).resolve().parent.parent
FASTA_PATH = REPO_ROOT / "utils" / "data" / "rhce_cdna_reference.fasta"

# Optional real-patient sample data (Sanger FASTAs split by exon amplicon).
# Tests that depend on it skip when the directory is absent.
PATIENT_ROOT = Path(r"C:\Users\ExPertComputer\Desktop\blood_group")
PATIENT_SAMPLES = ["Nampeung", "WHOreference 1", "WHOreference 5"]


@pytest.fixture(scope="module")
def analyzer() -> RHCEAnalyzer:
    a = RHCEAnalyzer()
    assert a.is_loaded(), "NM_020485.7 FASTA missing — cannot run RHCE tests"
    assert a.reference_kind == "cdna"
    return a


@pytest.fixture(scope="module")
def reference_seq(analyzer: RHCEAnalyzer) -> str:
    return analyzer.reference_seq


def _patch(seq: str, cdna_pos: int, new_base: str, offset: int) -> str:
    """Replace the base at cDNA position `cdna_pos` (1-based) with `new_base`."""
    idx = (cdna_pos - 1) + offset
    return seq[:idx] + new_base + seq[idx + 1 :]


def _patch_many(seq: str, mutations: dict, offset: int) -> str:
    out = seq
    for pos, base in mutations.items():
        out = _patch(out, pos, base, offset)
    return out


# --------------------------------------------------------------------------- #
# Reference sanity                                                            #
# --------------------------------------------------------------------------- #

def test_reference_loads_as_nm020485_cdna(analyzer: RHCEAnalyzer):
    assert analyzer.reference_kind == "cdna"
    # ATG context for RHCE: M-S-S-K-Y starts at the detected offset.
    ctx = analyzer.reference_seq[analyzer.cdna_offset : analyzer.cdna_offset + 15]
    assert ctx == "ATGAGCTCTAAGTAC", f"unexpected ATG context: {ctx}"


def test_file_bases_match_audit_note(analyzer: RHCEAnalyzer, reference_seq: str):
    """NM_020485.7 carries these specific bases at the diagnostic SNPs.

    The parity-inverted SNP table (c.178, c.203, c.1025) was built around these
    bases. If NCBI ever revises NM_020485.7 this test will fail and force a
    review of the inverted-parity setup.
    """
    expected = {
        "c.48G>C":   "G",
        "c.178A>C":  "C",  # file carries C-antigen base (inverted in table)
        "c.203G>A":  "A",  # file carries C-antigen base (inverted)
        "c.307C>T":  "C",
        "c.676G>C":  "G",
        "c.1025T>C": "C",  # file carries Cat-EII variant base (inverted)
        "c.1226A>G": "A",
    }
    for snp_name, expected_base in expected.items():
        info = RHCE_DIAGNOSTIC_SNPS[snp_name]
        idx = (info["cDNA_position"] - 1) + analyzer.cdna_offset
        assert reference_seq[idx] == expected_base, (
            f"{snp_name}: file has {reference_seq[idx]}, expected {expected_base}"
        )


# --------------------------------------------------------------------------- #
# Pure homozygous calls                                                       #
# --------------------------------------------------------------------------- #

def test_pure_c_allele_homozygous_gives_ccee(analyzer: RHCEAnalyzer, reference_seq: str):
    # Clean c-allele: revert the three inverted positions so they read as the
    # global wild-type, leaving c.48/c.307/c.676/c.1226 untouched. Also revert
    # c.1025 to T so the partial-E EII flag does not fire.
    clean_c = _patch_many(
        reference_seq,
        {178: "A", 203: "G", 1025: "T"},
        analyzer.cdna_offset,
    )
    result = analyzer.analyze([("clean_c_allele", clean_c)])

    assert result["phenotype"] == "ccee"
    assert result["c_e_call"]["genotype"] == "cc"
    assert result["big_E_call"]["genotype"] == "ee"
    assert result["allele_options"] == [
        {"haplotypes": "ce/ce", "isbt": "RHCE*01/RHCE*01"},
    ]
    # Every confirmatory C/c marker should agree.
    assert result["c_e_call"]["discordant_snps"] == []
    assert "c.178A>C" in result["c_e_call"]["supporting_snps"]
    assert "c.203G>A" in result["c_e_call"]["supporting_snps"]
    assert "c.307C>T" in result["c_e_call"]["supporting_snps"]
    # No partial-E marker fired.
    assert result["big_E_call"]["partial_markers"] == []


def test_pure_capital_C_homozygous_gives_CCEE_no_partial(
    analyzer: RHCEAnalyzer, reference_seq: str
):
    # Clean C-allele: c.48 G->C, c.307 C->T, c.676 G->C produce the C-antigen
    # + E-antigen primary calls. Leave c.178/c.203 at the file's C/A (already
    # the C-allele bases). Revert c.1025 to T so no partial-E flag.
    clean_C = _patch_many(
        reference_seq,
        {48: "C", 307: "T", 676: "C", 1025: "T"},
        analyzer.cdna_offset,
    )
    result = analyzer.analyze([("clean_C_allele", clean_C)])

    assert result["phenotype"] == "CCEE"
    assert result["c_e_call"]["genotype"] == "CC"
    assert result["big_E_call"]["genotype"] == "EE"
    assert result["c_e_call"]["discordant_snps"] == []
    assert result["big_E_call"]["partial_markers"] == []


# --------------------------------------------------------------------------- #
# Asian partial-E variants                                                    #
# --------------------------------------------------------------------------- #

def test_partial_E_category_EII_via_c1025_inverted(
    analyzer: RHCEAnalyzer, reference_seq: str
):
    # Clean c-allele setup BUT keep c.1025 at the file's C (= the partial-EII
    # variant base after parity swap). The inverted-parity trigger is
    # {het, hom_ref}, so hom_ref must fire the partial flag.
    seq = _patch_many(
        reference_seq,
        {178: "A", 203: "G"},        # clean confirmatories
        analyzer.cdna_offset,
    )  # c.1025 left as 'C' on purpose

    result = analyzer.analyze([("eii_carrier", seq)])

    assert result["c_e_call"]["genotype"] == "cc"
    assert result["big_E_call"]["genotype"] == "ee"
    partials = result["big_E_call"]["partial_markers"]
    assert any(p["snp"] == "c.1025T>C" for p in partials), partials
    assert "partial-E variant flagged" in result["phenotype"]
    assert "EII" in result["phenotype"] or "Category EII" in result["phenotype"]


def test_partial_E_category_EIII_via_c1226_variant(
    analyzer: RHCEAnalyzer, reference_seq: str
):
    # Clean c-allele setup + c.1226 A->G (the EIII variant base). c.1226 is NOT
    # parity-inverted, so the trigger set is {het, hom_alt}; hom_alt fires.
    seq = _patch_many(
        reference_seq,
        {178: "A", 203: "G", 1025: "T", 1226: "G"},
        analyzer.cdna_offset,
    )

    result = analyzer.analyze([("eiii_carrier", seq)])

    assert result["c_e_call"]["genotype"] == "cc"
    assert result["big_E_call"]["genotype"] == "ee"
    partials = result["big_E_call"]["partial_markers"]
    assert any(p["snp"] == "c.1226A>G" for p in partials), partials
    assert "partial-E variant flagged" in result["phenotype"]
    assert "EIII" in result["phenotype"] or "Category EIII" in result["phenotype"]


# --------------------------------------------------------------------------- #
# Heterozygous (IUPAC) decoding                                               #
# --------------------------------------------------------------------------- #

def test_iupac_het_at_primary_c48_gives_Cc(analyzer: RHCEAnalyzer, reference_seq: str):
    # S = G/C het at c.48 only; everything else clean c-allele.
    seq = _patch_many(
        reference_seq,
        {48: "S", 178: "A", 203: "G", 1025: "T"},
        analyzer.cdna_offset,
    )
    result = analyzer.analyze([("het_c48", seq)])

    assert result["c_e_call"]["genotype"] == "Cc"
    assert result["big_E_call"]["genotype"] == "ee"
    primary = result["snp_consensus"]["c.48G>C"]
    assert primary["consensus"] == "het"


def test_double_heterozygote_yields_both_phase_options(
    analyzer: RHCEAnalyzer, reference_seq: str
):
    # Het at every diagnostic position. Pick IUPAC codes so that each het is
    # exactly {file_base, the_OTHER_allele_base}:
    #   c.48  file=G, het G/C        -> S
    #   c.178 file=C, het C/A        -> M
    #   c.203 file=A, het A/G        -> R
    #   c.307 file=C, het C/T        -> Y
    #   c.676 file=G, het G/C        -> S
    # Leave c.1025 / c.1226 alone; this test is about phase ambiguity, not
    # partial-E, and the c.1025 partial flag is allowed but not asserted on.
    seq = _patch_many(
        reference_seq,
        {48: "S", 178: "M", 203: "R", 307: "Y", 676: "S", 1025: "T"},
        analyzer.cdna_offset,
    )

    result = analyzer.analyze([("double_het", seq)])

    assert result["c_e_call"]["genotype"] == "Cc"
    assert result["big_E_call"]["genotype"] == "Ee"
    # Sanger cannot phase: both haplotype pairs must be enumerated.
    haps = {opt["haplotypes"] for opt in result["allele_options"]}
    assert {"CE/ce", "Ce/cE"} <= haps, haps
    isbt = {opt["isbt"] for opt in result["allele_options"]}
    assert {"RHCE*04/RHCE*01", "RHCE*02/RHCE*03"} <= isbt, isbt


# --------------------------------------------------------------------------- #
# Robustness                                                                  #
# --------------------------------------------------------------------------- #

def test_random_low_identity_read_is_uncallable(analyzer: RHCEAnalyzer):
    rng = random.Random(1234)
    junk = "".join(rng.choices("ACGT", k=1000))
    result = analyzer.analyze_single_read(junk, read_id="junk")
    assert result["callable"] is False
    # With the local-window gate, individual SNP positions may be "covered"
    # (the local alignment reached them) but they should ALL come back as
    # no_call because random sequence cannot produce an 85%-identity window.
    for call in result["snp_calls"].values():
        if call.get("covered"):
            assert call.get("call") == "no_call", call
        else:
            assert call.get("covered") is False
    # And of course no SNP cleared the gate.
    assert result["trusted_snps"] == 0


def test_reverse_complement_input_matches_forward(
    analyzer: RHCEAnalyzer, reference_seq: str
):
    seq_fwd = _patch_many(
        reference_seq,
        {178: "A", 203: "G", 1025: "T"},
        analyzer.cdna_offset,
    )
    seq_rev = str(Seq(seq_fwd).reverse_complement())

    fwd = analyzer.analyze([("fwd", seq_fwd)])
    rev = analyzer.analyze([("rev", seq_rev)])

    assert fwd["phenotype"] == rev["phenotype"] == "ccee"
    # Reverse-complemented read must be detected on the reverse strand.
    assert rev["per_read_details"][0]["strand"] == "reverse"


def test_multi_read_concordance_yields_high_confidence(
    analyzer: RHCEAnalyzer, reference_seq: str
):
    seq = _patch_many(
        reference_seq,
        {178: "A", 203: "G", 1025: "T"},
        analyzer.cdna_offset,
    )
    result = analyzer.analyze([
        ("read_a", seq),
        ("read_b", seq),
        ("read_c", seq),
    ])
    assert result["snp_consensus"]["c.48G>C"]["confidence"] == "HIGH"
    assert result["snp_consensus"]["c.676G>C"]["confidence"] == "HIGH"
    assert result["reads_callable"] == 3


# --------------------------------------------------------------------------- #
# Phred Q-score gate at SNP column                                            #
# --------------------------------------------------------------------------- #

def test_phred_gate_blocks_call_when_snp_base_below_q30(
    analyzer: RHCEAnalyzer, reference_seq: str
):
    """Q-score gate on the c.48 SNP (G/C, primary C/c discriminator)."""
    info = RHCE_DIAGNOSTIC_SNPS["c.48G>C"]
    snp_idx = (info["cDNA_position"] - 1) + analyzer.cdna_offset
    quality = [40] * len(reference_seq)
    quality[snp_idx] = 15
    result = analyzer.analyze([("low_q_at_snp", reference_seq, quality)])
    primary = result["snp_consensus"]["c.48G>C"]
    assert primary["consensus"] == "no_call"
    snp_call = result["per_read_details"][0]["snp_calls"]["c.48G>C"]
    assert snp_call.get("phred_at_snp") == 15
    assert "Phred" in snp_call.get("reason", "")


def test_no_quality_supplied_skips_phred_gate_backward_compat(
    analyzer: RHCEAnalyzer, reference_seq: str
):
    """Legacy 2-tuple (no quality) must still work."""
    result = analyzer.analyze([("legacy_2tuple", reference_seq)])
    # Ensure the pipeline ran and at least one SNP was called.
    assert isinstance(result.get("phenotype"), str) and result["phenotype"]
    snp_call = result["per_read_details"][0]["snp_calls"]["c.48G>C"]
    # Phred field should not be set when quality wasn't passed.
    assert snp_call.get("phred_at_snp") is None


# --------------------------------------------------------------------------- #
# Real patient FASTAs (smoke test, skipped if data missing)                   #
# --------------------------------------------------------------------------- #

def _patient_fastas(sample: str) -> list:
    """Return list of (read_id, sequence) tuples for RHCE1/2/45 FASTAs."""
    out = []
    for exon_dir in ("RHCE1", "RHCE2", "RHCE45"):
        d = PATIENT_ROOT / sample / exon_dir
        if not d.is_dir():
            continue
        for fa in sorted(d.glob("*.fasta")):
            try:
                rec = next(SeqIO.parse(str(fa), "fasta"))
            except StopIteration:
                continue
            out.append((f"{sample}/{exon_dir}/{fa.stem}", str(rec.seq)))
    return out


@pytest.mark.skipif(
    not PATIENT_ROOT.exists(),
    reason=f"Real patient data not present at {PATIENT_ROOT}",
)
@pytest.mark.parametrize("sample", PATIENT_SAMPLES)
def test_real_patient_analyze_smoke(analyzer: RHCEAnalyzer, sample: str):
    reads = _patient_fastas(sample)
    if not reads:
        pytest.skip(f"No RHCE FASTAs found for sample {sample}")

    result = analyzer.analyze(reads)

    # Don't assert a specific phenotype — we don't have validated ground truth
    # for every patient. Just verify the analyzer produces a well-formed result.
    assert set(result) >= {
        "phenotype",
        "allele_options",
        "c_e_call",
        "big_E_call",
        "snp_consensus",
        "overall_confidence",
        "reads_total",
        "reads_callable",
        "per_read_details",
    }
    assert result["reads_total"] == len(reads)
    assert isinstance(result["phenotype"], str) and result["phenotype"]
    # Patient amplicons are genomic Sanger reads; whole-read identity often
    # drops below 85% because of introns. The per-SNP local-window gate is
    # what should rescue the call, so we assert on that rather than on the
    # legacy `reads_callable` count.
    assert result["reads_with_trusted_snps"] >= 1, (
        f"No SNP cleared the local-window gate for {sample}: "
        f"{[(r['read_id'], r['identity'], r.get('trusted_snps')) for r in result['per_read_details']]}"
    )
    # Regression check for the Option-2 local-window fix: a real RHCE1+2+45
    # amplicon set must resolve to a non-Indeterminate phenotype, not the
    # 'Indeterminate / NONE' we used to get when the whole-read 85% gate
    # threw out the RHCE45 read (which carries the primary E/e marker).
    assert result["phenotype"] != "Indeterminate", (
        f"{sample} still indeterminate after local-window fix: "
        f"reason={result['reason']}; "
        f"per-SNP={ {k: v['consensus'] for k, v in result['snp_consensus'].items()} }"
    )
    assert result["overall_confidence"] != "NONE", (
        f"{sample} overall_confidence still NONE: {result['reason']}"
    )
