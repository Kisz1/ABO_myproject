"""Unit tests for RHD path improvements:

  1. IUPAC-aware diagnostic SNP detection (heterozygous + N handling)
  2. Phred sliding-window end-trim + internal N-masking

These guard the clinical contract that:
  - A heterozygous c.1227 G/A (IUPAC 'R') is detected as Weak D type 4 / DEL
  - A homozygous c.1227 A is still detected as before
  - A no-call ('N') at the diagnostic position is suppressed
  - Low-quality 5'/3' ends are trimmed; remaining low-Q bases become 'N'
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from utils.rhd_analyzer import RHDAnalyzer  # noqa: E402
from main import quality_trim_and_mask  # noqa: E402


# --------------------------------------------------------------------------- #
# IUPAC-aware diagnostic SNP detection                                        #
# --------------------------------------------------------------------------- #

def _make_ref_and_query(pos_1based, ref_base, query_base, length=1500):
    ref = ['A'] * length
    ref[pos_1based - 1] = ref_base
    query = list(ref)
    query[pos_1based - 1] = query_base
    return ''.join(ref), ''.join(query)


def test_homozygous_alt_is_detected_with_zygosity_hom():
    analyzer = RHDAnalyzer()
    ref, query = _make_ref_and_query(1227, 'G', 'A')
    detected = analyzer._detect_diagnostic_snps(query, ref)
    assert 'c.1227G>A' in detected
    assert detected['c.1227G>A']['zygosity'] == 'hom'
    assert detected['c.1227G>A']['query_base'] == 'A'


def test_heterozygous_iupac_R_is_detected_with_zygosity_het():
    analyzer = RHDAnalyzer()
    ref, query = _make_ref_and_query(1227, 'G', 'R')  # R = A/G
    detected = analyzer._detect_diagnostic_snps(query, ref)
    assert 'c.1227G>A' in detected
    assert detected['c.1227G>A']['zygosity'] == 'het'
    assert detected['c.1227G>A']['query_base'] == 'R'


def test_N_at_diagnostic_position_is_no_call():
    analyzer = RHDAnalyzer()
    ref, query = _make_ref_and_query(1227, 'G', 'N')
    detected = analyzer._detect_diagnostic_snps(query, ref)
    assert 'c.1227G>A' not in detected, "N must be treated as no-call"


def test_unrelated_iupac_does_not_falsely_call():
    """Y = C/T contains neither ref G nor alt A → must not trigger c.1227G>A."""
    analyzer = RHDAnalyzer()
    ref, query = _make_ref_and_query(1227, 'G', 'Y')
    detected = analyzer._detect_diagnostic_snps(query, ref)
    assert 'c.1227G>A' not in detected


def test_reference_base_unchanged_is_not_called():
    analyzer = RHDAnalyzer()
    ref, query = _make_ref_and_query(1227, 'G', 'G')
    detected = analyzer._detect_diagnostic_snps(query, ref)
    assert 'c.1227G>A' not in detected


def test_phenotype_reason_includes_heterozygous_tag():
    """Heterozygous c.1227 G/A on a high-identity query → Weak D type 4 with '(heterozygous)' in reason."""
    analyzer = RHDAnalyzer()
    detected = {
        'c.1227G>A': {
            'exon': 9,
            'cDNA_position': 1227,
            'reference_base': 'G',
            'query_base': 'R',
            'zygosity': 'het',
            'alleles': ['RHD*01W.4'],
            'significance': 'test',
            'reference': 'test',
        }
    }
    decision = analyzer.determine_rhd_phenotype_snp_based(
        identity=98.0, diagnostic_snps=detected, query_seq='', variants=[])
    assert 'Weak D type 4' in decision['phenotype']
    assert decision['allele'] == 'RHD*01W.4'
    assert '(heterozygous)' in decision['reason']
    assert decision.get('zygosity') == 'het'


def test_phenotype_reason_includes_homozygous_tag():
    analyzer = RHDAnalyzer()
    detected = {
        'c.1227G>A': {
            'exon': 9,
            'cDNA_position': 1227,
            'reference_base': 'G',
            'query_base': 'A',
            'zygosity': 'hom',
            'alleles': ['RHD*01W.4'],
            'significance': 'test',
            'reference': 'test',
        }
    }
    decision = analyzer.determine_rhd_phenotype_snp_based(
        identity=98.0, diagnostic_snps=detected, query_seq='', variants=[])
    assert '(homozygous)' in decision['reason']
    assert decision.get('zygosity') == 'hom'


# --------------------------------------------------------------------------- #
# Quality trimming + N-masking                                                #
# --------------------------------------------------------------------------- #

def test_qtrim_high_quality_no_change():
    seq = 'ACGT' * 25
    qual = [40] * len(seq)
    left, right, masked, n_masked = quality_trim_and_mask(seq, qual, q_threshold=20, window=10)
    assert left == 0
    assert right == len(seq)
    assert masked == seq
    assert n_masked == 0


def test_qtrim_strips_low_quality_5p_end():
    seq = 'ACGT' * 25  # 100 bp
    qual = [5] * 12 + [40] * 88
    left, right, masked, n_masked = quality_trim_and_mask(seq, qual, q_threshold=20, window=10)
    assert left == 12
    assert right == 100
    assert n_masked == 0
    assert masked == seq[12:]


def test_qtrim_strips_low_quality_3p_end():
    seq = 'ACGT' * 25
    qual = [40] * 88 + [5] * 12
    left, right, masked, n_masked = quality_trim_and_mask(seq, qual, q_threshold=20, window=10)
    assert left == 0
    assert right == 88
    assert n_masked == 0


def test_qtrim_masks_internal_low_quality_as_N():
    seq = 'ACGTACGTACGT' * 10  # 120 bp
    qual = [40] * 120
    qual[60] = 5
    qual[61] = 5
    left, right, masked, n_masked = quality_trim_and_mask(seq, qual, q_threshold=20, window=10)
    assert left == 0
    assert right == 120
    assert n_masked == 2
    assert masked[60] == 'N'
    assert masked[61] == 'N'
    # Surrounding bases unchanged
    assert masked[59] == seq[59]
    assert masked[62] == seq[62]


def test_qtrim_all_low_quality_returns_empty():
    seq = 'ACGT' * 25
    qual = [5] * len(seq)
    left, right, masked, n_masked = quality_trim_and_mask(seq, qual, q_threshold=20, window=10)
    assert right <= left
    assert masked == ''


def test_qtrim_missing_quality_falls_back_to_permissive():
    """When quality length doesn't match seq length, treat as fully passing."""
    seq = 'ACGT' * 25
    left, right, masked, n_masked = quality_trim_and_mask(seq, [], q_threshold=20, window=10)
    assert left == 0
    assert right == len(seq)
    assert masked == seq
    assert n_masked == 0


# --------------------------------------------------------------------------- #
# End-to-end: heterozygous IUPAC propagates through analyze() and voting      #
# --------------------------------------------------------------------------- #

def test_heterozygous_iupac_flows_through_analyze():
    """A query with R at c.1227 should be reported with zygosity='het' in the
    full analyze() output (when identity is high enough to pass priority 4)."""
    analyzer = RHDAnalyzer()
    # Build a query that is mostly RHD456 reference but with R at cDNA pos 1227
    from utils.rhd_analyzer import RHD456_REFERENCE
    query = list(RHD456_REFERENCE)
    if 1226 < len(query) and query[1226] == 'G':
        query[1226] = 'R'
        query_seq = ''.join(query)
        result = analyzer.analyze(query_seq)
        # If the diagnostic SNP at 1227 fires, zygosity should be 'het'
        snps = result.get('diagnostic_snps', {})
        if 'c.1227G>A' in snps:
            assert snps['c.1227G>A']['zygosity'] == 'het'
            assert result.get('zygosity') == 'het'


# --------------------------------------------------------------------------- #
# Phred Q-score gate at SNP column (Q30, lab standard)                        #
# --------------------------------------------------------------------------- #

def test_phred_gate_blocks_snp_when_base_below_q30():
    """A SNP at the c.1227 position with the alt base 'A' but Phred Q15 must
    be silently suppressed — the analyzer does not include it in detected_snps
    because the base call is not reliable enough to act on clinically.
    """
    analyzer = RHDAnalyzer()
    ref, query = _make_ref_and_query(1227, 'G', 'A')
    # Build a quality list: all Q40 except c.1227 which is Q15.
    quality = [40] * len(query)
    quality[1227 - 1] = 15

    detected = analyzer._detect_diagnostic_snps(query, ref, query_quality=quality)
    # Below-Q30 Phred at the SNP position -> SNP not detected.
    assert 'c.1227G>A' not in detected, (
        f"Expected Q15 at c.1227 to suppress the SNP, got {detected}"
    )


def test_phred_gate_passes_snp_when_base_at_q30_or_above():
    """Q30 at the SNP column should pass the gate and detect the SNP normally,
    with the Phred score attached to the entry for audit."""
    analyzer = RHDAnalyzer()
    ref, query = _make_ref_and_query(1227, 'G', 'A')
    quality = [40] * len(query)
    quality[1227 - 1] = 30  # exactly at threshold

    detected = analyzer._detect_diagnostic_snps(query, ref, query_quality=quality)
    assert 'c.1227G>A' in detected
    assert detected['c.1227G>A']['zygosity'] == 'hom'
    assert detected['c.1227G>A']['phred_at_snp'] == 30


def test_no_quality_supplied_skips_phred_gate_backward_compat():
    """Legacy 2-arg call (no quality) must still detect the SNP — FASTA
    inputs have no Phred and shouldn't be penalised."""
    analyzer = RHDAnalyzer()
    ref, query = _make_ref_and_query(1227, 'G', 'A')
    detected = analyzer._detect_diagnostic_snps(query, ref)  # no quality
    assert 'c.1227G>A' in detected
    assert 'phred_at_snp' not in detected['c.1227G>A']


def test_analyze_multi_amplicon_accepts_3_tuples_with_quality():
    """The amplicon-voting entry point must accept 3-tuples (name, seq, qual)
    in addition to legacy 2-tuples without crashing — verifies the input
    contract change. The actual gate behaviour is tested at the
    _detect_diagnostic_snps level (the unit tests above), which is the
    correct layer to assert on since `analyze_multi_amplicon` requires real
    RHD-aligned reads to produce non-trivial diagnostic_snps output.
    """
    analyzer = RHDAnalyzer()
    _, query = _make_ref_and_query(1227, 'G', 'A', length=1500)
    quality = [40] * len(query)
    quality[1227 - 1] = 10  # low Phred at SNP

    # Both 2-tuple (legacy) and 3-tuple (new) inputs must produce well-formed
    # result dicts with the same top-level keys.
    result_2tuple = analyzer.analyze_multiple_amplicons([("legacy", query)])
    result_3tuple = analyzer.analyze_multiple_amplicons([("with_qual", query, quality)])

    for r in (result_2tuple, result_3tuple):
        assert set(r) >= {
            'amplicon_results', 'votes', 'final_verdict',
            'confidence', 'details', 'total_amplicons',
        }
        assert r['total_amplicons'] == 1


def test_legacy_2tuple_still_accepted_after_phred_refactor():
    """Sanity: a list of 2-tuples (the old FASTA-style API) must not crash
    after the 3-tuple support was added."""
    analyzer = RHDAnalyzer()
    _, query = _make_ref_and_query(1227, 'G', 'A', length=1500)
    result = analyzer.analyze_multiple_amplicons([
        ("a", query), ("b", query),
    ])
    assert result['total_amplicons'] == 2


# --------------------------------------------------------------------------- #
# Region picking + voting quality floor + RHD456 weighting                    #
# These guard the fix that stopped short exon 4-6 reads from being routed to  #
# the RHD1 reference (which made identity collapse to ~60% and flip them to   #
# RhD-, causing true RhD+ samples to be miscalled).                           #
# --------------------------------------------------------------------------- #

def test_region_picked_by_identity_not_length():
    """A short query that matches RHD456_REFERENCE near-perfectly must be
    routed to RHD456 even though its length (<1200 bp) would have triggered
    the old length-based router to send it to RHD1."""
    from utils.rhd_analyzer import RHD456_REFERENCE
    analyzer = RHDAnalyzer()
    # 600 bp slice from middle of RHD456 reference — clearly an exon 4-6 read.
    short_rhd456 = RHD456_REFERENCE[500:1100]
    result = analyzer.analyze(short_rhd456)
    assert result['region'] == 'RHD456', (
        f"600 bp RHD456 slice should pick RHD456, got {result['region']} "
        f"with identity {result['identity']}"
    )
    assert result['identity'] >= 95.0


def test_low_identity_amplicon_votes_inconclusive_not_rhd_minus():
    """An amplicon whose best identity to either reference is below
    MIN_AMPLICON_IDENTITY (noise floor) must vote Inconclusive, not RhD-.
    Garbage reads must not pollute the vote tally with false RhD- evidence."""
    analyzer = RHDAnalyzer()
    # Random-ish 500 bp sequence — won't align well to either reference.
    junk = ('ACGT' * 125)
    result = analyzer.analyze_multiple_amplicons([('junk', junk)])
    assert result['votes']['Inconclusive'] >= 1, (
        f"Junk read should vote Inconclusive, got {result['votes']}"
    )
    assert result['votes']['RhD-'] == 0


def test_rhd456_vote_outweighs_rhd1_when_they_disagree():
    """A single RHD456 deletion-level read should outweigh a single RHD1
    high-identity read in the final verdict. This encodes the biological
    fact that RHD1 can cross-amplify RHCE (giving false RHD presence),
    while RHD456 covers the discriminating antigenic determinants.
    Asserted via the weighted vote tally (RHD456 weight = 2, RHD1 = 1)."""
    from utils.rhd_analyzer import RHD1_REFERENCE, RHD456_REFERENCE, REGION_VOTE_WEIGHT
    assert REGION_VOTE_WEIGHT['RHD456'] > REGION_VOTE_WEIGHT['RHD1'], (
        "Regression: RHD456 must carry strictly more weight than RHD1"
    )
    analyzer = RHDAnalyzer()
    # A read that perfectly matches RHD1 only — picked region RHD1, vote RhD+.
    rhd1_hit = RHD1_REFERENCE[:800]
    # A read that matches RHD456 at ~82% identity (75-85% range -> RhD-,
    # well above the noise floor) by mutating every 6th base. The dominant
    # alignment is to RHD456, not RHD1.
    mutated = list(RHD456_REFERENCE[500:1500])
    for i in range(0, len(mutated), 6):
        mutated[i] = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}.get(mutated[i], 'A')
    rhd456_deletion_signal = ''.join(mutated)
    result = analyzer.analyze_multiple_amplicons([
        ('rhd1_hit', rhd1_hit),
        ('rhd456_del', rhd456_deletion_signal),
    ])
    # Sanity: confirm the second read was actually routed to RHD456 (so the
    # weighting actually kicks in). If region picking misroutes, the test
    # fails for the wrong reason — surface that here with a clearer message.
    regions = [a['region'] for a in result['amplicon_results']]
    assert 'RHD456' in regions, (
        f"Test setup failure: expected one amplicon routed to RHD456, "
        f"got regions={regions}"
    )
    assert result['final_verdict'].startswith('RhD-'), (
        f"RHD456 RhD- evidence should outweigh single RHD1 RhD+ vote; "
        f"got verdict={result['final_verdict']!r}, weighted_votes={result['votes']}"
    )


if __name__ == '__main__':
    import traceback
    tests = [v for k, v in dict(globals()).items() if k.startswith('test_')]
    failed = []
    for t in tests:
        try:
            t()
            print(f"[OK]   {t.__name__}")
        except AssertionError as e:
            failed.append(t.__name__)
            print(f"[FAIL] {t.__name__}: {e}")
        except Exception:
            failed.append(t.__name__)
            print(f"[ERROR] {t.__name__}:")
            traceback.print_exc()
    print(f"\n{len(tests) - len(failed)}/{len(tests)} passed")
    sys.exit(1 if failed else 0)
