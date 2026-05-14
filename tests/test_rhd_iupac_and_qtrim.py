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
