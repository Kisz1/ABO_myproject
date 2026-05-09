#!/usr/bin/env python
"""Test RHD analysis functionality"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from utils.rhd_analyzer import RHDAnalyzer

def test_rhd_analyzer():
    """Test RHD analyzer"""
    print("=" * 70)
    print("Testing RHD Analyzer")
    print("=" * 70)

    # Test with a sample sequence
    test_seq = (
        "ACTTCACCCTAAGGCTGGATCAGGATCCCCTCCAGGTTTTTACTAGAGCCAAACCCACATCTCCTTTCTCTTCTGCCACC"
        "CCCCCTTAAAATGCTTAGAAACACATAGATTTAAATACAAGTTCAAATGTAAGTAATTTCAACTGTGTAACTATGAGGAG"
        "TCAATTCTACGTGGGTCCTATCTGTATCCTCCCCAGGGCTCAGCTCCATTCTTTGCTTTCATTCATTCTCATTCAATACA"
        "TTGTTGTTAAGAGCTCACTGGGTGCCCTCTCTGTCATGTAGTAAGGTTTTAAAAAGAAAGCCTCTTCTGAGCTTCAGTTT"
        "CCTTATTTATAAAATAGGAATATTGATCTGTTCCTTGCTTTTCTTACAAGGATATGCTGAAGATGACTGAAGTACAGAGT"
    )

    print(f"\n[INFO] Test sequence length: {len(test_seq)} bp")
    print(f"[INFO] Testing with RHD reference data...")

    try:
        analyzer = RHDAnalyzer()
        print("[OK] RHDAnalyzer initialized")

        result = analyzer.analyze(test_seq)
        print("[OK] Analysis completed")

        print("\n" + "=" * 70)
        print("Analyzer Output:")
        print("=" * 70)

        for key, value in result.items():
            if key == 'variants':
                print(f"\n{key}: {len(value)} variants found")
                for v in value[:5]:
                    print(f"  - {v}")
                if len(value) > 5:
                    print(f"  ... and {len(value) - 5} more")
            else:
                print(f"{key}: {value}")

        # Check which fields are present
        print("\n" + "=" * 70)
        print("Field Analysis:")
        print("=" * 70)

        expected_fields = [
            'rhd_status',    # For RhD+ / RhD- determination
            'region',        # RHD1 or RHD456
            'reason',        # Why this decision was made
            'query_length',  # Length of input sequence
        ]

        print("\nExpected fields by main.py:")
        for field in expected_fields:
            if field in result:
                print(f"  [OK] {field}: present")
            else:
                print(f"  [MISSING] {field}: NOT PRESENT")

        print("\nActual fields returned:")
        for field in result.keys():
            print(f"  - {field}")

        return False if any(f not in result for f in expected_fields) else True

    except Exception as e:
        print(f"[FAIL] Error: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_rhd_analyzer()
    if not success:
        print("\n[WARNING] RHD Analyzer is missing required fields!")
        print("[INFO] The analyzer returns: variants, identity, score, strand")
        print("[INFO] But main.py expects: rhd_status, region, reason, query_length")
        print("[INFO] This needs to be fixed for RHD analysis to work properly.")
    sys.exit(0 if success else 1)
