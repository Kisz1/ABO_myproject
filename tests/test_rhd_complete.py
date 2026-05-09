#!/usr/bin/env python
"""Test complete RHD decision logic"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from utils.rhd_analyzer import RHDAnalyzer, RHD1_REFERENCE, RHD456_REFERENCE

def test_rhd_decision_logic():
    """Test RHD analyzer with complete decision logic"""
    print("=" * 70)
    print("RHD Analyzer - Complete Decision Logic Test")
    print("=" * 70)

    analyzer = RHDAnalyzer()

    # Test cases
    test_cases = [
        {
            'name': 'RhD+ case (RHD1, short sequence)',
            'sequence': RHD1_REFERENCE[:400],  # 400 bp < 800 bp threshold
            'expected_status': 'RhD+',
            'description': 'Short RHD1 amplicon indicates gene present'
        },
        {
            'name': 'RhD- case (RHD1, long sequence)',
            'sequence': RHD1_REFERENCE[:900],  # 900 bp >= 800 bp threshold
            'expected_status': 'RhD-',
            'description': 'Long RHD1 amplicon indicates gene absent'
        },
        {
            'name': 'RhD+ case (RHD456, high identity)',
            'sequence': RHD456_REFERENCE,  # 100% match
            'expected_status': 'RhD+',
            'description': 'High identity RHD456 indicates normal D antigen'
        },
        {
            'name': 'Perfect RHD1 match',
            'sequence': RHD1_REFERENCE,  # Full 951 bp
            'expected_status': 'RhD-',
            'description': 'Full RHD1 length indicates gene absence'
        },
    ]

    print("\n[Running Test Cases]\n")
    passed = 0
    failed = 0

    for i, test in enumerate(test_cases, 1):
        print(f"Test {i}: {test['name']}")
        print(f"  Description: {test['description']}")
        print(f"  Sequence length: {len(test['sequence'])} bp")

        result = analyzer.analyze(test['sequence'])

        print(f"  Region detected: {result['region']}")
        print(f"  Query length: {result['query_length']} bp")
        print(f"  Identity: {result['identity']}%")
        print(f"  RhD Status: {result['rhd_status']}")
        print(f"  Expected: {test['expected_status']}")
        print(f"  Reason: {result['reason']}")

        if result['rhd_status'] == test['expected_status']:
            print(f"  Result: [PASS]")
            passed += 1
        else:
            print(f"  Result: [FAIL]")
            failed += 1

        print()

    print("=" * 70)
    print(f"Test Summary: {passed} passed, {failed} failed out of {len(test_cases)}")
    print("=" * 70)

    # Verify field availability
    print("\n[Field Availability Check]\n")

    test_seq = RHD1_REFERENCE[:400]
    result = analyzer.analyze(test_seq)

    expected_fields = [
        'variants',
        'identity',
        'score',
        'strand',
        'query_length',
        'region',
        'rhd_status',
        'reason',
        'reference_description'
    ]

    print("Expected fields by main.py:")
    all_present = True
    for field in expected_fields:
        if field in result:
            print(f"  [OK] {field}: {result[field]}")
        else:
            print(f"  [MISSING] {field}")
            all_present = False

    print()
    return passed == len(test_cases) and all_present

if __name__ == "__main__":
    success = test_rhd_decision_logic()
    print(f"\n[Final Result] {'SUCCESS' if success else 'FAILURE'}")
    sys.exit(0 if success else 1)
