#!/usr/bin/env python
"""Simple test: RhD+ and RhD- determination"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from utils.rhd_analyzer import RHDAnalyzer, RHD1_REFERENCE

def test_rhd_phenotype():
    """Test RhD+ and RhD- phenotype determination"""
    print("=" * 70)
    print("RHD Phenotype Determination Test")
    print("=" * 70)

    analyzer = RHDAnalyzer()

    # Test cases for RhD+ and RhD-
    test_cases = [
        {
            'name': 'RhD+ (RHD1 < 800 bp)',
            'sequence': RHD1_REFERENCE[:400],
            'expected': 'RhD+',
        },
        {
            'name': 'RhD+ (RHD1 < 800 bp, 600 bp)',
            'sequence': RHD1_REFERENCE[:600],
            'expected': 'RhD+',
        },
        {
            'name': 'RhD+ (RHD1 < 800 bp, 799 bp)',
            'sequence': RHD1_REFERENCE[:799],
            'expected': 'RhD+',
        },
        {
            'name': 'RhD- (RHD1 >= 800 bp, 800 bp)',
            'sequence': RHD1_REFERENCE[:800],
            'expected': 'RhD-',
        },
        {
            'name': 'RhD- (RHD1 >= 800 bp, 900 bp)',
            'sequence': RHD1_REFERENCE[:900],
            'expected': 'RhD-',
        },
        {
            'name': 'RhD- (RHD1 full, 951 bp)',
            'sequence': RHD1_REFERENCE,
            'expected': 'RhD-',
        },
    ]

    print("\n[Testing RhD+ and RhD- Determination]\n")
    passed = 0
    failed = 0

    for test in test_cases:
        result = analyzer.analyze(test['sequence'])
        status = "[PASS]" if result['rhd_status'] == test['expected'] else "[FAIL]"

        if result['rhd_status'] == test['expected']:
            passed += 1
        else:
            failed += 1

        print(f"{status} {test['name']}")
        print(f"       Length: {result['query_length']} bp -> {result['rhd_status']} (expected {test['expected']})")
        print()

    # Summary
    print("=" * 70)
    print(f"Results: {passed} PASS, {failed} FAIL out of {len(test_cases)}")
    print("=" * 70)

    if failed == 0:
        print("\n[SUCCESS] RhD+ and RhD- determination is working correctly!")
        return True
    else:
        print(f"\n[WARNING] {failed} test(s) failed")
        return False

if __name__ == "__main__":
    success = test_rhd_phenotype()
    sys.exit(0 if success else 1)
