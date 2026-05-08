#!/usr/bin/env python
"""Test RHD analysis with real patient sample files from term paper data"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from Bio import SeqIO
from utils.rhd_analyzer import RHDAnalyzer

def test_real_rhd_samples():
    """Analyze real RHD1 and RHD456 samples from patient data"""
    print("=" * 80)
    print("REAL PATIENT RHD ANALYSIS TEST")
    print("=" * 80)

    test_data_path = Path(r"C:\Users\ExPertComputer\Downloads\ผล blood group NGS_term paper 2568 ใบปอ")

    # Test Ake patient RHD1 samples
    ake_rhd1_path = test_data_path / "Ake นศ.ป.เอก อ.น้ำผึ้ง ทำ uveitis" / "RHD1"
    ake_rhd456_path = test_data_path / "Ake นศ.ป.เอก อ.น้ำผึ้ง ทำ uveitis" / "RHD456"

    analyzer = RHDAnalyzer()

    print("\n[PATIENT: Ake นศ.ป.เอก อ.น้ำผึ้ง ทำ uveitis]")
    print("-" * 80)

    # Test RHD1 samples
    print("\n[RHD1 ANALYSIS] (< 800 bp = RhD+, >= 800 bp = RhD-)")
    print("-" * 80)

    if ake_rhd1_path.exists():
        rhd1_files = sorted(ake_rhd1_path.glob("*.fasta"))[:5]  # Test first 5 files

        for fasta_file in rhd1_files:
            try:
                record = SeqIO.read(str(fasta_file), "fasta")
                seq = str(record.seq)

                result = analyzer.analyze(seq)

                print(f"\nFile: {fasta_file.name}")
                print(f"  Sequence length: {result['query_length']} bp")
                print(f"  Region detected: {result['region']}")
                print(f"  Identity: {result['identity']:.1f}%")
                print(f"  RhD Status: {result['rhd_status']}")
                print(f"  Reason: {result['reason']}")
                if result['variants']:
                    print(f"  Variants: {len(result['variants'])} detected")

            except Exception as e:
                print(f"\n[ERROR] Failed to analyze {fasta_file.name}: {e}")
    else:
        print(f"\n[WARNING] RHD1 directory not found: {ake_rhd1_path}")

    # Test RHD456 samples
    print("\n\n[RHD456 ANALYSIS] (>= 90% identity = RhD+, < 90% = RhD-)")
    print("-" * 80)

    if ake_rhd456_path.exists():
        rhd456_files = sorted(ake_rhd456_path.glob("*.fasta"))[:5]

        for fasta_file in rhd456_files:
            try:
                record = SeqIO.read(str(fasta_file), "fasta")
                seq = str(record.seq)

                result = analyzer.analyze(seq)

                print(f"\nFile: {fasta_file.name}")
                print(f"  Sequence length: {result['query_length']} bp")
                print(f"  Region detected: {result['region']}")
                print(f"  Identity: {result['identity']:.1f}%")
                print(f"  RhD Status: {result['rhd_status']}")
                print(f"  Reason: {result['reason']}")
                if result['variants']:
                    print(f"  Variants: {len(result['variants'])} detected")

            except Exception as e:
                print(f"\n[ERROR] Failed to analyze {fasta_file.name}: {e}")
    else:
        print(f"\n[WARNING] RHD456 directory not found: {ake_rhd456_path}")

    # Summary
    print("\n\n" + "=" * 80)
    print("CLINICAL INTERPRETATION")
    print("=" * 80)
    print("\nRhD+ (Positive):")
    print("  - Patient has RHD gene present")
    print("  - Safe to transfuse as RhD+")
    print("  - If pregnant: monitor for alloimmunization")
    print("\nRhD- (Negative):")
    print("  - Patient lacks RHD gene (D antigen absent)")
    print("  - Must receive RhD- blood products")
    print("  - If pregnant: may need anti-D prophylaxis")

    print("\n" + "=" * 80)
    print("TEST COMPLETED")
    print("=" * 80)
    return True

if __name__ == "__main__":
    try:
        success = test_real_rhd_samples()
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"\n[ERROR] Test failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
