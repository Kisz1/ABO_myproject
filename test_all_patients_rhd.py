#!/usr/bin/env python
"""Test RHD analysis with all patient samples"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from Bio import SeqIO
from utils.rhd_analyzer import RHDAnalyzer

def test_all_patients():
    """Analyze RHD1 samples from all available patients"""
    print("=" * 80)
    print("MULTI-PATIENT RHD ANALYSIS TEST")
    print("=" * 80)

    test_data_path = Path(r"C:\Users\ExPertComputer\Downloads\ผล blood group NGS_term paper 2568 ใบปอ")

    if not test_data_path.exists():
        print(f"[ERROR] Test data path not found: {test_data_path}")
        return False

    analyzer = RHDAnalyzer()
    patient_dirs = [d for d in test_data_path.iterdir() if d.is_dir()]

    if not patient_dirs:
        print(f"[ERROR] No patient directories found in {test_data_path}")
        return False

    print(f"\n[INFO] Found {len(patient_dirs)} patient(s)\n")

    summary_results = {}

    for patient_dir in patient_dirs:
        patient_name = patient_dir.name
        rhd1_path = patient_dir / "RHD1"

        if not rhd1_path.exists():
            continue

        print(f"\n{'='*80}")
        print(f"PATIENT: {patient_name}")
        print(f"{'='*80}")

        rhd_files = sorted(rhd1_path.glob("*.fasta"))
        rhd_plus_count = 0
        rhd_minus_count = 0

        if not rhd_files:
            print(f"  [WARNING] No FASTA files found")
            continue

        print(f"\n  [RHD1 Amplicon Analysis] ({len(rhd_files)} files)")
        print(f"  {'-'*76}")

        for i, fasta_file in enumerate(rhd_files[:3], 1):  # Show first 3 files per patient
            try:
                record = SeqIO.read(str(fasta_file), "fasta")
                seq = str(record.seq)
                result = analyzer.analyze(seq)

                rhd_status = result['rhd_status']
                length = result['query_length']

                if 'RhD+' in rhd_status:
                    rhd_plus_count += 1
                    status_label = "[RhD+]"
                elif 'RhD-' in rhd_status:
                    rhd_minus_count += 1
                    status_label = "[RhD-]"
                else:
                    status_label = "[????]"

                print(f"\n  {i}. {fasta_file.name}")
                print(f"     Length: {length} bp -> {status_label}")
                print(f"     Identity: {result['identity']:.1f}%")

            except Exception as e:
                print(f"\n  {i}. {fasta_file.name}")
                print(f"     [ERROR] {str(e)[:50]}")

        total = len(rhd_files)
        print(f"\n  SUMMARY for {patient_name}:")
        print(f"    Total RHD1 files: {total}")
        print(f"    RhD+ results: {rhd_plus_count if rhd_plus_count > 0 else 'None analyzed'}")
        print(f"    RhD- results: {rhd_minus_count if rhd_minus_count > 0 else 'None analyzed'}")

        summary_results[patient_name] = {
            'total': total,
            'rhd_plus': rhd_plus_count,
            'rhd_minus': rhd_minus_count
        }

    # Final summary
    print(f"\n\n{'='*80}")
    print("FINAL SUMMARY")
    print(f"{'='*80}\n")

    for patient_name, results in summary_results.items():
        print(f"{patient_name}:")
        if results['rhd_plus'] > 0 and results['rhd_minus'] == 0:
            print(f"  RhD Status: RhD+ (POSITIVE - D antigen present)")
        elif results['rhd_minus'] > 0 and results['rhd_plus'] == 0:
            print(f"  RhD Status: RhD- (NEGATIVE - D antigen absent)")
        elif results['rhd_plus'] > 0 and results['rhd_minus'] > 0:
            print(f"  RhD Status: MIXED (variant or mixed genotype)")
        else:
            print(f"  RhD Status: INCONCLUSIVE")
        print()

    print(f"{'='*80}\n")
    return True

if __name__ == "__main__":
    try:
        success = test_all_patients()
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"\n[ERROR] Test failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
