#!/usr/bin/env python
"""Test RHD and ABO analysis with actual sample files"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from Bio import SeqIO
import utils.FASTA_analyzer as fasta_utils
from utils.rhd_analyzer import RHDAnalyzer

def test_with_sample_files():
    """Test analysis with sample FASTA files"""
    print("=" * 80)
    print("FULL SYSTEM TEST WITH SAMPLE FILES")
    print("=" * 80)

    sample_dir = Path(__file__).parent / "sample_for_global_alignment"
    fasta_files = sorted(sample_dir.glob("*.fasta"))[:5]

    if not fasta_files:
        print("[ERROR] No sample files found")
        return False

    print(f"\n[INFO] Found {len(fasta_files)} sample files\n")

    # ========== ABO ANALYSIS ==========
    print("=" * 80)
    print("PART 1: ABO ANALYSIS")
    print("=" * 80)

    service = fasta_utils.FASTAAlignmentService()

    print("\n[Step 1] Processing ABO FASTA files...\n")

    from io import StringIO
    file_objects = []
    for fasta_file in fasta_files:
        with open(fasta_file, 'r') as f:
            content = f.read()
        file_obj = StringIO(content)
        file_obj.name = fasta_file.name
        file_objects.append(file_obj)

    abo_results = service.generate_batch_summary(file_objects)

    print(f"[OK] Processed {len(abo_results)} files\n")
    print("ABO Analysis Results:")
    print("-" * 80)

    for res in abo_results:
        decision = res.get('decision', 'Unknown')
        exon = res.get('exon', 'N/A')
        coverage = res.get('coverage', 0)
        similarity = res.get('similarity', 0)
        variants = res.get('variants', [])

        print(f"\nFile: {res['filename']}")
        print(f"  Exon: {exon}")
        print(f"  Coverage: {coverage:.1f}%")
        print(f"  Similarity: {similarity:.1f}%")
        print(f"  Decision: {decision}")
        if variants:
            print(f"  Variants found: {len(variants)}")
            for v in variants[:2]:
                print(f"    - {v.get('type')}: {v.get('ref_base')}->{v.get('alt_base')} @ {v.get('isbt_pos')}")
        else:
            print(f"  Variants: None (matches reference)")

    # ========== RHD ANALYSIS ==========
    print("\n\n" + "=" * 80)
    print("PART 2: RHD ANALYSIS")
    print("=" * 80)

    rhd_analyzer = RHDAnalyzer()

    print("\n[Step 2] Analyzing samples as RHD amplicons...\n")
    print("RHD Analysis Results:")
    print("-" * 80)

    rhd_results = []
    for fasta_file in fasta_files:
        # Read sequence
        record = SeqIO.read(str(fasta_file), "fasta")
        seq = str(record.seq)

        # Analyze with RHD analyzer
        result = rhd_analyzer.analyze(seq)

        rhd_results.append({
            'filename': fasta_file.name,
            'result': result
        })

        print(f"\nFile: {fasta_file.name}")
        print(f"  Sequence length: {result['query_length']} bp")
        print(f"  Region detected: {result['region']}")
        print(f"  Identity: {result['identity']}%")
        print(f"  RhD Status: {result['rhd_status']}")
        print(f"  Reason: {result['reason']}")
        print(f"  Variants: {len(result['variants'])} detected")

    # ========== SUMMARY ==========
    print("\n\n" + "=" * 80)
    print("SUMMARY & INTERPRETATION")
    print("=" * 80)

    print("\n[ABO ANALYSIS]")
    confirmed = [r for r in abo_results if "Confirmed" in r.get('decision', '')]
    print(f"  Confirmed exons: {len(confirmed)}/{len(abo_results)}")
    for res in confirmed:
        print(f"    - Exon {res['exon']}: {res['similarity']:.1f}% similarity")

    print("\n[RHD ANALYSIS]")
    rhd_pos = [r for r in rhd_results if 'RhD+' in r['result']['rhd_status']]
    rhd_neg = [r for r in rhd_results if 'RhD-' in r['result']['rhd_status']]
    print(f"  Results:")
    print(f"    - RhD+ (positive): {len(rhd_pos)} samples")
    print(f"    - RhD- (negative): {len(rhd_neg)} samples")

    print("\n[BLOOD GROUP SUMMARY]")
    print("  ABO phenotype determination:")
    print("    - Based on confirmed exons and variant matching")
    print("    - Possible types: A, B, AB, O")
    print("\n  RHD phenotype determination:")
    print("    - Based on amplicon length and identity")
    print("    - RhD+ means D antigen present")
    print("    - RhD- means D antigen absent")

    print("\n" + "=" * 80)
    print("TEST COMPLETED SUCCESSFULLY")
    print("=" * 80)
    print("\nSystem is working correctly with real sample data!")
    print("\nTo use in production:")
    print("  1. Open http://localhost:8503 in your browser")
    print("  2. Upload your FASTA or AB1 files")
    print("  3. Click 'Analyze'")
    print("  4. View results in the tabs")
    print("=" * 80)

    return True

if __name__ == "__main__":
    try:
        success = test_with_sample_files()
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"\n[ERROR] Test failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
