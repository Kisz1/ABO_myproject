#!/usr/bin/env python
"""Test script for ABO analysis functionality"""

import sys
import os
from pathlib import Path

# Add utils to path
sys.path.insert(0, str(Path(__file__).parent))

import utils.FASTA_analyzer as fasta_utils

def test_abo_fasta_analysis():
    """Test ABO FASTA file analysis"""
    print("=" * 70)
    print("Testing ABO FASTA Analysis")
    print("=" * 70)

    # Initialize the service
    service = fasta_utils.FASTAAlignmentService()
    print("[OK] FASTAAlignmentService initialized")

    # Get sample files
    sample_dir = Path(__file__).parent / "sample_for_global_alignment"
    fasta_files = list(sample_dir.glob("*.fasta"))[:5]  # Test with first 5 files

    if not fasta_files:
        print("[FAIL] No FASTA sample files found")
        return False

    print(f"[OK] Found {len(fasta_files)} sample FASTA files")

    # Convert to file-like objects
    from io import StringIO
    file_objects = []
    for fasta_file in fasta_files:
        with open(fasta_file, 'r') as f:
            content = f.read()

        # Create a file-like object
        file_obj = StringIO(content)
        file_obj.name = fasta_file.name
        file_objects.append(file_obj)

    # Run analysis
    print("\nRunning batch analysis on sample files...")
    try:
        results = service.generate_batch_summary(file_objects)
        print(f"[OK] Analysis completed: {len(results)} files processed")

        # Display results
        print("\n" + "=" * 70)
        print("Analysis Results Summary:")
        print("=" * 70)

        for i, result in enumerate(results, 1):
            print(f"\n{i}. File: {result.get('filename', 'N/A')}")
            print(f"   Exon: {result.get('exon', 'N/A')}")
            print(f"   Coverage: {result.get('coverage', 0):.1f}%")
            print(f"   Similarity: {result.get('similarity', 0):.1f}%")
            print(f"   Decision: {result.get('decision', 'N/A')}")

            # Show variants if confirmed
            if "Confirmed" in result.get('decision', ""):
                variants = result.get('variants', [])
                if variants:
                    print(f"   Variants detected: {len(variants)}")
                    for v in variants[:3]:  # Show first 3
                        print(f"     - {v.get('type', 'N/A')}: {v.get('ref_base', '')}->{v.get('alt_base', '')} at pos {v.get('isbt_pos', 'N/A')}")
                    if len(variants) > 3:
                        print(f"     ... and {len(variants) - 3} more")
                else:
                    print("   No variants detected (matches reference)")

        # Identify confirmed exons
        confirmed = [r for r in results if "Confirmed" in r.get('decision', "")]
        print("\n" + "=" * 70)
        print(f"Confirmed Exons: {len(confirmed)}/{len(results)}")
        if confirmed:
            exons = [r['exon'] for r in confirmed]
            print(f"  - Exon {exons[0] if exons else 'N/A'}" +
                  (f" through {exons[-1]}" if len(exons) > 1 else ""))

        return True

    except Exception as e:
        print(f"[FAIL] Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_abo_fasta_analysis()
    sys.exit(0 if success else 1)
