#!/usr/bin/env python
"""Test ABO allele identification functionality"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import utils.FASTA_analyzer as fasta_utils
import utils.abo_identifier as abo_utils

def test_allele_identification():
    """Test ABO allele identification from FASTA analysis"""
    print("=" * 70)
    print("Testing ABO Allele Identification")
    print("=" * 70)

    # Run FASTA analysis
    service = fasta_utils.FASTAAlignmentService()
    sample_dir = Path(__file__).parent / "sample_for_global_alignment"
    fasta_files = list(sample_dir.glob("*.fasta"))[:5]

    if not fasta_files:
        print("[FAIL] No sample files found")
        return False

    from io import StringIO
    file_objects = []
    for fasta_file in fasta_files:
        with open(fasta_file, 'r') as f:
            content = f.read()
        file_obj = StringIO(content)
        file_obj.name = fasta_file.name
        file_objects.append(file_obj)

    print("[OK] Analyzing FASTA files...")
    try:
        results = service.generate_batch_summary(file_objects)
        confirmed_results = [r for r in results if "Confirmed" in r.get('decision', "")]

        print(f"[OK] {len(confirmed_results)} exons confirmed")

        if not confirmed_results:
            print("[INFO] No confirmed exons - cannot identify alleles")
            return True

        # Initialize ABOIdentifier
        abo_id = abo_utils.ABOIdentifier("ABO")
        print("[OK] ABOIdentifier initialized")

        # Test allele identification
        print("\n" + "=" * 70)
        print("Allele Identification Results:")
        print("=" * 70)

        # Group by exon
        files_by_exon = {}
        for res in confirmed_results:
            ex = res['exon_number']
            if ex not in files_by_exon:
                files_by_exon[ex] = []
            files_by_exon[ex].append(res)

        sorted_exons = sorted(files_by_exon.keys())
        print(f"\nDetected Exons: {', '.join([f'Exon {e}' for e in sorted_exons])}")

        # Test each exon individually
        print("\n[Individual Exon Analysis]")
        for res in confirmed_results:
            exon_num = res['exon_number']
            filename = res['filename']
            variants = res.get('variants', [])

            print(f"\n  Exon {exon_num} - {filename}")

            if variants:
                variant_list = [f"{v['type']}@{v['isbt_pos']}" for v in variants]
                print(f"    Variants: {', '.join(variant_list)}")

                # Identify alleles
                try:
                    preds, unk, var_names, node_map = abo_utils.identify_abo_alleles(
                        {'exon_alignments': [res]}, abo_identifier=abo_id)

                    if preds:
                        top_alleles = [list(d.keys())[0] for d in preds[:5]]
                        print(f"    Possible alleles: {', '.join(top_alleles)}" +
                              (f" (and {len(preds)-5} more)" if len(preds) > 5 else ""))
                    else:
                        print(f"    No matching alleles found")

                    if unk:
                        print(f"    Unknown variants: {len(unk)}")
                except Exception as e:
                    print(f"    Error identifying alleles: {e}")
            else:
                # Reference match
                print(f"    No variants (reference match)")
                all_alleles = abo_id.get_all_alleles()
                excluded = abo_id.get_alleles_with_variant_in_exon(exon_num)
                print(f"    Consistent with {len(all_alleles) - len(excluded)}/{len(all_alleles)} alleles")

        return True

    except Exception as e:
        print(f"[FAIL] Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_allele_identification()
    sys.exit(0 if success else 1)
