"""
Proper FASTA Alignment Service
Uses BioPython's pairwise2 module for gap-aware sequence alignment
to eliminate coordinate shift issues and false SNP cascades.
"""

import json
import glob
import os
from typing import Dict, List, Any 
import utils.referece_loader as rl

 
GLOBAL = 1000
DEBUG = False

try:
    from Bio.Align import PairwiseAligner
    from Bio.Seq import Seq
    from Bio import SeqIO
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    PairwiseAligner = None
    print("BioPython not available. Install with: pip install biopython")


class FASTAAlignmentService:
    """
    Service for proper sequence alignment using BioPython's gap-aware algorithms.
    Replaces simple position-by-position comparison to eliminate false SNPs from indels.
    """

    def __init__(self, gene="ABO"):
        """Initialize the alignment service"""
        if not BIOPYTHON_AVAILABLE:
            raise ImportError(
                "BioPython is required. Install with: pip install biopython")

        # Load NCBI ABO reference data
        if gene == "ABO":
            self.abo_reference = rl.ReferenceLoader().load_abo_reference()

        # Alignment parameters (optimized for DNA sequences)
        self.local_alignment_params = {
            'match_score': 2,        # Positive score for matches
            'mismatch_score': -1,    # Penalty for mismatches
            'gap_open_penalty': -2,  # Higher penalty for opening a gap
            'gap_extend_penalty': -0.5,  # Higher penalty for extending a gap
            'one_alignment_only': True  # Return only the best alignment
        }

        self.global_alignment_params = {
            'match_score': 2,
            'mismatch_score': -0.5,
            'gap_open_penalty': -2,
            'gap_extend_penalty': -0.5,
            'one_alignment_only': True
        }

        
    def getABO_ref(self, data_field):
        return self.abo_reference[data_field]

    def _extract_aligned_sequences(self, alignment):
        """Extract aligned sequences with gaps from BioPython alignment object"""
        alignment_str = str(alignment)
        lines = alignment_str.strip().split('\n')

        # if len(lines) >= 3:
        #     # Parse the alignment format:
        #     # target            0 ACGT--ACGT  8
        #     # query             0 ACGTTAACGT 10
        #     target_parts = lines[0].split()
        #     query_parts = lines[2].split()

        #     if len(target_parts) >= 3 and len(query_parts) >= 3:
        #         aligned_target = target_parts[2]  # The sequence with gaps
        #         aligned_query = query_parts[2]    # The sequence with gaps
        #         return aligned_query, aligned_target
        aligned_target = ""
        aligned_query = ""

        for i in range(0, len(lines), 4):

            target_parts = lines[i].split()
            query_parts = lines[i+2].split()

            if len(target_parts) >= 3 and len(query_parts) >= 3:
                aligned_target += target_parts[2]
                aligned_query += query_parts[2]

        return aligned_query, aligned_target

        # Fallback if parsing fails
       # return "", ""
 
    def align_sequence_to_exon(self, query_sequence: str, exon_number: int) -> Dict:
        """
        Perform proper local alignment between query sequence and specific exon.

        Args:
            query_sequence: The target sequence to align
            exon_number: Which ABO exon to align against (1-7)

        Returns:
            Dictionary with alignment results and variant information
        """
        if not self.abo_reference:
            return {'error': 'ABO reference not loaded'}

        # Get exon sequence from reference
        exon_data = None
        for exon in self.abo_reference['exons']:
            if exon['exon_number'] == exon_number:
                exon_data = exon
                break

        if not exon_data:
            return {'error': f'Exon {exon_number} not found'}

        exon_sequence = exon_data['sequence']
        

        # Create PairwiseAligner with Smith-Waterman (local) alignment settings
        if not BIOPYTHON_AVAILABLE or PairwiseAligner is None:
            return {'error': 'BioPython not available'}

        aligner = PairwiseAligner()  # type: ignore
        if len(exon_sequence) > GLOBAL:
            aligner.mode = 'global'
            aligner.match_score = self.global_alignment_params['match_score']
            aligner.mismatch_score = self.global_alignment_params['mismatch_score']
            aligner.open_gap_score = self.global_alignment_params['gap_open_penalty']
            aligner.extend_gap_score = self.global_alignment_params['gap_extend_penalty']
        else:
            aligner.mode = 'local'
            aligner.match_score = self.local_alignment_params['match_score']
            aligner.mismatch_score = self.local_alignment_params['mismatch_score']
            aligner.open_gap_score = self.local_alignment_params['gap_open_penalty']
            aligner.extend_gap_score = self.local_alignment_params['gap_extend_penalty']

        # Perform local alignment
        if len(query_sequence) > len(exon_sequence)*80/100:  # At least 80% length
            alignments = aligner.align(query_sequence, exon_sequence)

            # Check if there are too many alignments (indicating poor match)
            try:
                alignment_count = len(alignments)
                if alignment_count == 0:
                    return {'error': 'No significant alignment found'}
            except OverflowError:
                # Too many optimal alignments indicates sequences are very different
                # Get the best alignment
                return {'error': 'Sequences too divergent - too many possible alignments'}
            best_alignment = alignments[0]
            score = best_alignment.score  # type: ignore

            # Extract aligned sequences with gaps from the alignment object
            aligned_ref, aligned_query = self._extract_aligned_sequences(
                best_alignment)

            # Calculate alignment statistics
            alignment_length = len(aligned_query)
            matches = sum(1 for q, r in zip(aligned_query, aligned_ref)
                        if q == r and q != '-' and r != '-')
            similarity = matches / alignment_length if alignment_length > 0 else 0

            # Get alignment coordinates
            query_start = best_alignment.coordinates[1][0]  # type: ignore
            query_end = best_alignment.coordinates[1][-1]  # type: ignore
            ref_start = best_alignment.coordinates[0][0]  # type: ignore
            ref_end = best_alignment.coordinates[0][-1]  # type: ignore

            # Parse variants from the alignment
            variants = self._extract_variants_from_alignment(
                aligned_query, aligned_ref, query_start, ref_start, exon_data
            )

            return {
                'exon_number': exon_number,
                'alignment_score': score,
                'similarity': similarity,
                'alignment_length': alignment_length,
                'matches': matches,
                'query_start': query_start,
                'query_end': query_end,
                'ref_start': ref_start,
                'ref_end': ref_end,
                'aligned_query': aligned_query,
                'aligned_reference': aligned_ref,
                'variants': variants,
                'exon_info': {
                    'genomic_start': exon_data['start'],
                    'genomic_end': exon_data['end'],
                    'length': len(exon_sequence)
                }
            }
        else:
            return {'error': 'Query sequence shorter than exon sequence'}

    def _extract_variants_from_alignment(self, aligned_query: str, aligned_ref: str,
                                         query_start: int, ref_start: int, exon_data: Dict) -> List[Dict]:
        """
        Extract variants from a proper alignment, handling gaps correctly.

        Args:
            aligned_query: Query sequence with gaps
            aligned_ref: Reference sequence with gaps
            query_start: Start position in query
            exon_data: Exon information with genomic coordinates

        Returns:
            List of variants with proper coordinates
        """
        variants = []
        query_pos = query_start
        ref_pos = 0

        i = 0
        while i < len(aligned_query):
            query_base = aligned_query[i]
            ref_base = aligned_ref[i]

            if query_base == '-':
                # Deletion in query (insertion in reference)
                # Count consecutive deletions
                del_length = 0
                del_sequence = ""
                j = i
                while j < len(aligned_query) and aligned_query[j] == '-':
                    del_sequence += aligned_ref[j]
                    del_length += 1
                    j += 1

                genomic_pos = exon_data['start'] + ref_pos
                isbt_pos = exon_data['cds_start'] + ref_pos

                variants.append({
                    'type': 'deletion',
                    'query_pos': query_pos,
                    'ref_pos': ref_pos,
                    'genomic_pos': genomic_pos,
                    'isbt_pos': isbt_pos,
                    'length': del_length,
                    'deleted_sequence': del_sequence,
                    'change': f"del({del_sequence})",
                    'exon': exon_data['exon_number']
                })

                ref_pos += del_length
                i += del_length

            elif ref_base == '-':
                # Insertion in query (deletion in reference)
                # Count consecutive insertions
                ins_length = 0
                ins_sequence = ""
                j = i
                while j < len(aligned_query) and aligned_ref[j] == '-':
                    ins_sequence += aligned_query[j]
                    ins_length += 1
                    j += 1

                genomic_pos = exon_data['start'] + ref_pos
                isbt_pos = exon_data['cds_start'] + ref_pos

                variants.append({
                    'type': 'insertion',
                    'query_pos': query_pos,
                    'ref_pos': ref_pos,
                    'genomic_pos': genomic_pos,
                    'isbt_pos': isbt_pos,
                    'length': ins_length,
                    'inserted_sequence': ins_sequence,
                    'change': f"ins({ins_sequence})",
                    'exon': exon_data['exon_number']
                })

                query_pos += ins_length
                i += ins_length

            else:
                # Match or mismatch
                if query_base != ref_base:
                    # SNP
                    genomic_pos = exon_data['start'] + ref_pos
                    isbt_pos = exon_data['cds_start'] + ref_pos

                    variants.append({
                        'type': 'SNP',
                        'query_pos': query_pos,
                        'ref_pos': ref_pos,
                        'genomic_pos': genomic_pos,
                        'isbt_pos': isbt_pos,
                        'ref_base': ref_base,
                        'alt_base': query_base,
                        'change': f"{ref_base}>{query_base}",
                        'exon': exon_data['exon_number']
                    })

                query_pos += 1
                ref_pos += 1
                i += 1

        return variants

    def analyze_multi_exon_sequence(self, query_sequence: str,
                                    exon_combination: List[int]) -> Dict:
        """
        Analyze a sequence that may contain multiple exons using proper alignment.

        Args:
            query_sequence: The sequence to analyze
            exon_combination: List of exon numbers to try (e.g., [5, 6, 7])

        Returns:
            Analysis results with proper alignment-based variant calling
        """
        results = {
            'query_length': len(query_sequence),
            'exon_combination': exon_combination,
            'exon_alignments': [],
            'total_variants': 0,
            'variant_summary': {'SNPs': 0, 'insertions': 0, 'deletions': 0},
            'overall_similarity': 0.0
        }

        total_score = 0
        total_length = 0
        all_variants = []

        if DEBUG: print(f"\n=== PROPER ALIGNMENT ANALYSIS ===")
        if DEBUG: print(f"Query length: {len(query_sequence)} bp")
        if DEBUG: print(f"Analyzing exons: {exon_combination}")

        # Align to each exon separately
        for exon_num in exon_combination:
            if DEBUG: print(f"\nAligning to exon {exon_num}...")

            alignment_result = self.align_sequence_to_exon(
                query_sequence, exon_num)

            if 'error' in alignment_result:
                print(f"  Error: {alignment_result['error']}")
                continue

            results['exon_alignments'].append(alignment_result)

            # Accumulate statistics
            total_score += alignment_result['alignment_score']
            total_length += alignment_result['alignment_length']

            variants = alignment_result['variants']
            all_variants.extend(variants)

            # Count variant types
            for variant in variants:
                variant_type = variant['type']
                if variant_type == 'SNP':
                    results['variant_summary']['SNPs'] += 1
                elif variant_type == 'insertion':
                    results['variant_summary']['insertions'] += 1
                elif variant_type == 'deletion':
                    results['variant_summary']['deletions'] += 1

            if DEBUG: print(f"  Score: {alignment_result['alignment_score']:.1f}, "
                  f"Similarity: {alignment_result['similarity']:.3f}, "
                  f"Variants: {len(variants)}")

        # Calculate overall statistics
        #results['total_variants'] = len(all_variants)
        #results['all_variants'] = all_variants
        #results['overall_similarity'] = total_score / \
        #    total_length if total_length > 0 else 0
        if DEBUG:
            
            print(f"\n=== ALIGNMENT SUMMARY ===")
            print(f"Total variants: {results['total_variants']}")
            print(f"  SNPs: {results['variant_summary']['SNPs']}")
            print(f"  Insertions: {results['variant_summary']['insertions']}")
            print(f"  Deletions: {results['variant_summary']['deletions']}")
            print(f"Overall similarity: {results['overall_similarity']:.3f}")
       
        if results['exon_alignments']!=[]:
            return results
        else:
            return {'error': 'No exons were successfully aligned'}

    def create_coordinate_mapping_from_alignment(self, alignment_results: Dict) -> Dict:
        """
        Create coordinate mapping system from proper alignment results.
        Maps reporting coordinates to genomic coordinates while handling gaps.

        Args:
            alignment_results: Results from analyze_multi_exon_sequence

        Returns:
            Coordinate mapping compatible with existing system
        """
        coordinate_mapping = {
            'exon_boundaries': [],
            'total_positions': 0
        }

        current_reporting_pos = 1

        for alignment in alignment_results['exon_alignments']:
            exon_info = alignment['exon_info']
            exon_length = exon_info['length']

            boundary = {
                'exon_number': alignment['exon_number'],
                # Remove gaps
                'sequence': alignment['aligned_reference'].replace('-', ''),
                'reporting_start': current_reporting_pos,
                'reporting_end': current_reporting_pos + exon_length - 1,
                'genomic_start': exon_info['genomic_start'],
                'genomic_end': exon_info['genomic_end'],
                'cds_start': exon_info['cds_start'],
                'cds_end': exon_info['cds_end']
            }

            coordinate_mapping['exon_boundaries'].append(boundary)
            current_reporting_pos += exon_length

        coordinate_mapping['total_positions'] = current_reporting_pos - 1

        return coordinate_mapping

    def format_variants_for_reporting(self, variants: List[Dict],
                                      coordinate_mapping: Dict) -> List[Dict]:
        """
        Format variants for consistent reporting with coordinate mapping.

        Args:
            variants: Raw variants from alignment
            coordinate_mapping: Coordinate mapping system

        Returns:
            Formatted variants compatible with existing reporting system
        """
        formatted_variants = []

        for variant in variants:
            # Find the exon boundary for this variant
            exon_boundary = None
            for boundary in coordinate_mapping['exon_boundaries']:
                if boundary['exon_number'] == variant['exon']:
                    exon_boundary = boundary
                    break

            if not exon_boundary:
                continue

            # Calculate reporting coordinate
            offset_in_exon = variant['ref_pos']
            reporting_coord = exon_boundary['reporting_start'] + offset_in_exon

            formatted_variant = {
                'type': variant['type'],
                'change': variant['change'],
                'exon': variant['exon'],
                'coordinates': {
                    'reporting': reporting_coord,
                    'exon_relative': offset_in_exon + 1,  # 1-based
                    'genomic': variant['genomic_pos']
                }
            }

            # Add base information for SNPs
            if variant['type'] == 'SNP':
                formatted_variant['ref_base'] = variant['ref_base']
                formatted_variant['alt_base'] = variant['alt_base']

            formatted_variants.append(formatted_variant)

        formatted_variants.sort(key=lambda v: v['coordinates']['reporting'])

        return formatted_variants

    def align_to_genomic_reference(self, query_sequence: str) -> Dict:
        """
        Align a query sequence to the full genomic reference sequence.
        Checks both forward and reverse complement orientations.
        
        Args:
            query_sequence: The sequence to align
            
        Returns:
            Dictionary with alignment results including 'orientation'
        """
        if not self.abo_reference or 'full_sequence' not in self.abo_reference:
            return {'error': 'Genomic reference not available'}
            
        ref_sequence = self.abo_reference['full_sequence']
        
        # Create PairwiseAligner 
        if not BIOPYTHON_AVAILABLE or PairwiseAligner is None:
            return {'error': 'BioPython not available'}

        aligner = PairwiseAligner()  # type: ignore
        # Use local alignment to find the best matching subsequence
        aligner.mode = 'local' 
        aligner.match_score = self.local_alignment_params['match_score']
        aligner.mismatch_score = self.local_alignment_params['mismatch_score']
        aligner.open_gap_score = self.local_alignment_params['gap_open_penalty']
        aligner.extend_gap_score = self.local_alignment_params['gap_extend_penalty']
        
        # Helper to run alignment
        def run_single_alignment(seq, orient):
            alignments = aligner.align(seq, ref_sequence)
            try:
                if len(alignments) == 0:
                    return None
            except OverflowError:
                 return None # Treat overflow as no good alignment for now check
            return alignments[0]

        # 1. Forward Alignment
        best_fwd = run_single_alignment(query_sequence, 'forward')
        
        # 2. Reverse Complement Alignment
        try:
            rev_seq_obj = Seq(query_sequence).reverse_complement()
            rev_sequence = str(rev_seq_obj)
            best_rev = run_single_alignment(rev_sequence, 'reverse')
        except Exception as e:
            print(f"Warning: Failed to generate reverse complement: {e}")
            best_rev = None

        # Compare and select best
        final_alignment = None
        orientation = 'forward'
        aligned_seq_content = query_sequence # Sequence used for alignment

        score_fwd = best_fwd.score if best_fwd else -float('inf')
        score_rev = best_rev.score if best_rev else -float('inf')

        if best_rev and score_rev > score_fwd:
            final_alignment = best_rev
            orientation = 'reverse'
            aligned_seq_content = rev_sequence
        elif best_fwd:
            final_alignment = best_fwd
            orientation = 'forward'
            aligned_seq_content = query_sequence
        else:
             return {'error': 'No significant alignment found in either orientation'}

        score = final_alignment.score
        
        # Extract coordinates and sequences
        query_start = final_alignment.coordinates[0][0]
        query_end = final_alignment.coordinates[0][-1]
        ref_start = final_alignment.coordinates[1][0]
        ref_end = final_alignment.coordinates[1][-1]
        
        aligned_ref_seq, aligned_query_seq = self._extract_aligned_sequences(final_alignment)
        
        matches = sum(1 for q, r in zip(aligned_query_seq, aligned_ref_seq)
                    if q == r and q != '-' and r != '-')
        alignment_length = len(aligned_query_seq)
        similarity = matches / alignment_length if alignment_length > 0 else 0
        
        return {
            'alignment_score': score,
            'similarity': similarity,
            'orientation': orientation,
            'query_start': query_start,
            'query_end': query_end,
            'ref_start': ref_start,
            'ref_end': ref_end,
            'aligned_query': aligned_query_seq,
            'aligned_reference': aligned_ref_seq,
            'genomic_region': {
                'start': ref_start,
                'end': ref_end
            },
            'matches': matches,
            'alignment_length': alignment_length
        }

    
    def find_best_exon_match(self, query_sequence: str) -> Dict:
        """
        Align a query sequence against ALL exons to find the best match.
        Checks both forward and reverse complement orientations.

        Args:
            query_sequence: The sequence to align

        Returns:
            Dictionary with best exon match details
        """
        if not self.abo_reference or 'exons' not in self.abo_reference:
            return {'error': 'Exon data not available'}

        best_result = None
        highest_score = -float('inf')

        # Create PairwiseAligner (Local alignment)
        if not BIOPYTHON_AVAILABLE or PairwiseAligner is None:
            return {'error': 'BioPython not available'}

        aligner = PairwiseAligner() # type: ignore
        aligner.mode = 'local'
        aligner.match_score = self.local_alignment_params['match_score']
        aligner.mismatch_score = self.local_alignment_params['mismatch_score']
        aligner.open_gap_score = self.local_alignment_params['gap_open_penalty']
        aligner.extend_gap_score = self.local_alignment_params['gap_extend_penalty']

        # Pre-calculate reverse complement of query
        try:
            query_fwd = query_sequence
            query_rev = str(Seq(query_sequence).reverse_complement())
        except Exception as e:
            return {'error': f"Failed to generate reverse complement: {e}"}

        for exon in self.abo_reference['exons']:
            exon_num = exon['exon_number']
            exon_seq = exon['sequence']

            # Helper to check one orientation
            def check_orient(q_seq, orient_name):
                # Use simple alignment first to get score
                # Note: aligner.score() is faster than aligner.align() if we just want the score
                score = aligner.score(q_seq, exon_seq)
                return score

            # Check Forward
            score_fwd = check_orient(query_fwd, 'forward')
            
            # Check Reverse
            score_rev = check_orient(query_rev, 'reverse')

            # Determine best for this exon
            current_best_score = max(score_fwd, score_rev)
            current_orientation = 'forward' if score_fwd >= score_rev else 'reverse'
            current_seq = query_fwd if current_orientation == 'forward' else query_rev

            if current_best_score > highest_score:
                highest_score = current_best_score
                
                # generate full alignment object for the winner to get details
                alignments = aligner.align(current_seq, exon_seq)
                best_aln = alignments[0]
                
                aligned_ref_seq, aligned_query_seq = self._extract_aligned_sequences(best_aln)
                alignment_length = len(aligned_query_seq)
                matches = sum(1 for q, r in zip(aligned_query_seq, aligned_ref_seq) if q == r and q != '-' and r != '-')
                similarity = matches / alignment_length if alignment_length > 0 else 0

                best_result = {
                    'exon_number': exon_num,
                    'orientation': current_orientation,
                    'alignment_score': current_best_score,
                    'similarity': similarity,
                    'matches': matches,
                    'length': alignment_length,
                    'query_start': best_aln.coordinates[0][0],
                    'query_end': best_aln.coordinates[0][-1],
                    'ref_start': best_aln.coordinates[1][0], # Exon-relative
                    'ref_end': best_aln.coordinates[1][-1],   # Exon-relative
                    'aligned_query': aligned_query_seq,
                    'aligned_reference': aligned_ref_seq
                }

    def get_exon_coverage_stats(self, query_sequence: str) -> List[Dict]:
        """
        Calculate coverage AND similarity statistics for ALL exons.

        Args:
            query_sequence: The sequence to analyze

        Returns:
            List of dictionaries containing stats for each exon
        """
        stats = []
        if not self.abo_reference or 'exons' not in self.abo_reference:
            return stats

        # Create PairwiseAligner
        if not BIOPYTHON_AVAILABLE or PairwiseAligner is None:
            return stats

        aligner = PairwiseAligner() # type: ignore
        aligner.mode = 'local'
        aligner.match_score = self.local_alignment_params['match_score']
        aligner.mismatch_score = self.local_alignment_params['mismatch_score']
        aligner.open_gap_score = self.local_alignment_params['gap_open_penalty']
        aligner.extend_gap_score = self.local_alignment_params['gap_extend_penalty']

        # Pre-calculate reverse complement
        try:
            query_fwd = query_sequence
            query_rev = str(Seq(query_sequence).reverse_complement())
        except Exception:
            return stats

        for exon in self.abo_reference['exons']:
            exon_num = exon['exon_number']
            exon_seq = exon['sequence']
            exon_len = len(exon_seq)

            # Helper to get best alignment for this exon
            def get_best_alignment(q_seq):
                alignments = aligner.align(q_seq, exon_seq)
                try:
                    if not alignments: return None
                    return alignments[0]
                except OverflowError: return None

            aln_fwd = get_best_alignment(query_fwd)
            aln_rev = get_best_alignment(query_rev)

            score_fwd = aln_fwd.score if aln_fwd else -float('inf')
            score_rev = aln_rev.score if aln_rev else -float('inf')

            best_aln = None
            orientation = 'F'
            
            if score_rev > score_fwd:
                best_aln = aln_rev
                orientation = 'R'
            else:
                best_aln = aln_fwd
                orientation = 'F'

            coverage_pct = 0.0
            similarity_pct = 0.0
            
            if best_aln:
                # Calculate reference coverage
                ref_start = best_aln.coordinates[1][0]
                ref_end = best_aln.coordinates[1][-1]
                covered_len = ref_end - ref_start
                coverage_pct = (covered_len / exon_len) * 100
                
                # Calculate Similarity (Identity)
                # Note: BioPython alignment object handling
                aligned_ref_seq, aligned_query_seq = self._extract_aligned_sequences(best_aln)
                aln_len = len(aligned_query_seq)
                matches = sum(1 for q, r in zip(aligned_query_seq, aligned_ref_seq) if q == r and q != '-' and r != '-')
                similarity_pct = (matches / aln_len * 100) if aln_len > 0 else 0.0

                # Sanity check: if scores are exceedingly low, assume noise
                if best_aln.score < 20: 
                    coverage_pct = 0.0
                    similarity_pct = 0.0

            stats.append({
                'exon': exon_num,
                'coverage': coverage_pct,
                'similarity': similarity_pct,
                'orientation': orientation,
                'score': score_fwd if orientation == 'F' else score_rev
            })
            
        return stats

    def identify_variants(self, query_sequence: str) -> Dict:
        """
        Identify the best matching exon and call variants.
        Uses 2-step verification (Cov > 80%, Sim > 90%) to confirm exon first.
        
        Args:
            query_sequence: The sequence to analyze
            
        Returns:
            Dictionary with exon, orientation, and variants
        """
        if not query_sequence:
            return {'error': 'Empty sequence'}
            
        stats = self.get_exon_coverage_stats(query_sequence)
        
        # Filter for confirmed matches
        candidates = []
        for s in stats:
            if s['coverage'] >= 80.0 and s['similarity'] >= 90.0:
                candidates.append(s)
        
        if not candidates:
            return {'error': 'No exon confirmed (Requires >80% Cov, >90% Match)'}
        
        # Pick best candidate (highest score)
        best = max(candidates, key=lambda x: x['score'])
        
        exon_num = best['exon']
        orientation = best['orientation']
        
        # Prepare sequence for variant calling
        # align_sequence_to_exon expects sequence in Forward orientation relative to exon
        target_seq = query_sequence
        if orientation == 'R':
            try:
                target_seq = str(Seq(query_sequence).reverse_complement())
            except:
                return {'error': 'Failed to reverse complement'}
                
        # Call Variants
        # Note: align_sequence_to_exon performs alignment again, which is redundant 
        # but ensures we get the detailed variant breakdown structure.
        result = self.align_sequence_to_exon(target_seq, exon_num)
        
        if 'error' in result:
            return result
            
        result['orientation'] = orientation
        return result

    def generate_batch_summary(self, fasta_files: List[Any]) -> List[Dict]:
        """
        Process a batch of FASTA files (or Streamlit UploadedFiles) and generate a summary.
        
        Args:
            fasta_files: List of file paths (str) OR file-like objects (Streamlit UploadedFile)
            
        Returns:
            List of dictionaries containing:
            - filename
            - exon (if confirmed)
            - coverage
            - correctness (similarity)
            - decision (Confirmed/No Match)
            - ... details for visualization
        """
        summary_results = []
        
        for file_obj in fasta_files:
            filename = "Unknown"
            sequence = ""
            
            # Handle different input types
            try:
                if isinstance(file_obj, str): # File path
                    filename = os.path.basename(file_obj)
                    with open(file_obj, 'r') as f:
                        lines = f.readlines()
                        sequence = "".join([l.strip() for l in lines if not l.startswith(">")])
                else: # Streamlit UploadedFile or IO
                    filename = getattr(file_obj, "name", "Uploaded File")
                    # Read content
                    content = file_obj.getvalue()
                    if isinstance(content, bytes):
                        content = content.decode('utf-8')
                    # Parse simple FASTA
                    lines = content.splitlines()
                    sequence = "".join([l.strip() for l in lines if not l.startswith(">")])

                if not sequence:
                    summary_results.append({
                        'filename': filename,
                        'decision': "Empty Sequence",
                        'exon': "-",
                        'coverage': 0.0,
                        'similarity': 0.0,
                        'variants': []
                    })
                    continue

                # Run Robust Identification
                result = self.identify_variants(sequence)
                
                # Format Result for Table
                entry = {
                    'filename': filename,
                    'variants': []
                }
                
                if 'error' in result:
                    entry['decision'] = "No Match / Low Quality"
                    entry['exon'] = "-"
                    entry['coverage'] = 0.0
                    entry['similarity'] = 0.0
                else:
                    exon_num = result['exon_number']
                    # We need to retrieve the stats that led to this decision
                    # identify_variants calls get_exon_coverage_stats internally but doesn't return them all.
                    # However, since identify_variants SUCCEEDED, we know it passed the threshold.
                    # To get the exact numbers, we could re-run get_exon_coverage_stats or trust the internal check.
                    # For the table, let's re-fetch the exact stats for THIS exon.
                    
                    # Optimization: identify_variants essentially picked the best one.
                    # We can assume high confidence. But let's get the numbers for display.
                    stats = self.get_exon_coverage_stats(sequence)
                    best_stat = next((s for s in stats if s['exon'] == exon_num), None)
                    
                    cov = best_stat['coverage'] if best_stat else 0.0
                    sim = best_stat['similarity'] if best_stat else 0.0
                    
                    entry['decision'] = f"Exon {exon_num} Confirmed"
                    entry['exon'] = str(exon_num)
                    entry['coverage'] = cov
                    entry['similarity'] = sim
                    entry['variants'] = result.get('variants', [])
                    
                    # Merge full result data (like aligned sequences) into entry
                    entry.update(result)

                summary_results.append(entry)

            except Exception as e:
                summary_results.append({
                    'filename': filename,
                    'decision': f"Error: {str(e)}",
                    'exon': "-",
                    'coverage': 0.0,
                    'similarity': 0.0,
                    'variants': []
                })
                
        return summary_results

# Example usage and testing functions
 

# if __name__ == "__main__":
#     if BIOPYTHON_AVAILABLE:
#         test_proper_alignment_service()
#     else:
#         print("BioPython not available. Please install with: pip install biopython")
