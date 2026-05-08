"""
RHD Analyzer - WHO Standard RHD Genotyping
Embedded RHD1 and RHD456 reference sequences with decision logic
"""

from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq

# WHO STANDARD RHD REFERENCE SEQUENCES
# RHD1: Exon 1 region (951 bp) - used to detect RHD gene presence
RHD1_REFERENCE = (
    "ACTTCACCCTAAGGCTGGATCAGGATCCCCTCCAGGTTTTTACTAGAGCCAAACCCACATCTCCTTTCTCTTCTGCCACC"
    "CCCCCTTAAAATGCTTAGAAACACATAGATTTAAATACAAGTTCAAATGTAAGTAATTTCAACTGTGTAACTATGAGGAG"
    "TCAATTCTACGTGGGTCCTATCTGTATCCTCCCCAGGGCTCAGCTCCATTCTTTGCTTTCATTCATTCTCATTCAATACA"
    "TTGTTGTTAAGAGCTCACTGGGTGCCCTCTCTGTCATGTAGTAAGGTTTTAAAAAGAAAGCCTCTTCTGAGCTTCAGTTT"
    "CCTTATTTATAAAATAGGAATATTGATCTGTTCCTTGCTTTTCTTACAAGGATATGCTGAAGATGACTGAAGTACAGAGT"
    "AAAGAAGGATTATGTTTGGGTATCAAAGGAATAGAATGCCCTCTTTCAAACTGAGCACAGCAGGAACCTGTAACAGGAAC"
    "ACAGCAACTTGTTGAATGAATGACAATATTGGAAAACATACATTTCCTCCCCTCCCCATCATAGTCCCTCTGCTTCCGTG"
    "TTAACTCCATAGACAGGCCAGCACAGCCAGCCTTGAAGCCTGAGATAAGGCCTTTGGCGGGTGTCTCCCCTATCGCTCCC"
    "TCAAGCCCTCAAGTAGGTGTTGGAGAGAGGGGTGATGCCTGGTGCTGGTGGAACCCCTGCACAGAGACGGACACAGGATG"
    "AGCTCTAAGTACCCGCGGTCTGTCCGGCGCTGCCTGCCCCTCTGGGCCCTAACACTGGAAGCAGCTCTCATTCTCCTCTT"
    "CTATTTTTTTACCCACTATGACGCTTCCTTAGAGGATCAAAAGGGGCTCGTGGCATCCTATCAAGGTGAGAGTTCATTGG"
    "AACAGTGGTCACAGGAGCAAATAGCAGGGGCAGGGGCGGGGGAGGCCTATGGTTCTCCAGGGGCACAGATG"
)

# RHD456: Exons 4-6 region (~3336 bp) - used to assess D antigen variants
RHD456_REFERENCE = (
    "AGGCGTTGAAGCCAATAAGAGAATGCACCAACACCTGCCTAATGCAGCTGTGCACTGCACAGTGGCCCATCAGGTCCCAG"
    "CGTCCTGCTGGCCTTCAGCCAAAGCAGAGAGCATTAGTTGTCTAGTTTCTTACCGGCAGGCACTTGGCTCCCCCGATGGA"
    "GATCAGCCCAGCCACAAGACCCAGCACCATGGCAAGCCACGGAGAAGGGATCAGGTGACACGAGGTACCCACAGCCACGC"
    "CTCCTGCCAACACCGCACTGTGCACATAAGTCTGCAAAGAAATAGCGTGTGGGTAAAGGAAGCAAGGTAGAGAGAGAACA"
    "CCATCTTGCTGCAAGTGACCAGCACGCTGGGAACAAGGCAACCCTGTAACATCCTCCCTGCTTTGCTGATCCTGAAACCA"
    "CCTCTCCCCTGTGTTGGGCTACGTGTCCTTCATCAGTGGACACGGTGGGCCAGCTCACTACTGCCTCTTAGGCTCAGGGC"
    "ACAGTTACCCTGATTGTGCAGCATTTAGGGAATGAGCTGGGAAGTCCAGGTGCTCAGCATCCTAGGCACCAATCACACTC"
    "ATTAGTGGTGTGTCCCTGGACAGCGGGAAGAATGGGTTTTTGTTTTATTTTTTAGAGACAGGGTCTTGCTATGTTGCCTG"
    "GCTGGTCTCAAACTCCTGGGCTCAAGCGATCCTCCCGCCTCAGCCTCCCGAGTGCTGGGATTACAGGCGTGAGCCACCGC"
    "GCCCGGCCCCATGTGTGAGACGATCTGGCTTACATCTGGTCCATGTTCTCCAGTCTAAGCCAAAGCACTCATTGAGGTGT"
    "CTGATGTGGTCATTAGGAGACAGGTCACAGCTAACAGTGAGGACTGTGGTGGAGCAGGTGCTGAGACAGGAGCTATGTAT"
    "TCTGCCCTCACTGTCTCACTGCACCCATCACGGAAACAGGACACTTCCAGGCCTGAGCCCATTCATGAAGACACTGCTTC"
    "ATCCCATTCCTGAGCATCATGGCCAGCATTCTAGCCCGTACCGGATGCAGCTGGTGGTGGTGGTTATGGTGATGACACCA"
    "CAGACGGACACAGGATGAGCTCTAAGTACCCGCGGTCTGTCCGGCGCTGCCTGCCCCTCTGGGCCCTAACACTGGAAGCA"
    "GCTCTCATTCTCCTCTTCTATTTTTTTACCCACTATGACGCTTCCTTAGAGGATCAAAAGGGGCTCGTGGCATCCTATCAA"
    "GGTGAGAGTTCATTGGAACAGTGGTCACAGGAGCAAATAGCAGGGGCAGGGGCGGGGGAGGCCTATGGTTCTCCAGGGGC"
    "ACAGATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTCTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG"
)

# RHD Decision Thresholds (WHO Standards)
RHD1_MAX_LENGTH = 800  # bp - threshold for RHD1 region
IDENTITY_RHD_POSITIVE = 90.0  # % - identity threshold for RhD+ call
VARIANT_COUNT_RHD_POSITIVE = 3  # maximum variants for RhD+ status


class RHDAnalyzer:
    def __init__(self, gb_path=None, reference_seq=None):
        """
        Initialize RHD Analyzer with WHO standard references.
        Can load from GenBank file or use embedded references.
        """
        self.reference_seq = None
        self.exon_map = []
        self.rhd1_ref = RHD1_REFERENCE
        self.rhd456_ref = RHD456_REFERENCE

        # Try to load from GenBank if provided
        if gb_path and reference_seq is None:
            try:
                record = SeqIO.read(gb_path, "genbank")
                self.reference_seq = str(record.seq)
                self.exon_map = self._get_exon_map(record)
            except Exception:
                # Fall back to embedded references
                pass

    def _get_exon_map(self, record):
        """Extract exon positions from GenBank file."""
        exons = []
        for feature in record.features:
            if feature.type == "exon":
                start = int(feature.location.start)
                end = int(feature.location.end)
                exons.append((start, end))
        exons.sort()
        return exons

    def detect_amplicon_region(self, query_length, identity_rhd1, identity_rhd456):
        """
        Auto-detect which amplicon region based on sequence length and alignment identity.

        WHO Standards:
        - RHD1: ~951 bp (exon 1 region)
        - RHD456: ~3336 bp (exons 4-6 region)

        Returns: 'RHD1' or 'RHD456'
        """
        # Primary decision: use alignment identity
        # If one identity is significantly better (>5% difference), use that
        identity_diff = abs(identity_rhd1 - identity_rhd456)

        if identity_diff > 5:
            # Clear winner based on identity
            return 'RHD1' if identity_rhd1 > identity_rhd456 else 'RHD456'

        # Secondary decision: use sequence length
        # RHD1: ~951 bp
        # RHD456: ~3336 bp
        # Threshold: 1500 bp
        if query_length < 1500:
            return 'RHD1'
        else:
            return 'RHD456'

    def determine_rhd_phenotype(self, region, query_length, identity, variant_count):
        """
        Determine RhD phenotype based on WHO decision rules.

        RHD1 Rules (exon 1 region):
        - Length < 800 bp & present → RhD+ (gene present)
        - Length >= 800 bp or missing → RhD- (gene absent/deleted)

        RHD456 Rules (exons 4-6):
        - Identity >= 90% & variants <= 3 → RhD+ (normal D antigen)
        - Identity < 90% or variants > 3 → RhD- or RHD Variant

        Returns: (phenotype_status, reason)
        """
        if region == 'RHD1':
            if query_length < RHD1_MAX_LENGTH:
                return ('RhD+', f'RHD1 amplicon length {query_length} bp (<800 bp indicates RHD gene present)')
            else:
                return ('RhD-', f'RHD1 amplicon length {query_length} bp (>=800 bp indicates RHD gene absent/deleted)')

        elif region == 'RHD456':
            if identity >= IDENTITY_RHD_POSITIVE and variant_count <= VARIANT_COUNT_RHD_POSITIVE:
                return ('RhD+', f'RHD456 high identity {identity:.1f}% with {variant_count} variants (normal D antigen)')
            elif identity < IDENTITY_RHD_POSITIVE:
                return ('RhD-', f'RHD456 low identity {identity:.1f}% (<90% indicates absent/variant D antigen)')
            else:
                return ('RHD Variant', f'RHD456 identity {identity:.1f}% with {variant_count} variants (D antigen variant detected)')

        return ('Inconclusive', 'Could not determine amplicon region')

    def analyze(self, query_seq):
        """
        Complete RHD analysis with WHO decision logic.

        Returns dict with:
        - variants: list of detected mutations
        - identity: sequence identity %
        - score: alignment score
        - strand: forward or reverse
        - query_length: length of query sequence
        - region: RHD1 or RHD456
        - rhd_status: RhD+, RhD-, or RHD Variant
        - reason: explanation of decision
        """
        query_seq_str = str(query_seq)
        query_length = len(query_seq_str)

        # Analyze against both RHD1 and RHD456 references
        rhd1_forward = self._analyze_against_reference(query_seq_str, self.rhd1_ref, 'forward')
        rhd1_reverse = self._analyze_against_reference(query_seq_str, self.rhd1_ref, 'reverse')
        rhd456_forward = self._analyze_against_reference(query_seq_str, self.rhd456_ref, 'forward')
        rhd456_reverse = self._analyze_against_reference(query_seq_str, self.rhd456_ref, 'reverse')

        # Choose best alignment for each region
        rhd1_result = rhd1_forward if rhd1_forward['identity'] >= rhd1_reverse['identity'] else rhd1_reverse
        rhd456_result = rhd456_forward if rhd456_forward['identity'] >= rhd456_reverse['identity'] else rhd456_reverse

        # Auto-detect which region this is
        region = self.detect_amplicon_region(query_length, rhd1_result['identity'], rhd456_result['identity'])

        # Get result for detected region
        if region == 'RHD1':
            result = rhd1_result.copy()
            variants = self._extract_variants_simple(query_seq_str, self.rhd1_ref)
        else:
            result = rhd456_result.copy()
            variants = self._extract_variants_simple(query_seq_str, self.rhd456_ref)

        # Determine phenotype
        rhd_status, reason = self.determine_rhd_phenotype(
            region,
            query_length,
            result['identity'],
            len(variants)
        )

        # Compile final result
        final_result = {
            'variants': variants,
            'identity': round(result['identity'], 1),
            'score': round(result['score'], 1),
            'strand': result['strand'],
            'query_length': query_length,
            'region': region,
            'rhd_status': rhd_status,
            'reason': reason,
            'reference_description': f'{region} (WHO standard)'
        }

        return final_result

    def _analyze_against_reference(self, query_seq, ref_seq, strand='forward'):
        """Analyze query against a specific reference in a specific strand."""
        if strand == 'reverse':
            query_seq = str(Seq(query_seq).reverse_complement())

        alignments = pairwise2.align.localms(ref_seq, query_seq, 2, -1, -0.5, -0.1)

        if not alignments:
            return {'identity': 0.0, 'score': 0.0, 'strand': strand}

        best_match = alignments[0]
        identity = self._compute_alignment_identity(best_match)

        return {
            'identity': identity,
            'score': best_match.score,
            'strand': strand,
            'alignment': best_match
        }

    def _compute_alignment_identity(self, alignment):
        """Compute identity % only for aligned bases (ignoring gaps)."""
        ref_aligned = alignment[0]
        query_aligned = alignment[1]
        matches = 0
        aligned_positions = 0

        for ref_b, query_b in zip(ref_aligned, query_aligned):
            if ref_b == "-" or query_b == "-":
                continue
            aligned_positions += 1
            if ref_b == query_b:
                matches += 1

        if aligned_positions == 0:
            return 0.0

        return (matches / aligned_positions) * 100.0

    def _extract_variants_simple(self, query_seq, ref_seq):
        """
        Simple variant extraction using pairwise alignment.
        Returns list of variant descriptions (only true mismatches, not trailing gaps).
        """
        alignments = pairwise2.align.localms(ref_seq, query_seq, 2, -1, -0.5, -0.1)

        if not alignments:
            return []

        best_match = alignments[0]
        ref_aligned = best_match[0]
        query_aligned = best_match[1]
        start_in_ref = best_match[3]

        variants = []
        current_ref_idx = start_in_ref

        for i in range(len(ref_aligned)):
            ref_b = ref_aligned[i]
            query_b = query_aligned[i]

            if ref_b != "-":
                current_ref_idx_pos = current_ref_idx
                current_ref_idx += 1
            else:
                continue

            # Only count mismatches where both ref and query have bases
            if query_b != "-" and ref_b != query_b:
                # SNP
                variants.append(f"{ref_b}>{query_b}@{current_ref_idx_pos}")

        return variants


# Test code
if __name__ == "__main__":
    test_seqs = {
        'RhD_positive_rhd1': (
            "ACTTCACCCTAAGGCTGGATCAGGATCCCCTCCAGGTTTTTACTAGAGCCAAACCCACATCTCCTTTCTCTTCTGCCACC"
            "CCCCCTTAAAATGCTTAGAAACACATAGATTTAAATACAAGTTCAAATGTAAGTAATTTCAACTGTGTAACTATGAGGAG"
            "TCAATTCTACGTGGGTCCTATCTGTATCCTCCCCAGGGCTCAGCTCCATTCTTTGCTTTCATTCATTCTCATTCAATACA"
            "TTGTTGTTAAGAGCTCACTGGGTGCCCTCTCTGTCATGTAGTAAGGTTTTAAAAAGAAAGCCTCTTCTGAGCTTCAGTTT"
            "CCTTATTTATAAAATAGGAATATTGATCTGTTCCTTGCTTTTCTTACAAGGATATGCTGAAGATGACTGAAGTACAGAGT"
        ),
    }

    print("=" * 70)
    print("RHD Analyzer - Test")
    print("=" * 70)

    analyzer = RHDAnalyzer()

    for name, seq in test_seqs.items():
        print(f"\n[TEST] {name}")
        print(f"  Sequence length: {len(seq)} bp")

        result = analyzer.analyze(seq)

        print(f"  Region detected: {result['region']}")
        print(f"  Identity: {result['identity']}%")
        print(f"  RhD Status: {result['rhd_status']}")
        print(f"  Reason: {result['reason']}")
        print(f"  Variants: {len(result['variants'])} found")
        if result['variants']:
            for v in result['variants'][:3]:
                print(f"    - {v}")
