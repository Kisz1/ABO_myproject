"""
RHCE Analyzer - ISBT-Informed RHCE Analysis (C/c and E/e antigens)

References:
- ISBT 004 (RH system) v6.4
- NM_020485.7 (RHCE mRNA, cDNA reference)
- NG_009208.1 (RHCE genomic)

Diagnostic SNP scope (Asian-population focus, agreed in plan):
    Standard C/c + E/e:  c.48G>C, c.178A>C, c.203G>A, c.307C>T, c.676G>C
    Asian partial E:     c.1025T>C (Cat. EII), c.1226A>G (Cat. EIII)

----------------------------------------------------------------------------
NM_020485.7 INVERTED-PARITY EXCEPTION (audit note)
----------------------------------------------------------------------------
NM_020485.7 was cloned in the early 1990s from an RhCce individual, and the
deposited mRNA sequence therefore carries variant alleles at three positions
relative to the global HGVS consensus:

    c.178A>C   - NM_020485.7 has 'C' (C-antigen variant base)
    c.203G>A   - NM_020485.7 has 'A' (C-antigen variant base)
    c.1025T>C  - NM_020485.7 has 'C' (Category-EII-like variant base)

To keep using NM_020485.7 as the reference file (instead of switching to
NG_009208 or curating a synthetic c-allele consensus), we apply a per-SNP
Hard-coded Exception Mapping. For each of the three SNPs we:

  1. Swap ref_base / alt_base so the SNP's ref_base matches what the file
     actually carries at that position. The reference-base sanity check
     in _genotype_snp_at then passes naturally.
  2. Swap ref_call / alt_call (the antigen label that each base produces)
     so the biological antigen call stays correct: a query base equal to
     the swapped ref_base is interpreted as the C-antigen-carrying allele
     (or, for c.1025, the partial-E variant), not as wild-type.
  3. Flag the SNP with `parity_inverted: True`. The consensus voter
     (_call_antigen_axis and the Asian-partial-marker detector) reads
     this flag and FLIPS hom_ref <-> hom_alt before comparing the
     confirmatory marker's zygosity against the primary marker's. Without
     this flip, a pure c-allele patient would look "discordant" at c.178
     simply because A is now alt_base post-swap.

Net effect: every diagnostic SNP is callable against NM_020485.7, primary
markers (c.48, c.676) still drive the call, and confirmatory markers
contribute their evidence to the confidence score rather than being
auto-skipped with `ref_base mismatch`.

If you migrate to a clean c-allele reference, set parity_inverted=False on
these three SNPs and restore the original ref_base/alt_base + ref_call/alt_call.
"""

from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq


# ─── Reference file locations ────────────────────────────────────────────────
# User supplies one of these (see utils/data/RHCE_REFERENCE_README.md).
DATA_DIR = Path(__file__).resolve().parent / "data"
RHCE_GENBANK_PATH = DATA_DIR / "rhce_referance.gb"
RHCE_FASTA_PATH = DATA_DIR / "rhce_cdna_reference.fasta"


# ─── ISBT-defined diagnostic SNPs for C/c and E/e calling ────────────────────
# cDNA_position is 1-based per NM_020485 mRNA convention.
# When a cDNA reference is loaded, position N maps directly to ref_str[N-1].
# When a genomic reference (GenBank) is loaded, we translate cDNA -> genomic
# via the exon map (see RHCEAnalyzer._cdna_to_ref_index).
RHCE_DIAGNOSTIC_SNPS = {
    'c.48G>C': {
        'exon': 1,
        'cDNA_position': 48,
        'ref_base': 'G',
        'alt_base': 'C',
        'antigen_axis': 'C/c',
        'alt_call': 'C',          # alt allele expresses C antigen
        'ref_call': 'c',
        'role': 'primary',         # primary marker for C/c
        'significance': 'p.Ser16Cys - primary C/c discriminator',
        'reference': 'ISBT 004 v6.4',
    },
    'c.178A>C': {
        'exon': 2,
        'cDNA_position': 178,
        # NM_020485.7 carries 'C' here (the C-allele variant base). To pass
        # the reference-base sanity check we swap ref/alt and the antigen
        # calls; consensus voting flips zygosity via parity_inverted.
        'ref_base': 'C',
        'alt_base': 'A',
        'antigen_axis': 'C/c',
        'ref_call': 'C',
        'alt_call': 'c',
        'parity_inverted': True,
        'nm020485_note': "NM_020485.7 reference carries the C-antigen variant 'C' at this position (Cce-origin clone); parity inverted in the SNP table.",
        'role': 'confirmatory',
        'significance': 'C-haplotype confirmatory marker (parity inverted on NM_020485.7)',
        'reference': 'ISBT 004 v6.4',
    },
    'c.203G>A': {
        'exon': 2,
        'cDNA_position': 203,
        # NM_020485.7 carries 'A' here (C-allele variant base). Swap + parity flag.
        'ref_base': 'A',
        'alt_base': 'G',
        'antigen_axis': 'C/c',
        'ref_call': 'C',
        'alt_call': 'c',
        'parity_inverted': True,
        'nm020485_note': "NM_020485.7 reference carries the C-antigen variant 'A' at this position (Cce-origin clone); parity inverted in the SNP table.",
        'role': 'confirmatory',
        'significance': 'C-haplotype confirmatory marker (parity inverted on NM_020485.7)',
        'reference': 'ISBT 004 v6.4',
    },
    'c.307C>T': {
        'exon': 2,
        'cDNA_position': 307,
        'ref_base': 'C',
        'alt_base': 'T',
        'antigen_axis': 'C/c',
        'alt_call': 'C',
        'ref_call': 'c',
        'role': 'confirmatory',
        'significance': 'p.Leu103Ser - C-haplotype confirmatory marker',
        'reference': 'ISBT 004 v6.4',
    },
    'c.676G>C': {
        'exon': 5,
        'cDNA_position': 676,
        'ref_base': 'G',
        'alt_base': 'C',
        'antigen_axis': 'E/e',
        'alt_call': 'E',          # alt allele expresses E antigen
        'ref_call': 'e',
        'role': 'primary',
        'significance': 'p.Ala226Pro - primary E/e discriminator',
        'reference': 'ISBT 004 v6.4',
    },
    'c.1025T>C': {
        'exon': 7,
        'cDNA_position': 1025,
        # NM_020485.7 carries 'C' here (Category-EII-like variant base). Swap + parity flag.
        'ref_base': 'C',
        'alt_base': 'T',
        'antigen_axis': 'E/e',
        'ref_call': 'E-partial',
        'alt_call': 'e',
        'parity_inverted': True,
        'nm020485_note': "NM_020485.7 reference carries the partial-E variant 'C' at this position (Cce-origin clone); parity inverted in the SNP table.",
        'role': 'asian_partial',
        'significance': 'Category EII - Asian partial E variant (parity inverted on NM_020485.7)',
        'reference': 'ISBT 004 v6.4 (RHCE*02.EII)',
    },
    'c.1226A>G': {
        'exon': 8,
        'cDNA_position': 1226,
        'ref_base': 'A',
        'alt_base': 'G',
        'antigen_axis': 'E/e',
        'alt_call': 'E-partial',
        'ref_call': 'e',
        'role': 'asian_partial',
        'significance': 'Category EIII - Asian partial E variant',
        'reference': 'ISBT 004 v6.4 (RHCE*03.EIII)',
    },
}


# ─── IUPAC decode (same convention as rhd_analyzer) ──────────────────────────
IUPAC_DECODE = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC',
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG',
    'N': 'ACGT',
}


# Whole-read identity threshold. Kept for reporting (`callable_read` flag and
# `reads_callable` metric), but no longer gates SNP calling — see below.
MIN_IDENTITY_FOR_CALLING = 85.0

# Per-SNP local-window quality gate. The whole-read threshold is too coarse
# when a genomic Sanger amplicon is aligned against the cDNA reference: the
# embedded introns drag global identity down (e.g. 65%) even though the local
# alignment around each diagnostic SNP (which lies in an exon) is essentially
# perfect. We therefore evaluate identity in a window centered on each SNP
# and decide trust per-SNP rather than per-read.
SNP_LOCAL_WINDOW_HALF = 50          # ± columns around the SNP alignment column
SNP_LOCAL_WINDOW_MIN_BASES = 30     # minimum aligned (non-gap) bases required
MIN_LOCAL_IDENTITY_FOR_CALLING = 85.0  # window identity threshold to trust a SNP call

# Probe-based SNP calling. For each diagnostic SNP we extract a small
# ±SNP_PROBE_HALF bp window of the reference centered on the SNP and locally
# align JUST that probe against the read. This is dramatically more robust
# than reading the SNP off the whole-read alignment: a genomic Sanger
# amplicon that contains introns can be "smeared" across the entire cDNA
# reference by the whole-read aligner, putting intronic content under the
# SNP column. A small exonic probe matches the corresponding exonic region
# in the read at near-perfect identity regardless of surrounding introns.
SNP_PROBE_HALF = 50               # bp of reference flanking the SNP to use as probe
SNP_PROBE_MIN_BASES = 60          # minimum aligned (non-gap) bases for a usable probe hit
MIN_PHRED_AT_SNP = 30             # per-base Phred Q-score required AT the SNP column
                                  # when AB1 quality is available. Q30 = 1-in-1000
                                  # base error rate (lab standard).


class RHCEReferenceMissingError(RuntimeError):
    """Raised when no RHCE reference file is available."""


class RHCEAnalyzer:
    """
    Single-read RHCE analyzer.

    Loads either a GenBank reference (preferred — includes exon map for
    cDNA->genomic translation) or a cDNA FASTA (positions map directly to
    string indices).

    Use ``analyze_single_read(query_seq)`` to get per-SNP calls for one read.
    Multi-read consensus voting is layered on top in step 3.
    """

    def __init__(
        self,
        gb_path: Optional[str] = None,
        fasta_path: Optional[str] = None,
        reference_seq: Optional[str] = None,
    ):
        self.reference_seq: Optional[str] = None
        self.reference_kind: Optional[str] = None  # 'cdna' or 'genomic'
        self.exon_map: list = []                   # list of (start_genomic_0based, end_genomic_exclusive, exon_number)
        self.cds_start_genomic: Optional[int] = None  # 0-based start of CDS in reference
        self.cdna_offset: int = 0                  # cDNA pos 1 -> reference index = (1-1) + cdna_offset
        self._aligner = self._build_aligner()

        if reference_seq is not None:
            self.reference_seq = str(reference_seq).upper()
            self.reference_kind = 'cdna'
            self.cdna_offset = self._detect_cds_start_in_mrna(self.reference_seq)
            return

        gb = Path(gb_path) if gb_path else RHCE_GENBANK_PATH
        fa = Path(fasta_path) if fasta_path else RHCE_FASTA_PATH

        if gb.exists():
            self._load_genbank(gb)
        elif fa.exists():
            self._load_fasta(fa)
        # else: stay unloaded — analyze_single_read will raise with guidance.

    # ─── Reference loading ───────────────────────────────────────────────
    def _load_genbank(self, path: Path) -> None:
        record = SeqIO.read(str(path), "genbank")
        self.reference_seq = str(record.seq).upper()
        self.reference_kind = 'genomic'

        exons = []
        cds_start = None
        for feature in record.features:
            if feature.type == "exon":
                start = int(feature.location.start)
                end = int(feature.location.end)
                number = feature.qualifiers.get("number", [None])[0]
                try:
                    number = int(number) if number is not None else None
                except (TypeError, ValueError):
                    number = None
                exons.append((start, end, number))
            elif feature.type == "CDS" and cds_start is None:
                cds_start = int(feature.location.start)

        exons.sort(key=lambda e: e[0])
        # Backfill exon numbers if GenBank didn't carry them.
        for i, (s, e, n) in enumerate(exons):
            if n is None:
                exons[i] = (s, e, i + 1)
        self.exon_map = exons
        self.cds_start_genomic = cds_start

    def _load_fasta(self, path: Path) -> None:
        record = next(SeqIO.parse(str(path), "fasta"))
        self.reference_seq = str(record.seq).upper()
        self.reference_kind = 'cdna'
        # mRNA FASTAs include a 5' UTR before the ATG. HGVS cDNA numbering
        # starts at the A of ATG, so we need an offset for correct indexing.
        self.cdna_offset = self._detect_cds_start_in_mrna(self.reference_seq)

    @staticmethod
    def _detect_cds_start_in_mrna(seq: str) -> int:
        """Return the 0-based mRNA index of the ATG start codon for RHCE.

        Looks for the RHCE-specific start-codon context first
        (ATGAGCTCTAAGTAC = M-S-S-K-Y) to avoid false ATGs in the 5' UTR;
        falls back to the first ATG. Returns 0 if no ATG is found, which
        leaves cDNA-1 indexing unchanged (correct for non-mRNA inputs).
        """
        rhce_context = seq.find('ATGAGCTCTAAGTAC')
        if rhce_context >= 0:
            return rhce_context
        first_atg = seq.find('ATG')
        return first_atg if first_atg >= 0 else 0

    def is_loaded(self) -> bool:
        return self.reference_seq is not None

    def _require_reference(self) -> None:
        if not self.is_loaded():
            raise RHCEReferenceMissingError(
                "RHCE reference not found. Supply one of:\n"
                f"  - GenBank: {RHCE_GENBANK_PATH}\n"
                f"  - FASTA:   {RHCE_FASTA_PATH}\n"
                "See utils/data/RHCE_REFERENCE_README.md for download instructions."
            )

    # ─── cDNA-position translation ────────────────────────────────────────
    def _cdna_to_ref_index(self, cdna_pos: int) -> Optional[int]:
        """
        Translate a 1-based cDNA position to a 0-based index into the
        reference string.

        cDNA reference: trivially cdna_pos - 1.
        Genomic reference: walk the exon map, accumulating CDS length until
        we cover cdna_pos.
        """
        if self.reference_kind == 'cdna':
            idx = (cdna_pos - 1) + self.cdna_offset
            if 0 <= idx < len(self.reference_seq):
                return idx
            return None

        # Genomic: walk exons in transcript order. cDNA position 1 starts at
        # the CDS start (not the 5' UTR), so we skip bases before CDS.
        remaining = cdna_pos
        for start, end, _ in self.exon_map:
            exon_len = end - start
            exon_cds_start = start
            exon_cds_len = exon_len

            # Trim 5' UTR from the first relevant exon.
            if self.cds_start_genomic is not None and start < self.cds_start_genomic < end:
                exon_cds_start = self.cds_start_genomic
                exon_cds_len = end - self.cds_start_genomic
            elif self.cds_start_genomic is not None and end <= self.cds_start_genomic:
                continue  # entirely 5' UTR exon

            if remaining <= exon_cds_len:
                return exon_cds_start + (remaining - 1)
            remaining -= exon_cds_len

        return None

    # ─── Alignment ────────────────────────────────────────────────────────
    @staticmethod
    def _build_aligner() -> PairwiseAligner:
        """Configure a local aligner matching the pairwise2 localms(2,-1,-5,-0.5) setup."""
        aligner = PairwiseAligner()
        aligner.mode = 'local'
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -0.5
        return aligner

    def _align(self, query_seq: str, strand: str = 'forward'):
        seq = query_seq if strand == 'forward' else str(Seq(query_seq).reverse_complement())
        try:
            alignments = self._aligner.align(self.reference_seq, seq)
            best = alignments[0]
        except (ValueError, IndexError, OverflowError):
            return None
        return {
            'identity': self._identity(best),
            'score': best.score,
            'strand': strand,
            'alignment': best,
            'aligned_query': seq,
        }

    @staticmethod
    def _identity(alignment) -> float:
        ref_a, query_a = str(alignment[0]), str(alignment[1])
        matches = aligned = 0
        for r, q in zip(ref_a, query_a):
            if r == '-' or q == '-':
                continue
            aligned += 1
            if r == q:
                matches += 1
        return (matches / aligned * 100.0) if aligned else 0.0

    @staticmethod
    def _ref_index_to_alignment_column(alignment, ref_index_0based: int) -> Optional[int]:
        """Find the alignment column corresponding to a given reference index.

        PairwiseAligner returns aligned strings that span the full local
        alignment region (with gaps inserted to show indels). Column N of the
        aligned reference corresponds to original reference index =
        (count of non-gap chars in aligned_ref[0:N+1]) - 1 + start_offset,
        where start_offset is the first reference position covered by the
        local alignment (obtained via alignment.coordinates[0][0]).

        Returns None if the SNP position is outside the aligned window or
        falls in a query-gap region of the alignment.
        """
        ref_aligned = str(alignment[0])
        query_aligned = str(alignment[1])
        # For local alignments, the aligned strings start at coordinates[0][0]
        # on the reference, not at index 0. Adjust the cursor accordingly.
        try:
            cursor = int(alignment.coordinates[0][0])
        except (AttributeError, IndexError, TypeError):
            cursor = 0
        for col, ref_char in enumerate(ref_aligned):
            if ref_char == '-':
                continue
            if cursor == ref_index_0based:
                # Treat as uncovered when the query has a gap at this column
                # (the SNP position is outside the query's aligned window).
                if query_aligned[col] == '-':
                    return None
                return col
            cursor += 1
        return None

    @staticmethod
    def _query_pos_at_alignment_column(
        alignment, col: int, query_len: int, strand: str
    ) -> Optional[int]:
        """Map alignment column -> 0-based position in the original query
        sequence (before reverse-complement). Used for Phred-at-SNP lookup."""
        query_aligned = str(alignment[1])
        try:
            q_cursor = int(alignment.coordinates[1][0])
        except (AttributeError, IndexError, TypeError):
            q_cursor = 0
        for c, q_char in enumerate(query_aligned):
            if c == col:
                if q_char == '-':
                    return None
                if strand == 'reverse':
                    return query_len - 1 - q_cursor
                return q_cursor
            if q_char != '-':
                q_cursor += 1
        return None

    # ─── SNP genotyping ──────────────────────────────────────────────────
    @staticmethod
    def _local_window_identity(alignment, center_col: int) -> tuple:
        """Identity % in a ±SNP_LOCAL_WINDOW_HALF column window around center_col.

        Gap columns are excluded from both numerator and denominator (same
        convention as the whole-read _identity helper). Returns the tuple
        (identity_percent, num_aligned_bases) so callers can also enforce a
        minimum-coverage requirement before trusting the percentage.
        """
        ref_aligned = str(alignment[0])
        query_aligned = str(alignment[1])
        n = len(ref_aligned)
        lo = max(0, center_col - SNP_LOCAL_WINDOW_HALF)
        hi = min(n, center_col + SNP_LOCAL_WINDOW_HALF + 1)
        matches = total = 0
        for r, q in zip(ref_aligned[lo:hi], query_aligned[lo:hi]):
            if r == '-' or q == '-':
                continue
            total += 1
            if r.upper() == q.upper():
                matches += 1
        if total == 0:
            return 0.0, 0
        return (matches / total * 100.0), total

    def _call_snp_from_column(
        self,
        alignment,
        col: int,
        snp_info: dict,
        window_pct: float,
        window_n: int,
        query_seq: Optional[str] = None,
        query_quality: Optional[list] = None,
        strand: Optional[str] = None,
    ) -> dict:
        """Convert an aligned (ref_base, query_base) pair at `col` into a SNP call dict.

        Assumes the local-window quality check has already passed. Shared by
        the whole-read alignment path (`_genotype_snp_at`) and the per-SNP
        probe path (`_genotype_snp_via_probe`).

        When ``query_quality`` is provided, additionally gates the call by
        the per-base Phred Q-score at the SNP column. Skipped for FASTA
        inputs (no quality data).
        """
        ref_aligned = str(alignment[0])
        query_aligned = str(alignment[1])
        ref_b_at_col = ref_aligned[col].upper()
        query_b = query_aligned[col].upper()

        meta = {
            'local_window_identity': round(window_pct, 1),
            'local_window_bases': window_n,
        }

        # Phred Q-score gate at the SNP column.
        if query_quality is not None and query_seq is not None and strand is not None:
            q_pos = self._query_pos_at_alignment_column(
                alignment, col, len(query_seq), strand
            )
            if q_pos is not None and 0 <= q_pos < len(query_quality):
                phred = int(query_quality[q_pos])
                meta['phred_at_snp'] = phred
                if phred < MIN_PHRED_AT_SNP:
                    return {
                        'covered': True,
                        'call':    'no_call',
                        'reason':  (f'low Phred at SNP position '
                                    f'(Q={phred} < Q{MIN_PHRED_AT_SNP})'),
                        **meta,
                    }

        if ref_b_at_col != snp_info['ref_base']:
            return {
                'covered': True,
                'call': 'no_call',
                'reason': f"ref_base mismatch: reference has {ref_b_at_col}, expected {snp_info['ref_base']}",
                **meta,
            }

        if query_b in ('-', 'N'):
            return {'covered': True, 'call': 'no_call', 'reason': 'gap or N in query', **meta}

        decoded = IUPAC_DECODE.get(query_b, query_b)
        ref_b = snp_info['ref_base']
        alt_b = snp_info['alt_base']

        if len(decoded) == 1 and decoded == ref_b:
            zygosity, call = 'hom_ref', snp_info['ref_call']
        elif len(decoded) == 1 and decoded == alt_b:
            zygosity, call = 'hom_alt', snp_info['alt_call']
        elif len(decoded) == 2 and ref_b in decoded and alt_b in decoded:
            zygosity = 'het'
            call = f"{snp_info['ref_call']}/{snp_info['alt_call']}"
        else:
            return {
                'covered': True,
                'call': 'no_call',
                'reason': f'unexpected base {query_b} (decoded={decoded})',
                **meta,
            }

        return {
            'covered': True,
            'call': call,
            'zygosity': zygosity,
            'query_base': query_b,
            'ref_base': ref_b,
            'alt_base': alt_b,
            'antigen_axis': snp_info['antigen_axis'],
            'role': snp_info['role'],
            'exon': snp_info['exon'],
            'cDNA_position': snp_info['cDNA_position'],
            'significance': snp_info['significance'],
            **meta,
        }

    def _genotype_snp_at(
        self,
        alignment,
        ref_index_0based: int,
        snp_info: dict,
        query_seq: Optional[str] = None,
        query_quality: Optional[list] = None,
        strand: Optional[str] = None,
    ) -> Optional[dict]:
        """Whole-read-alignment SNP genotyper, gated by the local-window check."""
        col = self._ref_index_to_alignment_column(alignment, ref_index_0based)
        if col is None:
            return None

        window_pct, window_n = self._local_window_identity(alignment, col)
        if window_n < SNP_LOCAL_WINDOW_MIN_BASES:
            return {
                'covered': True,
                'call': 'no_call',
                'reason': f'local window too short ({window_n} aligned bases)',
                'local_window_identity': round(window_pct, 1),
                'local_window_bases': window_n,
            }
        if window_pct < MIN_LOCAL_IDENTITY_FOR_CALLING:
            return {
                'covered': True,
                'call': 'no_call',
                'reason': (
                    f'local window identity {window_pct:.1f}% '
                    f'< {MIN_LOCAL_IDENTITY_FOR_CALLING}%'
                ),
                'local_window_identity': round(window_pct, 1),
                'local_window_bases': window_n,
            }

        return self._call_snp_from_column(
            alignment, col, snp_info, window_pct, window_n,
            query_seq=query_seq,
            query_quality=query_quality,
            strand=strand,
        )

    @staticmethod
    def _probe_alignment_identity(alignment) -> tuple:
        """Identity % + aligned-base count over a (small) probe alignment.

        Same gap-excluding convention as _local_window_identity; computed over
        the FULL aligned span since the probe is ~100 bp by construction.
        """
        ref_aligned = str(alignment[0])
        query_aligned = str(alignment[1])
        matches = total = 0
        for r, q in zip(ref_aligned, query_aligned):
            if r == '-' or q == '-':
                continue
            total += 1
            if r.upper() == q.upper():
                matches += 1
        if total == 0:
            return 0.0, 0
        return (matches / total * 100.0), total

    def _genotype_snp_via_probe(
        self,
        query_seq: str,
        snp_info: dict,
        query_quality: Optional[list] = None,
    ) -> dict:
        """Per-SNP genotyping by aligning a small reference probe against the read.

        The whole-read alignment can be unreliable for multi-exon genomic
        Sanger amplicons (introns get smeared across the cDNA reference and
        push intronic content under the SNP column). A small reference probe
        (±SNP_PROBE_HALF bp around the SNP) is entirely exonic, so it locates
        the corresponding exonic region in the read at near-perfect identity
        regardless of surrounding introns.

        ``query_quality`` enables the Phred-Q30-at-SNP-column gate.
        """
        ref_index = self._cdna_to_ref_index(snp_info['cDNA_position'])
        if ref_index is None:
            return {'covered': False, 'reason': 'cDNA position outside reference'}

        lo = max(0, ref_index - SNP_PROBE_HALF)
        hi = min(len(self.reference_seq), ref_index + SNP_PROBE_HALF + 1)
        probe = self.reference_seq[lo:hi]
        snp_offset_in_probe = ref_index - lo

        best = None
        for strand_name, q in (
            ('forward', query_seq),
            ('reverse', str(Seq(query_seq).reverse_complement())),
        ):
            try:
                aln = self._aligner.align(probe, q)[0]
            except (ValueError, IndexError, OverflowError):
                continue
            ident, n = self._probe_alignment_identity(aln)
            if best is None or ident > best['identity']:
                best = {'alignment': aln, 'identity': ident, 'n': n, 'strand': strand_name}

        if best is None:
            return {'covered': False, 'reason': 'probe produced no alignment'}
        if best['n'] < SNP_PROBE_MIN_BASES:
            return {
                'covered': False,
                'reason': f'probe alignment too short ({best["n"]} aligned bases)',
                'local_window_identity': round(best['identity'], 1),
                'local_window_bases': best['n'],
            }
        if best['identity'] < MIN_LOCAL_IDENTITY_FOR_CALLING:
            return {
                'covered': True,
                'call': 'no_call',
                'reason': f'probe identity {best["identity"]:.1f}% < {MIN_LOCAL_IDENTITY_FOR_CALLING}%',
                'local_window_identity': round(best['identity'], 1),
                'local_window_bases': best['n'],
            }

        col = self._ref_index_to_alignment_column(best['alignment'], snp_offset_in_probe)
        if col is None:
            return {
                'covered': False,
                'reason': 'SNP position outside probe alignment or in query gap',
                'local_window_identity': round(best['identity'], 1),
                'local_window_bases': best['n'],
            }

        result = self._call_snp_from_column(
            best['alignment'], col, snp_info, best['identity'], best['n'],
            query_seq=query_seq,
            query_quality=query_quality,
            strand=best['strand'],
        )
        result['probe_strand'] = best['strand']
        return result

    # ─── Public single-read pipeline ─────────────────────────────────────
    def analyze_single_read(
        self,
        query_seq: str,
        read_id: Optional[str] = None,
        query_quality: Optional[list] = None,
    ) -> dict:
        """
        Analyze one read against the RHCE reference.

        ``query_quality`` is optional per-base Phred Q-scores. When provided,
        each SNP call is additionally gated by Q-score at the SNP column
        (Q30 default — lab standard). FASTA inputs pass None and skip the
        gate (legacy 2-tuple behaviour).

        Returns a dict shaped to be consumed by the multi-read consensus
        stage (step 3) and the UI:

            {
                'read_id': str or None,
                'query_length': int,
                'strand': 'forward' | 'reverse',
                'identity': float,
                'score': float,
                'reference_kind': 'cdna' | 'genomic',
                'snp_calls': {
                    'c.48G>C': {'covered': True, 'call': 'C/c', 'zygosity': 'het', ...},
                    'c.676G>C': {'covered': False},
                    ...
                },
                'callable': bool,   # True if identity high enough to trust SNP calls
                'reason': str,
            }
        """
        self._require_reference()

        query_str = str(query_seq).upper()
        query_length = len(query_str)

        forward = self._align(query_str, 'forward')
        reverse = self._align(query_str, 'reverse')

        if forward is None and reverse is None:
            return {
                'read_id': read_id,
                'query_length': query_length,
                'strand': None,
                'identity': 0.0,
                'score': 0.0,
                'reference_kind': self.reference_kind,
                'snp_calls': {},
                'callable': False,
                'reason': 'No alignment produced',
            }

        best = max(
            (a for a in (forward, reverse) if a is not None),
            key=lambda a: a['identity'],
        )

        # `callable_read` is now informational only — it reflects whether the
        # whole-read alignment hit the legacy 85% identity bar (which is hard
        # to meet for genomic Sanger amplicons against the cDNA reference).
        # Per-SNP calling uses the probe-based method below, which is robust
        # to introns and primer flanks because each probe is small and exonic.
        callable_read = best['identity'] >= MIN_IDENTITY_FOR_CALLING
        snp_calls: dict = {
            snp_name: self._genotype_snp_via_probe(
                query_str, snp_info, query_quality=query_quality,
            )
            for snp_name, snp_info in RHCE_DIAGNOSTIC_SNPS.items()
        }

        # Count how many SNPs cleared the per-SNP local-window gate. This is
        # the meaningful "trustworthy calls in this read" figure.
        trusted_snps = sum(
            1 for c in snp_calls.values()
            if c.get('covered') and c.get('call') and c.get('call') != 'no_call'
        )
        low_phred_at_snp = sum(
            1 for c in snp_calls.values()
            if (c.get('phred_at_snp') is not None
                and c.get('phred_at_snp') < MIN_PHRED_AT_SNP)
        )

        reason = (
            f"Read aligned at {best['identity']:.1f}% identity ({best['strand']}); "
            f"{trusted_snps} diagnostic SNP(s) cleared the local-window quality gate"
            + (f"; {low_phred_at_snp} blocked by Phred Q{MIN_PHRED_AT_SNP} gate"
               if low_phred_at_snp else "")
            + "."
        )

        return {
            'read_id': read_id,
            'query_length': query_length,
            'strand': best['strand'],
            'identity': round(best['identity'], 2),
            'score': round(best['score'], 2),
            'reference_kind': self.reference_kind,
            'snp_calls': snp_calls,
            'callable': callable_read,
            'trusted_snps': trusted_snps,
            'low_phred_at_snp': low_phred_at_snp,
            'reason': reason,
        }

    # ─── Multi-read consensus (step 3) ───────────────────────────────────
    def analyze(self, reads) -> dict:
        """
        International-standard multi-read RHCE pipeline.

        Accepts:
          - a single sequence string (treated as a single read)
          - a list of sequence strings
          - a list of (read_id, sequence) tuples

        Returns a consensus result with per-SNP voting, antigen-axis calls
        (C/c and E/e), combined genotype, possible ISBT haplotype pairs,
        and overall confidence.
        """
        self._require_reference()

        normalized = self._normalize_reads(reads)
        per_read = [
            self.analyze_single_read(seq, read_id=rid, query_quality=qual)
            for rid, seq, qual in normalized
        ]

        snp_consensus = self._consensus_snp_calls(per_read)
        c_axis = self._call_antigen_axis(snp_consensus, 'C/c')
        e_axis = self._call_antigen_axis(snp_consensus, 'E/e')
        phenotype_block = self._phenotype_and_alleles(c_axis, e_axis, snp_consensus)

        callable_reads = [r for r in per_read if r['callable']]
        reads_with_trusted_snps = [r for r in per_read if r.get('trusted_snps', 0) > 0]
        no_calls = [
            snp for snp, cons in snp_consensus.items()
            if cons['consensus'] == 'no_call' and RHCE_DIAGNOSTIC_SNPS[snp]['role'] == 'primary'
        ]

        # Overall confidence = min across both axes. The per-SNP local-window
        # gate already encodes trustworthiness, so we no longer downgrade to
        # NONE just because the legacy whole-read 85% gate rejected every read.
        # If both axes are NONE (e.g. both primaries no_call), _min_confidence
        # naturally returns NONE.
        overall = _min_confidence(c_axis['confidence'], e_axis['confidence'])

        return {
            'phenotype': phenotype_block['phenotype'],
            'allele_options': phenotype_block['allele_options'],
            'c_e_call': c_axis,
            'big_E_call': e_axis,
            'snp_consensus': snp_consensus,
            'no_calls_primary': no_calls,
            'overall_confidence': overall,
            'reads_total': len(per_read),
            'reads_callable': len(callable_reads),
            'reads_with_trusted_snps': len(reads_with_trusted_snps),
            'per_read_details': per_read,
            'reason': phenotype_block['reason'],
        }

    @staticmethod
    def _normalize_reads(reads) -> list:
        """Coerce input into ``[(read_id, sequence, quality_or_None), ...]``.

        Accepts: a single seq str, a list of seq strs, a list of (rid, seq)
        2-tuples (FASTA-style, no quality), or a list of (rid, seq, quality)
        3-tuples (AB1-style, per-base Phred Q-scores). 2-tuples are widened
        with quality=None for backward compat.
        """
        if isinstance(reads, str):
            return [('read_1', reads, None)]
        if isinstance(reads, (list, tuple)):
            out = []
            for i, item in enumerate(reads, 1):
                if isinstance(item, str):
                    out.append((f'read_{i}', item, None))
                elif isinstance(item, (list, tuple)) and len(item) == 2:
                    rid, seq = item
                    out.append((str(rid) if rid is not None else f'read_{i}',
                                seq, None))
                elif isinstance(item, (list, tuple)) and len(item) == 3:
                    rid, seq, qual = item
                    out.append((str(rid) if rid is not None else f'read_{i}',
                                seq, qual))
                else:
                    raise ValueError(
                        f"Unsupported read entry at index {i-1}: expected str, "
                        f"(id, seq), or (id, seq, quality)"
                    )
            return out
        raise ValueError("reads must be a string, list of strings, or list of (id, seq[, quality]) tuples")

    @staticmethod
    def _consensus_snp_calls(per_read_results: list) -> dict:
        """
        Vote on each diagnostic SNP across all reads.

        Rules:
          - ≥2 reads, all concordant on zygosity  -> HIGH confidence
          - ≥2 reads, majority agrees             -> MEDIUM (flagged discordant)
          - exactly 1 read                        -> MEDIUM (flagged single-read)
          - tie or no clear majority              -> LOW
          - 0 reads cover the SNP                 -> no_call
        """
        consensus: dict = {}

        for snp_name in RHCE_DIAGNOSTIC_SNPS:
            covering = []
            for read in per_read_results:
                call = read['snp_calls'].get(snp_name)
                if not call or not call.get('covered'):
                    continue
                if call.get('call') == 'no_call':
                    continue
                covering.append({
                    'read_id': read['read_id'],
                    'zygosity': call['zygosity'],
                    'call': call['call'],
                    'query_base': call.get('query_base'),
                })

            if not covering:
                consensus[snp_name] = {
                    'consensus': 'no_call',
                    'confidence': 'NONE',
                    'reads_covering': 0,
                    'votes': {},
                    'supporting_reads': [],
                    'discordant_reads': [],
                }
                continue

            votes: dict = {}
            for c in covering:
                votes[c['zygosity']] = votes.get(c['zygosity'], 0) + 1

            top_zygosity, top_count = max(votes.items(), key=lambda kv: kv[1])
            # Tie detection
            tied = [z for z, n in votes.items() if n == top_count]
            n_reads = len(covering)

            supporting = [c for c in covering if c['zygosity'] == top_zygosity]
            discordant = [c for c in covering if c['zygosity'] != top_zygosity]
            consensus_call = supporting[0]['call']

            if len(tied) > 1:
                confidence = 'LOW'
            elif n_reads >= 2 and top_count == n_reads:
                confidence = 'HIGH'
            elif n_reads >= 2 and top_count > n_reads / 2:
                confidence = 'MEDIUM'
            elif n_reads == 1:
                confidence = 'MEDIUM'
            else:
                confidence = 'LOW'

            consensus[snp_name] = {
                'consensus': top_zygosity,        # hom_ref | het | hom_alt
                'call': consensus_call,            # e.g. 'c', 'C', 'c/C'
                'confidence': confidence,
                'reads_covering': n_reads,
                'votes': votes,
                'supporting_reads': [c['read_id'] for c in supporting],
                'discordant_reads': [c['read_id'] for c in discordant],
            }

        return consensus

    @staticmethod
    def _call_antigen_axis(snp_consensus: dict, axis: str) -> dict:
        """
        Combine SNPs on one axis (C/c or E/e) into a single genotype call.

        Strategy:
          - Use the primary SNP if it has a non-no_call consensus.
          - Confirmatory SNPs should agree on zygosity; downgrade confidence
            on disagreement.
          - For E/e: detect Asian partial-E markers (c.1025/c.1226) and
            append a partial-E flag without changing the underlying e call.
        """
        primary_snps = [
            name for name, info in RHCE_DIAGNOSTIC_SNPS.items()
            if info['antigen_axis'] == axis and info['role'] == 'primary'
        ]
        confirmatory_snps = [
            name for name, info in RHCE_DIAGNOSTIC_SNPS.items()
            if info['antigen_axis'] == axis and info['role'] == 'confirmatory'
        ]
        asian_partial_snps = [
            name for name, info in RHCE_DIAGNOSTIC_SNPS.items()
            if info['antigen_axis'] == axis and info['role'] == 'asian_partial'
        ]

        primary_name = primary_snps[0] if primary_snps else None
        primary = snp_consensus.get(primary_name) if primary_name else None

        if not primary or primary['consensus'] == 'no_call':
            return {
                'axis': axis,
                'genotype': None,
                'confidence': 'NONE',
                'primary_snp': primary_name,
                'supporting_snps': [],
                'discordant_snps': [],
                'partial_markers': [],
                'reason': f'Primary marker {primary_name} not callable',
            }

        # Translate primary zygosity to a genotype label per axis.
        if axis == 'C/c':
            genotype = {'hom_ref': 'cc', 'het': 'Cc', 'hom_alt': 'CC'}[primary['consensus']]
        else:  # E/e
            genotype = {'hom_ref': 'ee', 'het': 'Ee', 'hom_alt': 'EE'}[primary['consensus']]

        supporting = [primary_name]
        discordant = []

        for conf_name in confirmatory_snps:
            conf = snp_consensus.get(conf_name, {})
            if conf.get('consensus') in (None, 'no_call'):
                continue
            conf_info = RHCE_DIAGNOSTIC_SNPS[conf_name]
            # Parity-inverted markers (NM_020485.7 exception) have their
            # ref/alt swapped vs the global HGVS reference, so a patient
            # who matches the primary's zygosity is biologically equivalent
            # to the FLIPPED zygosity at the inverted confirmatory marker.
            effective_conf_zyg = _flip_zygosity(conf['consensus']) \
                if conf_info.get('parity_inverted') else conf['consensus']
            if effective_conf_zyg == primary['consensus']:
                supporting.append(conf_name)
            else:
                discordant.append({
                    'snp': conf_name,
                    'expected': primary['consensus'],
                    'observed': conf['consensus'],
                    'parity_inverted': conf_info.get('parity_inverted', False),
                })

        partial_markers = []
        for ap_name in asian_partial_snps:
            ap = snp_consensus.get(ap_name, {})
            ap_info = RHCE_DIAGNOSTIC_SNPS[ap_name]
            # The "variant allele present" zygosities depend on parity:
            # - Normal: variant present when hom_alt or het.
            # - Inverted: variant present when hom_ref or het (because the
            #   file's ref_base IS the variant after the parity swap).
            if ap_info.get('parity_inverted'):
                trigger = {'het', 'hom_ref'}
            else:
                trigger = {'het', 'hom_alt'}
            if ap.get('consensus') in trigger:
                partial_markers.append({
                    'snp': ap_name,
                    'zygosity': ap['consensus'],
                    'significance': ap_info['significance'],
                })

        # Confidence: start from primary's, downgrade on confirmatory disagreement.
        confidence = primary['confidence']
        if discordant:
            confidence = _min_confidence(confidence, 'LOW')

        reason_parts = [f'Primary {primary_name} = {primary["consensus"]} ({primary["confidence"]})']
        if supporting[1:]:
            reason_parts.append(f'confirmed by {len(supporting) - 1} marker(s)')
        if discordant:
            reason_parts.append(f'{len(discordant)} confirmatory marker(s) disagree')
        if partial_markers:
            reason_parts.append(f'{len(partial_markers)} Asian partial-E marker(s) detected')

        return {
            'axis': axis,
            'genotype': genotype,
            'confidence': confidence,
            'primary_snp': primary_name,
            'supporting_snps': supporting,
            'discordant_snps': discordant,
            'partial_markers': partial_markers,
            'reason': '; '.join(reason_parts),
        }

    @staticmethod
    def _phenotype_and_alleles(c_axis: dict, e_axis: dict, snp_consensus: dict) -> dict:
        """
        Combine C/c and E/e axis calls into a phenotype label and possible
        ISBT haplotype pairs.

        Sanger genotyping cannot phase haplotypes, so for compound
        heterozygotes (CcEe) we report both possible haplotype pairs.
        """
        c_geno = c_axis['genotype']
        e_geno = e_axis['genotype']

        if c_geno is None or e_geno is None:
            return {
                'phenotype': 'Indeterminate',
                'allele_options': [],
                'reason': (
                    f"C/c: {c_axis['reason']} | E/e: {e_axis['reason']}"
                ),
            }

        phenotype = f'{c_geno}{e_geno}'  # e.g. 'CcEe', 'ccee', 'CCEE'

        # ISBT haplotype shortcodes for the 4 standard alleles.
        # RHCE*01 = ce, *02 = Ce, *03 = cE, *04 = CE
        c_alleles_per_strand = {
            'CC': ['C', 'C'],
            'Cc': ['C', 'c'],
            'cc': ['c', 'c'],
        }[c_geno]
        e_alleles_per_strand = {
            'EE': ['E', 'E'],
            'Ee': ['E', 'e'],
            'ee': ['e', 'e'],
        }[e_geno]

        # Enumerate haplotype pairs. For double-het only, both phase options
        # exist; for everything else, exactly one combination is possible.
        haplotype_to_isbt = {
            'ce': 'RHCE*01',
            'Ce': 'RHCE*02',
            'cE': 'RHCE*03',
            'CE': 'RHCE*04',
        }
        # Build the two strands; if double-het, also build the swap.
        strand_a_options = [
            c_alleles_per_strand[0] + e_alleles_per_strand[0],
            c_alleles_per_strand[0] + e_alleles_per_strand[1] if c_geno == 'Cc' and e_geno == 'Ee' else None,
        ]
        strand_b_options = [
            c_alleles_per_strand[1] + e_alleles_per_strand[1],
            c_alleles_per_strand[1] + e_alleles_per_strand[0] if c_geno == 'Cc' and e_geno == 'Ee' else None,
        ]

        allele_options = []
        seen = set()
        for sa, sb in zip(strand_a_options, strand_b_options):
            if sa is None or sb is None:
                continue
            pair_key = tuple(sorted([sa, sb]))
            if pair_key in seen:
                continue
            seen.add(pair_key)
            allele_options.append({
                'haplotypes': f'{sa}/{sb}',
                'isbt': f'{haplotype_to_isbt[sa]}/{haplotype_to_isbt[sb]}',
            })

        # Asian partial-E annotation.
        partial_note = ''
        if e_axis.get('partial_markers'):
            sigs = '; '.join(p['significance'] for p in e_axis['partial_markers'])
            partial_note = f' (partial-E variant flagged: {sigs})'

        reason = (
            f"C/c genotype {c_geno} ({c_axis['confidence']}), "
            f"E/e genotype {e_geno} ({e_axis['confidence']}){partial_note}"
        )

        return {
            'phenotype': phenotype + partial_note,
            'allele_options': allele_options,
            'reason': reason,
        }


# ─── Helpers ─────────────────────────────────────────────────────────────────
_CONFIDENCE_ORDER = {'NONE': 0, 'LOW': 1, 'MEDIUM': 2, 'HIGH': 3}


def _min_confidence(*levels: str) -> str:
    """Return the lowest confidence label among the inputs."""
    return min(levels, key=lambda lvl: _CONFIDENCE_ORDER.get(lvl, 0))


_ZYGOSITY_FLIP = {'hom_ref': 'hom_alt', 'hom_alt': 'hom_ref', 'het': 'het'}


def _flip_zygosity(zyg: str) -> str:
    """Flip hom_ref <-> hom_alt for parity-inverted SNPs; het stays het."""
    return _ZYGOSITY_FLIP.get(zyg, zyg)
