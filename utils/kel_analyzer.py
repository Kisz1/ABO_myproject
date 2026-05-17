"""
KEL Analyzer - ISBT-Informed Kell (K/k) Analysis

References:
- ISBT 006 (KEL system) v6.4
- NM_000420.4 (KEL mRNA, cDNA reference)
- NG_007492.1 (KEL genomic, optional)

Diagnostic SNP scope (Phase 1):
    K/k: c.578C>T (p.Thr193Met) - primary K1/K2 discriminator
         (KEL*01 = K = alt allele, KEL*02 = k = reference allele)

Phase 2 (drop additional entries into KEL_DIAGNOSTIC_SNPS when those
amplicons exist):
    Kp(a)/Kp(b): c.841C>T  (exon 8)
    Js(a)/Js(b): c.1790T>C (exon 17)

----------------------------------------------------------------------------
Reference-file audit note
----------------------------------------------------------------------------
NM_000420.4 SHOULD carry the k-allele (KEL*02) base at c.578, i.e. file_base
'C'. This is true for the publicly deposited record at the time of writing,
but K1 carriers are rare so it would also be possible for a future revision
to inherit a K-allele clone — in which case the parity-inversion treatment
from RHCEAnalyzer's c.178 entry should be applied here.

The `test_file_base_at_c578_matches_audit_note` regression test enforces this.
If it fails after a reference refresh, swap ref_base/alt_base and ref_call/
alt_call in the SNP table and add `parity_inverted: True`, mirroring RHCE.

Architecture: probe-based SNP calling (identical pattern to RHCEAnalyzer's
fix for genomic Sanger reads against a cDNA reference — the small reference
probe for each SNP is exon-local, so it lands on the right exonic region in
the read regardless of surrounding introns).
"""

from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq


# ─── Reference file locations ────────────────────────────────────────────────
DATA_DIR = Path(__file__).resolve().parent / "data"
KEL_GENBANK_PATH = DATA_DIR / "kel_referance.gb"
KEL_FASTA_PATH = DATA_DIR / "kel_cdna_reference.fasta"


# ─── ISBT-defined diagnostic SNPs for K/k calling ────────────────────────────
# cDNA_position is 1-based per NM_000420.4 mRNA convention.
KEL_DIAGNOSTIC_SNPS = {
    'c.578C>T': {
        'exon': 6,
        'cDNA_position': 578,
        'ref_base': 'C',
        'alt_base': 'T',
        'antigen_axis': 'K/k',
        'ref_call': 'k',           # reference allele expresses k (KEL*02)
        'alt_call': 'K',           # variant allele expresses K (KEL*01)
        'role': 'primary',
        'significance': 'p.Thr193Met - primary K/k discriminator',
        'reference': 'ISBT 006 v6.4',
    },
}


# ─── IUPAC decode (same convention as rhce_analyzer / rhd_analyzer) ──────────
IUPAC_DECODE = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC',
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG',
    'N': 'ACGT',
}


# ─── Quality gates (mirror rhce_analyzer for behavioural parity) ─────────────
MIN_IDENTITY_FOR_CALLING = 85.0          # whole-read identity (reporting only)
SNP_LOCAL_WINDOW_HALF = 50               # whole-read window check (legacy path)
SNP_LOCAL_WINDOW_MIN_BASES = 30
MIN_LOCAL_IDENTITY_FOR_CALLING = 85.0    # probe-identity threshold
SNP_PROBE_HALF = 50                      # probe = ±50 bp around each SNP
SNP_PROBE_MIN_BASES = 60                 # minimum aligned bases for a usable hit
MIN_PHRED_AT_SNP = 30                    # per-base Phred Q-score required AT the
                                         # SNP column when AB1 quality is available.
                                         # Q30 = 1-in-1000 base error rate (lab std).


class KELReferenceMissingError(RuntimeError):
    """Raised when no KEL reference file is available."""


class KELAnalyzer:
    """Single-axis (K/k) Kell blood group analyzer.

    Uses probe-based SNP calling (small reference window per SNP locally
    aligned against the read). This is robust to genomic Sanger amplicons
    containing introns that the cDNA reference doesn't have — exactly the
    issue that broke whole-read alignment on the RHCE45 amplicon.

    Use ``analyze_single_read(query_seq)`` for one read or ``analyze(reads)``
    for multi-read consensus voting.
    """

    def __init__(
        self,
        gb_path: Optional[str] = None,
        fasta_path: Optional[str] = None,
        reference_seq: Optional[str] = None,
    ):
        self.reference_seq: Optional[str] = None
        self.reference_kind: Optional[str] = None  # 'cdna' or 'genomic'
        self.exon_map: list = []
        self.cds_start_genomic: Optional[int] = None
        self.cdna_offset: int = 0
        self._aligner = self._build_aligner()

        if reference_seq is not None:
            self.reference_seq = str(reference_seq).upper()
            self.reference_kind = 'cdna'
            self.cdna_offset = self._detect_cds_start_in_mrna(self.reference_seq)
            return

        gb = Path(gb_path) if gb_path else KEL_GENBANK_PATH
        fa = Path(fasta_path) if fasta_path else KEL_FASTA_PATH

        if gb.exists():
            self._load_genbank(gb)
        elif fa.exists():
            self._load_fasta(fa)
        # else: stay unloaded — analyze_* raises with download guidance.

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
        for i, (s, e, n) in enumerate(exons):
            if n is None:
                exons[i] = (s, e, i + 1)
        self.exon_map = exons
        self.cds_start_genomic = cds_start

    def _load_fasta(self, path: Path) -> None:
        record = next(SeqIO.parse(str(path), "fasta"))
        self.reference_seq = str(record.seq).upper()
        self.reference_kind = 'cdna'
        self.cdna_offset = self._detect_cds_start_in_mrna(self.reference_seq)

    @staticmethod
    def _detect_cds_start_in_mrna(seq: str) -> int:
        """Return the 0-based mRNA index of the ATG start codon for KEL.

        KEL protein begins M-E-G-G-D-Q-T-E-Y-K (UniProt P23276), so the start
        codon context in the mRNA is `ATGGAAGGTGGGGAC`. Falls back to first
        ATG if the specific context is not found.
        """
        kel_context = seq.find('ATGGAAGGTGGGGAC')
        if kel_context >= 0:
            return kel_context
        first_atg = seq.find('ATG')
        return first_atg if first_atg >= 0 else 0

    def is_loaded(self) -> bool:
        return self.reference_seq is not None

    def _require_reference(self) -> None:
        if not self.is_loaded():
            raise KELReferenceMissingError(
                "KEL reference not found. Supply one of:\n"
                f"  - GenBank: {KEL_GENBANK_PATH}\n"
                f"  - FASTA:   {KEL_FASTA_PATH}\n"
                "See utils/data/KEL_REFERENCE_README.md for download instructions."
            )

    # ─── cDNA-position translation ────────────────────────────────────────
    def _cdna_to_ref_index(self, cdna_pos: int) -> Optional[int]:
        if self.reference_kind == 'cdna':
            idx = (cdna_pos - 1) + self.cdna_offset
            if 0 <= idx < len(self.reference_seq):
                return idx
            return None

        remaining = cdna_pos
        for start, end, _ in self.exon_map:
            exon_len = end - start
            exon_cds_start = start
            exon_cds_len = exon_len

            if self.cds_start_genomic is not None and start < self.cds_start_genomic < end:
                exon_cds_start = self.cds_start_genomic
                exon_cds_len = end - self.cds_start_genomic
            elif self.cds_start_genomic is not None and end <= self.cds_start_genomic:
                continue

            if remaining <= exon_cds_len:
                return exon_cds_start + (remaining - 1)
            remaining -= exon_cds_len

        return None

    # ─── Alignment ────────────────────────────────────────────────────────
    @staticmethod
    def _build_aligner() -> PairwiseAligner:
        """Local aligner with the same scoring as RHCE/RHD paths."""
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
        """Map a reference index to its column in a local PairwiseAligner alignment.

        Cursor starts at coordinates[0][0] (the local-alignment start on the
        reference) and increments for each non-gap reference character. Returns
        None if the SNP position is outside the aligned window or in a query gap.
        """
        ref_aligned = str(alignment[0])
        query_aligned = str(alignment[1])
        try:
            cursor = int(alignment.coordinates[0][0])
        except (AttributeError, IndexError, TypeError):
            cursor = 0
        for col, ref_char in enumerate(ref_aligned):
            if ref_char == '-':
                continue
            if cursor == ref_index_0based:
                if query_aligned[col] == '-':
                    return None
                return col
            cursor += 1
        return None

    @staticmethod
    def _query_pos_at_alignment_column(
        alignment, col: int, query_len: int, strand: str
    ) -> Optional[int]:
        """Map an alignment column to a 0-based index in the *original*
        query sequence (before reverse-complement, if applicable).

        Used to look up the Phred Q-score at the called SNP position when
        AB1 quality scores have been threaded through. Returns None if the
        column lands on a query gap or the column is outside the alignment.
        """
        query_aligned = str(alignment[1])
        try:
            q_cursor = int(alignment.coordinates[1][0])
        except (AttributeError, IndexError, TypeError):
            q_cursor = 0
        for c, q_char in enumerate(query_aligned):
            if c == col:
                if q_char == '-':
                    return None
                # q_cursor here is the position in the (possibly RC'd) query
                # that was passed to the aligner. For forward strand it
                # equals the position in the original query_seq.
                if strand == 'reverse':
                    return query_len - 1 - q_cursor
                return q_cursor
            if q_char != '-':
                q_cursor += 1
        return None

    # ─── SNP genotyping ──────────────────────────────────────────────────
    @staticmethod
    def _local_window_identity(alignment, center_col: int) -> tuple:
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

    @staticmethod
    def _probe_alignment_identity(alignment) -> tuple:
        """Identity over the FULL probe alignment (probe is ~100 bp by design)."""
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
        """Shared zygosity logic. Assumes the local-window check already passed.

        When ``query_quality`` is provided, additionally gates the call by
        the per-base Phred Q-score at the SNP column — refuses to call if
        the base is below MIN_PHRED_AT_SNP. Falls back to the legacy
        IUPAC-only path when quality is None (FASTA inputs).
        """
        ref_aligned = str(alignment[0])
        query_aligned = str(alignment[1])
        ref_b_at_col = ref_aligned[col].upper()
        query_b = query_aligned[col].upper()

        meta = {
            'local_window_identity': round(window_pct, 1),
            'local_window_bases': window_n,
        }

        # Phred Q-score gate at the SNP column (only applied if quality
        # is plumbed through — AB1 inputs preprocessed by the QC step).
        phred_at_snp: Optional[int] = None
        if query_quality is not None and query_seq is not None and strand is not None:
            q_pos = self._query_pos_at_alignment_column(
                alignment, col, len(query_seq), strand
            )
            if q_pos is not None and 0 <= q_pos < len(query_quality):
                phred_at_snp = int(query_quality[q_pos])
                meta['phred_at_snp'] = phred_at_snp
                if phred_at_snp < MIN_PHRED_AT_SNP:
                    return {
                        'covered': True,
                        'call': 'no_call',
                        'reason': (f'low Phred at SNP position '
                                   f'(Q={phred_at_snp} < Q{MIN_PHRED_AT_SNP})'),
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

    def _genotype_snp_via_probe(
        self,
        query_seq: str,
        snp_info: dict,
        query_quality: Optional[list] = None,
    ) -> dict:
        """Per-SNP genotyping by aligning a small reference probe against the read.

        ``query_quality`` is the per-base Phred Q-score list for ``query_seq``
        (e.g. ``trace['quality_trimmed']`` from process_rhd_ab1_files). When
        provided, the Q-score at the SNP column is gated against
        MIN_PHRED_AT_SNP. None for FASTA-only inputs (no quality data).
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
        """Analyze one read against the KEL reference.

        ``query_quality`` is optional per-base Phred Q-scores aligned to
        ``query_seq``. When provided, the analyzer additionally gates each
        SNP call by Q-score at the SNP column (default Q30, lab standard).
        FASTA inputs pass None and skip the Phred gate (legacy behaviour).

        Returns the same shape as RHCEAnalyzer.analyze_single_read so the UI
        layer can render it consistently.
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
                'trusted_snps': 0,
                'reason': 'No alignment produced',
            }

        best = max(
            (a for a in (forward, reverse) if a is not None),
            key=lambda a: a['identity'],
        )

        callable_read = best['identity'] >= MIN_IDENTITY_FOR_CALLING
        snp_calls: dict = {
            snp_name: self._genotype_snp_via_probe(
                query_str, snp_info, query_quality=query_quality,
            )
            for snp_name, snp_info in KEL_DIAGNOSTIC_SNPS.items()
        }
        trusted_snps = sum(
            1 for c in snp_calls.values()
            if c.get('covered') and c.get('call') and c.get('call') != 'no_call'
        )
        # Count SNPs gated out by low Phred (separate from probe-identity gate)
        low_phred_at_snp = sum(
            1 for c in snp_calls.values()
            if (c.get('phred_at_snp') is not None
                and c.get('phred_at_snp') < MIN_PHRED_AT_SNP)
        )

        reason = (
            f"Read aligned at {best['identity']:.1f}% identity ({best['strand']}); "
            f"{trusted_snps} diagnostic SNP(s) cleared the probe-identity gate"
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

    # ─── Multi-read consensus ────────────────────────────────────────────
    def analyze(self, reads) -> dict:
        """ISBT-style multi-read KEL pipeline.

        Accepts a single sequence string, a list of strings, a list of
        (read_id, sequence) tuples, or a list of (read_id, sequence,
        quality) triples (quality = per-base Phred Q-scores list, or None
        for FASTA inputs). Returns consensus result with antigen axis
        call, phenotype, ISBT haplotype pair, and confidence.
        """
        self._require_reference()

        normalized = self._normalize_reads(reads)
        per_read = [
            self.analyze_single_read(seq, read_id=rid, query_quality=qual)
            for rid, seq, qual in normalized
        ]

        snp_consensus = self._consensus_snp_calls(per_read)
        k_axis = self._call_antigen_axis(snp_consensus, 'K/k')
        phenotype_block = self._phenotype_and_alleles(k_axis)

        callable_reads = [r for r in per_read if r['callable']]
        reads_with_trusted_snps = [r for r in per_read if r.get('trusted_snps', 0) > 0]
        no_calls = [
            snp for snp, cons in snp_consensus.items()
            if cons['consensus'] == 'no_call' and KEL_DIAGNOSTIC_SNPS[snp]['role'] == 'primary'
        ]

        return {
            'phenotype': phenotype_block['phenotype'],
            'allele_options': phenotype_block['allele_options'],
            'k_axis_call': k_axis,
            'snp_consensus': snp_consensus,
            'no_calls_primary': no_calls,
            'overall_confidence': k_axis['confidence'],
            'reads_total': len(per_read),
            'reads_callable': len(callable_reads),
            'reads_with_trusted_snps': len(reads_with_trusted_snps),
            'per_read_details': per_read,
            'reason': phenotype_block['reason'],
        }

    @staticmethod
    def _normalize_reads(reads) -> list:
        """Normalize reads input to ``[(read_id, seq, quality_or_None), ...]``.

        Accepts: a single seq str, a list of seq strs, a list of (rid, seq)
        2-tuples (FASTA-style, no quality), or a list of (rid, seq, quality)
        3-tuples (AB1-style, per-base Phred Q-scores). 2-tuples are widened
        to 3-tuples with quality=None for backward compat.
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
        """Vote on each diagnostic SNP across all reads.

        Identical rules to RHCEAnalyzer:
          - ≥2 reads, all concordant on zygosity  -> HIGH
          - ≥2 reads, majority agrees             -> MEDIUM
          - exactly 1 read                        -> MEDIUM
          - tie or no clear majority              -> LOW
          - 0 reads cover the SNP                 -> no_call
        """
        consensus: dict = {}

        for snp_name in KEL_DIAGNOSTIC_SNPS:
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
                'consensus': top_zygosity,
                'call': consensus_call,
                'confidence': confidence,
                'reads_covering': n_reads,
                'votes': votes,
                'supporting_reads': [c['read_id'] for c in supporting],
                'discordant_reads': [c['read_id'] for c in discordant],
            }

        return consensus

    @staticmethod
    def _call_antigen_axis(snp_consensus: dict, axis: str) -> dict:
        """Combine SNPs on one axis (K/k) into a single genotype call.

        KEL Phase 1 has one primary SNP (c.578C>T) and no confirmatory or
        partial markers — so this is a thin wrapper around the primary SNP
        consensus. Kept structurally similar to RHCE's version so future
        confirmatory SNPs can drop straight in.
        """
        primary_snps = [
            name for name, info in KEL_DIAGNOSTIC_SNPS.items()
            if info['antigen_axis'] == axis and info['role'] == 'primary'
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
                'reason': f'Primary marker {primary_name} not callable',
            }

        # K/k genotype mapping.
        genotype = {'hom_ref': 'kk', 'het': 'Kk', 'hom_alt': 'KK'}[primary['consensus']]

        return {
            'axis': axis,
            'genotype': genotype,
            'confidence': primary['confidence'],
            'primary_snp': primary_name,
            'supporting_snps': [primary_name],
            'discordant_snps': [],
            'reason': f'Primary {primary_name} = {primary["consensus"]} ({primary["confidence"]})',
        }

    @staticmethod
    def _phenotype_and_alleles(k_axis: dict) -> dict:
        """Combine the K/k axis call into a phenotype label and ISBT haplotype pair.

        Single-axis system — exactly one haplotype option, no phase ambiguity.
        """
        k_geno = k_axis['genotype']

        if k_geno is None:
            return {
                'phenotype': 'Indeterminate',
                'allele_options': [],
                'reason': f"K/k: {k_axis['reason']}",
            }

        # Phenotype label uses the genetic-style genotype (matches RHCE).
        # Serological equivalents:  KK -> K+k-,  Kk -> K+k+,  kk -> K-k+
        allele_map = {
            'KK': ('K', 'K', 'KEL*01/KEL*01', 'K+k-'),
            'Kk': ('K', 'k', 'KEL*01/KEL*02', 'K+k+'),
            'kk': ('k', 'k', 'KEL*02/KEL*02', 'K-k+'),
        }
        a1, a2, isbt, serology = allele_map[k_geno]

        return {
            'phenotype': k_geno,
            'allele_options': [{
                'haplotypes': f'{a1}/{a2}',
                'isbt': isbt,
                'serology': serology,
            }],
            'reason': (
                f"K/k genotype {k_geno} ({k_axis['confidence']}); "
                f"serological equivalent {serology}"
            ),
        }
