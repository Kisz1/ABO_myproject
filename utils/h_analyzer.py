"""
H Analyzer - ISBT-Informed H Blood Group (FUT1 / Bombay) Analysis

References:
- ISBT 018 (H system) - FUT1 gene (alpha-1,2-fucosyltransferase)
- NM_000148.4 (FUT1 mRNA, cDNA reference)
- NG_007510.2 (FUT1 genomic, RefSeqGene; bundles FGF21 + FUT1 + IZUMO1 —
  the analyzer filters by gene name when loading)

Clinical context: FUT1 builds the H antigen, which is the precursor for
ABO antigens. Loss-of-function variants on both haplotypes produce the
**Bombay phenotype (Oh)** — RBCs lack H, A, and B antigens. Bombay
recipients REJECT all conventionally ABO-typed donor blood and require
Bombay-compatible units. Detecting Bombay carriers (Para-Bombay) is
clinically critical pre-transfusion.

Diagnostic SNP scope (Phase 1):
    c.725T>G  (p.Leu242Arg) - Indian Bombay; most common FUT1 null worldwide
    c.586C>T  (p.Gln196*)   - nonsense mutation; classic Bombay null
    c.460T>C  (p.Tyr154His) - FUT1*01M.01 / h2 weak allele (Para-Bombay)

All three are point mutations in the single FUT1 coding exon. Each was
verified against NG_007510.2 codons before being added to this list —
see tmp_h_verify.py in commit history for the verification grid.

----------------------------------------------------------------------------
Reference-file audit note
----------------------------------------------------------------------------
NG_007510.2 carries the wild-type FUT1 reference bases at all three SNP
sites (T, C, T respectively). No parity inversion needed. If a future
NCBI revision flips any of these, regression tests in tests/test_h_analyzer.py
will fail with a clear message pointing to this file.

----------------------------------------------------------------------------
Architecture
----------------------------------------------------------------------------
Probe-based SNP calling (same pattern as RHCE / KEL / FY / JK). All three
SNPs are coding-region positions, so works against either cDNA or genomic
reference; the analyzer prefers genomic for consistency with other systems.

No primary antigen-discriminator SNP — H is a binary system (H+ vs Bombay).
The phenotype-resolution function `_call_h_phenotype` aggregates the three
null SNPs into one of three states: H+, Para-Bombay potential, or Bombay (Oh).
"""

from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq


# ─── Reference file locations ────────────────────────────────────────────────
DATA_DIR = Path(__file__).resolve().parent / "data"
H_GENBANK_PATH = DATA_DIR / "h_referance.gb"
H_FASTA_PATH = DATA_DIR / "h_cdna_reference.fasta"


# ─── ISBT-defined diagnostic SNPs ────────────────────────────────────────────
# All three are 'cdna' position type (single coding exon in FUT1). All are
# 'modifier' role — H has no primary antigen discriminator, only nulls.
# `severity` distinguishes full nulls (Bombay) from weak alleles (Para-Bombay).
H_DIAGNOSTIC_SNPS = {
    'c.725T>G': {
        'position_type':    'cdna',
        'cDNA_position':    725,
        'exon':             4,
        'ref_base':         'T',
        'alt_base':         'G',
        'antigen_axis':     'expression',
        'ref_call':         'functional',
        'alt_call':         'silenced',
        'role':             'modifier',
        'severity':         'strong',   # full loss-of-function -> Bombay
        'significance':     ('p.Leu242Arg - Indian Bombay; most common FUT1 '
                             'loss-of-function variant worldwide'),
        'reference':        'ISBT 018 (FUT1*01N.04)',
    },
    'c.586C>T': {
        'position_type':    'cdna',
        'cDNA_position':    586,
        'exon':             4,
        'ref_base':         'C',
        'alt_base':         'T',
        'antigen_axis':     'expression',
        'ref_call':         'functional',
        'alt_call':         'silenced',
        'role':             'modifier',
        'severity':         'strong',   # nonsense -> truncation -> Bombay
        'significance':     ('p.Gln196* - nonsense mutation; truncates FUT1 '
                             'and produces classic Bombay null'),
        'reference':        'ISBT 018 (FUT1 nonsense; verified vs NG_007510.2)',
    },
    'c.460T>C': {
        'position_type':    'cdna',
        'cDNA_position':    460,
        'exon':             4,
        'ref_base':         'T',
        'alt_base':         'C',
        'antigen_axis':     'expression',
        'ref_call':         'functional',
        'alt_call':         'weak',
        'role':             'modifier',
        'severity':         'weak',     # partial activity -> Para-Bombay
        'significance':     ('p.Tyr154His - h2 weak allele; reduced FUT1 '
                             'activity, can produce Para-Bombay when '
                             'homozygous or compound heterozygous'),
        'reference':        'ISBT 018 (FUT1*01M.01)',
    },
}

# Names by which the FUT1 gene may be tagged in a GenBank record.
H_GENE_NAMES = ('FUT1', 'H', 'FT1', 'BCA')


# ─── IUPAC decode (same convention as other analyzers) ───────────────────────
IUPAC_DECODE = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC',
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG',
    'N': 'ACGT',
}


# ─── Quality gates (mirror jk/fy/kel for behavioural parity) ─────────────────
MIN_IDENTITY_FOR_CALLING = 85.0
SNP_LOCAL_WINDOW_HALF = 50
SNP_LOCAL_WINDOW_MIN_BASES = 30
MIN_LOCAL_IDENTITY_FOR_CALLING = 85.0
SNP_PROBE_HALF = 50
SNP_PROBE_MIN_BASES = 60
MIN_PHRED_AT_SNP = 30                    # per-base Phred Q-score required AT the
                                         # SNP column when AB1 quality is available.
                                         # Q30 = 1-in-1000 base error rate (lab std).


class HReferenceMissingError(RuntimeError):
    """Raised when no H (FUT1) reference file is available."""


class HAnalyzer:
    """H blood group (FUT1) analyzer covering three diagnostic null SNPs.

    Phase 1 scope is point-mutation Bombay variants only; deletions
    (c.547delC, c.880-882delTTC, etc.) are not yet handled.

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

        gb = Path(gb_path) if gb_path else H_GENBANK_PATH
        fa = Path(fasta_path) if fasta_path else H_FASTA_PATH

        if gb.exists():
            self._load_genbank(gb)
        elif fa.exists():
            self._load_fasta(fa)
        # else: stay unloaded — analyze_* raises with download guidance.

    # ─── Reference loading ───────────────────────────────────────────────
    # NG_007510.2 is a multi-gene RefSeqGene (FGF21 + FUT1 + IZUMO1). The
    # gene-name filter is critical — without it, the first CDS encountered
    # (FGF21 on the minus strand) becomes `cds_start_genomic` and every
    # cDNA SNP position is mis-mapped. Same fix pattern as FY's CADM3 case.
    @classmethod
    def _feature_is_h_gene(cls, feature) -> bool:
        gene = feature.qualifiers.get('gene', [''])[0] or ''
        return gene in H_GENE_NAMES

    def _load_genbank(self, path: Path) -> None:
        record = SeqIO.read(str(path), "genbank")
        self.reference_seq = str(record.seq).upper()
        self.reference_kind = 'genomic'

        exons = []
        cds_start = None
        any_h_tagged = any(
            self._feature_is_h_gene(f)
            for f in record.features
            if f.type in ('gene', 'CDS', 'exon', 'mRNA')
        )
        for feature in record.features:
            if any_h_tagged and not self._feature_is_h_gene(feature):
                continue
            if feature.type == "exon":
                start = int(feature.location.start)
                end = int(feature.location.end)
                number = feature.qualifiers.get("number", [None])[0]
                # FUT1 has alt-spliced exon-4 variants annotated as 4a/4b/4c/4d;
                # strip any non-numeric suffix when parsing the exon number.
                try:
                    import re
                    m = re.match(r'^\d+', str(number) if number is not None else '')
                    number = int(m.group()) if m else None
                except (TypeError, ValueError):
                    number = None
                exons.append((start, end, number))
            elif feature.type == "CDS" and cds_start is None:
                cds_start = int(feature.location.start)

        # Dedup by (start, end) so alt-spliced exon variants don't bloat the
        # walk-the-exons cDNA mapping. Not strictly needed for H since the
        # entire CDS sits inside one downstream exon, but kept for parity
        # with the JK analyzer.
        seen = set()
        deduped = []
        for s, e, n in sorted(exons, key=lambda x: x[0]):
            if (s, e) in seen:
                continue
            seen.add((s, e))
            deduped.append((s, e, n))
        for i, (s, e, n) in enumerate(deduped):
            if n is None:
                deduped[i] = (s, e, i + 1)
        self.exon_map = deduped
        self.cds_start_genomic = cds_start

    def _load_fasta(self, path: Path) -> None:
        record = next(SeqIO.parse(str(path), "fasta"))
        self.reference_seq = str(record.seq).upper()
        self.reference_kind = 'cdna'
        self.cdna_offset = self._detect_cds_start_in_mrna(self.reference_seq)

    @staticmethod
    def _detect_cds_start_in_mrna(seq: str) -> int:
        """Return the 0-based mRNA index of the ATG start codon for FUT1.

        FUT1 protein begins M-W-L-R-S (UniProt P19526), so the start-codon
        context in the mRNA is `ATGTGGCTCCGGAGC`. Falls back to first ATG.
        """
        fut1_context = seq.find('ATGTGGCTCCGGAGC')
        if fut1_context >= 0:
            return fut1_context
        first_atg = seq.find('ATG')
        return first_atg if first_atg >= 0 else 0

    def is_loaded(self) -> bool:
        return self.reference_seq is not None

    def _require_reference(self) -> None:
        if not self.is_loaded():
            raise HReferenceMissingError(
                "H (FUT1) reference not found. Supply one of:\n"
                f"  - GenBank: {H_GENBANK_PATH}\n"
                f"  - FASTA:   {H_FASTA_PATH}\n"
                "See utils/data/H_REFERENCE_README.md for download instructions."
            )

    # ─── Position translation ────────────────────────────────────────────
    def _position_to_ref_index(self, snp_info: dict) -> Optional[int]:
        """Map a SNP's cDNA position to a 0-based reference index.

        Phase 1 only handles 'cdna' positions — no intronic or promoter SNPs.
        """
        cdna_pos = int(snp_info['cDNA_position'])

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
    def _probe_alignment_identity(alignment) -> tuple:
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
        ref_aligned = str(alignment[0])
        query_aligned = str(alignment[1])
        ref_b_at_col = ref_aligned[col].upper()
        query_b = query_aligned[col].upper()

        meta = {
            'local_window_identity': round(window_pct, 1),
            'local_window_bases':    window_n,
        }

        # Phred Q-score gate at the SNP column (only when AB1 quality is plumbed).
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
                'covered':  True,
                'call':     'no_call',
                'reason':   (f"ref_base mismatch: reference has {ref_b_at_col}, "
                             f"expected {snp_info['ref_base']}"),
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
                'covered':  True,
                'call':     'no_call',
                'reason':   f'unexpected base {query_b} (decoded={decoded})',
                **meta,
            }

        return {
            'covered':       True,
            'call':          call,
            'zygosity':      zygosity,
            'query_base':    query_b,
            'ref_base':      ref_b,
            'alt_base':      alt_b,
            'antigen_axis':  snp_info['antigen_axis'],
            'role':          snp_info['role'],
            'severity':      snp_info.get('severity'),
            'exon':          snp_info.get('exon'),
            'cDNA_position': snp_info.get('cDNA_position'),
            'significance':  snp_info['significance'],
            **meta,
        }

    def _genotype_snp_via_probe(
        self,
        query_seq: str,
        snp_info: dict,
        query_quality: Optional[list] = None,
    ) -> dict:
        ref_index = self._position_to_ref_index(snp_info)
        if ref_index is None:
            return {'covered': False, 'reason': 'position outside reference'}

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
                'covered':                False,
                'reason':                 f'probe alignment too short ({best["n"]} aligned bases)',
                'local_window_identity':  round(best['identity'], 1),
                'local_window_bases':     best['n'],
            }
        if best['identity'] < MIN_LOCAL_IDENTITY_FOR_CALLING:
            return {
                'covered':                True,
                'call':                   'no_call',
                'reason':                 (f'probe identity {best["identity"]:.1f}% '
                                           f'< {MIN_LOCAL_IDENTITY_FOR_CALLING}%'),
                'local_window_identity':  round(best['identity'], 1),
                'local_window_bases':     best['n'],
            }

        col = self._ref_index_to_alignment_column(best['alignment'], snp_offset_in_probe)
        if col is None:
            return {
                'covered':                False,
                'reason':                 'SNP position outside probe alignment or in query gap',
                'local_window_identity':  round(best['identity'], 1),
                'local_window_bases':     best['n'],
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
        """Probe-only single-read analysis.

        Skips the full-read alignment that earlier versions used as a coarse
        identity gate (slow on 14 kb genomic reference for synthetic-patch
        tests). Per-SNP probe gate is the actual calling criterion; identity
        and strand are derived from the strongest probe call.

        ``query_quality`` enables the Phred-Q30-at-SNP-column gate.
        """
        self._require_reference()

        query_str = str(query_seq).upper()
        query_length = len(query_str)

        snp_calls: dict = {
            snp_name: self._genotype_snp_via_probe(
                query_str, snp_info, query_quality=query_quality,
            )
            for snp_name, snp_info in H_DIAGNOSTIC_SNPS.items()
        }

        # Derive identity / strand from the strongest probe call.
        best = None
        for call in snp_calls.values():
            ident = call.get('local_window_identity')
            if ident is None:
                continue
            if best is None or ident > best['identity']:
                best = {
                    'identity': ident,
                    'strand':   call.get('probe_strand'),
                }

        trusted_snps = sum(
            1 for c in snp_calls.values()
            if c.get('covered') and c.get('call') and c.get('call') != 'no_call'
        )
        low_phred_at_snp = sum(
            1 for c in snp_calls.values()
            if (c.get('phred_at_snp') is not None
                and c.get('phred_at_snp') < MIN_PHRED_AT_SNP)
        )

        if best is None:
            return {
                'read_id':         read_id,
                'query_length':    query_length,
                'strand':          None,
                'identity':        0.0,
                'score':           0.0,
                'reference_kind':  self.reference_kind,
                'snp_calls':       snp_calls,
                'callable':        False,
                'trusted_snps':    0,
                'low_phred_at_snp': low_phred_at_snp,
                'reason':          'No probe alignment produced for any SNP',
            }

        callable_read = (trusted_snps > 0
                         or best['identity'] >= MIN_IDENTITY_FOR_CALLING)

        reason = (
            f"Best probe identity {best['identity']:.1f}% ({best['strand']}); "
            f"{trusted_snps} diagnostic SNP(s) cleared the probe-identity gate"
            + (f"; {low_phred_at_snp} blocked by Phred Q{MIN_PHRED_AT_SNP} gate"
               if low_phred_at_snp else "")
            + "."
        )

        return {
            'read_id':         read_id,
            'query_length':    query_length,
            'strand':          best['strand'],
            'identity':        round(best['identity'], 2),
            'score':           0.0,  # full-read score no longer computed
            'reference_kind':  self.reference_kind,
            'snp_calls':       snp_calls,
            'callable':        callable_read,
            'trusted_snps':    trusted_snps,
            'low_phred_at_snp': low_phred_at_snp,
            'reason':          reason,
        }

    # ─── Multi-read consensus ────────────────────────────────────────────
    def analyze(self, reads) -> dict:
        self._require_reference()

        normalized = self._normalize_reads(reads)
        per_read = [
            self.analyze_single_read(seq, read_id=rid, query_quality=qual)
            for rid, seq, qual in normalized
        ]

        snp_consensus = self._consensus_snp_calls(per_read)
        h_state = self._call_h_phenotype(snp_consensus)
        phenotype_block = self._phenotype_and_alleles(h_state)

        callable_reads = [r for r in per_read if r['callable']]
        reads_with_trusted_snps = [r for r in per_read if r.get('trusted_snps', 0) > 0]

        return {
            'phenotype':               phenotype_block['phenotype'],
            'allele_options':          phenotype_block['allele_options'],
            'h_state_call':            h_state,
            'snp_consensus':           snp_consensus,
            'no_calls_primary':        [],   # H has no primary SNP
            'overall_confidence':      h_state['confidence'],
            'reads_total':             len(per_read),
            'reads_callable':          len(callable_reads),
            'reads_with_trusted_snps': len(reads_with_trusted_snps),
            'per_read_details':        per_read,
            'reason':                  phenotype_block['reason'],
        }

    @staticmethod
    def _normalize_reads(reads) -> list:
        """Normalize reads input to ``[(read_id, seq, quality_or_None), ...]``."""
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
        """Vote on each diagnostic SNP across reads. Same LOW/MEDIUM/HIGH
        gates as FY/JK/KEL — see those analyzers for the rule rationale."""
        consensus: dict = {}

        for snp_name in H_DIAGNOSTIC_SNPS:
            covering = []
            for read in per_read_results:
                call = read['snp_calls'].get(snp_name)
                if not call or not call.get('covered'):
                    continue
                if call.get('call') == 'no_call':
                    continue
                covering.append({
                    'read_id':    read['read_id'],
                    'zygosity':   call['zygosity'],
                    'call':       call['call'],
                    'query_base': call.get('query_base'),
                })

            if not covering:
                consensus[snp_name] = {
                    'consensus':        'no_call',
                    'confidence':       'NONE',
                    'reads_covering':   0,
                    'votes':            {},
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
                'consensus':        top_zygosity,
                'call':             consensus_call,
                'confidence':       confidence,
                'reads_covering':   n_reads,
                'votes':            votes,
                'supporting_reads': [c['read_id'] for c in supporting],
                'discordant_reads': [c['read_id'] for c in discordant],
            }

        return consensus

    @staticmethod
    def _call_h_phenotype(snp_consensus: dict) -> dict:
        """Aggregate the three null SNPs into one H-system phenotype state.

        Returns a dict with:
          state: 'h_positive' | 'para_bombay' | 'bombay' | 'indeterminate'
          confidence: HIGH / MEDIUM / LOW / NONE
          triggers_hom_alt: list of SNP entries that homozygously alt
          triggers_het: list of SNP entries that are heterozygous
          callable_ref: list of SNP entries that are homozygous reference
          any_strong_hom_alt: whether a 'strong' (full-null) variant is homozygous
          any_strong_het: whether a 'strong' variant is heterozygous

        State rules:
          - any strong null hom_alt          -> 'bombay'      (Oh, both haplotypes null)
          - any weak null hom_alt (only)     -> 'para_bombay' (h2/h2 weak)
          - any null het (strong or weak)    -> 'para_bombay' (one haplotype affected)
          - all callable nulls are hom_ref   -> 'h_positive'
          - none callable                    -> 'indeterminate'
        """
        triggers_hom_alt = []
        triggers_het = []
        callable_ref = []
        any_strong_hom_alt = False
        any_strong_het = False

        for snp_name, info in H_DIAGNOSTIC_SNPS.items():
            cons = snp_consensus.get(snp_name)
            if not cons or cons['consensus'] in (None, 'no_call'):
                continue
            entry = {'snp': snp_name, 'severity': info.get('severity'),
                     'confidence': cons['confidence'],
                     'significance': info['significance']}
            if cons['consensus'] == 'hom_alt':
                triggers_hom_alt.append(entry)
                if info.get('severity') == 'strong':
                    any_strong_hom_alt = True
            elif cons['consensus'] == 'het':
                triggers_het.append(entry)
                if info.get('severity') == 'strong':
                    any_strong_het = True
            elif cons['consensus'] == 'hom_ref':
                callable_ref.append(entry)

        if any_strong_hom_alt:
            state = 'bombay'
            confidence = next(t['confidence'] for t in triggers_hom_alt
                              if t['severity'] == 'strong')
        elif triggers_hom_alt:
            # weak-only hom_alt -> Para-Bombay (h2/h2)
            state = 'para_bombay'
            confidence = triggers_hom_alt[0]['confidence']
        elif triggers_het:
            state = 'para_bombay'
            confidence = triggers_het[0]['confidence']
        elif callable_ref:
            state = 'h_positive'
            confs = [c['confidence'] for c in callable_ref]
            confidence = ('HIGH' if 'HIGH' in confs
                          else 'MEDIUM' if 'MEDIUM' in confs
                          else 'LOW')
        else:
            state = 'indeterminate'
            confidence = 'NONE'

        return {
            'state':              state,
            'confidence':         confidence,
            'triggers_hom_alt':   triggers_hom_alt,
            'triggers_het':       triggers_het,
            'callable_ref':       callable_ref,
            'any_strong_hom_alt': any_strong_hom_alt,
            'any_strong_het':     any_strong_het,
        }

    @staticmethod
    def _phenotype_and_alleles(h_state: dict) -> dict:
        """Map the combined H state to a phenotype label and ISBT haplotype(s)."""
        state = h_state['state']

        if state == 'indeterminate':
            return {
                'phenotype':      'Indeterminate',
                'allele_options': [],
                'reason':         'No diagnostic SNPs were callable on this run.',
            }

        if state == 'h_positive':
            return {
                'phenotype':      'H+',
                'allele_options': [{
                    'haplotypes': 'FUT1*A/FUT1*A',
                    'isbt':       'FUT1*01/FUT1*01',
                    'serology':   'H+ (normal H antigen expression on RBCs)',
                }],
                'reason': ('No silencing variants detected; FUT1 reference '
                           'on both haplotypes.'),
            }

        if state == 'bombay':
            triggers = ', '.join(t['snp'] for t in h_state['triggers_hom_alt'])
            return {
                'phenotype':      'Bombay (Oh) — H-negative',
                'allele_options': [{
                    'haplotypes': 'FUT1*0/FUT1*0',
                    'isbt':       'FUT1*null / FUT1*null',
                    'serology':   ('Bombay (Oh) — RBCs lack H, A, and B antigens. '
                                   f'Trigger(s): {triggers}. '
                                   '**Rejects all standard ABO-typed donor blood; '
                                   'requires Bombay-compatible units.**'),
                }],
                'reason': (f'Loss-of-function variant(s) homozygous: {triggers}; '
                           'FUT1 protein non-functional on both haplotypes.'),
            }

        # para_bombay
        het_triggers = ', '.join(t['snp'] for t in h_state['triggers_het'])
        weak_hom = [t for t in h_state['triggers_hom_alt']
                    if t['severity'] == 'weak']
        weak_str = ', '.join(t['snp'] for t in weak_hom)

        if weak_hom and not h_state['triggers_het']:
            # weak hom_alt (e.g. c.460 C/C) — h2/h2 Para-Bombay
            return {
                'phenotype':      'Para-Bombay (h2 homozygous)',
                'allele_options': [{
                    'haplotypes': 'FUT1*M/FUT1*M',
                    'isbt':       'FUT1*01M.01 / FUT1*01M.01',
                    'serology':   ('Para-Bombay — weak H expression. '
                                   f'Both haplotypes carry the weak allele ({weak_str}). '
                                   'May serotype as H-weak/negative; cross-match cautiously.'),
                }],
                'reason': (f'Weak allele(s) homozygous: {weak_str}; '
                           'reduced but residual FUT1 activity expected.'),
            }

        # null het(s) -> Para-Bombay potential (one haplotype affected)
        return {
            'phenotype':      'Para-Bombay potential',
            'allele_options': [{
                'haplotypes': f'FUT1*A / FUT1*affected({het_triggers})',
                'isbt':       'FUT1*01 / FUT1*null-or-weak',
                'serology':   ('Heterozygous carrier — one functional FUT1 + one '
                               f'affected haplotype ({het_triggers}). Typically H+ '
                               'serologically; weak/variable expression possible.'),
            }],
            'reason': (f'Null/weak variant(s) heterozygous: {het_triggers}. '
                       'One haplotype affected; other expresses normal FUT1.'),
        }
