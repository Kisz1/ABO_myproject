"""
JK Analyzer - ISBT-Informed Kidd (Jka/Jkb/Jk(a-b-)) Analysis

References:
- ISBT 009 (JK system) - SLC14A1 gene (formerly JK / UT-B urea transporter)
- NM_015865.7 (SLC14A1 mRNA, cDNA reference)
- NG_011775.2 (SLC14A1 genomic, includes introns — REQUIRED for c.342-1G>A)

Diagnostic SNP scope (Phase 1):
    JK*A/JK*B:    c.838G>A (p.Asp280Asn) - primary JKA/JKB discriminator
                  JK*01 = JK*A (Jk(a+)), JK*02 = JK*B (Jk(b+))
    JK*Null#1:    c.342-1G>A in the intron-5 splice acceptor
                  Abolishes splicing -> Polynesian / Pacific-Islander Jk(a-b-).
                  *Requires a genomic reference; the intronic position is not
                  in the cDNA.*
    JK*Null#2:    c.871T>C (p.Ser291Pro) - Asian / Finnish Jk(a-b-) point
                  mutation. Encoded protein is unstable / non-functional.

----------------------------------------------------------------------------
Reference-file audit note
----------------------------------------------------------------------------
NG_011775.2 SHOULD carry the JK*A allele base ('G') at c.838 — JK*A is the
more common allele globally and the standard reference deposit. The
regression test
    tests/test_jk_analyzer.py::test_file_base_at_c838_matches_audit_note
enforces this. If a future NCBI revision inherits a JK*B clone (base 'A'),
swap ref_base/alt_base and ref_call/alt_call in the c.838 entry and set
`parity_inverted: True`, mirroring the RHCE c.178 treatment.

----------------------------------------------------------------------------
Architecture
----------------------------------------------------------------------------
Probe-based SNP calling (same pattern as RHCE / KEL / FY): a small reference
probe (±50 bp) around each SNP is locally aligned against the read. This is
intron-tolerant for cDNA-vs-genomic mismatch and lets the splice-site SNP
be addressed via the same code path as the coding SNPs.
"""

from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq


# ─── Reference file locations ────────────────────────────────────────────────
DATA_DIR = Path(__file__).resolve().parent / "data"
JK_GENBANK_PATH = DATA_DIR / "jk_referance.gb"
JK_FASTA_PATH = DATA_DIR / "jk_cdna_reference.fasta"


# ─── ISBT-defined diagnostic SNPs ────────────────────────────────────────────
# `position_type` is one of:
#   'cdna'     : 1-based mRNA position (works against either reference kind)
#   'intronic' : `cDNA_anchor` + `intronic_offset` (HGVS c.<N><offset> form);
#                requires a genomic reference.
JK_DIAGNOSTIC_SNPS = {
    'c.838G>A': {
        'position_type':    'cdna',
        'cDNA_position':    838,
        'exon':             9,
        # Parity inversion: NG_011775.2 deposits the JK*B allele (Asn280, base
        # 'A') at c.838, not the textbook JK*A. So the locally-observed ref
        # base is 'A' and the variant base is 'G', even though the textbook
        # SNP name is "c.838G>A". Same trick as RHCE c.178 — see audit note.
        'ref_base':         'A',
        'alt_base':         'G',
        'antigen_axis':     'A/B',
        'ref_call':         'B',    # the file's reference base encodes JK*B
        'alt_call':         'A',    # the file's variant base encodes JK*A
        'role':             'primary',
        'significance':     'p.Asp280Asn - primary JKA/JKB discriminator',
        'reference':        'ISBT 009 (JK*01 / JK*02)',
        'parity_inverted':  True,
    },
    'c.342-1G>A': {
        'position_type':    'intronic',
        'cDNA_anchor':      342,    # first base of exon 6
        'intronic_offset':  -1,     # 1 bp 5' of c.342 (last base of intron 5)
        'ref_base':         'G',
        'alt_base':         'A',
        'antigen_axis':     'expression',
        'ref_call':         'expressed',
        'alt_call':         'silenced',
        'role':             'modifier',
        'significance':     ('Intron-5 splice acceptor abolished — Polynesian '
                             '/ Pacific-Islander Jk(a-b-) when homozygous'),
        'reference':        'ISBT 009 (JK*01N.05 / JK*02N.04)',
        'requires_reference_kind': 'genomic',
    },
    'c.871T>C': {
        'position_type':    'cdna',
        'cDNA_position':    871,
        'exon':             9,
        'ref_base':         'T',
        'alt_base':         'C',
        'antigen_axis':     'expression',
        'ref_call':         'expressed',
        'alt_call':         'silenced',
        'role':             'modifier',
        'significance':     ('p.Ser291Pro — Asian / Finnish Jk(a-b-) point '
                             'mutation; protein is non-functional'),
        'reference':        'ISBT 009 (JK*01N.06 / JK*02N.05)',
    },
}

# Names by which the SLC14A1 (Kidd) gene may be tagged in a GenBank record.
JK_GENE_NAMES = ('SLC14A1', 'JK', 'JK1', 'UT-B', 'HUT11')


# ─── IUPAC decode (same convention as fy / kel / rhce / rhd analyzers) ───────
IUPAC_DECODE = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC',
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG',
    'N': 'ACGT',
}


# ─── Quality gates (mirror fy_analyzer for behavioural parity) ───────────────
MIN_IDENTITY_FOR_CALLING = 85.0
SNP_LOCAL_WINDOW_HALF = 50
SNP_LOCAL_WINDOW_MIN_BASES = 30
MIN_LOCAL_IDENTITY_FOR_CALLING = 85.0
SNP_PROBE_HALF = 50
SNP_PROBE_MIN_BASES = 60
MIN_PHRED_AT_SNP = 30                    # per-base Phred Q-score required AT the
                                         # SNP column when AB1 quality is available.
                                         # Q30 = 1-in-1000 base error rate (lab std).


class JKReferenceMissingError(RuntimeError):
    """Raised when no JK reference file is available."""


class JKAnalyzer:
    """Kidd (JK) blood group analyzer covering three diagnostic SNPs.

    See module docstring for the SNP scope and the rationale for needing a
    genomic reference if the intronic c.342-1G>A SNP is to be callable.

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

        gb = Path(gb_path) if gb_path else JK_GENBANK_PATH
        fa = Path(fasta_path) if fasta_path else JK_FASTA_PATH

        if gb.exists():
            self._load_genbank(gb)
        elif fa.exists():
            self._load_fasta(fa)
        # else: stay unloaded — analyze_* raises with download guidance.

    # ─── Reference loading ───────────────────────────────────────────────
    # NG_011775.2 is currently a single-gene record (SLC14A1 only) — the
    # gene-name filter is unnecessary today but kept for parity with FY in
    # case a future RefSeqGene revision bundles a neighbour.
    @classmethod
    def _feature_is_jk_gene(cls, feature) -> bool:
        gene = feature.qualifiers.get('gene', [''])[0] or ''
        return gene in JK_GENE_NAMES

    def _load_genbank(self, path: Path) -> None:
        record = SeqIO.read(str(path), "genbank")
        self.reference_seq = str(record.seq).upper()
        self.reference_kind = 'genomic'

        exons = []
        cds_start = None
        any_jk_tagged = any(
            self._feature_is_jk_gene(f)
            for f in record.features
            if f.type in ('gene', 'CDS', 'exon', 'mRNA')
        )
        for feature in record.features:
            if any_jk_tagged and not self._feature_is_jk_gene(feature):
                continue
            if feature.type == "exon":
                start = int(feature.location.start)
                end = int(feature.location.end)
                number = feature.qualifiers.get("number", [None])[0]
                # SLC14A1 has alt-spliced exon-3 variants annotated as 3a/3b/3c;
                # parse the int prefix and skip any that don't fit the
                # canonical numbering when collisions occur.
                try:
                    number = int(str(number).split('a')[0]
                                 .split('b')[0].split('c')[0]) if number else None
                except (TypeError, ValueError):
                    number = None
                exons.append((start, end, number))
            elif feature.type == "CDS" and cds_start is None:
                cds_start = int(feature.location.start)

        # De-duplicate exon entries by (start, end) — alt-spliced variants
        # would otherwise inflate the cDNA-position walk.
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
        """Return the 0-based mRNA index of the ATG start codon for SLC14A1.

        SLC14A1 (UT-B) protein begins M-E-D-S-P-T-M-V-R (UniProt Q13336), so
        the start-codon context in the mRNA is `ATGGAGGACAGCCCC`. Falls back
        to first ATG if the specific context is not found.
        """
        jk_context = seq.find('ATGGAGGACAGCCCC')
        if jk_context >= 0:
            return jk_context
        first_atg = seq.find('ATG')
        return first_atg if first_atg >= 0 else 0

    def is_loaded(self) -> bool:
        return self.reference_seq is not None

    def _require_reference(self) -> None:
        if not self.is_loaded():
            raise JKReferenceMissingError(
                "JK reference not found. Supply one of:\n"
                f"  - GenBank (genomic, REQUIRED for c.342-1 splice SNP): {JK_GENBANK_PATH}\n"
                f"  - FASTA (cDNA, c.342-1 will be uncoverable):  {JK_FASTA_PATH}\n"
                "See utils/data/JK_REFERENCE_README.md for download instructions."
            )

    # ─── Position translation ────────────────────────────────────────────
    def _cdna_pos_to_ref_index(self, cdna_pos: int) -> Optional[int]:
        """Translate a 1-based cDNA position to a 0-based reference index."""
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

    def _position_to_ref_index(self, snp_info: dict) -> Optional[int]:
        """Map a SNP's annotated position to a 0-based reference index.

        Handles three position types: 'cdna' (KEL/RHCE-style), 'intronic'
        (cDNA-anchored offset, requires genomic reference), and 'promoter'
        (cds-anchored offset, requires genomic reference; kept for symmetry
        with FYAnalyzer).
        """
        ptype = snp_info.get('position_type', 'cdna')

        if ptype == 'intronic':
            if self.reference_kind != 'genomic':
                return None
            anchor_idx = self._cdna_pos_to_ref_index(int(snp_info['cDNA_anchor']))
            if anchor_idx is None:
                return None
            # SLC14A1 is on the + strand of NG_011775.2. For - strand genes
            # this would need to be reversed; revisit when needed.
            idx = anchor_idx + int(snp_info['intronic_offset'])
            if 0 <= idx < len(self.reference_seq):
                return idx
            return None

        if ptype == 'promoter':
            if self.reference_kind != 'genomic' or self.cds_start_genomic is None:
                return None
            idx = self.cds_start_genomic + int(snp_info['promoter_offset'])
            if 0 <= idx < len(self.reference_seq):
                return idx
            return None

        return self._cdna_pos_to_ref_index(int(snp_info['cDNA_position']))

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
            'exon':          snp_info.get('exon'),
            'cDNA_position': snp_info.get('cDNA_position') or snp_info.get('cDNA_anchor'),
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
            ptype = snp_info.get('position_type', 'cdna')
            if ptype == 'intronic':
                return {
                    'covered': False,
                    'reason': ('intronic SNP requires a genomic reference; '
                               'this run is using a cDNA-only reference'),
                }
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

        Skips the full-read alignment against the 35 kb genomic reference
        that earlier versions used as a coarse identity gate. The per-SNP
        probe gate (MIN_LOCAL_IDENTITY_FOR_CALLING over a ±50 bp probe) is
        the actual calling criterion; the full-read alignment added nothing
        beyond perf cost — synthetic-patch tests passing the whole 35 kb
        reference back as the read were O(N²) per call.

        identity / strand are derived from the strongest probe call —
        same shape, same UI, ~100× faster on long synthetic inputs.
        """
        self._require_reference()

        query_str = str(query_seq).upper()
        query_length = len(query_str)

        snp_calls: dict = {
            snp_name: self._genotype_snp_via_probe(
                query_str, snp_info, query_quality=query_quality,
            )
            for snp_name, snp_info in JK_DIAGNOSTIC_SNPS.items()
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
        a_b_axis = self._call_antigen_axis(snp_consensus, 'A/B')
        null_state = self._call_combined_null(snp_consensus)
        phenotype_block = self._phenotype_and_alleles(a_b_axis, null_state)

        callable_reads = [r for r in per_read if r['callable']]
        reads_with_trusted_snps = [r for r in per_read if r.get('trusted_snps', 0) > 0]
        no_calls = [
            snp for snp, cons in snp_consensus.items()
            if cons['consensus'] == 'no_call' and JK_DIAGNOSTIC_SNPS[snp]['role'] == 'primary'
        ]

        return {
            'phenotype':               phenotype_block['phenotype'],
            'allele_options':          phenotype_block['allele_options'],
            'a_b_call':                a_b_axis,
            'null_call':               null_state,
            'snp_consensus':           snp_consensus,
            'no_calls_primary':        no_calls,
            'overall_confidence':      a_b_axis['confidence'],
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
        """Vote on each diagnostic SNP across reads. Identical rules to
        FY/KEL — see those analyzers for the LOW/MEDIUM/HIGH gates."""
        consensus: dict = {}

        for snp_name in JK_DIAGNOSTIC_SNPS:
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
    def _call_antigen_axis(snp_consensus: dict, axis: str) -> dict:
        """Combine the primary SNP on one axis into a single genotype call.

        Parity-safe: derives the genotype string from the SNP entry's
        `ref_call` / `alt_call` rather than assuming hom_ref == 'AA'. This
        is required because the c.838G>A entry is parity-inverted (the
        deposited NG_011775.2 reference carries JK*B, not JK*A).
        """
        primary_snps = [
            name for name, info in JK_DIAGNOSTIC_SNPS.items()
            if info['antigen_axis'] == axis and info['role'] == 'primary'
        ]
        primary_name = primary_snps[0] if primary_snps else None
        primary = snp_consensus.get(primary_name) if primary_name else None

        if not primary or primary['consensus'] == 'no_call':
            return {
                'axis':            axis,
                'genotype':        None,
                'confidence':      'NONE',
                'primary_snp':     primary_name,
                'supporting_snps': [],
                'discordant_snps': [],
                'reason':          f'Primary marker {primary_name} not callable',
            }

        info = JK_DIAGNOSTIC_SNPS[primary_name]
        ref_call = info['ref_call']
        alt_call = info['alt_call']
        # Het is always sorted alphabetically so 'AB' is the single het
        # genotype string regardless of which allele the reference carries.
        genotype = {
            'hom_ref': ref_call * 2,
            'het':     ''.join(sorted([ref_call, alt_call])),
            'hom_alt': alt_call * 2,
        }[primary['consensus']]

        return {
            'axis':            axis,
            'genotype':        genotype,
            'confidence':      primary['confidence'],
            'primary_snp':     primary_name,
            'supporting_snps': [primary_name],
            'discordant_snps': [],
            'reason':          f'Primary {primary_name} = {primary["consensus"]} ({primary["confidence"]})',
        }

    @staticmethod
    def _call_combined_null(snp_consensus: dict) -> dict:
        """Combine the two independent JK-null SNPs into a single null state.

        - 'both_silenced' : EITHER null SNP is hom_alt (one variant alone is
                            sufficient to produce Jk(a-b-) when homozygous).
        - 'one_silenced'  : EITHER null SNP is heterozygous (one haplotype
                            silenced; phase ambiguity if A/B is also het).
        - 'no_silencing'  : both null SNPs are hom_ref (or both uncallable).

        When both nulls are heterozygous on the same sample (very rare —
        compound silencing), confidence is dropped to LOW and the trigger
        list records both.
        """
        triggers = []
        worst_confidence = None
        any_hom_alt = False
        any_het = False

        for snp_name, info in JK_DIAGNOSTIC_SNPS.items():
            if info['antigen_axis'] != 'expression':
                continue
            cons = snp_consensus.get(snp_name)
            if not cons or cons['consensus'] in (None, 'no_call'):
                continue
            if cons['consensus'] == 'hom_alt':
                any_hom_alt = True
                triggers.append({'snp': snp_name, 'state': 'hom_alt',
                                 'significance': info['significance']})
            elif cons['consensus'] == 'het':
                any_het = True
                triggers.append({'snp': snp_name, 'state': 'het',
                                 'significance': info['significance']})
            # Track the lowest confidence among contributing nulls.
            if cons['consensus'] in ('hom_alt', 'het'):
                if (worst_confidence is None
                        or cons['confidence'] in ('LOW', 'NONE')):
                    worst_confidence = cons['confidence']

        if any_hom_alt:
            state = 'both_silenced'
        elif any_het:
            state = 'one_silenced'
        else:
            state = 'no_silencing'

        # Compound het across both nulls is unusual — hedge confidence.
        compound_het = sum(1 for t in triggers if t['state'] == 'het') >= 2
        if compound_het and worst_confidence != 'NONE':
            worst_confidence = 'LOW'

        return {
            'state':         state,
            'triggers':      triggers,
            'compound_het':  compound_het,
            'confidence':    worst_confidence or 'NONE',
        }

    @staticmethod
    def _phenotype_and_alleles(a_b_axis: dict, null_state: dict) -> dict:
        """Combine A/B genotype + combined null state into a phenotype label
        and ISBT haplotype option(s).

        Phase-ambiguity handling: when the null state is 'one_silenced' AND
        A/B is heterozygous, the silencing variant could be in cis with JK*A
        or with JK*B — Sanger genotyping cannot phase the two. Both options
        are returned.
        """
        ab_geno = a_b_axis['genotype']

        if ab_geno is None:
            return {
                'phenotype':      'Indeterminate',
                'allele_options': [],
                'reason':         f"A/B: {a_b_axis['reason']}",
            }

        nstate = null_state['state']
        triggers = null_state.get('triggers', [])
        trigger_names = ', '.join(t['snp'] for t in triggers) or '-'

        # 1) Both haplotypes silenced — Jk(a-b-) regardless of A/B genotype.
        if nstate == 'both_silenced':
            return {
                'phenotype':      'Jk(a-b-)',
                'allele_options': [{
                    'haplotypes': 'JK*0/JK*0',
                    'isbt':       'JK*null / JK*null',
                    'serology':   ('Jk(a-b-) — both haplotypes silenced '
                                   f'(trigger: {trigger_names})'),
                }],
                'reason': (f'Null variant(s) homozygous: {trigger_names}; '
                           'both haplotypes silenced regardless of A/B genotype.'),
            }

        # 2) No silencing (or null SNPs uncallable) — A/B determines phenotype.
        if nstate == 'no_silencing':
            ab_map = {
                'AA': ('JK*A', 'JK*A', 'JK*01/JK*01', 'Jk(a+b-)'),
                'AB': ('JK*A', 'JK*B', 'JK*01/JK*02', 'Jk(a+b+)'),
                'BB': ('JK*B', 'JK*B', 'JK*02/JK*02', 'Jk(a-b+)'),
            }
            a1, a2, isbt, serology = ab_map[ab_geno]
            return {
                'phenotype':      serology,
                'allele_options': [{
                    'haplotypes': f'{a1}/{a2}',
                    'isbt':       isbt,
                    'serology':   serology,
                }],
                'reason': (f'A/B = {ab_geno} ({a_b_axis["confidence"]}); '
                           'no silencing variants detected.'),
            }

        # 3) One haplotype silenced — phase ambiguous when A/B is het.
        if ab_geno == 'AA':
            return {
                'phenotype':      'Jk(a+b-)',
                'allele_options': [{
                    'haplotypes': f'JK*A/JK*A,null({trigger_names})',
                    'isbt':       'JK*01 / JK*01N',
                    'serology':   ('Jk(a+b-) — one JK*A silenced; net '
                                   'Jka with reduced copy number'),
                }],
                'reason': ('A/B homozygous JK*A; one null variant heterozygous '
                           f'({trigger_names}) — one haplotype silenced.'),
            }

        if ab_geno == 'BB':
            return {
                'phenotype':      'Jk(a-b+)',
                'allele_options': [{
                    'haplotypes': f'JK*B/JK*B,null({trigger_names})',
                    'isbt':       'JK*02 / JK*02N',
                    'serology':   ('Jk(a-b+) — one JK*B silenced; net '
                                   'Jkb with reduced copy number'),
                }],
                'reason': ('A/B homozygous JK*B; one null variant heterozygous '
                           f'({trigger_names}) — one haplotype silenced.'),
            }

        # ab_geno == 'AB' with a null het -> phase ambiguous
        return {
            'phenotype':      'Jk(a+b-) or Jk(a-b+) (phase ambiguous)',
            'allele_options': [
                {
                    'haplotypes': f'JK*A,null({trigger_names}) / JK*B',
                    'isbt':       'JK*01N / JK*02',
                    'serology':   'Jk(a-b+) — null on JK*A haplotype',
                },
                {
                    'haplotypes': f'JK*A / JK*B,null({trigger_names})',
                    'isbt':       'JK*01 / JK*02N',
                    'serology':   'Jk(a+b-) — null on JK*B haplotype',
                },
            ],
            'reason': ('A/B heterozygous and a null variant heterozygous '
                       f'({trigger_names}): the null could be in cis with '
                       'JK*A or JK*B. Sanger genotyping cannot phase the '
                       'two — both options listed.'),
        }
