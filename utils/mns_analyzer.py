"""
MNS Analyzer - ISBT-Informed MNS Blood Group (GYPA + GYPB) Analysis

References:
- ISBT 002 (MNS system) - GYPA (M/N antigens) + GYPB (S/s antigens)
- NM_002099.7 / NG_007470.3 (GYPA, glycophorin A)
- NM_002100.5 / NG_007483.3 (GYPB, glycophorin B)

Diagnostic SNP scope (Phase 1, 2 SNPs across 2 genes):
    M / N:   GYPA c.59C>T  (p.Ser20Leu) — primary M/N discriminator
    S / s:   GYPB c.143C>T (p.Thr48Met) — primary S/s discriminator
             *parity-inverted*: NG_007483.3 deposits the s allele

Why only 2 SNPs in Phase 1:
    The M/N locus also has changes at codon 24 (positions c.71+c.72), but
    those require haplotype-pair calling — c.71G>A *alone* against the
    GYPA reference (codon GGT) gives Asp, not Glu, because the N allele
    actually carries c.71G>A AND c.72T>A together (GGT->GAA). Single-SNP
    calling at c.71 would be wrong, so codon 24 is deferred to Phase 2
    where a haplotype-aware caller can be added.

----------------------------------------------------------------------------
Reference-file audit notes
----------------------------------------------------------------------------
GYPA (NG_007470.3): SHOULD carry the M allele base 'C' at c.59. This is
the assumption in MNS_DIAGNOSTIC_SNPS. The regression test
`test_file_base_at_c59_matches_audit_note` enforces this.

GYPB (NG_007483.3): SHOULD carry the s allele base 'C' at c.143. The
SNP entry is parity-inverted (ref_call='s', alt_call='S') because the
deposit is the s allele. If a future revision flips to the S allele
('T' at c.143), un-invert the entry: swap ref_base/alt_base and
ref_call/alt_call, and remove `parity_inverted: True`. Same pattern
as JK c.838 and RHCE c.178.

----------------------------------------------------------------------------
Architecture
----------------------------------------------------------------------------
Dual-gene system: the analyzer loads BOTH GYPA and GYPB references in
__init__ and tracks per-gene CDS / exon maps in `self.gene_refs`. Each
diagnostic SNP carries a `gene` field; `_position_to_ref_index` and
`_genotype_snp_via_probe` dispatch to the correct reference per SNP.

Probe-based per-SNP calling (same pattern as RHCE / KEL / FY / JK / H)
mitigates the GYPA/GYPB ~97% paralogy: each ±50 bp probe is anchored
in gene-specific context (UTR-adjacent for c.59, intron-3 boundary for
c.143), so cross-paralog mis-alignment is unlikely.
"""

from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq


# ─── Reference file locations ────────────────────────────────────────────────
DATA_DIR = Path(__file__).resolve().parent / "data"
GYPA_GENBANK_PATH = DATA_DIR / "mns_gypa_referance.gb"
GYPB_GENBANK_PATH = DATA_DIR / "mns_gypb_referance.gb"
# (no FASTA fallback for MNS; both SNPs need the genomic coordinates)


# ─── ISBT-defined diagnostic SNPs ────────────────────────────────────────────
# `gene` field is mandatory and selects which gene's reference to use for
# position translation and probe slicing.
MNS_DIAGNOSTIC_SNPS = {
    'c.59C>T (GYPA)': {
        'gene':             'GYPA',
        'position_type':    'cdna',
        'cDNA_position':    59,
        'exon':             2,
        'ref_base':         'C',
        'alt_base':         'T',
        'antigen_axis':     'M/N',
        'ref_call':         'M',
        'alt_call':         'N',
        'role':             'primary',
        'significance':     'p.Ser20Leu - primary M/N discriminator',
        'reference':        'ISBT 002 (GYPA*M / GYPA*N)',
    },
    'c.143C>T (GYPB)': {
        'gene':             'GYPB',
        'position_type':    'cdna',
        'cDNA_position':    143,
        'exon':             3,
        # Parity inversion: NG_007483.3 deposits the s allele (Thr48, ACG)
        # at this position. The textbook SNP name is "c.143T>C (S>s)" but
        # the file's ref base is C and the file's variant base is T.
        # Same trick as JK c.838 and RHCE c.178.
        'ref_base':         'C',
        'alt_base':         'T',
        'antigen_axis':     'S/s',
        'ref_call':         's',  # the file's reference base encodes s
        'alt_call':         'S',  # the file's variant base encodes S (Met48)
        'role':             'primary',
        'significance':     'p.Thr48Met - primary S/s discriminator',
        'reference':        'ISBT 002 (GYPB*S / GYPB*s)',
        'parity_inverted':  True,
    },
}

# Names by which each MNS gene may be tagged in a GenBank record. Both
# RefSeqGene records are single-gene at time of writing, but the filter
# is kept for parity with FY/H multi-gene loaders.
GYPA_GENE_NAMES = ('GYPA', 'MN', 'GPA')
GYPB_GENE_NAMES = ('GYPB', 'SS', 'GPB')


# ─── IUPAC decode (same convention as other analyzers) ───────────────────────
IUPAC_DECODE = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC',
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG',
    'N': 'ACGT',
}


# ─── Quality gates (mirror jk/fy/h/kel for behavioural parity) ───────────────
MIN_IDENTITY_FOR_CALLING = 85.0
SNP_LOCAL_WINDOW_HALF = 50
SNP_LOCAL_WINDOW_MIN_BASES = 30
MIN_LOCAL_IDENTITY_FOR_CALLING = 85.0
SNP_PROBE_HALF = 50
SNP_PROBE_MIN_BASES = 60
MIN_PHRED_AT_SNP = 30                    # per-base Phred Q-score required AT the
                                         # SNP column when AB1 quality is available.
                                         # Q30 = 1-in-1000 base error rate (lab std).


class MNSReferenceMissingError(RuntimeError):
    """Raised when one or both of the GYPA / GYPB reference files are absent."""


class MNSAnalyzer:
    """MNS blood group analyzer covering M/N (GYPA) and S/s (GYPB) axes.

    Phase 1 covers two SNPs across two genes. Codon 24 (M/N confirmatory),
    U antigen detection, and Miltenberger hybrid alleles are out of scope
    (would need haplotype-pair calling and/or hybrid-breakpoint detection).

    Use ``analyze_single_read(query_seq)`` for one read or ``analyze(reads)``
    for multi-read consensus voting.
    """

    def __init__(
        self,
        gypa_gb_path: Optional[str] = None,
        gypb_gb_path: Optional[str] = None,
    ):
        self.gene_refs: dict = {}   # {'GYPA': {seq, kind, exon_map, cds_start, ...}}
        self._aligner = self._build_aligner()

        gypa = Path(gypa_gb_path) if gypa_gb_path else GYPA_GENBANK_PATH
        gypb = Path(gypb_gb_path) if gypb_gb_path else GYPB_GENBANK_PATH

        if gypa.exists():
            self._load_genbank(gypa, 'GYPA', GYPA_GENE_NAMES)
        if gypb.exists():
            self._load_genbank(gypb, 'GYPB', GYPB_GENE_NAMES)
        # else: stay unloaded — analyze_* raises with download guidance.

    # ─── Reference loading ───────────────────────────────────────────────
    @staticmethod
    def _feature_is_target_gene(feature, gene_names) -> bool:
        gene = feature.qualifiers.get('gene', [''])[0] or ''
        return gene in gene_names

    def _load_genbank(self, path: Path, gene_key: str, gene_names) -> None:
        record = SeqIO.read(str(path), "genbank")
        seq = str(record.seq).upper()

        exons = []
        cds_start = None
        any_target_tagged = any(
            self._feature_is_target_gene(f, gene_names)
            for f in record.features
            if f.type in ('gene', 'CDS', 'exon', 'mRNA')
        )
        for feature in record.features:
            if any_target_tagged and not self._feature_is_target_gene(feature, gene_names):
                continue
            if feature.type == "exon":
                start = int(feature.location.start)
                end = int(feature.location.end)
                number = feature.qualifiers.get("number", [None])[0]
                try:
                    import re
                    m = re.match(r'^\d+', str(number) if number is not None else '')
                    number = int(m.group()) if m else None
                except (TypeError, ValueError):
                    number = None
                exons.append((start, end, number))
            elif feature.type == "CDS" and cds_start is None:
                cds_start = int(feature.location.start)

        # Dedup by (start, end) — alt-spliced exon variants would otherwise
        # bloat the cDNA-position walk.
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

        self.gene_refs[gene_key] = {
            'seq':                seq,
            'kind':               'genomic',
            'cds_start_genomic':  cds_start,
            'exon_map':           deduped,
        }

    def is_loaded(self) -> bool:
        # Both genes must be loaded for MNS analysis to be meaningful.
        return ('GYPA' in self.gene_refs and 'GYPB' in self.gene_refs)

    # Backward-compat properties so simple callers can probe a "primary"
    # reference (use GYPA — that's the M/N gene people usually mean first).
    @property
    def reference_seq(self) -> Optional[str]:
        return self.gene_refs.get('GYPA', {}).get('seq')

    @property
    def reference_kind(self) -> Optional[str]:
        return self.gene_refs.get('GYPA', {}).get('kind')

    def _require_reference(self) -> None:
        if not self.is_loaded():
            missing = []
            if 'GYPA' not in self.gene_refs:
                missing.append(f'GYPA: {GYPA_GENBANK_PATH}')
            if 'GYPB' not in self.gene_refs:
                missing.append(f'GYPB: {GYPB_GENBANK_PATH}')
            raise MNSReferenceMissingError(
                "MNS reference not found. Both GYPA and GYPB are required.\n"
                "Missing: " + "; ".join(missing) + "\n"
                "See utils/data/MNS_REFERENCE_README.md for download instructions."
            )

    # ─── Position translation ────────────────────────────────────────────
    def _position_to_ref_index(self, snp_info: dict) -> Optional[int]:
        gene = snp_info['gene']
        gene_data = self.gene_refs.get(gene)
        if not gene_data:
            return None

        cdna_pos = int(snp_info['cDNA_position'])
        cds_start = gene_data.get('cds_start_genomic')
        exon_map = gene_data.get('exon_map') or []

        remaining = cdna_pos
        for start, end, _ in exon_map:
            exon_len = end - start
            exon_cds_start = start
            exon_cds_len = exon_len

            if cds_start is not None and start < cds_start < end:
                exon_cds_start = cds_start
                exon_cds_len = end - cds_start
            elif cds_start is not None and end <= cds_start:
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

    def _align_against_gene(self, query_seq: str, gene: str, strand: str = 'forward'):
        ref_seq = self.gene_refs[gene]['seq']
        seq = query_seq if strand == 'forward' else str(Seq(query_seq).reverse_complement())
        try:
            alignments = self._aligner.align(ref_seq, seq)
            best = alignments[0]
        except (ValueError, IndexError, OverflowError):
            return None
        return {
            'identity': self._identity(best),
            'score':    best.score,
            'strand':   strand,
            'gene':     gene,
            'alignment': best,
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
            'gene':          snp_info['gene'],
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
        """Per-SNP probe-based calling. Probe is sliced from the SNP's gene
        reference (GYPA or GYPB), so paralog mis-alignment is mitigated."""
        gene = snp_info['gene']
        gene_data = self.gene_refs.get(gene)
        if not gene_data:
            return {'covered': False, 'reason': f'reference for {gene} not loaded'}

        ref_index = self._position_to_ref_index(snp_info)
        if ref_index is None:
            return {'covered': False, 'reason': 'position outside reference'}

        ref_seq = gene_data['seq']
        lo = max(0, ref_index - SNP_PROBE_HALF)
        hi = min(len(ref_seq), ref_index + SNP_PROBE_HALF + 1)
        probe = ref_seq[lo:hi]
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
        """Analyze one read using probe-based per-SNP calling only.

        Skips the full-read alignment against each gene reference that the
        other analyzers use as a coarse "callable" gate. With dual ~38 kb
        references, full-read alignment of long sequences (e.g. patched
        full-gene synthetic inputs) was O(N×M) per read per gene — minutes
        per call. The per-SNP probe gate (MIN_LOCAL_IDENTITY_FOR_CALLING,
        default 85% over a ±50 bp window) already filters out off-target
        reads, so the coarse gate adds nothing beyond performance cost.

        identity / strand / best_gene are derived from the probe-call with
        the highest local-window identity instead — same shape, same UI.

        ``query_quality`` enables the Phred-Q30-at-SNP-column gate.
        """
        self._require_reference()

        query_str = str(query_seq).upper()
        query_length = len(query_str)

        snp_calls: dict = {
            snp_name: self._genotype_snp_via_probe(
                query_str, snp_info, query_quality=query_quality,
            )
            for snp_name, snp_info in MNS_DIAGNOSTIC_SNPS.items()
        }

        # Derive identity / strand / best_gene from the strongest probe call,
        # and per-gene best identity across that gene's SNPs.
        best = None
        per_gene_identity: dict = {}
        for snp_name, call in snp_calls.items():
            ident = call.get('local_window_identity')
            if ident is None:
                continue
            gene = MNS_DIAGNOSTIC_SNPS[snp_name]['gene']
            if ident > per_gene_identity.get(gene, -1.0):
                per_gene_identity[gene] = round(ident, 2)
            if best is None or ident > best['identity']:
                best = {
                    'identity': ident,
                    'gene':     gene,
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
                'read_id':           read_id,
                'query_length':      query_length,
                'strand':            None,
                'identity':          0.0,
                'best_gene':         None,
                'per_gene_identity': {},
                'snp_calls':         snp_calls,
                'callable':          False,
                'trusted_snps':      0,
                'low_phred_at_snp':  low_phred_at_snp,
                'reason':            'No probe alignment produced for any SNP',
            }

        # callable = the read carried at least one SNP through the probe gate
        # OR has high overall probe identity (a relevant but maybe-uncalled read).
        callable_read = (trusted_snps > 0
                         or best['identity'] >= MIN_IDENTITY_FOR_CALLING)

        per_gene_str = ", ".join(f"{g}={i}%" for g, i in per_gene_identity.items())
        reason = (
            f"Best probe identity {best['identity']:.1f}% on {best['gene']} "
            f"({best['strand']}); per-gene best {per_gene_str}; "
            f"{trusted_snps} diagnostic SNP(s) cleared the probe-identity gate"
            + (f"; {low_phred_at_snp} blocked by Phred Q{MIN_PHRED_AT_SNP} gate"
               if low_phred_at_snp else "")
            + "."
        )

        return {
            'read_id':           read_id,
            'query_length':      query_length,
            'strand':            best['strand'],
            'identity':          round(best['identity'], 2),
            'best_gene':         best['gene'],
            'per_gene_identity': per_gene_identity,
            'snp_calls':         snp_calls,
            'callable':          callable_read,
            'trusted_snps':      trusted_snps,
            'low_phred_at_snp':  low_phred_at_snp,
            'reason':            reason,
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
        m_n_axis = self._call_antigen_axis(snp_consensus, 'M/N')
        s_s_axis = self._call_antigen_axis(snp_consensus, 'S/s')
        phenotype_block = self._phenotype_and_alleles(m_n_axis, s_s_axis)

        callable_reads = [r for r in per_read if r['callable']]
        reads_with_trusted_snps = [r for r in per_read if r.get('trusted_snps', 0) > 0]
        no_calls = [
            snp for snp, cons in snp_consensus.items()
            if cons['consensus'] == 'no_call' and MNS_DIAGNOSTIC_SNPS[snp]['role'] == 'primary'
        ]

        # Overall confidence = lower of the two axis confidences (NONE>LOW>MEDIUM>HIGH).
        order = ['NONE', 'LOW', 'MEDIUM', 'HIGH']
        confs = [c['confidence'] for c in (m_n_axis, s_s_axis) if c.get('confidence')]
        overall = min(confs, key=order.index) if confs else 'NONE'

        return {
            'phenotype':               phenotype_block['phenotype'],
            'allele_options':          phenotype_block['allele_options'],
            'm_n_call':                m_n_axis,
            's_s_call':                s_s_axis,
            'snp_consensus':           snp_consensus,
            'no_calls_primary':        no_calls,
            'overall_confidence':      overall,
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
        gates as JK/FY/H — see those analyzers for the rule rationale."""
        consensus: dict = {}

        for snp_name in MNS_DIAGNOSTIC_SNPS:
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
        is required because the GYPB c.143 entry is parity-inverted (the
        deposited reference carries s, not S).
        """
        primary_snps = [
            name for name, info in MNS_DIAGNOSTIC_SNPS.items()
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
                'reason':          f'Primary marker {primary_name} not callable',
            }

        info = MNS_DIAGNOSTIC_SNPS[primary_name]
        ref_call = info['ref_call']
        alt_call = info['alt_call']
        # Het uses plain alphabetical sort so the canonical genotype string
        # is 'MN' (uppercase) and 'Ss' (capital first, then lowercase —
        # ASCII puts 'S' < 's'). A `key=str.upper` sort would collapse 'S'
        # and 's' to equal-keys, and stable sort would then preserve the
        # parity-inverted input order [s, S] -> 'sS' (wrong).
        sorted_pair = ''.join(sorted([ref_call, alt_call]))
        genotype = {
            'hom_ref': ref_call * 2,
            'het':     sorted_pair,
            'hom_alt': alt_call * 2,
        }[primary['consensus']]

        return {
            'axis':            axis,
            'genotype':        genotype,
            'confidence':      primary['confidence'],
            'primary_snp':     primary_name,
            'reason':          f'Primary {primary_name} = {primary["consensus"]} ({primary["confidence"]})',
        }

    @staticmethod
    def _phenotype_and_alleles(m_n_axis: dict, s_s_axis: dict) -> dict:
        """Combine M/N + S/s genotypes into MNS phenotype + ISBT haplotype(s).

        Phase-ambiguity handling: when BOTH axes are heterozygous (MN + Ss),
        the four MNS haplotypes (MS, Ms, NS, Ns) can pair in two distinct
        ways: MS/Ns or Ms/NS. Sanger genotyping cannot phase the two —
        both options are returned.
        """
        mn = m_n_axis['genotype']
        ss = s_s_axis['genotype']

        if mn is None and ss is None:
            return {
                'phenotype':      'Indeterminate',
                'allele_options': [],
                'reason':         f"M/N: {m_n_axis['reason']}; S/s: {s_s_axis['reason']}",
            }

        # Build the serology phenotype label (e.g. "M+N-S+s+")
        m_pos = mn is not None and 'M' in mn
        n_pos = mn is not None and 'N' in mn
        s_pos = ss is not None and 'S' in ss
        s_low_pos = ss is not None and 's' in ss

        mn_label = (f"M{'+' if m_pos else '-'}N{'+' if n_pos else '-'}"
                    if mn else "M?N?")
        ss_label = (f"S{'+' if s_pos else '-'}s{'+' if s_low_pos else '-'}"
                    if ss else "S?s?")
        phenotype = f"{mn_label}{ss_label}"

        # Determine haplotype options
        if mn is None or ss is None:
            # One axis missing — emit one option with the known half
            haps = []
            if mn and ss:
                pass
            else:
                known_axis = mn or ss
                haps.append({
                    'haplotypes': f'{known_axis} (other axis indeterminate)',
                    'isbt':       f'GYPA*{known_axis}' if mn else f'GYPB*{known_axis}',
                    'serology':   phenotype,
                })
            return {
                'phenotype':      phenotype,
                'allele_options': haps,
                'reason': (f'M/N {mn or "n/a"}, S/s {ss or "n/a"}; '
                           f'one axis missing limits haplotype assignment.'),
            }

        # Both axes called. Build haplotype string for non-double-het cases.
        if mn != 'MN' or ss != 'Ss':
            # Only one axis (or neither) is heterozygous → unique haplotype pair.
            mn_alleles = (['M', 'M'] if mn == 'MM'
                          else ['N', 'N'] if mn == 'NN'
                          else ['M', 'N'])  # mn == 'MN'
            ss_alleles = (['S', 'S'] if ss == 'SS'
                          else ['s', 's'] if ss == 'ss'
                          else ['S', 's'])  # ss == 'Ss'

            hap1 = f"{mn_alleles[0]}{ss_alleles[0]}"
            hap2 = f"{mn_alleles[1]}{ss_alleles[1]}"
            return {
                'phenotype':      phenotype,
                'allele_options': [{
                    'haplotypes': f'{hap1}/{hap2}',
                    'isbt':       f'GYPA*{mn} GYPB*{ss}',
                    'serology':   phenotype,
                }],
                'reason': (f'M/N = {mn} ({m_n_axis["confidence"]}); '
                           f'S/s = {ss} ({s_s_axis["confidence"]}). '
                           'Single haplotype assignment.'),
            }

        # Double-heterozygous → two phase-ambiguous haplotype pairs.
        return {
            'phenotype':      phenotype,
            'allele_options': [
                {
                    'haplotypes': 'MS/Ns',
                    'isbt':       'GYPA*M GYPB*S / GYPA*N GYPB*s',
                    'serology':   phenotype + ' (option 1)',
                },
                {
                    'haplotypes': 'Ms/NS',
                    'isbt':       'GYPA*M GYPB*s / GYPA*N GYPB*S',
                    'serology':   phenotype + ' (option 2)',
                },
            ],
            'reason': ('Both M/N and S/s are heterozygous: the four MNS '
                       'haplotypes (MS, Ms, NS, Ns) can pair as MS/Ns or '
                       'Ms/NS. Sanger genotyping cannot phase the two — '
                       'both options listed.'),
        }
