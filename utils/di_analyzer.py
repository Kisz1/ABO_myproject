"""
DI Analyzer - ISBT-Informed Diego Blood Group (Di(a)/Di(b)) Analysis

References:
- ISBT 010 (DI system) - SLC4A1 gene (also AE1 / Band 3)
- NM_000342.4 (SLC4A1 mRNA, cDNA reference)
- NG_007498.1 (SLC4A1 genomic, RefSeqGene LRG_803)

Diagnostic SNP scope (Phase 1):
    Di(a)/Di(b):  c.2561T>C (p.Pro854Leu) — primary discriminator
                  Di(a)  = T at c.2561 (Leu854); ~10% in indigenous Americans
                                                  and East Asians
                  Di(b)  = C at c.2561 (Pro854); ~99% worldwide

Located in SLC4A1 exon 19. The lab amplicon `DI1819` (exons 18-19)
covers this position and is the standard Phase-1 genotyping target.

Phase 2 (deferred — different exon, no patient amplicon coverage yet):
    Wr(a)/Wr(b): c.1972G>A (p.Glu658Lys) — exon 16
                 Wr(a) extremely rare; Wr(b) common.

----------------------------------------------------------------------------
Reference-file audit note
----------------------------------------------------------------------------
NG_007498.1 deposits the Di-b allele (`C` at c.2561, codon 854 = CCG = Pro)
— the more common allele globally. So the SNP entry below is
*parity-inverted*: the file's "ref" base is `C` (Di-b) and the file's
"alt" base is `T` (Di-a), even though the textbook SNP name is
"c.2561T>C". Same trick as JK c.838 and MNS GYPB c.143.

The regression test
`tests/test_di_analyzer.py::test_file_base_at_c2561_matches_audit_note`
enforces this. If a future NCBI revision flips to a Di-a clone (`T` at
c.2561), un-invert the entry and remove `parity_inverted: True`.

----------------------------------------------------------------------------
Architecture
----------------------------------------------------------------------------
Single-SNP single-gene system — closest in shape to KELAnalyzer. Probe-based
SNP calling (same pattern as RHCE / KEL / FY / JK / H / MNS); the per-SNP
probe is anchored in SLC4A1 exon-19 context, and the probe gate is the
sole "callable" criterion (matching the MNS optimisation).
"""

from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq


# ─── Reference file locations ────────────────────────────────────────────────
DATA_DIR = Path(__file__).resolve().parent / "data"
DI_GENBANK_PATH = DATA_DIR / "di_referance.gb"
DI_FASTA_PATH = DATA_DIR / "di_cdna_reference.fasta"


# ─── ISBT-defined diagnostic SNPs ────────────────────────────────────────────
DI_DIAGNOSTIC_SNPS = {
    'c.2561T>C': {
        'position_type':    'cdna',
        'cDNA_position':    2561,
        'exon':             19,
        # Parity inversion: NG_007498.1 deposits the Di-b allele (Pro854,
        # base 'C') at c.2561. Textbook SNP name keeps the rare-first
        # convention "c.2561T>C", but the file's ref base is C and the
        # file's variant base is T. Same pattern as JK c.838 and MNS GYPB c.143.
        'ref_base':         'C',
        'alt_base':         'T',
        'antigen_axis':     'Dia/Dib',
        'ref_call':         'b',     # the file's reference base encodes Di-b
        'alt_call':         'a',     # the file's variant base encodes Di-a
        'role':             'primary',
        'significance':     'p.Pro854Leu - primary Di(a)/Di(b) discriminator',
        'reference':        'ISBT 010 (DI*01 = Di-b, DI*02 = Di-a)',
        'parity_inverted':  True,
    },
}

# Names by which the SLC4A1 (Diego) gene may be tagged in a GenBank record.
# NG_007498.1 is single-gene at time of writing; filter kept for parity with
# FY / H / MNS multi-gene loaders.
DI_GENE_NAMES = ('SLC4A1', 'DI', 'AE1', 'EPB3', 'BAND3')


# ─── IUPAC decode (same convention as other analyzers) ───────────────────────
IUPAC_DECODE = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC',
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG',
    'N': 'ACGT',
}


# ─── Quality gates (mirror jk/fy/h/mns for behavioural parity) ───────────────
MIN_IDENTITY_FOR_CALLING = 85.0
MIN_LOCAL_IDENTITY_FOR_CALLING = 85.0
SNP_PROBE_HALF = 50
SNP_PROBE_MIN_BASES = 60
MIN_PHRED_AT_SNP = 30                    # per-base Phred Q-score required AT the
                                         # SNP column when AB1 quality is available.
                                         # Q30 = 1-in-1000 base error rate (lab std).


class DIReferenceMissingError(RuntimeError):
    """Raised when no Diego (SLC4A1) reference file is available."""


class DIAnalyzer:
    """Diego (Di) blood group analyzer covering the primary Di(a)/Di(b) SNP.

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
        self.reference_kind: Optional[str] = None
        self.exon_map: list = []
        self.cds_start_genomic: Optional[int] = None
        self.cdna_offset: int = 0
        self._aligner = self._build_aligner()

        if reference_seq is not None:
            self.reference_seq = str(reference_seq).upper()
            self.reference_kind = 'cdna'
            self.cdna_offset = self._detect_cds_start_in_mrna(self.reference_seq)
            return

        gb = Path(gb_path) if gb_path else DI_GENBANK_PATH
        fa = Path(fasta_path) if fasta_path else DI_FASTA_PATH

        if gb.exists():
            self._load_genbank(gb)
        elif fa.exists():
            self._load_fasta(fa)
        # else: stay unloaded — analyze_* raises with download guidance.

    # ─── Reference loading ───────────────────────────────────────────────
    @classmethod
    def _feature_is_di_gene(cls, feature) -> bool:
        gene = feature.qualifiers.get('gene', [''])[0] or ''
        return gene in DI_GENE_NAMES

    def _load_genbank(self, path: Path) -> None:
        record = SeqIO.read(str(path), "genbank")
        self.reference_seq = str(record.seq).upper()
        self.reference_kind = 'genomic'

        exons = []
        cds_start = None
        any_di_tagged = any(
            self._feature_is_di_gene(f)
            for f in record.features
            if f.type in ('gene', 'CDS', 'exon', 'mRNA')
        )
        for feature in record.features:
            if any_di_tagged and not self._feature_is_di_gene(feature):
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
        """Return the 0-based mRNA index of the ATG start codon for SLC4A1.

        SLC4A1 (Band 3 / AE1) protein begins M-E-E-L-Q (UniProt P02730), so
        the start-codon context is `ATGGAGGAGCTGCAG`. Falls back to first ATG.
        """
        di_context = seq.find('ATGGAGGAGCTGCAG')
        if di_context >= 0:
            return di_context
        first_atg = seq.find('ATG')
        return first_atg if first_atg >= 0 else 0

    def is_loaded(self) -> bool:
        return self.reference_seq is not None

    def _require_reference(self) -> None:
        if not self.is_loaded():
            raise DIReferenceMissingError(
                "Diego (SLC4A1) reference not found. Supply one of:\n"
                f"  - GenBank: {DI_GENBANK_PATH}\n"
                f"  - FASTA:   {DI_FASTA_PATH}\n"
                "See utils/data/DI_REFERENCE_README.md for download instructions."
            )

    # ─── Position translation ────────────────────────────────────────────
    def _position_to_ref_index(self, snp_info: dict) -> Optional[int]:
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
        """Probe-only single-read analysis (no full-read alignment gate).

        Single-SNP system, so identity / strand / callable are derived from
        the one probe's call. Mirrors the MNS optimisation — full-read
        alignment against a 27 kb reference is wasteful for long synthetic
        inputs (patched references in tests) and adds nothing beyond the
        probe gate's filtering.

        ``query_quality`` enables the Phred-Q30-at-SNP-column gate (lab std).
        """
        self._require_reference()

        query_str = str(query_seq).upper()
        query_length = len(query_str)

        snp_calls: dict = {
            snp_name: self._genotype_snp_via_probe(
                query_str, snp_info, query_quality=query_quality,
            )
            for snp_name, snp_info in DI_DIAGNOSTIC_SNPS.items()
        }

        # Identity / strand from the (single) probe call.
        primary = next(iter(snp_calls.values()))
        identity = primary.get('local_window_identity', 0.0) or 0.0
        strand = primary.get('probe_strand')

        trusted_snps = sum(
            1 for c in snp_calls.values()
            if c.get('covered') and c.get('call') and c.get('call') != 'no_call'
        )
        low_phred_at_snp = sum(
            1 for c in snp_calls.values()
            if (c.get('phred_at_snp') is not None
                and c.get('phred_at_snp') < MIN_PHRED_AT_SNP)
        )
        callable_read = (trusted_snps > 0
                         or identity >= MIN_IDENTITY_FOR_CALLING)

        reason = (
            f"Probe identity {identity:.1f}% ({strand}); "
            f"{trusted_snps} diagnostic SNP(s) cleared the probe-identity gate"
            + (f"; {low_phred_at_snp} blocked by Phred Q{MIN_PHRED_AT_SNP} gate"
               if low_phred_at_snp else "")
            + "."
        )

        return {
            'read_id':         read_id,
            'query_length':    query_length,
            'strand':          strand,
            'identity':        round(identity, 2),
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
        di_axis = self._call_antigen_axis(snp_consensus, 'Dia/Dib')
        phenotype_block = self._phenotype_and_alleles(di_axis)

        callable_reads = [r for r in per_read if r['callable']]
        reads_with_trusted_snps = [r for r in per_read if r.get('trusted_snps', 0) > 0]
        no_calls = [
            snp for snp, cons in snp_consensus.items()
            if cons['consensus'] == 'no_call' and DI_DIAGNOSTIC_SNPS[snp]['role'] == 'primary'
        ]

        return {
            'phenotype':               phenotype_block['phenotype'],
            'allele_options':          phenotype_block['allele_options'],
            'di_axis_call':            di_axis,
            'snp_consensus':           snp_consensus,
            'no_calls_primary':        no_calls,
            'overall_confidence':      di_axis['confidence'],
            'reads_total':             len(per_read),
            'reads_callable':          len(callable_reads),
            'reads_with_trusted_snps': len(reads_with_trusted_snps),
            'per_read_details':        per_read,
            'reason':                  phenotype_block['reason'],
        }

    @staticmethod
    def _normalize_reads(reads) -> list:
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
        gates as MNS/JK/FY/H — see those analyzers for the rule rationale."""
        consensus: dict = {}

        for snp_name in DI_DIAGNOSTIC_SNPS:
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

        Parity-safe: the genotype string is built from the SNP entry's
        `ref_call` / `alt_call` rather than assuming hom_ref == 'aa'. Plain
        alphabetical sort produces 'ab' for het (capital before lowercase
        in ASCII; the Diego entry uses lowercase 'a' and 'b' so a < b).
        """
        primary_snps = [
            name for name, info in DI_DIAGNOSTIC_SNPS.items()
            if info['antigen_axis'] == axis and info['role'] == 'primary'
        ]
        primary_name = primary_snps[0] if primary_snps else None
        primary = snp_consensus.get(primary_name) if primary_name else None

        if not primary or primary['consensus'] == 'no_call':
            return {
                'axis':        axis,
                'genotype':    None,
                'confidence':  'NONE',
                'primary_snp': primary_name,
                'reason':      f'Primary marker {primary_name} not callable',
            }

        info = DI_DIAGNOSTIC_SNPS[primary_name]
        ref_call = info['ref_call']
        alt_call = info['alt_call']
        sorted_pair = ''.join(sorted([ref_call, alt_call]))
        genotype = {
            'hom_ref': ref_call * 2,
            'het':     sorted_pair,
            'hom_alt': alt_call * 2,
        }[primary['consensus']]

        return {
            'axis':        axis,
            'genotype':    genotype,
            'confidence':  primary['confidence'],
            'primary_snp': primary_name,
            'reason':      f'Primary {primary_name} = {primary["consensus"]} ({primary["confidence"]})',
        }

    @staticmethod
    def _phenotype_and_alleles(di_axis: dict) -> dict:
        """Map the Di-axis genotype to a phenotype label and ISBT haplotype.

        Single-axis system — exactly one haplotype option, no phase ambiguity.
        """
        geno = di_axis['genotype']

        if geno is None:
            return {
                'phenotype':      'Indeterminate',
                'allele_options': [],
                'reason':         f"Di(a)/Di(b): {di_axis['reason']}",
            }

        # Genotype mapping with parity-inverted reference (ref=b, alt=a).
        # Phenotypes:
        #   bb -> Di(a-b+)
        #   ab -> Di(a+b+)
        #   aa -> Di(a+b-)
        ab_map = {
            'bb': ('Di-b', 'Di-b', 'DI*01/DI*01', 'Di(a-b+)'),
            'ab': ('Di-a', 'Di-b', 'DI*02/DI*01', 'Di(a+b+)'),
            'aa': ('Di-a', 'Di-a', 'DI*02/DI*02', 'Di(a+b-)'),
        }
        a1, a2, isbt, serology = ab_map[geno]

        return {
            'phenotype':      serology,
            'allele_options': [{
                'haplotypes': f'{a1}/{a2}',
                'isbt':       isbt,
                'serology':   serology,
            }],
            'reason': (f'Di(a)/Di(b) genotype = {geno} ({di_axis["confidence"]}); '
                       f'serological equivalent {serology}.'),
        }
