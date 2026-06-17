"""File ingestion + sequence processing, decoupled from Streamlit UI.

Extracted from main.py so the AB1/FASTA processors are unit-testable
without importing streamlit.
"""

import io

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

import utils.ab1_analyzer as ab1_utils
import utils.FASTA_analyzer as fasta_utils


# ── Pipeline het-calling / quality defaults ──────────────────────────────────
# Stricter heterozygote ratio than the permissive analyzer default
# (ab1_analyzer.HET_RATIO_DEFAULT = 0.10): clinical SNP calling demands a
# clearer secondary peak (>=30% of the major) before reporting a 2nd allele,
# trading sensitivity for specificity. Provenance + measured stability on the
# n=9 RHD panel: docs/THRESHOLDS.md §Heterozygosity.
HET_RATIO_PIPELINE = 0.30
# Phred sliding-window end-trim defaults (quality_trim_and_mask): Q20 = 1-in-100
# base error (standard Sanger trim floor); 10 bp window smooths single dips.
QTRIM_Q_THRESHOLD = 20
QTRIM_WINDOW = 10


def _record_error(errors, file_obj, exc):
    """Append a (filename, message) failure record to ``errors`` if a list was
    provided. Lets callers surface unreadable/corrupt uploads instead of having
    them silently dropped. No-op when ``errors`` is None (back-compat)."""
    if errors is None:
        return
    name = getattr(file_obj, 'name', 'unknown')
    errors.append({'filename': name, 'error': f"{type(exc).__name__}: {exc}"})


def process_ab1_files(fwd_ab1_files, exons_ref, threshold_ratio=HET_RATIO_PIPELINE, errors=None):
    """
    Process AB1 files for RHD analysis.

    For RHD analysis with multiple amplicons:
    - Process EACH file separately (don't merge)
    - Return list of individual traces
    - Voting system uses all amplicons to determine RhD+/RhD-
    """
    ab1_service = ab1_utils.AB1Analyzer()
    results = []
    all_hets = []

    # Process each AB1 file separately for multi-amplicon analysis
    for ab1_file in fwd_ab1_files:
        try:
            trace = ab1_service.read_ab1_trace(ab1_file)
            if not trace:
                continue

            # For single file or each file in batch: reverse chromatogram
            reversed_trace = ab1_service.reverse_chromatogram(trace)
            normalized_trace = ab1_service.normalize_trace_per_channel(reversed_trace)

            # Detect heterozygotes in this file
            raw_hets = ab1_service.detect_hetero(reversed_trace, ratio=threshold_ratio)
            hets = []
            for position, top_bases in raw_hets:
                major_base, major_signal = top_bases[0]
                minor_base, minor_signal = top_bases[1]
                ratio = minor_signal / (major_signal + 1e-6)
                hets.append({
                    "position": position,
                    "ref_base": major_base,
                    "alt_base": minor_base,
                    "ratio": ratio,
                    "major_signal": major_signal,
                    "minor_signal": minor_signal,
                    "top_bases": top_bases
                })

            all_hets.extend(hets)

            # Store trace with filename for RHD analysis identification
            reversed_trace['filename'] = ab1_file.name
            results.append(reversed_trace)

        except Exception as e:
            _record_error(errors, ab1_file, e)
            continue

    if results:
        return results, all_hets if all_hets else None

    return None, None

_BASES_TO_IUPAC = {
    frozenset('AG'): 'R', frozenset('CT'): 'Y',
    frozenset('CG'): 'S', frozenset('AT'): 'W',
    frozenset('GT'): 'K', frozenset('AC'): 'M',
}

def quality_trim_and_mask(seq, qual, q_threshold=QTRIM_Q_THRESHOLD, window=QTRIM_WINDOW):
    """Sliding-window end-trim + internal-N-mask.

    Returns (left, right, masked_seq, n_masked) where:
      - ``seq[left:right]`` is the trimmed slice with internal low-Q bases
        replaced by 'N' in ``masked_seq``
      - ``n_masked`` counts the internal N substitutions
    Returns (0, 0, '', 0) if the read is entirely below threshold.
    """
    n = len(seq)
    if n == 0:
        return 0, 0, '', 0
    if len(qual) != n:
        qual = [q_threshold] * n

    w = min(window, n)

    left = 0
    while left + w <= n:
        if all(q >= q_threshold for q in qual[left:left + w]):
            break
        left += 1
    else:
        return 0, 0, '', 0  # no clean 5' window found

    right = n
    while right - w >= left:
        if all(q >= q_threshold for q in qual[right - w:right]):
            break
        right -= 1
    else:
        return 0, 0, '', 0  # no clean 3' window found

    if right <= left:
        return 0, 0, '', 0

    trimmed = list(seq[left:right])
    n_masked = 0
    for i, q in enumerate(qual[left:right]):
        if q < q_threshold:
            trimmed[i] = 'N'
            n_masked += 1

    return left, right, ''.join(trimmed), n_masked

def process_rhd_ab1_files(rhd_ab1_files, q_threshold=QTRIM_Q_THRESHOLD, window=QTRIM_WINDOW,
                          het_ratio=HET_RATIO_PIPELINE, errors=None):
    """RHD-only AB1 processor with Phred trimming + signal-based het encoding.

    Diverges from process_ab1_files (used by ABO) in three ways:
      1. Captures Phred quality scores from the AB1 file.
      2. Sliding-window end-trims and N-masks low-Q bases before SNP analysis.
      3. Uses chromatogram signal (detect_hetero) to find heterozygous positions
         and encodes them as IUPAC codes in the basecalled sequence so the
         downstream RHDAnalyzer can detect heterozygous diagnostic SNPs.

    Returns (results_list, hets_list_or_None) shaped like process_ab1_files so
    the existing call site can swap in without other changes.
    """
    if not rhd_ab1_files:
        return None, None

    ab1_service = ab1_utils.AB1Analyzer()
    results = []
    all_hets = []

    for ab1_file in rhd_ab1_files:
        try:
            trace = ab1_service.read_ab1_trace_with_quality(ab1_file)
            if not trace:
                continue

            seq = trace['seq']
            qual = list(trace.get('quality_scores') or [])
            n = len(seq)
            if n == 0:
                continue

            left, right, masked_str, n_masked = quality_trim_and_mask(
                seq, qual, q_threshold=q_threshold, window=window)
            if right <= left:
                continue  # entirely low-quality read
            trimmed_seq = list(masked_str)

            # 3. Signal-based het detection on the original (untrimmed) trace.
            #    Normalize first so detect_hetero's thresholds are meaningful.
            normalized = ab1_service.normalize_trace_per_channel(trace)
            raw_hets = ab1_service.detect_hetero(normalized, ratio=het_ratio)

            # 4. Encode hets as IUPAC in the trimmed sequence. detect_hetero
            #    yields (sample_idx, top_bases) where sample_idx is a peak
            #    position in PLOC2 sample-space. Map back to base index.
            pos_array = np.asarray(trace['pos'])
            n_hets_encoded = 0
            for sample_idx, top_bases in raw_hets:
                base_idx_arr = np.where(pos_array == sample_idx)[0]
                if len(base_idx_arr) == 0:
                    continue
                base_idx = int(base_idx_arr[0])
                trimmed_idx = base_idx - left
                if trimmed_idx < 0 or trimmed_idx >= len(trimmed_seq):
                    continue
                if trimmed_seq[trimmed_idx] == 'N':
                    continue  # quality wins over signal at masked positions
                major_base = top_bases[0][0]
                minor_base = top_bases[1][0]
                iupac = _BASES_TO_IUPAC.get(frozenset([major_base, minor_base]))
                if not iupac:
                    continue
                trimmed_seq[trimmed_idx] = iupac
                n_hets_encoded += 1
                all_hets.append({
                    'position': base_idx,
                    'ref_base': major_base,
                    'alt_base': minor_base,
                    'ratio': top_bases[1][1] / (top_bases[0][1] + 1e-6),
                    'iupac': iupac,
                    'filename': getattr(ab1_file, 'name', 'unknown'),
                })

            trace['seq'] = ''.join(trimmed_seq)
            trace['filename'] = getattr(ab1_file, 'name', 'unknown')
            # Per-base Phred quality aligned to the trimmed sequence. Lets
            # downstream analyzers apply a stricter Q-score gate at each
            # SNP column without re-doing AB1 parsing.
            trace['quality_trimmed'] = list(qual[left:right])
            trace['qc'] = {
                'q_threshold': q_threshold,
                'window': window,
                'trimmed_5p': left,
                'trimmed_3p': n - right,
                'masked_internal': n_masked,
                'het_positions_encoded': n_hets_encoded,
                'final_length': len(trace['seq']),
                'original_length': n,
            }
            results.append(trace)

        except Exception as e:
            _record_error(errors, ab1_file, e)
            continue

    if results:
        return results, (all_hets if all_hets else None)
    return None, None

def process_rhd_fasta_files(rhd_fasta_files, errors=None):
    """
    Process RHD FASTA files for RHD multi-amplicon analysis.

    Reads each FASTA file, extracts the sequence, and packages it as a
    dict trace ({'seq': ..., 'filename': ...}) so analyze_rhd_multifactor
    can process it identically to AB1 traces.
    """
    if not rhd_fasta_files:
        return None

    traces = []
    for fasta_file in rhd_fasta_files:
        try:
            content = fasta_file.read()
            if isinstance(content, bytes):
                content = content.decode('utf-8')
            fasta_text = io.StringIO(content)
            records = list(SeqIO.parse(fasta_text, "fasta"))
            if not records:
                continue
            seq = str(records[0].seq).upper()
            if len(seq) < 50:
                continue
            traces.append({
                'seq': seq,
                'filename': fasta_file.name,
                'source': 'fasta',
            })
        except Exception as e:
            _record_error(errors, fasta_file, e)
            continue

    return traces if traces else None

def process_fasta_file(fasta_file, exon_start=0, exon_end=0):
    # Convert Streamlit uploaded file to text mode for BioPython

    # Read the file content and convert to string
    fasta_content = fasta_file.read()
    if isinstance(fasta_content, bytes):
        fasta_content = fasta_content.decode('utf-8')

    # Create a text StringIO object for BioPython
    fasta_text = io.StringIO(fasta_content)

    # Parse the FASTA file
    fasta_records = list(SeqIO.parse(fasta_text, "fasta"))
    if not fasta_records:
        raise ValueError("No sequences found in FASTA file")

    # Use the first sequence
    first_record = fasta_records[0]
    service = fasta_utils.FASTAAlignmentService()
    fwd_seq = first_record.seq
    rev_seq = fwd_seq.reverse_complement()

    if exon_end > exon_start > 0:
        fwd_fasta_analysis = service.analyze_multi_exon_sequence(
            fwd_seq, list(range(exon_start, exon_end+1)))
        rev_fasta_analysis = service.analyze_multi_exon_sequence(
            rev_seq, list(range(exon_start, exon_end+1)))
    elif exon_start == 0:
        fwd_fasta_analysis = service.analyze_multi_exon_sequence(
            fwd_seq, list(range(1, 8)))
        rev_fasta_analysis = service.analyze_multi_exon_sequence(
            rev_seq, list(range(1, 8)))
    else:
        fwd_fasta_analysis = {}
        rev_fasta_analysis = {}
    strand = {"forward": fwd_fasta_analysis,
              "reverse": rev_fasta_analysis, "none": "none"}
    fwd_similarities = {}
    rev_similarities = {}
    if 'error' in fwd_fasta_analysis or 'error' in rev_fasta_analysis:
        return {}, []
    for i in fwd_fasta_analysis['exon_alignments']:
        fwd_similarities[i['exon_number']] = i['similarity']

    for i in rev_fasta_analysis['exon_alignments']:
        rev_similarities[i['exon_number']] = i['similarity']

    exon_comparison = {}
    for exon_num in fwd_similarities.keys():
        fwd_sim = fwd_similarities.get(exon_num, 0)
        rev_sim = rev_similarities.get(exon_num, 0)
        if fwd_sim > rev_sim and fwd_sim > 0.:
            exon_comparison[exon_num] = {
                "winner": "forward", "similarity": fwd_sim}
        elif rev_sim > fwd_sim and rev_sim > 0.9:
            exon_comparison[exon_num] = {
                "winner": "reverse", "similarity": rev_sim}
        else:
            exon_comparison[exon_num] = {
                "winner": "tie", "similarity": fwd_sim}

    count_forward_wins = sum(
        1 for v in exon_comparison.values() if v['winner'] == 'forward')
    count_reverse_wins = sum(
        1 for v in exon_comparison.values() if v['winner'] == 'reverse')
    count_ties = sum(1 for v in exon_comparison.values()
                     if v['winner'] == 'tie')
    summary = {
        "forward_wins": count_forward_wins,
        "reverse_wins": count_reverse_wins,
        "ties": count_ties
    }
    if summary["forward_wins"] > summary["reverse_wins"]:

        best_match = "forward"

    elif summary["reverse_wins"] > summary["forward_wins"]:
        best_match = "reverse"
    else:
        best_match = "forward"

    selected_strand = strand[best_match]
    aboRef = service.getABO_ref("exons")

    exons_ref = []
    exon_combination = []
    for i in selected_strand['exon_alignments']:
        if i['similarity'] > 0.9:
            x = i['exon_number']-1
            exon = {}
            exon['exon'] = i['exon_number']

            exon['ref_start'] = i['ref_start']
            exon['ref_end'] = i['ref_end']

            exon['cds_start'] = aboRef[x]['cds_start']
            exon['cds_end'] = aboRef[x]['cds_end']
            exons_ref.append(exon)
            exon_combination.append(i['exon_number'])

    filtered_exons = [exon for exon in selected_strand['exon_alignments']
                      if exon['exon_number'] in [e['exon'] for e in exons_ref]]
    selected_strand['exon_alignments'] = filtered_exons
    selected_strand['exon_combination'] = exon_combination

    total_variants = 0
    SNPs = 0
    insertions = 0
    deletions = 0
    for exon in selected_strand['exon_alignments']:
        total_variants += len(exon['variants'])
        for var in exon['variants']:
            if var['type'] == 'SNP':
                SNPs += 1
            elif var['type'] == 'insertion':
                insertions += 1
            elif var['type'] == 'deletion':
                deletions += 1
    selected_strand['total_variants'] = total_variants
    selected_strand['variant_summary'] = {
        "SNPs": SNPs,
        "insertions": insertions,
        "deletions": deletions
    }

    return selected_strand, exons_ref
