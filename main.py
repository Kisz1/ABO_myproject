import streamlit as st
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import io

import utils.FASTA_analyzer as fasta_utils
import utils.ab1_analyzer as ab1_utils
import utils.abo_identifier as abo_utils
from utils.rhd_analyzer import RHDAnalyzer
from utils.rhce_analyzer import RHCEAnalyzer, RHCEReferenceMissingError
from utils.kel_analyzer import KELAnalyzer, KELReferenceMissingError
from utils.fy_analyzer import FYAnalyzer, FYReferenceMissingError
from utils.jk_analyzer import JKAnalyzer, JKReferenceMissingError
from utils.h_analyzer import HAnalyzer, HReferenceMissingError
from utils.mns_analyzer import MNSAnalyzer, MNSReferenceMissingError
from utils.di_analyzer import DIAnalyzer, DIReferenceMissingError
from utils.isbt_handler import ISBTDataHandler
from utils.bloodgroup import get_system, get_available_system_keys
from utils.bloodgroup.router import route_filename

import plotly.graph_objects as go
import itertools

st.set_page_config(layout="wide")

RHD_REFERENCE_OPTIONS = [
    "Built-in positive reference",
    "Built-in negative reference",
    "Upload custom reference"
]


IUPAC_CODES = {
    'A':	'A',
    'C':	'C',
    'G':	'G',
    'T':	'T',
    'R':	'A or G',
    'Y':	'C or T',
    'S':	'G or C',
    'W':	'A or T',
    'K':	'G or T',
    'M':	'A or C',
    'B':	'C or G or T',
    'D':	'A or G or T',
    'H':	'A or C or T',
    'V':	'A or C or G',
    'N':	'A or C or G or T'}


def plot_chromatogram_plotly_old(trace, base_width=2):
    """
    Create an interactive chromatogram plot with Plotly.

    Parameters
    ----------
    trace : dict
        Dictionary containing "A", "C", "G", "T", "pos", and "seq".
    base_width : float
        Controls horizontal scaling (pixels per base).
    """
    seq_len = len(trace["seq"])
    width_px = max(800, int(seq_len * base_width))

    fig = go.Figure()
    x = list(range(seq_len))

    # Plot each base channel
    fig.add_trace(go.Scatter(
        y=trace["A"], x=x, mode="lines", name="A", line=dict(color="green", width=1)))
    fig.add_trace(go.Scatter(
        y=trace["C"], x=x, mode="lines", name="C", line=dict(color="blue", width=1)))
    fig.add_trace(go.Scatter(
        y=trace["G"], x=x, mode="lines", name="G", line=dict(color="black", width=1)))
    fig.add_trace(go.Scatter(
        y=trace["T"], x=x, mode="lines", name="T", line=dict(color="red", width=1)))

    # Optional: base labels (every Nth base)
    step = max(1, seq_len // 100)
    labels = [trace["seq"][i] if i % step == 0 else "" for i in range(seq_len)]
    fig.add_trace(go.Scatter(
        x=x,
        y=[0]*seq_len,
        text=labels,
        mode="text",
        textposition="top center",
        hoverinfo="skip",
        showlegend=False
    ))

    fig.update_layout(
        title=f"ABO Exon: {trace['exon']} (length={seq_len})",
        width=width_px,
        height=400,
        xaxis_title="Base Index",
        yaxis_title="Signal Intensity",
        hovermode="x unified",
        template="plotly_white",
        legend=dict(orientation="h", y=-0.2),
        margin=dict(l=40, r=20, t=40, b=40)
    )

    return fig


def plot_chromatogram_plotly(trace, base_width=2, hetero_sites=None,
                             cds_start=None, cds_end=None):
    """
    Interactive chromatogram plot with exon coordinate support and heterozygous tooltips.

    Parameters
    ----------
    trace : dict
        Must contain "A", "C", "G", "T", "pos", and "seq".
    base_width : float
        Horizontal scaling (pixels per base).
    hetero_sites : list[tuple[int, dict]]
        List of (base_index, base_intensity_dict) for heterozygous peaks.
    cds_start, cds_end : int, optional
        Genomic coordinates for the exon region (for x-axis labels).
    """
    seq_len = len(trace["seq"])
    width_px = max(800, int(seq_len * base_width))

    # Define coordinate offset
    start_offset = cds_start if cds_start is not None else 0
    x = [start_offset + i for i in range(seq_len)]

    fig = go.Figure()

    # Base channels
    fig.add_trace(go.Scatter(
        x=x, y=trace["A"], mode="lines", name="A", line=dict(color="green", width=1)))
    fig.add_trace(go.Scatter(
        x=x, y=trace["C"], mode="lines", name="C", line=dict(color="blue", width=1)))
    fig.add_trace(go.Scatter(
        x=x, y=trace["G"], mode="lines", name="G", line=dict(color="black", width=1)))
    fig.add_trace(go.Scatter(
        x=x, y=trace["T"], mode="lines", name="T", line=dict(color="red", width=1)))

    # Optional base labels
    step = max(1, seq_len // 100)
    labels = [trace["seq"][i] if i % step == 0 else "" for i in range(seq_len)]
    fig.add_trace(go.Scatter(
        x=x,
        y=[0]*seq_len,
        text=labels,
        mode="text",
        textposition="top center",
        hoverinfo="skip",
        showlegend=False
    ))

    # Heterozygous markers with tooltip
    if hetero_sites:
        ymax = max(max(trace["A"]), max(trace["C"]),
                   max(trace["G"]), max(trace["T"]))
        for base_idx, intensities in hetero_sites:
            pos = start_offset + base_idx
            bases_sorted = sorted(intensities.items(),
                                  key=lambda kv: kv[1], reverse=True)
            top = bases_sorted[0]
            second = bases_sorted[1] if len(bases_sorted) > 1 else ("", 0)
            ratio = f"{second[0]}:{second[1]} / {top[0]}:{top[1]}  (ratio={second[1]/(top[1]+1e-6):.2f})"
            hover_text = "<br>".join(
                [f"{b}: {v}" for b, v in intensities.items()]) + f"<br>{ratio}"
            fig.add_trace(go.Scatter(
                x=[pos, pos],
                y=[0, ymax],
                mode="lines",
                line=dict(color="orange", width=1, dash="dot"),
                name=f"Hetero {pos}",
                hovertemplate=f"Pos {pos}<br>{hover_text}<extra></extra>"
            ))

    exon_label = trace.get('exon')
    title_label = f"ABO Exon {exon_label}" if exon_label is not None else "Raw Trace"
    fig.update_layout(
        title=f"Chromatogram – {title_label} (len={seq_len})",
        width=width_px,
        height=400,
        xaxis_title="Genomic Position" if cds_start else "Base Index",
        yaxis_title="Signal Intensity",
        hovermode="x unified",
        template="plotly_white",
        legend=dict(orientation="h", y=-0.2),
        margin=dict(l=40, r=20, t=40, b=40)
    )

    return fig


def display_alignment_with_snps(aligned_query, aligned_reference, cds_start=None, cds_end=None, variants=None, exon_number=None, unique_id=None):
    """
    Display aligned sequences in code boxes with copy functionality and red highlighting for variants.
    """

    if exon_number:
        st.write(f"#### 🧬 Exon {exon_number} Alignment")
    else:
        st.write("#### 🧬 Sequence Alignment")

    # Build sequences with HTML highlighting for differences
    ref_html = "REF: "
    query_html = "QRY: "

    for ref_base, query_base in zip(aligned_reference, aligned_query):
        if ref_base != query_base:
            # Red highlighting for differences
            ref_html += f'<span style="color: red; font-weight: bold; background-color: #ffeeee;">{ref_base}</span>'
            query_html += f'<span style="color: red; font-weight: bold; background-color: #ffeeee;">{query_base}</span>'
        else:
            ref_html += ref_base
            query_html += query_base

    # Display sequences with highlighting
    st.write("**Reference Sequence:**")
    st.markdown(f'<div style="background-color: #f8f9fa; padding: 10px; border-radius: 5px; border: 1px solid #e9ecef; font-family: monospace; font-size: 14px; position: relative;">{ref_html}</div>',
                unsafe_allow_html=True)

    # Also provide plain text version for easy copying
    ref_plain = f"REF: {aligned_reference}"
    with st.container():
        key_suffix = f"{unique_id}_{exon_number}" if unique_id else f"{exon_number}_{str(aligned_reference)[:10]}"
        if st.checkbox("Show plain text reference", key=f"align_ref_{key_suffix}"):
            st.code(ref_plain, language=None)

    st.write("**Query Sequence:**")
    st.markdown(f'<div style="background-color: #f8f9fa; padding: 10px; border-radius: 5px; border: 1px solid #e9ecef; font-family: monospace; font-size: 14px; position: relative;">{query_html}</div>',
                unsafe_allow_html=True)

    # Also provide plain text version for easy copying
    query_plain = f"QRY: {aligned_query}"
    with st.container():
        key_suffix = f"{unique_id}_{exon_number}" if unique_id else f"{exon_number}_{str(aligned_query)[:10]}"
        if st.checkbox("Show plain text query", key=f"align_query_{key_suffix}"):
            st.code(query_plain, language=None)

    # Analyze differences and show summary
    differences = []
    snps = 0
    insertions = 0
    deletions = 0

    for i, (ref_base, query_base) in enumerate(zip(aligned_reference, aligned_query)):
        if ref_base != query_base:
            pos = i + 1
            if ref_base == '-':
                insertions += 1
                differences.append(f"Pos {pos}: Insertion ({query_base})")
            elif query_base == '-':
                deletions += 1
                differences.append(f"Pos {pos}: Deletion ({ref_base})")
            else:
                snps += 1
                differences.append(f"Pos {pos}: SNP ({ref_base}→{query_base})")

    # Display summary statistics
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Total Differences", len(differences))
    with col2:
        st.metric("SNPs", snps)
    with col3:
        st.metric("Insertions", insertions)
    with col4:
        st.metric("Deletions", deletions)


def display_detailed_alignment_table(aligned_query, aligned_reference, variants=None, cds_start=None, cds_end=None):
    """
    Display detailed alignment in table format with position-by-position comparison.
    """

    st.write("#### 📊 Detailed Position-by-Position Analysis")

    # Prepare data for table
    positions = []
    ref_bases = []
    query_bases = []
    match_status = []
    in_cds = []
    variant_info = []

    # Create variant lookup by position
    variant_lookup = {}
    if variants:
        for var in variants:
            pos = var.get('alignment_pos', var.get(
                'ref_pos', var.get('position')))
            if pos is not None:
                variant_lookup[pos] = var

    for i, (ref_base, query_base) in enumerate(zip(aligned_reference, aligned_query)):
        pos = i + 1
        positions.append(pos)
        ref_bases.append(ref_base)
        query_bases.append(query_base)

        # Determine match status
        if ref_base == query_base:
            match_status.append("✅ Match")
        elif ref_base == '-':
            match_status.append("🔹 Insertion")
        elif query_base == '-':
            match_status.append("🔸 Deletion")
        else:
            match_status.append("❌ SNP")

        # Check if in CDS
        if cds_start and cds_end and cds_start <= pos <= cds_end:
            in_cds.append("Yes")
        else:
            in_cds.append("No")

        # Add variant information
        if pos in variant_lookup:
            var = variant_lookup[pos]
            variant_info.append(
                f"ISBT: {var.get('isbt_pos', 'N/A')}, Type: {var.get('type', 'N/A')}")
        else:
            variant_info.append("-")

    # Show only positions with differences (SNPs, insertions, deletions)
    show_all = st.checkbox("Show all positions", value=False)

    if not show_all:
        # Filter to show only differences
        filtered_data = []
        for i in range(len(positions)):
            if match_status[i] != "✅ Match":
                filtered_data.append({
                    "Position": positions[i],
                    "Reference": ref_bases[i],
                    "Query": query_bases[i],
                    "Status": match_status[i],
                    "In CDS": in_cds[i],
                    "Variant Info": variant_info[i]
                })

        if filtered_data:
            df = pd.DataFrame(filtered_data)
            st.write(
                f"**Showing {len(filtered_data)} positions with differences:**")
            st.dataframe(df, width='stretch', height=400)
        else:
            st.success("No differences found - sequences are identical!")
    else:
        # Show all positions
        all_data = {
            "Position": positions,
            "Reference": ref_bases,
            "Query": query_bases,
            "Status": match_status,
            "In CDS": in_cds,
            "Variant Info": variant_info
        }
        df = pd.DataFrame(all_data)
        st.write(f"**Showing all {len(positions)} positions:**")
        st.dataframe(df, width='stretch', height=400)


def process_ab1_files(fwd_ab1_files, exons_ref, threshold_ratio=0.3):
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
            continue

    if results:
        return results, all_hets if all_hets else None

    return None, None


_BASES_TO_IUPAC = {
    frozenset('AG'): 'R', frozenset('CT'): 'Y',
    frozenset('CG'): 'S', frozenset('AT'): 'W',
    frozenset('GT'): 'K', frozenset('AC'): 'M',
}


def quality_trim_and_mask(seq, qual, q_threshold=20, window=10):
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


def process_rhd_ab1_files(rhd_ab1_files, q_threshold=20, window=10, het_ratio=0.30):
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

        except Exception:
            continue

    if results:
        return results, (all_hets if all_hets else None)
    return None, None


def process_rhd_fasta_files(rhd_fasta_files):
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
        except Exception:
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


def analyze_rhd_multifactor(ab1_traces):
    """
    Analyze RHD using RHDAnalyzer with embedded WHO references.
    Auto-detects amplicon region and applies correct decision logic.
    """
    analyzer = RHDAnalyzer()
    amplicon_results = []

    for i, trace in enumerate(ab1_traces):
        if not isinstance(trace, dict):
            continue

        query_seq = trace.get('seq', '')
        if not query_seq or len(query_seq) < 50:
            continue

        r = analyzer.analyze(query_seq)

        status_map = {
            "RhD+":         "RhD+ (D positive)",
            "RhD-":         "RhD- (D negative)",
            "RHD Variant":  "RHD Variant - confirm required",
            "Inconclusive": "Inconclusive - confirm required",
        }
        phenotype = status_map.get(r['rhd_status'], r['rhd_status'])
        region    = r['region']
        identity  = r['identity'] or 0.0
        variants  = r['variants']

        result = {
            'amplicon_index':           i,
            'query_length':             r['query_length'],
            'reference_length':         951 if region == "RHD1" else 3336,
            'matched_region':           region,
            'reference_description':    f"{region}_ref (WHO standard)",
            'matched_identity':         round(identity, 1),
            'identity':                 round(identity, 1),
            'matched_score':            0.0,
            'variant_count':            len(variants),
            'variants':                 variants,
            'score_1':                  round(identity, 1) if region == "RHD1"  else 0.0,
            'score_456':                round(identity, 1) if region == "RHD456" else 0.0,
            'phenotype':                phenotype,
            'reason':                   r['reason'],
            'rule':                     f"{region} rule (strand: {r['strand']})",
            'strand':                   r['strand'],
            'region_assignment_reason': f"{region} auto-detected by alignment",
        }
        amplicon_results.append(result)

    return amplicon_results


def consolidate_rhd_results(amplicon_results):
    """
    Consolidate per-amplicon RHD results into a final verdict.
    If multiple amplicons: check agreement. If single: use that result.
    """
    if not amplicon_results:
        return None
    
    if len(amplicon_results) == 1:
        return amplicon_results[0]
    
    phenotypes = [r.get('phenotype', 'Unknown') for r in amplicon_results]
    positive_count = sum(1 for p in phenotypes if 'RhD+' in p)
    negative_count = sum(1 for p in phenotypes if 'RhD-' in p)
    
    if positive_count == len(amplicon_results):
        result = amplicon_results[0].copy()
        result['phenotype'] = 'RhD+ (confirmed)'
        result['reason'] = 'All amplicons indicate RHD gene present'
        result['multi_amplicon'] = True
    elif negative_count == len(amplicon_results):
        result = amplicon_results[0].copy()
        result['phenotype'] = 'RhD- (confirmed)'
        result['reason'] = 'All amplicons indicate RHD gene absent'
        result['multi_amplicon'] = True
    else:
        result = {
            'phenotype': 'Inconclusive',
            'reason': 'Results inconsistent across amplicons',
            'multi_amplicon': True,
            'amplicon_count': len(amplicon_results)
        }
    
    return result

    types_list = {'SNP': 'alt_base', 'insertion': 'inserted_sequence',
                  'deletion': 'deleted_sequence'}
    het_variants = []
    var_nodes = []
    unknown = []
    variant_base = i.get(types_list[types], "")
    possible_bases = IUPAC_CODES.get(variant_base, variant_base).split(" or ")
    ref_base = i.get('ref_base')

    for base in possible_bases:
        if not base:
            continue

        field = types_list[i['type']]
        var_node = None
        if i['type'] == 'deletion':
            var_node = abo_identifier.get_variant_node(
                i['isbt_pos'], base, "", variant_base)
        elif i['type'] == 'insertion':
            var_node = abo_identifier.get_variant_node(
                i['isbt_pos'], "", base, variant_base)
        else:
            var_node = abo_identifier.get_variant_node(
                i['isbt_pos'], ref_base, base, variant_base)

        if var_node is not None:
            var_nodes.append(var_node)
            het_variants.append(i[field])
        else:
            if ref_base is None or base != ref_base:
                unknown.append(i)
    return var_nodes, het_variants, unknown


def get_display_base(base):
    """Convert IUPAC code to display string."""
    if base in IUPAC_CODES:
        base_display = IUPAC_CODES[base]
        if base == base_display:
            return base
        else:
            return f"{base} ({base_display})"
    return base

def handle_IUPAC_codes(abo_identifier, i, types):
    types_list = {
        'SNP':       'alt_base',
        'insertion': 'inserted_sequence',
        'deletion':  'deleted_sequence'
    }
    het_variants = []
    var_nodes    = []
    unknown      = []

    variant_base   = i.get(types_list[types], "")
    possible_bases = IUPAC_CODES.get(variant_base, variant_base).split(" or ")
    ref_base       = i.get('ref_base')

    for base in possible_bases:
        if not base:
            continue

        var_node = None
        if i['type'] == 'deletion':
            var_node = abo_identifier.get_variant_node(
                i['isbt_pos'], base, "", variant_base)
        elif i['type'] == 'insertion':
            var_node = abo_identifier.get_variant_node(
                i['isbt_pos'], "", base, variant_base)
        else:
            var_node = abo_identifier.get_variant_node(
                i['isbt_pos'], ref_base, base, variant_base)

        field = types_list[i['type']]
        if var_node is not None:
            var_nodes.append(var_node)
            het_variants.append(i[field])
        else:
            if ref_base is None or base != ref_base:
                unknown.append(i)

    return var_nodes, het_variants, unknown


def identify_abo_alleles(FASTA_variant_list, abo_identifier=None):
    if abo_identifier is None:
        abo_identifier = abo_utils.ABOIdentifier("ABO")
    var_nodes = []
    het_variants = []
    unknown = []
    for exon in FASTA_variant_list['exon_alignments']:
        variants = exon['variants']
        for i in variants:
            if i['type'] == 'insertion':
                var_node, het_var, unk = handle_IUPAC_codes(
                    abo_identifier, i, 'insertion')
            elif i['type'] == 'deletion':
                var_node, het_var, unk = handle_IUPAC_codes(
                    abo_identifier, i, 'deletion')
            else:
                var_node, het_var, unk = handle_IUPAC_codes(
                    abo_identifier, i, 'SNP')
            var_nodes.extend(var_node)
            het_variants.extend(het_var)
            unknown.extend(unk)
    alleles = []
    node_iupac_map = {node[0]: node[2] for node in var_nodes}
    print(node_iupac_map)
    for node_name, node_data, iupac_code in var_nodes:
        # The identify_alleles method expects a list of tuples (node_name, node_data)
        # So we pass the current node as a list containing one tuple
        allele = abo_identifier.identify_alleles([(node_name, node_data)])
        # Use extend to add elements from the list
        alleles.append({node_name: allele})

    # Extract the lists of alleles from the 'alleles' list of dictionaries
    allele_lists = [list(d.values())[0] for d in alleles]

    # Find the intersection of all allele lists
    if allele_lists:
        # Start with the first list
        common_alleles = set(allele_lists[0])
        # Intersect with the remaining lists
        for allele_list in allele_lists[1:]:
            common_alleles.intersection_update(allele_list)
    else:
        common_alleles = set()
    allele_variants_list = []
    for i in common_alleles:
        v = abo_identifier.get_variants_for_allele(i)
        av_list = []
        for j in v:
            gene, location, change = j[0].split("_")
            exon = abo_identifier.get_exon(location)
            av_list.append({"name": j[0], "exon": exon, "location": int(
                location), "change": change})  # Convert location to int for sorting

        av_list.sort(key=lambda x: x['location'])

        allele_variants_list.append({i: av_list})
    allele_variants_list.sort(key=lambda x: list(x.keys())[
                              0])  # Sort by allele name

    variants_name = [x[0] for x in var_nodes]

    unknown_alleles_to_display = []
    for u in unknown:
        item = {}
        exon = abo_identifier.get_exon(u.get('isbt_pos'))
        item['isbt_pos'] = u.get('isbt_pos')
        if exon is None:
            item['exon'] = 'N/A'
        item['type'] = u.get('type')
        if u.get('type') == 'deletion':
            item['ref_base'] = get_display_base(
                u.get('deleted_sequence', 'N/A')),
            item['alt_base'] = '-',
        elif u.get('type') == 'insertion':
            item['ref_base'] = '-',
            item['alt_base'] = get_display_base(
                u.get('inserted_sequence', 'N/A')),

        else:
            item['ref_base'] = get_display_base(u.get('ref_base', 'N/A')),
            item['alt_base'] = get_display_base(u.get('alt_base', 'N/A')),

        unknown_alleles_to_display.append(item)
    return allele_variants_list, unknown_alleles_to_display, variants_name, node_iupac_map


def get_display_iupac_change(change, iupac_code):
    change_list = change.split(">")
    if len(change_list) > 1:
        ref_base = change_list[0]
        return f"{ref_base}>{iupac_code}"
    else:
        change = change[:-1] + f"({iupac_code})"
        return change


st.title("🧬 ABO blood group analysis")

# Upload section
# with st.sidebar:

# ─── Unified auto-detect upload ──────────────────────────────────────────
# One uploader that routes each file to its blood-group system by filename
# pattern (RHCE, RHD, KEL56, FY12, JK78, MIA234, DI1819, etc.). Routed
# reads are appended to the existing per-system input lists below, so the
# downstream pipeline is unchanged. Manual per-system uploaders remain as
# an override path.
st.sidebar.markdown("### 🚀 Auto-detect Upload (all systems)")
unified_ab1_files = st.sidebar.file_uploader(
    "Drop AB1 files of any blood-group amplicon",
    type=["ab1"], accept_multiple_files=True,
    key="unified_ab1",
    help=("Filenames are matched against amplicon-naming conventions "
          "(RHCE45, RHD1, KEL56, FY12, JK78, MIA234, DI1819, ...). "
          "Unrecognized names are flagged 'unrouted' and skipped."))
unified_fasta_files = st.sidebar.file_uploader(
    "Drop FASTA files of any blood-group amplicon",
    type=["fasta", "fa", "fas"], accept_multiple_files=True,
    key="unified_fasta",
    help=("Same routing rules as the AB1 uploader above. Files for ABO/RHD/"
          "RHCE/Kell/Duffy/Kidd/H/MNS/Diego are auto-detected and dispatched."))

st.sidebar.markdown("---")
st.sidebar.markdown("### 🅰️🅱️ ABO Inputs")
fwd_ab1 = st.sidebar.file_uploader(
    "Upload ABO AB1 file",
    type=["ab1"], accept_multiple_files=True,
    help="AB1 chromatogram for ABO analysis (chromatogram tab + heterozygote detection).")

fasta_files = st.sidebar.file_uploader(
    "Upload ABO FASTA file",
    type=["fasta", "fa", "fas"], accept_multiple_files=True,
    help="Exon-specific FASTA aligned to ABO references.")

st.sidebar.markdown("### 🩸 RHD Inputs")
rhd_ab1_files = st.sidebar.file_uploader(
    "Upload RHD AB1 file",
    type=["ab1"], accept_multiple_files=True,
    help="AB1 chromatogram for RHD multi-amplicon analysis.")

rhd_fasta_files = st.sidebar.file_uploader(
    "Upload RHD FASTA file",
    type=["fasta", "fa", "fas"], accept_multiple_files=True,
    help="RHD FASTA (e.g. RHD1, RHD456). Routed to the RHD analyzer with multi-amplicon voting.")

st.sidebar.markdown("### 🩸 RHCE Inputs")
rhce_ab1_files = st.sidebar.file_uploader(
    "Upload RHCE AB1 file(s)",
    type=["ab1"], accept_multiple_files=True,
    help="AB1 chromatograms covering RHCE exons 1, 2, 5 (C/c + E/e diagnostic regions).")

rhce_fasta_files = st.sidebar.file_uploader(
    "Upload RHCE FASTA file(s)",
    type=["fasta", "fa", "fas"], accept_multiple_files=True,
    help="Multiple FASTA reads recommended for international-standard multi-read consensus voting.")

st.sidebar.markdown("### 🩸 Kell Inputs")
kel_ab1_files = st.sidebar.file_uploader(
    "Upload Kell AB1 file(s)",
    type=["ab1"], accept_multiple_files=True,
    help="AB1 chromatograms covering KEL exon 6 (K/k diagnostic region, c.578C>T).")

kel_fasta_files = st.sidebar.file_uploader(
    "Upload Kell FASTA file(s)",
    type=["fasta", "fa", "fas"], accept_multiple_files=True,
    help="FASTA reads covering KEL exon 6. Multiple reads enable consensus voting.")

st.sidebar.markdown("### 🩸 Duffy Inputs")
fy_ab1_files = st.sidebar.file_uploader(
    "Upload Duffy AB1 file(s)",
    type=["ab1"], accept_multiple_files=True,
    help=("AB1 chromatograms covering Duffy exon 2 (c.125G>A FYA/FYB, c.265C>T "
          "FY*X) and/or the promoter (-67T>C GATA, requires genomic amplicon)."))

fy_fasta_files = st.sidebar.file_uploader(
    "Upload Duffy FASTA file(s)",
    type=["fasta", "fa", "fas"], accept_multiple_files=True,
    help=("FASTA reads covering the Duffy diagnostic SNPs. Multiple reads enable "
          "consensus voting. Promoter SNP (-67T>C) needs a genomic amplicon."))

st.sidebar.markdown("### 🩸 Kidd Inputs")
jk_ab1_files = st.sidebar.file_uploader(
    "Upload Kidd AB1 file(s)",
    type=["ab1"], accept_multiple_files=True,
    help=("AB1 chromatograms covering Kidd exon 9 (c.838G>A JKA/JKB, c.871T>C "
          "Asian-null) and/or intron 5 (c.342-1G>A Polynesian-null splice, "
          "requires genomic amplicon)."))

jk_fasta_files = st.sidebar.file_uploader(
    "Upload Kidd FASTA file(s)",
    type=["fasta", "fa", "fas"], accept_multiple_files=True,
    help=("FASTA reads covering the Kidd diagnostic SNPs. Multiple reads enable "
          "consensus voting. Splice SNP (c.342-1) needs a genomic amplicon."))

st.sidebar.markdown("### 🩸 H Inputs")
h_ab1_files = st.sidebar.file_uploader(
    "Upload H AB1 file(s)",
    type=["ab1"], accept_multiple_files=True,
    help=("AB1 chromatograms covering FUT1 exon 4 (c.460T>C h2, c.586C>T "
          "nonsense, c.725T>G Indian Bombay). Detects Bombay (Oh) and "
          "Para-Bombay carriers."))

h_fasta_files = st.sidebar.file_uploader(
    "Upload H FASTA file(s)",
    type=["fasta", "fa", "fas"], accept_multiple_files=True,
    help=("FASTA reads covering the FUT1 diagnostic SNPs. Multiple reads "
          "enable consensus voting."))

st.sidebar.markdown("### 🩸 MNS Inputs")
mns_ab1_files = st.sidebar.file_uploader(
    "Upload MNS AB1 file(s)",
    type=["ab1"], accept_multiple_files=True,
    help=("AB1 chromatograms covering GYPA exon 2 (c.59C>T, M/N) and/or "
          "GYPB exon 3 (c.143C>T, S/s). Drop reads from either gene — the "
          "analyzer dispatches per-SNP."))

mns_fasta_files = st.sidebar.file_uploader(
    "Upload MNS FASTA file(s)",
    type=["fasta", "fa", "fas"], accept_multiple_files=True,
    help=("FASTA reads for GYPA and/or GYPB. Multiple reads enable "
          "consensus voting; both genes can be uploaded together."))

st.sidebar.markdown("### 🩸 Diego Inputs")
di_ab1_files = st.sidebar.file_uploader(
    "Upload Diego AB1 file(s)",
    type=["ab1"], accept_multiple_files=True,
    help=("AB1 chromatograms covering SLC4A1 exon 19 (c.2561T>C, "
          "Di(a)/Di(b) discriminator). Lab amplicon: DI1819."))

di_fasta_files = st.sidebar.file_uploader(
    "Upload Diego FASTA file(s)",
    type=["fasta", "fa", "fas"], accept_multiple_files=True,
    help=("FASTA reads covering SLC4A1 c.2561. Multiple reads enable "
          "consensus voting."))

exon_start = st.sidebar.number_input(
    "Exon start (optional)", min_value=0, value=0)
exon_end = st.sidebar.number_input(
    "Exon end (optional, 0 = full length)", min_value=0, value=0)

threshold_ratio = st.sidebar.slider(
    "Heterozygosity threshold ratio", 0.1, 0.9, 0.3, 0.05)

analyze_button = st.sidebar.button("Analyze")

st.sidebar.markdown("---")
st.sidebar.caption("Version 2.0")


exons_ref = []


def get_cds(exon_number):
    for exon in exons_ref:
        if exon['exon'] == exon_number:
            return exon['cds_start'], exon['cds_end']
    return None, None


def generate_final_blood_group_summary(robust_summary, processed_AB1, hets, isbt_handler, has_fasta, rhd_reference_seq=None, rhce_result=None, kel_result=None, fy_result=None, jk_result=None, h_result=None, mns_result=None, di_result=None):
    """
    Generate a comprehensive final blood group summary from all analyzed systems.

    Args:
        robust_summary: ABO analysis results from FASTA files
        processed_AB1: AB1 processing results (for RHD)
        hets: Heterozygote data
        isbt_handler: ISBT handler for phenotype suggestions
        has_fasta: bool indicating whether FASTA files were uploaded
        rhce_result: optional dict from RHCEAnalyzer.analyze() with consensus result
        kel_result: optional dict from KELAnalyzer.analyze() with consensus result
        fy_result: optional dict from FYAnalyzer.analyze() with consensus result
        jk_result: optional dict from JKAnalyzer.analyze() with consensus result
        h_result:  optional dict from HAnalyzer.analyze()  with consensus result

    Returns:
        Dictionary with final blood group summary
    """
    summary = {
        'abo': {'status': 'not_analyzed', 'phenotype': None, 'alleles': [], 'variants': [], 'message': None},
        'rhd': {'status': 'not_analyzed', 'phenotype': None, 'alleles': [], 'variants': [], 'identity': None, 'query_length': None, 'reference_length': None, 'note': None},
        'rhce': {'status': 'not_analyzed', 'phenotype': None, 'c_e': None, 'big_E': None, 'allele_options': [], 'confidence': None, 'reads_callable': 0, 'reads_total': 0, 'partial_markers': [], 'message': None},
        'kel': {'status': 'not_analyzed', 'phenotype': None, 'k_axis': None, 'allele_options': [], 'confidence': None, 'reads_callable': 0, 'reads_total': 0, 'message': None},
        'fy':  {'status': 'not_analyzed', 'phenotype': None, 'a_b': None, 'gata': None, 'fy_x': None, 'allele_options': [], 'confidence': None, 'reads_callable': 0, 'reads_total': 0, 'phase_ambiguous': False, 'message': None},
        'jk':  {'status': 'not_analyzed', 'phenotype': None, 'a_b': None, 'null_state': None, 'null_triggers': [], 'allele_options': [], 'confidence': None, 'reads_callable': 0, 'reads_total': 0, 'phase_ambiguous': False, 'message': None},
        'h':   {'status': 'not_analyzed', 'phenotype': None, 'state': None, 'triggers_hom_alt': [], 'triggers_het': [], 'allele_options': [], 'confidence': None, 'reads_callable': 0, 'reads_total': 0, 'bombay_warning': False, 'message': None},
        'mns': {'status': 'not_analyzed', 'phenotype': None, 'm_n': None, 's_s': None, 'allele_options': [], 'confidence': None, 'reads_callable': 0, 'reads_total': 0, 'phase_ambiguous': False, 'message': None},
        'di':  {'status': 'not_analyzed', 'phenotype': None, 'genotype': None, 'allele_options': [], 'confidence': None, 'reads_callable': 0, 'reads_total': 0, 'message': None},
    }

    # ABO Analysis Summary
    if not has_fasta:
        summary['abo']['status'] = 'missing_data'
        summary['abo']['message'] = 'ABO analysis requires FASTA file upload.'
    elif robust_summary and isinstance(robust_summary, list) and len(robust_summary) > 0:
        confirmed_results = [r for r in robust_summary if "Confirmed" in r.get('decision', "")]

        if confirmed_results and len(confirmed_results) > 0:
            # Extract all variants from confirmed exons
            all_variants = []
            for result in confirmed_results:
                if result.get('variants'):
                    for variant in result['variants']:
                        if variant.get('isbt_pos'):
                            all_variants.append(f"c.{variant['isbt_pos']}{variant.get('type', '')}")

            if all_variants and len(all_variants) > 0:
                # Get ABO phenotype suggestions
                abo_suggestion = isbt_handler.suggest_blood_group_from_variants(all_variants, "ABO")
                if abo_suggestion['status'] == 'matched':
                    summary['abo']['phenotype'] = abo_suggestion['phenotypes'][0] if abo_suggestion['phenotypes'] else None
                    summary['abo']['alleles'] = [allele.get('name', '') for allele in abo_suggestion['matches']]
                summary['abo']['status'] = 'analyzed'
                summary['abo']['variants'] = all_variants
            else:
                summary['abo']['status'] = 'no_variants'
                summary['abo']['phenotype'] = 'Reference'  # No variants = reference phenotype
        else:
            # Results exist but none are confirmed
            summary['abo']['status'] = 'not_confirmed'
            summary['abo']['phenotype'] = 'Awaiting confirmation'
            summary['abo']['message'] = 'No confirmed exon matches were found in the FASTA results.'
    else:
        summary['abo']['status'] = 'failed'
        summary['abo']['message'] = 'No ABO FASTA results were produced.'

    # RHD Analysis Summary - Multi-Factor Decision Tree
    if processed_AB1 and isinstance(processed_AB1, list) and len(processed_AB1) > 0:
        try:
            amplicon_results = analyze_rhd_multifactor(processed_AB1)
            consolidated = consolidate_rhd_results(amplicon_results)
            
            if consolidated:
                summary['rhd']['status'] = 'analyzed'
                summary['rhd']['phenotype'] = consolidated.get('phenotype')
                summary['rhd']['reason'] = consolidated.get('reason')
                summary['rhd']['amplicon_results'] = amplicon_results
                summary['rhd']['variants'] = consolidated.get('variants', [])
                if len(amplicon_results) == 1:
                    r = amplicon_results[0]
                    summary['rhd']['identity'] = r.get('identity')
                    summary['rhd']['query_length'] = r.get('query_length')
                    summary['rhd']['reference_length'] = r.get('reference_length')
                    summary['rhd']['rule_applied'] = r.get('rule')
            else:
                summary['rhd']['status'] = 'no_sequence'
                summary['rhd']['phenotype'] = 'No valid sequences'
                
        except Exception as e:
            st.warning(f"⚠️ RHD analysis error: {str(e)}")
            summary['rhd']['status'] = 'error'
            summary['rhd']['phenotype'] = 'Analysis failed'

    # RHCE Analysis Summary - multi-read consensus
    if rhce_result is not None:
        if rhce_result.get('phenotype') == 'Indeterminate' or rhce_result.get('reads_callable', 0) == 0:
            summary['rhce']['status'] = 'indeterminate'
            summary['rhce']['phenotype'] = 'Indeterminate'
            summary['rhce']['message'] = rhce_result.get('reason')
        else:
            summary['rhce']['status'] = 'analyzed'
            summary['rhce']['phenotype'] = rhce_result.get('phenotype')
            summary['rhce']['c_e'] = (rhce_result.get('c_e_call') or {}).get('genotype')
            summary['rhce']['big_E'] = (rhce_result.get('big_E_call') or {}).get('genotype')
            summary['rhce']['allele_options'] = rhce_result.get('allele_options', [])
            summary['rhce']['confidence'] = rhce_result.get('overall_confidence')
            summary['rhce']['reads_callable'] = rhce_result.get('reads_callable', 0)
            summary['rhce']['reads_total'] = rhce_result.get('reads_total', 0)
            summary['rhce']['partial_markers'] = (rhce_result.get('big_E_call') or {}).get('partial_markers', [])

    # KEL Analysis Summary - multi-read consensus
    if kel_result is not None:
        if kel_result.get('phenotype') == 'Indeterminate' or kel_result.get('reads_with_trusted_snps', 0) == 0:
            summary['kel']['status'] = 'indeterminate'
            summary['kel']['phenotype'] = 'Indeterminate'
            summary['kel']['message'] = kel_result.get('reason')
        else:
            summary['kel']['status'] = 'analyzed'
            summary['kel']['phenotype'] = kel_result.get('phenotype')
            summary['kel']['k_axis'] = (kel_result.get('k_axis_call') or {}).get('genotype')
            summary['kel']['allele_options'] = kel_result.get('allele_options', [])
            summary['kel']['confidence'] = kel_result.get('overall_confidence')
            summary['kel']['reads_callable'] = kel_result.get('reads_callable', 0)
            summary['kel']['reads_total'] = kel_result.get('reads_total', 0)

    # FY (Duffy) Analysis Summary - multi-read consensus
    if fy_result is not None:
        if fy_result.get('phenotype') == 'Indeterminate' or fy_result.get('reads_with_trusted_snps', 0) == 0:
            summary['fy']['status'] = 'indeterminate'
            summary['fy']['phenotype'] = 'Indeterminate'
            summary['fy']['message'] = fy_result.get('reason')
        else:
            summary['fy']['status'] = 'analyzed'
            summary['fy']['phenotype'] = fy_result.get('phenotype')
            summary['fy']['a_b'] = (fy_result.get('a_b_call') or {}).get('genotype')
            summary['fy']['gata'] = (fy_result.get('gata_call') or {}).get('consensus')
            summary['fy']['fy_x'] = (fy_result.get('fy_x_call') or {}).get('consensus')
            summary['fy']['allele_options'] = fy_result.get('allele_options', [])
            summary['fy']['confidence'] = fy_result.get('overall_confidence')
            summary['fy']['reads_callable'] = fy_result.get('reads_callable', 0)
            summary['fy']['reads_total'] = fy_result.get('reads_total', 0)
            summary['fy']['phase_ambiguous'] = len(fy_result.get('allele_options', [])) > 1

    # JK (Kidd) Analysis Summary - multi-read consensus
    if jk_result is not None:
        if jk_result.get('phenotype') == 'Indeterminate' or jk_result.get('reads_with_trusted_snps', 0) == 0:
            summary['jk']['status'] = 'indeterminate'
            summary['jk']['phenotype'] = 'Indeterminate'
            summary['jk']['message'] = jk_result.get('reason')
        else:
            summary['jk']['status'] = 'analyzed'
            summary['jk']['phenotype'] = jk_result.get('phenotype')
            summary['jk']['a_b'] = (jk_result.get('a_b_call') or {}).get('genotype')
            summary['jk']['null_state'] = (jk_result.get('null_call') or {}).get('state')
            summary['jk']['null_triggers'] = (jk_result.get('null_call') or {}).get('triggers', [])
            summary['jk']['allele_options'] = jk_result.get('allele_options', [])
            summary['jk']['confidence'] = jk_result.get('overall_confidence')
            summary['jk']['reads_callable'] = jk_result.get('reads_callable', 0)
            summary['jk']['reads_total'] = jk_result.get('reads_total', 0)
            summary['jk']['phase_ambiguous'] = len(jk_result.get('allele_options', [])) > 1

    # H (FUT1 / Bombay) Analysis Summary - multi-read consensus
    if h_result is not None:
        if h_result.get('phenotype') == 'Indeterminate' or h_result.get('reads_with_trusted_snps', 0) == 0:
            summary['h']['status'] = 'indeterminate'
            summary['h']['phenotype'] = 'Indeterminate'
            summary['h']['message'] = h_result.get('reason')
        else:
            summary['h']['status'] = 'analyzed'
            summary['h']['phenotype'] = h_result.get('phenotype')
            summary['h']['state'] = (h_result.get('h_state_call') or {}).get('state')
            summary['h']['triggers_hom_alt'] = (h_result.get('h_state_call') or {}).get('triggers_hom_alt', [])
            summary['h']['triggers_het'] = (h_result.get('h_state_call') or {}).get('triggers_het', [])
            summary['h']['allele_options'] = h_result.get('allele_options', [])
            summary['h']['confidence'] = h_result.get('overall_confidence')
            summary['h']['reads_callable'] = h_result.get('reads_callable', 0)
            summary['h']['reads_total'] = h_result.get('reads_total', 0)
            # Bombay is a transfusion-safety flag — surface it explicitly.
            summary['h']['bombay_warning'] = summary['h']['state'] == 'bombay'

    # MNS (GYPA + GYPB) Analysis Summary - multi-read consensus
    if mns_result is not None:
        if mns_result.get('phenotype') == 'Indeterminate' or mns_result.get('reads_with_trusted_snps', 0) == 0:
            summary['mns']['status'] = 'indeterminate'
            summary['mns']['phenotype'] = 'Indeterminate'
            summary['mns']['message'] = mns_result.get('reason')
        else:
            summary['mns']['status'] = 'analyzed'
            summary['mns']['phenotype'] = mns_result.get('phenotype')
            summary['mns']['m_n'] = (mns_result.get('m_n_call') or {}).get('genotype')
            summary['mns']['s_s'] = (mns_result.get('s_s_call') or {}).get('genotype')
            summary['mns']['allele_options'] = mns_result.get('allele_options', [])
            summary['mns']['confidence'] = mns_result.get('overall_confidence')
            summary['mns']['reads_callable'] = mns_result.get('reads_callable', 0)
            summary['mns']['reads_total'] = mns_result.get('reads_total', 0)
            summary['mns']['phase_ambiguous'] = len(mns_result.get('allele_options', [])) > 1

    # DI (Diego / SLC4A1) Analysis Summary - multi-read consensus
    if di_result is not None:
        if di_result.get('phenotype') == 'Indeterminate' or di_result.get('reads_with_trusted_snps', 0) == 0:
            summary['di']['status'] = 'indeterminate'
            summary['di']['phenotype'] = 'Indeterminate'
            summary['di']['message'] = di_result.get('reason')
        else:
            summary['di']['status'] = 'analyzed'
            summary['di']['phenotype'] = di_result.get('phenotype')
            summary['di']['genotype'] = (di_result.get('di_axis_call') or {}).get('genotype')
            summary['di']['allele_options'] = di_result.get('allele_options', [])
            summary['di']['confidence'] = di_result.get('overall_confidence')
            summary['di']['reads_callable'] = di_result.get('reads_callable', 0)
            summary['di']['reads_total'] = di_result.get('reads_total', 0)

    return summary


def display_final_blood_group_result(summary):
    """
    Display the final blood group result prominently at the top of analysis results.

    Args:
        summary: Blood group summary from generate_final_blood_group_summary
    """
    st.markdown("---")
    st.markdown("## 🩸 **FINAL BLOOD GROUP RESULT**")

    # Create a prominent display
    result_parts = []

    # ABO result
    if summary['abo']['status'] == 'analyzed':
        phenotype = summary['abo']['phenotype'] or 'Variants detected'
        result_parts.append(f"**ABO: {phenotype}**")
    elif summary['abo']['status'] == 'no_variants':
        result_parts.append("**ABO: Reference**")
    else:
        result_parts.append("*ABO: Not analyzed*")

    # RHD result
    if summary['rhd']['status'] == 'analyzed' and summary['rhd']['phenotype']:
        result_parts.append(f"**RHD: {summary['rhd']['phenotype']}**")
    else:
        result_parts.append("*RHD: Not analyzed*")

    # RHCE result
    if summary['rhce']['status'] == 'analyzed' and summary['rhce']['phenotype']:
        conf = summary['rhce'].get('confidence') or ''
        conf_tag = f" ({conf})" if conf else ''
        result_parts.append(f"**RHCE: {summary['rhce']['phenotype']}{conf_tag}**")
    elif summary['rhce']['status'] == 'indeterminate':
        result_parts.append("**RHCE: Indeterminate**")
    else:
        result_parts.append("*RHCE: Not analyzed*")

    # KEL result
    if summary['kel']['status'] == 'analyzed' and summary['kel']['phenotype']:
        conf = summary['kel'].get('confidence') or ''
        conf_tag = f" ({conf})" if conf else ''
        result_parts.append(f"**Kell: {summary['kel']['phenotype']}{conf_tag}**")
    elif summary['kel']['status'] == 'indeterminate':
        result_parts.append("**Kell: Indeterminate**")
    else:
        result_parts.append("*Kell: Not analyzed*")

    # FY (Duffy) result
    if summary['fy']['status'] == 'analyzed' and summary['fy']['phenotype']:
        conf = summary['fy'].get('confidence') or ''
        conf_tag = f" ({conf})" if conf else ''
        amb = " ⚠" if summary['fy'].get('phase_ambiguous') else ''
        result_parts.append(f"**Duffy: {summary['fy']['phenotype']}{conf_tag}{amb}**")
    elif summary['fy']['status'] == 'indeterminate':
        result_parts.append("**Duffy: Indeterminate**")
    else:
        result_parts.append("*Duffy: Not analyzed*")

    # JK (Kidd) result
    if summary['jk']['status'] == 'analyzed' and summary['jk']['phenotype']:
        conf = summary['jk'].get('confidence') or ''
        conf_tag = f" ({conf})" if conf else ''
        amb = " ⚠" if summary['jk'].get('phase_ambiguous') else ''
        result_parts.append(f"**Kidd: {summary['jk']['phenotype']}{conf_tag}{amb}**")
    elif summary['jk']['status'] == 'indeterminate':
        result_parts.append("**Kidd: Indeterminate**")
    else:
        result_parts.append("*Kidd: Not analyzed*")

    # H (FUT1 / Bombay) result — Bombay gets a loud 🚨 because it's a
    # transfusion-safety issue (recipient rejects standard ABO blood).
    if summary['h']['status'] == 'analyzed' and summary['h']['phenotype']:
        conf = summary['h'].get('confidence') or ''
        conf_tag = f" ({conf})" if conf else ''
        bombay = " 🚨" if summary['h'].get('bombay_warning') else ''
        result_parts.append(f"**H: {summary['h']['phenotype']}{conf_tag}{bombay}**")
    elif summary['h']['status'] == 'indeterminate':
        result_parts.append("**H: Indeterminate**")
    else:
        result_parts.append("*H: Not analyzed*")

    # MNS result — phase-ambiguity (double-het MN/Ss) gets the ⚠ flag.
    if summary['mns']['status'] == 'analyzed' and summary['mns']['phenotype']:
        conf = summary['mns'].get('confidence') or ''
        conf_tag = f" ({conf})" if conf else ''
        amb = " ⚠" if summary['mns'].get('phase_ambiguous') else ''
        result_parts.append(f"**MNS: {summary['mns']['phenotype']}{conf_tag}{amb}**")
    elif summary['mns']['status'] == 'indeterminate':
        result_parts.append("**MNS: Indeterminate**")
    else:
        result_parts.append("*MNS: Not analyzed*")

    # Diego result — single-axis system, no phase ambiguity flag.
    if summary['di']['status'] == 'analyzed' and summary['di']['phenotype']:
        conf = summary['di'].get('confidence') or ''
        conf_tag = f" ({conf})" if conf else ''
        result_parts.append(f"**Diego: {summary['di']['phenotype']}{conf_tag}**")
    elif summary['di']['status'] == 'indeterminate':
        result_parts.append("**Diego: Indeterminate**")
    else:
        result_parts.append("*Diego: Not analyzed*")

    # Display the result
    final_result = " | ".join(result_parts)

    # Create a styled box for the result
    st.markdown(f"""
    <div style="
        background-color: #f0f8ff;
        border: 3px solid #4CAF50;
        border-radius: 10px;
        padding: 20px;
        margin: 10px 0;
        text-align: center;
        font-size: 24px;
        font-weight: bold;
        color: #2E7D32;
    ">
        {final_result}
    </div>
    """, unsafe_allow_html=True)

    # Additional details in expandable section
    with st.expander("📋 Detailed Analysis Results", expanded=False):
        col1, col2, col3, col4, col5, col6, col7, col8, col9 = st.columns(9)

        with col1:
            st.markdown("### ABO System")
            if summary['abo']['status'] == 'analyzed':
                if summary['abo']['phenotype']:
                    st.success(f"Phenotype: {summary['abo']['phenotype']}")
                if summary['abo']['alleles']:
                    st.info(f"Matching alleles: {', '.join(summary['abo']['alleles'])}")
                if summary['abo']['variants']:
                    st.write(f"**Variants detected:** {', '.join(summary['abo']['variants'])}")
            elif summary['abo']['status'] == 'no_variants':
                st.info("No variants detected - Reference phenotype")
            elif summary['abo']['status'] == 'missing_data':
                st.warning("ABO analysis not performed: FASTA file required for ABO.")
                if summary['abo']['message']:
                    st.write(summary['abo']['message'])
            elif summary['abo']['status'] == 'failed':
                st.error("ABO analysis failed to produce results.")
                if summary['abo']['message']:
                    st.write(summary['abo']['message'])
            else:
                st.warning("ABO analysis not performed")

        with col2:
            st.markdown("### RHD System")
            if summary['rhd']['status'] == 'analyzed':
                if summary['rhd']['phenotype']:
                    st.success(f"Phenotype: {summary['rhd']['phenotype']}")
                if summary['rhd'].get('identity') is not None:
                    st.write(f"Identity: {summary['rhd']['identity']}%")
                if summary['rhd'].get('query_length') is not None and summary['rhd'].get('reference_length') is not None:
                    st.write(f"Sequence length: {summary['rhd']['query_length']} bp vs reference {summary['rhd']['reference_length']} bp")
                if summary['rhd']['alleles']:
                    st.info(f"Matching alleles: {', '.join(summary['rhd']['alleles'])}")
                if summary['rhd']['variants']:
                    st.write(f"**Variants detected:** {', '.join(summary['rhd']['variants'])}")
                else:
                    st.write("**Variants:** None (Wild-type/Reference)")
                if summary['rhd'].get('note'):
                    st.warning(summary['rhd']['note'])
            else:
                st.warning("RHD analysis not performed")

        with col3:
            st.markdown("### RHCE System")
            if summary['rhce']['status'] == 'analyzed':
                st.success(f"Phenotype: {summary['rhce']['phenotype']}")
                if summary['rhce'].get('c_e') or summary['rhce'].get('big_E'):
                    st.write(f"**C/c:** {summary['rhce'].get('c_e') or '-'}  |  "
                             f"**E/e:** {summary['rhce'].get('big_E') or '-'}")
                if summary['rhce'].get('confidence'):
                    st.write(f"**Confidence:** {summary['rhce']['confidence']}")
                if summary['rhce'].get('reads_total'):
                    st.write(f"**Reads:** {summary['rhce']['reads_callable']} callable / "
                             f"{summary['rhce']['reads_total']} total")
                if summary['rhce']['allele_options']:
                    options = ", ".join(o['isbt'] for o in summary['rhce']['allele_options'])
                    st.info(f"ISBT haplotype options: {options}")
                if summary['rhce']['partial_markers']:
                    sigs = "; ".join(p['significance'] for p in summary['rhce']['partial_markers'])
                    st.warning(f"Asian partial-E flagged: {sigs}")
            elif summary['rhce']['status'] == 'indeterminate':
                st.warning("RHCE indeterminate (no callable reads or primary markers missed)")
                if summary['rhce'].get('message'):
                    st.caption(summary['rhce']['message'])
            else:
                st.warning("RHCE analysis not performed")

        with col4:
            st.markdown("### Kell System")
            if summary['kel']['status'] == 'analyzed':
                st.success(f"Phenotype: {summary['kel']['phenotype']}")
                if summary['kel'].get('k_axis'):
                    st.write(f"**K/k:** {summary['kel']['k_axis']}")
                if summary['kel'].get('confidence'):
                    st.write(f"**Confidence:** {summary['kel']['confidence']}")
                if summary['kel'].get('reads_total'):
                    st.write(f"**Reads:** {summary['kel']['reads_callable']} callable / "
                             f"{summary['kel']['reads_total']} total")
                if summary['kel']['allele_options']:
                    opt = summary['kel']['allele_options'][0]
                    serology = opt.get('serology', '')
                    serology_tag = f" — {serology}" if serology else ''
                    st.info(f"ISBT haplotype: {opt['isbt']}{serology_tag}")
            elif summary['kel']['status'] == 'indeterminate':
                st.warning("Kell indeterminate (primary marker c.578 not callable)")
                if summary['kel'].get('message'):
                    st.caption(summary['kel']['message'])
            else:
                st.warning("Kell analysis not performed")

        with col5:
            st.markdown("### Duffy System")
            if summary['fy']['status'] == 'analyzed':
                st.success(f"Phenotype: {summary['fy']['phenotype']}")
                if summary['fy'].get('a_b'):
                    st.write(f"**FY*A/FY*B:** {summary['fy']['a_b']}")
                if summary['fy'].get('gata'):
                    st.write(f"**GATA -67:** {summary['fy']['gata']}")
                if summary['fy'].get('fy_x'):
                    st.write(f"**FY*X (c.265):** {summary['fy']['fy_x']}")
                if summary['fy'].get('confidence'):
                    st.write(f"**Confidence:** {summary['fy']['confidence']}")
                if summary['fy'].get('reads_total'):
                    st.write(f"**Reads:** {summary['fy']['reads_callable']} callable / "
                             f"{summary['fy']['reads_total']} total")
                if summary['fy']['allele_options']:
                    if summary['fy'].get('phase_ambiguous'):
                        st.warning("Phase-ambiguous — multiple haplotype options:")
                    for opt in summary['fy']['allele_options']:
                        serology = opt.get('serology', '')
                        serology_tag = f" — {serology}" if serology else ''
                        st.info(f"ISBT haplotype: {opt['isbt']}{serology_tag}")
            elif summary['fy']['status'] == 'indeterminate':
                st.warning("Duffy indeterminate (primary marker c.125 not callable)")
                if summary['fy'].get('message'):
                    st.caption(summary['fy']['message'])
            else:
                st.warning("Duffy analysis not performed")

        with col6:
            st.markdown("### Kidd System")
            if summary['jk']['status'] == 'analyzed':
                st.success(f"Phenotype: {summary['jk']['phenotype']}")
                if summary['jk'].get('a_b'):
                    st.write(f"**JK*A/JK*B:** {summary['jk']['a_b']}")
                if summary['jk'].get('null_state'):
                    st.write(f"**Null state:** {summary['jk']['null_state']}")
                if summary['jk'].get('null_triggers'):
                    triggers = ", ".join(t['snp'] for t in summary['jk']['null_triggers'])
                    st.caption(f"Triggers: {triggers}")
                if summary['jk'].get('confidence'):
                    st.write(f"**Confidence:** {summary['jk']['confidence']}")
                if summary['jk'].get('reads_total'):
                    st.write(f"**Reads:** {summary['jk']['reads_callable']} callable / "
                             f"{summary['jk']['reads_total']} total")
                if summary['jk']['allele_options']:
                    if summary['jk'].get('phase_ambiguous'):
                        st.warning("Phase-ambiguous — multiple haplotype options:")
                    for opt in summary['jk']['allele_options']:
                        serology = opt.get('serology', '')
                        serology_tag = f" — {serology}" if serology else ''
                        st.info(f"ISBT haplotype: {opt['isbt']}{serology_tag}")
            elif summary['jk']['status'] == 'indeterminate':
                st.warning("Kidd indeterminate (primary marker c.838 not callable)")
                if summary['jk'].get('message'):
                    st.caption(summary['jk']['message'])
            else:
                st.warning("Kidd analysis not performed")

        with col7:
            st.markdown("### H System")
            if summary['h']['status'] == 'analyzed':
                if summary['h'].get('bombay_warning'):
                    st.error(f"🚨 {summary['h']['phenotype']}")
                else:
                    st.success(f"Phenotype: {summary['h']['phenotype']}")
                if summary['h'].get('state'):
                    st.write(f"**State:** {summary['h']['state']}")
                if summary['h'].get('triggers_hom_alt'):
                    triggers = ", ".join(t['snp'] for t in summary['h']['triggers_hom_alt'])
                    st.caption(f"Hom-alt: {triggers}")
                if summary['h'].get('triggers_het'):
                    triggers = ", ".join(t['snp'] for t in summary['h']['triggers_het'])
                    st.caption(f"Het: {triggers}")
                if summary['h'].get('confidence'):
                    st.write(f"**Confidence:** {summary['h']['confidence']}")
                if summary['h'].get('reads_total'):
                    st.write(f"**Reads:** {summary['h']['reads_callable']} callable / "
                             f"{summary['h']['reads_total']} total")
                if summary['h']['allele_options']:
                    for opt in summary['h']['allele_options']:
                        serology = opt.get('serology', '')
                        serology_tag = f" — {serology}" if serology else ''
                        if summary['h'].get('bombay_warning'):
                            st.error(f"ISBT haplotype: {opt['isbt']}{serology_tag}")
                        else:
                            st.info(f"ISBT haplotype: {opt['isbt']}{serology_tag}")
            elif summary['h']['status'] == 'indeterminate':
                st.warning("H indeterminate (no diagnostic SNPs callable)")
                if summary['h'].get('message'):
                    st.caption(summary['h']['message'])
            else:
                st.warning("H analysis not performed")

        with col8:
            st.markdown("### MNS System")
            if summary['mns']['status'] == 'analyzed':
                st.success(f"Phenotype: {summary['mns']['phenotype']}")
                if summary['mns'].get('m_n'):
                    st.write(f"**M/N:** {summary['mns']['m_n']}")
                if summary['mns'].get('s_s'):
                    st.write(f"**S/s:** {summary['mns']['s_s']}")
                if summary['mns'].get('confidence'):
                    st.write(f"**Confidence:** {summary['mns']['confidence']}")
                if summary['mns'].get('reads_total'):
                    st.write(f"**Reads:** {summary['mns']['reads_callable']} callable / "
                             f"{summary['mns']['reads_total']} total")
                if summary['mns']['allele_options']:
                    if summary['mns'].get('phase_ambiguous'):
                        st.warning("Phase-ambiguous — multiple haplotype options:")
                    for opt in summary['mns']['allele_options']:
                        serology = opt.get('serology', '')
                        serology_tag = f" — {serology}" if serology else ''
                        st.info(f"ISBT haplotype: {opt['isbt']}{serology_tag}")
            elif summary['mns']['status'] == 'indeterminate':
                st.warning("MNS indeterminate (no primary marker callable)")
                if summary['mns'].get('message'):
                    st.caption(summary['mns']['message'])
            else:
                st.warning("MNS analysis not performed")

        with col9:
            st.markdown("### Diego System")
            if summary['di']['status'] == 'analyzed':
                st.success(f"Phenotype: {summary['di']['phenotype']}")
                if summary['di'].get('genotype'):
                    st.write(f"**Di(a)/Di(b):** {summary['di']['genotype']}")
                if summary['di'].get('confidence'):
                    st.write(f"**Confidence:** {summary['di']['confidence']}")
                if summary['di'].get('reads_total'):
                    st.write(f"**Reads:** {summary['di']['reads_callable']} callable / "
                             f"{summary['di']['reads_total']} total")
                if summary['di']['allele_options']:
                    opt = summary['di']['allele_options'][0]
                    serology = opt.get('serology', '')
                    serology_tag = f" — {serology}" if serology else ''
                    st.info(f"ISBT haplotype: {opt['isbt']}{serology_tag}")
            elif summary['di']['status'] == 'indeterminate':
                st.warning("Diego indeterminate (primary marker c.2561 not callable)")
                if summary['di'].get('message'):
                    st.caption(summary['di']['message'])
            else:
                st.warning("Diego analysis not performed")


# --- Main Panel ---
st.title("Genetic Analysis Dashboard")

# How to Use the App section
with st.expander("📖 How to Use the App", expanded=False):
    st.markdown("""
    ### Quick Start Guide
    
    **1. Upload your sample file**
    - Drag & drop your .ab1 (Sanger) or .fasta file in the sidebar
    - For RHD analysis, upload a single AB1 or FASTA file
    - For ABO analysis, upload multiple exon-specific FASTA files
    
    **2. Select blood group system**
    - Choose from the left sidebar: ABO, RHD, RHCE, Kell, Duffy, Kidd, H
    - Each system has specific upload requirements
    
    **3. Configure analysis settings**
    - Adjust heterozygosity threshold if needed (default: 0.3)
    - Set exon ranges for targeted analysis (optional)
    
    **4. Run analysis**
    - Click the "Analyze" button to start processing
    - The app automatically aligns sequences, detects variants, and reports blood group systems
    
    **5. View results**
    - Check chromatogram analysis for heterozygotes
    - Review exon-based SNPs and variant details
    - See automatic blood group suggestions
    - Browse ISBT allele database information
    
    **6. Export data**
    - Download results as CSV files for further analysis
    """)

# Tabs for displaying results
tab1, tab2, tab3, tab4 = st.tabs([
    "Chromatogram Check for Heterozygotes",
    "Exon-based SNP",
    "Allele Prediction",
    "Teams & References"
])

with tab4:
    st.markdown("### 📚 References")
    st.markdown("""
            Sirikul C, Wita R, Anukul N.,
                "Assessment of a New ABO Blood Group Genotyping Online Platform."
                , J Hemato Transfus Med. Vol. 35  Supplement 2025. ISSN 2985-2404 (online).")
    """)
    st.markdown("""
    S. Prananpaeng, T. Thaiyanto, R. Wita and N. Anukul, 
    "Whole-Exome Sequencing (WES) Analysis for ABO Subgroups Identification," 
    *2023 20th International Joint Conference on Computer Science and Software Engineering (JCSSE)*, 
    2023, pp. 264-269, url: [https://ieeexplore.ieee.org/document/10202117/](https://ieeexplore.ieee.org/document/10202117/).
    """)
    st.markdown("""
    R. Wita, S. Somhom, J. Chawachat, A. Thongratsameethong, N. Anukul and C. Sirikul, 
    "DNA Sequencing Analysis Framework for ABO Genotyping and ABO Discrepancy Resolution," 
    *2021 18th International Conference on Electrical Engineering/Electronics, Computer, 
    Telecommunications and Information Technology (ECTI-CON)*, Chiang Mai, Thailand, 2021, 
    pp. 913-916, doi: [10.1109/ECTI-CON51831.2021.9454861](https://ieeexplore.ieee.org/document/9454861).
    """)
    st.markdown("""
            Anukul N, Wita R, Leetrakool N, Sirikul C, Veeraphan N, Wongchai S. 
             "Two novel alleles on Fucosyltransferase 2 from northern Thai paraBombay 
             family and computational prediction on mutation effect."
              Transfusion. 2021;1–11. [https://doi.org/10.1111/trf.16646](https://doi.org/10.1111/trf.16646)
             """)

    st.subheader("👥 Teams")
    st.markdown("#### InnoGeHLA Lab")
    st.markdown("""
The InnoGeHLA (**Inno**vation **Ge**nomics **HLA**) Lab at Chiang Mai University is a collaborative research group 
                that brings together expertise in genomics, bioinformatics, and transfusion science. 
                The lab’s mission is to blend molecular biology with computational innovation to enhance 
                precision medicine and genetic diagnostics, with a special focus on blood group and HLA 
                genotyping.

The team is led by **Asst. Prof. Nampeung Anukul**, whose work centers on blood group genetics and transfusion 
                science, and **Miss Chonthicha Sirikul**, who specializes in molecular diagnostics and hematology.
                 Joining them is **Asst. Prof. Ratsameetip Wita**, a computer scientist from the Faculty of 
                Science who explores the use of artificial intelligence and bioinformatics to interpret 
                complex genetic data.

Together, they form a dynamic and interdisciplinary team that bridges the gap between biomedical research 
                and computational science—working to turn genomic insights into practical tools that 
                benefit both laboratories and patients.
                """)
    st.markdown("##### Lab members")

    # CSS for circular images
    st.markdown("""
        <style>
        .member-card {
            text-align: center;
            padding: 20px;
        }
        
        .member-name {
            font-weight: bold;
            font-size: 18px;
            margin-top: 15px;
            color: #333;
            text-align: center;
        }
        
        .member-position {
            font-size: 14px;
            color: #666;
            margin-top: 5px;
            line-height: 1.6;
            text-align: center;
        }
        
        /* Center the image within its container */
        .element-container:has(div[data-testid="stImage"]) {
            display: flex;
            justify-content: center;
            align-items: center;
        }
        
        div[data-testid="stImage"] {
            text-align: center;
        }
        
        /* Remove width fit-content from Streamlit's emotion cache class */
        .st-emotion-cache-p75nl5 {
            width: 100% !important;
        }
        
        /* Make Streamlit images circular */
        div[data-testid="stImage"] img {
            border-radius: 50%;
            width: 200px;
            height: 200px;
            object-fit: cover;
            object-position: center;
            box-shadow: 0 4px 8px rgba(0, 0, 0, 0.2);
            display: inline-block;
        }
        </style>
    """, unsafe_allow_html=True)

    

    # Create three columns for team members
    col1, col2, col3 = st.columns(3)

    with col1:
        st.markdown('<div class="member-card">', unsafe_allow_html=True)
        try:
            from PIL import Image
            img = Image.open("utils/img/NP.jpg")
            st.image(img, width=200)
        except Exception as e:
            st.image("https://via.placeholder.com/200", width=200)
        st.markdown('<div class="member-name">Nampeung Anukul</div>',
                    unsafe_allow_html=True)
        st.markdown('''<div class="member-position">
            Division of
Blood Transfusion Science,<br>Faculty of Associated Medical Sciences<br>
            Chiang Mai University<br>
            nampeung.a@cmu.ac.th
        </div>''', unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)

    with col2:
        st.markdown('<div class="member-card">', unsafe_allow_html=True)
        try:
            from PIL import Image
            img = Image.open("utils/img/CS.jpg")
            st.image(img, width=200)
        except Exception as e:
            st.image("https://via.placeholder.com/200", width=200)
        st.markdown(
            '<div class="member-name">Chonthicha Sirikul</div>', unsafe_allow_html=True)
        st.markdown('''<div class="member-position">
            Division of
Blood Transfusion Science,<br> Faculty of Associated Medical Sciences<br>
            Chiang Mai University<br>
            chonthicha.sir@cmu.ac.th
        </div>''', unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)

    with col3:
        st.markdown('<div class="member-card">', unsafe_allow_html=True)
        try:
            from PIL import Image
            img = Image.open("utils/img/RW.jpg")
            st.image(img, width=200)
        except Exception as e:
            st.image("https://via.placeholder.com/200", width=200)
        st.markdown('<div class="member-name">Ratsameetip Wita</div>',
                    unsafe_allow_html=True)
        st.markdown('''<div class="member-position">
            Department of Computer Science,<br> Faculty of Science<br>
            Chiang Mai University<br>
            ratsameetip.wit@cmu.ac.th
        </div>''', unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)
    # Student members table - read from CSV file

    st.markdown("##### Student Contributors")
    import pandas as pd
    try:
        student_data = pd.read_csv("utils/data/students.csv")
        st.dataframe(student_data, hide_index=True)
    except FileNotFoundError:
        st.warning(
            "Student data file not found. Please add utils/data/students.csv")
    except Exception as e:
        st.error(f"Error loading student data: {e}")

if analyze_button:
    if (not fwd_ab1 and not fasta_files and not rhd_fasta_files and not rhd_ab1_files
            and not rhce_ab1_files and not rhce_fasta_files
            and not kel_ab1_files and not kel_fasta_files
            and not fy_ab1_files and not fy_fasta_files
            and not jk_ab1_files and not jk_fasta_files
            and not h_ab1_files and not h_fasta_files
            and not mns_ab1_files and not mns_fasta_files
            and not di_ab1_files and not di_fasta_files
            and not unified_ab1_files and not unified_fasta_files):
        st.warning("Please upload at least one file before analyzing.")
    else:
        status_container = st.empty()
        status_container.success("Files uploaded successfully! Starting analysis...")

        # ─── Unified auto-detect upload: filename-based routing ─────────
        # Each unified file is routed by its filename to a blood-group
        # system, then appended to that system's per-system upload list so
        # the existing downstream processing handles it transparently. The
        # routing_decisions log is rendered above the results banner.
        routing_decisions: list = []   # [{filename, kind, system_or_None}]
        _routed: dict = {
            'ABO': {'ab1': [], 'fasta': []},
            'RHD': {'ab1': [], 'fasta': []},
            'RHCE': {'ab1': [], 'fasta': []},
            'KEL': {'ab1': [], 'fasta': []},
            'FY':  {'ab1': [], 'fasta': []},
            'JK':  {'ab1': [], 'fasta': []},
            'H':   {'ab1': [], 'fasta': []},
            'MNS': {'ab1': [], 'fasta': []},
            'DI':  {'ab1': [], 'fasta': []},
        }
        for _f in (unified_ab1_files or []):
            _sys = route_filename(_f.name)
            routing_decisions.append({
                'filename': _f.name, 'kind': 'ab1', 'system': _sys,
            })
            if _sys in _routed:
                _routed[_sys]['ab1'].append(_f)
        for _f in (unified_fasta_files or []):
            _sys = route_filename(_f.name)
            routing_decisions.append({
                'filename': _f.name, 'kind': 'fasta', 'system': _sys,
            })
            if _sys in _routed:
                _routed[_sys]['fasta'].append(_f)

        # Merge routed files into the per-system uploader variables so the
        # existing downstream pipeline (process_*_files, analyzers, summary,
        # banner, detail-view) doesn't need to know unified upload exists.
        fwd_ab1        = list(fwd_ab1        or []) + _routed['ABO']['ab1']
        fasta_files    = list(fasta_files    or []) + _routed['ABO']['fasta']
        rhd_ab1_files  = list(rhd_ab1_files  or []) + _routed['RHD']['ab1']
        rhd_fasta_files = list(rhd_fasta_files or []) + _routed['RHD']['fasta']
        rhce_ab1_files = list(rhce_ab1_files or []) + _routed['RHCE']['ab1']
        rhce_fasta_files = list(rhce_fasta_files or []) + _routed['RHCE']['fasta']
        kel_ab1_files  = list(kel_ab1_files  or []) + _routed['KEL']['ab1']
        kel_fasta_files = list(kel_fasta_files or []) + _routed['KEL']['fasta']
        fy_ab1_files   = list(fy_ab1_files   or []) + _routed['FY']['ab1']
        fy_fasta_files = list(fy_fasta_files or []) + _routed['FY']['fasta']
        jk_ab1_files   = list(jk_ab1_files   or []) + _routed['JK']['ab1']
        jk_fasta_files = list(jk_fasta_files or []) + _routed['JK']['fasta']
        h_ab1_files    = list(h_ab1_files    or []) + _routed['H']['ab1']
        h_fasta_files  = list(h_fasta_files  or []) + _routed['H']['fasta']
        mns_ab1_files  = list(mns_ab1_files  or []) + _routed['MNS']['ab1']
        mns_fasta_files = list(mns_fasta_files or []) + _routed['MNS']['fasta']
        di_ab1_files   = list(di_ab1_files   or []) + _routed['DI']['ab1']
        di_fasta_files = list(di_fasta_files or []) + _routed['DI']['fasta']

        # ABO AB1 channel -> chromatogram + heterozygote detection
        processed_AB1, hets = process_ab1_files(
            fwd_ab1, [], threshold_ratio) if fwd_ab1 else (None, None)

        # Build synthetic FASTA file objects from AB1 basecalled sequences so
        # they flow through the same variant/allele pipeline as user FASTAs.
        ab1_derived_fastas = []
        if processed_AB1:
            for trace in processed_AB1:
                seq = trace.get('seq', '')
                if not seq or len(seq) < 50:
                    continue
                fname = trace.get('filename', 'ab1_seq') + ".fasta"
                fasta_text = f">{fname}\n{seq}\n".encode('utf-8')
                synth = io.BytesIO(fasta_text)
                synth.name = fname
                ab1_derived_fastas.append(synth)

        # Combine user FASTAs and AB1-derived FASTAs for the ABO pipeline
        abo_inputs_for_fasta = []
        if fasta_files:
            abo_inputs_for_fasta.extend(fasta_files)
        if ab1_derived_fastas:
            abo_inputs_for_fasta.extend(ab1_derived_fastas)

        # --- Robust Processing Logic ---
        robust_summary = []
        service = None
        if abo_inputs_for_fasta:
            service = fasta_utils.FASTAAlignmentService()
            robust_summary = service.generate_batch_summary(abo_inputs_for_fasta)

        confirmed_results = [r for r in robust_summary if "Confirmed" in r.get('decision', "")]

        # Build exons_ref from confirmed results so AB1 exon labels can map to CDS coords
        if confirmed_results and service is not None:
            aboRef = service.getABO_ref("exons")
            exons_ref = []
            for r in confirmed_results:
                exon_num = r['exon_number']
                x = exon_num - 1
                if x < len(aboRef):
                    exons_ref.append({
                        'exon': exon_num,
                        'ref_start': r.get('ref_start', 0),
                        'ref_end': r.get('ref_end', 0),
                        'cds_start': aboRef[x]['cds_start'],
                        'cds_end': aboRef[x]['cds_end']
                    })

        # RHD AB1 channel -> RHD analyzer only (no chromatogram tab).
        # Uses RHD-specific processor with Phred quality trimming and
        # signal-based heterozygous IUPAC encoding for diagnostic SNPs.
        rhd_ab1_traces, _ = process_rhd_ab1_files(
            rhd_ab1_files, het_ratio=threshold_ratio) if rhd_ab1_files else (None, None)

        # RHD FASTA channel -> RHD analyzer
        rhd_fasta_traces = process_rhd_fasta_files(rhd_fasta_files) if rhd_fasta_files else None

        # Combine RHD inputs (channel-based, no filename filtering)
        rhd_input_traces = []
        if rhd_ab1_traces:
            rhd_input_traces.extend(rhd_ab1_traces)
        if rhd_fasta_traces:
            rhd_input_traces.extend(rhd_fasta_traces)

        # RHCE inputs - reuse the RHD AB1/FASTA processors (same trace shape).
        rhce_ab1_traces, _ = process_rhd_ab1_files(
            rhce_ab1_files, het_ratio=threshold_ratio) if rhce_ab1_files else (None, None)
        rhce_fasta_traces = process_rhd_fasta_files(rhce_fasta_files) if rhce_fasta_files else None

        rhce_input_traces = []
        if rhce_ab1_traces:
            rhce_input_traces.extend(rhce_ab1_traces)
        if rhce_fasta_traces:
            rhce_input_traces.extend(rhce_fasta_traces)

        # Run RHCE consensus once so both the headline banner and the
        # detailed RHCE section share the same result. None if no inputs
        # or if the reference is missing.
        rhce_result = None
        rhce_reads_for_summary = []
        rhce_error_message = None
        if rhce_input_traces:
            for i, trace in enumerate(rhce_input_traces):
                if isinstance(trace, dict) and 'seq' in trace:
                    seq = trace.get('seq', '')
                    if seq and len(seq) > 50:
                        rid = trace.get('filename') or f"RHCE_read_{i+1}"
                        qual = trace.get('quality_trimmed')
                        rhce_reads_for_summary.append((rid, seq, qual))
            if rhce_reads_for_summary:
                try:
                    rhce_result = RHCEAnalyzer().analyze(rhce_reads_for_summary)
                except RHCEReferenceMissingError as e:
                    rhce_error_message = str(e)

        # KEL inputs — reuse the RHD AB1/FASTA processors (same trace shape).
        kel_ab1_traces, _ = process_rhd_ab1_files(
            kel_ab1_files, het_ratio=threshold_ratio) if kel_ab1_files else (None, None)
        kel_fasta_traces = process_rhd_fasta_files(kel_fasta_files) if kel_fasta_files else None

        kel_input_traces = []
        if kel_ab1_traces:
            kel_input_traces.extend(kel_ab1_traces)
        if kel_fasta_traces:
            kel_input_traces.extend(kel_fasta_traces)

        kel_result = None
        kel_reads_for_summary = []
        kel_error_message = None
        if kel_input_traces:
            for i, trace in enumerate(kel_input_traces):
                if isinstance(trace, dict) and 'seq' in trace:
                    seq = trace.get('seq', '')
                    if seq and len(seq) > 50:
                        rid = trace.get('filename') or f"KEL_read_{i+1}"
                        # Per-base Phred quality (AB1 inputs) or None (FASTA)
                        qual = trace.get('quality_trimmed')
                        kel_reads_for_summary.append((rid, seq, qual))
            if kel_reads_for_summary:
                try:
                    kel_result = KELAnalyzer().analyze(kel_reads_for_summary)
                except KELReferenceMissingError as e:
                    kel_error_message = str(e)

        # FY (Duffy) inputs — reuse the RHD AB1/FASTA processors (same trace shape).
        fy_ab1_traces, _ = process_rhd_ab1_files(
            fy_ab1_files, het_ratio=threshold_ratio) if fy_ab1_files else (None, None)
        fy_fasta_traces = process_rhd_fasta_files(fy_fasta_files) if fy_fasta_files else None

        fy_input_traces = []
        if fy_ab1_traces:
            fy_input_traces.extend(fy_ab1_traces)
        if fy_fasta_traces:
            fy_input_traces.extend(fy_fasta_traces)

        fy_result = None
        fy_reads_for_summary = []
        fy_error_message = None
        if fy_input_traces:
            for i, trace in enumerate(fy_input_traces):
                if isinstance(trace, dict) and 'seq' in trace:
                    seq = trace.get('seq', '')
                    if seq and len(seq) > 50:
                        rid = trace.get('filename') or f"FY_read_{i+1}"
                        qual = trace.get('quality_trimmed')
                        fy_reads_for_summary.append((rid, seq, qual))
            if fy_reads_for_summary:
                try:
                    fy_result = FYAnalyzer().analyze(fy_reads_for_summary)
                except FYReferenceMissingError as e:
                    fy_error_message = str(e)

        # JK (Kidd) inputs — reuse the RHD AB1/FASTA processors (same trace shape).
        jk_ab1_traces, _ = process_rhd_ab1_files(
            jk_ab1_files, het_ratio=threshold_ratio) if jk_ab1_files else (None, None)
        jk_fasta_traces = process_rhd_fasta_files(jk_fasta_files) if jk_fasta_files else None

        jk_input_traces = []
        if jk_ab1_traces:
            jk_input_traces.extend(jk_ab1_traces)
        if jk_fasta_traces:
            jk_input_traces.extend(jk_fasta_traces)

        jk_result = None
        jk_reads_for_summary = []
        jk_error_message = None
        if jk_input_traces:
            for i, trace in enumerate(jk_input_traces):
                if isinstance(trace, dict) and 'seq' in trace:
                    seq = trace.get('seq', '')
                    if seq and len(seq) > 50:
                        rid = trace.get('filename') or f"JK_read_{i+1}"
                        qual = trace.get('quality_trimmed')
                        jk_reads_for_summary.append((rid, seq, qual))
            if jk_reads_for_summary:
                try:
                    jk_result = JKAnalyzer().analyze(jk_reads_for_summary)
                except JKReferenceMissingError as e:
                    jk_error_message = str(e)

        # H (FUT1 / Bombay) inputs — reuse the RHD AB1/FASTA processors.
        h_ab1_traces, _ = process_rhd_ab1_files(
            h_ab1_files, het_ratio=threshold_ratio) if h_ab1_files else (None, None)
        h_fasta_traces = process_rhd_fasta_files(h_fasta_files) if h_fasta_files else None

        h_input_traces = []
        if h_ab1_traces:
            h_input_traces.extend(h_ab1_traces)
        if h_fasta_traces:
            h_input_traces.extend(h_fasta_traces)

        h_result = None
        h_reads_for_summary = []
        h_error_message = None
        if h_input_traces:
            for i, trace in enumerate(h_input_traces):
                if isinstance(trace, dict) and 'seq' in trace:
                    seq = trace.get('seq', '')
                    if seq and len(seq) > 50:
                        rid = trace.get('filename') or f"H_read_{i+1}"
                        qual = trace.get('quality_trimmed')
                        h_reads_for_summary.append((rid, seq, qual))
            if h_reads_for_summary:
                try:
                    h_result = HAnalyzer().analyze(h_reads_for_summary)
                except HReferenceMissingError as e:
                    h_error_message = str(e)

        # MNS (GYPA + GYPB) inputs — reuse the RHD AB1/FASTA processors.
        mns_ab1_traces, _ = process_rhd_ab1_files(
            mns_ab1_files, het_ratio=threshold_ratio) if mns_ab1_files else (None, None)
        mns_fasta_traces = process_rhd_fasta_files(mns_fasta_files) if mns_fasta_files else None

        mns_input_traces = []
        if mns_ab1_traces:
            mns_input_traces.extend(mns_ab1_traces)
        if mns_fasta_traces:
            mns_input_traces.extend(mns_fasta_traces)

        mns_result = None
        mns_reads_for_summary = []
        mns_error_message = None
        if mns_input_traces:
            for i, trace in enumerate(mns_input_traces):
                if isinstance(trace, dict) and 'seq' in trace:
                    seq = trace.get('seq', '')
                    if seq and len(seq) > 50:
                        rid = trace.get('filename') or f"MNS_read_{i+1}"
                        qual = trace.get('quality_trimmed')
                        mns_reads_for_summary.append((rid, seq, qual))
            if mns_reads_for_summary:
                try:
                    mns_result = MNSAnalyzer().analyze(mns_reads_for_summary)
                except MNSReferenceMissingError as e:
                    mns_error_message = str(e)

        # Diego (SLC4A1) inputs — reuse the RHD AB1/FASTA processors.
        di_ab1_traces, _ = process_rhd_ab1_files(
            di_ab1_files, het_ratio=threshold_ratio) if di_ab1_files else (None, None)
        di_fasta_traces = process_rhd_fasta_files(di_fasta_files) if di_fasta_files else None

        di_input_traces = []
        if di_ab1_traces:
            di_input_traces.extend(di_ab1_traces)
        if di_fasta_traces:
            di_input_traces.extend(di_fasta_traces)

        di_result = None
        di_reads_for_summary = []
        di_error_message = None
        if di_input_traces:
            for i, trace in enumerate(di_input_traces):
                if isinstance(trace, dict) and 'seq' in trace:
                    seq = trace.get('seq', '')
                    if seq and len(seq) > 50:
                        rid = trace.get('filename') or f"DI_read_{i+1}"
                        qual = trace.get('quality_trimmed')
                        di_reads_for_summary.append((rid, seq, qual))
            if di_reads_for_summary:
                try:
                    di_result = DIAnalyzer().analyze(di_reads_for_summary)
                except DIReferenceMissingError as e:
                    di_error_message = str(e)

        status_container.empty()

        # === Routing decisions from the unified auto-detect upload ===
        # Only rendered when the unified upload was actually used. Shows
        # which system each file was dispatched to, and flags any files
        # whose names didn't match a known amplicon convention.
        if routing_decisions:
            n_routed = sum(1 for r in routing_decisions if r['system'])
            n_unrouted = sum(1 for r in routing_decisions if not r['system'])
            st.markdown(f"### 🚀 Auto-detect routing  "
                        f"({n_routed} routed, {n_unrouted} unrouted)")
            _routing_rows = [{
                'File':   r['filename'],
                'Kind':   r['kind'].upper(),
                'Routed to': r['system'] or '— (unrouted, skipped)',
                'Status': '✓ Routed' if r['system'] else '⚠ Name not recognised',
            } for r in routing_decisions]
            st.dataframe(pd.DataFrame(_routing_rows),
                         hide_index=True, use_container_width=True)
            if n_unrouted:
                st.caption(
                    "Unrouted files were not analyzed. Rename them to include "
                    "a recognised pattern (RHCE, RHD, ABO, KEL, FY, JK, MIA, "
                    "DI, FUT1) or upload via the system-specific uploaders."
                )

        # === FINAL BLOOD GROUP RESULT banner (consolidated) ===
        try:
            _isbt_handler = ISBTDataHandler()
            _summary = generate_final_blood_group_summary(
                robust_summary,
                processed_AB1,
                hets,
                _isbt_handler,
                has_fasta=bool(fasta_files),
                rhce_result=rhce_result,
                kel_result=kel_result,
                fy_result=fy_result,
                jk_result=jk_result,
                h_result=h_result,
                mns_result=mns_result,
                di_result=di_result,
            )
            display_final_blood_group_result(_summary)
        except Exception as banner_exc:
            st.warning(f"Could not render consolidated summary: {banner_exc}")

        with tab1:
            st.subheader("Chromatogram Check for Heterozygotes")
            st.write("### Heterozygote Positions Detected: ",
                     len(hets) if hets else 0)
            if fasta_files and not fwd_ab1:
                st.info("FASTA results are in the **Exon-based SNP** and **Allele Prediction** tabs.")
            if processed_AB1:
                def _het_to_intensities(h):
                    tb = h.get('top_bases')
                    if tb:
                        return dict(tb)
                    if 'major_signal' in h and 'minor_signal' in h:
                        return {h['ref_base']: h['major_signal'], h['alt_base']: h['minor_signal']}
                    return {h['ref_base']: 1, h['alt_base']: 1}
                hetero_sites_plot = [(h['position'], _het_to_intensities(h)) for h in hets] if hets else None
                for i in processed_AB1:  # type: ignore
                    exon_num = i.get('exon')
                    cds_start, cds_end = get_cds(exon_num) if exon_num is not None else (None, None)
                    fig = plot_chromatogram_plotly(
                        i, base_width=2, cds_start=cds_start, cds_end=cds_end, hetero_sites=hetero_sites_plot)
                    st.plotly_chart(fig, use_container_width=True)
            elif not fwd_ab1:
                st.info("Upload an ABO AB1 file to view the chromatogram.")

            if hets:
                st.write("### Detected Heterozygous Positions:")
                het_data = {
                    "Position": [h['position'] for h in hets],
                    "Top Base": [h['ref_base'] for h in hets],
                    "Alt Base": [h['alt_base'] for h in hets],
                    "Ratio": [round(h['ratio'], 2) for h in hets]
                }
                het_df = pd.DataFrame(het_data)
                st.table(het_df)

        with tab2:
            st.subheader("Exon-based SNP Analysis (Robust)")

            if not robust_summary:
                 st.info("No FASTA files processed.")
            else:
                 # --- 1. Summary Table ---
                 st.write("### 📄 Analysis Summary Table")
                 
                 summary_data = []
                 for res in robust_summary:
                     summary_data.append({
                         "Filename": res['filename'],
                         "Exon": res['exon'],
                         "Coverage (%)": f"{res['coverage']:.1f}%",
                         "Correctness (%)": f"{res['similarity']:.1f}%",
                         "Final Decision": res['decision']
                     })
                 
                 st.dataframe(pd.DataFrame(summary_data), use_container_width=True)
                 
                 # --- 2. Detailed Alignments for Confirmed Exons ---
                 st.write("### 🔍 Detailed Alignments (Confirmed Only)")
                 
                 found_confirmed = False
                 for res in robust_summary:
                     if "Confirmed" in res['decision']:
                         found_confirmed = True
                         with st.expander(f"✅ {res['filename']} - Exon {res['exon']} Details"):
                             st.write(f"**Orientation**: {res.get('orientation', 'N/A')}")
                             st.write(f"**Score**: {res.get('score', 0)}")
                             
                             # Use existing display function
                             # We need to construct 'aligned_query' and 'aligned_reference' from result if available
                             # identify_variants calls align_sequence_to_exon, which returns aligned_query/ref
                             if 'aligned_query' in res and 'aligned_reference' in res:
                                 display_alignment_with_snps(
                                     aligned_query=res['aligned_query'],
                                     aligned_reference=res['aligned_reference'],
                                     cds_start=res.get('cds_start'), # Might be None if not passed, but alignment handles it
                                     cds_end=res.get('cds_end'),
                                     variants=res.get('variants', []),
                                     exon_number=res['exon_number'],
                                     unique_id=res['filename']
                                 )
                                 
                                 # Variant Table for this file
                                 if res.get('variants'):
                                     st.write("**Variants Identified:**")
                                     v_data = []
                                     for v in res['variants']:
                                         v_data.append({
                                             "ISBT Pos": v.get('isbt_pos'),
                                             "Type": v.get('type'),
                                             "Change": f"{v.get('ref_base', '')}>{v.get('alt_base', '')}" if v.get('type') == 'SNP' else f"{v.get('type')}"
                                         })
                                     st.dataframe(pd.DataFrame(v_data), hide_index=True)
                                 else:
                                     st.success("No variants found (Perfect Match to Reference)")
                             else:
                                 st.warning("Alignment details not available.")
                                 
                 if not found_confirmed:
                     st.warning("No exons were confidently confirmed in the uploaded files.")

        with tab3:
            st.subheader("Allele Prediction & Strand Analysis")
            st.write("Analysis of potential genotype combinations based on exon segregation.")
            st.markdown(
                "**Reference:** [ISBT ABO Alleles Table](https://www.isbtweb.org/resource/001aboalleles.html)")

            # Color code legend
            st.markdown("""
            #### 🎨 Color Code Legend:
            """)

            legend_html = """
            <table style="border-collapse: collapse; margin-bottom: 20px;">
                <tr>
                    <td style="background-color: #FEC98F; padding: 8px 12px; border: 2px solid #999; font-weight: bold;">Heterozygous Variant</td>
                    <td style="padding: 8px 12px; border: 2px solid #999;">Variant detected with IUPAC ambiguity code (heterozygous)</td>
                </tr>
                <tr>
                    <td style="background-color: #E0FFCC; padding: 8px 12px; border: 2px solid #999; font-weight: bold;">Found in Tested Exon</td>
                    <td style="padding: 8px 12px; border: 2px solid #999;">Variant confirmed in the tested exon</td>
                </tr>
                <tr>
                    <td style="background-color: #FFD6C9; padding: 8px 12px; border: 2px solid #999; font-weight: bold;">Not Found in Tested Exon</td>
                    <td style="padding: 8px 12px; border: 2px solid #999;">Variant expected but not detected in the tested exon</td>
                </tr>
            </table>
            """
            st.markdown(legend_html, unsafe_allow_html=True)

            if not confirmed_results:
                 st.info("No confirmed exons to analyze.")
            else:
                 abo_id = abo_utils.ABOIdentifier("ABO")

                 # --- 1. Group Files by Exon ---
                 files_by_exon = {}
                 for res in confirmed_results:
                     ex = res['exon_number']
                     if ex not in files_by_exon:
                         files_by_exon[ex] = []
                     files_by_exon[ex].append(res)
                 
                 sorted_exons = sorted(files_by_exon.keys())
                 st.write(f"**Detected Exons:** {', '.join([f'Exon {e} ({len(files_by_exon[e])} files)' for e in sorted_exons])}")

                 # --- 2. Generate Strand Combinations ---
                 # Cartesian product of files from each exon
                 # data_lists = [files_by_exon[e] for e in sorted_exons]
                 # combinations = list(itertools.product(*data_lists)) # This works if we have at least one file for each exon
                 
                 # Optimization: Only combine relevant exons? 
                 # For now, we combine ALL available exons.
                 
                 strand_candidates = list(itertools.product(*[files_by_exon[e] for e in sorted_exons]))
                 
                 st.write(f"### 🧬 Strand Combinations (Potential Genotypes)")
                 st.write(f"Testing {len(strand_candidates)} possible combinations of files...")
                 
                 valid_strands_found = 0
                 
                 for i, combo in enumerate(strand_candidates):
                     # combo is a tuple of file_results, one per exon
                     combo_name = " + ".join([f"{c['filename']} (E{c['exon_number']})" for c in combo])
                     
                     # Check if we should display this combo
                     # Logic: Run identification on the Union of variants
                     
                     preds, unk, var_names, node_map = identify_abo_alleles(
                         {'exon_alignments': list(combo)}, abo_identifier=abo_id)
                     
                     if preds:
                         valid_strands_found += 1
                         with st.expander(f"✅ Combination {i+1}: matches **{', '.join(list(preds[0].keys()))}** ..."):
                             st.write(f"**Files:** {combo_name}")
                             
                             # Display Prediction Table
                             html_string = "<table>"
                             html_string += "<tr><th>Allele</th><th>Variant Name</th><th>Exon</th><th>Location</th><th>Change</th></tr>"
                             
                             tested_exons_in_combo = [c['exon_number'] for c in combo]

                             for allele_data in preds:
                                for allele_name, allele_variants in allele_data.items():
                                    num_variants = len(allele_variants)
                                    html_string += f"<tr style='border-top: 3px solid #999;'><td rowspan='{num_variants}'><b>{allele_name}</b></td>"
                                    
                                    for idx, variant in enumerate(allele_variants):
                                        hets = False
                                        if variant['name'] in node_map:
                                            iupac_code = node_map[variant['name']]
                                            if iupac_code not in ['A', 'T', 'C', 'G']:
                                                hets = True
                                        
                                        if idx > 0:
                                            html_string += "<tr>"
                                        
                                        style = "background-color: #EFECE6;"
                                        
                                        # Check if variant exon is in our tested combo
                                        if variant['exon'] in tested_exons_in_combo:
                                            if variant['name'] in var_names:
                                                if hets:
                                                     style = "background-color: #FEC98F;"
                                                else:
                                                     style = "background-color: #E0FFCC;"
                                            else:
                                                 style = "background-color: #FFD6C9;"
                                        
                                        html_string += f"<td style='{style}'>{variant['name']}</td><td style='{style}'>{variant['exon']}</td><td style='{style}'>{variant['location']}</td><td style='{style}'>{variant['change']}</td></tr>"
                             
                             html_string += "</table>"
                             st.markdown(html_string, unsafe_allow_html=True)

                 if valid_strands_found == 0:
                     st.warning("No consistent alleles found for any combination of files.")

                 st.write("---")
                 st.write("### 📂 Individual File Analysis (Exon-specific Support)")
                 
                 for res in confirmed_results:
                     filename = res['filename']
                     exon_num = res['exon_number']
                     variants = res['variants']
                     
                     st.markdown(f"#### 📄 {filename} (Exon {exon_num})")
                     
                     if variants:
                         # Identify based on variants
                         preds, unk, var_names, node_map = identify_abo_alleles({'exon_alignments': [res]}, abo_identifier=abo_id)
                         
                         if preds:
                             # Simplified view for single file
                             top_alleles = [list(d.keys())[0] for d in preds[:5]]
                             st.success(f"Supports: {', '.join(top_alleles)}" + ("..." if len(preds)>5 else ""))
                         else:
                            st.warning("Variants detected but no known allele combination matches exactly.")

                         if unk:
                             st.write("**Unknown Variants:**")
                             st.dataframe(pd.DataFrame(unk))

                     else:
                         # Reference Matching Case
                         st.success("✅ No variants detected (Matches Reference)")
                         
                         all_alleles = abo_id.get_all_alleles()
                         excluded = abo_id.get_alleles_with_variant_in_exon(exon_num)
                         msg = f"Consistent with {len(all_alleles)} alleles (Total) - {len(excluded)} (Excluded) = {len(all_alleles)-len(excluded)} Candidates."
                         st.write(msg)
                     st.write("---")

        # === DETAILED SYSTEM ANALYSIS ===
        # Only one blood-group system's detailed view is shown at a time so
        # that results don't overlap. The selectbox lists only systems whose
        # inputs were provided this run.
        detail_systems = []
        if rhd_input_traces:
            detail_systems.append("RHD")
        if rhce_input_traces:
            detail_systems.append("RHCE")
        if kel_input_traces:
            detail_systems.append("Kell")
        if fy_input_traces:
            detail_systems.append("Duffy")
        if jk_input_traces:
            detail_systems.append("Kidd")
        if h_input_traces:
            detail_systems.append("H")
        if mns_input_traces:
            detail_systems.append("MNS")
        if di_input_traces:
            detail_systems.append("Diego")

        if detail_systems:
            st.markdown("---")
            st.header("🩸 Detailed System Analysis")
            selected_system = st.selectbox(
                "Show details for:",
                detail_systems,
                key="detailed_system_selector",
            )

            if selected_system == "RHD":
                st.subheader("RHD Analysis Results")

                # Extract sequences from AB1 traces and RHD FASTA traces for RHD analysis
                # Each trace is a separate amplicon for voting
                rhd_sequences = []
                qc_by_name = {}
                for i, trace in enumerate(rhd_input_traces):
                    if isinstance(trace, dict) and 'seq' in trace:
                        seq = trace.get('seq', '')
                        if seq and len(seq) > 50:
                            filename = trace.get('filename', '')
                            if not filename:
                                filename = f"Amplicon_{i+1}"
                            # Per-base Phred (AB1) or None (FASTA). The RHD
                            # analyzer's analyze_multiple_amplicons accepts
                            # 3-tuples and applies Q30 at SNP positions.
                            qual = trace.get('quality_trimmed')
                            rhd_sequences.append((filename, seq, qual))
                            if 'qc' in trace:
                                qc_by_name[filename] = trace['qc']

                if rhd_sequences:
                    analyzer = RHDAnalyzer()
                    voting_result = analyzer.analyze_multiple_amplicons(rhd_sequences)

                    for r in voting_result['amplicon_results']:
                        qc = qc_by_name.get(r['name'])
                        if qc:
                            r['qc'] = qc

                    st.subheader(f"Multi-Amplicon Analysis ({voting_result['total_amplicons']} amplicon(s))")

                    st.write("### Individual Amplicon Results:")

                    table_data = []
                    for result in voting_result['amplicon_results']:
                        qc = result.get('qc') or {}
                        if qc:
                            qc_summary = (
                                f"5':{qc.get('trimmed_5p', 0)} "
                                f"3':{qc.get('trimmed_3p', 0)} "
                                f"N:{qc.get('masked_internal', 0)} "
                                f"het:{qc.get('het_positions_encoded', 0)}"
                            )
                        else:
                            qc_summary = '-'
                        zyg = result.get('zygosity')
                        table_data.append({
                            'File': result['name'],
                            'Length (bp)': result['length'],
                            'Region': result['region'],
                            'Identity (%)': f"{result['identity']:.1f}%",
                            'Variants': result['variants'],
                            'ISBT Allele': result.get('allele') or '-',
                            'Zygosity': zyg.upper() if zyg else '-',
                            'QC (trim/mask/het)': qc_summary,
                            'Vote': result['vote'],
                        })

                    df_results = pd.DataFrame(table_data)
                    st.dataframe(df_results, use_container_width=True, hide_index=True)

                    col1, col2, col3, col4 = st.columns(4)
                    with col1:
                        st.metric("RhD+ Votes", voting_result['votes']['RhD+'])
                    with col2:
                        st.metric("RhD- Votes", voting_result['votes']['RhD-'])
                    with col3:
                        st.metric("Inconclusive", voting_result['votes']['Inconclusive'])
                    with col4:
                        st.metric("Confidence", voting_result['confidence'])

                    st.markdown("---")
                    verdict = voting_result['final_verdict']
                    confidence = voting_result['confidence']

                    if 'RhD+' in verdict:
                        st.success(f"### ✅ Final Result: {verdict}")
                    elif 'RhD-' in verdict:
                        st.warning(f"### ⚠️ Final Result: {verdict}")
                    else:
                        st.info(f"### ℹ️ Final Result: {verdict}")

                    st.write(f"**Confidence Level:** {confidence}")
                    st.write(f"**Reasoning:** {voting_result['details']}")

                    with st.expander("📋 Detailed Amplicon Analysis"):
                        for result in voting_result['amplicon_results']:
                            with st.expander(f"{result['name']} - {result['region']} ({result['length']}bp)"):
                                st.write(f"**Region:** {result['region']}")
                                st.write(f"**Length:** {result['length']} bp")
                                st.write(f"**Identity:** {result['identity']:.1f}%")
                                st.write(f"**Variants:** {result['variants']}")
                                st.write(f"**ISBT Allele:** {result.get('allele') or '-'}")
                                st.write(f"**Phenotype:** {result.get('rhd_status', '-')}")
                                st.write(f"**Vote:** {result['vote']}")
                                if result.get('zygosity'):
                                    st.write(f"**Zygosity (diagnostic SNP):** {result['zygosity'].upper()}")
                                st.write(f"**Reason:** {result['reason']}")
                                if result.get('mechanism'):
                                    st.write(f"**Mechanism:** {result['mechanism']}")
                                if result.get('serology'):
                                    st.warning(f"**Serology recommendation:** {result['serology']}")
                                if result.get('note'):
                                    st.info(f"**Note:** {result['note']}")
                                qc = result.get('qc')
                                if qc:
                                    st.write(
                                        f"**QC (Phred Q≥{qc.get('q_threshold')}):** "
                                        f"trimmed {qc.get('trimmed_5p', 0)} bp from 5', "
                                        f"{qc.get('trimmed_3p', 0)} bp from 3'; "
                                        f"masked {qc.get('masked_internal', 0)} internal low-Q bases as N; "
                                        f"encoded {qc.get('het_positions_encoded', 0)} heterozygous position(s) as IUPAC; "
                                        f"final length {qc.get('final_length', 0)} bp "
                                        f"(original {qc.get('original_length', 0)} bp)."
                                    )
                                snps = result.get('diagnostic_snps') or {}
                                if snps:
                                    st.write("**Diagnostic SNPs detected:**")
                                    snp_rows = []
                                    for snp_name, info in snps.items():
                                        snp_rows.append({
                                            'SNP': snp_name,
                                            'cDNA pos': info.get('cDNA_position'),
                                            'Ref': info.get('reference_base'),
                                            'Query': info.get('query_base'),
                                            'Zygosity': (info.get('zygosity') or '-').upper(),
                                            'Significance': info.get('significance'),
                                        })
                                    st.dataframe(pd.DataFrame(snp_rows), hide_index=True, use_container_width=True)
                else:
                    st.info("AB1 files processed but no valid sequences found for RHD analysis.")

            elif selected_system == "RHCE":
                st.subheader("RHCE Analysis Results (C/c + E/e)")

                if not rhce_reads_for_summary:
                    st.info("RHCE files processed but no valid sequences found.")
                elif rhce_error_message is not None:
                    st.error("RHCE reference data is missing.")
                    st.code(rhce_error_message)
                    st.info("See **utils/data/RHCE_REFERENCE_README.md** for download instructions.")
                else:
                    if rhce_result is not None:
                        # Headline result
                        conf = rhce_result['overall_confidence']
                        pheno = rhce_result['phenotype']
                        if conf == 'HIGH':
                            st.success(f"### ✅ Phenotype: **{pheno}**  (confidence: {conf})")
                        elif conf == 'MEDIUM':
                            st.info(f"### Phenotype: **{pheno}**  (confidence: {conf})")
                        elif conf == 'LOW':
                            st.warning(f"### ⚠️ Phenotype: **{pheno}**  (confidence: {conf})")
                        else:
                            st.error(f"### Phenotype: **{pheno}**  (confidence: {conf})")

                        st.caption(rhce_result['reason'])

                        # Read coverage metrics
                        col1, col2, col3, col4 = st.columns(4)
                        with col1:
                            st.metric("Reads (callable / total)",
                                      f"{rhce_result['reads_callable']} / {rhce_result['reads_total']}")
                        with col2:
                            st.metric("C/c call", rhce_result['c_e_call']['genotype'] or "–")
                        with col3:
                            st.metric("E/e call", rhce_result['big_E_call']['genotype'] or "–")
                        with col4:
                            st.metric("Overall confidence", conf)

                        # ISBT haplotype interpretation
                        if rhce_result['allele_options']:
                            st.markdown("#### Possible ISBT haplotype pairs")
                            st.caption(
                                "Sanger genotyping cannot phase haplotypes. "
                                "For compound heterozygotes both phase options are listed."
                            )
                            for opt in rhce_result['allele_options']:
                                st.write(f"- **{opt['haplotypes']}** → `{opt['isbt']}`")

                        # Per-SNP consensus table
                        st.markdown("#### Per-SNP consensus")
                        snp_rows = []
                        for snp_name, cons in rhce_result['snp_consensus'].items():
                            votes_str = ", ".join(f"{z}:{n}" for z, n in cons['votes'].items()) or "–"
                            snp_rows.append({
                                'SNP': snp_name,
                                'Consensus': cons['consensus'],
                                'Call': cons.get('call', '–'),
                                'Confidence': cons['confidence'],
                                'Reads covering': cons['reads_covering'],
                                'Votes': votes_str,
                                'Discordant reads': ', '.join(cons['discordant_reads']) or '–',
                            })
                        st.dataframe(pd.DataFrame(snp_rows),
                                     hide_index=True, use_container_width=True)

                        # Partial-E flags
                        partial = rhce_result['big_E_call'].get('partial_markers') or []
                        if partial:
                            st.warning(
                                "**Asian partial-E variant(s) detected:** "
                                + "; ".join(f"{p['snp']} ({p['zygosity']}) - {p['significance']}" for p in partial)
                            )

                        # Per-read drill-down
                        with st.expander("📋 Per-read details"):
                            for read in rhce_result['per_read_details']:
                                st.markdown(f"**{read['read_id']}** — "
                                            f"strand={read['strand']}, "
                                            f"identity={read['identity']}%, "
                                            f"length={read['query_length']} bp, "
                                            f"callable={read['callable']}")
                                st.caption(read['reason'])
                                read_rows = []
                                for snp_name, call in read['snp_calls'].items():
                                    if call.get('covered'):
                                        read_rows.append({
                                            'SNP': snp_name,
                                            'Call': call.get('call'),
                                            'Zygosity': call.get('zygosity', '–'),
                                            'Query base': call.get('query_base', '–'),
                                        })
                                    else:
                                        read_rows.append({
                                            'SNP': snp_name,
                                            'Call': 'not covered',
                                            'Zygosity': '–',
                                            'Query base': '–',
                                        })
                                st.dataframe(pd.DataFrame(read_rows),
                                             hide_index=True, use_container_width=True)

            elif selected_system == "Kell":
                # Mirrors the RHCE detail layout exactly: phenotype banner →
                # read-coverage metrics → ISBT haplotype → per-SNP consensus
                # table → per-read drill-down.
                st.subheader("Kell Analysis Results (K/k)")

                if not kel_reads_for_summary:
                    st.info("Kell files processed but no valid sequences found.")
                elif kel_error_message is not None:
                    st.error("KEL reference data is missing.")
                    st.code(kel_error_message)
                    st.info("See **utils/data/KEL_REFERENCE_README.md** for download instructions.")
                else:
                    if kel_result is not None:
                        conf = kel_result['overall_confidence']
                        pheno = kel_result['phenotype']
                        if conf == 'HIGH':
                            st.success(f"### ✅ Phenotype: **{pheno}**  (confidence: {conf})")
                        elif conf == 'MEDIUM':
                            st.info(f"### Phenotype: **{pheno}**  (confidence: {conf})")
                        elif conf == 'LOW':
                            st.warning(f"### ⚠️ Phenotype: **{pheno}**  (confidence: {conf})")
                        else:
                            st.error(f"### Phenotype: **{pheno}**  (confidence: {conf})")

                        st.caption(kel_result['reason'])

                        col1, col2, col3 = st.columns(3)
                        with col1:
                            st.metric("Reads (callable / total)",
                                      f"{kel_result['reads_callable']} / {kel_result['reads_total']}")
                        with col2:
                            st.metric("K/k call",
                                      (kel_result.get('k_axis_call') or {}).get('genotype') or "–")
                        with col3:
                            st.metric("Overall confidence", conf)

                        if kel_result['allele_options']:
                            st.markdown("#### ISBT haplotype")
                            st.caption(
                                "Single-axis (K/k) call — exactly one haplotype pair, "
                                "no phase ambiguity."
                            )
                            for opt in kel_result['allele_options']:
                                serology = opt.get('serology', '')
                                serology_tag = f" — {serology}" if serology else ''
                                st.write(f"- **{opt['haplotypes']}** → `{opt['isbt']}`{serology_tag}")

                        st.markdown("#### Per-SNP consensus")
                        snp_rows = []
                        for snp_name, cons in kel_result['snp_consensus'].items():
                            votes_str = ", ".join(f"{z}:{n}" for z, n in cons['votes'].items()) or "–"
                            snp_rows.append({
                                'SNP': snp_name,
                                'Consensus': cons['consensus'],
                                'Call': cons.get('call', '–'),
                                'Confidence': cons['confidence'],
                                'Reads covering': cons['reads_covering'],
                                'Votes': votes_str,
                                'Discordant reads': ', '.join(cons['discordant_reads']) or '–',
                            })
                        st.dataframe(pd.DataFrame(snp_rows),
                                     hide_index=True, use_container_width=True)

                        with st.expander("📋 Per-read details"):
                            for read in kel_result['per_read_details']:
                                st.markdown(f"**{read['read_id']}** — "
                                            f"strand={read['strand']}, "
                                            f"identity={read['identity']}%, "
                                            f"length={read['query_length']} bp, "
                                            f"callable={read['callable']}")
                                st.caption(read['reason'])
                                read_rows = []
                                for snp_name, call in read['snp_calls'].items():
                                    if call.get('covered'):
                                        read_rows.append({
                                            'SNP': snp_name,
                                            'Call': call.get('call'),
                                            'Zygosity': call.get('zygosity', '–'),
                                            'Query base': call.get('query_base', '–'),
                                        })
                                    else:
                                        read_rows.append({
                                            'SNP': snp_name,
                                            'Call': 'not covered',
                                            'Zygosity': '–',
                                            'Query base': '–',
                                        })
                                st.dataframe(pd.DataFrame(read_rows),
                                             hide_index=True, use_container_width=True)

            elif selected_system == "Duffy":
                # Mirrors the RHCE/KEL detail layout: phenotype banner ->
                # coverage metrics -> ISBT haplotype(s) -> per-SNP consensus
                # table -> per-read drill-down expander. Adds a phase-ambiguity
                # warning when GATA + A/B both heterozygous.
                st.subheader("Duffy Analysis Results (FY)")

                if not fy_reads_for_summary:
                    st.info("Duffy files processed but no valid sequences found.")
                elif fy_error_message is not None:
                    st.error("FY reference data is missing.")
                    st.code(fy_error_message)
                    st.info("See **utils/data/FY_REFERENCE_README.md** for download instructions.")
                else:
                    if fy_result is not None:
                        conf = fy_result['overall_confidence']
                        pheno = fy_result['phenotype']
                        if conf == 'HIGH':
                            st.success(f"### ✅ Phenotype: **{pheno}**  (confidence: {conf})")
                        elif conf == 'MEDIUM':
                            st.info(f"### Phenotype: **{pheno}**  (confidence: {conf})")
                        elif conf == 'LOW':
                            st.warning(f"### ⚠️ Phenotype: **{pheno}**  (confidence: {conf})")
                        else:
                            st.error(f"### Phenotype: **{pheno}**  (confidence: {conf})")

                        st.caption(fy_result['reason'])

                        col1, col2, col3, col4, col5 = st.columns(5)
                        with col1:
                            st.metric("Reads (callable / total)",
                                      f"{fy_result['reads_callable']} / {fy_result['reads_total']}")
                        with col2:
                            st.metric("FY*A/FY*B",
                                      (fy_result.get('a_b_call') or {}).get('genotype') or "–")
                        with col3:
                            gata = (fy_result.get('gata_call') or {}).get('consensus')
                            st.metric("GATA -67",
                                      {'hom_ref': 'T/T (ref)', 'het': 'T/C',
                                       'hom_alt': 'C/C (silenced)'}.get(gata, "–"))
                        with col4:
                            fy_x = (fy_result.get('fy_x_call') or {}).get('consensus')
                            st.metric("FY*X (c.265)",
                                      {'hom_ref': 'C/C (ref)', 'het': 'C/T',
                                       'hom_alt': 'T/T (weak)'}.get(fy_x, "–"))
                        with col5:
                            st.metric("Overall confidence", conf)

                        # GATA / cDNA-reference warning
                        if (fy_result.get('gata_call') or {}).get('confidence') == 'NONE':
                            st.warning(
                                "GATA -67T>C is **not callable** in this run. "
                                "If a genomic amplicon covering the ACKR1 "
                                "promoter is available, use a GenBank reference "
                                "(see FY_REFERENCE_README.md) to detect "
                                "Fy(a-b-) erythroid-silencing."
                            )

                        if fy_result['allele_options']:
                            phase_ambiguous = len(fy_result['allele_options']) > 1
                            if phase_ambiguous:
                                st.markdown("#### Possible ISBT haplotype pairs (phase ambiguous)")
                                st.caption(
                                    "GATA-null heterozygous + FY*A/FY*B heterozygous "
                                    "→ Sanger genotyping cannot phase which allele the "
                                    "GATA mutation is in cis with. Both options shown."
                                )
                            else:
                                st.markdown("#### ISBT haplotype")
                            for opt in fy_result['allele_options']:
                                serology = opt.get('serology', '')
                                serology_tag = f" — {serology}" if serology else ''
                                st.write(f"- **{opt['haplotypes']}** → `{opt['isbt']}`{serology_tag}")

                        st.markdown("#### Per-SNP consensus")
                        snp_rows = []
                        for snp_name, cons in fy_result['snp_consensus'].items():
                            votes_str = ", ".join(f"{z}:{n}" for z, n in cons['votes'].items()) or "–"
                            snp_rows.append({
                                'SNP': snp_name,
                                'Consensus': cons['consensus'],
                                'Call': cons.get('call', '–'),
                                'Confidence': cons['confidence'],
                                'Reads covering': cons['reads_covering'],
                                'Votes': votes_str,
                                'Discordant reads': ', '.join(cons['discordant_reads']) or '–',
                            })
                        st.dataframe(pd.DataFrame(snp_rows),
                                     hide_index=True, use_container_width=True)

                        with st.expander("📋 Per-read details"):
                            for read in fy_result['per_read_details']:
                                st.markdown(f"**{read['read_id']}** — "
                                            f"strand={read['strand']}, "
                                            f"identity={read['identity']}%, "
                                            f"length={read['query_length']} bp, "
                                            f"callable={read['callable']}")
                                st.caption(read['reason'])
                                read_rows = []
                                for snp_name, call in read['snp_calls'].items():
                                    if call.get('covered'):
                                        read_rows.append({
                                            'SNP': snp_name,
                                            'Call': call.get('call'),
                                            'Zygosity': call.get('zygosity', '–'),
                                            'Query base': call.get('query_base', '–'),
                                        })
                                    else:
                                        read_rows.append({
                                            'SNP': snp_name,
                                            'Call': 'not covered',
                                            'Zygosity': '–',
                                            'Query base': '–',
                                        })
                                st.dataframe(pd.DataFrame(read_rows),
                                             hide_index=True, use_container_width=True)

            elif selected_system == "Kidd":
                # Mirrors the RHCE/KEL/Duffy detail layout: phenotype banner ->
                # coverage metrics -> ISBT haplotype(s) -> per-SNP consensus ->
                # per-read drill-down. Includes phase-ambiguity warning and a
                # splice-uncallable warning when the c.342-1 SNP wasn't reached.
                st.subheader("Kidd Analysis Results (JK)")

                if not jk_reads_for_summary:
                    st.info("Kidd files processed but no valid sequences found.")
                elif jk_error_message is not None:
                    st.error("JK reference data is missing.")
                    st.code(jk_error_message)
                    st.info("See **utils/data/JK_REFERENCE_README.md** for download instructions.")
                else:
                    if jk_result is not None:
                        conf = jk_result['overall_confidence']
                        pheno = jk_result['phenotype']
                        if conf == 'HIGH':
                            st.success(f"### ✅ Phenotype: **{pheno}**  (confidence: {conf})")
                        elif conf == 'MEDIUM':
                            st.info(f"### Phenotype: **{pheno}**  (confidence: {conf})")
                        elif conf == 'LOW':
                            st.warning(f"### ⚠️ Phenotype: **{pheno}**  (confidence: {conf})")
                        else:
                            st.error(f"### Phenotype: **{pheno}**  (confidence: {conf})")

                        st.caption(jk_result['reason'])

                        col1, col2, col3, col4, col5 = st.columns(5)
                        with col1:
                            st.metric("Reads (callable / total)",
                                      f"{jk_result['reads_callable']} / {jk_result['reads_total']}")
                        with col2:
                            st.metric("JK*A/JK*B",
                                      (jk_result.get('a_b_call') or {}).get('genotype') or "–")
                        with col3:
                            splice = (jk_result.get('snp_consensus') or {}).get('c.342-1G>A', {}).get('consensus')
                            st.metric("Splice (c.342-1)",
                                      {'hom_ref': 'G/G (ref)', 'het': 'G/A',
                                       'hom_alt': 'A/A (silenced)'}.get(splice, "–"))
                        with col4:
                            asian = (jk_result.get('snp_consensus') or {}).get('c.871T>C', {}).get('consensus')
                            st.metric("Asian null (c.871)",
                                      {'hom_ref': 'T/T (ref)', 'het': 'T/C',
                                       'hom_alt': 'C/C (silenced)'}.get(asian, "–"))
                        with col5:
                            st.metric("Overall confidence", conf)

                        # Splice / cDNA-reference warning
                        splice_cons = (jk_result.get('snp_consensus') or {}).get('c.342-1G>A', {})
                        if splice_cons.get('confidence') == 'NONE':
                            st.warning(
                                "Splice-acceptor SNP **c.342-1G>A is not callable** "
                                "in this run. If a genomic amplicon covering intron "
                                "5 of SLC14A1 is available, use a GenBank reference "
                                "(see JK_REFERENCE_README.md) to detect Polynesian "
                                "Jk(a-b-)."
                            )

                        # Compound-het flag
                        if (jk_result.get('null_call') or {}).get('compound_het'):
                            st.warning(
                                "**Compound heterozygous null variants detected** "
                                "(both c.342-1 and c.871 het in the same sample). "
                                "This is rare — confidence has been lowered; "
                                "recommend chromatogram review and serology."
                            )

                        if jk_result['allele_options']:
                            phase_ambiguous = len(jk_result['allele_options']) > 1
                            if phase_ambiguous:
                                st.markdown("#### Possible ISBT haplotype pairs (phase ambiguous)")
                                st.caption(
                                    "Null variant heterozygous + JK*A/JK*B heterozygous "
                                    "→ Sanger genotyping cannot phase which allele the "
                                    "null mutation is in cis with. Both options shown."
                                )
                            else:
                                st.markdown("#### ISBT haplotype")
                            for opt in jk_result['allele_options']:
                                serology = opt.get('serology', '')
                                serology_tag = f" — {serology}" if serology else ''
                                st.write(f"- **{opt['haplotypes']}** → `{opt['isbt']}`{serology_tag}")

                        st.markdown("#### Per-SNP consensus")
                        snp_rows = []
                        for snp_name, cons in jk_result['snp_consensus'].items():
                            votes_str = ", ".join(f"{z}:{n}" for z, n in cons['votes'].items()) or "–"
                            snp_rows.append({
                                'SNP': snp_name,
                                'Consensus': cons['consensus'],
                                'Call': cons.get('call', '–'),
                                'Confidence': cons['confidence'],
                                'Reads covering': cons['reads_covering'],
                                'Votes': votes_str,
                                'Discordant reads': ', '.join(cons['discordant_reads']) or '–',
                            })
                        st.dataframe(pd.DataFrame(snp_rows),
                                     hide_index=True, use_container_width=True)

                        with st.expander("📋 Per-read details"):
                            for read in jk_result['per_read_details']:
                                st.markdown(f"**{read['read_id']}** — "
                                            f"strand={read['strand']}, "
                                            f"identity={read['identity']}%, "
                                            f"length={read['query_length']} bp, "
                                            f"callable={read['callable']}")
                                st.caption(read['reason'])
                                read_rows = []
                                for snp_name, call in read['snp_calls'].items():
                                    if call.get('covered'):
                                        read_rows.append({
                                            'SNP': snp_name,
                                            'Call': call.get('call'),
                                            'Zygosity': call.get('zygosity', '–'),
                                            'Query base': call.get('query_base', '–'),
                                        })
                                    else:
                                        read_rows.append({
                                            'SNP': snp_name,
                                            'Call': 'not covered',
                                            'Zygosity': '–',
                                            'Query base': '–',
                                        })
                                st.dataframe(pd.DataFrame(read_rows),
                                             hide_index=True, use_container_width=True)

            elif selected_system == "H":
                # Mirrors the RHCE/KEL/Duffy/Kidd detail layout, but the H
                # system has no primary antigen-discriminator SNP and adds
                # a transfusion-safety alert when the Bombay (Oh) state is
                # detected (recipient rejects standard ABO blood).
                st.subheader("H Analysis Results (FUT1)")

                if not h_reads_for_summary:
                    st.info("H files processed but no valid sequences found.")
                elif h_error_message is not None:
                    st.error("H (FUT1) reference data is missing.")
                    st.code(h_error_message)
                    st.info("See **utils/data/H_REFERENCE_README.md** for download instructions.")
                else:
                    if h_result is not None:
                        conf = h_result['overall_confidence']
                        pheno = h_result['phenotype']
                        state = (h_result.get('h_state_call') or {}).get('state')

                        # Bombay gets the loud red banner regardless of
                        # confidence — it's a transfusion-safety alert.
                        if state == 'bombay':
                            st.error(f"### 🚨 **{pheno}** (confidence: {conf})")
                            st.error(
                                "**TRANSFUSION ALERT** — Bombay (Oh) recipients "
                                "have antibodies against the H antigen and **reject "
                                "all standard ABO-typed donor blood**. "
                                "Bombay-compatible units must be sourced."
                            )
                        elif conf == 'HIGH':
                            st.success(f"### ✅ Phenotype: **{pheno}**  (confidence: {conf})")
                        elif conf == 'MEDIUM':
                            st.info(f"### Phenotype: **{pheno}**  (confidence: {conf})")
                        elif conf == 'LOW':
                            st.warning(f"### ⚠️ Phenotype: **{pheno}**  (confidence: {conf})")
                        else:
                            st.error(f"### Phenotype: **{pheno}**  (confidence: {conf})")

                        st.caption(h_result['reason'])

                        col1, col2, col3, col4, col5 = st.columns(5)
                        with col1:
                            st.metric("Reads (callable / total)",
                                      f"{h_result['reads_callable']} / {h_result['reads_total']}")
                        with col2:
                            c725 = (h_result.get('snp_consensus') or {}).get('c.725T>G', {}).get('consensus')
                            st.metric("c.725 (Indian Bombay)",
                                      {'hom_ref': 'T/T (ref)', 'het': 'T/G',
                                       'hom_alt': 'G/G (Bombay)'}.get(c725, "–"))
                        with col3:
                            c586 = (h_result.get('snp_consensus') or {}).get('c.586C>T', {}).get('consensus')
                            st.metric("c.586 (nonsense)",
                                      {'hom_ref': 'C/C (ref)', 'het': 'C/T',
                                       'hom_alt': 'T/T (Bombay)'}.get(c586, "–"))
                        with col4:
                            c460 = (h_result.get('snp_consensus') or {}).get('c.460T>C', {}).get('consensus')
                            st.metric("c.460 (h2 weak)",
                                      {'hom_ref': 'T/T (ref)', 'het': 'T/C',
                                       'hom_alt': 'C/C (weak)'}.get(c460, "–"))
                        with col5:
                            st.metric("Overall confidence", conf)

                        if h_result['allele_options']:
                            st.markdown("#### ISBT haplotype")
                            for opt in h_result['allele_options']:
                                serology = opt.get('serology', '')
                                serology_tag = f" — {serology}" if serology else ''
                                st.write(f"- **{opt['haplotypes']}** → `{opt['isbt']}`{serology_tag}")

                        st.markdown("#### Per-SNP consensus")
                        snp_rows = []
                        for snp_name, cons in h_result['snp_consensus'].items():
                            votes_str = ", ".join(f"{z}:{n}" for z, n in cons['votes'].items()) or "–"
                            snp_rows.append({
                                'SNP': snp_name,
                                'Consensus': cons['consensus'],
                                'Call': cons.get('call', '–'),
                                'Confidence': cons['confidence'],
                                'Reads covering': cons['reads_covering'],
                                'Votes': votes_str,
                                'Discordant reads': ', '.join(cons['discordant_reads']) or '–',
                            })
                        st.dataframe(pd.DataFrame(snp_rows),
                                     hide_index=True, use_container_width=True)

                        with st.expander("📋 Per-read details"):
                            for read in h_result['per_read_details']:
                                st.markdown(f"**{read['read_id']}** — "
                                            f"strand={read['strand']}, "
                                            f"identity={read['identity']}%, "
                                            f"length={read['query_length']} bp, "
                                            f"callable={read['callable']}")
                                st.caption(read['reason'])
                                read_rows = []
                                for snp_name, call in read['snp_calls'].items():
                                    if call.get('covered'):
                                        read_rows.append({
                                            'SNP': snp_name,
                                            'Call': call.get('call'),
                                            'Zygosity': call.get('zygosity', '–'),
                                            'Query base': call.get('query_base', '–'),
                                        })
                                    else:
                                        read_rows.append({
                                            'SNP': snp_name,
                                            'Call': 'not covered',
                                            'Zygosity': '–',
                                            'Query base': '–',
                                        })
                                st.dataframe(pd.DataFrame(read_rows),
                                             hide_index=True, use_container_width=True)

            elif selected_system == "MNS":
                # Mirrors the established detail layout: phenotype banner ->
                # coverage metrics -> ISBT haplotype(s) -> per-SNP consensus
                # -> per-read drill-down. Adds per-gene identity reporting
                # and a phase-ambiguity warning for the MN+Ss double-het case.
                st.subheader("MNS Analysis Results (GYPA + GYPB)")

                if not mns_reads_for_summary:
                    st.info("MNS files processed but no valid sequences found.")
                elif mns_error_message is not None:
                    st.error("MNS reference data is missing.")
                    st.code(mns_error_message)
                    st.info("See **utils/data/MNS_REFERENCE_README.md** for download instructions.")
                else:
                    if mns_result is not None:
                        conf = mns_result['overall_confidence']
                        pheno = mns_result['phenotype']
                        if conf == 'HIGH':
                            st.success(f"### ✅ Phenotype: **{pheno}**  (confidence: {conf})")
                        elif conf == 'MEDIUM':
                            st.info(f"### Phenotype: **{pheno}**  (confidence: {conf})")
                        elif conf == 'LOW':
                            st.warning(f"### ⚠️ Phenotype: **{pheno}**  (confidence: {conf})")
                        else:
                            st.error(f"### Phenotype: **{pheno}**  (confidence: {conf})")

                        st.caption(mns_result['reason'])

                        col1, col2, col3, col4 = st.columns(4)
                        with col1:
                            st.metric("Reads (callable / total)",
                                      f"{mns_result['reads_callable']} / {mns_result['reads_total']}")
                        with col2:
                            st.metric("M/N (GYPA)",
                                      (mns_result.get('m_n_call') or {}).get('genotype') or "–")
                        with col3:
                            st.metric("S/s (GYPB)",
                                      (mns_result.get('s_s_call') or {}).get('genotype') or "–")
                        with col4:
                            st.metric("Overall confidence", conf)

                        if mns_result['allele_options']:
                            phase_ambiguous = len(mns_result['allele_options']) > 1
                            if phase_ambiguous:
                                st.markdown("#### Possible ISBT haplotype pairs (phase ambiguous)")
                                st.caption(
                                    "Both M/N and S/s heterozygous → the four "
                                    "MNS haplotypes (MS, Ms, NS, Ns) can pair as "
                                    "MS/Ns or Ms/NS. Sanger genotyping cannot "
                                    "phase the two — both options shown."
                                )
                            else:
                                st.markdown("#### ISBT haplotype")
                            for opt in mns_result['allele_options']:
                                serology = opt.get('serology', '')
                                serology_tag = f" — {serology}" if serology else ''
                                st.write(f"- **{opt['haplotypes']}** → `{opt['isbt']}`{serology_tag}")

                        st.markdown("#### Per-SNP consensus")
                        snp_rows = []
                        for snp_name, cons in mns_result['snp_consensus'].items():
                            votes_str = ", ".join(f"{z}:{n}" for z, n in cons['votes'].items()) or "–"
                            snp_rows.append({
                                'SNP': snp_name,
                                'Consensus': cons['consensus'],
                                'Call': cons.get('call', '–'),
                                'Confidence': cons['confidence'],
                                'Reads covering': cons['reads_covering'],
                                'Votes': votes_str,
                                'Discordant reads': ', '.join(cons['discordant_reads']) or '–',
                            })
                        st.dataframe(pd.DataFrame(snp_rows),
                                     hide_index=True, use_container_width=True)

                        with st.expander("📋 Per-read details"):
                            for read in mns_result['per_read_details']:
                                per_gene = read.get('per_gene_identity') or {}
                                per_gene_str = (
                                    "  ".join(f"{g}={i}%" for g, i in per_gene.items())
                                    if per_gene else "—"
                                )
                                st.markdown(
                                    f"**{read['read_id']}** — "
                                    f"best={read.get('best_gene', '?')} "
                                    f"strand={read['strand']}, "
                                    f"identity={read['identity']}%, "
                                    f"length={read['query_length']} bp, "
                                    f"callable={read['callable']}"
                                )
                                st.caption(f"Per-gene identity: {per_gene_str}")
                                st.caption(read['reason'])
                                read_rows = []
                                for snp_name, call in read['snp_calls'].items():
                                    if call.get('covered'):
                                        read_rows.append({
                                            'SNP': snp_name,
                                            'Gene': call.get('gene', '–'),
                                            'Call': call.get('call'),
                                            'Zygosity': call.get('zygosity', '–'),
                                            'Query base': call.get('query_base', '–'),
                                        })
                                    else:
                                        read_rows.append({
                                            'SNP': snp_name,
                                            'Gene': '–',
                                            'Call': 'not covered',
                                            'Zygosity': '–',
                                            'Query base': '–',
                                        })
                                st.dataframe(pd.DataFrame(read_rows),
                                             hide_index=True, use_container_width=True)

            elif selected_system == "Diego":
                # Mirrors the established detail layout. Single-axis system
                # (Di(a)/Di(b)) so the metric strip is compact and there's
                # no phase-ambiguity branch.
                st.subheader("Diego Analysis Results (DI / SLC4A1)")

                if not di_reads_for_summary:
                    st.info("Diego files processed but no valid sequences found.")
                elif di_error_message is not None:
                    st.error("Diego (SLC4A1) reference data is missing.")
                    st.code(di_error_message)
                    st.info("See **utils/data/DI_REFERENCE_README.md** for download instructions.")
                else:
                    if di_result is not None:
                        conf = di_result['overall_confidence']
                        pheno = di_result['phenotype']
                        if conf == 'HIGH':
                            st.success(f"### ✅ Phenotype: **{pheno}**  (confidence: {conf})")
                        elif conf == 'MEDIUM':
                            st.info(f"### Phenotype: **{pheno}**  (confidence: {conf})")
                        elif conf == 'LOW':
                            st.warning(f"### ⚠️ Phenotype: **{pheno}**  (confidence: {conf})")
                        else:
                            st.error(f"### Phenotype: **{pheno}**  (confidence: {conf})")

                        st.caption(di_result['reason'])

                        col1, col2, col3 = st.columns(3)
                        with col1:
                            st.metric("Reads (callable / total)",
                                      f"{di_result['reads_callable']} / {di_result['reads_total']}")
                        with col2:
                            st.metric("Di(a)/Di(b) call",
                                      (di_result.get('di_axis_call') or {}).get('genotype') or "–")
                        with col3:
                            st.metric("Overall confidence", conf)

                        if di_result['allele_options']:
                            st.markdown("#### ISBT haplotype")
                            for opt in di_result['allele_options']:
                                serology = opt.get('serology', '')
                                serology_tag = f" — {serology}" if serology else ''
                                st.write(f"- **{opt['haplotypes']}** → `{opt['isbt']}`{serology_tag}")

                        st.markdown("#### Per-SNP consensus")
                        snp_rows = []
                        for snp_name, cons in di_result['snp_consensus'].items():
                            votes_str = ", ".join(f"{z}:{n}" for z, n in cons['votes'].items()) or "–"
                            snp_rows.append({
                                'SNP': snp_name,
                                'Consensus': cons['consensus'],
                                'Call': cons.get('call', '–'),
                                'Confidence': cons['confidence'],
                                'Reads covering': cons['reads_covering'],
                                'Votes': votes_str,
                                'Discordant reads': ', '.join(cons['discordant_reads']) or '–',
                            })
                        st.dataframe(pd.DataFrame(snp_rows),
                                     hide_index=True, use_container_width=True)

                        with st.expander("📋 Per-read details"):
                            for read in di_result['per_read_details']:
                                st.markdown(f"**{read['read_id']}** — "
                                            f"strand={read['strand']}, "
                                            f"identity={read['identity']}%, "
                                            f"length={read['query_length']} bp, "
                                            f"callable={read['callable']}")
                                st.caption(read['reason'])
                                read_rows = []
                                for snp_name, call in read['snp_calls'].items():
                                    if call.get('covered'):
                                        read_rows.append({
                                            'SNP': snp_name,
                                            'Call': call.get('call'),
                                            'Zygosity': call.get('zygosity', '–'),
                                            'Query base': call.get('query_base', '–'),
                                        })
                                    else:
                                        read_rows.append({
                                            'SNP': snp_name,
                                            'Call': 'not covered',
                                            'Zygosity': '–',
                                            'Query base': '–',
                                        })
                                st.dataframe(pd.DataFrame(read_rows),
                                             hide_index=True, use_container_width=True)


else:
    st.info("Upload files and click **Analyze** to start.")
