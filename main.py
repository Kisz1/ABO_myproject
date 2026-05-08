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
from utils.rhd_analyzer import RHDAnalyzer, RHD1_REFERENCE, RHD456_REFERENCE
from utils.isbt_handler import ISBTDataHandler
from utils.bloodgroup import get_system, get_available_system_keys

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

    mapping_service = fasta_utils.FASTAAlignmentService()
    traces = None
    if len(fwd_ab1_files) > 1:
        for i in fwd_ab1_files:
            ab1_service = ab1_utils.AB1Analyzer()
            trace = ab1_service.read_ab1_trace(i)
            trace = ab1_service.merge_overlap(
                traces, trace) if traces else trace
            traces = trace

    else:
        ab1_service = ab1_utils.AB1Analyzer()
        traces = ab1_service.read_ab1_trace(fwd_ab1_files[0])

    if traces:
        merged_reverse = ab1_service.reverse_chromatogram(traces)
        merged_norm = ab1_service.normalize_trace_per_channel(merged_reverse)
        raw_hets = ab1_service.detect_hetero(merged_reverse, ratio=threshold_ratio)
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

        # For AB1 files: Return the full traces as-is for RHD analysis
        # (AB1 files are typically RHD sequences, not ABO exon-specific)
        # We'll extract the sequence for RHD variant analysis later
        results = [merged_reverse]  # Wrap in list to maintain consistency
        return results, hets

    return None, None


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
            'reference_length':         len(RHD1_REFERENCE) if region == "RHD1" else len(RHD456_REFERENCE),
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
fwd_ab1 = st.sidebar.file_uploader("Upload  AB1 file", type=[
    "ab1"], accept_multiple_files=True, help="You can upload multiple files for batch processing.")

fasta_files = st.sidebar.file_uploader(
    "Upload exon-specific FASTA", type=["fasta", "fa", "fas"], accept_multiple_files=True)
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


def generate_final_blood_group_summary(robust_summary, processed_AB1, hets, isbt_handler, has_fasta, rhd_reference_seq=None):
    """
    Generate a comprehensive final blood group summary from all analyzed systems.

    Args:
        robust_summary: ABO analysis results from FASTA files
        processed_AB1: AB1 processing results (for RHD)
        hets: Heterozygote data
        isbt_handler: ISBT handler for phenotype suggestions
        has_fasta: bool indicating whether FASTA files were uploaded

    Returns:
        Dictionary with final blood group summary
    """
    summary = {
        'abo': {'status': 'not_analyzed', 'phenotype': None, 'alleles': [], 'variants': [], 'message': None},
        'rhd': {'status': 'not_analyzed', 'phenotype': None, 'alleles': [], 'variants': [], 'identity': None, 'query_length': None, 'reference_length': None, 'note': None}
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
        col1, col2 = st.columns(2)

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
    if not fwd_ab1 and not fasta_files:
        st.warning("Please upload at least one file before analyzing.")
    else:
        status_container = st.empty()
        status_container.success("Files uploaded successfully! Starting analysis...")

        # --- Robust Processing Logic ---
        robust_summary = []
        if fasta_files:
            service = fasta_utils.FASTAAlignmentService()
            robust_summary = service.generate_batch_summary(fasta_files)

        confirmed_results = [r for r in robust_summary if "Confirmed" in r.get('decision', "")]

        # Build exons_ref from confirmed FASTA results so AB1 exon extraction has coordinates
        if confirmed_results and fasta_files:
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

        processed_AB1, hets = process_ab1_files(
            fwd_ab1, exons_ref, threshold_ratio) if fwd_ab1 else (None, None)

        status_container.empty()

        with tab1:
            st.subheader("Chromatogram Check for Heterozygotes")
            st.write("### Heterozygote Positions Detected: ",
                     len(hets) if hets else 0)
            if fasta_files and not fwd_ab1:
                st.info("FASTA results are in the **Exon-based SNP** and **Allele Prediction** tabs.")
            if processed_AB1:
                hetero_sites_plot = [(h[0], dict(h[1])) for h in hets] if hets else None
                for i in processed_AB1:  # type: ignore
                    x = i['exon']
                    cds_start, cds_end = get_cds(x)
                    fig = plot_chromatogram_plotly(
                        i, base_width=2, cds_start=cds_start, cds_end=cds_end, hetero_sites=hetero_sites_plot)
                    st.plotly_chart(fig, use_container_width=True)
            elif not fwd_ab1:
                st.info("Upload an AB1 file to view the chromatogram.")

            if hets:
                st.write("### Detected Heterozygous Positions:")
                het_data = {
                    "Position": [h[0] for h in hets],
                    "Top Base": [h[1][0][0] for h in hets],
                    "Alt Base": [h[1][1][0] if len(h[1]) > 1 else "" for h in hets],
                    "Ratio": [round(h[1][1][1] / (h[1][0][1] + 1e-6), 2) if len(h[1]) > 1 else 0 for h in hets]
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


else:
    st.info("Upload files and click **Analyze** to start.")
