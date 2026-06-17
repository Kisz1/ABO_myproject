"""RHD multi-amplicon voting/consolidation logic, decoupled from Streamlit UI.

Extracted from main.py for unit-testability.
"""


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
            'reference_length':         {'RHD1': 951, 'RHD456': 3314, 'RHD7': 934, 'RHD9': 874}.get(region, 0),
            'matched_region':           region,
            'reference_description':    f"{region}_ref (WHO standard)",
            'matched_identity':         round(identity, 1),
            'identity':                 round(identity, 1),
            'matched_score':            0.0,
            'variant_count':            len(variants),
            'variants':                 variants,
            'score_1':                  round(identity, 1) if region == "RHD1"   else 0.0,
            'score_456':                round(identity, 1) if region == "RHD456" else 0.0,
            'score_7':                  round(identity, 1) if region == "RHD7"   else 0.0,
            'score_9':                  round(identity, 1) if region == "RHD9"   else 0.0,
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
