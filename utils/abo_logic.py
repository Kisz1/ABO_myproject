"""ABO allele identification logic, decoupled from Streamlit UI.

Extracted from main.py for unit-testability.
"""

import utils.abo_identifier as abo_utils


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
