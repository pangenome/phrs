#!/usr/bin/env python3

import json
import csv

# Load g:Profiler results for protein-coding genes only
with open('gprofiler_results_coding_only.json', 'r') as f:
    data = json.load(f)

print("=== g:Profiler GO Enrichment Analysis Results (Protein-Coding Genes Only) ===")

if 'result' not in data:
    print("No results found in the response")
    print("Full response:", data)
    exit(1)

results = data['result']
meta = data['meta']

query_genes = meta['query_metadata']['queries']['query_1']
print(f"Query info: {len(query_genes)} protein-coding genes analyzed")
print(f"Genes: {', '.join(query_genes)}")
print(f"Genome background sources: {meta['query_metadata']['sources']}")
print()

# Separate results by source
go_bp_results = []
go_mf_results = []
kegg_results = []

for result in results:
    if result['source'] == 'GO:BP':
        go_bp_results.append(result)
    elif result['source'] == 'GO:MF':
        go_mf_results.append(result)
    elif result['source'] == 'KEGG':
        kegg_results.append(result)

print(f"Enrichment Results Summary:")
print(f"- GO:BP (Biological Process): {len(go_bp_results)} terms")
print(f"- GO:MF (Molecular Function): {len(go_mf_results)} terms")
print(f"- KEGG pathways: {len(kegg_results)} terms")
print()

# Save GO:BP results to CSV
if go_bp_results:
    fieldnames = ['native', 'name', 'description', 'p_value', 'intersection_size', 'query_size', 'term_size', 'precision', 'recall']

    with open('phr_coding_only_GO_BP_enrichment.csv', 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for result in go_bp_results:
            row = {key: result.get(key, '') for key in fieldnames}
            writer.writerow(row)

    print("GO:BP Enrichment Results:")
    for i, result in enumerate(go_bp_results, 1):
        print(f"  {i}. {result['name']} (p={result['p_value']:.2e})")
        print(f"     Genes: {result['intersection_size']}/{result['query_size']} (precision: {result['precision']:.1%})")
    print()
else:
    # Create empty file
    with open('phr_coding_only_GO_BP_enrichment.csv', 'w', newline='') as csvfile:
        csvfile.write("")
    print("No GO:BP enrichment found\n")

# Save GO:MF results to CSV
if go_mf_results:
    fieldnames = ['native', 'name', 'description', 'p_value', 'intersection_size', 'query_size', 'term_size', 'precision', 'recall']

    with open('phr_coding_only_GO_MF_enrichment.csv', 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for result in go_mf_results:
            row = {key: result.get(key, '') for key in fieldnames}
            writer.writerow(row)

    print("GO:MF Enrichment Results:")
    for i, result in enumerate(go_mf_results, 1):
        print(f"  {i}. {result['name']} (p={result['p_value']:.2e})")
        print(f"     Genes: {result['intersection_size']}/{result['query_size']} (precision: {result['precision']:.1%})")
    print()
else:
    # Create empty file
    with open('phr_coding_only_GO_MF_enrichment.csv', 'w', newline='') as csvfile:
        csvfile.write("")
    print("No GO:MF enrichment found\n")

# Get gene mapping info
print("Gene Mapping Information:")
if 'genes_metadata' in meta:
    mapping_data = meta['genes_metadata']['query']['query_1']['mapping'] if 'query' in meta['genes_metadata'] else {}
    failed_genes = meta['genes_metadata'].get('failed', [])
    ambiguous_genes = meta['genes_metadata'].get('ambiguous', {})

    print(f"Successfully mapped genes: {len(mapping_data)}")
    if failed_genes:
        print(f"Failed to map: {', '.join(failed_genes)}")
    if ambiguous_genes:
        print(f"Ambiguous mappings: {', '.join(ambiguous_genes.keys())}")
    print(f"Total genes in analysis: {len(mapping_data)}")