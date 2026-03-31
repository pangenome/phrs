#!/usr/bin/env python3

import json
import csv

# Load g:Profiler results for filtered analysis (excluding acrocentric PHRs)
with open('gprofiler_results_no_acro.json', 'r') as f:
    data = json.load(f)

print("=== g:Profiler GO Enrichment Analysis Results (No Acrocentric PHRs) ===")

if 'result' not in data:
    print("No results found in the response")
    print("Full response:", data)
    exit(1)

results = data['result']
meta = data['meta']

query_genes = meta['query_metadata']['queries']['query_1']
print(f"Query info: {len(query_genes)} genes analyzed (excluding acrocentric PHRs)")
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

print(f"Found {len(go_bp_results)} significant GO Biological Process terms")
print(f"Found {len(go_mf_results)} significant GO Molecular Function terms")
print(f"Found {len(kegg_results)} significant KEGG pathways")
print()

# Save GO:BP results to CSV
if go_bp_results:
    print("Top 10 GO Biological Process enrichments:")
    with open('phr_no_acro_GO_BP_enrichment.csv', 'w', newline='') as csvfile:
        fieldnames = ['term_id', 'term_name', 'p_value', 'adjusted_p_value', 'term_size', 'query_size', 'intersection_size', 'effective_domain_size', 'source']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for i, result in enumerate(sorted(go_bp_results, key=lambda x: x['p_value'])[:10]):
            print(f"{i+1}. {result['name']} (p={result['p_value']:.2e})")
            writer.writerow({
                'term_id': result['native'],
                'term_name': result['name'],
                'p_value': result['p_value'],
                'adjusted_p_value': result['p_value'],  # g:Profiler already provides adjusted p-values
                'term_size': result['term_size'],
                'query_size': result['query_size'],
                'intersection_size': result['intersection_size'],
                'effective_domain_size': result['effective_domain_size'],
                'source': result['source']
            })

        # Write all results
        for result in go_bp_results:
            writer.writerow({
                'term_id': result['native'],
                'term_name': result['name'],
                'p_value': result['p_value'],
                'adjusted_p_value': result['p_value'],
                'term_size': result['term_size'],
                'query_size': result['query_size'],
                'intersection_size': result['intersection_size'],
                'effective_domain_size': result['effective_domain_size'],
                'source': result['source']
            })

# Save GO:MF results to CSV
if go_mf_results:
    print("\nTop 10 GO Molecular Function enrichments:")
    with open('phr_no_acro_GO_MF_enrichment.csv', 'w', newline='') as csvfile:
        fieldnames = ['term_id', 'term_name', 'p_value', 'adjusted_p_value', 'term_size', 'query_size', 'intersection_size', 'effective_domain_size', 'source']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for i, result in enumerate(sorted(go_mf_results, key=lambda x: x['p_value'])[:10]):
            print(f"{i+1}. {result['name']} (p={result['p_value']:.2e})")
            writer.writerow({
                'term_id': result['native'],
                'term_name': result['name'],
                'p_value': result['p_value'],
                'adjusted_p_value': result['p_value'],
                'term_size': result['term_size'],
                'query_size': result['query_size'],
                'intersection_size': result['intersection_size'],
                'effective_domain_size': result['effective_domain_size'],
                'source': result['source']
            })

        # Write all results
        for result in go_mf_results:
            writer.writerow({
                'term_id': result['native'],
                'term_name': result['name'],
                'p_value': result['p_value'],
                'adjusted_p_value': result['p_value'],
                'term_size': result['term_size'],
                'query_size': result['query_size'],
                'intersection_size': result['intersection_size'],
                'effective_domain_size': result['effective_domain_size'],
                'source': result['source']
            })

# Save KEGG results to CSV
if kegg_results:
    print(f"\nTop 10 KEGG pathway enrichments:")
    with open('phr_no_acro_KEGG_enrichment.csv', 'w', newline='') as csvfile:
        fieldnames = ['term_id', 'term_name', 'p_value', 'adjusted_p_value', 'term_size', 'query_size', 'intersection_size', 'effective_domain_size', 'source']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for i, result in enumerate(sorted(kegg_results, key=lambda x: x['p_value'])[:10]):
            print(f"{i+1}. {result['name']} (p={result['p_value']:.2e})")
            writer.writerow({
                'term_id': result['native'],
                'term_name': result['name'],
                'p_value': result['p_value'],
                'adjusted_p_value': result['p_value'],
                'term_size': result['term_size'],
                'query_size': result['query_size'],
                'intersection_size': result['intersection_size'],
                'effective_domain_size': result['effective_domain_size'],
                'source': result['source']
            })

        # Write all results
        for result in kegg_results:
            writer.writerow({
                'term_id': result['native'],
                'term_name': result['name'],
                'p_value': result['p_value'],
                'adjusted_p_value': result['p_value'],
                'term_size': result['term_size'],
                'query_size': result['query_size'],
                'intersection_size': result['intersection_size'],
                'effective_domain_size': result['effective_domain_size'],
                'source': result['source']
            })

print(f"\nGene mapping info:")
query_genes = meta['query_metadata']['queries']['query_1']
print(f"Total genes submitted: 220 (excluding acrocentric PHRs)")
print(f"Genes recognized: {len(query_genes)}")
print(f"Genes lost: {220 - len(query_genes)}")

print(f"\nResults saved:")
print(f"- phr_no_acro_GO_BP_enrichment.csv ({len(go_bp_results)} terms)")
print(f"- phr_no_acro_GO_MF_enrichment.csv ({len(go_mf_results)} terms)")
print(f"- phr_no_acro_KEGG_enrichment.csv ({len(kegg_results)} pathways)")