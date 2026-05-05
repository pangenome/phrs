#!/usr/bin/env python3
"""
Run GO enrichment analysis on protein-coding genes only from PHR intervals
"""

import requests
import json
import pandas as pd
import sys
import time

def load_gene_names(filename):
    """Load gene names from file, one per line"""
    with open(filename, 'r') as f:
        genes = [line.strip() for line in f if line.strip()]
    return genes

def run_gprofiler_enrichment(gene_list, organism="hsapiens"):
    """Run g:Profiler enrichment analysis"""

    # Prepare the request
    request_data = {
        "organism": organism,
        "query": gene_list,
        "sources": ["GO:BP", "GO:MF", "KEGG"],
        "user_threshold": 0.05,
        "significance_threshold_method": "fdr"
    }

    # Save request for reference
    with open('gprofiler_request_coding_only.json', 'w') as f:
        json.dump(request_data, f, indent=2)

    print(f"Running g:Profiler analysis on {len(gene_list)} protein-coding genes...")
    print(f"Gene list: {', '.join(gene_list)}")

    # Make the request
    response = requests.post('https://biit.cs.ut.ee/gprofiler/api/gost/profile/',
                           json=request_data)

    if response.status_code == 200:
        results = response.json()

        # Save full results
        with open('gprofiler_results_coding_only.json', 'w') as f:
            json.dump(results, f, indent=2)

        return results
    else:
        print(f"Error: {response.status_code}")
        print(response.text)
        return None

def process_results(results):
    """Process g:Profiler results and save to CSV"""

    if not results or 'result' not in results:
        print("No results returned from g:Profiler")
        return None, None

    result_data = results['result']

    if not result_data:
        print("No significant enrichment found!")
        # Create empty CSVs to match expected output format
        empty_df = pd.DataFrame()
        empty_df.to_csv('phr_coding_only_GO_BP_enrichment.csv', index=False)
        empty_df.to_csv('phr_coding_only_GO_MF_enrichment.csv', index=False)
        return None, None

    # Convert to DataFrame
    df = pd.DataFrame(result_data)

    # Separate by source
    go_bp_results = df[df['source'] == 'GO:BP'].copy() if 'source' in df.columns else pd.DataFrame()
    go_mf_results = df[df['source'] == 'GO:MF'].copy() if 'source' in df.columns else pd.DataFrame()

    # Save results
    go_bp_results.to_csv('phr_coding_only_GO_BP_enrichment.csv', index=False)
    go_mf_results.to_csv('phr_coding_only_GO_MF_enrichment.csv', index=False)

    # Print summary
    print(f"\nEnrichment Results Summary:")
    print(f"GO:BP terms found: {len(go_bp_results)}")
    print(f"GO:MF terms found: {len(go_mf_results)}")

    if len(go_bp_results) > 0:
        print(f"\nTop GO:BP terms:")
        for _, row in go_bp_results.head().iterrows():
            print(f"  {row['term_name']} (p={row['p_value']:.2e})")

    if len(go_mf_results) > 0:
        print(f"\nTop GO:MF terms:")
        for _, row in go_mf_results.head().iterrows():
            print(f"  {row['term_name']} (p={row['p_value']:.2e})")

    return go_bp_results, go_mf_results

def main():
    # Load protein-coding gene names
    gene_names = load_gene_names('phrs.no_acro.coding_gene_names.txt')

    print(f"Loaded {len(gene_names)} protein-coding gene names from PHR intervals")
    print(f"Genes: {', '.join(gene_names)}")

    # Run enrichment
    results = run_gprofiler_enrichment(gene_names)

    if results:
        go_bp, go_mf = process_results(results)
        print("\nEnrichment analysis complete!")
        print("Results saved to:")
        print("- phr_coding_only_GO_BP_enrichment.csv")
        print("- phr_coding_only_GO_MF_enrichment.csv")
        print("- gprofiler_request_coding_only.json")
        print("- gprofiler_results_coding_only.json")
    else:
        print("Enrichment analysis failed!")
        return 1

    return 0

if __name__ == "__main__":
    sys.exit(main())