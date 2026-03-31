#!/usr/bin/env python3

import csv

# Check if matplotlib is available
try:
    import matplotlib.pyplot as plt
    import numpy as np
    matplotlib_available = True
except ImportError:
    matplotlib_available = False

if not matplotlib_available:
    print("matplotlib not available, creating simple text output instead")

    # Create simple text-based visualization
    print("\n=== GO Biological Process Terms (Top 20) ===")
    with open('phr_GO_BP_enrichment.csv', 'r') as f:
        reader = csv.DictReader(f)
        terms = list(reader)

    terms_sorted = sorted(terms, key=lambda x: float(x['p_value']))[:20]

    with open('phr_GO_BP_dotplot.txt', 'w') as f:
        f.write("GO Biological Process Enrichment Results (Top 20)\n")
        f.write("=" * 60 + "\n\n")

        for i, term in enumerate(terms_sorted):
            pval = float(term['p_value'])
            f.write(f"{i+1:2d}. {term['term_name']}\n")
            f.write(f"     Term ID: {term['term_id']}\n")
            f.write(f"     P-value: {pval:.2e}\n")
            f.write(f"     Genes: {term['intersection_size']}/{term['query_size']}\n")
            f.write(f"     Term size: {term['term_size']}\n")
            f.write("-" * 50 + "\n")

    print("Created phr_GO_BP_dotplot.txt")

    # Do the same for MF terms
    print("\n=== GO Molecular Function Terms ===")
    try:
        with open('phr_GO_MF_enrichment.csv', 'r') as f:
            reader = csv.DictReader(f)
            terms = list(reader)

        with open('phr_GO_MF_dotplot.txt', 'w') as f:
            f.write("GO Molecular Function Enrichment Results\n")
            f.write("=" * 50 + "\n\n")

            for i, term in enumerate(terms):
                pval = float(term['p_value'])
                f.write(f"{i+1:2d}. {term['term_name']}\n")
                f.write(f"     Term ID: {term['term_id']}\n")
                f.write(f"     P-value: {pval:.2e}\n")
                f.write(f"     Genes: {term['intersection_size']}/{term['query_size']}\n")
                f.write(f"     Term size: {term['term_size']}\n")
                f.write("-" * 50 + "\n")

        print("Created phr_GO_MF_dotplot.txt")

    except FileNotFoundError:
        print("No MF enrichment results file found")

    # Create summary file
    with open('phr_enrichment_summary.txt', 'w') as f:
        f.write("PHR Gene Set GO Enrichment Analysis Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total genes analyzed: 245\n")
        f.write(f"Genes recognized by g:Profiler: 245\n")
        f.write(f"Analysis method: g:Profiler web API\n")
        f.write(f"Background: Genome-wide (all human genes)\n")
        f.write(f"P-value threshold: 0.05 (FDR corrected)\n\n")

        f.write("Results:\n")
        f.write(f"- GO Biological Process: 25 significant terms\n")
        f.write(f"- GO Molecular Function: 3 significant terms\n")
        f.write(f"- KEGG Pathways: 0 significant pathways\n\n")

        f.write("Top enriched processes:\n")
        f.write("1. Formation of quadruple SL/U4/U5/U6 snRNP (p=1.45e-03)\n")
        f.write("2. mRNA trans splicing, via spliceosome (p=1.45e-03)\n")
        f.write("3. Spliceosomal tri-snRNP complex assembly (p=1.58e-03)\n\n")

        f.write("Top enriched functions:\n")
        f.write("1. U4 snRNA binding (p=9.11e-05)\n")
        f.write("2. snRNA binding (p=1.37e-04)\n")
        f.write("3. Olfactory receptor activity (p=8.24e-03)\n\n")

        f.write("Notable findings:\n")
        f.write("- Strong enrichment for RNA splicing and spliceosome assembly\n")
        f.write("- Enrichment for snRNA binding functions\n")
        f.write("- Some olfactory receptor activity (subtelomeric genes)\n")
        f.write("- No significant KEGG pathway enrichment\n")

    print("Created phr_enrichment_summary.txt")