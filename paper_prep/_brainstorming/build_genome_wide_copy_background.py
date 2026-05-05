#!/usr/bin/env python3

import gzip
import csv
import re
from collections import defaultdict

def parse_gff_attributes(attr_string):
    """Parse GFF3 attributes string into dictionary."""
    attrs = {}
    for item in attr_string.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            attrs[key] = value
    return attrs

def extract_gene_name(attributes):
    """Extract standardized gene name from GFF attributes."""
    # Look for Name first, then gene, then other identifiers
    if 'Name' in attributes:
        return attributes['Name']
    elif 'gene' in attributes:
        return attributes['gene']
    elif 'gene_name' in attributes:
        return attributes['gene_name']
    elif 'ID' in attributes:
        # For gene features, ID might be the gene name
        return attributes['ID']
    return None

def build_genome_wide_copy_counts():
    """Build genome-wide copy counts for all genes from GFF3 file."""
    print("Processing genome-wide GFF3 file...")

    gene_counts = defaultdict(int)
    gene_info = {}

    # Open the compressed GFF3 file
    with gzip.open('chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz', 'rt') as f:
        for line_num, line in enumerate(f):
            if line_num % 100000 == 0:
                print(f"Processed {line_num} lines...")

            # Skip comments and headers
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) != 9:
                continue

            seqname, source, feature_type, start, end, score, strand, frame, attributes = parts

            # Only process gene features
            if feature_type != 'gene':
                continue

            # Parse attributes
            attrs = parse_gff_attributes(attributes)
            gene_name = extract_gene_name(attrs)

            if gene_name:
                # Clean up gene name (remove version numbers, etc.)
                clean_name = re.sub(r'\.\d+$', '', gene_name)  # Remove .1, .2 etc.
                gene_counts[clean_name] += 1

                # Store additional info
                if clean_name not in gene_info:
                    gene_info[clean_name] = {
                        'biotype': attrs.get('gene_biotype', 'unknown'),
                        'description': attrs.get('description', ''),
                        'chromosomes': set(),
                        'arms': set()
                    }

                # Add chromosome and arm info
                gene_info[clean_name]['chromosomes'].add(seqname)

                # Determine chromosome arm (simplified)
                chromosome = seqname.replace('chr', '')
                if 'centromere' in attrs.get('description', '').lower():
                    arm = 'centromere'
                else:
                    # This is simplified - in reality you'd need centromere coordinates
                    arm = 'unknown'
                gene_info[clean_name]['arms'].add(arm)

    print(f"Found {len(gene_counts)} unique gene names")
    print(f"Total gene copies: {sum(gene_counts.values())}")

    # Convert to list of dictionaries
    genome_wide_data = []
    for gene_name, count in gene_counts.items():
        info = gene_info[gene_name]
        genome_wide_data.append({
            'gene_name': gene_name,
            'genome_wide_copies': count,
            'gene_biotype': info['biotype'],
            'chromosomes': ','.join(sorted(info['chromosomes'])),
            'num_chromosomes': len(info['chromosomes'])
        })

    return genome_wide_data

def main():
    """Main function to build genome-wide copy background."""
    print("=== Building Genome-Wide Copy Background ===")

    # Build genome-wide copy counts
    genome_wide_data = build_genome_wide_copy_counts()

    # Load PHR gene copy summary
    print("\nLoading PHR gene copy summary...")
    phr_genes = {}
    with open('gene_copy_summary.csv', 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            phr_genes[row['gene_name']] = {
                'copies_in_phrs': int(row['total_copies']),
                'gene_biotype': row['gene_biotype']
            }

    print(f"PHR genes: {len(phr_genes)}")
    print(f"Genome-wide genes: {len(genome_wide_data)}")

    # Create comprehensive background by merging
    print("\nCreating comprehensive copy background...")
    comprehensive_data = []

    # Create lookup dict for genome-wide data
    genome_wide_dict = {item['gene_name']: item for item in genome_wide_data}

    # Get all unique gene names
    all_genes = set(phr_genes.keys()) | set(genome_wide_dict.keys())

    for gene_name in all_genes:
        phr_info = phr_genes.get(gene_name, {})
        genome_info = genome_wide_dict.get(gene_name, {})

        comprehensive_data.append({
            'gene_name': gene_name,
            'copies_in_phrs': phr_info.get('copies_in_phrs', 0),
            'genome_wide_copies': genome_info.get('genome_wide_copies', 0),
            'gene_biotype': phr_info.get('gene_biotype', genome_info.get('gene_biotype', 'unknown')),
            'num_chromosomes': genome_info.get('num_chromosomes', 0),
            'chromosomes': genome_info.get('chromosomes', ''),
            'in_phrs': gene_name in phr_genes,
            'in_genome': gene_name in genome_wide_dict
        })

    # Save results
    print("\nSaving results...")

    # Save comprehensive background
    with open('comprehensive_copy_background.csv', 'w', newline='') as f:
        if comprehensive_data:
            writer = csv.DictWriter(f, fieldnames=comprehensive_data[0].keys())
            writer.writeheader()
            writer.writerows(comprehensive_data)

    # Save genome-wide data
    with open('genome_wide_gene_copies.csv', 'w', newline='') as f:
        if genome_wide_data:
            writer = csv.DictWriter(f, fieldnames=genome_wide_data[0].keys())
            writer.writeheader()
            writer.writerows(genome_wide_data)

    # Print summary statistics
    print("\n=== Summary Statistics ===")
    print(f"Total unique genes in genome: {len(genome_wide_data)}")
    print(f"Total copies in genome: {sum(item['genome_wide_copies'] for item in genome_wide_data)}")
    print(f"Genes with PHR copies: {sum(1 for item in comprehensive_data if item['in_phrs'])}")
    print(f"Total PHR copies: {sum(item['copies_in_phrs'] for item in comprehensive_data)}")
    print(f"Genes in both PHRs and genome-wide: {sum(1 for item in comprehensive_data if item['in_phrs'] and item['in_genome'])}")

    # Show some examples
    print("\nExample genes with high copy numbers:")
    high_copy_genes = [item for item in comprehensive_data if item['genome_wide_copies'] > 10]
    high_copy_genes.sort(key=lambda x: x['genome_wide_copies'], reverse=True)

    for gene in high_copy_genes[:10]:
        print(f"{gene['gene_name']}: PHR={gene['copies_in_phrs']}, Genome={gene['genome_wide_copies']}, Type={gene['gene_biotype']}")

    print("\n=== Background Building Complete ===")

if __name__ == "__main__":
    main()