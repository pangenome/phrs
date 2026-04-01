#!/bin/bash

# Run GO enrichment analysis using g:Profiler web API for protein-coding genes only
# Read gene list
GENE_LIST=$(cat phrs.no_acro.coding_gene_names.txt | tr '\n' ' ')

echo "Running GO enrichment analysis for $(wc -l < phrs.no_acro.coding_gene_names.txt) protein-coding genes..."
echo "Gene list: $(cat phrs.no_acro.coding_gene_names.txt | tr '\n' ' ')"

# Create JSON payload for g:Profiler
cat > gprofiler_request_coding_only.json << EOF
{
  "organism": "hsapiens",
  "query": [$(cat phrs.no_acro.coding_gene_names.txt | sed 's/^/"/; s/$/",/' | tr -d '\n' | sed 's/,$//')],
  "sources": ["GO:BP", "GO:MF", "KEGG"],
  "user_threshold": 0.05,
  "significance_threshold_method": "fdr"
}
EOF

echo "Submitting query to g:Profiler..."

# Submit to g:Profiler API
curl -X POST \
  -H "Content-Type: application/json" \
  -d @gprofiler_request_coding_only.json \
  "https://biit.cs.ut.ee/gprofiler/api/gost/profile/" \
  -o gprofiler_results_coding_only.json

if [ $? -eq 0 ]; then
  echo "g:Profiler analysis completed successfully"
  echo "Results saved to gprofiler_results_coding_only.json"
else
  echo "Error: g:Profiler API call failed"
  exit 1
fi