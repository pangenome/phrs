#!/bin/bash

# Run GO enrichment analysis using g:Profiler web API
# Read gene list
GENE_LIST=$(cat gene_list_for_gprofiler.txt | tr '\n' ' ')

echo "Running GO enrichment analysis for $(wc -l < gene_list_for_gprofiler.txt) genes..."
echo "Gene list sample: $(head -5 gene_list_for_gprofiler.txt | tr '\n' ' ')..."

# Create JSON payload for g:Profiler
cat > gprofiler_request.json << EOF
{
  "organism": "hsapiens",
  "query": [$(cat gene_list_for_gprofiler.txt | sed 's/^/"/; s/$/",/' | tr -d '\n' | sed 's/,$//')],
  "sources": ["GO:BP", "GO:MF", "KEGG"],
  "user_threshold": 0.05,
  "significance_threshold_method": "fdr"
}
EOF

echo "Submitting query to g:Profiler..."

# Submit to g:Profiler API
curl -X POST \
  -H "Content-Type: application/json" \
  -d @gprofiler_request.json \
  "https://biit.cs.ut.ee/gprofiler/api/gost/profile/" \
  -o gprofiler_results.json

if [ $? -eq 0 ]; then
  echo "g:Profiler analysis completed successfully"
  echo "Results saved to gprofiler_results.json"
else
  echo "Error: g:Profiler API call failed"
  exit 1
fi