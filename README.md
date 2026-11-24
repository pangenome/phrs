# Concerted evolution of human subtelomeres

Code, figure inputs and manuscript source for the analysis of inter-chromosomal sequence sharing at human
chromosome ends across the 465 near-complete assemblies of the Human Pangenome Reference Consortium
release 2 (232 individuals plus CHM13), with three-dimensional contact maps and a three-generation pedigree.

This is a research repository, not a software package. Most analyses read large intermediate files
(all-against-all alignments, similarity matrices, contact maps) from the cluster path
`/moosefs/guarracino/HPRCv2/PHR_III/`; the tables needed to redraw the figures are in `data/`.

## Layout

| path | content |
| --- | --- |
| `submission/` | the manuscript: `paper.tex` (single file), `bibliography.bib`, `fig/` (main and Extended Data figures), `compile.sh` |
| `paper_prep/figures/<figure>/` | one directory per figure: the script that draws it, `sources.tsv` (inputs), `caption.md` and the rendered PDF/PNG |
| `data/` | figure and analysis inputs: PHR calls per flank, CHM13 PHR intervals, similarity matrix, community assignments, per-pair similarity and contact tables, CHM13 gene annotation |
| `scripts/phr-detection/` | contig classification and the PHR caller (`find-multichr-regions-incremental.py`) |
| `scripts/phr_flank5mb/` | PHR calls on 5 Mb flanks (extraction and all-vs-one alignment jobs) |
| `scripts/cladistics/` | duplicon copy trees (`duplicon_monophyly.py`, `plot_duplicon_trees.py`), pangenome gene families, enrichment controls |
| `scripts/conversion/` | gene-conversion tests: GENECONV wrapper, the cohort scans in 2 kb bins (`conversion_scan.py`) and at site level (`conversion_sites.py`), read and trio checks, known-answer tests (`test_*.py`) |
| `scripts/ci/` | divergence and confidence-interval analyses (paralogous versus allelic divergence, Fig. 5 inputs) and the audit of numbers quoted in the response |
| `scripts/hic/`, `scripts/human/`, `scripts/mouse/` | sequence-versus-contact analyses for Hi-C, Pore-C, CiFi and mouse meiotic Hi-C |
| `scripts/popgen/` | population differentiation (Hudson F<sub>ST</sub>) |
| `scripts/pedigree/` | the WashU pedigree analyses |
| `response/` | the point-by-point response to the referees and the notes behind each answer |
| `end-to-end-report/` | a 14-section report of the whole analysis (`end-to-end-report/README.md` is the table of contents) |
| `slides/` | conference decks |

## Running

Build the manuscript (needs `pdflatex` and `bibtex` only):

```
cd submission
bash compile.sh
```

Redraw a figure from its directory, for example Fig. 5 (Python, matplotlib). It reads its three inputs from
`/moosefs`, or from a local copy of them named by `FIG5_INPUT_DIR`:

```
cd paper_prep/figures/fig5_concerted
FIG5_INPUT_DIR=/path/to/inputs python3 figure_fig5_concerted.py
```

The R figure scripts run under guix; the manifest path is in each script's header comment.
The known-answer tests of the conversion scans need no cluster data:

```
python3 scripts/conversion/test_conversion_scan.py
python3 scripts/conversion/test_conversion_sites.py
```
