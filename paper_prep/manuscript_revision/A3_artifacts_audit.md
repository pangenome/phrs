# A3 Artifacts Audit: construction paths, notebook citations, and internal pointers

Date: 2026-06-17
Task: manuscript-revision-a3

## Scope and search record

Audited public-facing manuscript files and figure-support files for construction artifacts that should not appear in a submission package:

- `submission/paper.tex`
- `submission/bibliography.bib`
- `submission/README.md`, `submission/BUILD_LOG.md`, `submission/notes/`
- `paper_prep/figures/*/caption.md`
- `paper_prep/figures/*/sources.tsv`

Search classes run:

- `/moosefs`, `/home/`, `Dropbox`, `Desktop`
- Markdown line pointers and draft pointers: `.md:[0-9]`, `.md line`, `L[0-9]`, `NATURE_DRAFT`, `REVISION_LOG`, `PEER_REVIEW`
- Pipeline report names and lab notebooks: `end-to-end-report`, `SURVEY_[0-9]`, `paper_prep/surveys`, `paper_prep/_brainstorming`
- Local scripts as citations or public provenance: `scripts/`, `submission/scripts`, `README`, `BUILD_LOG`
- Notebook/provenance vocabulary: `lab notebook`, `notebook`, `source-of-truth`, `absolute path`, `on disk`, `local copy`, `local snapshot`

## Summary

The active manuscript has one direct private-path leak in Data and Code Availability and several Methods sentences that cite local script paths as if they were public evidence. It also has one explicit internal Markdown line pointer in the Methods. These should be fixed before submission.

`submission/bibliography.bib` does not show `/moosefs`, line-number pointers, internal pipeline report citations, or local script-path citations under the searched patterns. It contains ordinary bibliographic notes and public URLs only for this artifact class.

The figure `sources.tsv` files contain many `/moosefs` paths and `SURVEY_*` back-references. These are useful internal build manifests and can remain in the repo, but they should not be submitted as public source-data notes without conversion to public repo-relative artifacts, accession identifiers, or literature citations.

## Actionable manuscript leaks

### A3-1. Private moosefs roots in Data and Code Availability

- File/line: `submission/paper.tex:641-642`
- Current text:

```tex
GitHub \texttt{ekg/phrs}; on-disk roots
\texttt{/moosefs/guarracino/HPRCv2/PHR\_III/} and \texttt{/moosefs/erikg/phrs/};
```

- Problem: Private filesystem roots are not public accessions and cannot be used by reviewers or readers. This is the highest-priority construction artifact because it is in the submission's Data and Code Availability section.
- Exact replacement:

```tex
Analysis code, manuscript source, figure-generation scripts, and vendored figure inputs are available at \url{https://github.com/ekg/phrs}. HPRC v2 assemblies, CHM13v2.0, WashU pedigree assemblies, CEPH1463 assemblies, and published 3D-genome datasets are available from their cited public releases and accessions; accession-level identifiers are listed with the corresponding citations and source-data tables.
```

- Follow-up note: If source-data tables are created for submission, each table should use public dataset names/citations rather than `/moosefs` roots. Examples: HPRC Data Release 2 / HPRC v2 citation for haplotype assemblies, T2T-CHM13 for the reference, Cechova et al. for WashU T2T pedigree assemblies, Porubsky et al. for CEPH1463, Tan 2018 for Dip-C, Xu 2025 for sperm scHi-C, Zuo 2021 / GEO accessions already noted in the internal report for mouse meiotic Hi-C, and Lalli 2025 for the T2T-CHM13 recombination map.

### A3-2. Local script path cited as bootstrap-method evidence

- File/line: `submission/paper.tex:517-521`
- Current text:

```tex
Support from a 1,000-replicate character-level bootstrap
(\texttt{scripts/cladistics/char\_bootstrap\_d\_m9.R}): each replicate
resamples PHR flanks with replacement, rebuilds the arm-level distance from
cached per-PHR cross-chromosome involvement (no PGGB re-run; $\approx 108$ s
single-core), and recomputes NJ + UPGMA.
```

- Problem: The local script path is useful for reproducibility but reads like a lab-notebook citation. The runtime detail is a construction note, not Methods content.
- Exact replacement:

```tex
Support was estimated by a 1,000-replicate character-level bootstrap: each replicate resampled PHR flanks with replacement, rebuilt the arm-level distance matrix from per-PHR cross-chromosome involvement, and recomputed NJ and UPGMA trees. The implementation is included in the public analysis-code repository.
```

### A3-3. Local strict-MAPQ script path in Methods

- File/line: `submission/paper.tex:548-550`
- Current text:

```tex
A strict-MAPQ re-binning from upstream \texttt{.allValidPairs} (both mates
MAPQ $\geq 30$; \texttt{scripts/hic/mapq\_strict\_d\_peerq1.py}) reproduces the
flanking B/W within Poisson noise while PHR-internal contact collapses to a
```

- Problem: The local script path is not a public source and the wording implies a completed numerical control. If the control is retained, cite it as a repository implementation and avoid private provenance.
- Exact replacement:

```tex
A strict-MAPQ re-binning from upstream valid-pair files (both mates MAPQ $\geq 30$) reproduced the flanking B/W within Poisson noise while PHR-internal contact collapsed to a noise floor.
```

- If the strict-MAPQ table is not available in a public source-data file, safer replacement:

```tex
A strict-MAPQ re-binning workflow for upstream valid-pair files (both mates MAPQ $\geq 30$) is included in the public analysis-code repository; the unique-sequence flanking analysis provides the reported artifact-control statistic.
```

### A3-4. Local script directory cited in reproducibility subsection

- File/line: `submission/paper.tex:572-577`
- Current text:

```tex
We re-aligned all Hi-C, Pore-C, CiFi and Dip-C data from raw FASTQ against the
corresponding T2T reference (CHM13v2.0 or matched HPRC v2 haplotype) with
multi-mappers retained.
Re-alignment scripts and command lines are at \texttt{scripts/hic-realign/}.
Anyone attempting to reproduce these results from existing deposited processed
files at default parameters will see no signal.
```

- Problem: `scripts/hic-realign/` is a repo-relative construction pointer. It is acceptable in Code Availability, but in the Methods it should be phrased as public repository content and paired with public raw-data accessions.
- Exact replacement:

```tex
We re-aligned all Hi-C, Pore-C, CiFi and Dip-C data from raw FASTQ against the corresponding T2T reference (CHM13v2.0 or matched HPRC v2 haplotype) with multi-mappers retained. The re-alignment workflows and command lines are included in the public analysis-code repository. Deposited processed files generated with default MAPQ filters do not preserve the inter-arm signal.
```

### A3-5. Internal Markdown report line pointer in Methods

- File/line: `submission/paper.tex:595-598`
- Current text:

```tex
PBMC Dip-C negative control: 18 cells with PHR coordinates projected from
CHM13 to hg19 via impg, yielding $W/B = 0.983$, $p = 0.305$
(\texttt{05\_hic\_validation.md} \S PBMC, L455-469).
```

- Problem: `05_hic_validation.md` and `L455-469` are internal lab-notebook pointers. A submitted paper should cite the public Dip-C source and report the statistic without a private line reference.
- Exact replacement:

```tex
PBMC Dip-C negative control: 18 cells from the public Dip-C dataset \citeMeth{Tan2018}, with PHR coordinates projected from CHM13 to hg19 via impg, yielding $W/B = 0.983$ and $p = 0.305$.
```

### A3-6. Local mouse Mantel script path in Methods

- File/line: `submission/paper.tex:615-616`
- Current text:

```tex
\texttt{vegan::mantel} cross-check within $\pm 0.02$
(\texttt{scripts/mouse/mantel\_d\_m5.\{py,R\}}) \citeMeth{Zuo2021}.
```

- Problem: A local script path is being used as a method citation. Keep the public Zuo citation and describe the cross-check generally.
- Exact replacement:

```tex
\texttt{vegan::mantel} cross-check within $\pm 0.02$ \citeMeth{Zuo2021}; the Python and R implementations are included in the public analysis-code repository.
```

### A3-7. Local pedigree-classification script path in Methods

- File/line: `submission/paper.tex:620-624`
- Current text:

```tex
\texttt{odgi untangle nth-best=1} per flank; high-quality filter
\texttt{min\_score} $\geq 0.8$ with $500\ \mathrm{bp} \leq \mathrm{size} \leq 100\ \mathrm{kb}$;
within-Leiden filter as credibility constraint; pattern classification via
\texttt{scripts/pedigree/analyze-pedigree-recombination.py}
\citeMeth{Cechova2025,Porubsky2025}.
```

- Problem: The script path is a construction pointer. The public scientific sources are the pedigree assembly papers/releases; code location belongs in Code Availability.
- Exact replacement:

```tex
\texttt{odgi untangle nth-best=1} per flank; high-quality filter
\texttt{min\_score} $\geq 0.8$ with $500\ \mathrm{bp} \leq \mathrm{size} \leq 100\ \mathrm{kb}$;
within-Leiden filter as credibility constraint; pattern classes were assigned from patch geometry and parent/child haplotype consistency
\citeMeth{Cechova2025,Porubsky2025}. The classification implementation is included in the public analysis-code repository.
```

### A3-8. Local Monte Carlo script path in Methods

- File/line: `submission/paper.tex:628-631`
- Current text:

```tex
The within-Leiden fraction was tested against a permutation null preserving
per-arm patch-count marginals ($B = 10000$; per-community BH correction;
\texttt{scripts/pedigree/monte\_carlo\_null\_d\_m4.py}); Wilson 95\% CI on
494/538 is 89.2--93.9\%.
```

- Problem: Same local-script-as-citation issue.
- Exact replacement:

```tex
The within-Leiden fraction was tested against a permutation null preserving per-arm patch-count marginals ($B = 10000$; per-community BH correction); Wilson 95\% CI on 494/538 is 89.2--93.9\%. The permutation implementation is included in the public analysis-code repository.
```

### A3-9. Script path in Data and Code Availability

- File/line: `submission/paper.tex:646-648`
- Current text:

```tex
Re-alignment scripts and command lines for Hi-C/Pore-C/CiFi/Dip-C are at
\texttt{scripts/hic-realign/}; deposited processed files (default MAPQ
$\geq 30$) are insufficient to reproduce the inter-arm signal.
```

- Problem: This is less severe than A3-4 because Code Availability can refer to repository paths. However, it should anchor the path to the public GitHub repository and raw-data accessions.
- Exact replacement:

```tex
The public repository includes re-alignment workflows and command lines for Hi-C, Pore-C, CiFi and Dip-C. These workflows start from the cited public raw-read datasets; deposited processed files generated with default MAPQ $\geq 30$ filters are insufficient to reproduce the inter-arm signal.
```

### A3-10. Internal figure-regeneration comment in manuscript source

- File/line: `submission/paper.tex:39`
- Current text:

```tex
% Figures are regenerated by submission/scripts/figures/ (see its README.md).
```

- Problem: This is a LaTeX comment and will not print in the manuscript PDF. It is not a reader-facing leak, but it is still embedded in the submission source. If raw `.tex` is uploaded, it exposes internal construction structure.
- Exact replacement if source cleanliness is required:

```tex
% Figure-generation workflows are provided in the public analysis-code repository.
```

## Bibliography audit

No actionable construction-artifact leaks found in `submission/bibliography.bib` under the searched patterns.

Non-issues:

- `submission/bibliography.bib:691` contains `https://github.com/pangenome/impg`, a public URL and appropriate software reference.
- Bibliography `note` fields contain scientific context such as "co-localisation" and "used in this manuscript"; these are not internal path or notebook leaks.

## Figure captions in `paper_prep/figures`

These captions are not currently the active `submission/paper.tex` legends, but they are likely upstream sources for figure legends. They should not be exported unchanged where they cite surveys or internal source files.

### A3-11. Survey citations in upstream captions

- Files/lines:
  - `paper_prep/figures/fig1/caption.md:21`
  - `paper_prep/figures/fig2/caption.md:6`
  - `paper_prep/figures/fig3/caption.md:6`
  - `paper_prep/figures/fig3/caption.md:20`
  - `paper_prep/figures/fig4/caption.md:16`
  - `paper_prep/figures/fig4/caption.md:22`
  - `paper_prep/figures/fig4/caption.md:28`
  - `paper_prep/figures/ed2/caption.md:10`
  - `paper_prep/figures/ed8/caption.md:3`
- Problem: `SURVEY_04`, `SURVEY_07`, `SURVEY_08`, `SURVEY_09`, and `SURVEY_14` are internal evidence inventories, not public citations.
- Replacement rule: Remove the survey back-reference from any public caption. If a citation is needed, cite the relevant paper or the public data source in the manuscript legend or Methods.
- Exact examples:
  - Replace "`SURVEY_14 §1.6`" with "cross-assembler CEPH1463 features defined as described in Methods".
  - Replace "`SURVEY_09 §1.3`" with "the known RPE-1 t(X;10) rearrangement \cite{hic3d_Volpe2025RPE1}" if that panel is retained.
  - Replace "`SURVEY_08 §1.7`" with "mouse meiotic Hi-C from Zuo et al. \cite{Zuo2021}".
  - Replace "`SURVEY_07 §1.2`" with "flanking unique-sequence control described in Methods".
  - Replace "`SURVEY_04 §1.1`" with "allele-vs-paralog distances computed from the pangenome Jaccard matrix as described in Methods".

## Figure source manifests: internal docs that should remain internal

The `paper_prep/figures/*/sources.tsv` files are internal build manifests. They correctly record absolute `/moosefs` inputs and survey back-references needed to rebuild figures in the lab environment. They should remain internal unless converted to public source-data manifests.

### A3-12. Private `/moosefs` paths in figure source manifests

- Files/lines:
  - `paper_prep/figures/fig1/sources.tsv:6-10`
  - `paper_prep/figures/fig2/sources.tsv:2-10`
  - `paper_prep/figures/fig3/sources.tsv:2-25`
  - `paper_prep/figures/fig4/sources.tsv:2-3`, `:7`, `:9-20`
  - `paper_prep/figures/ed1/sources.tsv:3-5`
  - `paper_prep/figures/ed2/sources.tsv:2-4`, `:6-9`
  - `paper_prep/figures/ed3/sources.tsv:2-6`
  - `paper_prep/figures/ed4/sources.tsv:6`
  - `paper_prep/figures/ed5/sources.tsv:2-6`
  - `paper_prep/figures/ed8/sources.tsv:5-8`
- Recommendation: Do not submit these manifests as-is. If a public source-data manifest is required, create a separate export with columns such as `panel`, `public_source`, `public_accession_or_citation`, `repo_relative_processed_file`, and `note`.
- Replacement pattern:

```tsv
panel	public_source	public_accession_or_citation	repo_relative_processed_file	note
Fig4A	WashU T2T pedigree assemblies	Cechova et al. 2025; public release cited in bibliography	data/<public-or-vendored-table>.tsv	Inter-chromosomal patch summary used for the plotted statistic.
Fig4C	Mouse meiotic Hi-C	Zuo et al. 2021; GEO accessions listed in source-data table	data/<public-or-vendored-table>.tsv	Stage-resolved 20 kb per-PHR-pair contact summaries.
```

### A3-13. Internal survey and report paths in figure source manifests

- Files/lines:
  - `paper_prep/figures/fig1/sources.tsv:2`, `:4`, `:6`, `:9`
  - `paper_prep/figures/fig2/sources.tsv:2-10`
  - `paper_prep/figures/fig4/sources.tsv:4-6`, `:8`, `:22-24`
  - `paper_prep/figures/ed1/sources.tsv:6`
  - `paper_prep/figures/ed3/sources.tsv:7`
  - `paper_prep/figures/ed4/sources.tsv:7-8`
  - `paper_prep/figures/ed5/sources.tsv:7`
  - `paper_prep/figures/ed8/sources.tsv:2-4`, `:9-10`
- Problem: `SURVEY_*`, `paper_prep/surveys/*.md`, and `end-to-end-report/...` are internal notebooks. They are appropriate for provenance inside this repo but should not be treated as source-data citations.
- Recommendation: Keep the existing internal manifests for rebuilds. For submission-facing manifests, replace internal survey paths with public citations or public repository-relative generated inputs.
- Exact examples:
  - `paper_prep/surveys/SURVEY_14_pedigree_recombination.md` -> "Cechova et al. 2025 WashU T2T pedigree assemblies; Porubsky et al. 2025 CEPH1463 assemblies; processed patch summary in public repository."
  - `paper_prep/surveys/SURVEY_08_mouse.md` -> "Zuo et al. 2021 mouse meiotic Hi-C; B6 and CAST T2T assemblies cited in Methods; processed per-stage tables in public repository."
  - `end-to-end-report/pedigree-plots/washu/*.pdf` -> "Rendered pedigree untangle panels vendored in `submission/fig/MainFigures/` or regenerated from public repository workflows."

## Submission support docs: internal and should not be bundled as manuscript evidence

These files are under `submission/` but are not part of the active manuscript text. They should remain internal or be omitted from any journal upload unless rewritten.

### A3-14. `submission/BUILD_LOG.md`

- Lines:
  - `submission/BUILD_LOG.md:68`
  - `submission/BUILD_LOG.md:70`
  - `submission/BUILD_LOG.md:143`
- Current artifacts include a local `/home/guarracino/...` source, "local copies", and `NATURE_DRAFT_v6.md`.
- Recommendation: Do not upload `BUILD_LOG.md`. If a build record is required, replace with a minimal public statement: "The manuscript was built from `submission/paper.tex` with the bundled Springer Nature class and bibliography files; figure assets are included under `submission/fig/`."

### A3-15. `submission/README.md`

- Lines:
  - `submission/README.md:130-160`
  - `submission/README.md:183`
- Current artifacts include stale build-log notes, `scripts/figures/README.md`, and `paper_prep/synthesis/NATURE_DRAFT_v6.md`.
- Recommendation: Keep as an internal developer README. Do not use it as a submission Methods/Data Availability source. If public repository documentation is needed, rewrite around stable public paths and remove stale/frozen-draft references.

### A3-16. `submission/notes/20260521_pptx-slide-image-provenance.md`

- Lines:
  - `submission/notes/20260521_pptx-slide-image-provenance.md:4`
  - `:26`, `:39`, `:53`, `:78`, `:92-104`, `:120`
- Current artifacts include `/home/guarracino/...`, `/moosefs/...`, and slide-build provenance.
- Recommendation: Do not bundle this note with the manuscript. It is an internal provenance note for slide-derived placeholder figures, not a public source-data manifest.

## Internal docs that should remain internal, not fixed by A3

The broad search found many `/moosefs`, `end-to-end-report`, `SURVEY_*`, and script-path references in `paper_prep/synthesis/`, `paper_prep/surveys/`, and `end-to-end-report/report/`. These are expected for internal lab notebooks, reviewer-response analyses, and build documentation. They should not be mechanically deleted.

Rule for downstream A fan-in:

- Edit `submission/paper.tex` for A3-1 through A3-10.
- Leave `submission/bibliography.bib` unchanged for this artifact class.
- Do not edit `paper_prep/figures/*/sources.tsv` unless creating a separate public export. The existing files are internal rebuild manifests.
- Do not submit `submission/BUILD_LOG.md`, `submission/README.md`, or `submission/notes/` as manuscript evidence without rewriting.

## Validation checklist

- Artifact exists: this file.
- Searches covered `/moosefs`, Markdown line pointers, pipeline report names, survey names, and local scripts as citations.
- Manuscript text leaks are separated from internal repo docs and figure build manifests.
- Exact replacement wording is provided where the fix is obvious.
