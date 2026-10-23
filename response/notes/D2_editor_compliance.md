# D2 — Editor compliance draft (text only)

> Provenance for this note: see `response/PROVENANCE.md` (tool versions: `response/notes/tool_versions_and_pins.md`).

Authority: the editor's letter at the top of `PEER_REVIEW.md` (lines 1–33). Every requirement
below is quoted or paraphrased from that letter only; nothing here is inferred from a simulated
review (that file was not opened) or invented.

This note drafts the internal tracking row **"Editor compliance"** / item 17 in
`response/POINT_BY_POINT.md` (table at :47) and the mandatory list at `POINT_BY_POINT.md`:198–201.
That file itself is NOT edited here (outside my write scope); the mapping is for the parent to
apply.

Statuses used throughout:

- **[DRAFTED]** — text supplied here, ready to paste.
- **[CHECKED]** — I verified the current state of the repo/manuscript and recorded the evidence.
- **[AUTHOR ACTION]** — a factual or administrative step that only an author can take.
- **[AUTHOR TO SUPPLY]** — a fact (name, number, wording) that is not in this repository and
  that I will not invent.

---

## 0. What the letter asks for, in its own words

| # | Requirement (`PEER_REVIEW.md` line) | Status today [CHECKED] |
| --- | --- | --- |
| 1 | "must include a Data availability statement … conditions of access to the 'minimum dataset'" (:21) | present, but as a **combined** `\subsection*{Data and code availability}` inside Methods (`submission/paper.tex`:889–910) |
| 2 | "a Code availability statement must be included … as a separate section after the Data availability statement but before the references … deposited in a DOI-minting repository such as Zenodo … and cited in the reference list" (:23) | not separate, not deposited, not cited in the reference list |
| 3 | "the Methods section must include an ethics statement … name of the committee that approved the study … guidelines and regulations … informed consent was obtained" (:27) | **absent** — `grep -n -i "consent\|ethic\|IRB\|approval" submission/paper.tex` returns only an affiliation line (:1113) |
| 4 | "Extended Data should be assembled into a maximum of 10 A4 size, multi-panelled display items, submitted as individual files in .jpg, .tif or .eps format only" (:31) | 3 items (limit is fine) but all three are **.pdf** |
| 5 | "Source data underlying graphs should be provided in spreadsheet format (.xls, .xlsx, or .csv) … one file per figure" (:33) | none submitted as a Source Data file |
| 6 | "include the final submission word count … within your covering letter" (:19) | not computed in any repo document I could find |
| 7 | "we also strongly suggest that your revised manuscript has tracked changes" (:15) | none |
| 8 | "All of the checklists provided with the current submission (Reporting summary, and Code and software checklist (if applicable)) should be updated" (:20) | no Reporting Summary or Code/software checklist file exists in this repo |
| 9 | "fill out a Third party rights table and upload this with your revised manuscript" (:35) | none |
| 10 | LENGTH: "~2,500–4,300 words" (:19) | main text 3,330 words under the counting rule in §4.3 → within range |

---

## 1. Data availability statement — draft

[CHECKED] The two URLs and the BioProject accession already used by the manuscript are
`https://github.com/human-pangenomics/hprc_intermediate_assembly/tree/main/data_tables`
(`submission/paper.tex`:892), `https://humanpangenome.org/hprc-data-release-2/`
(`submission/bibliography.bib`:739, with "BioProject PRJNA730823" in its note and in the
Acknowledgements at `submission/paper.tex`:1165), and
`https://github.com/biomonika/washu-pedigree` (`submission/paper.tex`:905).

[DRAFTED] — *Data availability.* (to replace `submission/paper.tex`:889–910, split from code)

> **Data availability**
>
> The haplotype-resolved assemblies analysed in this study are the Human Pangenome Reference
> Consortium (HPRC) v2 release: 232 individuals from five superpopulations, 464
> haplotype-resolved assemblies, released as HPRC Data Release 2
> (<https://humanpangenome.org/hprc-data-release-2/>) with the underlying sequencing and
> assembly data tables at
> <https://github.com/human-pangenomics/hprc_intermediate_assembly/tree/main/data_tables>
> (BioProject PRJNA730823). The complete telomere-to-telomere reference T2T-CHM13v2.0 is the
> assembly of ref. [Nurk2022]. The three-generation telomere-to-telomere pedigree assemblies
> from Washington University are available at <https://github.com/biomonika/washu-pedigree>
> [Cechova2025].
>
> The **minimum dataset** needed to interpret, verify and extend the conclusions — the
> per-flank pseudo-homolog region (PHR) calls, the arm-level similarity matrix and community
> assignments, the per-pair sequence-versus-contact tables, and the pedigree PHR calls — is
> provided as plain tab-separated tables in the code repository below and archived at
> [DOI — AUTHOR ACTION], specifically:
>
> | table | content |
> | --- | --- |
> | `data/all-vs-all.1Mb.p95.id95.len.tsv` | one row per telomere-anchored 500 kb flank for all 465 assemblies; called PHR interval, partner chromosome set and partner arm set |
> | `data/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv` | arm-level (41 × 41) subtelomeric similarity matrix underlying Fig. 2B |
> | `data/hprcv2.1Mb.subtelo.arm-leiden.communities.tsv`, `…arm-leiden-k15.assignments.tsv` | arm-level Leiden communities (C1–C15) and per-arm assignments, as used in Fig. 2B and Extended Data Fig. 3 |
> | `data/fig5_washu.all-vs-all.1Mb.p95.id95.len.tsv` | pedigree whole-genome flank comparison underlying Fig. 5 |
> | `data/human_CHM13_hic_*_seqlevel.tsv`, `data/human_HG002_{hic,porec,cifi}_*_seqlevel.tsv`, `data/mouse_{leptotene,zygotene,pachytene,diplotene}_phr_*_seqlevel.tsv` | per-PHR-pair similarity and length-normalized contact values underlying Fig. 4 and Extended Data Fig. 1 |
> | `scripts/cladistics/phr_by_superpopulation.tsv` | PHR count, PHR length and community membership by superpopulation |
> | `scripts/popgen/matched_fst_d_m6.tsv` | arm × community-block Hudson *F*<sub>ST</sub> against published continental baselines |
> | `submission/supporting_material/*.tsv` | Extended Data Fig. 2 source tables and the Fig. 5 per-call margin table |
>
> The chromosome-contact and annotation datasets we re-analysed are public: HPRC v2 Hi-C
> (CHM13, HG002 and four further samples), CHM13 Arima Hi-C, and the CHM13v2.0 RefSeq/Liftoff
> v5.2 annotation are in the HPRC data tables and S3 bucket
> (<https://github.com/human-pangenomics/hprc_intermediate_assembly>); HG002 Pore-C is the
> UCSC/GIAB R10.4.1 nanopore Pore-C release in the same bucket; HG002 CiFi is ENA run
> ERR14656265 in BioProject PRJEB83708 [hic3d_cifi2025]; mouse meiotic Hi-C is GEO GSE155638
> (runs SRR12376885–SRR12376892) [Zuo2021]; the B6 and CAST T2T assemblies are in ENA
> BioProject PRJEB47108 (GCA_964188535 and GCA_964188545) [Francis2025]. No new human genomic
> data were generated in this study.

Notes on the draft, all [CHECKED]:

- The working URLs are the ones the current manuscript already uses; I did not invent any new
  one. The four re-analysis datasets (Pore-C, CiFi, mouse Hi-C, B6/CAST) are carried over
  verbatim from `submission/paper.tex`:896–906 (the paragraph beginning "The
  chromosome-contact and annotation datasets we reanalyzed are public"), not restated from
  memory.
- Table paths in the table above all exist (`ls data/` — 63 entries; per-file presence
  confirmed for each path named). Exception and **[AUTHOR ACTION]**: the *per-window* IMPG
  table from which the PHRs were called (`all-vs-all.1Mb.p95.id95.len.per_window.tsv`) is
  **not** vendored in this repo (`ls data | grep per_window` returns nothing); it is only on
  the `/moosefs` mount. Either deposit it too, or state that the derived per-flank table is the
  deposited minimum dataset — an author's decision about size versus completeness.
- **[AUTHOR TO SUPPLY]** the archival DOI (Zenodo/Figshare/Dryad) and the access conditions
  wording that Nature wants ("conditions of access to the minimum dataset"). If the deposit is
  public and unrestricted, say so explicitly.
- The existing text says "Analysis code, manuscript source, figure-generation scripts, and
  figure input files are available at …" — that sentence belongs in Code availability, not
  here; the split in §2 does that.

---

## 2. Code availability statement — draft

[CHECKED] `submission/paper.tex`:893–894 currently reads "Analysis code, manuscript source,
figure-generation scripts, and figure input files are available at
\url{https://github.com/pangenome/phrs}." [CHECKED] The repository's own git remote
(`.git/config`, read as a file; no git command run) is `git@github.com:ekg/phrs.git`, so the
public org/name in the statement and where the DOI will actually point must be reconciled —
**[AUTHOR ACTION]**.

[DRAFTED]

> **Code availability**
>
> All analysis code, manuscript source, figure-generation scripts and figure input files are
> available at <https://github.com/pangenome/phrs>, archived at [DOI — AUTHOR ACTION]. The
> pipeline is organised by analysis stage: `scripts/phr-detection/` (flank extraction and
> the incremental multi-chromosome region caller), `scripts/popgen/` (arm-level Jaccard
> similarity, community detection, the Hudson *F*<sub>ST</sub> block-jackknife comparison),
> `scripts/hic/` (the PHR-pair versus contact analyses for Hi-C, Pore-C and CiFi),
> `scripts/mouse/` (mouse meiotic Hi-C), `scripts/pedigree/` (pedigree direct-alignment
> windows) and `scripts/cladistics/` (tree resampling controls and the descriptive tables).
> Only standard Python 3 and R are required; no proprietary software is used. The code is
> released under [LICENCE — AUTHOR ACTION].

[DEPOSIT AS AUTHOR ACTION] — the letter is explicit: "Code should be deposited in a
DOI-minting repository such as Zenodo, Gigantum or Code Ocean and cited in the reference
list" (`PEER_REVIEW.md`:23). Concretely, still to do:

1. mint the DOI (Zenodo GitHub integration on a tagged release is the shortest route);
2. add a `LICENSE` file — **[CHECKED]** there is no `LICENSE`/`COPYING` at the repo root
   today, and the letter asks for "a license approved by the open source initiative";
3. add the DOI to the reference list as a numbered citation (the letter's "cited in the
   reference list");
4. **[CHECKED]** there is no `VERSIONS.md` in the repo (`find -name "VERSIONS*"` → empty),
   although `paper_prep/synthesis/SCRIPT_INVENTORY.md` expects one with pinned commits. Pin
   the analysis commits before minting so the DOI resolves to exactly the analysed code.
5. reconcile the two repository URLs (§2, first [CHECKED] note) and state any access
   restriction; the letter asks the statement to indicate "whether and how the code … can be
   accessed, including any restrictions to access".

---

## 3. Ethics statement — draft, with the parts we can attest to separated out

The letter (`PEER_REVIEW.md`:27) requires, in Methods: the **name of the approving committee**,
a confirmation that the research followed **all relevant guidelines and regulations**, and a
confirmation that **informed consent was obtained from all participants** (or the details of
who granted an exemption).

[DRAFTED] — attestable part, i.e. what the manuscript itself can state without new facts:

> **Ethics statement**
>
> This study did not generate new human-participant data. All human sequence data analysed
> here are from previously published, publicly released assembly datasets: the HPRC v2
> release of 232 individuals (HPRC Data Release 2, BioProject PRJNA730823; ref. [hprc v2]),
> the T2T-CHM13v2.0 reference assembly [Nurk2022], and the three-generation
> telomere-to-telomere pedigree assemblies from Washington University [Cechova2025]
> (<https://github.com/biomonika/washu-pedigree>). Each release was made under the consent
> and data-use terms of the original study that produced it, and the assemblies we used are
> public and de-identified; we re-analysed the released assemblies only, with no access to
> individual-level identifiable data. [COMMITTEE / APPROVAL NUMBERS / CONSENT WORDING —
> AUTHOR TO SUPPLY]

[CHECKED] Why the bracketed part is not filled in: `grep -rn -i "informed consent|ethics|IRB|
de-identif"` over the repo's `.md`/`.tex`/`.bib` files returns no ethics text of any kind.
There is no approval number, no committee name and no consent wording anywhere in this
repository, so any that appeared here would be fabricated.

[AUTHOR TO SUPPLY] — the items the letter names explicitly, with where they must come from:

| item | who has it | why I cannot supply it |
| --- | --- | --- |
| Name of the committee that approved the **HPRC** component (e.g. the relevant IRB/Ethics Committee of the consortium member institution(s)) | HPRC / the consortium member institutions that collected the samples | no such name in the repo |
| Name of the committee that approved the **WashU pedigree** component, and its approval/protocol number | the Washington University authors of the pedigree study [Cechova2025] | no such name in the repo |
| Confirmation that the research was "performed in accordance with all relevant guidelines and regulations", and which guidelines (e.g. Declaration of Helsinki) | corresponding authors | a claim only authors can make |
| Consent wording: what participants consented to, including consent for public data release, or the exemption details if consent/approval was waived | HPRC and the pedigree study's authors | a claim only authors can make |
| Whether any sample is from a population or community with additional governance (e.g. the HPRC's sample-governance statements) | HPRC | not in the repo |

[DRAFTED] — the two "not applicable" statements, which are safe because the paper has neither
embryos/gametes/stem cells nor animal work of its own (it re-analyses mouse Hi-C, GEO
GSE155638):

> Human embryos, gametes and stem cells were not used in this study. No animals were used in
> this study; the mouse data analysed are previously published meiotic Hi-C maps (GEO
> GSE155638) [Zuo2021].

[CHECKED] Placement: the letter requires this **in the Methods section**
(`PEER_REVIEW.md`:27), i.e. a `\subsection*{Ethics statement}` inside `\section*{Methods}`
(`submission/paper.tex`:598 onward), not in the Acknowledgements.

---

## 4. Format checklist — remaining items, with the current state verified

### 4.1 Extended Data: format and count

[CHECKED] There are exactly **3** Extended Data items, all as PDF, in
`submission/fig/ExtendedDataFigures/`:

| item | current file | needs |
| --- | --- | --- |
| ED Fig. 1 (contact-assay companions) | `ED_Fig1_human_contacts.pdf` | `.tif` or `.eps` (letter allows `.jpg` too) at print quality |
| ED Fig. 2 (copy-number-aware enrichment heatmap) | `ED_Fig2_copy_number_ontology_heatmap.pdf` | `.tif`/`.eps` |
| ED Fig. 3 (module maps / architecture) | `ED_Fig3_module_maps.pdf` | `.tif`/`.eps` |

- Count is within the "maximum of 10 A4 size, multi-panelled display items".
- Legends are already "listed sequentially at the end of the main text, not in the Extended
  Data files" (`submission/paper.tex`:912–976), which matches the letter.
- **[AUTHOR ACTION]** rasterise/convert the three PDFs to `.tif` (or `.eps`) at the quality of
  the print figures; the LaTeX build can keep the PDFs while the submission files are `.tif`.
- **[AUTHOR TO SUPPLY]** decision on two orphan panels that exist in the same directory and are
  referenced nowhere in the manuscript: `ED_Fig2a_PAN027_maternal_recombination.png` and
  `ED_Fig2b_PAN028_maternal_recombination.png`. Either they become an Extended Data item (the
  count allows up to 10) or they should be dropped from the submission set.

### 4.2 Source Data: one file per figure

[CHECKED] No Source Data file is submitted today; `submission/supporting_material/` holds
Extended Data source tables and the Fig. 5 margin table only. The materials that would become
per-figure Source Data files, and whether the underlying table exists [CHECKED]:

| figure | source table(s) available in repo | state |
| --- | --- | --- |
| Fig. 1 | `data/all-vs-all.1Mb.p95.id95.len.tsv` (PHR lengths/positions, panel B), `data/hprc25272*.tsv.gz` (per-window partner counts, panel A) | **[AUTHOR ACTION]** assemble one `.xlsx`/`.csv` per figure |
| Fig. 2 | `data/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`, `data/hprcv2.1Mb.subtelo.arm-leiden.communities.tsv`, `submission/fig/MainFigures/arm_order_community.tsv`, `arm_order_tree.tsv`, `community_blocks.tsv`, `data/fig2a_node_community.tsv.gz` | **[AUTHOR ACTION]** assemble |
| Fig. 3 | `submission/fig/MainFigures/Fig3a/b/c` panels drawn from `scripts/cladistics/arm_module_*.tsv`, `data/fig2a_node_community.tsv.gz` | **[AUTHOR ACTION]** assemble |
| Fig. 4 | `data/human_CHM13_hic_*_seqlevel.tsv`, `data/human_HG002_{hic,porec,cifi}_*_seqlevel.tsv`, `data/hg002_porec_contact_matrix_20000bp.tsv`, `data/mouse_{leptotene,zygotene,pachytene,diplotene}_phr_*_seqlevel.tsv` | **[AUTHOR ACTION]** assemble |
| Fig. 5 | `data/fig5_washu.all-vs-all.1Mb.p95.id95.len.tsv`, `data/fig5_PAN027pat_vs_PAN011_joint.class_winners.impg_similarity.tsv.gz` (and the PAN028mat/PAN027mat/PAN010 counterparts), `submission/supporting_material/Fig5_margin_table.tsv`, `Fig5_pedigree_supporting_material.tsv`, `submission/fig/MainFigures/Fig5D_highlights.tsv` | **[AUTHOR ACTION]** assemble |

- The letter's format: "spreadsheet format (.xls, .xlsx, or .csv) … Submit one file per figure
  … and indicate the relevant figure number(s) in the file description." The formats above are
  all TSV; the deposit needs XLSX/CSV per figure.
- The letter's caveat: "Source data should not be submitted as Extended Data."

### 4.3 Word count for the cover letter

[CHECKED] Counted with the rule: take the main text of `submission/paper.tex` (line 95, the
first paragraph of the Introduction, through line 421, the last line of the Discussion; line
422 is blank and line 423 is `\subsection*{Acknowledgments}`), strip `\cite{}`, `\ref{}`, `\label{}`, `\url{}`
and `$…$` math, then `wc -w`:

| block | lines | words |
| --- | --- | --- |
| main text (Introduction → Discussion) | 95–421 | **3,330** |
| abstract | 65–91 | 228 |
| Methods | 598–910 | 3,132 |
| figure legends | 451–597 | 913 |
| Extended Data legends | 912–976 | 393 |

- The letter's window is "~2,500–4,300 words" excluding abstract, methods and references
  (`PEER_REVIEW.md`:19), so 3,330 is inside the window by ~370 words of headroom before the
  upper bound — the **cover-letter number is 3,330** under this rule, and the exact rule must
  be restated in the letter so the editor can reproduce it.
- **[AUTHOR ACTION]** recompute after the final revision (my number will be stale then), and
  state in the cover letter that the count excludes abstract, Methods, references and figure
  legends.
- Caveat on my number: the stripping rule above is crude (macro arguments become single
  tokens), so treat 3,330 as ±1 %. It is a same-day measurement of the current file, not a
  quoted figure.

### 4.4 Tracked-changes version

[AUTHOR ACTION] The letter "strongly suggests" a tracked-changes version
(`PEER_REVIEW.md`:15). For a LaTeX manuscript the practical route is `latexdiff` between the
submitted and revised `paper.tex`, producing a separate tracked-changes PDF; it needs the
original submitted `.tex`, which is not in this worktree.

### 4.5 Updated Reporting Summary and Code/software checklist

[AUTHOR TO SUPPLY] **[CHECKED]** no Reporting Summary or Code and software checklist exists
anywhere in this repository (`find . -iname "*reporting*" -o -iname "*code*checklist*" -o
-iname "*software*checklist*"` → empty). These are filled-in web forms uploaded with the
submission; the reported statistics must be updated to reflect the revisions (new *n*, the
added analyses, the *F*<sub>ST</sub> jackknife intervals, the new descriptive superpopulation
table).

### 4.6 Third party rights table

[AUTHOR TO SUPPLY] The letter requires the table with the revised manuscript
(`PEER_REVIEW.md`:35). Candidate third-party items to enumerate, from what the manuscript
contains: the CHM13v2.0 RefSeq/Liftoff v5.2 annotation, GEO GSE155638 mouse Hi-C, ENA
PRJEB47108 B6/CAST assemblies, ENA PRJEB83708 CiFi, and the HPRC v2 data tables. All are
re-analysed public data rather than reproduced display items, which is the point to state
explicitly in the table.

### 4.7 Also worth confirming while restructuring (not an editor item in the list above)

[CHECKED] `submission/paper.tex` currently has **one** section,
`\subsection*{Data and code availability}` (:889), placed inside Methods. The letter wants
Data availability, then **Code availability as its own section "after the Data availability
statement but before the references"** (:23). The restructure is **[AUTHOR ACTION]** in
`paper.tex`, which is outside my write scope.

---

## 5. What this note does NOT establish

1. No statement here has been placed into `submission/paper.tex` — I may not modify existing
   files. §1–§3 are drafts for an author to paste, and §4 is a checklist.
2. **[AUTHOR TO SUPPLY]** items are exactly that: committee names, approval numbers, consent
   wording, the archival DOI, the licence choice, the Reporting Summary and the third-party
   rights table. Nothing in this note substitutes for them.
3. The data availability draft names only URLs and accessions already present in the
   manuscript or the bibliography. I did not verify that any of them still resolves (no
   network check was run), and I did not verify the Zenodo/Figshare/Dryad policy text the
   editor links to.
4. The word count 3,330 is a measurement of the current `paper.tex` under a stated crude rule,
   not a recomputation of Nature's own counter, and it will be stale after the revision.
5. No claim is made that the remaining Extended Data conversions, Source Data assembly or the
   ethics wording will satisfy the editor — that is her call.
