# Slide Table: GO/Functional Context With Labels

Use as a compact methods-contrast table. It should not replace the report-backed HPRCv2 gene-family table.

| Label | Counted unit | Representative terms | Numeric context | Slide use | Caveat |
|---|---|---|---|---|---|
| Historical standard GO BP ORA | Unique gene symbols; copy-collapsed | Formation of quadruple SL/U4/U5/U6 snRNP; mRNA trans splicing; spliceosomal tri-snRNP assembly | `query_size=23`; `intersection_size=8`; adjusted p about 0.00145-0.00158 | Method contrast: standard ORA surfaces RNA/splicing-like terms when repeated copies are collapsed. | Historical/brainstorming; not HPRCv2 community-family testing; not copy-number-aware. |
| Historical standard GO MF ORA | Unique gene symbols; copy-collapsed | U4 snRNA binding; snRNA binding; olfactory receptor activity | `query_size=18`; intersection 8 for snRNA terms and 4 for olfactory receptor activity; adjusted p 9.11e-05 to 0.00824 | Context for why OR biology appears in older GO output. | Duplicate source rows are present; not canonical proof. |
| Historical coding-only GO BP ORA | Deduplicated coding-only query | Sensory perception of smell and related chemical-stimulus terms | `query_size=9`; `intersection_size=3`; p = 0.04008 | Context for olfactory biology in a small coding-only list. | No adjusted p column; small query; not copy-number-aware. |
| Historical coding-only GO MF ORA | Deduplicated coding-only query | Olfactory receptor activity; GTP binding; GTPase activity; G protein-coupled receptor activity | `query_size=9`; intersection 3 for OR/GTP terms; p = 0.02931 to 0.03903 | Context for OR and PAR1/C5 gene anchors. | No adjusted p column; not canonical HPRCv2 community analysis. |
| Original copy-weighted ORA | Gene copies against copy background | Regulation of transcription; immune response; cytokine activity; GTP binding; GTPase activity | Fold enrichment = 1 and adjusted p = 1 for listed rows | Mention only to show an abandoned/older copy-weighted attempt. | Exploratory; not useful as a biology claim. |
| Improved copy-weighted ORA exploratory rows | Gene copies against genome-wide background | Sensory perception of smell; olfactory receptor activity; GTP binding; GTPase activity; regulation of transcription | Fold enrichment 598.166 for smell/olfactory, 309.396 for GTP, 928.188 for transcription | Use only as caveated motivation for copy-aware methods. | Exploratory and statistically caveated; smell/olfactory rows list IL9R/IL9RP genes, not OR4F; validation flagged anti-conservative behavior and failed FDR control. |

Source rows:

- `paper_prep/_brainstorming/phr_GO_BP_enrichment.csv` rows 2-6.
- `paper_prep/_brainstorming/phr_GO_MF_enrichment.csv` rows 2-4.
- `paper_prep/_brainstorming/phr_coding_only_GO_BP_enrichment.csv` rows 2-7.
- `paper_prep/_brainstorming/phr_coding_only_GO_MF_enrichment.csv` rows 2, 4, 8, and 9.
- `paper_prep/_brainstorming/phr_copy_weighted_enrichment.csv` rows 2-6.
- `paper_prep/_brainstorming/improved_copy_weighted_enrichment.csv` rows 2-6.
- Caveat: `paper_prep/_brainstorming/copy_number_weighted_ora_methodology_synthesis.md:13-19`, `:67-104`, and `:388-418`.

Caption suggestion:

> Older GO/ORA tables are useful context, not final evidence. Standard/coding-only ORA collapses copy number; exploratory copy-weighted rows are statistically caveated. The report-backed slide should lead with HPRCv2 community-arm gene-family support.
