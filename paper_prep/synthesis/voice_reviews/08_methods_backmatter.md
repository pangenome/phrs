# Methods and back-matter consistency audit

## Verdict

Minor-to-moderate revise before implementation.

The Methods mostly contain the support needed for the current figure spine, but
they also preserve several reviewer-era or over-detailed analyses that pull the
paper away from the house abstract voice. The main implementation risk is not
missing technical support for the figures; it is mismatch between the compressed
BoG/Nature story and Methods/back-matter promises that imply a broader reviewer-
response manuscript.

The implementation pass should first decide which non-figured 3D controls remain
part of the current story. After that, Methods can be shortened and made more
consistent with the body. Do not polish text that should be removed.

## Keep / Remove / Relocate Table

| Lines | Current content | Recommendation | Reason |
| --- | --- | --- | --- |
| `submission/paper.tex:413-420` | Sample selection and reference frame | Keep, revise wording | Required for 466/465+CHM13 promise. `233 HPRC v2 v1.1 individuals` at line 415 is awkward; check intended HPRC v2 naming. |
| `submission/paper.tex:422-433` | Sample exclusions | Keep, possibly compress | Supports flank census and explains GRCh38/Y/PAR and chr18_q chimera exclusions. Technical but relevant. |
| `submission/paper.tex:435-440` | Telomere-anchored flank extraction | Keep | Required support for 18,827 flanks, 48 arms and 500 kb truncation caveat. |
| `submission/paper.tex:442-448` | wfmash all-vs-all alignment | Keep, tighten | Supports implicit graph method and 95% identity. Consider moving panmixia/transitive-closure rationale into the following subsection to avoid repetition. |
| `submission/paper.tex:450-458` | IMPG transitive closure and 11.6% sampling | Keep | Required for abstract's `approximately 12%` and no chromosomal partitioning. |
| `submission/paper.tex:460-465` | PHR detection | Keep | Required support for 15,668 PHRs and 41/48 arms. |
| `submission/paper.tex:467-473` | Pangenome graph and Jaccard similarity | Keep | Required for Fig. 2A and sequence-level similarity. |
| `submission/paper.tex:475-478` | Arm-level distance matrix | Keep | Required for arm-level community analyses. |
| `submission/paper.tex:480-487` | Leiden community detection | Keep | Required for 15 arm-level and 50 sequence-level communities. |
| `submission/paper.tex:489-492` | UPGMA dendrogram | Keep if Fig. 2 legend retains UPGMA | Required only if Fig. 2 continues to use UPGMA ordering/robustness language. |
| `submission/paper.tex:494-508` | Neighbour-joining tree and bootstrap | Keep but qualify | Supports ordering/known relationship recovery, but the low support for the six-q-arm group could distract from the abstract unless framed as an ordering-control limitation. |
| `submission/paper.tex:510-542` | Hi-C, Pore-C and CiFi pipeline | Keep, revise caution | Required for Fig. 4/ED1 and 3D evidence. Align with 3D memo: controls argue against a simple multi-mapping artifact, but do not estimate exact PHR-internal contact magnitude. |
| `submission/paper.tex:544-554` | Reproducibility from deposited data | Keep if 3D reanalysis remains central | Useful, but the final sentence is unusually direct. It fits Methods better than Results. |
| `submission/paper.tex:556-565` | Exclusion controls | Remove or relocate unless explicitly retained as Methods-only support | This reads like reviewer-era 3D control material and has no current figure. Keep only if body/legend explicitly needs it. |
| `submission/paper.tex:566-576` | Single-cell 3D and controls | Decide before editing | Body lines `submission/paper.tex:218-220` mention CiFi, Dip-C and sperm scHi-C, so some support is currently required. PBMC and `S_all` controls at lines 572-575 look like reviewer-era drift unless retained deliberately. |
| `submission/paper.tex:577-593` | Mouse pipeline | Keep | Required for Fig. 4c and abstract's mouse zygotene claim. |
| `submission/paper.tex:595-608` | Pedigree odgi-untangle | Keep, tune evidence language | Required for Fig. 5 and abstract's cautious pedigree sentence. |
| `submission/paper.tex:610-614` | Software versions | Keep | Necessary back matter. |
| `submission/paper.tex:616-625` | Data and code availability | Keep, revise after final scope | Align re-alignment scripts and deposited-file caveat with whichever 3D datasets remain. |
| `submission/paper.tex:627-640` | Limitations | Revise substantially | Contains reviewer-era/stale items. The FST CI at lines 636-638 is inconsistent with the instruction that popgen/FST was cut. |
| `submission/paper.tex:647-665` | ED1 legend | Keep | ED1 is current CHM13 Hi-C replicate and supports Fig. 4a. |

## Required Support Map

| Claim or figure promise | Body / legend location | Required Methods support |
| --- | --- | --- |
| 466 near-complete assemblies; 465 HPRC v2 haplotypes plus CHM13 | Abstract `submission/paper.tex:53`; Results `submission/paper.tex:111-113` | Sample selection `submission/paper.tex:413-420`; check wording at line 415. |
| Telomere-anchored 500 kb flanks and 18,827 flanks | Results `submission/paper.tex:111-113`; Fig. 1/PHR extent | Flank extraction `submission/paper.tex:435-440`; sample exclusions `submission/paper.tex:422-433`. |
| No chromosomal partitioning; implicit pangenome graph | Abstract `submission/paper.tex:53`; Results `submission/paper.tex:110-125` | wfmash/IMPG support `submission/paper.tex:442-458`. |
| 15,668 PHRs on 41/48 arms | Results `submission/paper.tex:129-134`; abstract `submission/paper.tex:53` | PHR detection `submission/paper.tex:460-465`. |
| 3.51 Mb outside acrocentric short arms and PARs | Results `submission/paper.tex:148-150`; Fig. 1 legend `submission/paper.tex:304-307` | Flank extraction and PHR detection `submission/paper.tex:435-465`; limitations should retain the 500 kb truncation caveat at `submission/paper.tex:633-635`. |
| Fig. 2 connected graph and arm-level similarity/community structure | Results `submission/paper.tex:155-173`; Fig. 2 legend `submission/paper.tex:311-332` | Pangenome/Jaccard, distance matrix, Leiden, UPGMA and NJ support `submission/paper.tex:467-508`. |
| Known systems plus linked q-arm group | Abstract `submission/paper.tex:53`; Results `submission/paper.tex:164-170` | Community methods `submission/paper.tex:480-508`; note low bootstrap support at `submission/paper.tex:503-507` must not be over-read as disproving the descriptive community. |
| Fig. 4 human Pore-C and CHM13 Hi-C association | Results `submission/paper.tex:211-217`; Fig. 4/ED1 legends `submission/paper.tex:361-386`, `submission/paper.tex:647-665` | 3D pipeline and pointwise Spearman methods `submission/paper.tex:510-542`; reproducibility `submission/paper.tex:544-554`. |
| Mouse zygotene bouquet peak | Abstract `submission/paper.tex:53`; Results `submission/paper.tex:205-210`; Fig. 4c legend `submission/paper.tex:379-386` | Mouse pipeline `submission/paper.tex:577-593`. |
| Other 3D datasets: CiFi, Dip-C, sperm scHi-C | Results `submission/paper.tex:218-220` | Single-cell/3D support `submission/paper.tex:566-576`; if body removes this sentence, much of this Methods subsection can go. |
| Pedigree patch enrichment and cautious exchange evidence | Abstract `submission/paper.tex:53`; Results `submission/paper.tex:229-243`; Fig. 5 legend `submission/paper.tex:390-406` | Pedigree methods `submission/paper.tex:595-608`. |
| Data/code reproducibility | Data availability `submission/paper.tex:616-625` | Ensure paths and `scripts/hic-realign/` claims match final retained analyses. |

## Drift Audit

- FST/popgen: `submission/paper.tex:636-638` mentions
  `$F_{\mathrm{ST}}$ block-jackknife CI`, but the body no longer carries a
  population-genetic result and the project instruction says popgen/FST was cut
  from body and Methods. Remove this sentence or narrow it to currently reported
  confidence intervals.

- CEPH1463 and RPE-1: not in Methods lines `408-667`, but Fig. 5 legend
  `submission/paper.tex:403-404` says cross-assembler validation and RPE-1 are
  described in the main text. The pedigree memo flags this as an internal
  consistency problem. Remove the legend promise unless those analyses are
  deliberately restored.

- Single-cell controls: `submission/paper.tex:566-576` mixes required support
  for body-mentioned Dip-C/sperm scHi-C with controls that feel reviewer-era
  (`PBMC Dip-C negative control`, `$S_{\mathrm{all}}$ negative control`). Decide
  whether the body keeps the other-contact-map sentence at
  `submission/paper.tex:218-220`; then cut or retain only the needed Methods.

- Exclusion controls: `submission/paper.tex:556-565` reports a multi-resolution
  Mantel exclusion walk with no current figure. This looks like reviewer-era
  control material. Remove unless the 3D section or response package explicitly
  needs it.

- Gene enrichment: the current Methods scope does not include gene enrichment,
  which is good. Do not reintroduce it in this manuscript pass.

- Within-community heterogeneity: no obvious Methods block remains, which is
  good. Do not reintroduce.

- 14-test 3D forest/controls: the exclusion-control and single-cell-control
  blocks may be remnants of this broader reviewer-era validation package. Keep
  the pieces required for Fig. 4/ED1; avoid retaining a control catalog that the
  body does not use.

- "Observed directly" pedigree language: Methods `submission/paper.tex:595-608`
  is mostly fine because it describes the untangle procedure and null. The
  stronger risk is in body/legend wording; keep Methods technical and let the
  Results use "evidence consistent with" for biological exchange.

## Integration Risks Before Implementation

1. Do not launch a broad prose edit until the 3D scope is decided. If the body
   keeps CiFi, Dip-C and sperm scHi-C, Methods need a concise support paragraph;
   if the body follows the tight figure-only story, much of `submission/paper.tex:566-576`
   should be removed.

2. The limitations paragraph is currently the highest-risk back-matter drift.
   It mixes valid limitations (somatic human Hi-C, no germline LAD/Lamin B1,
   12% sampling threshold, 500 kb flank truncation, small Hi-C N) with cut or
   stale material (`F_ST`, low-callability cM/Mb). Rewrite after the body scope
   is settled.

3. Keep the Methods support for multi-mapping cautious. The 3D review recommends
   language such as "supports" or "argues against a simple multi-mapping
   artifact" rather than "confirms"; Methods lines `516-526` should be aligned
   with that.

4. Align spelling globally during implementation. Methods still use
   `inter-chromosomal` at `submission/paper.tex:531`, `submission/paper.tex:534`
   and `submission/paper.tex:536`; the house style is `interchromosomal` in
   prose.

5. Figure legends should not promise analyses absent from the body. The Fig. 5
   legend issue at `submission/paper.tex:403-404` should be fixed in the same
   pass as pedigree/conclusion edits.

6. Preserve the single-file LaTeX constraint. These reviews identify edits, but
   implementation should remain scoped to `submission/paper.tex` and any
   synthesis notes intentionally updated to match.

## Implementation Priority

1. Resolve obvious inconsistencies: Fig. 5 CEPH/RPE-1 promise; FST limitation;
   `interchromosomal` spelling.
2. Decide 3D support scope: keep or cut non-figured CiFi/Dip-C/sperm/control
   material.
3. Rewrite the limitations paragraph after the above decisions.
4. Then perform the section-level voice edits from memos 01-07.
