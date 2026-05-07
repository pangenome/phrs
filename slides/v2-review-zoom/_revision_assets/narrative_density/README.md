# Narrative Density Plan for Review Zoom Deck

Task: `review-zoom-density-methods-narrative`
Date: 2026-05-06 UTC
Scope: narrative and integration plan only. This file does not edit the Typst
deck source.

## Executive Recommendation

Treat the current `slides/v2-review-zoom` deck as a figure-review asset pool,
not as the linear talk deck. The current deck has 26 pages because agent-951
split composite figures into readable zoom pages. That was the right review
mode, but the next revision should route only one page from each analytical
cluster into the core talk path and move the rest to backup.

Recommended structure:

1. Core talk path: 14-15 pages, with four tiny methodology cards.
2. Backup/revision path: all dense raw crops, legacy labels, alternate tree/MDS
   views, deep 3D controls, and detailed gene-background slidelets.
3. Speaker notes carry method caveats; slides carry only one claim each.

The methodology should not be repeated on every result slide. Add short
interstitials at the four places where the audience changes analytical mode:

1. Sequence similarity/community detection before the heatmap/tree/MDS block.
2. MDS/population dispersion before any superpopulation claim.
3. 3D contact validation before the Hi-C/Pore-C/Dip-C block.
4. Gene annotation/enrichment before the DUX4/OR4F/TAR1 block.

## Source Evidence Used

- Provenance audit: the current review zoom deck is the agent-951 render
  (`10bee88`, squash `4862ec7`), while the earlier zoom deck came from
  agent-878 and used stale absolute worktree paths. See
  `slides/v2-review-zoom/_revision_assets/git_provenance/README.md:13`,
  `:17`, `:32`, and `:92`.
- Provenance audit: most current zoom assets are copied/cropped artifacts, not
  reproducible new renders. The table flags specific crop/script risks for
  slides 03, 04, 06, 07, 08, 09, 10, 11, 12, 13, and 14. See
  `slides/v2-review-zoom/_revision_assets/git_provenance/README.md:105`,
  `:106`, `:107`, `:109`, `:111`, `:113`, `:114`, `:115`, `:118`,
  `:120`, `:121`, `:124`, and `:125`.
- Current deck map: `zoom_review_deck.typ` currently lists 26 review pages,
  including split pages 03a/03b, 06a/06b, 07a/07b, 08a/08b,
  09a/09b/09c, 10a/10b, 12a/12b, 13a/13b, and 14a/14b/14c. See
  `slides/v2-review-zoom/_typst/zoom_review_deck.typ:109`,
  `:122`, `:131`, `:140`, `:149`, `:158`, `:167`, `:179`, `:188`,
  `:197`, `:206`, `:215`, `:224`, `:233`, `:242`, `:251`, `:260`,
  `:269`, `:278`, `:287`, `:296`, `:305`, `:314`, `:323`, `:332`,
  and `:341`.
- v2 narrative plan: the talk-level arc is methods -> empirical observation
  -> mechanism -> proof, with slide 14 explicitly compressible if time runs
  tight. See `slides/v2/SLIDES_v2_PLAN.md:56` and `:130`.
- v2 narrative plan: the MDS/PCA terminology drift is already identified and
  the recommended low-friction talk choice is MDS/PCoA relabeling unless a real
  PCA is generated. See `slides/v2/SLIDES_v2_PLAN.md:120` and
  `slides/v2-review-zoom/_revision_assets/git_provenance/README.md:212`.
- Fanout topics: Erik's numbered feedback is represented by the active
  review-zoom tasks. The relevant clusters are ER/connectivity, HPRCv2
  karyogram, length/clades, heatmap/tree labels, tree completeness, MDS
  dispersion, MDS/community consolidation, Hi-C methods, Hi-C visual redesign,
  and gene background/enrichment. See `.wg/graph.jsonl:343`, `:460`,
  `:441`, `:475`, `:799`, `:917`, `:929`, `:451`, `:807`, `:137`,
  `:581`, and `:349`.

## Narrative Rules for the Fan-in Task

1. One slide, one job. A slide can show a visual or define a method, but should
   not also carry caveats, backup statistics, and source history.
2. Core slides should be the talk. Backup slides should answer review/Q&A
   questions.
3. Do not keep a page in core only because a fanout task produced an asset.
   Every fanout output should be either core, backup, or provenance-only.
4. Use "MDS/PCoA" consistently unless a new true PCA asset is produced.
5. Avoid "canonical review asset" wording in audience-facing titles. Provenance
   belongs in `REVISION_NOTES.md`, not in slide titles.
6. Do not reintroduce agent-878 absolute paths. Any revised crop should have a
   recorded source path or generation script.

## Revised Core Talk Path

Target: 14-15 talk pages. Use "optional" rows only if timing allows or if the
fanout output is materially stronger than the current asset.

| Core order | Source / revision source | Proposed title | One-line takeaway | Speaker-note prompt | Integration notes |
|---:|---|---|---|---|---|
| 01 | Current 01 | Concerted evolution and unorthodox recombination of human subtelomeres | HPRC v2 makes population-scale subtelomeric exchange visible. | "This is a companion-scale view of chromosome ends, not a single-reference subtelomere catalog." | Keep title slide light. Do not add method bullets here. |
| 02 | Collapse current 02/03a/03b plus ER fanout | Method: an implicit graph over chromosome ends | Pairwise alignments are the graph; IMPG queries walk the interval forest by transitive closure. | "No GFA is built. The PAF edges are enough; interval trees make them queryable." | Use one compact workflow visual with a small ER inset. Move standalone interval-tree and full ER pages to backup. |
| 03 | 04 karyogram fanout or current 04/05 | Inter-chromosomal matches concentrate at subtelomeres | The genome-wide signal is end-biased and PAR-scale, not a genome-wide alignment artifact. | "Once we can query without chromosomal partitioning, the answer is concentrated at chromosome ends." | Prefer HPRCv2 karyogram if fanout recommends it; otherwise current 04 followed by 05 as backup. |
| 04 | 06 length redesign fanout | PHRs are long, uneven, and clade-structured | PHR length is not random noise; outlier arms map to named biological communities. | "This slide replaces the four-quadrant density page with a talk-readable summary of scale and outliers." | Core should be a ranked or community summary, not current dense 06a. |
| 05 | New interstitial | How sequence sharing becomes communities | Jaccard similarity, Leiden/community detection, and trees are three views of the same sequence-sharing structure. | "Before the heatmap: all points and trees in the next block come from the same PHR similarity matrix." | Tiny method card. 20-30 seconds. No equations beyond "Jaccard distance -> community/tree." |
| 06 | 07a heatmap/tree/p/q fanout | Arm-level similarity recovers named clades | The arm-by-arm map names the clades that drive the rest of the talk. | "Name only the recognisable groups: PAR1/PAR2, acrocentric p, DUX4/D4Z4, 10p-18p, major q-arm group." | Use improved 07a if available. Tree-left/p-q labels help. Keep 07b as backup unless it is clearer than the heatmap. |
| 07 | 09 MDS/community fanout; current 09c leading | Sequence-level MDS shows the same communities | The 41-arm structure persists when unfolded to individual PHR sequences. | "This is not a separate result from the heatmap; it is the same similarity system viewed at sequence level." | Lead with 09c. Integrate 09b labels into the plot if possible. Retire 09a from core. |
| 08 | New interstitial | How to read MDS and population spread | MDS is a distance display; population dispersion is a summary statistic, not an ancestry axis. | "The unit, sample imbalance, and dimensionality matter. We are quantifying spread, not claiming a clean population classifier." | Put before 08b only if a population-dispersion slide remains core. Otherwise make it backup speaker note. |
| 09 | Optional: 08b dispersion fanout | Population dispersion on the PHR MDS | If retained, AFR-wide spread is quantified with one simple metric and stated limitations. | "Say the metric once, then the caveat: unequal sample sizes and MDS dimensions limit interpretation." | Optional core. Drop first if time is tight or if the metric is not yet defensible. |
| 10 | New interstitial | How sequence communities are tested in 3D | Independent contact maps ask whether sequence-community arms contact each other more than non-community arms. | "The next plots are validation, not another way to call communities." | Tiny method card before 10a. Define W/B and Mantel in plain language. |
| 11 | Hi-C visual redesign fanout; current 10a | Sequence communities co-localize in 3D | The contact matrix is the visual proof that sequence-similar arms meet in 3D. | "Rows and columns are arms ordered by sequence community; warm blocks mean stronger within-community contact." | Use a square 1:1 contact matrix. Keep labels and color scale legible. |
| 12 | Hi-C redesign fanout; current 10b | The 3D signal survives confound exclusions | Removing acrocentric and sex-arm drivers strengthens rather than destroys the sequence-contact relationship. | "Above y=x means the association got stronger after exclusion; this is the confound defense." | Core if slide is simplified. Otherwise make it backup and state one robustness sentence on slide 11. |
| 13 | Current 13a/b, preferably consolidated | Pedigrees show exchange events directly | Family-level exchange events fall where the population graph predicts. | "This is the proof slide: population communities predict actual inherited inter-chromosomal patches." | Use one readable pedigree lead, not both top and bottom halves in core. Put detailed halves in backup. |
| 14 | New gene method/interstitial plus 14 fanouts | Gene cargo is counted copy-aware | DUX4, OR4F, and TAR1 are gene/repeat cargo riding the PHR exchange network; count copies, not just symbols. | "This is a catalog/enrichment slide, not a claim that these genes cause the exchange." | Optional if time is tight. If kept, make it one consolidated slide, not three slidelets. |
| 15 | Current 15 plus optional ED8 loop | Closing: sequence sharing, 3D proximity, exchange | Subtelomeres concertedly evolve through a feedback loop linking similarity, contact, and exchange. | "End on the thesis sentence. Do not enumerate every clade again." | Consider using the unused ED8 feedback-loop asset as a backup visual, not a dense core figure. |

### Timing Cut

If this must land inside a 15-minute talk without rushing, drop core row 09
first and present population dispersion as a backup/Q&A slide. Drop row 14
second if the gene-background fanout is not compact. The v2 plan already flags
slide 14 as the cleanest timing cut (`slides/v2/SLIDES_v2_PLAN.md:130`).

## Backup / Revision Slide Grouping

These pages should remain in the revised zoom deck as backup/revision pages or
appendix pages, not in the default talk path.

| Backup group | Current pages / fanout sources | Backup title pattern | Why backup, not core |
|---|---|---|---|
| Method details | 02, 03a, 03b | Backup: interval-tree schematic; Backup: full IMPG workflow; Backup: ER threshold plot | Useful for methods Q&A, but three separate pages stall the opening. Core should collapse them into one method slide. |
| Genome-wide alternatives | Current 04, 05; HPRCv2 karyogram fanout | Backup: manuscript Fig 1 heatmap; Backup: HPRCv2 karyogram; Backup: interchromosomal similarity count view | Keep whichever genome-wide visual is not chosen for core. Do not show all in sequence. |
| Length density | Current 06a; revised 06 alternatives | Backup: full per-arm length-density facets | Current 06a is accurate but too dense for live narration. Use it when someone asks to inspect all arms. |
| Clade callout details | Current 06b; revised clade summary | Backup: named clade callouts | Core should state 3-5 named groups; backup can list all outliers and introvert arms. |
| Tree completeness | Current 07b; 07b fanout rooted/unrooted candidates | Backup: full NJ/tree layout | The tree is a validation/inspection view. Core should usually use the heatmap unless the new tree is clearer. |
| MDS variants | 08a, 08b, 09a, 09b, 09c alternates | Backup: chromosome-colored MDS; Backup: legacy PCA-labeled raster; Backup: clade legend | Avoid the current 08/09 duplication. One labeled community MDS is core; the rest supports Q&A. |
| 3D deep controls | 10b if not core, 11, 12a, 12b | Backup: Mantel exclusions; Backup: single-cell 3D; Backup: mouse meiotic Hi-C; Backup: stage trajectory | Strong material, but too many assays in the main path dilute the argument. Use one core 3D result plus one robustness check. |
| Pedigree detail | 13a, 13b | Backup: full pedigree ribbon top/bottom | The proof matters, but the split view is a review affordance. Core needs one guided image. |
| Gene background | 14a, 14b, 14c; 14 gene-background and OR4F enrichment fanouts | Backup: DUX4/D4Z4 context; Backup: OR4F gradient; Backup: TAR1 repeat context; Backup: copy-aware enrichment table | Core can have one "gene cargo" synthesis slide. Detailed gene/repeat explanations belong in backup. |
| Closing evidence | Current unused `s15_ed8.png` / ED8 feedback loop | Backup: mechanistic feedback-loop schematic | Useful if the closing text slide feels too abstract, but not necessary if the closing thesis is spoken cleanly. |

## Tiny Methodology Interstitial Specs

Each interstitial should look like a transition card, not a methods slide from a
paper. Maximum three visual elements, maximum one numeric fact unless it is the
point of the card.

| Insert before | Title | On-slide content | Takeaway | Speaker-note prompt |
|---|---|---|---|---|
| Before arm heatmap/tree/MDS block | How sequence sharing becomes communities | `PHR alignments -> Jaccard distance -> Leiden / tree / MDS` | All community visuals derive from one sequence-similarity matrix. | "Do not treat the heatmap, tree, and MDS as independent evidence; they are complementary views of the same sequence-sharing graph." |
| Before population-colored MDS or dispersion | How to read MDS and population spread | `MDS/PCoA of Jaccard distances`; `dispersion = spread around group centroid`; `limits: unequal N, sequence-level units` | MDS shows geometry of PHR similarity; dispersion quantifies spread without overclaiming ancestry. | "Say MDS/PCoA, not PCA, unless a true PCA asset is generated. Name the unit of analysis." |
| Before Hi-C/Pore-C block | How 3D contact validates communities | `same sequence communities`; `independent contact maps`; `W/B + Mantel + exclusions` | 3D contact is orthogonal validation of the sequence communities. | "The community labels are frozen from sequence. The contact data ask whether those labels predict physical proximity." |
| Before gene block | How gene cargo is counted | `intersect PHRs with annotation`; `count gene copies/families`; `compare to copy-aware background` | Gene cargo should be framed as copy-aware catalog/enrichment, not one-symbol ORA. | "Do not imply DUX4/OR4F/TAR1 cause the PHR structure. They are cargo and markers carried by repeated exchange." |

## Exact Slide Title Changes

These are title-level changes only. The fan-in task owns any Typst edits.

| Current page | Current title | Proposed title | Core / backup |
|---|---|---|---|
| 01 | Title - review focus page | Concerted evolution and unorthodox recombination of human subtelomeres | Core |
| 02 | Implicit interval tree | Backup: interval trees make alignments queryable | Backup |
| 03a | IMPG workflow | Method: an implicit graph over chromosome ends | Core, collapsed with 02/03b |
| 03b | Erdos-Renyi callout | Why 12% sampling is enough for genome-wide closure | Core inset or backup, depending on ER fanout |
| 04 | Fig 1 panel a - genome-wide heatmap | Inter-chromosomal matches concentrate at subtelomeres | Core or backup, depending on karyogram decision |
| 05 | Interchromosomal similarities | Backup: genome-wide count view of interchromosomal sharing | Backup unless 04 is replaced |
| 06a | Length distribution - full faceted split view | Backup: full per-arm PHR length distributions | Backup |
| 06b | Length distribution clade callouts | PHR length outliers define named clades | Core if redesigned |
| 07a | Fig 1 panel c - arm heatmap | Arm-level similarity recovers named clades | Core |
| 07b | NJ tree with named clades | Backup: tree confirms the same clades | Backup by default |
| 08a | MDS colored by chromosome | Backup: MDS/PCoA colored by chromosome | Backup |
| 08b | MDS colored by superpopulation | Population dispersion on the PHR MDS | Optional core |
| 09a | PCA communities | Backup: legacy PCA-labeled community raster | Backup or retire |
| 09b | Clade legend | Backup: named-clade legend | Backup, or integrate into 09c |
| 09c | Community-colored MDS | Sequence-level MDS shows the same communities | Core |
| 10a | Fig 3 panel a - Hi-C/Pore-C contact matrix | Sequence communities co-localize in 3D | Core |
| 10b | Mantel exclusions | Confound check: the 3D signal survives exclusions | Core if simplified; backup if dense |
| 11 | Fig 3 panel c - single-cell 3D | Backup: single-cell 3D confirms the direction | Backup |
| 12a | Fig 4 panel d - mouse meiotic Hi-C | Backup: meiotic 3D contact peaks at zygotene | Backup |
| 12b | Meiotic stage trajectory inset | Backup: stage trajectory for meiotic contact | Backup |
| 13a | Pedigree exchange figure - top half | Pedigree proof: exchange events fall in predicted communities | Core, ideally consolidated |
| 13b | Pedigree exchange figure - bottom half | Backup: detailed pedigree exchange events | Backup |
| 14a | Gene biology - DUX4/D4Z4 | Gene cargo: DUX4/D4Z4 in the PHR exchange network | Backup/detail |
| 14b | Gene biology - OR4F | Gene cargo: OR4F family copies mark subtelomeric exchange | Backup/detail or part of consolidated core gene slide |
| 14c | Gene biology - TAR1 | Gene cargo: TAR1 as subtelomeric repeat context | Backup/detail |
| 15 | Closing - review focus page | Closing: sequence sharing, 3D proximity, exchange | Core |

## Speaker-Note Prompts by Block

### Opening / Method

- Keep the method promise narrow: "I need one trick to survey all chromosome
  ends without picking a reference or partitioning by chromosome."
- Say "alignment edges are the graph" once, then move on.
- If the ER card is shown: "The 12% figure is the sampling density that reaches
  full alignment after wfmash prefiltering; it is far above the connectivity
  threshold, so closure is not a sparse disconnected walk."

### Empirical Signal

- Open slide 03 with: "Now that we can query the graph, what does it see?"
- Use PAR-scale language early, but do not list all clades yet.
- On length/clade slide, name only the groups that recur later: PARs,
  acrocentric p-arms, DUX4/D4Z4 4q/10q, 10p/18p, and the large q-arm group.

### Communities / MDS

- Before heatmap: "Same matrix, three views."
- On 07a: "This is the arm-level map."
- On 09c: "This is the individual-sequence map."
- If 08b remains core: "Population spread is a secondary readout; the primary
  claim is still community structure."

### 3D Validation

- Before 10a: "Now freeze the sequence communities and test them against
  independent 3D contact data."
- On 10a: "The picture is the matrix; the statistic quantifies the blockiness."
- On 10b: "The exclusion result is the defense against the obvious confound."
- Keep single-cell and mouse material as "we also tested this" backup unless
  the talk specifically needs a meiotic bridge before the pedigree slide.

### Pedigree / Proof

- Say "prediction" and "direct exchange" explicitly.
- Avoid reading all event counts on screen. The headline is that exchange
  events fall where population communities predicted.

### Gene Cargo / Closing

- For gene cargo, use "copy-aware catalog" language before "enrichment."
- Close on the loop: sequence sharing creates 3D opportunity; 3D opportunity
  enables exchange; exchange renews sequence sharing.

## Integration Checklist

- [ ] Start fan-in from `zoom_review_deck.typ` but route pages through the
      core/backup grouping above instead of preserving current 26-page order.
- [ ] Add `REVISION_NOTES.md` rows for every changed major asset, with source
      path, script, external URL, or commit. This is required because many
      current PNGs are copied/cropped assets with weak reproducibility.
- [ ] Do not copy any Typst snippets from `slides/v2-zoom/_typst/zoom_deck.typ`
      unless prior-agent absolute worktree paths are removed.
- [ ] Resolve MDS/PCA wording globally: slide titles, axis labels, footers, and
      speaker notes should say MDS/PCoA unless a real PCA artifact is generated.
- [ ] Choose exactly one core visual from each repeated cluster:
      genome-wide signal, length/clade summary, arm-level community map,
      sequence-level MDS, 3D contact, pedigree, gene cargo.
- [ ] Put current dense pages 06a, 09a/09b, 11, 12a/12b, 13b, and 14a/14b/14c
      into backup unless a fanout task produces a simpler replacement.
- [ ] Add four tiny method cards or equivalent interstitials:
      sequence/community, MDS/dispersion, 3D validation, gene cargo/enrichment.
- [ ] Keep each method card visually sparse: no more than three terms and one
      short speaker prompt.
- [ ] For slide 03/ER and slide 04/karyogram, record external repo commit or URL
      in revision notes; the provenance audit already identifies these as
      external hooks.
- [ ] For revised crop assets, include either a generation script, a source PDF
      plus crop geometry, or a clear "manual crop" note in revision notes.
- [ ] Render the revised deck as a new v2 PDF and PNG sequence, preserving the
      current zoom deck as reference until the fan-in task validates the new
      render.

## Non-goals

- This note does not choose final images from the still-running fanout tasks.
- This note does not edit `slides/v2-review-zoom/_typst/zoom_review_deck.typ`.
- This note does not decide whether slide 14 belongs in the final timed talk;
  it provides a compact route if Erik wants the gene-background material.
