# Speaker Notes

## Overall Talk Track

This section is a short interpretation bridge, not a methods appendix. The
message is that subtelomeric PHRs are a copy-rich biological system, so unique
gene-symbol ORA can understate families that have expanded through duplicated
sequence. The production v5 claims should come from HPRCv2 arm/community
support tables and the OR4F gradient. The older genome-wide copy-weighted ORA
tables can motivate the method, but they should not be treated as validated
statistical results.

## Slide 1 Notes

Start by contrasting questions. Standard ORA asks, "is this gene symbol in my
set?" Copy-aware enrichment asks, "how much repeated PHR sequence is occupied
by this family?" That distinction matters here because the region being tested
is enriched for duplicated subtelomeric families, pseudogenes, lncRNAs, and
retrogenes. A family copied across many arms can look modest in a unique-symbol
list while dominating the actual annotation landscape.

Useful phrasing:

"For these regions, copy number is not a nuisance variable. It is part of the
biology we are trying to summarize."

## Slide 2 Notes

OR4F is the cleanest visual entry point. OR4F genes are olfactory receptor
family members frequently found near chromosome ends. The HPRCv2 OR4F table
shows 5,023 annotations across 16 arms, split between coding and pseudogene
copies. The wide pseudogene-fraction range, from chr7p to chr15q, is the visual
evidence for arm-specific turnover.

Be precise about significance. The OR family is broadly present across
communities, with OR4F5 and OR4F8P each found on 14 arms, but the HPRCv2
community-family Fisher test does not make this a BH-significant enrichment
claim. The C3 OR row is a suggestive qualitative signal with q ~= 0.071.

Avoid saying:

- "OR4F is significantly enriched after correction."
- "The copy-weighted ORA proves olfactory enrichment."

The parked copy-weighted ORA table has a smell/olfactory row, but the current
inventory warns that the improved script maps that row to IL9R/IL9RP genes
rather than OR4F. Keep the OR4F claim anchored to the HPRCv2 OR4F gradient and
community support.

## Slide 3 Notes

This slide broadens the biology beyond OR4F. RPL23A pseudogene copies,
DDX11L, WASH, and FAM138 are better described as a recurrent subtelomeric
duplicon backbone than as one protein-coding pathway. RPL23AP45, DDX11L16,
FAM138D, and WASH6P are useful examples because they recur across many
communities and arms in the HPRCv2 tables.

Use brief background, not functional overreach. In these annotation tables,
DDX11L and FAM138 behave as subtelomeric noncoding families. WASH and RPL23A
entries are mostly duplicated retrogene/pseudogene-style annotations. The
slide claim is repeated subtelomeric architecture, not a single active pathway.

For GTP-related language, stay at the anchor-gene level. GTP binding and
GTPase activity are plausible biological labels because the old copy-weighted
ORA rows are driven by GTPBP6 and IQSEC3. GTPBP6 is a PAR1 GTP-binding gene.
IQSEC3 is an ARF guanine nucleotide-exchange factor, so it is related to
GTPase regulation rather than being the same kind of direct GTP-binding signal.
Use this as a hypothesis or biological hook, not a validated GO result.

Suggested phrasing:

"The GTP signal is interesting because it points back to named genes we can
locate in the community tables, especially GTPBP6 in PAR1 and IQSEC3 in C5.
But the 309x copy-weighted ORA number is still a reanalysis target, not a
production statistic."

## Slide 4 Notes

This is the honesty slide. It should prevent three overclaims.

First, community assignment is not the same thing as a called CHM13 interval.
The HPRCv2 community tables assign arms in the pangenome analysis, while
`chm13.phrs.bed` is a CHM13 reference interval extract. Community-assigned
arms without called CHM13 PHR rows include C5 `chr6_p`, C7 `chr13_p`, C14
`chrY_q`, and C15 `chrY_p`. If a downstream renderer draws CHM13 intervals,
those arms must either be omitted from that interval-specific graphic or
captioned as community arms without called CHM13 interval rows.

Second, the HPRCv2 Fisher screen is a useful statistical check but not a
positive corrected enrichment result. There were 116 community-family tests and
no BH-significant rows. The right wording is "presence pattern", "copy/support
signal", or "duplicon architecture", not "significant enrichment" unless a
future validated analysis supports that exact claim.

Third, the older weighted-hypergeometric ORA stream is exploratory. Prior
validation found non-uniform null p-values, inflated Type I error under
gene-level sampling, and failed BH/FDR behavior. Before using its p-values or
large fold enrichments as slide truth, rerun the analysis with an empirical
null, permutation test, or another calibrated copy-number-aware method; check
the weighted parameter scaling as part of that validation; and show the
standard-vs-copy-weighted contrast side by side.

TAR1 can be mentioned only if the presenter wants a repeat-landscape bridge.
It is useful context because TAR1 is telomere-proximal and C2 has high TAR1
density, but TAR1 is not a gene family and should not be mixed into the
gene-enrichment claim.

## Safe Closing Line

"The take-home is architectural: copy-rich subtelomeric regions carry repeated
gene-family signals that a unique-symbol ORA can flatten. For v5, we show the
table-backed copy/support patterns and keep formal copy-weighted ORA as an
explicitly caveated reanalysis target."
