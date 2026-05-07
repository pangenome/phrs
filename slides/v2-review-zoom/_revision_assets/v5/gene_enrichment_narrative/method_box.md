# Method Box

## Slide-Ready Box

**Copy-number-aware enrichment**

Standard ORA asks whether each gene name is present. That is useful, but it
collapses a copied subtelomeric family to one symbol.

For PHRs, the copy-aware question is different:

- Which gene families recur across duplicated PHR sequence?
- How many gene instances, arms, and HPR communities carry the family?
- Does the signal remain visible when we compare standard unique-symbol
  results with copy/support-aware counts?

In v5, the main evidence is HPRCv2 arm-level community support plus the OR4F
per-arm copy gradient. Older genome-wide copy-weighted ORA results are
exploratory until rerun with calibrated permutation or empirical-null
validation. Do not transfer their p=0 rows or large fold-enrichment values to
the deck without checking the weighted parameters and FDR behavior.

## One-Sentence Speaker Version

"Standard ORA treats OR4F5 as one gene name; a copy-aware view asks whether
OR4F5 and related families are repeatedly occupying subtelomeric sequence
across arms and communities."

## Wording Guardrail

Say "community-arm support" or "HPRCv2 pangenome support" when discussing
families across arms. Do not imply that every community arm has a called CHM13
PHR interval; C5 `chr6_p`, C7 `chr13_p`, C14 `chrY_q`, and C15 `chrY_p` are
community-assigned arms without rows in `chm13.phrs.bed`.
