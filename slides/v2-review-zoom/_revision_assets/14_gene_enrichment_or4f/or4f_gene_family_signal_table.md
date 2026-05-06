# Candidate slide table: OR4F and high-copy subtelomeric gene families

| Signal | Slide-ready metric | Caveat |
|---|---|---|
| OR4F decay gradient | 5,023 OR4F annotations across 16 arms; 62.1% pseudogene overall | Use as canonical slide metric; source is HPRCv2 or4f_pseudogene_fraction.csv |
| OR family community signal | 10 OR4F genes in 7 communities; OR4F5 on 14 arms; OR4F8P on 14 arms | Presence pattern, not BH-significant enrichment: C3 OR q=0.071 |
| High-copy coding families | Recovered ORA table: OR4F=72, DUX4/FRG2/FRG2B=54, IL9R/IL9RP=58 copies | Historical Erik copy-aware ORA; parked under paper_prep/_brainstorming, use only as context |
| D4Z4 / DUX4L community | C1 chr4q/chr10q carries D4Z4; section 9 reports 22 DUX4L pseudogenes specific to C1 | Use as linked biological comparator, not as an OR4F substitute |
| Subtelomeric duplicon backbone | RPL23AP45 10 communities / 21 arms; SEPTIN14P22 9 communities / 22 arms; DDX11L16 9 communities / 20 arms | Canonical community-enrichment qualitative signal; Fisher tests do not survive BH |

Recommended one-line slide claim:

> OR4F is the clean visual: 5,023 annotations across 16 subtelomeric arms form a 11.1% to 99.8% pseudogenization gradient, while the broader gene-family analysis shows the same PHRs are dominated by copy-rich duplicon backbones rather than statistically significant one-family enrichments after BH correction.

Source note: the OR4F gradient comes from `/moosefs/guarracino/HPRCv2/PHR_III/plots/or4f_pseudogene_fraction.csv`; the old copy-aware ORA numbers are recovered from `paper_prep/_brainstorming/gene_copy_summary.csv` and should be treated as contextual because that directory is explicitly parked as noncanonical.
