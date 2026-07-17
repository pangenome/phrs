# V7 community functional attribution

## Interpretation

This is a **post-inference attribution** of the 209 exact midpoint `PRIMARY_SUPPORTED` V7 term rows. It introduces no ontology test, no functional-class p-value, and no change to any V7 support decision. The inferential object remains each exact `(collection, relation, term_id)` row.

Physical-copy burden is counted at coordinate resolution: N distinct coordinate copies contribute N units, even when they share a gene name, source gene, family, locus family, arm community, or display class. Ontology rows are redundant annotations and must not be interpreted as independent biological systems. The class labels below are post-inference display summaries only.

The primary universe contains 187 ontology-eligible midpoint PHR copies; 186 distinct copies contribute to at least one supported exact term. Any-overlap-supported rows are retained only as a labeled sensitivity layer in `COMMUNITY_EXACT_TERM_ATTRIBUTION.tsv` and are not used to define the headline classes.

## Clean community–function mappings

| Display class | Copy burden by community | Pattern | Representative exact supported terms |
|---|---|---|---|
| DUX4 / ZGA / transcription / nuclear envelope / cell cycle | C1: 65 | C1 carries 65/65 | negative regulation of G0 to G1 transition (GO:0070317, direct); Zygotic genome activation (ZGA) (R-HSA-9819196, direct); RNA polymerase II transcription regulatory region sequence-specific DNA binding (GO:0000977, direct) |
| WASH / endosomal actin / exocyst | C3: 3, C5: 3, C8: 1, C9: 1, C12: 1, C14: 1 | distributed across 6 communities | WASH complex (GO:0071203, direct); exocyst (GO:0000145, direct); Arp2/3 complex-mediated actin nucleation (GO:0034314, direct) |
| DDX11 / helicase / chromosome | C3: 3, C5: 3, C8: 1, C9: 1, C12: 1, C14: 1 | distributed across 6 communities | 5'-3' DNA helicase activity (GO:0043139, direct); G-quadruplex unwinding activity (GO:0160225, direct); establishment of sister chromatid cohesion (GO:0034085, direct) |
| RPL23A / ribosomal / nucleolar | C1: 2, C3: 2, C6: 4, C11: 1, C12: 2 | distributed across 5 communities | nucleolus (GO:0005730, direct); rRNA binding (GO:0019843, direct) |
| SEPTIN14 / septin / cytokinesis | C3: 12, C5: 2, C11: 4, C12: 2 | C3 carries 12/20 | septin complex (GO:0031105, direct); septin ring (GO:0005940, direct); cytoskeleton-dependent cytokinesis (GO:0061640, direct) |
| WBP1L / CXCL12 signaling | C3: 6, C8: 2, C11: 6 | distributed across 3 communities | CXCL12-activated CXCR4 signaling pathway (GO:0038160, direct); chemokine (C-X-C motif) ligand 12 signaling pathway (GO:0038146, ancestor); C-X-C chemokine receptor CXCR4 signaling pathway (GO:0038159, ancestor) |

C1 is the cleanest mapping: all DUX4-source contributors fall in C1, tying that community to exact ZGA, transcription-regulatory, nuclear-envelope, and cell-cycle signals. WASHC1 and DDX11 signatures are distributed across multiple communities rather than constituting a single-community system. RPL23A and SEPTIN14 likewise span communities, with their exact nucleolar/rRNA-binding and septin/cytokinesis terms retained as separate source-defined displays.

WBP1L-source copies robustly carry exact CXCL12/CXCR4 signaling rows. The supported exact term `protein localization to perinuclear region of cytoplasm` has a SEPTIN14 contributor signature, so it remains in the SEPTIN14 class rather than being conflated with WBP1L/CXCL12 signaling.

## Community synopsis

The repeated contributor-row burden below is supplied only as an audit of annotation incidence; because one copy can contribute to many related direct and ancestor terms, it is ontology-redundant. Unique-copy burden is the interpretable physical-copy summary.

| Community | Eligible PHR copies | Unique supported copies | Redundant supported term-copy rows | Dominant display class(es) | Contributing source symbols |
|---|---:|---:|---:|---|---|
| C1 | 75 | 75 | 4784 | DUX4 / ZGA / transcription / nuclear envelope / cell cycle | AGGF1, CLUH, DUX4, FRG2, FRG2B, RARRES2, RPL23A |
| C2 | 3 | 3 | 40 | no defined display class | IL9R, TUBB8, TUBB8B |
| C3 | 46 | 45 | 1509 | SEPTIN14 / septin / cytokinesis | CIC, DDX11, GTF2I, IL9R, MIR1302-2, OR4C13, OR4F17, OR4F21, OR4F5, RPL23A, SEPTIN14, TUBB, WASHC1, WBP1L |
| C4 | 0 | 0 | 0 | no defined display class | none |
| C5 | 10 | 10 | 518 | DDX11 / helicase / chromosome; WASH / endosomal actin / exocyst | CIC, DDX11, IQSEC3, SEPTIN14, SEPTIN14P20, WASHC1 |
| C6 | 4 | 4 | 88 | RPL23A / ribosomal / nucleolar | RPL23A |
| C7 | 2 | 2 | 32 | no defined display class | SNX18 |
| C8 | 4 | 4 | 165 | WBP1L / CXCL12 signaling | DDX11, WASHC1, WBP1L |
| C9 | 3 | 3 | 146 | DDX11 / helicase / chromosome; WASH / endosomal actin / exocyst | DDX11, IL9R, WASHC1 |
| C10 | 0 | 0 | 0 | no defined display class | none |
| C11 | 17 | 17 | 430 | WBP1L / CXCL12 signaling | CIC, GTF2I, OR4F29, OR4F3, RPL23A, SEPTIN14, WBP1L |
| C12 | 8 | 8 | 313 | RPL23A / ribosomal / nucleolar; SEPTIN14 / septin / cytokinesis | CIC, DDX11, MIR1302-2, RPL23A, SEPTIN14, WASHC1 |
| C13 | 0 | 0 | 0 | no defined display class | none |
| C14 | 9 | 9 | 220 | DDX11 / helicase / chromosome; WASH / endosomal actin / exocyst | AMD1, DDX11, DPH3, ELOC, IL9R, TRPC6, VAMP7, WASHC1 |
| C15 | 6 | 6 | 95 | no defined display class | FABP5, GTPBP6, KRT18, PLCXD1, PPP2R3B, SHOX |

## Boundary and negative findings

No broad immune or olfactory term is a headline-supported exact midpoint enrichment in V7. 
Source and class labels describe annotation-bearing CHM13 physical copies, not expression, protein activity, dosage, retained pseudogene function, biological independence of neighboring copies, or population prevalence. The any-overlap layer is a sensitivity analysis and should not replace the midpoint headline.

## Reproducibility

Run `python3 build_community_attribution.py` from any directory, then `python3 -m unittest -v test_community_attribution.py`. JSON arrays are used for copy IDs because the atomic V7 copy identifiers themselves contain pipe characters. `OUTPUT_MANIFEST.sha256.tsv` records the exact files covered by release hashes.
