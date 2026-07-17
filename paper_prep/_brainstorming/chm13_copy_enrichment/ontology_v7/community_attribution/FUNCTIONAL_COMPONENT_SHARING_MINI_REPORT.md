# Functional component sharing across subtelomeric arms

This note is an interpretive layer over the V7 copy-number-aware ontology enrichment and the arm-level projection in `ARM_FUNCTIONAL_CLASS_PROJECTION.tsv`. It does not add a new enrichment test. The statistical calls remain the exact V7 term rows; this report asks how the copy-number-bearing components are arranged across chromosome ends and PHR communities.

## General pattern

The subtelomeric system looks like a mosaic, not a set of chromosome ends that are uniformly homologous from telomere inward. Individual arms carry mixtures of shared components. Some components are large and concentrated, while others appear as smaller, repeated markers spread across several otherwise distinct communities. This explains why similarity clustering can place arms into broad neighborhoods while pairwise Jaccard similarity remains modest: arms can share one or a few blocks without sharing most of their PHR content.

The useful mental model is three nested layers:

- **PHR communities** are sequence-neighborhoods among chromosome arms.
- **Gene/source components** are copy-number-bearing pieces within those neighborhoods.
- **Ontology display classes** summarize the functional annotations carried by specific source components, without treating redundant ontology rows as independent biology.

In this view, communities are not pure functional compartments. They are assortments of components. The same component can bridge several communities, and one community can carry several unrelated components.

## Enriched display components by chromosome arm

| Component | Copy burden and arm placement | Community pattern | Interpretation |
|---|---|---|---|
| DUX4 / ZGA / transcription / nuclear envelope / cell cycle | 4q C1: 33; 10q C1: 32 | Restricted to C1 | A high-copy, community-private cassette. This is the clearest case where one community is dominated by one component. |
| WASH / endosomal actin / exocyst | 3q C3: 1; 9p C5: 1; 11p C3: 1; 12p C5: 1; 15q C8: 1; 16p C9: 1; 19p C3: 1; 20p C12: 1; 20q C5: 1; Xq C14: 1 | Distributed across C3, C5, C8, C9, C12, C14 | A small repeated component that crosses community boundaries. |
| DDX11 / helicase / chromosome | 3q C3: 1; 9p C5: 1; 11p C3: 1; 12p C5: 1; 15q C8: 1; 16p C9: 1; 19p C3: 1; 20p C12: 1; 20q C5: 1; Xq C14: 1 | Distributed across C3, C5, C8, C9, C12, C14 | Coextensive with the WASH component in the CHM13 projection, suggesting a recurrent WASH/DDX11-bearing subtelomeric module. |
| SEPTIN14 / septin / cytokinesis | 1p C11: 2; 2q C12: 2; 3q C3: 2; 5q C11: 2; 7p C3: 2; 9q C3: 2; 11p C3: 2; 16q C3: 2; 19p C3: 2; 20q C5: 2 | Mainly C3, with bridges into C11, C12, and C5 | A repeated two-copy pattern across several arms. This looks like a component that helps connect the C3 and C11/C12/C5 neighborhoods without making those arms globally identical. |
| WBP1L / CXCL12 signaling | 1p C11: 2; 3q C3: 2; 5q C11: 2; 6q C11: 2; 11p C3: 2; 15q C8: 2; 19p C3: 2 | C11 and C3, plus C8 | A paired-copy component that specifically links the C11 and C3 neighborhoods, with an additional copy pair on 15q. |
| RPL23A / ribosomal / nucleolar | 1q C6: 1; 2q C12: 1; 4q C1: 1; 5q C11: 1; 9q C3: 1; 10q C1: 1; 16q C3: 1; 17q C6: 1; 19q C6: 1; 20p C12: 1; 21q C6: 1 | Scattered across C1, C3, C6, C11, C12 | A small, broadly dispersed marker component. It contributes functional signal but is not the basis of one clean community. |

## Source-level components that help explain the mosaic

Not all recurrent sources define one of the six display classes. Some still help explain why arms share pieces without becoming highly similar overall.

| Source component | Arm placement | Community pattern | Note |
|---|---|---|---|
| CIC | 1p C11; 2q C12; 3q C3; 5q C11; 7p C3; 11p C3; 16q C3; 19p C3; 20q C5 | Mostly C3/C11, with C12/C5 | Tracks part of the same broad C3-C11-C12-C5 subtelomeric mixture, but is not a display class here. |
| GTF2I | 1p C11; 3q C3; 5q C11; 7p C3; 9q C3; 11p C3; 16q C3; 19p C3 | C3 and C11 | Another repeated C3/C11-associated component, helping explain partial sharing between these communities. |
| TUBB8/TUBB8B | 10p C2; 18p C2 | C2 | A compact 10p/18p component, consistent with C2 being a small two-arm community rather than a broad distributed module. |
| SNX18 | 21p C7 | C7 | A small acrocentric-p-arm signal in the C7 acrocentric community. |
| OR4F/OR4G-like copies | 1p C11; 3q C3; 5q C11; 11p C3; 19p C3 | C3 and C11 | Sparse olfactory-receptor-family copies are present, but broad olfactory ontology terms are not headline-supported V7 enrichments. |
| Xp/Yp-associated sources such as SHOX, PPP2R3B, PLCXD1, GTPBP6, FABP5 | Xp C15 in the CHM13 coordinate projection | C15 | These are supported copies but not part of the six enriched display classes. In this CHM13-centered analysis, Yp/Yq have no ontology-eligible PHR copies. |
| Xq-associated sources such as TRPC6, ELOC, DPH3, VAMP7, AMD1 | Xq C14 in the CHM13 coordinate projection | C14 | Xq carries both the WASH/DDX11 component and additional supported source-level copies; Yq has no ontology-eligible PHR copy in this CHM13 projection. |

## Community-level intuition

- **C1 (4q/10q)** is the strongest single-component case. It is mostly a DUX4/D4Z4-associated system with a few additional small source components. This is the best example of high copy number driving a clear functional signal.
- **C3 (3q/7p/9q/11p/16q/19p)** is a multi-component hub. It contains SEPTIN14, WBP1L, WASH, DDX11, RPL23A, CIC, GTF2I, and sparse OR4F/OR4G-like copies. It is less a single block than a recurring mixture of subtelomeric components.
- **C5 (6p/9p/12p/20q)** shares the WASH/DDX11 module and touches the SEPTIN14 component, but not the full C3/C11 mixture.
- **C11 (1p/5q/6q/8p)** is tied to C3 mainly through WBP1L, SEPTIN14, CIC, and GTF2I-like components. Its arms do not need high whole-region Jaccard similarity to share these source components.
- **C6 (1q/13q/17q/19q/21q/22q)** is mostly represented here by scattered RPL23A copies, not by the larger WASH/DDX11 or SEPTIN14/WBP1L modules.
- **C7 (13p/14p/15p/21p/22p)** is the acrocentric p-arm community, but this non-acrocentric-focused display only picks up a small SNX18-supported signal on 21p and no headline display class.
- **C14/C15 (sex-chromosome communities)** show asymmetry in this CHM13-coordinate projection: Xq has WASH/DDX11 plus other supported sources, Xp has supported source copies outside the six display classes, while Yp/Yq have no ontology-eligible PHR copies here.

## Working interpretation

The repeated components suggest a system of subtelomeric exchange in which particular cassettes, not whole chromosome ends, recur across different arms. Some cassettes are community-private and high-copy, like DUX4 in C1. Others are small and widely distributed, like WASH/DDX11 or RPL23A. Still others bridge selected neighborhoods, like SEPTIN14 and WBP1L between C3 and C11/C12/C5/C8.

This gives a coherent explanation for the mixed similarity structure: the PHR compartment is organized by recurrent, partially overlapping components. Communities capture preferred combinations of those components, while low Jaccard similarity reflects that most arms carry only a subset of the shared pieces plus arm-specific sequence.
