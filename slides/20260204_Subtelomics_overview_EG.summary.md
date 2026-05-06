Source: 

# Slide 1 — Title

**Title:** Inter-chromosomal subtelomeric relationships

**Visible text:**

* Andrea Guarracino, Erik Garrison
* Update — 2026/02/04

**Description:**
A clean title slide introducing a talk/update about inter-chromosomal relationships in subtelomeric regions.

---

# Slide 2 — Implicit interval tree

**Title:** Implicit interval tree

**Visible text:**

* Li and Rong, 2020
* Alignment information
* CIGAR string
* Target range
* explicit interval tree
* implicit interval tree
* Node Legend
* `[Start, End) Index, MaxEnd`
* Level 0, Level 1, Level 2, Level 3

**Interval list shown:**

```text
100 150
130 200
170 300
180 250
200 250
250 350
270 300
300 320
350 450
390 420
```

**Description:**
The slide explains how genomic alignment intervals can be represented as an implicit interval tree. At the top is a horizontal interval plot with ranges labeled by target position. Below it, an explicit tree diagram shows interval nodes with metadata such as start/end, index, and max end. A highlighted node corresponds to target range `[300, 320)`, with an arrow indicating associated alignment information. On the right, the same intervals are shown as a compact ordered list representing the implicit interval tree.

---

# Slide 3 — IMPG workflow

**Title:** IMPG: IMplicit Pangenome Graphs

**Visible text:**

* Alignment `[start, end)` → intervals → implicit interval tree `(per target)`
* `wfmash`
* `https://github.com/waveygang/wfmash`
* `https://github.com/pangenome/impg`
* We perform the alignment.
* We index the alignment.
* We create an interval forest, with an interval tree for each sequence.
* sequence 1
* sequence 2
* sequence 3
* sequence 4
* sequence 5
* sequence N

**Description:**
This slide shows the IMPG pipeline. On the left is a wave image representing `wfmash`, with a GitHub link underneath. In the center is a dotplot-like alignment matrix, labeled as the alignment indexing step. On the right is a stack of interval-tree diagrams representing an interval forest, where each sequence has its own interval tree. The workflow proceeds from alignment intervals to indexed implicit interval trees per target sequence.

---

# Slide 4 — Querying the pangenome

**Title:** Querying the pangenome — HPRCv2

**Visible text:**

* Genome-wide alignment identity patterns `(100kb windows)`
* Average identity per matching chromosome `(range: 0.980–1.000)`
* Avg Identity
* 1.000, 0.995, 0.990, 0.985, 0.980
* Target chromosome
* Position `(Mbp)`
* chr18 q-arm `(subtelomeric)`
* Subtelomeric regions feature patterns that have been described but lack precise quantification.
* Slide number: 4

**Description:**
The main visual is a genome-wide plot of average alignment identity across chromosomes in 100 kb windows. Chromosomes are listed on the y-axis and genomic position in megabases on the x-axis. A red-to-black color scale indicates average identity, ranging from 0.980 to 1.000. An inset zooms into the subtelomeric q-arm of chromosome 18, showing dense matching patterns near the chromosome end. The right side emphasizes the motivation: subtelomeric patterns are known qualitatively but need more precise quantification.

---

# Slide 5 — Interchromosomal similarities

**Title:** Interchromosomal similarities — HPRCv2

**Visible text:**

* Number of unique chromosomes per region across CHM13 `(100kb windows)`
* Number of chromosomes
* Position `(Mbp)`
* Region

  * CEN
  * PAR
  * PHR
  * XTR

**Description:**
The slide contains a chromosome-by-chromosome plot across CHM13 using 100 kb windows. The x-axis is genomic position in megabases, and the y-axis tracks chromosomes/arms. Orange traces show the number of unique chromosomes matching each region. Colored background highlights identify annotated region types: CEN, PAR, PHR, and XTR. The plot highlights strong interchromosomal similarity signals, especially near subtelomeric and special repetitive regions.

---

# Slide 6 — Length distributions of interchromosomal matches

**Title:** Length of regions with inter-chromosomal matches

**Visible text:**

* Region length `(kb)`
* Count
* Facets labeled by chromosome and arm, including p-arm and q-arm panels
* Multiple panels include sample counts, e.g. `n=428`, `n=446`, etc.

**Description:**
This slide shows many small histograms faceted by chromosome arm. Blue histograms appear for p-arm panels and orange histograms for q-arm panels. Some panels are shaded pink, indicating missing or unavailable arms/regions. The x-axis is region length in kilobases, and the y-axis is count. The distributions show how long interchromosomal matching regions are across different chromosome arms, with some arms showing concentrated peaks and others broader distributions.

---

# Slide 7 — All-vs-all heatmap

**Title:** All-vs-all — Heatmap

**Visible text:**

* Missing introvert arms
* 2p, 3p, 5p, 8q, 11q, 14q
* Arm

  * p
  * q
* Color scale: 0.2, 0.6, 0.8, 1

**Description:**
The slide shows a clustered all-vs-all heatmap comparing chromosome arms. Dendrograms appear along the top and left, grouping arms by similarity. A red/blue side annotation marks p versus q arms. The strongest similarities appear as a red/orange diagonal, with additional off-diagonal blocks suggesting clusters of related subtelomeric regions. The upper-right note lists missing “introvert arms”: 2p, 3p, 5p, 8q, 11q, and 14q.

---

# Slide 8 — PCA by chromosome

**Title:** All-vs-all — PCA colored by chromosome

**Visible text:**

* Missing introvert arms
* 2p, 3p, 5p, 8q, 11q, 14q
* Chromosome

  * chr1 through chr22
  * chrX
  * chrY
* Arm

  * p
  * q
* Dimension 1 `(16.05%)`
* Dimension 2 `(11.2%)`

**Description:**
This slide shows a PCA scatterplot of chromosome-arm similarity patterns. Points are colored by chromosome, and point shapes distinguish p arms from q arms. Many points are labeled with chromosome arm and sample/count annotations. The scatter forms several visible clusters, suggesting groups of chromosome arms with similar subtelomeric similarity profiles. The same missing introvert arms are noted in the upper-right corner.

---

# Slide 9 — PCA by superpopulation

**Title:** All-vs-all — PCA colored by superpopulation

**Visible text:**

* Missing introvert arms
* 2p, 3p, 5p, 8q, 11q, 14q
* Superpopulation

  * AFR
  * AMR
  * EAS
  * EUR
  * SAS
* Arm

  * p
  * q
* Dimension 1 `(16.05%)`
* Dimension 2 `(11.2%)`

**Description:**
This slide uses the same PCA layout as the previous slide, but colors points by superpopulation instead of chromosome. The legend includes AFR, AMR, EAS, EUR, and SAS. The same p/q arm shapes are used. Points from different superpopulations appear mixed across clusters, suggesting that the major PCA structure may be driven more by chromosome-arm/community relationships than by superpopulation alone.

---

# Slide 10 — PCA by community

**Title:** All-vs-all — PCA colored by community

**Visible text:**

* Missing introvert arms
* 2p, 3p, 5p, 8q, 11q, 14q
* PARs-driven
* PHRs-driven
* DUX4-driven?

**Community list:**

```text
C1: 4q, 10q
C2: 10p, 18p
C3: 3q, 11p, 15q, 19p
C4: 7q, 12q
C5: 6p, 9p, 12p, 20q
C6: 1q, 13q, 17q, 19q, 21q, 22q
C7: 13p, 14p, 15p, 21p, 22p
C8: 16p
C9: 7p, 9q, 16q
C10: 17p
C11: 1p, 5q, 6q, 8p
C12: 2q, 20p
C13: 4p
C14: Xq, Yq
C15: 18q (only 1 case), Xp, Yp
```

**Description:**
The final slide shows the PCA scatterplot colored by inferred community. Labeled boxes identify communities C1 through C15 directly on the plot. Several broader regions are annotated with interpretive labels: “PARs-driven” near the lower-left/lower-center, “PHRs-driven” near the center, and “DUX4-driven?” near the lower-right/right-center. A legend maps each community to a color and number of arms. The right side lists each community and its associated chromosome arms, with C15 noting that 18q appears in only one case.
