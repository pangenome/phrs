# Tan et al. 2026 — verified extract from the full text

Tan K-T, Ong RJX, Kee BBR, Yap RKH, Ng AJT, Rajapaksha GK, Garrity-Janger M, Lin Q, Chong JR,
Kyi CT, HPRC, Meyerson M, Li H. "Widespread structural variations at human chromosome ends."
bioRxiv 2026. doi:10.64898/2026.09.03.748401. PMC13684503, PMID 42922562. CC BY-NC 4.0.

**Source of this extract:** the bioRxiv full-text PDF, retrieved through Exa's crawler (our own host is
IP-rate-limited to 429 by bioRxiv; see the retrieval record at the end). **Coverage: front matter,
Abstract, Introduction, all Results, Discussion — text available up to ~50,000 characters, which is
where the copy truncates (mid-Discussion). The METHODS SECTION WAS NOT RETRIEVED**, so aligner
choices, identity/length thresholds, divergence metric and any statistical parameters remain unknown.
Line numbers below are the preprint's own marginal numbering.

## 1. Dataset (verified)

- 1,670 assemblies from **five consortia** → "**860 unique haplotypes**" (L111-114), six ancestry
  groups (AFR, AMR, EAS, EUR, MENA, SAS) (L97-100).
- Subtelomere is operationalised as "the **500 kb** immediately proximal to the telomeric repeat
  arrays at each chromosome end" (L60-61) — identical operational definition to our windows.
- 860 assemblies × 41 arms = **35,150 subtelomeric sequences**, mean 40.9 per haploid assembly
  (L132-133).
- Arm assignment required "consistent assignment to the same arm across **at least four of five**
  reference sequences" (L115-116).
- Completeness: HPRC1 only 52% complete; HPRC2 >95% complete with >1,000 kb, "ensuring the 500 kb
  subtelomeric region ... was fully captured" (L123-129). **Only 12% of expected acrocentric contigs
  could be assigned** (L137-139) — acrocentric p arms are effectively unusable, as in our analysis.
  Sex chromosomes handled by accounting for PAR1/PAR2 (L140-143).

## 2. Structural-variation method (partially verified — Methods not retrieved)

From Results only: contigs were compared to CHM13 and "sequence divergence" computed (L149-151);
Stong et al. paralogy blocks were mapped to CHM13 subtelomeres and then to all pangenome contigs
(L159-163); contigs clustered "into major haplotypes differing by **at least two paralogy blocks of
≥3 kb**" (L170-172), with a second round of clustering for minor haplotypes (L186-188). Blocks were
present in 83% (40/48) of arms (L166-167); three or more major haplotypes in 29% (14/48) of arms,
11p most variable with seven (L172-174).

## 3. DATING — THE PRIORITY QUESTION IS NOW ANSWERED: they do NOT date by molecular clock

**There is no mutation rate, no molecular clock, no coalescent model and no calibrated age for any
block anywhere in Results or Discussion.** Their dating is inference from two lines of evidence:

**(a) Variance partitioning (L254-261):** "we used a variance-partitioning model with ancestry and
haplotype as explanatory variables ... Haplotype explained more of the divergence than ancestry ...
indicating that divergence is structured by arm and haplotype lineage rather than by ancestry,
consistent with ... most cross-arm block duplications are ancient, predating human population
divergence."

**(b) Great-ape presence/absence and block phylogenies (L264-302):** 45/50 blocks retained in Pan
troglodytes, only 6/50 in Gorilla gorilla (L270-271); human-only blocks 19a, 19b and 41 (L272-273);
"most subtelomeric paralogy blocks predate human speciation, being present across multiple great-ape
lineages" (L280-281).

The only "dated" numbers they quote are **literature species-divergence times used as context, not as
block ages**: Pongo species 0.96 vs Pan species 1.6 mya (L274-275); Gorilla 10.7 vs human 19.1 mya
(L279).

**Consequence for us:** their "ancient" is a relative statement (older than human population
divergence; shared with great apes). Our molecular-clock bound — d = 2·mu·T giving 5-12 Myr as an
upper bound, ~2-5 Myr at mu = 3e-8 — is therefore **quantitatively new relative to this paper**, not a
replication of it. We can say they establish antiquity qualitatively and we put a number on it.

## 4. Concerted evolution / homogenisation (verified, and very close to our framing)

They classify blocks into three evolutionary categories (L284-302), and the middle one is exactly our
concerted-evolution argument:

- "**Block 12**'s topology matches the species tree ... suggesting the block had a post-speciation
  expansion with **limited interchromosomal exchange**" (L286-290).
- "**Block 8** ... without arm clustering ... indicating post-speciation expansion followed by
  **frequent interchromosomal exchange that homogenized its sequence across arms**" (L290-293).
- "**Block 36** ... bipartite topology ... a lineage-specific expansion after the human Pan/Pongo
  split" (L293-298).
- Conclusion: "these patterns show markedly distinct evolutionary trajectories across blocks, with
  **interchromosomal exchange frequency varying drastically**" (L300-302).
- Block 20, present across many arms, "showed poor stratification at both levels, suggesting greater
  sequence homogeneity and **more frequent inter-arm exchange**" (L250-252).

They also quantify block copy number and arm specificity: **Block 40 median CN 94, restricted to 4q
and 10q** — the same pair as our community C1; **Block 2 median CN 17 across 23 arms**, length
18-26 kb (L231-236).

## 5. Telomere length (verified)

Simulated reads from terminal 50 kb mapped to CHM13 69%, pangenome 94%, self-assembly 99% (L430-433).
Real HPRC2 telomeric reads: pangenome concordance median 100% (mean 87.2%) vs CHM13 median 91.7%
(mean 66.6%) (L456-459). Relative telomere length spans -875 bp (20q) to +534 bp (4p) (L461-462).
Haplotype-telomere-length association: only 1p (q=0.078), 5q (q=0.090), 6q (q=0.087) at relaxed
FDR<0.1 — "suggestive rather than establish" (L468-472).

## 6. Citation of Guarracino et al. 2026 — STILL UNVERIFIED (references not retrieved)

The truncated copy ends in the Discussion, so the reference list is unavailable. Nothing in Results or
Discussion discusses our preprint. Unresolved, not disproven.

## 7. Their two headline exchange results — direct comparanda for our work

**(a) 4q/10q D4Z4 exchange (L411-423).** "Among 4q-S1 contigs, **7/275 (2.5%)** showed evidence of
10q-S1-derived sequence displacing the functional PAS at their terminal ends ... Conversely, **3/514
(0.58%)** 10q-S1 contigs carried 4q-S1-derived sequence introducing a functional PAS." Orthogonal
check: "we independently confirmed a detected translocation using **HiFi long reads and the
corresponding Hifiasm assembly of the same sample**" (L418-420).

This is the same locus as our C1 (4q/10q) candidate, the same read-level validation strategy we
applied at 9q/3q, and a comparable frequency to our 2.3% of arms carrying a strong candidate.

**(b) Arm-misassignment map (L437-448) — IMPORTANT FOR OUR 9q/3q CLAIM.** "Reads from 11q and 13q
mapped to CHM13 with 95% median accuracy, while **8p and 9q reads mapped with 0% median accuracy**.
These reads were misassigned to other arms sharing the same subtelomeric blocks (**8p reads
misassigned to 6q; 9q reads misassigned to 11p and 3q**)."

So the population-wide data say 9q terminal sequence is **routinely misassigned to 3q and 11p** —
independent confirmation that 9q and 3q carry strongly homologous terminal sequence, but also the
reason to be careful: the 9q/3q homology is a **general population feature**, not something unique to
the pedigree. Our pedigree result is narrower and still valid — it concerns *which paternal haplotype*
carries the module (DAD chr9:h2 yes, chr9:h1 no) and its transmission through two meioses at
99.96-99.97% — but it must be described as one sampled carrier of a common architecture, not as a
private rearrangement.

## 8. Their positional claim, which bears on our telomere-distance profile (L178-182)

"All-to-all comparisons between representative contigs of each major haplotype ... showed that
**differences between haplotypes of the same arm lay largely at the telomeric end, while centromeric
ends were more likely to be homologous**."

That is a telomere-proximity effect at the *haplotype* level, obtained by contig comparison rather than
by our distance binning. It is consistent in direction with what we see structurally, and it is not
subject to the allelic-denominator failure that invalidated our first 10 kb of bins — so their
framing is the one to cite for the positional gradient.

## Retrieval record (complete; do not repeat from this host)

- bioRxiv PDF via our own tools: **HTTP 429 at the IP level**, with and without a browser user agent,
  on both `v1.full.pdf` and `biorxiv/early/...` routes.
- Europe PMC `fullTextXML` for PMC13684503: front matter + *"The license terms selected by the
  author(s) for this preprint version do not permit archiving in PMC."*
- OpenAlex W7211892909: `oa_status: green`, single OA PDF = the blocked bioRxiv URL.
- **Exa web search with `includeContent` WORKED** and returned the PDF text (truncated at 50,000
  chars, ending mid-Discussion). This is the route to use again if the Methods are needed.
- Delegated Codex run: could not write files (`bwrap: loopback: Failed RTM_NEWADDR`); use the
  `sxgc.*` user lanes on this host.

**Remaining gap:** the Methods section (aligner, thresholds, divergence metric) and the reference list
(citation check). Re-running Exa against the same URL may only return the same truncated slice; the
most reliable fix is still a browser download of the PDF.
