# Slide Text

## Slide 1: Copy Number Changes The Enrichment Question

Subtitle: Standard ORA sees symbols; subtelomeric biology repeats copies.

- Standard ORA counts each gene symbol once, so a gene family copied across many PHR arms can still enter as one hit.
- PHR/subtelomeric duplicons are copy-rich by design: repeated genes, pseudogenes, lncRNAs, retrogenes, and nearby repeats are part of the biology.
- Copy-number-aware interpretation asks which families accumulate across gene instances, arms, and communities, not only which unique names appear.
- Keep the scopes separate: community copy/support patterns are slide-worthy, while formal copy-weighted ORA needs validated statistics.

## Slide 2: OR4F Is The Clean Copy-Rich Signal

Subtitle: Olfactory receptor copies expose subtelomeric gene turnover.

- OR4F is a subtelomeric olfactory receptor family; coding and pseudogene copies mark ongoing birth, decay, and duplication at chromosome ends.
- HPRCv2 source data show 5,023 OR4F annotations across 16 arms: 3,117 pseudogene and 1,906 coding annotations.
- The arm gradient is wide: 62.1 percent pseudogene overall, ranging from 11.1 percent at chr7p to 99.8 percent at chr15q.
- Community support is broad but qualitative: 10 OR4F genes span 7 communities, and OR4F5/OR4F8P each appear on 14 arms; the C3 OR Fisher row has q ~= 0.071.

## Slide 3: Duplicon Backbone And GTP Anchors

Subtitle: The repeated architecture extends beyond olfactory receptors.

- RPL23A pseudogenes, DDX11L, WASH, and FAM138 mark a recurring subtelomeric duplicon backbone across many communities and arms.
- Example hub signals are highly copied across the arm-level table: RPL23AP45 spans 10 communities/21 arms, DDX11L16 spans 9/20, FAM138D spans 9/17, and WASH6P spans 8/21.
- C5 carries a DDX11L/WASH/FAM138 module with IQSEC3 as a gene-level anchor; C15/PAR1 includes SHOX and the GTP-binding gene GTPBP6.
- Older copy-weighted ORA also flagged GTP binding/GTPase through IQSEC3/GTPBP6, but those fold enrichments remain exploratory until recalibrated.

## Slide 4: Boundaries For The Claim

Subtitle: Community support is not the same as a called interval or validated GO result.

- Some community-assigned arms lack called CHM13 PHR intervals in `chm13.phrs.bed`: C5 `chr6_p`, C7 `chr13_p`, C14 `chrY_q`, and C15 `chrY_p`.
- The HPRCv2 family Fisher screen tested 116 community-family rows; none survives BH correction, so present family signals as qualitative architecture.
- The parked weighted-hypergeometric ORA is exploratory because p-value calibration and FDR control failed under gene-level sampling.
- TAR1 is useful repeat context, especially near telomeres and in C2, but it is not a gene-enrichment result.
