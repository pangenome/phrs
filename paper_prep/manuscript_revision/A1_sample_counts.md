# A1 Sample Counts Audit

Date: 2026-06-17

## Decision

Use this count triple for the active manuscript:

| Count type | Confirmed count | Scope |
| --- | ---: | --- |
| HPRC Release 2 individuals used by the PHR analysis | 232 | Distinct non-reference sample IDs in the local PanSN-to-FASTA index, including HG002 |
| HPRC haplotype assemblies used by the PHR analysis | 464 | Two haplotype assemblies for each of the 232 HPRC individuals |
| Total assemblies in the PHR analysis | 465 | 464 HPRC haplotypes plus CHM13v2.0; GRCh38 is present in the wider assembly index/mirror but is not present in the PHR flank/region tables |

The manuscript should not describe the PHR analysis as 233 individuals, 465 HPRC haplotypes, or 466 analyzed assemblies. Those numbers conflate CHM13 with HPRC haplotypes, or conflate the wider Release 2 assembly/index mirror with the PHR-specific input set.

## Source Paths Inspected

Local sources:

| Source | Evidence |
| --- | --- |
| `submission/paper.tex` | Active manuscript text and figure legends. Count-sensitive occurrences are listed below with line numbers. |
| `AGENTS.md` | Project guide currently says "233 HPRCv2 samples (465 haplotypes...)"; this is inconsistent with the PHR inputs and should be updated with `CLAUDE.md` in lock-step by a later doc-maintenance task. |
| `CLAUDE.md` | Same count text as `AGENTS.md`; included because the files are required to stay in lock-step. |
| `subtelomeric_analysis_report.md` | Contains multiple 233/465 prose occurrences. Some source inventory text also labels 465 FASTA files as HPRCv2 assemblies while including CHM13 in that count. |
| `/moosefs/guarracino/HPRCv2/PHR_III/pansn_to_fasta.tsv` | 466 PanSN entries: 464 HPRC haplotypes, CHM13, and GRCh38. Every non-reference HPRC sample ID has exactly two entries. |
| `/moosefs/pangenomes/HPRCv2/*.fa.gz` | 466 FASTA files in the local mirror: 464 HPRC haplotype assemblies plus CHM13 plus GRCh38. |
| `/moosefs/guarracino/HPRCv2/PHR_III/hprcv2.1Mb.telo_500kb_trimmed.fa.gz.fai` | 18,827 flank records from 465 haplotype prefixes: 464 HPRC haplotypes plus CHM13; no GRCh38 prefix. |
| `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv` | PHR region table contains 465 haplotype prefixes: 464 HPRC haplotypes plus CHM13; no GRCh38 prefix. |
| `/moosefs/guarracino/HPRCv2/PHR_III/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.fai` | 15,668 signal-bearing PHR sequence records downstream of the same analyzed set. |
| `/moosefs/guarracino/HPRCv2/PHR_III/liftoff_genes_hprc_r2_v1.0.index.csv` | 462 indexed annotation rows plus header; annotation coverage evidence only, not the PHR sequence-analysis denominator. |
| `/moosefs/guarracino/HPRCv2/PHR_III/repeat_masker_bed_hprc_r2_v1.0.index.csv` | 462 indexed annotation rows plus header; annotation coverage evidence only, not the PHR sequence-analysis denominator. |
| `/moosefs/guarracino/HPRCv2/data/hprc-sequence-production.tsv` | Broader production metadata table. It has many samples outside the PHR assembly set and lacks several local PanSN IDs; do not use it as the assembly count source. |

Discoverable release metadata:

| Source | Evidence |
| --- | --- |
| HPRC Data Release 2 page, `https://humanpangenome.org/hprc-data-release-2` | States that Release 2 incorporates genomes from 232 individuals. |
| HPRC Release 2 assembly repository, `https://github.com/human-pangenomics/hprc_intermediate_assembly` | Documents that Release 2 assembly metadata are in `data_tables/` and that the assembly index includes reference-level assemblies CHM13, GRCh38, and HG002 to standardize analysis versions. |
| HPRC Release 2 assembly index, `https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/main/data_tables/assemblies_release2_v1.0.index.csv` | Contains 466 assembly rows: 462 haplotypes for 231 non-reference HPRC/HPP sample IDs, 2 HG002 reference haplotypes, CHM13, and GRCh38. This is consistent with 232 HPRC individuals and explains why the wider release/index total is not the PHR analyzed total. |
| NIH HPRC Assemblies Release 2 mirror, `https://hpc.nih.gov/refdb/dbview.php?id=1176` | Lists 466 relevant files for Release 2, including CHM13 and GRCh38 alongside HPRC assemblies; consistent with the local mirror but not the PHR analysis denominator. |

## Count Checks

Commands run from `/moosefs/erikg/phrs`:

```sh
awk -F'\t' '{print $1}' /moosefs/guarracino/HPRCv2/PHR_III/pansn_to_fasta.tsv |
  awk -F'#' '{print $1}' | sort | uniq -c |
  awk 'BEGIN{hprc_ind=0; hprc_hap=0; total=0; chm13=0; grch38=0}
       {total+=$1;
        if ($2=="CHM13") chm13+=$1;
        else if ($2=="GRCh38") grch38+=$1;
        else {hprc_ind++; hprc_hap+=$1; if ($1!=2) print "NON_DIPLOID",$2,$1}}
       END{print "hprc_individuals",hprc_ind;
           print "hprc_haplotypes",hprc_hap;
           print "chm13",chm13;
           print "grch38",grch38;
           print "total_pansn_entries",total}'
```

Observed:

```text
hprc_individuals 232
hprc_haplotypes 464
chm13 1
grch38 1
total_pansn_entries 466
```

PHR flank index check:

```sh
awk -F'\t' '{split($1,a,"#"); print a[1]"#"a[2]}' \
  /moosefs/guarracino/HPRCv2/PHR_III/hprcv2.1Mb.telo_500kb_trimmed.fa.gz.fai |
  sort -u |
  awk -F'#' 'BEGIN{h=0;c=0;g=0}
              {if($1=="CHM13") c++; else if($1=="GRCh38") g++; else h++}
              END{print "hprc_hap_prefixes_with_flanks",h;
                  print "chm13_hap_prefixes_with_flanks",c;
                  print "grch38_hap_prefixes_with_flanks",g;
                  print "total_hap_prefixes_with_flanks",h+c+g}'
```

Observed:

```text
hprc_hap_prefixes_with_flanks 464
chm13_hap_prefixes_with_flanks 1
grch38_hap_prefixes_with_flanks 0
total_hap_prefixes_with_flanks 465
```

PHR region table check:

```sh
awk -F'\t' 'NR>1{split($1,a,"#"); print a[1]"#"a[2]}' \
  /moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv |
  sort -u |
  awk -F'#' 'BEGIN{h=0;c=0;g=0}
              {if($1=="CHM13") c++; else if($1=="GRCh38") g++; else h++}
              END{print "hprc_hap_prefixes_region_table",h;
                  print "chm13_hap_prefixes_region_table",c;
                  print "grch38_hap_prefixes_region_table",g;
                  print "total_hap_prefixes_region_table",h+c+g}'
```

Observed:

```text
hprc_hap_prefixes_region_table 464
chm13_hap_prefixes_region_table 1
grch38_hap_prefixes_region_table 0
total_hap_prefixes_region_table 465
```

## Affected `submission/paper.tex` Locations

These are the active manuscript occurrences that need mechanical replacement. Line numbers refer to the file state audited on 2026-06-17.

| Line(s) | Current text | Issue | Proposed replacement |
| ---: | --- | --- | --- |
| 55 | `466 near-complete assemblies (465 HPRC v2 haplotypes together with CHM13)` | Off by one in both total assemblies and HPRC haplotypes; 465 is total analyzed assemblies after adding CHM13. | `465 near-complete assemblies (464 HPRC v2 haplotypes from 232 individuals, together with CHM13)` |
| 101-102 | `HPRC v2) provides near-complete, haplotype-resolved assemblies for 233 individuals from five superpopulations` | Release 2 and local PHR inputs support 232 individuals. | `HPRC v2) provides near-complete, haplotype-resolved assemblies for 232 individuals from five superpopulations` |
| 113-114 | `From the 233 individuals we used 465 HPRC v2 haplotype assemblies together with CHM13v2.0` | Should be 232 HPRC individuals and 464 HPRC haplotype assemblies; CHM13 makes total 465. | `From 232 HPRC v2 individuals we used 464 haplotype assemblies together with CHM13v2.0` |
| 308-309 | `Interchromosomal homology across 465 HPRC v2 haplotypes plus CHM13.` | Figure legend counts CHM13 in addition to 465 HPRC haplotypes, implying 466 analyzed assemblies. | `Interchromosomal homology across 464 HPRC v2 haplotypes plus CHM13.` |
| 434-435 | `233 HPRC v2 v1.1 individuals contributed haplotype-resolved assemblies, and CHM13v2.0 was added as the closed-reference anchor.` | Should be 232 HPRC individuals. Also "v1.1" is potentially confusing because the release/index uses mixed assembly file versions, including v1.0.1 and v1.1.0 for corrected haplotype swaps; the manuscript can simply say HPRC v2. | `232 HPRC v2 individuals contributed haplotype-resolved assemblies, and CHM13v2.0 was added as the closed-reference anchor.` |
| 436-437 | `The arm-flank census uses 465 haplotypes (233 samples \times 2, minus 1 because CHM13 contributes a single haploid sequence).` | Formula is wrong: CHM13 is not one of the HPRC samples. The census uses 464 HPRC haplotypes plus one haploid CHM13 assembly. | `The arm-flank census uses 465 assemblies: 464 HPRC haplotype assemblies (232 individuals \times 2) plus the haploid CHM13 assembly.` |
| 438-439 | `Superpopulation breakdown: AFR 67, EAS 52, AMR 44, SAS 37, EUR 33` | These values sum to 233 and therefore cannot describe the confirmed PHR assembly set. Needs recomputation from an exact sample-to-superpopulation table before keeping any breakdown. | Replace with either no breakdown sentence, or: `Superpopulation labels were used only for downstream population summaries and were assigned at the HPRC individual level.` |

No other `submission/paper.tex` lines with `HPRC`, `haplotypes`, `assemblies`, `samples`, `individuals`, or the competing counts `232/233/464/465/466` require replacement for this count reconciliation.

## Exact Proposed Wording

Abstract:

```tex
Here we measure this sharing across 465 near-complete assemblies (464 HPRC v2 haplotypes from 232 individuals, together with CHM13) using an implicit pangenome graph that samples approximately 12\% of haplotype-pair comparisons and queries transitive relationships without chromosomal or positional priors.
```

Introduction/background:

```tex
Consortium (HPRC v2) provides near-complete, haplotype-resolved assemblies for
232 individuals from five superpopulations \cite{hprc_hprcv2_2025}.
```

Results, implicit pangenome graph section:

```tex
From 232 HPRC v2 individuals we used 464 haplotype assemblies together with
CHM13v2.0, extracted 18,827 telomere-anchored 500 kb flanks across the 48
chromosome arms, and aligned them all-against-all with wfmash at 95\% minimum
identity \cite{Guarracino2023}, without restricting alignment to homologous
chromosomes.
```

Figure 1 legend:

```tex
\textbf{(A)} Interchromosomal homology across 464 HPRC v2 haplotypes plus
CHM13.
```

Methods, sample selection and reference frame:

```tex
232 HPRC v2 individuals contributed haplotype-resolved assemblies, and
CHM13v2.0 was added as the closed-reference anchor.
The arm-flank census uses 465 assemblies: 464 HPRC haplotype assemblies
(232 individuals $\times$ 2) plus the haploid CHM13 assembly.
```

If the superpopulation breakdown must be retained, recompute it from the exact 232-sample assembly set before editing. Do not carry forward the current `AFR 67, EAS 52, AMR 44, SAS 37, EUR 33` sentence because it sums to 233.

## Notes for Downstream Fan-in

- Do not edit `submission/paper.tex` from this A1 task; the downstream A fan-in owns mechanical manuscript edits.
- Treat "HPRC haplotypes" and "total analyzed assemblies" as separate denominators:
  - HPRC haplotypes: 464.
  - Analyzed assemblies after adding CHM13: 465.
  - Wider local/public assembly inventory including GRCh38: 466, but GRCh38 is not in the PHR flank/region tables.
- The public Release 2 assembly index includes CHM13, GRCh38, and HG002 reference-level assemblies to standardize versions. This is why a raw "assembly index rows" count is not the same thing as "HPRC haplotypes analyzed for PHRs."
- `AGENTS.md`, `CLAUDE.md`, and `subtelomeric_analysis_report.md` also contain 233/465 language, but this task requested a decision artifact and no direct edits to those files.
