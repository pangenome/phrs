# A2 Bibliography Duplication and Citation Audit

Date: 2026-06-17
Task: `manuscript-revision-a2`

## Scope

Inspected:

- `submission/paper.tex`
- `submission/bibliography.bib`
- build state in `submission/`
- revision prompt note in `wg_manuscript_revision_prompt.md` lines 121-122

The active manuscript uses `multibib` with a main bibliography and a Methods
bibliography:

- `submission/paper.tex` lines 16-20 define `\newcites{Meth}{Methods References}`.
- `submission/paper.tex` line 293 prints the main bibliography.
- `submission/paper.tex` lines 665-666 print the Methods bibliography.
- No `\citeSupp{...}` commands are present, so the Supplementary bibliography is
  currently empty.

The `.bib` database itself does not contain duplicate BibTeX records by DOI or
normalized title in the inspected file. The duplication is caused by citing the
same BibTeX key in both `\cite{...}` and `\citeMeth{...}`, which prints the same
paper twice under different reference numbers.

## Build Artifact Status

No pre-existing `.aux`, `.bbl`, or `.blg` artifacts were present in
`submission/` at audit start. I attempted `make` in `submission/` to recover exact
current numbering, but the local TeX environment failed on the first pass because
`geometry.sty` is missing:

```text
! LaTeX Error: File `geometry.sty' not found.
```

The current numbers below are therefore reconstructed mechanically from citation
order in `submission/paper.tex`. This matches the prompt's numbering pattern:
39 main references followed by 28 Methods references, with Methods numbers
continuing at 40 because `resetlabels` is not used.

## Duplicate Main/Methods Reference Groups

These entries are cited once in the main bibliography and again in the Methods
bibliography. The `Methods printed no.` column is the apparent global reference
number if the manuscript compiles with the current `multibib` setup.

| BibTeX key | Main no. | Methods local no. | Methods printed no. | Main citation lines | Methods citation lines | Paper |
|---|---:|---:|---:|---|---|---|
| `Nurk2022` | 16 | 2 | 41 | 92, 100 | 439 | The complete sequence of a human genome |
| `hprc_hprcv2_2025` | 18 | 1 | 40 | 102 | 439 | HPRC Data Release 2 |
| `Guarracino2023` | 19 | 4 | 43 | 116 | 464 | Recombination between heterologous human acrocentric chromosomes |
| `Garrison2024pggb` | 20 | 8 | 47 | 119 | 490 | Building pangenome graphs |
| `Francis2025` | 28 | 23 | 62 | 210 | 603 | Complete genome assemblies of two mouse subspecies |
| `Zuo2021` | 29 | 24 | 63 | 210 | 616 | Stage-resolved Hi-C analyses of mouse meiosis |
| `Ulahannan2019` | 30 | 13 | 52 | 230 | 556 | Pore-C concatemer sequencing |
| `hic3d_cifi2025` | 31 | 14 | 53 | 232 | 556 | CiFi long-read chromosome conformation capture |
| `Tan2018` | 32 | 18 | 57 | 232 | 591 | Dip-C single-cell genome structures |
| `Xu2025` | 33 | 19 | 58 | 233 | 594 | Human sperm single-cell 3D genome organization |
| `Cechova2025` | 34 | 25 | 64 | 242 | 624 | Telomere-to-telomere assembly of a multi-generation human pedigree |

Prompt A-2 listed the same core duplicates as examples: HPRC release 18/40,
T2T-CHM13 16/41, acrocentric recombination 19/43, mouse assemblies 28/62,
Zuo Hi-C 29/63, pedigree 34/64, Pore-C 30/52, CiFi 31/53, Dip-C 32/57, and sperm
scHi-C 33/58. The audit additionally flags `Garrison2024pggb` as a duplicated
main/Methods entry at 20/47.

## Wrong or Mismatched Citation Targets

### Definite mismatch: IMPG transitive closure cites ODGI and Minigraph-Cactus

Location:

- `submission/paper.tex` lines 469-473, subsection
  `Implicit pangenome graph and IMPG transitive closure`

Current command:

```tex
\citeMeth{pangenome_graphs_impg_GuarracinoHeumos2022,pangenome_graphs_impg_Hickey2024}
```

Current target papers:

- `pangenome_graphs_impg_GuarracinoHeumos2022` is ODGI, not IMPG
  (`submission/bibliography.bib` lines 662-672).
- `pangenome_graphs_impg_Hickey2024` is Minigraph-Cactus, not IMPG
  (`submission/bibliography.bib` lines 675-684).

Correct target for IMPG:

- `pangenome_graphs_impg_IMPG2023`, currently main no. 22, defined at
  `submission/bibliography.bib` lines 687-692.

Safe patch:

- Replace the Methods citation at `submission/paper.tex` line 473 with
  `\cite{pangenome_graphs_impg_IMPG2023}` after the bibliography merge, or with
  `\citeMeth{pangenome_graphs_impg_IMPG2023}` only if the two-list structure is
  temporarily retained.
- Do not cite `pangenome_graphs_impg_Hickey2024` for this sentence unless a new
  sentence is added that explicitly discusses Minigraph-Cactus.

### Related odgi placement issue

Location:

- `submission/paper.tex` lines 486-492, subsection
  `Pangenome graph and Jaccard similarity`

Current command:

```tex
\citeMeth{Garrison2024pggb,andreace2023pangenome,heumos2024nfcore}
```

The sentence explicitly mentions `odgi similarity`, `odgi layout`, and
`odgi draw`, but it does not cite the ODGI paper. The ODGI key is currently
misused in the IMPG subsection above.

Safe patch:

- Move `pangenome_graphs_impg_GuarracinoHeumos2022` from the IMPG subsection to
  this pangenome-graph/Jaccard subsection, e.g. cite
  `\cite{Garrison2024pggb,pangenome_graphs_impg_GuarracinoHeumos2022,andreace2023pangenome,heumos2024nfcore}`.
- Keep `pangenome_graphs_impg_Hickey2024` out of this command unless the text
  explicitly needs a Minigraph-Cactus contrast.

### Likely stray citation: PHR detection cites VG toolkit

Location:

- `submission/paper.tex` lines 479-484, subsection `PHR detection`

Current command:

```tex
\citeMeth{Garrison2018}
```

Current target paper:

- `Garrison2018`: "Variation graph toolkit improves read mapping by representing
  genetic variation in the reference" (VG toolkit).

The cited sentence reports IMPG sliding-window detection thresholds and PHR
counts, not VG read mapping. This appears to resolve to the wrong paper or to be
a leftover citation.

Safe patch options:

- If a citation is desired for the sliding-window query method, cite
  `pangenome_graphs_impg_IMPG2023`.
- If the sentence is reporting only this manuscript's output counts, remove the
  citation rather than replacing it.

## Affected Citation Commands

All Methods citations that would need conversion or review in a single-list
patch:

- `submission/paper.tex` line 439:
  `\citeMeth{hprc_hprcv2_2025,Nurk2022,logsdon2025hgsvc}`
- `submission/paper.tex` line 464:
  `\citeMeth{Guarracino2023}`
- `submission/paper.tex` line 473:
  `\citeMeth{pangenome_graphs_impg_GuarracinoHeumos2022,pangenome_graphs_impg_Hickey2024}`
- `submission/paper.tex` line 484:
  `\citeMeth{Garrison2018}`
- `submission/paper.tex` line 490:
  `\citeMeth{Garrison2024pggb,andreace2023pangenome,heumos2024nfcore}`
- `submission/paper.tex` line 556:
  `\citeMeth{hic3d_dixon2012,hic3d_imakaev2012,Ulahannan2019,hic3d_cifi2025}`
- `submission/paper.tex` line 587:
  `\citeMeth{hic3d_alavattam2019,hic3d_wolff2018,hic3d_deshpande2022}`
- `submission/paper.tex` line 591:
  `\citeMeth{Tan2018}`
- `submission/paper.tex` line 594:
  `\citeMeth{Xu2025,hic3d_scnanoHiC2023,hic3d_scnanoHiC2_2025,hic3d_kitamura2025}`
- `submission/paper.tex` line 603:
  `\citeMeth{Francis2025}`
- `submission/paper.tex` line 616:
  `\citeMeth{Zuo2021}`
- `submission/paper.tex` line 624:
  `\citeMeth{Cechova2025,Porubsky2025}`
- `submission/paper.tex` line 627:
  `\citeMeth{pedigree_Porsborg2025primaterecom,pedigree_Joseph2024PRDM9indep}`

No `\citeSupp{...}` commands were found.

## Safe Merge and Renumber Plan

Recommended integration patch for Cluster A fan-in:

1. Convert the manuscript to one reference list.
   - In `submission/paper.tex` lines 16-20, remove the `multibib` setup or stop
     defining/using the `Meth` bibliography.
   - Replace every `\citeMeth{...}` command listed above with `\cite{...}`.
   - Keep `submission/bibliography.bib` intact during this mechanical pass unless
     unused entries are deliberately pruned later.

2. Move the single bibliography to the final reference-list location.
   - Remove the current main `\bibliography{bibliography}` at
     `submission/paper.tex` line 293, where it appears before tables, figures,
     Methods, and Extended Data legends.
   - Replace the Methods bibliography block at `submission/paper.tex` lines
     665-666 with the single bibliography:

     ```tex
     \bibliographystyle{unsrt}
     \bibliography{bibliography}
     ```

   - This keeps all references in one numbered sequence and avoids printing a
     reference list before later Methods citations.

3. Repair the IMPG/odgi/Minigraph-Cactus mismatch while converting citations.
   - At `submission/paper.tex` line 473, cite
     `pangenome_graphs_impg_IMPG2023`.
   - At `submission/paper.tex` line 490, add
     `pangenome_graphs_impg_GuarracinoHeumos2022` to the odgi/Jaccard citation
     group.
   - Remove `pangenome_graphs_impg_Hickey2024` from the IMPG sentence. Retain it
     only if another sentence explicitly discusses Minigraph-Cactus.

4. Resolve the likely stray VG citation.
   - At `submission/paper.tex` line 484, either remove `Garrison2018` or replace
     it with `pangenome_graphs_impg_IMPG2023`, depending on whether the final
     sentence is meant to cite the method or merely report the present study's
     output counts.

5. Rebuild and check numbering.
   - Run `cd submission && make clean && make` in an environment with the required
     LaTeX packages.
   - Confirm only one `.bbl` is generated for the active reference list.
   - Confirm the duplicate pairs above no longer print twice.
   - Confirm the IMPG transitive-closure sentence resolves to the IMPG entry, not
     to ODGI or Minigraph-Cactus.

## Files and Sections for Integration

Primary integration file:

- `submission/paper.tex`

Specific sections:

- Preamble bibliography setup: lines 16-20.
- Main bibliography placement: line 293.
- Methods `Sample selection and reference frame`: line 439.
- Methods `wfmash all-vs-all alignment`: line 464.
- Methods `Implicit pangenome graph and IMPG transitive closure`: lines 469-473.
- Methods `PHR detection`: lines 479-484.
- Methods `Pangenome graph and Jaccard similarity`: lines 486-492.
- Methods 3D/pedigree citation commands: lines 556, 587, 591, 594, 603, 616,
  624, and 627.
- Methods bibliography block: lines 665-666.

Reference database for target-key verification:

- `submission/bibliography.bib` lines 662-672: ODGI key
  `pangenome_graphs_impg_GuarracinoHeumos2022`.
- `submission/bibliography.bib` lines 675-684: Minigraph-Cactus key
  `pangenome_graphs_impg_Hickey2024`.
- `submission/bibliography.bib` lines 687-692: IMPG key
  `pangenome_graphs_impg_IMPG2023`.

