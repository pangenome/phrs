# A4 Guarracino Acrocentric Great-Apes DOI Check

Date: 2026-06-17

## Question

Verify whether the local `acrocentric_Guarracino2025ape` bibliography entry is real and correctly formatted:

- title: `Origin and evolution of acrocentric chromosomes in human and great apes`
- DOI prefix/stem: `10.64898/2025.12.22.696095`
- expected date: December 2025

Do not edit the bibliography directly.

## Current Local Entry

Source: `submission/bibliography.bib`

```bibtex
@misc{acrocentric_Guarracino2025ape,
  author  = {Guarracino, Andrea and others},
  title   = {Origin and evolution of acrocentric chromosomes in human and great apes},
  howpublished = {bioRxiv},
  year    = {2025},
  doi     = {10.64898/2025.12.22.696095},
  url     = {https://www.biorxiv.org/content/10.64898/2025.12.22.696095v1},
  note    = {Cross-species analysis of acrocentric chromosome architecture in
             chimpanzee and bonobo; confirmed chr14 inversion is human-specific,
             making SST1-mediated 13;14 and 14;21 Robertsonian fusions human-unique.}
}
```

## Verification Sources

Network lookup was available and used after checking the local bibliography.

- bioRxiv record: <https://www.biorxiv.org/content/10.64898/2025.12.22.696095v1>
- DOI resolver: <https://doi.org/10.64898/2025.12.22.696095>
- Crossref metadata API: <https://api.crossref.org/works/10.64898/2025.12.22.696095>
- PubMed record: <https://pubmed.ncbi.nlm.nih.gov/41497641/>
- NCBI citation API used for structured citation text: <https://api.ncbi.nlm.nih.gov/lit/ctxp/v1/pubmed/?format=citation&id=41497641>
- bioRxiv December 2025 archive entry: <https://connect.biorxiv.org/archive/2025-12>

## Findings

The entry is real. DOI `10.64898/2025.12.22.696095` resolves and Crossref lists it as a bioRxiv/openRxiv posted-content preprint with title `Origin and evolution of acrocentric chromosomes in human and great apes`.

The date should be treated as posted/published `2025-12-23`, not simply "December 2025" and not December 22. The `2025.12.22` component is part of the bioRxiv DOI/article identifier, but Crossref, PubMed, NCBI citation output, and the bioRxiv archive all report posting/publication on December 23, 2025.

The current local author field is not correctly formatted for this paper. Andrea Guarracino is an author, but he is not the first author. The first author is Steven J. Solar. PubMed/NCBI citation output gives:

`Solar SJ, Hebbar P, de Lima LG, Sweeten A, Rhie A, Potapova T, de Gennaro L, Guarracino A, Kim J, Pickett BD, Paten B, Wilson MA, Koren S, Garrison E, Eichler EE, Ventura M, Gerton JL, Phillippy AM. Origin and evolution of acrocentric chromosomes in human and great apes. bioRxiv [Preprint]. 2025 Dec 23:2025.12.22.696095. doi: 10.64898/2025.12.22.696095. PMID: 41497641; PMCID: PMC12767296.`

The local key `acrocentric_Guarracino2025ape` is therefore misleading if keys are expected to track first author. Keeping the key is acceptable only as a mechanical no-citation-edit choice; a future bibliography cleanup should rename it to something like `acrocentric_Solar2025ape` and update cite commands.

## Recommended Correction

Use the existing key if the A fan-in wants a low-risk bibliography-only replacement without editing manuscript cite commands:

```bibtex
@misc{acrocentric_Guarracino2025ape,
  author       = {Solar, Steven J. and Hebbar, Prajna and de Lima, Leonardo Gomes and Sweeten, Alex and Rhie, Arang and Potapova, Tamara and de Gennaro, Luciana and Guarracino, Andrea and Kim, Juhyun and Pickett, Brandon D. and Paten, Benedict and Wilson, Melissa A. and Koren, Sergey and Garrison, Erik and Eichler, Evan E. and Ventura, Mario and Gerton, Jennifer L. and Phillippy, Adam M.},
  title        = {Origin and evolution of acrocentric chromosomes in human and great apes},
  howpublished = {bioRxiv preprint},
  year         = {2025},
  month        = dec,
  doi          = {10.64898/2025.12.22.696095},
  url          = {https://www.biorxiv.org/content/10.64898/2025.12.22.696095v1},
  note         = {Posted December 23, 2025. PMID: 41497641; PMCID: PMC12767296.}
}
```

If the fan-in is allowed to update cite keys, prefer this semantically corrected key instead:

```bibtex
@misc{acrocentric_Solar2025ape,
  author       = {Solar, Steven J. and Hebbar, Prajna and de Lima, Leonardo Gomes and Sweeten, Alex and Rhie, Arang and Potapova, Tamara and de Gennaro, Luciana and Guarracino, Andrea and Kim, Juhyun and Pickett, Brandon D. and Paten, Benedict and Wilson, Melissa A. and Koren, Sergey and Garrison, Erik and Eichler, Evan E. and Ventura, Mario and Gerton, Jennifer L. and Phillippy, Adam M.},
  title        = {Origin and evolution of acrocentric chromosomes in human and great apes},
  howpublished = {bioRxiv preprint},
  year         = {2025},
  month        = dec,
  doi          = {10.64898/2025.12.22.696095},
  url          = {https://www.biorxiv.org/content/10.64898/2025.12.22.696095v1},
  note         = {Posted December 23, 2025. PMID: 41497641; PMCID: PMC12767296.}
}
```

## Recommendation

Do not remove this citation. It is a real December 2025 bioRxiv preprint with the DOI currently used locally.

For manuscript revision A fan-in, replace the local entry's author/date/citation metadata with the corrected BibTeX above. Use the existing key for the smallest mechanical change, or rename to `acrocentric_Solar2025ape` only if the fan-in also updates all `\cite{acrocentric_Guarracino2025ape}` occurrences.
