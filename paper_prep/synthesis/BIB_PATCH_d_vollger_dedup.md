# BIB_PATCH: d-bib-doi — Vollger2023 vs concerted_evolution_nahr_Vollger2023

## Entry 1: `Vollger2023`

```bibtex
@article{Vollger2023,
  author  = {Vollger, Mitchell R. and Guitart, Xavi and Dishuck, Philip C.
             and Mercuri, Ludovica and Harvey, William T. and Gershman, Ariel and
             Diekhans, Mark and Sulovari, Arvis and Munson, Katherine M. and
             Lewis, Alexandra P. and Hoekzema, Kendra and Porubsky, David and
             Li, Ruiyang and Nurk, Sergey and Koren, Sergey and Miga, Karen H.
             and Phillippy, Adam M. and Timp, Winston and Ventura, Mario and
             Eichler, Evan E.},
  title   = {Segmental duplications and their variation in a complete human
             genome},
  journal = {Science},
  year    = {2023},
  volume  = {376},
  number  = {6588},
  pages   = {eabj6965},
  doi     = {10.1126/science.abj6965},
  note    = {T2T-CHM13 segmental duplication catalogue; segmental duplications
             at chromosome ends and acrocentric short arms reach >50\% per-arm
             coverage; complements Bailey 2002 and Eichler 2001 frameworks
             at T2T scale}
}
```

**DOI:** `10.1126/science.abj6965`
**Journal:** Science, vol. 376, no. 6588, eabj6965

## Entry 2: `concerted_evolution_nahr_Vollger2023`

```bibtex
@article{concerted_evolution_nahr_Vollger2023,
  author    = {Vollger, Mitchell R. and Dishuck, Philip C. and Harvey, William T. and DeWitt, William S. and Guitart, Xavier and Goldberg, Michael E. and Rozanski, Allison N. and others},
  title     = {Increased mutation and gene conversion within human segmental duplications},
  journal   = {Nature},
  year      = {2023},
  volume    = {617},
  number    = {7960},
  pages     = {325--334},
  doi       = {10.1038/s41586-023-05895-y},
  url       = {https://www.nature.com/articles/s41586-023-05895-y}
}
```

**DOI:** `10.1038/s41586-023-05895-y`
**Journal:** Nature, vol. 617, no. 7960, pp. 325–334

## DOI Comparison

| Key | DOI | Journal | Title |
|-----|-----|---------|-------|
| `Vollger2023` | `10.1126/science.abj6965` | Science 2023 | Segmental duplications and their variation in a complete human genome |
| `concerted_evolution_nahr_Vollger2023` | `10.1038/s41586-023-05895-y` | Nature 2023 | Increased mutation and gene conversion within human segmental duplications |

DOIs are **different**. These are two distinct peer-reviewed papers, both by Vollger et al. in 2023, published in different journals (Science vs Nature) with different subjects (T2T segdup catalogue vs mutation/gene conversion rates within segdups).

## Verdict: NOT-DUPLICATE

Both entries must be retained. They describe different aspects of segmental duplications:
- `Vollger2023` (Science): T2T-CHM13 segmental duplication catalogue — structural inventory.
- `concerted_evolution_nahr_Vollger2023` (Nature): elevated mutation and gene conversion rates inside segdups — evolutionary dynamics.

## Usage Counts

Searched in `paper_prep/synthesis/` and `paper_prep/lit_review/`:

| Key | Occurrences |
|-----|-------------|
| `Vollger2023` (bare, excluding the `concerted_evolution_nahr_` variant) | 39 |
| `concerted_evolution_nahr_Vollger2023` | 33 |

## Suggested Key Rename for `Vollger2023`

The bare key `Vollger2023` is ambiguous now that two Vollger 2023 papers exist in the bib. Recommend renaming it to `segdup_t2t_Vollger2023` in a future consolidation task to match the descriptive-key convention already used for the Nature entry.

If the consolidation task implements this rename, the sed script is:

```bash
# Rename Vollger2023 -> segdup_t2t_Vollger2023 everywhere in synthesis/ and lit_review/
# WARNING: must NOT match concerted_evolution_nahr_Vollger2023 — use word-boundary anchors
find paper_prep/synthesis/ paper_prep/lit_review/ -type f \
  | xargs sed -i 's/\bconcerted_evolution_nahr_Vollger2023\b/PLACEHOLDER_NAHR/g; s/\bVollger2023\b/segdup_t2t_Vollger2023/g; s/PLACEHOLDER_NAHR/concerted_evolution_nahr_Vollger2023/g'
```

This two-pass approach (protect the longer key first, rename the shorter, restore) avoids partial matches.

No sed script is required for deduplication since the entries are NOT duplicates.
