# BIB_MERGE_v6_LOG — consolidation of bib hygiene patches into REFERENCES_v6.bib

Task: `bib-merge-v6`
Agent: `agent-224`
Date: 2026-05-18

## Source patches reviewed

| Patch file | Verdict in patch | Action taken |
|---|---|---|
| `BIB_PATCH_d_vollger_dedup.md` | NOT-DUPLICATE (Vollger2023 = Science 10.1126/science.abj6965; concerted_evolution_nahr_Vollger2023 = Nature 10.1038/s41586-023-05895-y) | No edit. Both entries retained as in v5. No sed script run. The patch suggests a future rename Vollger2023 -> segdup_t2t_Vollger2023; deferred (out of scope of this merge). |
| `BIB_PATCH_d_smith1976.bib` (`smith1976crossover`) | Patch added a stanza for Smith G P 1976 Science 191(4227):528-535, DOI 10.1126/science.1251186 | **SKIPPED** — paper already present in REFERENCES_v5.bib under bibkey `Smith1976` with identical DOI. Adding `smith1976crossover` would create a duplicate-paper-under-different-key (the d-bib2-add agent verified only that the *key* was unique, not that the paper was absent). |
| `BIB_PATCH_d_cech2004.bib` (`cech2004chromosomeend`) | New stanza for Cech 2004 Cell 116(2):273-279, DOI 10.1016/s0092-8674(04)00038-8 | ADDED. Confirmed absent from v5 (no matching DOI or PMID). |
| `BIB_PATCH_f20_dux4_cancer.bib` (`dux4_cancer_Chew2019`) | New stanza for Chew et al. 2019 Dev Cell 50(5):658-671, DOI 10.1016/j.devcel.2019.06.011 | ADDED. Confirmed absent from v5 (no matching DOI or PMID). |

## Pre-merge duplicate detection

`grep` against REFERENCES_v5.bib for each patch's DOI:

```
DOI 10.1126/science.1251186 (smith1976crossover)
  -> @article{Smith1976, ... doi = {10.1126/science.1251186}}  HIT (duplicate paper)

DOI 10.1016/s0092-8674(04)00038-8 (cech2004chromosomeend)
  -> no hits

DOI 10.1016/j.devcel.2019.06.011 (dux4_cancer_Chew2019)
  -> no hits

DOI 10.1126/science.abj6965 (Vollger2023)
  -> @article{Vollger2023, ... doi = {10.1126/science.abj6965}}  retained
DOI 10.1038/s41586-023-05895-y (concerted_evolution_nahr_Vollger2023)
  -> @article{concerted_evolution_nahr_Vollger2023, ... doi = {10.1038/s41586-023-05895-y}}  retained
```

## Edits applied

1. Copied `REFERENCES_v5.bib` (372 entries) -> `REFERENCES_v6.bib`.
2. Appended `cech2004chromosomeend` stanza from `BIB_PATCH_d_cech2004.bib`.
3. Appended `dux4_cancer_Chew2019` stanza from `BIB_PATCH_f20_dux4_cancer.bib`.
4. Re-sorted entries alphabetically (case-insensitive, `LC_ALL=en_US.UTF-8 sort -f`) by bibkey, matching v5's sort convention. This places:
   - `cech2004chromosomeend` between `Caspersson1970` and `Cechova2025` (positions 68-70 in v6).
   - `dux4_cancer_Chew2019` between `DuretGaltier2009` and `dux4_d4z4_fshd_degreef2010` (positions 115-117 in v6).
5. Rewrote header block to reflect v6 provenance.
6. No sed mass-rename run against drafts: Vollger NOT-DUPLICATE so the rename script in `BIB_PATCH_d_vollger_dedup.md` is intentionally NOT executed.

## Entry count reconciliation

| Step | Count |
|---|---|
| REFERENCES_v5.bib (baseline) | 372 |
| + cech2004chromosomeend | +1 |
| + dux4_cancer_Chew2019 | +1 |
| + smith1976crossover (would have been) | +0 (skipped, duplicate) |
| - Vollger dedup | -0 (NOT-DUPLICATE) |
| **REFERENCES_v6.bib total** | **374** |

Task description predicted 374 or 375; 374 is within the predicted range (the lower bound corresponded to "Vollger applied + Smith added"; here neither Vollger dedup applied nor Smith newly added, but Smith's de-facto duplicate makes 374 the right count).

## Validation results

```
$ python3 -c 'import pybtex.database; db = pybtex.database.parse_file("paper_prep/synthesis/REFERENCES_v6.bib"); print(len(db.entries))'
374

$ grep '^@' paper_prep/synthesis/REFERENCES_v6.bib | awk -F'[{,]' '{print $2}' | sort | uniq -d
(empty — zero duplicate bibkeys)

$ grep '^@' paper_prep/synthesis/REFERENCES_v6.bib | wc -l
374
```

## RENDERED_REFERENCES_v6.md

`RENDERED_REFERENCES_v5.md` (73 numbered entries) carried forward verbatim to `RENDERED_REFERENCES_v6.md`. Header pointers updated (`source_bibtex`, `draft`, `notes`). No renumbering needed:
- `cech2004chromosomeend` and `dux4_cancer_Chew2019` are not cited anywhere in `NATURE_DRAFT_v5.md` (verified by grep — zero hits).
- Vollger dedup did not execute (NOT-DUPLICATE).
- `smith1976crossover` was not added, so no impact on existing `Smith1976` (entry [52] in v5 rendered).

If `NATURE_DRAFT_v6` introduces citations to the 2 new keys, the renderer should append them with new numbers (>73) and re-sort.

## Files touched in this commit

- `paper_prep/synthesis/REFERENCES_v6.bib` (new)
- `paper_prep/synthesis/RENDERED_REFERENCES_v6.md` (new)
- `paper_prep/synthesis/BIB_MERGE_v6_LOG.md` (new — this file)

No `v1-v5` historical drafts modified.

## Open follow-ups (not actioned here)

- **Vollger key disambiguation**: the bare `Vollger2023` key is now ambiguous because two distinct Vollger 2023 papers (Science vs Nature) coexist in the bib. `BIB_PATCH_d_vollger_dedup.md` recommends renaming `Vollger2023` -> `segdup_t2t_Vollger2023` and provides a safe two-pass sed. Not done here because it requires touching synthesis/ and lit_review/ drafts beyond this task's scope. Recommend a separate `bib-rename-vollger-v6` task if the next paper draft uses both keys.
- **`Smith1976` note enrichment**: existing `Smith1976` stanza lacks a `note` field. The `smith1976crossover` patch contributed useful context ("foundational paper establishing unequal crossover ... cited in P12 alongside Dover 1982"). If the v6 draft needs that semantic anchor, a follow-up could add a `note` field to the existing `Smith1976` entry rather than re-adding the duplicate.
