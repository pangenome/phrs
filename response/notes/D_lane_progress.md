# Lane D progress notes (superpopulation descriptives + editor compliance)

> Provenance for this note: see `response/PROVENANCE.md` (tool versions: `response/notes/tool_versions_and_pins.md`).

Created: start of resumed lane D (prior attempt killed at 30 min).

## commands run

## numbers

## what this does NOT establish

### Step 1 — inputs located (all pre-existing; nothing recomputed yet)

- `scripts/cladistics/phr_by_superpopulation.py` (31,281 B) and
  `scripts/cladistics/phr_by_superpopulation.tsv` (11,855 B, 113 lines = header
  + 112 data rows) both exist and are readable.
- TSV sections present, from `cut -f1 | sort | uniq -c`: cohort / superpopulation /
  community / reference (verified below).
- Superpopulation mapping file (repo-local, NOT moosefs):
  `slides/v2-review-zoom/_revision_assets/v9/labels_superpop_stats_polish/nearest_same_superpop_mds_distances.tsv`
  (4,180,581 B, 15,669 lines incl. header; columns include Name, Sample,
  Haplotype, Chromosome, Arm, Superpopulation).
- Its provenance manifest:
  `slides/.../labels_superpop_stats_polish/source_manifest.tsv` lists
  `seq_assignment_tsv = /moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv`
  (sha256 17bb574d…), i.e. the repo copy comes from that upstream assignment table.
- Both moosefs cross-check inputs exist and are readable *today*:
  `/moosefs/guarracino/HPRCv2/data/hprc-sequence-production.tsv` (15,435 B,
  header `ChildID N Sex Subpopulation Superpopulation Production Year`) and
  `/moosefs/.../similarity/hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv`
  (1,410,827 B).
- F_ST comparison script (the differentiation analysis, NOT this one):
  `scripts/popgen/matched_fst_from_blocks.py` (2,666 B) -> writes
  `scripts/popgen/matched_fst_d_m6.tsv`; reuses constants/jackknife from
  `scripts/popgen/matched_fst_d_m6.py`.

## commands run
- `cat scripts/cladistics/phr_by_superpopulation.tsv` (full read)
- `head -120 scripts/cladistics/phr_by_superpopulation.py` (docstring/provenance)
- `grep -n -i "superpop\|moosefs\|labels" scripts/cladistics/phr_by_superpopulation.py`
- `ls -la slides/v2-review-zoom/_revision_assets/v9/labels_superpop_stats_polish/`
- `cat slides/.../labels_superpop_stats_polish/source_manifest.tsv`
- `ls -la /moosefs/guarracino/HPRCv2/data/hprc-sequence-production.tsv` + `head -2`
- `ls -la /moosefs/.../hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv`
- `grep -n -i superpopulation PEER_REVIEW.md` -> line 97 = referee #1 minor point 11
- `head -40 scripts/popgen/matched_fst_from_blocks.py`
- `sed -n 1,80p PEER_REVIEW.md` (editor letter)

### Step 2 — D1 written; reproduction gate passed

- `response/notes/D1_superpopulation_descriptives.md` written (Tables 1-3 all copied from
  the TSV; no cell recomputed by hand).
- Reran the script to a scratch path ONLY (`--out /tmp/erikg-laneD/phr_by_superpopulation.rerun.tsv`)
  so the in-repo TSV was untouched: `md5sum` identical
  (ac0d6116abda54ebcaf159014416ba10 both), `diff` -> IDENTICAL.
- Captured its stdout (gate lines + the two permutation checks) verbatim into D1 §6.
  Permutation stdout: Cramers V=0.0407 chi2=14886869.9 df=56 p_perm=0.0000 (20000 perms,
  seed 20260201); Kruskal-Wallis H=16.245 df=4 p_perm=0.0024.
- Internal-consistency awk gates on the TSV: superpop sums == cohort row exactly
  (232 / 464 / 15630 / 2245670000); per-population community frac_phrs and frac_bp each
  sum to 1.0000 +- 0.0002; 112 data rows.
- Corrected three numbers after checking them: C7 bp share range is 15.2% (AMR) to 19.6%
  (EUR), not "16.1-19.6%"; C7 500kb-censoring share is 91.6-95.8%; C3/C6 track cohort within
  0.7 pp not 0.9; CHM13 reference row n_arms = 38 (not a garbled cell).

## numbers
- TSV md5 `ac0d6116abda54ebcaf159014416ba10`, 113 lines (112 data rows), 11,855 B.
- cohort: 232 individuals / 464 haplotypes / 41 arms / 15,630 PHRs / 2,245,670,000 bp.
- AFR 67/134/4549/656470000; AMR 44/88/2974/428005000; EAS 52/104/3438/480285000;
  EUR 32/64/2168/319460000; SAS 37/74/2501/361450000.
- CHM13 reference section: 1/1/38 arms/38 PHRs/6,840,000 bp, excluded from aggregates.
- label-table cross-check: 0 conflicts vs both moosefs sources; 232/232 covered by the
  per-PHR assignment table; 29 individuals absent from hprc-sequence-production.tsv.

### Step 3 — D2 fact-gathering (editor compliance)

Verified for D2 (each is a command I ran, or a file:line I read):
- editor letter = `PEER_REVIEW.md`:1-33 (DATA AND CODE AVAILABILITY STATEMENTS at :21 and :23,
  ETHICS STATEMENT :27, EXTENDED DATA :31, SOURCE DATA (GRAPHS) :33, LENGTH :19,
  REPRODUCIBILITY :20, THIRD PARTY RIGHTS :35).
- existing combined section: `submission/paper.tex`:889-910 `\subsection*{Data and code availability}`
  (inside Methods; URLs at :892, :894, :899, :905). No ethics statement anywhere in the tex
  (`grep -n -i "consent\|ethic\|IRB\|approval" submission/paper.tex` -> only an affiliation line).
- HPRC v2 release page + BioProject: `submission/bibliography.bib`:734-740
  (`https://humanpangenome.org/hprc-data-release-2/`, note cites BioProject PRJNA730823);
  PRJNA730823 also in `submission/paper.tex`:1165. Pedigree: bib :415-424
  (bioRxiv doi 10.64898/2025.12.14.693655) + `https://github.com/biomonika/washu-pedigree`
  at paper.tex:905.
- Extended Data: exactly 3 items, all `.pdf` (`submission/fig/ExtendedDataFigures/`:
  ED_Fig1_human_contacts.pdf, ED_Fig2_copy_number_ontology_heatmap.pdf,
  ED_Fig3_module_maps.pdf), legends in paper.tex:912-976 (correctly at end of main text,
  not in the ED files). Two unreferenced panels exist in the same folder
  (ED_Fig2a_PAN027_maternal_recombination.png, ED_Fig2b_PAN028_maternal_recombination.png)
  and appear nowhere in the tex.
- `front data/`: 63 entries; `data/all-vs-all.1Mb.p95.id95.len.tsv` present;
  **`all-vs-all.1Mb.p95.id95.len.per_window.tsv` is NOT in the repo** (`ls data | grep per_window`
  -> nothing, exit 1), matching response/notes/R5_...md.
- No REPORTING SUMMARY / Code-and-software checklist file anywhere in the repo
  (`find . -iname "*reporting*" ...` -> empty).
- No LICENSE/COPYING at repo root, no VERSIONS.md (`find -name "VERSIONS*"` -> empty).
- `git remote origin` (read from `.git/config`, no git command run) = `git@github.com:ekg/phrs.git`,
  while paper.tex:894 states `https://github.com/pangenome/phrs`.
- word counts, crude LaTeX-stripped rule (sed strip \cite/\ref/\label/\url/$math$, then wc -w):
  main text lines 95-422 = **3,330**; abstract 65-91 = 228; Methods 598-910 = 3,132;
  figure legends 451-597 = 913; ED legends 912-976 = 393.

### Step 4 — D2 written; both deliverables complete

- `response/notes/D2_editor_compliance.md` written: §0 requirement table (10 editor items mapped to
  PEER_REVIEW.md line numbers), §1 Data availability draft (with the minimum-dataset table naming
  every repo table), §2 Code availability draft + Zenodo/DOI/licence/VERSIONS author actions,
  §3 ethics statement split into [DRAFTED attestable] / [AUTHOR TO SUPPLY] table, §4 format
  checklist (ED .tif/.eps, Source Data per figure Fig1-5, word count, tracked changes,
  Reporting Summary, third-party rights), §5 what it does not establish.
- Every path named in D2 exists: 18-file existence loop + 23 glob hits, all OK
  (`ls`/`[ -e ]` checks). NOTE: the task brief spells `scripts/popogen/`; the real directory is
  `scripts/popgen/` (verified `find scripts -maxdepth 1 -type d`). Both notes use the real path.
- Orphan ED panels found: `ED_Fig2a_PAN027_maternal_recombination.png`,
  `ED_Fig2b_PAN028_maternal_recombination.png` exist in the ED folder and appear nowhere in
  paper.tex -> raised as AMBIGUITY in D2 §4.1 (AUTHOR TO SUPPLY).
- Line-number corrections made after checking: paper.tex 896-906 (re-analysis paragraph),
  95-421 (main text for word count; 422 is blank), ED legends 912-976.

## numbers (final)
- main text 3,330 words (lines 95-421, crude LaTeX-strip rule); abstract 228; Methods 3,132;
  figure legends 913; ED legends 393.
- ED items: 3, all .pdf (limit 10).
- data/ : 63 entries; `all-vs-all.1Mb.p95.id95.len.per_window.tsv` NOT in repo (exit 1 on grep).
- repo has no LICENSE, no COPYING, no VERSIONS.md, no Reporting Summary / code checklist.

## what this does NOT establish
- Nothing was placed into submission/paper.tex (read-only for this lane).
- Ethics facts (committee, approval numbers, consent wording), the archival DOI, the licence
  choice, Reporting Summary and third-party rights table are all AUTHOR TO SUPPLY.
- D1 is descriptive only; differentiation is the matched F_ST table (scripts/popgen/).
- No git operations were run (explicit instruction), so these notes are uncommitted.
