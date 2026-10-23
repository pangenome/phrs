# M8 — Referee #1 major point 8: the 500 kb flank right-censors the PHR length distribution

Status: IN PROGRESS. Append-only log; each numbered step records commands run and what
they measured. Numbers are copied from command output, never estimated.

Repo: /home/erikg/phrs (read-only task: only the three deliverable files are created).
Work dir for intermediates: `/scratch/erikg/m8-censoring/` (own directory; not a
running job's scratch).

## Deliverables
- `scripts/cladistics/flank_censoring_sensitivity.py`
- `scripts/cladistics/flank_censoring_sensitivity.tsv`
- `response/notes/M8_flank_censoring.md` (this file)

## Authority and scope
- `PEER_REVIEW.md` Referee #1 major point 8 (read at `PEER_REVIEW.md:87`).
- `response/POINT_BY_POINT.md:129-132` (R1-M8, L0, OPEN).
- `scripts/phr-detection/PIPELINE.md` (stage 2 flank extraction, stage 3 wfmash,
  stage 4 caller), `scripts/phr-detection/find-multichr-regions-incremental.py`
  (caller semantics: telomere-anchored walk, 5 kb windows, `--max`,
  `--min-consecutive-stop 4`, region = first→last passing window).
- `submission/paper.tex:153-165`, `:382`, `:470` (claims under review). Not modified.

## Log

- [step 0] notes skeleton created.
- [step 1] **Inputs inventoried** (all reads succeed; no moosefs hang).
  - `data/all-vs-all.1Mb.p95.id95.len.tsv`: 18,827 rows (one per telomere-anchored
    500 kb flank) + header. Columns `seq, arm, self_chr, region_start, region_end,
    chrs_involved, arms_involved`. `region_*` are 0-based offsets inside the 500 kb
    flank; `.` = no PHR called. 15,668 rows carry a PHR, 3,158 rows are `.`.
  - `/moosefs/guarracino/HPRCv2/PHR_III/hprcv2.1Mb.telo_500kb_trimmed.fa.gz.fai`:
    18,827 flanks, **every one exactly 500,000 bp** (verified: `Counter` of fai
    field 2 = `{500000: 18827}`). 8,771 `parm` + 10,056 `qarm`.
  - `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.per_window.tsv`
    (93,149,310 B, md5 not needed): 521,456 window rows, 18,827 flanks,
    `max dist_from_telomere = 495000`, `max window_end = 500000`. **There is no
    per-window information beyond the 500 kb flank in the published run.**
  - `find` for human 2 Mb / 4 Mb / 5 Mb flank sets: the only such directories are
    `mouse_T2T/{subtelo,flanking,similarity,community_analysis}_{2,4}Mb`. **No human
    2–5 Mb flank alignment exists.** The 500 kb run's own PAFs
    (`PHR_III/all-vs-all.1Mb.p95/`, 37,655 files) have query length 500000
    (e.g. `HG00126#2#CM090128.1:134736978-135236977_chr10_qarm 500000 243000 261000 …`),
    i.e. the alignments themselves were computed on the trimmed 500 kb sequences, so
    the flank cannot be extended from the existing alignment set.
  - `hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz` is *not* 1 Mb flanks: its fai has 15,668
    entries whose lengths are the PHR lengths (e.g. `CHM13#0#chr10:134574995-134754994_chr10_qarm 180000`).
  - Deeper pre-existing similarity source that *is* genome-wide (i.e. not flank-limited):
    `/moosefs/guarracino/HPRCv2/aln-ba8599f934840b4b1a836e83fee938d9c0bc5e77/`
    = 25,274 `*.p95.Pinf.aln.paf.gz` (297 GB) whole-assembly-vs-whole-assembly PAFs
    (`head /scratch/guarracino/hprcv2-paf-list.txt` is the same list under a dead
    `/lizardfs` prefix; the `/moosefs` copy is present and readable). 25,274 ≈ 2 ×
    12,636 unordered pairs of 465 assemblies = the paper's "~12 % of haplotype-pair
    comparisons" survey set. **This is the only deeper source; it was used as the
    second aligner regime.**
  - `data/hprc25272.CHM13.w100kb-xm5-id098-l5k.tsv.gz` (31,202 windows, 100 kb,
    id ≥ 0.98, min output length 5 kb, *transitive*) and
    `data/hprc25272-wf.CHM13.100kb-xm5-id098-l50000.tsv.gz`: the authors' own
    genome-wide per-window tables, used as a third regime.
  - wfmash v0.23.0 (the paper's binary, `/moosefs/guarracino/pggb_wfmash023/...`) does
    not run on this host: `error while loading shared libraries: libgsl.so.27`. Only
    libgsl.so.25 is installed. Local `wfmash` is `v0.24.2-16-g4a35bfc5`.
    One timing probe of the published arrangement (500 kb flank set as queries,
    `--target-prefix CHM13#0#chr10:134254995-134754994_chr10_qarm`, `-t 96`):
    88,781 PAF lines in **6:38 wall, 21.1 GB peak RSS**. At that rate the original
    18,827-job all-vs-all is ~2,070 h single-stream on this host, and a 2 Mb
    re-extraction of the **1,087** saturated flanks (`1,087 × ~26 min` at ~2
    concurrent jobs) is ~230 h. **A 2–5 Mb wfmash re-run is not feasible in this
    task window**, so the extension is measured from the pre-existing genome-wide
    alignment set instead (step 5+), and bounded from the published per-window data
    (step 4). This is recorded as a limitation, not papered over.
