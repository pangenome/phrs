# Tool versions and pins

Purpose: one place that records **which binary produced which result**. Two reasons this exists:
(1) the Code availability statement the editor requires; (2) we have already been bitten once by a
version living only in prose (sweepga), and once by a stale-index/CLI drift (impg).

Status: living document. Update it whenever a new tool or version enters a result that we quote.

## Alignment / pangenome graph tools

| tool | pinned in Methods | what is actually installed | note |
|---|---|---|---|
| **wfmash** | v0.23.0 | — | used for the all-vs-all flank alignments that seed the PHR call |
| **sweepga** | **v0.1.1** | `sweepga 0.1.1` from `github.com/pangenome/SweepGA#02e7ed0d` | **matches the pin exactly**; candidate calls at FastGA frequency 32; validation passes at 512–8192 with `-s200` |
| **FastGA** | not separately pinned (used via sweepga) | `~/.cargo/bin/FastGA`, no version string in the binary | noted in Methods as sweepga's bundled aligner; for our own validations we record the flags (`-f512/-f8192 -s200 -c70 -i0.6`) rather than a version, which is imperfect and should be stated as such |
| **pggb** | v0.7.4 (commit 498c5d7) with seqwish v0.7.11, smoothxg v0.8.2, odgi v0.9.2, GFAffix v0.2.1, vg v1.62.0 | — | the graph that produced the published 15,668-PHR call |
| **IMPG** | **v0.4.1** | `impg 0.5.0` = `v0.5.0-552-gbe61a4e` (dev tree, 552 commits past the tag) | **CLI has diverged**: the dev build requires `-d/--merge-distance` and `--target-range/--target-bed`, and demands `--alignment-files` even in index mode, so it cannot run the published command lines. Pinned build: `/home/erikg/impg-v0.4.1/` (worktree at tag `v0.4.1` = `8dca5b6`, version 0.4.1), release binary at `target/release/impg` |

## Data products and their provenance

| product | path | provenance | caveat |
|---|---|---|---|
| IMPG index (all-vs-all, whole flank set) | `/moosefs/guarracino/HPRCv2/PHR_III/hprc25272.bidirectional.impg` (37 GB) | built Jan 2026; magic `IMPGIDX2` (format v2) from the 18,827 all-vs-all PAFs | the PAF list (`all-vs-all.1Mb.p95.paf.list`) points at dead `/scratch` paths, so **the index probably cannot be rebuilt**; it is bounded by the 95% identity threshold used for the PHR call, so it cannot speak about copies below 95% |
| HPRC r2 hifi assemblies | `/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/...` and S3 | HPRC release 2 | — |
| PCLAI local ancestry | S3 `release2_pclai_beds`, three coordinate releases: `asm_coord`, `chm13_coord`, `grch38_coord` | HPRC r2 annotation workgroup | **use `asm_coord`/`chm13_coord`, not `grch38`**: GRCh38 has less subtelomeric sequence, so grch38-coordinate coverage (100% of the terminal window) overstates resolvability; chm13-coord retains ~55% of a flank, asm-coord covers 83–100% of the terminal 500 kb at confidence 885–1000. Only HG06807 (PAN027) is in r2 — the rest of the pedigree is absent |
| RepeatMasker per haplotype | index: `/moosefs/.../repeat_masker_bed_hprc_r2_v1.0.index.csv` (462 haplotypes) | HPRC r2 annotation; S3 → HTTPS | this replaces the minor-satellite *proxy* used in the mouse analysis; for human satellite/repeat masking we should use these real annotations |
| mouse meiotic Hi-C | `data/mouse_<stage>_phr_{20000,50000}bp_seqlevel.tsv` | Zuo 2021 (B6 spermatocytes) | satellite annotation for mouse is partial: major satellite real (Dfam `GSAT_MM`), minor satellite proxied |

## Rules we apply

1. One version per comparison. Never mix the dev impg and the pinned impg inside one analysis.
2. Any number we quote states the binary that produced it (path + `--version` + `git describe` where
   relevant).
3. A version that lives only in prose is a bug: it goes in this file.
4. If a tool cannot read our data product in its pinned version, we say so and choose explicitly —
   we do not silently switch and describe it as the pinned analysis.
