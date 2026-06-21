# Whole-genome joint-parent sweepGA for Fig5 pedigree events

This package is the corrected direct sweepGA experiment for the WashU Fig5
pedigree events. It uses full whole-genome assembly FASTAs from:

- `/moosefs/pangenomes/washu_pedigree/PAN010.fa.gz`
- `/moosefs/pangenomes/washu_pedigree/PAN011.fa.gz`
- `/moosefs/pangenomes/washu_pedigree/PAN027.fa.gz`
- `/moosefs/pangenomes/washu_pedigree/PAN028.fa.gz`

The canonical `/moosefs/pangenomes/washu_pedigree/*.fa.gz` copies returned
BGZF/gzip input/output errors on Slurm compute nodes during retry. The executed
workflow therefore stages equivalent full WashU v1.1 diploid assemblies under
`/moosefs/erikg/phrs/recovery/fig5-whole-genome-joint-parent-sweepga/`, rebuilt
from the public WashU v1.1 URLs with consistent `SAMPLE#hap#chr` headers,
bgzip compression, and `samtools faidx` indexes. `config/comparisons.tsv` points
all four samples at those recovered full-genome paths.

It does not use
`/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/washu.1Mb.telo_500kb_trimmed.fa.gz`.
The earlier 500 kb telomeric-window packages are retained only as historical
diagnostics because their PAF coordinates are local to trimmed subtelomeric
windows and cannot answer the whole-genome recombination-path question.

## Sequence naming and haplotype choices

The canonical source `.fai` indexes use sample/haplotype-prefixed chromosome
records such as `PAN027#1#chr1.maternal` and `PAN027#2#chr1.paternal`. The
recovered full-genome copies use the consistent records requested for this
retry, for example `PAN027#1#chr1` and `PAN027#2#chr1`. The preparation step
selects all full-chromosome records matching the requested haplotype prefix and
records the first/last selected sequence plus counts in
`summaries/input_manifest.tsv`.

Required joint comparisons:

| comparison_id | query haplotype | combined parental target |
| --- | --- | --- |
| `PAN027pat_vs_PAN011_joint` | `PAN027#2` paternal full-genome records | `PAN011#1` + `PAN011#2` full-genome records collapsed to `PAN011#joint` |
| `PAN027mat_vs_PAN010_joint` | `PAN027#1` maternal full-genome records | `PAN010#1` + `PAN010#2` full-genome records collapsed to `PAN010#joint` |
| `PAN028mat_vs_PAN027_joint` | `PAN028#1` maternal full-genome records | `PAN027#1` + `PAN027#2` full-genome records collapsed to `PAN027#joint` |

The target headers are collapsed to one joint parent group before sweepGA so
`1:1`, `1:many`, `2:many`, and `4:many` choose jointly across both parental
haplotypes. With the recovered `SAMPLE#hap#chr` naming, the preparation script
keeps target record names unique by moving the original haplotype into the
joint contig name, e.g. `PAN011#1#chr1` becomes `PAN011#joint#h1_chr1` and
`PAN011#2#chr1` becomes `PAN011#joint#h2_chr1`.

## Why chopping precedes filtering

The raw whole-genome sweepGA output is preserved as `raw_paf/*.paf.gz`, but the
joint filters are run on chopped whole-genome-derived PAF fragments in
`chopped_paf/*.paf.gz`. This correction is intentional: the similarity/path
metric can be distorted when raw alignment intervals merge too far together.
Chopping bounds each raw PAF row into deterministic query-axis fragments before
running the sweepGA mapping filters.

No `rustybam` executable was available on `PATH` in this worktree, so this
package uses `scripts/chop_paf.py`, a deterministic splitter. The default
parameters are:

- `PAF_CHOP_LENGTH=500000`
- `PAF_CHOP_OVERLAP=0`

For each raw PAF row longer than the chop length, the script emits consecutive
query-axis fragments no longer than 500 kb, linearly interpolates target
coordinates, and scales match/alignment counts. The exact command, parameters,
raw input, chopped output, and record counts are written to
`summaries/chop_manifest.tsv`.

## Rerun commands

Run from the repository root. Heavy extraction/alignment/chopping/filtering is
submitted through Slurm. Extraction, manifest creation, and chopping use normal
node-local staging as needed. sweepGA/FastGA raw-alignment graph/database
temporaries and chopped-PAF filter temporaries are explicitly placed under
`/dev/shm` by default via `SWEEPGA_DEVSHM_BASE=/dev/shm`, with cleanup traps
removing each job scratch directory at exit. The raw alignment runner also
copies the prepared full whole-genome query and joint-parent target FASTAs into
the per-job `/dev/shm` directory before invoking FastGA, because FastGA/FAtoGDB
creates source-adjacent database files. This is full-genome input staging, not
a telomeric-window, arm-only, or per-chromosome substitute.

```bash
bash paper_prep/_brainstorming/pedigree_whole_genome_sweepga_joint_parent/scripts/submit_prepare_inputs.sh
bash paper_prep/_brainstorming/pedigree_whole_genome_sweepga_joint_parent/scripts/submit_raw_many_many.sh
bash paper_prep/_brainstorming/pedigree_whole_genome_sweepga_joint_parent/scripts/submit_chop_matrix.sh
bash paper_prep/_brainstorming/pedigree_whole_genome_sweepga_joint_parent/scripts/submit_filter_matrix.sh
bash paper_prep/_brainstorming/pedigree_whole_genome_sweepga_joint_parent/scripts/validate_outputs.sh
```

The stages may also be run under explicit Slurm dependencies if rerunning the
whole workflow from scratch.

## Outputs

First-class alignment artifacts:

- `raw_paf/*.paf.gz`: raw whole-genome `many:many` sweepGA outputs.
- `chopped_paf/*.paf.gz`: chopped fragments derived from the raw whole-genome
  PAFs and used as filter input.
- `filtered_paf/*.paf.gz`: chopped-input joint filters for `many:many`,
  `1:1`, `1:many`, `2:many`, and `4:many`.

Summary/provenance tables:

- `summaries/input_manifest.tsv`
- `summaries/recovery_manifest.tsv`
- `summaries/slurm_jobs.tsv`
- `summaries/chop_manifest.tsv`
- `summaries/filter_manifest.tsv`
- `summaries/paf_file_summary.tsv`
- `summaries/output_file_manifest.tsv`

Strict `1:1` is retained as a diagnostic control only. The chopped raw
`many:many` layer and chopped `4:many` layer are the likely evidence layers for
downstream review. This task intentionally does not generate a final Fig5
schematic; evidence selection is left to the downstream review task.

## Slurm provenance

The first preparation submission, job `1704277`, failed before touching data
because Slurm staged the shell script under `/var/spool/slurmd` and the run
script resolved helper paths relative to that staged copy. The run scripts were
patched to use the exported repository `PACKAGE_DIR`.

The current Slurm job tables are in `summaries/prepare_slurm_jobs.tsv`,
`summaries/slurm_jobs.tsv`, `summaries/chop_slurm_jobs.tsv`, and
`summaries/filter_slurm_jobs.tsv` as applicable.

## Current execution status

The first attempt could not read the canonical PAN011 full assembly from
`/moosefs/pangenomes/washu_pedigree/PAN011.fa.gz` on Slurm compute nodes.
Retry preparation then exposed the same failure class for canonical PAN027.
This completed workflow rebuilds/stages full WashU v1.1 bgzip+faidx copies for
all four required samples in
`/moosefs/erikg/phrs/recovery/fig5-whole-genome-joint-parent-sweepga/` and uses
those recovered full-genome sources for every comparison. No 500 kb
telomeric-window fallback, per-chromosome-only fallback, or arm/window
substitute is used.

Completed Slurm stages:

- Recovery: jobs `1704287` and `1704288`, recorded in
  `summaries/recovery_slurm_jobs.tsv` and `summaries/recovery_manifest.tsv`.
- Full-genome FASTA preparation: job `1704290`, recorded in
  `summaries/prepare_slurm_jobs.tsv` and `summaries/input_manifest.tsv`.
- Raw whole-genome `many:many` sweepGA alignments: jobs `1704307`-`1704309`,
  recorded in `summaries/slurm_jobs.tsv`.
- Deterministic PAF chopping: job `1704310`, recorded in
  `summaries/chop_slurm_jobs.tsv` and `summaries/chop_manifest.tsv`.
- Chopped-input joint filters: job `1704311`, recorded in
  `summaries/filter_slurm_jobs.tsv` and `summaries/filter_manifest.tsv`.

`scripts/validate_outputs.sh` passed after these jobs. Raw, chopped, and
filtered PAF gzip integrity checks passed for all required comparisons.
