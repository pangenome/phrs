# PAN011 full-assembly recovery

Task `prereq-restore-readable` found that the canonical MooseFS PAN011 files are
not readable from the login node or from Slurm node `octopus07`:

```text
/moosefs/pangenomes/washu_pedigree/PAN011.fa.gz
/moosefs/pangenomes/washu_pedigree/PAN011.fa.gz.fai
```

Both returned `Input/output error` during direct `gzip -t` and `.fai` reads.
The directory is owned by `guarracino:guarracino` and is not writable by this
agent account, so the files could not be replaced in place without an owner or
administrator repair.

The recovered full whole-genome PAN011 assembly is:

```text
/moosefs/erikg/phrs/recovery/prereq-restore-readable/PAN011.fa.gz
```

This is not the 500 kb telomeric-window FASTA. It was rebuilt from the public
WashU pedigree v1.1 diploid assembly:

```text
https://public.gi.ucsc.edu/~mcechova/pedigree/assemblies/v1.1/assembly.v1.1.PAN011.diploid.fa
```

The upstream headers have the form `PAN011.chr1.haplotype1`. During recovery
they were converted to the local PanSN-style names expected by downstream
workflow code:

```text
PAN011#1#chr1
PAN011#2#chr1
```

The rebuilt FASTA is bgzip-compressed and indexed with `samtools faidx`; its
`.fai` should show 46 full-chromosome records, beginning with `PAN011#1#chr1`
and including the corresponding `PAN011#2#chr*` records.

To use the recovered full assembly with the joint-parent preparation code, set
`WASHU_PEDIGREE_FASTA`:

```bash
WASHU_PEDIGREE_FASTA=/moosefs/erikg/phrs/recovery/prereq-restore-readable/PAN011.fa.gz \
  python3 paper_prep/_brainstorming/pedigree_direct_sweepga_joint_parent/scripts/prepare_inputs.py
```

That override is intentionally explicit. Without it,
`prepare_inputs.py` keeps its historical default telomeric-window FASTA for the
old window-level package.

Validation commands for Slurm compute nodes:

```bash
gzip -t /moosefs/erikg/phrs/recovery/prereq-restore-readable/PAN011.fa.gz
head /moosefs/erikg/phrs/recovery/prereq-restore-readable/PAN011.fa.gz.fai
```

If a future owner/admin repair can replace the canonical MooseFS files, copy the
recovered `PAN011.fa.gz`, `PAN011.fa.gz.fai`, and `PAN011.fa.gz.gzi` to
`/moosefs/pangenomes/washu_pedigree/` and rerun those same checks against the
canonical path on `octopus07`.
