# Pedigree Native odgi BEDPE/PAF Merge-Distance Decision Record

Primary runnable: `scripts/pedigree/run_untangle_native_merge_tracts.sbatch`.
Parser: `scripts/pedigree/untangle_native_merge_tracts.py`.

Large native odgi intermediates are outside git under:

- `/moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm`

## Were native BEDPE/PAF outputs generated on Slurm?

Yes. Slurm job `1703959` generated the valid rerun outputs in
`/moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm` using 24 CPUs, 96G, and the `workers` partition. The first
worker also ran a direct head-node pass before the Slurm-only constraint was
issued; that earlier pass is not treated as the valid provenance. The committed
sbatch script is now the natural runnable and directly runs both BEDPE and PAF
forms of `odgi untangle`, including `-p` for PAF.

The Slurm grid was three WashU child-parent comparisons, merge distances `0`
and `1000`, `-e 50000`, `-j 0.8`, and `-n 4`.
The default compact summary parses the `m1000` files because that is sufficient
for the BEDPE/PAF/sweepGA comparison; pass `--merge-dist 0 1000` for a full
native merge-distance audit. The parser projects odgi's emitted `nb`/`nth.best`
fields to analysis top-N values 1, 2, and 4 without rerunning heavy work.

Command manifest:

```bash
/home/erikg/.guix-profile/bin/odgi untangle -t ${SLURM_CPUS_PER_TASK:-24} -i /moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/pggb/washu.1Mb.telo_500kb_trimmed.fa.gz.6e0e250.11fba48.13f423a.smooth.final.og -Q /moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm/queries_PAN027_vs_PAN010.txt -R /moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle/targets_PAN027_vs_PAN010.txt -e 50000 -m 1000 -j 0.8 -n 4 | gzip -c > /moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm/PAN027_vs_PAN010.e50000.m1000.j0.8.n4.bedpe.gz
/home/erikg/.guix-profile/bin/odgi untangle -t ${SLURM_CPUS_PER_TASK:-24} -i /moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/pggb/washu.1Mb.telo_500kb_trimmed.fa.gz.6e0e250.11fba48.13f423a.smooth.final.og -Q /moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm/queries_PAN027_vs_PAN010.txt -R /moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle/targets_PAN027_vs_PAN010.txt -e 50000 -m 1000 -j 0.8 -n 4 -p | gzip -c > /moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm/PAN027_vs_PAN010.e50000.m1000.j0.8.n4.paf.gz
/home/erikg/.guix-profile/bin/odgi untangle -t ${SLURM_CPUS_PER_TASK:-24} -i /moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/pggb/washu.1Mb.telo_500kb_trimmed.fa.gz.6e0e250.11fba48.13f423a.smooth.final.og -Q /moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm/queries_PAN027_vs_PAN011.txt -R /moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle/targets_PAN027_vs_PAN011.txt -e 50000 -m 1000 -j 0.8 -n 4 | gzip -c > /moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm/PAN027_vs_PAN011.e50000.m1000.j0.8.n4.bedpe.gz
/home/erikg/.guix-profile/bin/odgi untangle -t ${SLURM_CPUS_PER_TASK:-24} -i /moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/pggb/washu.1Mb.telo_500kb_trimmed.fa.gz.6e0e250.11fba48.13f423a.smooth.final.og -Q /moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm/queries_PAN027_vs_PAN011.txt -R /moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle/targets_PAN027_vs_PAN011.txt -e 50000 -m 1000 -j 0.8 -n 4 -p | gzip -c > /moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm/PAN027_vs_PAN011.e50000.m1000.j0.8.n4.paf.gz
/home/erikg/.guix-profile/bin/odgi untangle -t ${SLURM_CPUS_PER_TASK:-24} -i /moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/pggb/washu.1Mb.telo_500kb_trimmed.fa.gz.6e0e250.11fba48.13f423a.smooth.final.og -Q /moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm/queries_PAN028_vs_PAN027.txt -R /moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle/targets_PAN028_vs_PAN027.txt -e 50000 -m 1000 -j 0.8 -n 4 | gzip -c > /moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm/PAN028_vs_PAN027.e50000.m1000.j0.8.n4.bedpe.gz
/home/erikg/.guix-profile/bin/odgi untangle -t ${SLURM_CPUS_PER_TASK:-24} -i /moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/pggb/washu.1Mb.telo_500kb_trimmed.fa.gz.6e0e250.11fba48.13f423a.smooth.final.og -Q /moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm/queries_PAN028_vs_PAN027.txt -R /moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle/targets_PAN028_vs_PAN027.txt -e 50000 -m 1000 -j 0.8 -n 4 -p | gzip -c > /moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm/PAN028_vs_PAN027.e50000.m1000.j0.8.n4.paf.gz
```

## Does sweepGA accept/filter odgi-emitted PAF directly?

Yes for an uncompressed odgi-emitted PAF. The installed binary is
`/home/erikg/.cargo/bin/sweepga` reporting `sweepga 0.1.1`, from `/moosefs/erikg/sweepga` commit
`018e4ce49d2c125820e0ac50dc5feaa02d423683`. The required commit for this task is `018e4ce49d2c125820e0ac50dc5feaa02d423683`.

Minimal representative test used
`PAN028_vs_PAN027.e50000.m1000.j0.8.n4.paf.gz` from the Slurm output. Passing
the `.paf.gz` path directly failed with `invalid BGZF header`; after plain
decompression to `.paf`, sweepGA accepted and filtered the native odgi PAF:

- `/home/erikg/.cargo/bin/sweepga --num-mappings 1:many --scaffold-jump 0 --output-file /moosefs/erikg/phrs/pedigree_native_untangle_agent2556_sweepga_test/PAN028_vs_PAN027.e50000.m1000.j0.8.n4.uncompressed.1_many.paf /moosefs/erikg/phrs/pedigree_native_untangle_agent2556_sweepga_test/PAN028_vs_PAN027.e50000.m1000.j0.8.n4.native.paf`: accepted, 4827 PAF rows ([INFO ] [sweepga::done] total:0.1s cpu:0.6x rss:12.80 MB disk_peak:0 bytes disk_written:0 bytes; gzip path rejected: invalid BGZF header)
- `/home/erikg/.cargo/bin/sweepga --num-mappings 2:many --scaffold-jump 0 --output-file /moosefs/erikg/phrs/pedigree_native_untangle_agent2556_sweepga_test/PAN028_vs_PAN027.e50000.m1000.j0.8.n4.uncompressed.2_many.paf /moosefs/erikg/phrs/pedigree_native_untangle_agent2556_sweepga_test/PAN028_vs_PAN027.e50000.m1000.j0.8.n4.native.paf`: accepted, 6998 PAF rows ([INFO ] [sweepga::done] total:0.1s cpu:0.7x rss:13.40 MB disk_peak:0 bytes disk_written:0 bytes; gzip path rejected: invalid BGZF header)
- `/home/erikg/.cargo/bin/sweepga --num-mappings 4:many --scaffold-jump 0 --output-file /moosefs/erikg/phrs/pedigree_native_untangle_agent2556_sweepga_test/PAN028_vs_PAN027.e50000.m1000.j0.8.n4.uncompressed.4_many.paf /moosefs/erikg/phrs/pedigree_native_untangle_agent2556_sweepga_test/PAN028_vs_PAN027.e50000.m1000.j0.8.n4.native.paf`: accepted, 8742 PAF rows ([INFO ] [sweepga::done] total:0.1s cpu:0.6x rss:17.61 MB disk_peak:0 bytes disk_written:0 bytes; gzip path rejected: invalid BGZF header)

The compact comparison is in `scripts/pedigree/untangle_native_merge_summary.tsv`.

## Does native merge or sweepGA filtering justify a later manuscript edit?

No. Native BEDPE/PAF output and sweepGA filtering establish a cleaner
provenance path, but they do not clearly improve tract calls enough to justify
a later manuscript edit from this task alone. For orientation only, native PAF
`m1000/n4` has 1994 parser-counted
interchromosomal tracts and native BEDPE `m1000/n4` has
1979. These are methods/provenance counts, not
a conversion-vs-crossover mechanism claim.
