# Fig5 homolog-vs-interchrom telomeric zoom panels

Zoom panels for the PAN027 paternal vs PAN011 father `sweepga_f32` 10:10
pre-IMPG depth-filtered class-winner scan.

The generator keeps every 2 kb query window in the first or last 500 kb of a
query chromosome where the best interchromosomal IMPG similarity ranks above the
best same-chromosome/homolog similarity. Arms are ranked by winning base pairs;
the rendered panel includes the strongest arms and explicitly retains PAR1
(`chrXp`) and the distal `chr9q`/`chr3q` candidate.

Outputs:

- `fig5_homolog_vs_interchrom_zoom_panels.pdf`
- `fig5_homolog_vs_interchrom_zoom_panels.png`
- `fig5_homolog_vs_interchrom_zoom_panels.svg`
- `telomeric_interchrom_winner_segments.tsv`
- `telomeric_interchrom_winner_arm_summary.tsv`
- `artifact_manifest.tsv`

Regenerate from the repository root:

```bash
python3 paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels/scripts/plot_homolog_vs_interchrom_zoom_panels.py
Rscript paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels/scripts/render_homolog_vs_interchrom_zoom_panels.R
```

The script first looks for the input class-winner table in this worktree under
`paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/outputs/`.
If that large upstream output is absent, it falls back to the shared mirror at
`/moosefs/erikg/phrs/.../outputs/`.
