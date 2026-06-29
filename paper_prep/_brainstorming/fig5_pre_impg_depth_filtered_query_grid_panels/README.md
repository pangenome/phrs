# Fig5 pre-IMPG depth-filtered query-grid panels

Query-grid style visualization of the completed one-node pre-IMPG pilot:

- input: `fig5_pre_impg_depth_filtered_similarity` class-winner IMPG similarity
  TSVs, with one best same-chromosome and one best interchromosomal match per
  2 kb query window
- comparison: `PAN027pat_vs_PAN011_joint`
- windows: PAR1 control and PAN027 chr9q -> chr3q candidate from the existing
  Fig5 query-grid panel config
- rows: SweepGA f32 `1:1`, `4:4`, and `10:10` mapping bases after pre-IMPG
  query-window depth filtering; each window is colored by whether the
  same-chromosome/homolog or interchromosomal candidate wins

Build:

```bash
bash paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_query_grid_panels/scripts/make_panels.sh
```
