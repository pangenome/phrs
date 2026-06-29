# Fig5 pre-IMPG depth-filtered query-grid panels

Query-grid style visualization of the completed one-node pre-IMPG pilot:

- input: `fig5_pre_impg_depth_filtered_similarity` top-N IMPG similarity TSVs
- comparison: `PAN027pat_vs_PAN011_joint`
- windows: PAR1 control and PAN027 chr9q -> chr3q candidate from the existing
  Fig5 query-grid panel config
- rows: SweepGA f32 `1:1`, `4:4`, and `10:10` mapping bases after pre-IMPG
  query-window depth filtering

Build:

```bash
bash paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_query_grid_panels/scripts/make_panels.sh
```

