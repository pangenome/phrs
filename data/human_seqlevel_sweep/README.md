# human_seqlevel_sweep

Per-PHR-pair sequence-vs-3D tables for the human Fig 4 panels, at every Hi-C /
Pore-C / CiFi resolution (5 / 10 / 20 / 50 / 100 kb). Same schema as the mouse
sweep (`data/mouse_meiosis_sweep/seqlevel/`): one row per inter-chromosomal PHR
sequence pair, columns

```
seq_a seq_b chr_a arm_a chr_b arm_b size_a size_b jaccard hic_contact_raw hic_contact_norm hic_bins
```

`hic_contact_norm = hic_contact_raw / hic_bins` = mean balanced contact per
bin-pair (a contact *density*, NOT a cooltools O/E).

Datasets (each at 5/10/20/50/100 kb):
- `human_HG002_porec_*bp_seqlevel.tsv`  — HG002 Pore-C (Fig 4a), 2830 pairs
- `human_HG002_hic_*bp_seqlevel.tsv`    — HG002 Hi-C
- `human_CHM13_hic_*bp_seqlevel.tsv`    — CHM13 Hi-C (ED 1), 652 pairs @50 kb
- `human_HG002_cifi_*bp_seqlevel.tsv`   — HG002 CiFi (very sparse, ~90% zero)

The 50 kb HG002 Pore-C and CHM13 Hi-C files match the panels already vendored at
repo root (`data/human_*_50000bp_seqlevel.tsv`).

## Source / re-fetch (needs /moosefs)

Upstream dir: `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_free/`
Pipeline: `/moosefs/guarracino/HPRCv2/scripts/community/sequence_hic_correlation.py`
(multi-resolution `.mcool` coolers in `.../PHR_III/HiC/<sample>/`).

```
B=/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_free
for ds in HG002_porec HG002_hic CHM13_hic HG002_cifi; do
  scp "HOST:$B/human_${ds}_[0-9]*bp_seqlevel.tsv" data/human_seqlevel_sweep/
done
```

## Resolution evaluation

`scripts/human/human_seqlevel_resolution.R` — cross-PHR ρ(Jaccard, contact) per
dataset × resolution with the PHR-node bootstrap. The human coupling is
**resolution-invariant** and solid at 50 kb (Pore-C ρ≈0.38, Hi-C ≈0.66,
CHM13 ≈0.72; CiFi weak due to sparsity). 50 kb is the canonical paper setting
(see `paper_prep/_brainstorming/mouse_bouquet_reanalysis/NOTE.md`).
