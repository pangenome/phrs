// Snippet for fan-in only. Do not apply in this task.
//
// Goal: remove the detached bottom legend from slide 13b while keeping the
// detailed pedigree panels. The event labels are already printed directly in
// the tracks, so the bottom legend is not needed for talk-speed reading.

#figure-slide(
  "13b",
  "Backup: detailed pedigree exchange events",
  "../_revision_assets/v8/typography_legend_cleanup/slide13b_pedigree_bottom_no_unused_legend.png",
  source: "v8/typography_legend_cleanup/crop_png_top.py; materialized crop of s13_pedigree_bottom.png with unused bottom legend removed; event labels are direct in-panel annotations",
)
