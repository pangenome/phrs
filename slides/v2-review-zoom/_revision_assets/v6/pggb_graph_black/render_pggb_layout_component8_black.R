#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
cmd_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "render_pggb_layout_component8_black.R"

value_after <- function(flag, default = NA_character_) {
  hit <- which(args == flag)
  if (length(hit) == 0 || hit[length(hit)] == length(args)) {
    return(default)
  }
  args[hit[length(hit)] + 1]
}

display_path <- function(path) {
  abs_path <- normalizePath(path, mustWork = FALSE)
  cwd <- normalizePath(getwd(), mustWork = TRUE)
  cwd_prefix <- paste0(cwd, .Platform$file.sep)
  if (startsWith(abs_path, cwd_prefix)) {
    return(sub(cwd_prefix, "", abs_path, fixed = TRUE))
  }
  abs_path
}

shell_command <- function(parts) {
  paste(shQuote(parts), collapse = " ")
}

layout_tsv <- value_after(
  "--layout-tsv",
  "/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.og.lay.tsv"
)
out_png <- value_after("--out", "pggb_graph_2d_black.png")
component_id <- as.integer(value_after("--component", "8"))
render_log <- value_after("--render-log", "render_log.tsv")
point_color <- value_after("--point-color", "#111111")
point_alpha <- as.numeric(value_after("--point-alpha", "0.30"))
point_cex <- as.numeric(value_after("--point-cex", "0.12"))
border_color <- value_after("--border-color", "#b8c0cc")
background_color <- value_after("--background-color", "white")

if (!file.exists(layout_tsv)) {
  stop("Layout TSV does not exist: ", layout_tsv)
}
if (is.na(point_alpha) || point_alpha <= 0 || point_alpha > 1) {
  stop("--point-alpha must be in (0, 1]")
}
if (is.na(point_cex) || point_cex <= 0) {
  stop("--point-cex must be positive")
}

layout <- read.delim(layout_tsv, sep = "\t", stringsAsFactors = FALSE)
required_cols <- c("idx", "X", "Y", "component")
missing_cols <- setdiff(required_cols, names(layout))
if (length(missing_cols) > 0) {
  stop("Layout TSV is missing required columns: ", paste(missing_cols, collapse = ", "))
}

component <- layout[layout$component == component_id, required_cols]
if (nrow(component) == 0) {
  stop("No rows found for component ", component_id)
}

component_counts <- sort(table(layout$component), decreasing = TRUE)
component_rank <- match(as.character(component_id), names(component_counts))
status <- if (component_rank == 1) "main_component" else paste0("component_rank_", component_rank)

plot_x <- component$Y
plot_y <- component$X

png(out_png, width = 1920, height = 1080, res = 144, bg = background_color, type = "cairo")
par(mar = c(0.18, 0.18, 0.18, 0.18), xaxs = "i", yaxs = "i", family = "sans")
plot(
  plot_x,
  plot_y,
  pch = 16,
  cex = point_cex,
  col = grDevices::adjustcolor(point_color, alpha.f = point_alpha),
  axes = FALSE,
  ann = FALSE,
  asp = 1
)
usr <- par("usr")
rect(usr[1], usr[3], usr[2], usr[4], border = border_color, lwd = 2)
legend(
  "bottomleft",
  legend = c(
    sprintf("PGGB ODGI layout main component (component %s)", component_id),
    sprintf("%s layout nodes; rotated 90 degrees for 16:9 slide fit", format(nrow(component), big.mark = ","))
  ),
  bty = "n",
  text.col = "#1f1f1f",
  cex = 0.82
)
invisible(dev.off())

png_info <- file.info(out_png)

render_record <- data.frame(
  key = c(
    "render_task",
    "derived_from_task",
    "render_mode",
    "status",
    "layout_tsv",
    "source_odgi_graph",
    "source_odgi_layout_binary",
    "prior_slide_asset",
    "inspected_existing_odgi_draw_png",
    "inspected_existing_odgi_viz_png",
    "component",
    "component_nodes",
    "component_rank_by_layout_node_count",
    "x_range",
    "y_range",
    "orientation",
    "palette",
    "background_color",
    "point_color",
    "point_alpha",
    "point_cex",
    "border_color",
    "output_png",
    "output_dimensions",
    "output_bytes",
    "slurm_job_id",
    "render_command"
  ),
  value = c(
    "review-zoom-v6-pggb-graph-black",
    "review-zoom-v5-pggb-gfalook-2d-render",
    "rerender_existing_odgi_layout_tsv_component_with_dark_neutral_palette",
    status,
    normalizePath(layout_tsv, mustWork = TRUE),
    "/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.og",
    "/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.og.lay",
    "slides/v2-review-zoom/_revision_assets/v5/pggb_graph_odgi/pggb_graph_2d.png",
    "/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.og.lay.draw.png (167x1000; too sparse/narrow for slide use)",
    "/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.og.viz_multiqc.png (1968x157164; too tall/raw for slide use)",
    as.character(component_id),
    as.character(nrow(component)),
    as.character(component_rank),
    paste(range(component$X), collapse = "-"),
    paste(range(component$Y), collapse = "-"),
    "ODGI X/Y coordinates plotted as Y/X to rotate the tall layout into a 16:9 slide frame; no re-layout",
    "dark neutral graph marks on a light background for projection readability; no blue graph strokes",
    background_color,
    point_color,
    as.character(point_alpha),
    as.character(point_cex),
    border_color,
    display_path(out_png),
    "1920x1080",
    as.character(png_info$size),
    "not_used_no_heavy_ODGI_extraction_or_render_submitted",
    shell_command(c("Rscript", display_path(script_path), args))
  ),
  stringsAsFactors = FALSE
)

write.table(render_record, render_log, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Wrote ", out_png, "\n", sep = "")
cat("Wrote ", render_log, "\n", sep = "")
