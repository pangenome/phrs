#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
cmd_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "render_pggb_layout_component8.R"

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

layout_tsv <- value_after(
  "--layout-tsv",
  "/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.og.lay.tsv"
)
out_png <- value_after("--out", "pggb_graph_2d.png")
component_id <- as.integer(value_after("--component", "8"))
render_log <- value_after("--render-log", "render_log.tsv")

if (!file.exists(layout_tsv)) {
  stop("Layout TSV does not exist: ", layout_tsv)
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

png(out_png, width = 1920, height = 1080, res = 144, bg = "white", type = "cairo")
par(mar = c(0.18, 0.18, 0.18, 0.18), xaxs = "i", yaxs = "i", family = "sans")
plot(
  plot_x,
  plot_y,
  pch = 16,
  cex = 0.11,
  col = grDevices::adjustcolor("#163b69", alpha.f = 0.16),
  axes = FALSE,
  ann = FALSE,
  asp = 1
)
usr <- par("usr")
rect(usr[1], usr[3], usr[2], usr[4], border = "#dce8f7", lwd = 2)
legend(
  "bottomleft",
  legend = c(
    sprintf("PGGB ODGI layout main component (component %s)", component_id),
    sprintf("%s layout nodes; rotated 90 degrees for 16:9 slide fit", format(nrow(component), big.mark = ","))
  ),
  bty = "n",
  text.col = "#222222",
  cex = 0.82
)
invisible(dev.off())

png_info <- file.info(out_png)

render_record <- data.frame(
  key = c(
    "status",
    "layout_tsv",
    "source_odgi_graph",
    "source_odgi_layout_binary",
    "inspected_existing_odgi_draw_png",
    "inspected_existing_odgi_viz_png",
    "component",
    "component_nodes",
    "component_rank_by_layout_node_count",
    "x_range",
    "y_range",
    "orientation",
    "output_png",
    "output_dimensions",
    "output_bytes",
    "slurm_job_id",
    "render_command"
  ),
  value = c(
    status,
    normalizePath(layout_tsv, mustWork = TRUE),
    "/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.og",
    "/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.og.lay",
    "/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.og.lay.draw.png (167x1000; too sparse/narrow for slide use)",
    "/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.og.viz_multiqc.png (1968x157164; too tall/raw for slide use)",
    as.character(component_id),
    as.character(nrow(component)),
    as.character(component_rank),
    paste(range(component$X), collapse = "-"),
    paste(range(component$Y), collapse = "-"),
    "ODGI X/Y coordinates plotted as Y/X to rotate the tall layout into a 16:9 slide frame; no re-layout",
    display_path(out_png),
    "1920x1080",
    as.character(png_info$size),
    "not_used_no_heavy_ODGI_extraction_or_render_submitted",
    paste(c("Rscript", display_path(script_path), args), collapse = " ")
  ),
  stringsAsFactors = FALSE
)

write.table(render_record, render_log, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Wrote ", out_png, "\n", sep = "")
cat("Wrote ", render_log, "\n", sep = "")
