options(scipen = 10000)

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(ggnewscale)
library(viridisLite)

# Resolve repo root from this script's own location (submission/scripts/figures/)
# so it runs from anywhere and reads the vendored data/ inputs.
.args <- commandArgs(trailingOnly = FALSE)
.f    <- sub("^--file=", "", .args[grep("^--file=", .args)])
root  <- if (length(.f)) normalizePath(file.path(dirname(.f), "..", "..", "..")) else getwd()
out_dir <- file.path(root, "submission/fig/MainFigures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Drop only trailing windows with no inter-chromosomal signal (the empty telomeric
# tip at the distal end); keep leading windows for context.
trim_empty_windows <- function(df) {
  pk <- df$position[df$haplo_count > 0]
  if (!length(pk)) return(df)
  df[df$position <= max(pk), , drop = FALSE]
}

# Parameters for images
width <- 16
height <- 10
dpi <- 300

# Read the data
prefix <- 'hprc25272' # hprc25272 or hprc7524
suffix <- 'wf' # wf or fg
window_size <- '100kb'
l_size <- '50000'
num_haplo <- 466
num_sample <- 234
data <- read_tsv(file.path(root, "data/hprc25272-wf.CHM13.100kb-xm5-id098-l50000.tsv.gz"))

# Parse the chroms-num_haplotypes column to extract chromosome information
parse_chroms_column <- function(chroms_str) {
  if (is.na(chroms_str) || chroms_str == "") {
    return(list(num_chromosomes = 0, chromosomes = NA))
  }

  # Split by comma
  pairs <- strsplit(chroms_str, ",")[[1]]

  # Extract chromosome names (everything before the last hyphen)
  chrom_names <- sapply(pairs, function(x) {
    parts <- strsplit(x, "-")[[1]]
    # Join all but the last part (in case chromosome name contains hyphen)
    paste(parts[-length(parts)], collapse = "-")
  })

  # Remove any whitespace
  chrom_names <- trimws(chrom_names)

  return(list(
    num_chromosomes = length(unique(chrom_names)),
    chromosomes = paste(unique(chrom_names), collapse = ",")
  ))
}

# Apply parsing to the entire dataset
parsed_data <- lapply(data$`chroms-num_haplotypes`, parse_chroms_column)
data$num_chromosomes <- sapply(parsed_data, function(x) x$num_chromosomes)
data$chromosomes <- sapply(parsed_data, function(x) x$chromosomes)
rm(parsed_data)

# Extract chromosome information from the chrom column
data <- data %>%
  mutate(
    chromosome = gsub("CHM13#0#(chr[^\\s]+).*", "\\1", chrom),
    # Extract just the number/letter from the chromosome
    chrom_num = gsub("chr([0-9XYM]+)", "\\1", chromosome)
  ) %>%
  # Create numeric position for sorting (X, Y, M at the end)
  left_join(
    tibble(
      chrom_num = c(as.character(1:22), "X", "Y", "M"),
      chrom_order = c(1:22, 23, 24, 25)
    ),
    by = "chrom_num"
  )

# Convert start and end positions to Mbp, also keep midpoint for heatmaps
data <- data %>%
  mutate(
    start_mbp = start / 1e6,
    end_mbp = end / 1e6,
    position = (start + end) / 2 / 1e6  # Keep for heatmaps
  )

# Reorder the factor levels to ensure chromosomes are in the correct order
data$chromosome <- factor(data$chromosome, levels = c(paste0("chr", 1:22), "chrX", "chrY"))
data <- data %>% filter(chromosome != "chrM")

# Optional: Set BED file path here (set to NULL if no BED file)
bed_file_path <- file.path(root, 'data/chm13-annotations.bed')  # Change this to your BED file path

# Function to read and process BED file
read_bed_regions <- function(bed_path) {
  if (!is.null(bed_path) && file.exists(bed_path)) {
    bed_data <- read_tsv(bed_path,
                         col_names = c("chromosome", "start", "end", "name", "score", "strand"),
                         col_types = cols(.default = "c", start = "d", end = "d"))

    # Convert positions to Mbp to match the plot scale
    bed_data <- bed_data %>%
      mutate(
        start_mbp = start / 1e6,
        end_mbp = end / 1e6,
        # Ensure chromosome names match the format in the main data
        chromosome = ifelse(grepl("^chr", chromosome), chromosome, paste0("chr", chromosome)),
        # Create a factor for the name column for consistent coloring
        name = as.factor(name)
      )

    return(bed_data)
  } else {
    return(NULL)
  }
}

# Read the BED file
bed_regions <- read_bed_regions(bed_file_path) %>%
  filter(name != 'PHR-sex') %>%
  mutate(name = if_else(name == 'PHR-acro', 'PHR', name)) %>%
  mutate(name = if_else(name %in% c('PAR1', 'PAR2'), 'PAR', name)) %>%
  mutate(name = if_else(name %in% c('XTR1', 'XTR2'), 'XTR', name)) %>%
  mutate(name = if_else(name == 'Centromere', 'CEN', name))

# Create color palette for BED regions if they exist
if (!is.null(bed_regions)) {
  cat("\nBED regions loaded successfully!\n")
  cat("Number of regions:", nrow(bed_regions), "\n")
  cat("Chromosomes covered:", paste(unique(bed_regions$chromosome), collapse = ", "), "\n\n")

  unique_labels <- unique(bed_regions$name)
  n_labels <- length(unique_labels)

  # Explicit colors per annotation type
  bed_colors <- c(
    "CEN"  = "black",
    "PAR"  = "#E69F00",   # strong orange
    "PHR"  = "#009E73",   # strong teal
    "XTR"  = "#0072B2"    # strong blue
  )
  # Keep only labels present in the data
  bed_colors <- bed_colors[names(bed_colors) %in% unique_labels]

  cat("\nBED region color mapping:\n")
  for (i in 1:length(bed_colors)) {
    cat(names(bed_colors)[i], ":", bed_colors[i], "\n")
  }

  # If BED regions exist, ensure they have the same factor levels
  bed_regions$chromosome <- factor(bed_regions$chromosome, levels = levels(data$chromosome))
} else {
  cat("\nNo BED file loaded. To add region annotations, set bed_file_path to your BED file.\n\n")
}

#===============================================================================
# Inter-chromosomal mapping karyogram
#======================================
# Parse chroms-num_haplotypes to get inter-chromosomal (non-self) mapping info
# For each window, identify which OTHER chromosomes have alignments there

parse_interchrom <- function(data) {
  all_chroms <- paste0("chr", c(1:22, "X", "Y", "M"))

  # Initialize matrix: rows = windows, cols = chromosomes
  count_matrix <- matrix(0L,
                         nrow = nrow(data),
                         ncol = length(all_chroms),
                         dimnames = list(NULL, all_chroms))

  for (i in 1:nrow(data)) {
    raw <- data$`chroms-num_haplotypes`[i]
    if (is.na(raw) || raw == "") next

    pairs <- strsplit(raw, ",")[[1]]
    for (pair in pairs) {
      parts <- strsplit(trimws(pair), "-")[[1]]
      if (length(parts) >= 2) {
        chr_name <- paste(parts[-length(parts)], collapse = "-")
        count <- as.numeric(parts[length(parts)])
        if (chr_name %in% all_chroms && !is.na(count)) {
          count_matrix[i, chr_name] <- count
        }
      }
    }
  }

  return(count_matrix)
}

cat("Parsing inter-chromosomal mapping data...\n")
interchrom_counts <- parse_interchrom(data)

# For each window, zero out the self-chromosome column and compute:
#   - num_other_chroms: how many OTHER chromosomes have alignments
#   - which_other_chroms: comma-separated list of those chromosomes
#   - total_other_haplotypes: sum of haplotypes from other chromosomes
all_chroms <- paste0("chr", c(1:22, "X", "Y", "M"))

data$num_other_chroms <- 0L
data$total_other_haplotypes <- 0L
data$which_other_chroms <- NA_character_

for (i in 1:nrow(data)) {
  self_chr <- data$chromosome[i]
  counts_i <- interchrom_counts[i, ]
  # Zero out self
  if (!is.na(self_chr) && self_chr %in% all_chroms) {
    counts_i[self_chr] <- 0
  }
  nonzero <- counts_i[counts_i > 0]
  data$num_other_chroms[i] <- length(nonzero)
  data$total_other_haplotypes[i] <- sum(nonzero)
  if (length(nonzero) > 0) {
    data$which_other_chroms[i] <- paste(names(nonzero), collapse = ",")
  }
}

cat("Inter-chromosomal parsing complete.\n")
cat("  Windows with inter-chrom mappings:", sum(data$num_other_chroms > 0), "/", nrow(data), "\n")
cat("  Max other chromosomes in a single window:", max(data$num_other_chroms), "\n\n")

# Prepare karyogram data (used by multiple plots below)
karyogram_chrom_levels <- rev(c(paste0("chr", 1:22), "chrX", "chrY"))

karyogram_data <- data %>%
  filter(!is.na(chromosome)) %>%
  mutate(chromosome = factor(chromosome, levels = karyogram_chrom_levels))

max_other <- max(karyogram_data$num_other_chroms, na.rm = TRUE)

# Prepare BED annotations for karyogram-style plots (chromosomes on Y axis)
# Map each annotation region to the numeric Y position of its chromosome
bed_karyogram <- NULL
if (!is.null(bed_regions)) {
  bed_karyogram <- bed_regions %>%
    filter(!is.na(chromosome)) %>%
    mutate(chromosome = factor(chromosome, levels = karyogram_chrom_levels)) %>%
    filter(!is.na(chromosome)) %>%
    mutate(chrom_y = as.numeric(chromosome),
           # pad tiny regions (PAR/PHR are <1 Mb) to a minimum drawn width
           pad = pmax(0, (2.0 - (end_mbp - start_mbp)) / 2),
           xs  = pmax(0, start_mbp - pad), xe = end_mbp + pad)
}

# --- G-band ideogram (CHM13 cytobands) --------------------------------------
gie_cols <- c(gneg = "#f4f4f4", gpos25 = "#cfcfcf", gpos50 = "#9e9e9e",
              gpos75 = "#6b6b6b", gpos100 = "#383838", acen = "#c9636b",
              gvar = "#b3b3b3", stalk = "#8fa0c4")
bands <- read_tsv(file.path(root, "data/chm13v2.0_cytobands_allchrs.bed"),
                  col_names = c("chromosome", "start", "end", "band", "gieStain"),
                  col_types = cols(.default = "c", start = "d", end = "d")) %>%
  filter(chromosome %in% c(paste0("chr", 1:22), "chrX", "chrY")) %>%
  mutate(chromosome = factor(chromosome, levels = karyogram_chrom_levels),
         chrom_y   = as.numeric(chromosome),
         start_mbp = start / 1e6, end_mbp = end / 1e6,
         gieStain  = factor(gieStain, levels = names(gie_cols)))



# --- Plot 1: Karyogram colored by number of inter-chromosomal mappings ---
# Each CHM13 chromosome as a compact horizontal bar, colored by how many OTHER
# chromosomes have alignments mapping to that window.

p_karyogram_count <- ggplot(karyogram_data)

if (!is.null(bed_karyogram)) {
  p_karyogram_count <- p_karyogram_count +
    geom_rect(data = bed_karyogram,
              aes(xmin = start_mbp, xmax = end_mbp,
                  ymin = chrom_y - 0.44, ymax = chrom_y + 0.12,
                  fill = name),
              inherit.aes = FALSE, alpha = 0.6) +
    scale_fill_manual(values = bed_colors, name = "Region",
                      guide = guide_legend(override.aes = list(alpha = 1),
                                           keywidth = unit(0.5, "cm"),
                                           keyheight = unit(0.5, "cm"))) +
    new_scale_fill()
}

p_karyogram_count <- p_karyogram_count +
  geom_rect(aes(xmin = start_mbp, xmax = end_mbp,
                ymin = as.numeric(chromosome) - 0.4,
                ymax = as.numeric(chromosome) + 0.4,
                fill = num_other_chroms),
            linewidth = 0) +
  scale_fill_gradientn(
    colors = c("grey95", "#fee08b", "#fc8d59", "#d73027", "#67001f"),
    values = scales::rescale(c(0, 1, 3, 8, max(max_other, 10))),
    limits = c(0, max_other),
    breaks = scales::pretty_breaks(n = 5),
    name = "# other\nchromosomes"
  ) +
  scale_x_continuous(breaks = seq(0, 300, 50), expand = c(0.01, 0)) +
  scale_y_continuous(
    breaks = 1:length(levels(karyogram_data$chromosome)),
    labels = levels(karyogram_data$chromosome),
    expand = expansion(add = c(0.3, 0.7))
  ) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Position (Mbp)",
    y = NULL,
    caption = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 13),
    plot.caption = element_text(size = 9, hjust = 1, color = "gray40"),
    plot.margin = margin(26, 2, 2, 2, "pt"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.15),
    legend.position = "right",
    legend.key.height = unit(1.2, "cm"),
    legend.key.width = unit(0.6, "cm"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.spacing.y = unit(0.3, "cm")
  )

# --- Plot 1b: Same karyogram with chr4q-end rainbow inset ---
library(patchwork)

# Last 400 kbp of chr4q (last 4 windows — where the inter-chr signal is)
chr4_end <- max(data$end_mbp[data$chromosome == "chr10"], na.rm = TRUE)
chr4_inset_start <- chr4_end - 0.4  # 0.4 Mbp = 4 windows of 100kb

# Highlight box on chr4 in the main plot
chr4_y <- as.numeric(factor("chr10", levels = karyogram_chrom_levels))
p_karyogram_count_main <- p_karyogram_count +
  annotate("rect", xmin = chr4_inset_start, xmax = chr4_end,
           ymin = chr4_y - 0.5, ymax = chr4_y + 0.5,
           color = "black", fill = NA, linewidth = 0.6, linetype = "solid")

# Build heatmap data for inset: x = window, y = source chromosome, fill = haplotype count
# This shows both WHICH chromosomes and HOW MANY haplotypes support each inter-chr mapping
inset_rows <- which(data$chromosome == "chr10" &
                    data$start_mbp >= chr4_inset_start &
                    data$end_mbp <= chr4_end)

inset_heatmap <- cbind(
  data[inset_rows, ] %>% select(start_mbp, end_mbp),
  interchrom_counts[inset_rows, ]
) %>%
  as.data.frame() %>%
  pivot_longer(cols = all_of(all_chroms),
               names_to = "source_chr",
               values_to = "haplo_count") %>%
  # Remove self (chr4)
  filter(source_chr != "chr10") %>%
  mutate(
    source_chr = factor(source_chr, levels = rev(all_chroms[all_chroms != "chr10"])),
    position = (start_mbp + end_mbp) / 2
  )

# Only keep chromosomes that have at least one non-zero entry
chroms_with_signal <- inset_heatmap %>%
  group_by(source_chr) %>%
  summarize(total = sum(haplo_count), .groups = "drop") %>%
  filter(total > 0) %>%
  pull(source_chr)

# Re-level factor to only include chromosomes with signal (in natural order)
chroms_with_signal_ordered <- rev(all_chroms[all_chroms != "chr10" & all_chroms %in% chroms_with_signal])

inset_heatmap_filtered <- inset_heatmap %>%
  filter(source_chr %in% chroms_with_signal) %>%
  mutate(source_chr = factor(source_chr, levels = chroms_with_signal_ordered))

# Tile width for proper rendering
tile_w <- if (nrow(inset_heatmap_filtered) > 0) {
  median(diff(sort(unique(inset_heatmap_filtered$position))), na.rm = TRUE)
} else { 0.1 }

# Add text color column (white on dark tiles, black on light)
inset_heatmap_filtered <- inset_heatmap_filtered %>%
  mutate(text_color = ifelse(haplo_count > max(haplo_count, na.rm = TRUE) * 0.45, "white", "grey30"))
inset_heatmap_filtered <- trim_empty_windows(inset_heatmap_filtered)
i1_lo <- min(inset_heatmap_filtered$start_mbp); i1_hi <- max(inset_heatmap_filtered$end_mbp)

p_inset <- ggplot(inset_heatmap_filtered,
                  aes(x = position, y = source_chr, fill = haplo_count)) +
  geom_tile(width = tile_w * 1.05, height = 1) +
  geom_text(data = inset_heatmap_filtered %>% filter(haplo_count > 0),
            aes(label = haplo_count, color = text_color),
            size = 2.5, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_gradientn(
    colors = c("white", "#deebf7", "#9ecae1", "#4292c6", "#08306b"),
    limits = c(0, max(inset_heatmap_filtered$haplo_count, na.rm = TRUE)),
    breaks = scales::pretty_breaks(n = 3),
    guide = "none"
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4), expand = c(0.01, 0)) +
  coord_cartesian(xlim = c(i1_lo, i1_hi)) +
  labs(x = "Mbp", y = NULL,
       title = paste0("chr10 end  ", round(i1_lo, 1),
                      "\u2013", round(i1_hi, 1), " Mbp")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    axis.title.x = element_text(size = 9),
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(3, 5, 3, 5, "pt"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.position = "left",
    legend.justification = c(0, 0.5),  # vertically centered
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.3, "cm"),
    legend.background = element_rect(fill = "transparent", color = NA)
  )

# --- Second inset: chr5 q-end (last 400 kbp) -- strongest community + pedigree
# case: shares with chr6 (C11 OR4F 5q/6q community) and chr1 (pedigree chr5q/chr1p).
chr5_end_mb <- max(data$end_mbp[data$chromosome == "chr5"], na.rm = TRUE)
chr1_inset_start <- chr5_end_mb - 0.4
chr1_inset_end <- chr5_end_mb

# Add highlight box on chr5 (q-end)
chr1_y <- as.numeric(factor("chr5", levels = karyogram_chrom_levels))
p_karyogram_count_main <- p_karyogram_count_main +
  annotate("rect", xmin = chr1_inset_start, xmax = chr1_inset_end,
           ymin = chr1_y - 0.5, ymax = chr1_y + 0.5,
           color = "black", fill = NA, linewidth = 0.6, linetype = "solid")

# Build heatmap for the chr5 q-end window
inset2_rows <- which(data$chromosome == "chr5" &
                     data$start_mbp >= chr1_inset_start &
                     data$end_mbp <= chr1_inset_end)

inset2_heatmap <- cbind(
  data[inset2_rows, ] %>% select(start_mbp, end_mbp),
  interchrom_counts[inset2_rows, ]
) %>%
  as.data.frame() %>%
  pivot_longer(cols = all_of(all_chroms),
               names_to = "source_chr",
               values_to = "haplo_count") %>%
  filter(source_chr != "chr5") %>%
  mutate(
    source_chr = factor(source_chr, levels = rev(all_chroms[all_chroms != "chr5"])),
    position = (start_mbp + end_mbp) / 2
  )

# Only keep chromosomes with signal
chroms_with_signal2 <- inset2_heatmap %>%
  group_by(source_chr) %>%
  summarize(total = sum(haplo_count), .groups = "drop") %>%
  filter(total > 0) %>%
  pull(source_chr)

chroms_with_signal2_ordered <- rev(all_chroms[all_chroms != "chr5" & all_chroms %in% chroms_with_signal2])

inset2_filtered <- inset2_heatmap %>%
  filter(source_chr %in% chroms_with_signal2) %>%
  mutate(source_chr = factor(source_chr, levels = chroms_with_signal2_ordered))

tile_w2 <- if (nrow(inset2_filtered) > 0) {
  median(diff(sort(unique(inset2_filtered$position))), na.rm = TRUE)
} else { 0.1 }

inset2_filtered <- inset2_filtered %>%
  mutate(text_color = ifelse(haplo_count > max(haplo_count, na.rm = TRUE) * 0.45, "white", "grey30"))
inset2_filtered <- trim_empty_windows(inset2_filtered)
i2_lo <- min(inset2_filtered$start_mbp); i2_hi <- max(inset2_filtered$end_mbp)

p_inset2 <- ggplot(inset2_filtered,
                   aes(x = position, y = source_chr, fill = haplo_count)) +
  geom_tile(width = tile_w2 * 1.05, height = 1) +
  geom_text(data = inset2_filtered %>% filter(haplo_count > 0),
            aes(label = haplo_count, color = text_color),
            size = 2.5, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_gradientn(
    colors = c("white", "#deebf7", "#9ecae1", "#4292c6", "#08306b"),
    limits = c(0, max(inset2_filtered$haplo_count, na.rm = TRUE)),
    breaks = scales::pretty_breaks(n = 3),
    guide = "none"
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4), expand = c(0.01, 0)) +
  coord_cartesian(xlim = c(i2_lo, i2_hi)) +
  labs(x = "Mbp", y = NULL,
       title = paste0("chr5 end  ", round(i2_lo, 1), "\u2013", round(i2_hi, 1), " Mbp")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    axis.title.x = element_text(size = 9),
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(3, 5, 3, 5, "pt"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.position = "left",
    legend.justification = c(0, 0.5),  # vertically centered
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.3, "cm"),
    legend.background = element_rect(fill = "transparent", color = NA)
  )

# --- Third inset: chr13 acrocentric p arm (0-5 Mbp) ---
acro_chr <- "chr21"
acro_inset_start <- 0
acro_inset_end <- 5

# Add highlight box on acrocentric chromosome
acro_y <- as.numeric(factor(acro_chr, levels = karyogram_chrom_levels))
p_karyogram_count_main <- p_karyogram_count_main +
  annotate("rect", xmin = acro_inset_start, xmax = acro_inset_end,
           ymin = acro_y - 0.5, ymax = acro_y + 0.5,
           color = "black", fill = NA, linewidth = 0.6, linetype = "solid")

# Build heatmap for acrocentric p arm
inset3_rows <- which(data$chromosome == acro_chr &
                     data$start_mbp >= acro_inset_start &
                     data$end_mbp <= acro_inset_end)

inset3_heatmap <- cbind(
  data[inset3_rows, ] %>% select(start_mbp, end_mbp),
  interchrom_counts[inset3_rows, ]
) %>%
  as.data.frame() %>%
  pivot_longer(cols = all_of(all_chroms),
               names_to = "source_chr",
               values_to = "haplo_count") %>%
  filter(source_chr != acro_chr) %>%
  mutate(
    source_chr = factor(source_chr, levels = rev(all_chroms[all_chroms != acro_chr])),
    position = (start_mbp + end_mbp) / 2
  )

chroms_with_signal3 <- inset3_heatmap %>%
  group_by(source_chr) %>%
  summarize(total = sum(haplo_count), .groups = "drop") %>%
  filter(total > 0) %>%
  pull(source_chr)

chroms_with_signal3_ordered <- rev(all_chroms[all_chroms != acro_chr & all_chroms %in% chroms_with_signal3])

inset3_filtered <- inset3_heatmap %>%
  filter(source_chr %in% chroms_with_signal3) %>%
  mutate(source_chr = factor(source_chr, levels = chroms_with_signal3_ordered))
inset3_filtered <- trim_empty_windows(inset3_filtered)
i3_lo <- min(inset3_filtered$start_mbp); i3_hi <- max(inset3_filtered$end_mbp)

tile_w3 <- if (nrow(inset3_filtered) > 0) {
  median(diff(sort(unique(inset3_filtered$position))), na.rm = TRUE)
} else { 0.1 }

p_inset3 <- ggplot(inset3_filtered,
                   aes(x = position, y = source_chr, fill = haplo_count)) +
  geom_tile(width = tile_w3 * 1.05, height = 1) +
  scale_fill_gradientn(
    colors = c("white", "#deebf7", "#9ecae1", "#4292c6", "#08306b"),
    limits = c(0, max(inset3_filtered$haplo_count, na.rm = TRUE)),
    breaks = scales::pretty_breaks(n = 3),
    guide = "none"
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4), expand = c(0.01, 0)) +
  coord_cartesian(xlim = c(i3_lo, i3_hi)) +
  labs(x = "Mbp", y = NULL,
       title = paste0(acro_chr, " p arm  ", round(i3_lo, 1), "\u2013", round(i3_hi, 1), " Mbp")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    axis.title.x = element_text(size = 9),
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(3, 5, 3, 5, "pt"),
    legend.position = "none"
  )


# Place all three insets
# chr1 top-right, chr4q bottom-right, chr13 middle-left of the right column
n_chr_inset1 <- length(chroms_with_signal_ordered)   # chr4q
n_chr_inset2 <- length(chroms_with_signal2_ordered)  # chr1
total_chr <- n_chr_inset1 + n_chr_inset2
available_height <- 0.89
gap <- 0.02
chr4_extra <- 0.03
inset1_height <- available_height * n_chr_inset1 / total_chr - gap / 2 + chr4_extra
inset2_height <- available_height * n_chr_inset2 / total_chr - gap / 2 - chr4_extra

cat("Inset chr4q: ", n_chr_inset1, "chroms, height =", round(inset1_height, 3), "\n")
cat("Inset chr1:  ", n_chr_inset2, "chroms, height =", round(inset2_height, 3), "\n")

# Right column (chr1 top, chr4q bottom): narrower than before, moved right
right_left <- 0.78
right_right <- 0.97

# Vertical base offset for all three insets (lower = further down)
inset_base_y <- -0.01
inset2_yshift <- 0.0     # chr5-end inset: no shift (its earlier position)
inset3_yshift <- 0.10    # push the acrocentric (chr14 p-arm) inset lower

# Inset3 (chr13): placed to the LEFT of the right column, centered vertically
# between inset1 (top) and inset2 (bottom)
n_chr_inset3 <- length(chroms_with_signal3_ordered)
inset3_height <- 0.46
inset3_top <- inset_base_y + inset1_height + gap / 2 + inset3_height / 2 - inset3_yshift
inset3_bottom <- inset3_top - inset3_height

p_karyogram_count_with_inset <- p_karyogram_count_main +
  inset_element(p_inset2,
                left = right_left, right = right_right,
                bottom = inset_base_y + inset1_height + gap - inset2_yshift,
                top = inset_base_y + inset1_height + gap + inset2_height - inset2_yshift) +
  inset_element(p_inset,
                left = right_left, right = right_right,
                bottom = inset_base_y,
                top = inset_base_y + inset1_height) +
  inset_element(p_inset3,
                left = 0.57, right = 0.755,
                bottom = inset3_bottom, top = inset3_top)

# Sequential-gradient variant (not used for the manuscript figure) -> scratch.
ggsave(
  filename = file.path(tempdir(), "p_interchrom_karyogram_count_inset.png"),
  plot = p_karyogram_count_with_inset,
  width = 16, height = 8.5, dpi = dpi, units = "in", bg = "white"
)


# --- Plot 2: Karyogram with rainbow-colored count ---
# Same compact bar layout as Plot 1, but the number of other chromosomes is
# shown as a discrete rainbow palette (spectral) instead of a sequential gradient.

# Bin num_other_chroms into discrete categories for rainbow coloring
count_breaks <- c(0, 1, 2, 3, 4, 5, 10, 15, Inf)
count_labels <- c("0", "1", "2", "3", "4", "5", "6-10", "11+")
# Adjust breaks/labels to match actual max, collapsing bins above max_other
if (max_other <= 5) {
  count_breaks <- c(0:max_other, max_other + 1)
  count_labels <- as.character(0:max_other)
} else if (max_other < 11) {
  count_breaks <- c(0, 1, 2, 3, 4, 5, 10, max_other + 1)
  count_labels <- c("0", "1", "2", "3", "4", "5", paste0("6-", max_other))
}

karyogram_data$count_bin <- cut(
  karyogram_data$num_other_chroms,
  breaks = count_breaks,
  labels = count_labels,
  include.lowest = TRUE, right = FALSE
)

# Turbo palette: spans full visible spectrum (blue→cyan→green→yellow→orange→red)
# Each bin lands in a distinct hue family, maximising distinguishability
n_bins <- length(count_labels)
# Viridis instead of turbo, to differ from the companion coverage paper, with the
# low-count steps pulled apart (plain viridis(7) compresses 1/2/3 in the dark end).
nb   <- n_bins - 1
vpos <- if (nb == 7) c(0.00, 0.34, 0.50, 0.62, 0.74, 0.87, 1.00) else seq(0, 1, length.out = nb)
count_rainbow <- setNames(
  c("grey95", viridis(101, option = "C")[round(vpos * 100) + 1]),
  count_labels
)

p_karyogram_count_rainbow <- ggplot(karyogram_data) +
  # G-band ideogram as a track above each chromosome bar
  geom_rect(data = bands,
            aes(xmin = start_mbp, xmax = end_mbp,
                ymin = chrom_y + 0.13, ymax = chrom_y + 0.40, fill = gieStain),
            inherit.aes = FALSE) +
  scale_fill_manual(values = gie_cols, guide = "none") +
  new_scale_fill()

if (!is.null(bed_karyogram)) {
  p_karyogram_count_rainbow <- p_karyogram_count_rainbow +
    geom_rect(data = bed_karyogram,
              aes(xmin = xs, xmax = xe,
                  ymin = chrom_y - 0.40, ymax = chrom_y - 0.28,
                  fill = name),
              inherit.aes = FALSE) +
    scale_fill_manual(values = bed_colors, name = "Region",
                      guide = guide_legend(keywidth = unit(0.5, "cm"),
                                           keyheight = unit(0.5, "cm"))) +
    new_scale_fill()
}

p_karyogram_count_rainbow <- p_karyogram_count_rainbow +
  geom_rect(aes(xmin = start_mbp, xmax = end_mbp,
                ymin = as.numeric(chromosome) - 0.25,
                ymax = as.numeric(chromosome) + 0.09,
                fill = count_bin),
            linewidth = 0) +
  scale_fill_manual(values = count_rainbow, name = "# other\nchromosomes",
                    guide = guide_legend(keywidth = unit(0.6, "cm"),
                                         keyheight = unit(0.6, "cm"))) +
  scale_x_continuous(breaks = seq(0, 300, 50), expand = c(0.01, 0)) +
  scale_y_continuous(
    breaks = 1:length(levels(karyogram_data$chromosome)),
    labels = levels(karyogram_data$chromosome),
    expand = expansion(add = c(0.3, 0.7))
  ) +
  labs(
    title = NULL, subtitle = NULL,
    x = "Position (Mbp)",
    y = NULL,
    caption = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 13),
    plot.caption = element_text(size = 9, hjust = 1, color = "gray40"),
    plot.margin = margin(26, 2, 2, 2, "pt"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.15),
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.spacing.y = unit(0.3, "cm")
  )

# --- Rainbow version with the same three insets ---
# Build rainbow variants of the insets (same data, different color scale)

# chr4q end — rainbow
inset_heatmap_filtered_r <- inset_heatmap_filtered %>%
  mutate(count_bin = cut(haplo_count, breaks = c(-Inf, 0, 1, 5, 20, 100, 300, Inf),
                         labels = c("0", "1", "2-5", "6-20", "21-100", "101-300", "300+"),
                         include.lowest = TRUE))

# Reuse exact same p_inset code but for rainbow, we just keep the same blue gradient
# since the user said the only difference is the main plot's # other chromosomes coloring.
# Insets stay as blue heatmaps in both plot variants.

# Apply highlight boxes to the rainbow main plot (same as count_main)
p_karyogram_count_rainbow_main <- p_karyogram_count_rainbow

p_karyogram_count_rainbow_with_inset <- p_karyogram_count_rainbow_main +
  inset_element(p_inset2,
                left = right_left, right = right_right,
                bottom = inset_base_y + inset1_height + gap - inset2_yshift,
                top = inset_base_y + inset1_height + gap + inset2_height - inset2_yshift) +
  inset_element(p_inset,
                left = right_left, right = right_right,
                bottom = inset_base_y,
                top = inset_base_y + inset1_height) +
  inset_element(p_inset3,
                left = 0.57, right = 0.755,
                bottom = inset3_bottom, top = inset3_top)

# Manuscript Figure 1a (viridis, chr1/chr14/chr10 insets).
ggsave(
  filename = file.path(out_dir, "Fig1a_genomewide.png"),
  plot = p_karyogram_count_rainbow_with_inset,
  width = 16, height = 8.5, dpi = dpi, units = "in", bg = "white"
)
cat("Saved ", file.path(out_dir, "Fig1a_genomewide.png"), "\n")
ggsave(
  filename = file.path(out_dir, "Fig1a_genomewide.pdf"),
  plot = p_karyogram_count_rainbow_with_inset,
  width = 16, height = 8.5, units = "in", bg = "white",
  device = cairo_pdf
)
cat("Saved ", file.path(out_dir, "Fig1a_genomewide.pdf"), "\n")
