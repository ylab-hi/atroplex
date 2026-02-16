#!/usr/bin/env Rscript
#' Manuscript-ready visualization of SQANTI3 classification results.
#'
#' Auto-discovers samples from a SQANTI3 output directory containing
#' subfolders with *_classification.txt files. Each figure is saved as
#' a separate PDF+PNG. ISM-specific figures go into an ism/ subfolder.
#'
#' Usage:
#'   Rscript visualize_sqanti3.R <sqanti3_dir> [-o output_dir]

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(scales)
  library(patchwork)
  library(pheatmap)
})

# ── Configuration ────────────────────────────────────────────────────────────

CATEGORY_ORDER <- c(
  "full-splice_match", "incomplete-splice_match",
  "novel_in_catalog", "novel_not_in_catalog",
  "genic", "genic_intron", "antisense", "intergenic", "fusion"
)

CATEGORY_SHORT <- c(
  "full-splice_match"       = "FSM",
  "incomplete-splice_match" = "ISM",
  "novel_in_catalog"        = "NIC",
  "novel_not_in_catalog"    = "NNC",
  "genic"                   = "Genic",
  "genic_intron"            = "Genic intron",
  "antisense"               = "Antisense",
  "intergenic"              = "Intergenic",
  "fusion"                  = "Fusion"
)

CATEGORY_COLORS <- c(
  "FSM"          = "#2166AC",
  "ISM"          = "#67A9CF",
  "NIC"          = "#D6604D",
  "NNC"          = "#F4A582",
  "Genic"        = "#B2ABD2",
  "Genic intron" = "#878787",
  "Antisense"    = "#FDB863",
  "Intergenic"   = "#E0E0E0",
  "Fusion"       = "#1B7837"
)

ISM_SUBCAT_COLORS <- c(
  "3prime fragment"  = "#66c2a5",
  "5prime fragment"  = "#fc8d62",
  "internal fragment" = "#8da0cb",
  "mono-exon"        = "#e78ac3",
  "intron retention" = "#a6d854"
)

theme_manuscript <- theme_bw(base_size = 8) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    legend.key.size = unit(0.35, "cm"),
    plot.title = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6)
  )

# ── Helpers ──────────────────────────────────────────────────────────────────

discover_samples <- function(sqanti_dir) {
  dirs <- list.dirs(sqanti_dir, recursive = FALSE, full.names = FALSE)
  samples <- dirs[file.exists(file.path(sqanti_dir, dirs,
                                        paste0(dirs, "_classification.txt")))]
  sort(samples)
}

load_classifications <- function(sqanti_dir, sample_names) {
  map_dfr(sample_names, function(s) {
    path <- file.path(sqanti_dir, s, paste0(s, "_classification.txt"))
    read_tsv(path, col_types = cols(.default = "c"), show_col_types = FALSE) |>
      mutate(
        sample = s,
        exons = as.integer(exons),
        ref_exons = as.integer(ref_exons),
        length = as.integer(length),
        ref_length = as.integer(ref_length),
        diff_to_TSS = as.numeric(diff_to_TSS),
        diff_to_TTS = as.numeric(diff_to_TTS)
      )
  })
}

save_fig <- function(p, path, w, h) {
  ggsave(paste0(path, ".pdf"), p, width = w, height = h, units = "in", dpi = 300)
  ggsave(paste0(path, ".png"), p, width = w, height = h, units = "in", dpi = 300)
}

# ── Figure: Structural category counts ───────────────────────────────────────

fig_category_counts <- function(data, sample_names, out_dir) {
  cats <- intersect(CATEGORY_ORDER, unique(data$structural_category))
  counts <- data |>
    filter(structural_category %in% cats) |>
    mutate(
      sample = factor(sample, levels = sample_names),
      cat_short = factor(CATEGORY_SHORT[structural_category],
                         levels = CATEGORY_SHORT[cats])
    ) |>
    count(sample, cat_short)

  p <- ggplot(counts, aes(x = sample, y = n, fill = cat_short)) +
    geom_col(width = 0.8, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = CATEGORY_COLORS, name = NULL) +
    scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k")) +
    labs(title = "SQANTI3 structural categories per sample",
         y = "Number of transcripts", x = NULL) +
    theme_manuscript

  save_fig(p, file.path(out_dir, "category_counts"),
           w = max(7, length(sample_names) * 0.45), h = 3.5)
  message("  category_counts")
}

# ── Figure: Structural category proportions ──────────────────────────────────

fig_category_proportions <- function(data, sample_names, out_dir) {
  cats <- intersect(CATEGORY_ORDER, unique(data$structural_category))
  counts <- data |>
    filter(structural_category %in% cats) |>
    mutate(
      sample = factor(sample, levels = sample_names),
      cat_short = factor(CATEGORY_SHORT[structural_category],
                         levels = CATEGORY_SHORT[cats])
    ) |>
    count(sample, cat_short)

  p <- ggplot(counts, aes(x = sample, y = n, fill = cat_short)) +
    geom_col(position = "fill", width = 0.8, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = CATEGORY_COLORS, name = NULL) +
    scale_y_continuous(labels = percent_format()) +
    labs(title = "SQANTI3 structural category proportions",
         y = "Proportion of transcripts", x = NULL) +
    theme_manuscript

  save_fig(p, file.path(out_dir, "category_proportions"),
           w = max(7, length(sample_names) * 0.45), h = 3.5)
  message("  category_proportions")
}

# ── Figure: Subcategory breakdown (FSM, ISM, NIC, NNC) ──────────────────────

fig_subcategory_breakdown <- function(data, sample_names, out_dir) {
  cats_of_interest <- c(
    "full-splice_match"       = "FSM",
    "incomplete-splice_match" = "ISM",
    "novel_in_catalog"        = "NIC",
    "novel_not_in_catalog"    = "NNC"
  )
  cats_of_interest <- cats_of_interest[names(cats_of_interest) %in%
                                         unique(data$structural_category)]
  if (length(cats_of_interest) == 0) return(invisible())

  plots <- imap(cats_of_interest, function(short, cat) {
    sub <- data |>
      filter(structural_category == cat) |>
      mutate(sample = factor(sample, levels = sample_names)) |>
      count(sample, subcategory)
    sub_order <- sub |> count(subcategory, wt = n) |> arrange(desc(n))
    sub$subcategory <- factor(sub$subcategory, levels = sub_order$subcategory)

    ggplot(sub, aes(x = sample, y = n, fill = subcategory)) +
      geom_col(position = "fill", width = 0.8, colour = "white", linewidth = 0.15) +
      scale_fill_brewer(palette = "Set3", name = NULL) +
      scale_y_continuous(labels = percent_format()) +
      labs(title = paste(short, "subcategories"), y = NULL, x = NULL) +
      theme_manuscript +
      theme(legend.text = element_text(size = 5),
            legend.key.size = unit(0.25, "cm"),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4))
  })

  p <- wrap_plots(plots, nrow = 1) +
    plot_annotation(title = "SQANTI3 subcategory breakdown",
                    theme = theme(plot.title = element_text(face = "bold", size = 11)))
  save_fig(p, file.path(out_dir, "subcategories"),
           w = length(cats_of_interest) * 3.5, h = 4)
  message("  subcategories")
}

# ── Figure: Sample similarity heatmap ────────────────────────────────────────

fig_sample_heatmap <- function(data, sample_names, out_dir) {
  sample_sets <- data |>
    distinct(sample, associated_gene, structural_category) |>
    mutate(key = paste(associated_gene, structural_category, sep = "::")) |>
    split(~sample) |>
    map(~ .x$key)

  n <- length(sample_names)
  jaccard_mat <- matrix(0, n, n, dimnames = list(sample_names, sample_names))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      si <- sample_sets[[sample_names[i]]]
      sj <- sample_sets[[sample_names[j]]]
      inter <- length(intersect(si, sj))
      uni <- length(union(si, sj))
      jaccard_mat[i, j] <- if (uni > 0) inter / uni else 0
    }
  }

  size <- max(5, n * 0.4)
  pdf(file.path(out_dir, "sample_heatmap.pdf"), width = size, height = size)
  pheatmap(jaccard_mat,
           clustering_method = "average",
           color = colorRampPalette(c("white", "#FDAE61", "#D73027"))(100),
           display_numbers = TRUE, number_format = "%.2f",
           fontsize_number = max(4, 7 - n %/% 8), fontsize = 7,
           main = "Sample similarity (gene x category Jaccard)")
  dev.off()

  png(file.path(out_dir, "sample_heatmap.png"),
      width = size, height = size, units = "in", res = 300)
  pheatmap(jaccard_mat,
           clustering_method = "average",
           color = colorRampPalette(c("white", "#FDAE61", "#D73027"))(100),
           display_numbers = TRUE, number_format = "%.2f",
           fontsize_number = max(4, 7 - n %/% 8), fontsize = 7,
           main = "Sample similarity (gene x category Jaccard)")
  dev.off()
  message("  sample_heatmap")
}

# ── Figure: Exon count distributions ─────────────────────────────────────────

fig_exon_distributions <- function(data, sample_names, out_dir) {
  cats <- c("full-splice_match", "incomplete-splice_match",
            "novel_in_catalog", "novel_not_in_catalog")
  cats <- cats[cats %in% unique(data$structural_category)]
  if (length(cats) == 0) return(invisible())

  sub <- data |>
    filter(structural_category %in% cats) |>
    mutate(
      sample = factor(sample, levels = sample_names),
      cat_short = factor(CATEGORY_SHORT[structural_category],
                         levels = CATEGORY_SHORT[cats]),
      exons_capped = pmin(exons, 30)
    )

  p <- ggplot(sub, aes(x = sample, y = exons_capped, fill = sample)) +
    geom_violin(scale = "width", linewidth = 0.2, show.legend = FALSE) +
    stat_summary(fun = median, geom = "point", size = 0.5, colour = "black") +
    facet_wrap(~cat_short, nrow = 1) +
    scale_fill_manual(values = rep(hue_pal()(min(length(sample_names), 20)),
                                   length.out = length(sample_names))) +
    labs(title = "Exon count distributions by structural category",
         y = "Exon count", x = NULL) +
    theme_manuscript +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4))

  save_fig(p, file.path(out_dir, "exon_distributions"),
           w = length(cats) * 3.5, h = 3.5)
  message("  exon_distributions")
}

# ══════════════════════════════════════════════════════════════════════════════
# ISM figures (output to ism/ subfolder)
# ══════════════════════════════════════════════════════════════════════════════

# ── ISM: Subcategory counts per sample ───────────────────────────────────────

fig_ism_subcategories <- function(ism, sample_names, ism_dir) {
  sub_order <- ism |> count(subcategory) |> arrange(desc(n)) |> pull(subcategory)
  p <- ism |>
    mutate(sample = factor(sample, levels = sample_names),
           subcategory = factor(subcategory, levels = sub_order)) |>
    count(sample, subcategory) |>
    ggplot(aes(x = sample, y = n, fill = subcategory)) +
    geom_col(width = 0.8, colour = "white", linewidth = 0.15) +
    scale_fill_brewer(palette = "Set2", name = NULL) +
    scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k")) +
    labs(title = "ISM subcategories per sample", y = "ISM count", x = NULL) +
    theme_manuscript
  save_fig(p, file.path(ism_dir, "ism_subcategories"), w = 7, h = 3.5)
  message("  ism/ism_subcategories")
}

# ── ISM: Exon retention ratio ────────────────────────────────────────────────

fig_ism_exon_retention <- function(ism, ism_dir) {
  ism_multi <- ism |>
    filter(exons > 1, ref_exons > 0) |>
    mutate(exon_ratio = pmin(exons / ref_exons, 1.5))

  p <- ism_multi |>
    filter(subcategory %in% c("3prime_fragment", "5prime_fragment", "internal_fragment")) |>
    mutate(subcategory = gsub("_", " ", subcategory)) |>
    ggplot(aes(x = exon_ratio, colour = subcategory, fill = subcategory)) +
    geom_density(alpha = 0.15, linewidth = 0.6) +
    scale_colour_manual(values = ISM_SUBCAT_COLORS, name = NULL) +
    scale_fill_manual(values = ISM_SUBCAT_COLORS, name = NULL) +
    labs(title = "Exon retention ratio (multi-exon ISM)",
         x = "Observed / reference exon count", y = "Density") +
    theme_manuscript +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  save_fig(p, file.path(ism_dir, "ism_exon_retention"), w = 5, h = 3.5)
  message("  ism/ism_exon_retention")
}

# ── ISM: Missing exons by subcategory ────────────────────────────────────────

fig_ism_missing_exons <- function(ism, ism_dir) {
  ism_with_ref <- ism |>
    filter(ref_exons > 0, exons > 0) |>
    mutate(
      missing_exons = ref_exons - exons,
      subcategory_label = gsub("_", " ", subcategory)
    ) |>
    filter(missing_exons >= 0)

  subcats <- c("3prime_fragment", "5prime_fragment", "mono-exon",
               "internal_fragment", "intron_retention")
  subcats_present <- subcats[subcats %in% unique(ism_with_ref$subcategory)]

  # (A) Missing exons distribution per subcategory
  pa <- ism_with_ref |>
    filter(subcategory %in% subcats_present) |>
    mutate(
      missing_capped = pmin(missing_exons, 20),
      subcategory_label = factor(gsub("_", " ", subcategory),
                                 levels = gsub("_", " ", subcats_present))
    ) |>
    ggplot(aes(x = missing_capped, fill = subcategory_label)) +
    geom_histogram(binwidth = 1, colour = "white", linewidth = 0.15) +
    facet_wrap(~subcategory_label, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = ISM_SUBCAT_COLORS, guide = "none") +
    scale_x_continuous(breaks = seq(0, 20, 5)) +
    scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k")) +
    labs(title = "Number of missing exons per ISM transcript",
         x = "Missing exons (ref_exons - observed exons)", y = "Count") +
    theme_manuscript +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

  # (B) Missing exons vs reference complexity (scatter/hex)
  pb <- ism_with_ref |>
    filter(subcategory %in% c("3prime_fragment", "5prime_fragment")) |>
    mutate(subcategory_label = gsub("_", " ", subcategory)) |>
    ggplot(aes(x = ref_exons, y = missing_exons)) +
    geom_hex(bins = 40) +
    scale_fill_viridis_c(trans = "log10", name = "Count") +
    geom_abline(slope = 0.5, intercept = 0, linetype = "dashed",
                colour = "red", alpha = 0.5) +
    facet_wrap(~subcategory_label) +
    labs(title = "Missing exons vs reference transcript complexity",
         x = "Reference exon count", y = "Missing exons") +
    coord_cartesian(xlim = c(0, 40), ylim = c(0, 30)) +
    theme_manuscript +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

  # (C) Per-sample: median missing exons
  pc <- ism_with_ref |>
    filter(subcategory %in% c("3prime_fragment", "5prime_fragment")) |>
    mutate(subcategory_label = gsub("_", " ", subcategory)) |>
    group_by(sample, subcategory_label) |>
    summarise(
      median_missing = median(missing_exons),
      q25 = quantile(missing_exons, 0.25),
      q75 = quantile(missing_exons, 0.75),
      .groups = "drop"
    ) |>
    ggplot(aes(x = sample, y = median_missing, fill = subcategory_label)) +
    geom_col(position = "dodge", width = 0.7, colour = "white", linewidth = 0.2) +
    geom_errorbar(aes(ymin = q25, ymax = q75),
                  position = position_dodge(width = 0.7), width = 0.3, linewidth = 0.3) +
    scale_fill_manual(values = ISM_SUBCAT_COLORS, name = NULL) +
    labs(title = "Median missing exons per sample (IQR bars)",
         y = "Missing exons", x = NULL) +
    theme_manuscript

  p <- pa / (pb | pc) +
    plot_annotation(
      title = "ISM: missing exon analysis",
      theme = theme(plot.title = element_text(face = "bold", size = 12))
    )
  save_fig(p, file.path(ism_dir, "ism_missing_exons"), w = 12, h = 8)
  message("  ism/ism_missing_exons")
}

# ── ISM: ISM-to-FSM ratio per sample ────────────────────────────────────────

fig_ism_fsm_ratio <- function(ism, fsm, sample_names, ism_dir) {
  ratio_df <- ism |> count(sample, name = "n_ism") |>
    left_join(fsm |> count(sample, name = "n_fsm"), by = "sample") |>
    mutate(ratio = n_ism / n_fsm,
           sample = factor(sample, levels = sample_names),
           above = ratio >= 1)

  p <- ggplot(ratio_df, aes(x = sample, y = ratio, fill = above)) +
    geom_col(width = 0.8, show.legend = FALSE, colour = "white", linewidth = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
    scale_fill_manual(values = c("TRUE" = "#D6604D", "FALSE" = "#67A9CF")) +
    labs(title = "ISM-to-FSM ratio per sample", y = "ISM / FSM ratio", x = NULL) +
    theme_manuscript
  save_fig(p, file.path(ism_dir, "ism_fsm_ratio"), w = 7, h = 3.5)
  message("  ism/ism_fsm_ratio")
}

# ── ISM: Redundancy per reference isoform ────────────────────────────────────

fig_ism_redundancy <- function(ism, ism_dir) {
  ism_per_ref <- ism |>
    group_by(associated_transcript) |>
    summarise(n_talon = n_distinct(isoform), .groups = "drop") |>
    mutate(bin = cut(n_talon, breaks = c(1, 2, 3, 5, 10, 20, 50, 100, 500),
                     right = FALSE, include.lowest = TRUE))

  p <- ism_per_ref |>
    count(bin) |>
    ggplot(aes(x = bin, y = n)) +
    geom_col(fill = "#67A9CF", colour = "white", linewidth = 0.2) +
    scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k")) +
    labs(title = "ISM redundancy: TALON IDs per reference isoform",
         x = "Unique TALON IDs per reference isoform",
         y = "Reference isoforms") +
    theme_manuscript +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6))
  save_fig(p, file.path(ism_dir, "ism_redundancy"), w = 5, h = 3.5)
  message("  ism/ism_redundancy")
}

# ── ISM: Fragment truncation offsets ─────────────────────────────────────────

fig_ism_truncation_offsets <- function(ism, ism_dir) {
  offsets <- bind_rows(
    ism |> filter(subcategory == "3prime_fragment") |>
      transmute(offset = diff_to_TSS, type = "3' frag: TSS offset"),
    ism |> filter(subcategory == "5prime_fragment") |>
      transmute(offset = diff_to_TTS, type = "5' frag: TTS offset")
  ) |>
    filter(!is.na(offset), abs(offset) < 50000)

  p <- ggplot(offsets, aes(x = offset, colour = type, fill = type)) +
    geom_density(alpha = 0.15, linewidth = 0.6) +
    scale_colour_manual(values = c("3' frag: TSS offset" = "#66c2a5",
                                   "5' frag: TTS offset" = "#fc8d62"), name = NULL) +
    scale_fill_manual(values = c("3' frag: TSS offset" = "#66c2a5",
                                 "5' frag: TTS offset" = "#fc8d62"), name = NULL) +
    xlim(-20000, 5000) +
    labs(title = "Fragment truncation offsets",
         x = "Offset from reference (bp)", y = "Density") +
    theme_manuscript +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  save_fig(p, file.path(ism_dir, "ism_truncation_offsets"), w = 5, h = 3.5)
  message("  ism/ism_truncation_offsets")
}

# ── ISM: Length completeness by subcategory ──────────────────────────────────

fig_ism_length_completeness <- function(ism, ism_dir) {
  ism_len <- ism |>
    filter(ref_length > 0) |>
    mutate(length_ratio = length / ref_length) |>
    filter(length_ratio >= 0, length_ratio <= 2)
  subcats_main <- c("3prime_fragment", "5prime_fragment", "mono-exon",
                     "internal_fragment", "intron_retention")
  subcats_present <- subcats_main[subcats_main %in% unique(ism_len$subcategory)]

  p <- ism_len |>
    filter(subcategory %in% subcats_present) |>
    mutate(subcategory = factor(subcategory, levels = subcats_present),
           subcat_label = gsub("_", "\n", subcategory)) |>
    ggplot(aes(x = subcat_label, y = length_ratio, fill = subcategory)) +
    geom_violin(linewidth = 0.2, show.legend = FALSE) +
    stat_summary(fun = median, geom = "crossbar", width = 0.4, linewidth = 0.3) +
    geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.3) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "ISM length completeness by subcategory",
         x = NULL, y = "Length ratio (observed / reference)") +
    theme_manuscript +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 6))
  save_fig(p, file.path(ism_dir, "ism_length_completeness"), w = 5, h = 3.5)
  message("  ism/ism_length_completeness")
}

# ── Main ─────────────────────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript visualize_sqanti3.R <sqanti3_dir> [-o output_dir]\n")
  quit(status = 1)
}

sqanti_dir <- normalizePath(args[1])
out_dir <- if ("-o" %in% args) args[which(args == "-o") + 1] else
  file.path(sqanti_dir, "figures")
ism_dir <- file.path(out_dir, "ism")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ism_dir, recursive = TRUE, showWarnings = FALSE)

sample_names <- discover_samples(sqanti_dir)
if (length(sample_names) == 0) {
  stop("No SQANTI3 classification files found in ", sqanti_dir, "/\n",
       "Expected: <dir>/<sample>/<sample>_classification.txt")
}

message("Found ", length(sample_names), " samples in ", sqanti_dir, "/")
for (s in sample_names) message("  ", s)

message("\nLoading classifications...")
data <- load_classifications(sqanti_dir, sample_names)
message("  ", format(nrow(data), big.mark = ","), " transcripts total")

message("\nGenerating figures in ", out_dir, "/")
fig_category_counts(data, sample_names, out_dir)
fig_category_proportions(data, sample_names, out_dir)
fig_subcategory_breakdown(data, sample_names, out_dir)
fig_sample_heatmap(data, sample_names, out_dir)
fig_exon_distributions(data, sample_names, out_dir)

message("\nGenerating ISM figures in ", ism_dir, "/")
ism <- data |> filter(structural_category == "incomplete-splice_match")
fsm <- data |> filter(structural_category == "full-splice_match")
if (nrow(ism) > 0) {
  fig_ism_subcategories(ism, sample_names, ism_dir)
  fig_ism_exon_retention(ism, ism_dir)
  fig_ism_missing_exons(ism, ism_dir)
  fig_ism_fsm_ratio(ism, fsm, sample_names, ism_dir)
  fig_ism_redundancy(ism, ism_dir)
  fig_ism_truncation_offsets(ism, ism_dir)
  fig_ism_length_completeness(ism, ism_dir)
} else {
  message("  No ISM transcripts found, skipping ISM figures")
}

message("\nDone. ", length(sample_names), " samples.")