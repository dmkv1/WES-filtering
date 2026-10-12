## plot_snv_filters.R
## Diagnostic plots for filter_SNVs.R, sharing one FILTER-status colour scheme:
##   1. filters_SNV.png: per-sample stacked bar of variant counts by FILTER status (see
##      filter_SNVs.R header for the full list), log10 x-axis, faceted by metadata.csv's
##      group1 (whatever column config.yaml's group1 names). Two bars per sample: top, the
##      quality-clean calls (PASS_coding, PASS_other, low_priority_consequence); bottom, all
##      other calls. Red dashed line at snv_qc_threshold: a sample whose top bar ends right of
##      it is high_SNV_count.
##   2. af_SNV_<Disease>.png: per-sample tumor AF of every variant (quasirandom dots)
##      coloured by FILTER status, faceted by patient ID, one file per Disease (af_SNV.png for none). Dashed
##      lines at AF 0.2 and 0.5; N = variant count of the sample.
##
##   - FLAGS_large_gene / FLAGS_polymorphic / FLAGS_repetitive are pooled into one FLAGS
##     bucket.
##   - FilterMutectCalls labels kept by WES-snakemake (keep_filter_labels, germline and
##     strand_bias pass-through) are pooled: any set containing germline -> M2_germline,
##     strand_bias alone -> M2_strand_bias, other sets of clustered_events / haplotype /
##     strand_bias -> M2_clustered_haplotype.
##   - PASS is split into PASS_coding / PASS_other via is_coding_consequence(). PASS is
##     already allowlist-restricted to coding + 3' UTR, so PASS_other is predominantly
##     3' UTR, not an unfiltered noncoding bucket.
##
## Bar label: sID, calling mode (TN/TO), capture kit, group2 when set, and " [excluded]"
## for samples with QC == "excluded".
##
## Inputs:
##   config.yaml                              -- its directory is outdir
##   <outdir>/analysis-SNVs/combined_snvs.qs2 -- filter_SNVs.R output (requires ./runall --load)
##
## Outputs:
##   <outdir>/analysis-SNVs/filters_SNV.png
##   <outdir>/analysis-SNVs/af_SNV_<Disease>.png (af_SNV.png without a Disease)
##
## Usage:
##   Rscript plot_snv_filters.R --config results/<date>_WES-filtering/config.yaml

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggbeeswarm)
  library(ggh4x)
  library(scales)
  library(yaml)
  library(qs2)
})

.script_path <- if (sys.nframe() > 0 && !is.null(sys.frame(1)$ofile)) {
  sys.frame(1)$ofile
} else {
  sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
}
.this_dir <- normalizePath(dirname(.script_path))
rm(.script_path)

source(file.path(.this_dir, "snv_consequence_helpers.R"))
source(file.path(.this_dir, "..", "config.R"))

cfg <- read_wes_config(config_arg("plot_snv_filters.R"))
config_path <- cfg$config_path

outdir    <- dirname(config_path)
plots_dir <- file.path(outdir, "analysis-SNVs")

snvs <- as.data.table(qs_read(file.path(outdir, "analysis-SNVs", "combined_snvs.qs2")))

# ── Per-sample label: "sID  calling_mode  capture_kit  [group2]  [excluded]" ────────

calling_mode_abbrev <- c("tumor-normal" = "TN", "tumor-only" = "TO")

snvs[, sample_label := trimws(paste(
  sID,
  calling_mode_abbrev[calling_mode],
  sub("\\+UTR$", "", capture_kit),
  if ("group2" %in% names(snvs)) ifelse(is.na(group2), "", as.character(group2)) else "",
  ifelse(QC %in% "excluded", "[excluded]", "")
))]

# ── Bucket FILTER into the plot's display categories ────────────────────────────────

snvs[, display_status := dplyr::case_when(
  startsWith(FILTER, "FLAGS_")                             ~ "FLAGS",
  grepl("(^|;)germline(;|$)", FILTER)                      ~ "M2_germline",
  FILTER == "strand_bias"                                  ~ "M2_strand_bias",
  grepl("^(clustered_events|haplotype|strand_bias)(;(clustered_events|haplotype|strand_bias))*$",
        FILTER)                                            ~ "M2_clustered_haplotype",
  FILTER == "PASS" & is_coding_consequence(Consequence)    ~ "PASS_coding",
  FILTER == "PASS"                                         ~ "PASS_other",
  TRUE                                                      ~ FILTER
)]

status_levels <- c("M2_germline", "M2_strand_bias", "M2_clustered_haplotype", "FLAGS", "low_alt_reads", "low_AF", "tumor_only_low_AF",
                    "germline_ctrl", "germline_gnomad", "low_priority_consequence",
                    "PASS_other", "PASS_coding")
status_colors <- c(
  M2_germline              = "#5E3C99",
  M2_strand_bias           = "#B2ABD2",
  M2_clustered_haplotype   = "#1B9E77",
  FLAGS                    = "#AAAAAA",
  low_alt_reads            = "#F4A261",
  low_AF                   = "#E76F51",
  tumor_only_low_AF        = "#F1C40F",
  germline_ctrl            = "#C1121F",
  germline_gnomad          = "#780000",
  low_priority_consequence = "#8D99AE",
  PASS_other               = "#90C6E8",
  PASS_coding              = "#2471A3"
)

unknown_status <- setdiff(unique(snvs$display_status), status_levels)
if (length(unknown_status) > 0) {
  stop("FILTER value(s) not recognized by plot_snv_filters.R's status_levels/status_colors: ",
       paste(unknown_status, collapse = ", "),
       " -- update this script if filter_SNVs.R gained a new FILTER label.", call. = FALSE)
}

# filter_SNVs.R's .snv_qc_threshold: more quality-clean calls than this -> high_SNV_count
snv_qc_threshold <- 2000L
quality_clean <- c("PASS_coding", "PASS_other", "low_priority_consequence")

cat("=== plot_snv_filters.R ===\n")
cat("Patients:", uniqueN(snvs$ID), " Samples:", uniqueN(snvs$sID),
    " Variants:", nrow(snvs), "\n")
print(table(snvs$display_status))

# Order rows within each group1 facet by ascending total variant count (top = most)
sample_order <- snvs[, .N, by = .(sample_label, group1)][order(group1, N), sample_label]
snvs[, sample_f := factor(sample_label, levels = sample_order)]
snvs[, display_status_f := factor(display_status, levels = rev(status_levels))]
snvs[, group1_f := factor(group1, levels = sort(unique(group1)))]

# Pre-aggregate with explicit cumulative xmin/xmax for log-scale stacked bars
# (geom_bar starts from x=0, which is invalid on a log scale)
bar_dt <- snvs[, .N, by = .(sample_f, group1_f, display_status_f)]
bar_dt[, display_status_f := factor(display_status_f, levels = status_levels)]
bar_dt[, clean := display_status_f %in% quality_clean]
setorder(bar_dt, sample_f, display_status_f)
bar_dt[, xmax := cumsum(N),              by = .(sample_f, clean)]
bar_dt[, xmin := shift(xmax, fill = 1L), by = .(sample_f, clean)]
# Quality-clean bar above, the rest below, within the sample's row
bar_dt[, ymid := as.numeric(sample_f) + fifelse(clean, 0.2, -0.2)]

p <- ggplot(bar_dt) +
  geom_rect(aes(
    xmin = xmin, xmax = xmax,
    ymin = ymid - 0.18,
    ymax = ymid + 0.18,
    fill = display_status_f
  )) +
  geom_vline(xintercept = c(30, 100, 300, 1000),
             linetype = "dashed", color = "grey30", linewidth = 0.35) +
  geom_vline(xintercept = snv_qc_threshold,
             linetype = "dashed", color = "red", linewidth = 0.6) +
  scale_x_log10(
    breaks = c(10, 30, 100, 300, 1000, 3000),
    labels = scales::comma
  ) +
  scale_y_continuous(
    breaks = seq_along(levels(snvs$sample_f)),
    labels = levels(snvs$sample_f),
    expand = expansion(add = 0.5)
  ) +
  scale_fill_manual(
    name   = "Filter status",
    values = status_colors,
    breaks = rev(status_levels)
  ) +
  facet_wrap2(~ group1_f, ncol = 1, scales = "free_y", strip.position = "top",
              axes = "x") +
  force_panelsizes(rows = snvs[, uniqueN(sample_label), by = group1_f][
                               match(levels(snvs$group1_f), group1_f), V1]) +
  labs(y = NULL, x = "Variant count (log10)") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.y        = element_text(size = 8),
    axis.text.x        = element_text(size = 8),
    axis.title.x       = element_text(size = 12),
    strip.text.x       = element_text(angle = 0, size = 12),
    legend.position    = "bottom",
    plot.margin        = margin(5.5, 60, 5.5, 5.5),
    panel.grid.major.y = element_blank()
  )

out_plot <- file.path(plots_dir, "filters_SNV.png")
n_bars <- uniqueN(snvs$sample_label)
ggsave(out_plot, p, width = 8, height = min(40, 3 + 0.18 * n_bars), dpi = 300)
cat("\nWrote", out_plot, "\n")

# ── AF dot plots: one file per Disease, faceted by patient ID ───────────────────────

snvs[, sID_label := paste0(sID, ifelse(QC %in% "excluded", " [excluded]", ""))]
snvs[, display_status_f := factor(display_status, levels = status_levels)]

# Samples without a Disease go to af_SNV.png
snvs[, Disease_file := fifelse(is.na(Disease), "", as.character(Disease))]
for (dis in sort(unique(snvs$Disease_file))) {
  df <- snvs[Disease_file == dis]
  setorder(df, ID, sID)
  df[, sID_f := factor(sID_label, levels = unique(sID_label))]
  # PASS drawn last, on top
  df <- df[order(display_status_f %in% c("PASS_coding", "PASS_other"))]
  n_labels <- df[, .(label = sprintf("N=%d", .N)), by = .(ID, sID_f)]

  p_af <- ggplot(df, aes(x = sID_f, y = AF_tumor, color = display_status_f)) +
    geom_quasirandom(size = 0.5, alpha = 0.55, bandwidth = 0.4) +
    geom_hline(yintercept = c(0.2, 0.5), linetype = "dashed", colour = "grey50",
               linewidth = 0.35) +
    geom_text(data = n_labels, aes(x = sID_f, y = 1.01, label = label),
              inherit.aes = FALSE, size = 3, fontface = "bold", vjust = -1,
              colour = "grey20") +
    facet_grid(~ ID, scales = "free_x", space = "free_x") +
    scale_color_manual(name = "Filter status", values = status_colors,
                       breaks = rev(status_levels)) +
    scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25),
                       labels = percent_format(accuracy = 1)) +
    labs(x = NULL, y = "Tumor AF") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
      strip.text      = element_text(face = "bold", size = 8),
      legend.position = "right",
      panel.spacing   = unit(0.3, "lines")
    )

  out_af <- file.path(plots_dir, if (nzchar(dis)) {
    sprintf("af_SNV_%s.png", gsub("[^A-Za-z0-9]", "_", dis))
  } else "af_SNV.png")
  ggsave(out_af, p_af, width = max(6.67, uniqueN(df$sID) * 0.75), height = 4.6,
         dpi = 300, limitsize = FALSE)
  cat("Wrote", out_af, "\n")
}
