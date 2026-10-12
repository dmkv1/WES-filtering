## plot_cn0_rel_depth.R
## Diagnostic plot for the artifact_zero_cov threshold in filter_CNVs.R: distribution
## of rel_depth (depth / sample diploid-baseline depth) across all cn == 0 segments of
## samples not marked QC == "excluded" in metadata.csv, with the configured
## artifact_rel_depth_min cutoff marked. Same baseline/rel_depth logic as
## filter_CNVs.R section 5, kept standalone so it can be re-run without the
## full filter pipeline's side effects (log file, metadata.csv rewrite).
##
## Inputs:
##   config.yaml  — results_dir, metadata_dir, cnv_filter_params$artifact_rel_depth_min
##
## Output:
##   <outdir>/analysis-CNVs/cn0_rel_depth_distribution.png
##
## Usage:
##   Rscript plot_cn0_rel_depth.R --config results/<date>_WES-filtering/config.yaml

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(yaml)
})

.script_path <- if (sys.nframe() > 0 && !is.null(sys.frame(1)$ofile)) {
  sys.frame(1)$ofile
} else {
  sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
}
.this_dir <- normalizePath(dirname(.script_path))
rm(.script_path)
source(file.path(.this_dir, "..", "config.R"))

cfg <- read_wes_config(config_arg("plot_cn0_rel_depth.R"))

combined_dir <- file.path(normalizePath(cfg$results_dir, mustWork = TRUE), "combined")
outdir <- dirname(cfg$config_path)
metadata_dir <- cfg$metadata_dir
cnv_outdir <- file.path(outdir, "analysis-CNVs")
dir.create(cnv_outdir, showWarnings = FALSE, recursive = TRUE)

artifact_rel_depth_min <- as.numeric(cfg$cnv_filter_params$artifact_rel_depth_min)

metadata <- read.csv(file.path(metadata_dir, "metadata.csv"), stringsAsFactors = FALSE)
cnv <- as.data.frame(fread(file.path(combined_dir, "combined_cnvs.tsv"),
                            sep = "\t", na.strings = c("", "NA")))
.included <- if ("QC" %in% names(metadata)) metadata$sID[!metadata$QC %in% "excluded"] else metadata$sID
cat(sprintf("Samples: %d in metadata, %d excluded\n", nrow(metadata),
            nrow(metadata) - length(.included)))
cnv <- cnv %>% filter(sID %in% .included)

baseline <- cnv %>%
  group_by(sID) %>%
  summarise(
    n_diploid      = sum(cn == 2 & !is.na(cn)),
    diploid_depth  = median(depth[cn == 2], na.rm = TRUE),
    fallback_depth = median(depth, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(diploid_depth = ifelse(n_diploid == 0 | is.na(diploid_depth),
                                fallback_depth, diploid_depth))

cnv <- cnv %>%
  left_join(baseline %>% select(sID, diploid_depth), by = "sID") %>%
  mutate(rel_depth = depth / diploid_depth)

zero <- cnv %>% filter(cn == 0, !is.na(rel_depth))
n_flag <- sum(zero$rel_depth >= artifact_rel_depth_min & !is.na(zero$nbaf) & zero$nbaf > 0)

cat(sprintf("cn == 0 segments: %d\n", nrow(zero)))
cat(sprintf("artifact_rel_depth_min: %.2f\n", artifact_rel_depth_min))
cat(sprintf("Flagged artifact_zero_cov (rel_depth >= cutoff & nbaf > 0): %d\n", n_flag))
cat(sprintf("Retained as genuine homozygous deletion: %d\n", nrow(zero) - n_flag))

p <- ggplot(zero, aes(x = pmin(rel_depth, 1.5))) +
  geom_histogram(binwidth = 0.02, boundary = 0, colour = "gray20",
                 fill = "gray70", linewidth = 0.2) +
  geom_vline(xintercept = artifact_rel_depth_min, colour = "firebrick",
             linewidth = 0.6, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 1.5), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "rel_depth (segment depth / sample diploid baseline)", y = "cn = 0 segments") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())

png_path <- file.path(cnv_outdir, "cn0_rel_depth_distribution.png")
ggsave(png_path, p, width = 6, height = 4, dpi = 300)
cat("Wrote", png_path, "\n")
