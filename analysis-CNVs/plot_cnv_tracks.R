## plot_cnv_tracks.R
## Genome-wide CNV tracks, the linear counterpart of model-comparisons' Circos rings: one row
## per sample (QC == "PASS" in metadata.csv), chr1-22, X, Y on the x axis. Panel widths
## scale with sqrt(chromosome length), so the small chromosomes stay readable.
##
## Copy number is relative to the sample's expected copy number at that locus (cn_normal,
## as in filter_CNVs.R section 6): the resolved ploidy on autosomes and female chrX,
## ploidy / 2 on male chrX/chrY. Female chrY is not drawn.
##   - PASS segments of combined_cnvs.qs2 (filter_CNVs.R keeps only cn != cn_normal) are
##     filled by cn - cn_normal, clamped to -2..+2, with the Circos palette.
##   - Segments at cn_normal come from WES-snakemake's combined_cnvs.tsv (filter_CNVs.R
##     drops them) and are drawn gray90 underneath.
##   - Regions without a segment (gaps, filtered segments) show the gray70 background.
## A thin cytoband row (UCSC gieStain, config `cytoband`) sits above the first sample.
##
## Inputs:
##   config.yaml  — results_dir, metadata_dir, cytoband
##   <outdir>/analysis-CNVs/combined_cnvs.qs2
##
## Output:
##   <outdir>/analysis-CNVs/cnv_tracks.png
##
## Usage:
##   Rscript plot_cnv_tracks.R --config results/<date>_WES-filtering/config.yaml

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggh4x)
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

cfg <- read_wes_config(config_arg("plot_cnv_tracks.R"))

combined_dir <- file.path(normalizePath(cfg$results_dir, mustWork = TRUE), "combined")
outdir <- dirname(cfg$config_path)
cnv_outdir <- file.path(outdir, "analysis-CNVs")

# GRCh38 chromosome lengths
chr_len <- data.table(
  chromosome = paste0("chr", c(1:22, "X", "Y")),
  len = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
          145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
          101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468,
          156040895, 57227415))

# Relative CN buckets and colours (model-comparisons' cnv_palette; neutral = its gray90
# at-ploidy baseline)
rel_levels <- c("≤-2", "-1", "0", "+1", "≥+2")
rel_colors <- setNames(c("#2D6FAB", "#92AFCF", "gray90", "#E6AAAA", "#EB3B25"), rel_levels)

filt <- as.data.table(qs2::qs_read(file.path(cnv_outdir, "combined_cnvs.qs2")))
metadata <- fread(file.path(cfg$metadata_dir, "metadata.csv"))
samples <- metadata[QC == "PASS"][order(ID, sID), sID]
sample_info <- unique(filt[sID %in% samples, .(sID, ploidy, gender)])
sample_info[, label := sprintf("%s (%sN)", sID, ploidy)]
labels <- sample_info[match(samples, sID), label]
rows <- c("", labels)   # "" = cytoband row

cols <- c("sID", "chromosome", "start", "end", "cn")
prep <- function(d) {
  d <- merge(d[sID %in% samples & chromosome %in% chr_len$chromosome], sample_info, by = "sID")
  d[, cn_normal := fcase(
    !chromosome %in% c("chrX", "chrY"), as.numeric(ploidy),
    chromosome == "chrX" & gender == "f", as.numeric(ploidy),
    gender == "m", pmax(1, round(ploidy / 2)),
    default = NA_real_)]
  d <- d[!is.na(cn) & !is.na(cn_normal) & !(chromosome == "chrY" & cn_normal == 0)]
  d[, `:=`(rel = factor(rel_levels[pmin(pmax(cn - cn_normal, -2), 2) + 3], levels = rel_levels),
           label = factor(label, levels = rows),
           chromosome = factor(chromosome, levels = chr_len$chromosome))]
}
altered <- prep(filt[FILTER == "PASS", ..cols])[cn != cn_normal]
neutral <- prep(fread(file.path(combined_dir, "combined_cnvs.tsv"), select = cols))[cn == cn_normal]
cat(sprintf("Samples: %d; segments: %d altered (PASS), %d at expected cn\n",
            length(samples), nrow(altered), nrow(neutral)))

band_colors <- c(gneg = "white", gpos25 = "#C8C8C8", gpos50 = "#969696", gpos75 = "#646464",
                 gpos100 = "black", gvar = "#DCDCDC", stalk = "#647FA4", acen = "#C0392B")
bands <- fread(cfg$cytoband, col.names = c("chromosome", "start", "end", "band", "stain"))
bands <- bands[chromosome %in% chr_len$chromosome]
bands[, `:=`(label = factor("", levels = rows),
             chromosome = factor(chromosome, levels = chr_len$chromosome),
             col = band_colors[stain])]

# invisible points pin each chromosome panel to its full length
bounds <- chr_len[, .(label = factor(labels[1], levels = rows),
                      chromosome = factor(chromosome, levels = chr_len$chromosome),
                      x = c(0, len)), by = .(chr = chromosome)]

g <- ggplot() +
  geom_blank(data = bounds, aes(x = x)) +
  geom_rect(data = bands, aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = I(col)),
            colour = NA, linewidth = 0) +
  geom_rect(data = neutral, aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = rel),
            colour = NA, linewidth = 0) +
  geom_rect(data = altered, aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = rel),
            colour = NA, linewidth = 0) +
  scale_fill_manual("Relative CN", values = rel_colors, breaks = rev(rel_levels), drop = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(label ~ chromosome, switch = "y", scales = "free_x",
             labeller = labeller(chromosome = function(x) sub("^chr", "", x))) +
  force_panelsizes(cols = sqrt(chr_len$len), rows = c(0.25, rep(1, length(labels)))) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(
    panel.background = element_rect(fill = "gray70"),
    panel.grid = element_blank(),
    panel.spacing = unit(0, "lines"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    strip.text.y.left = element_text(angle = 0),
    strip.text.x = element_text(size = 10)
  )

out <- file.path(cnv_outdir, "cnv_tracks.png")
ggsave(out, g, width = 15, height = 0.35 * length(samples) + 1.1, dpi = 300, bg = "white")
cat("Wrote", out, "\n")
