## filter_CNVs.R
## Reads the WES-snakemake combined CNV table (results/combined/combined_cnvs.tsv,
## CNVkit segments, one row per segment per sample, sID as join key), joins project
## metadata, drops segments at the sample's own resolved ploidy (not always 2 —
## see section 3b/6), and labels the remaining segments in a FILTER column.
##
## Inputs:
##   <metadata_dir>/metadata.csv          — patient/sample metadata (Disease, ...)
##   <outdir>/<flags_genelist>             — genes to flag; config.yaml's flags_genelist
##                                          names the file (default "flags_genes.csv")
##   <combined_dir>/combined_cnvs.tsv     — WES-snakemake combined CNV table
##   config.yaml                          — cnv_filter_params block sets the artifact thresholds
##
## Outputs:
##   <outdir>/analysis-CNVs/combined_cnvs.qs2  — non-diploid segments with a FILTER column
##
## FILTER labels written:
##   PASS                    — survives all filters
##   artifact_zero_cov       — called cn = 0 while coverage and heterozygous SNPs are
##                             still present (cnv_filter_params)
##   low_probe_small_segment — few CNVkit probes AND small (cnv_filter_params); large
##                             low-probe segments are kept — see section 8
##   FLAGS_all_genes         — every gene in the segment is on the flags list
##
## Segments at the sample's expected normal copy number are dropped.
## The baseline is the sample's own resolved ploidy (joined in from combined_qc.tsv, 
## section 3b — cnvkit's --ploidy, PureCN-fit or samplesheet-known):
## on autosomes that is `ploidy`; on chrX/chrY it also depends on gender
## (male chrX/chrY are normally hemizygous, i.e. ploidy/2; female chrY is normally
## absent, cn == 0) — see section 6.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(qs2)
})

# Sourced by load_WES_results.R, which resolves every path and the config; this script
# has no standalone defaults.
.required <- c("rootdir", "cfg", "combined_dir", "outdir", "metadata_dir")
.missing  <- .required[!vapply(.required, exists, logical(1))]
if (length(.missing) > 0) {
  stop("missing: ", paste(.missing, collapse = ", "),
       " — run via load_WES_results.R --config <config.yaml>", call. = FALSE)
}
rm(.required, .missing)
cnv_outdir <- file.path(outdir, "analysis-CNVs")
dir.create(cnv_outdir, showWarnings = FALSE, recursive = TRUE)

log_con <- file(file.path(cnv_outdir, "filter_CNVs.log"), open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

# ============================================================================
# 1. Load metadata
# ============================================================================

metadata <- read.csv(file.path(metadata_dir, "metadata.csv"), stringsAsFactors = FALSE)

cat("Metadata rows:", nrow(metadata), "\n")
cat("group1:", paste(names(table(metadata$group1)),
    table(metadata$group1), sep = "=", collapse = ", "), "\n")

# ============================================================================
# 2. Load combined CNV table
# ============================================================================

cnv_path <- file.path(combined_dir, "combined_cnvs.tsv")
cnv_combined <- as.data.frame(fread(cnv_path, sep = "\t", na.strings = c("", "NA")))

cat("\nLoaded", cnv_path, "\n")
cat("Combined CNV table:", nrow(cnv_combined), "segments from",
    n_distinct(cnv_combined$ID), "patients\n")

# ============================================================================
# 3. Left-join project metadata (Disease, calling_mode, group1, group2 if set, QC)
# ============================================================================

meta_add <- metadata %>% select(sID, Disease, calling_mode, group1, any_of("group2"), QC)
cnv_combined <- cnv_combined %>% left_join(meta_add, by = "sID")

n_no_meta <- sum(!cnv_combined$sID %in% metadata$sID)
if (n_no_meta > 0) {
  cat("Segments from samples not in metadata (dropped):", n_no_meta, "\n")
  cat("  Missing sIDs:", paste(unique(cnv_combined$sID[!cnv_combined$sID %in% metadata$sID]),
      collapse = ", "), "\n")
  cnv_combined <- cnv_combined %>% filter(sID %in% metadata$sID)
}

# ============================================================================
# 3b. Join resolved ploidy (the value actually passed to cnvkit --ploidy)
# ============================================================================
# combined_qc.tsv is written by the same WES-snakemake run as combined_cnvs.tsv,
# from the resolve_purity_source sidecar (purecn.smk) — the single source of
# truth cnvkit_call itself used, so this can't drift from what actually produced
# `cn`. Falls back to 2 (the pipeline's own default) for any sID absent from it.

qc_path <- file.path(combined_dir, "combined_qc.tsv")
qc_add <- as.data.frame(fread(qc_path, sep = "\t", na.strings = c("", "NA"))) %>%
  select(sID, ploidy, ploidy_source)
cnv_combined <- cnv_combined %>% left_join(qc_add, by = "sID")

n_no_ploidy <- sum(is.na(cnv_combined$ploidy))
if (n_no_ploidy > 0) {
  warning(n_no_ploidy, " segments have no resolved ploidy (sID missing from ",
          "combined_qc.tsv); normal-ploidy baseline falls back to diploid (2) ",
          "for these: ", paste(unique(cnv_combined$sID[is.na(cnv_combined$ploidy)]),
          collapse = ", "), call. = FALSE)
}
cnv_combined$ploidy[is.na(cnv_combined$ploidy)] <- 2L
cnv_combined$ploidy <- as.integer(round(cnv_combined$ploidy))

# Every segment starts as PASS; sections 5 and 8 relabel the ones they catch.
cnv_combined$FILTER <- "PASS"

# Reorder: metadata first, then CNV columns in logical order
meta_cols <- intersect(c("group1", "group2", "Disease", "ID", "sID", "sample_type", "gender",
                         "capture_kit", "tumor_fraction", "ploidy", "ploidy_source",
                         "calling_mode", "QC"),
                       names(cnv_combined))
cnv_col_order <- c(
  "chromosome", "start", "end", "gene", "FILTER",
  "cn", "cn1", "cn2",
  "log2", "baf",
  "ci_lo", "ci_hi",
  "depth", "probes", "nbaf", "weight", "p_ttest"
)
remaining <- setdiff(names(cnv_combined), c(meta_cols, cnv_col_order))
cnv_combined <- cnv_combined[, c(
  meta_cols,
  intersect(cnv_col_order, names(cnv_combined)),
  remaining
)]

# ============================================================================
# 4. Assess statistical columns before filtering
# ============================================================================

cat("\n=== Statistical assessment (all segments, pre-filter) ===\n")
cat("cn distribution:\n")
print(table(cnv_combined$cn, useNA = "ifany"))

cat("\nprobes distribution (quantiles):\n")
print(quantile(cnv_combined$probes, probs = c(0, .05, .10, .25, .50, .75, .90, .95, 1),
               na.rm = TRUE))

cat("\np_ttest distribution (quantiles):\n")
print(quantile(cnv_combined$p_ttest, probs = c(0, .05, .25, .50, .75, .95, .99, 1),
               na.rm = TRUE))

cat("\nweight distribution (quantiles):\n")
print(quantile(cnv_combined$weight, probs = c(0, .05, .25, .50, .75, .95, .99, 1),
               na.rm = TRUE))

# ============================================================================
# 5. Artifact segments — cn = 0 called where coverage is still present
# ============================================================================
# A homozygous deletion must have lost its coverage. CNVkit nevertheless calls
# cn = 0 over off-target and low-mappability regions whose underlying sequencing
# coverage is untouched, often at a log2 far below any biological ratio but by no
# means always. Such segments span tens of Mb and recur at identical coordinates
# across unrelated patients, which makes them score as concordant between a
# patient and their models.
#
# The entry criterion is therefore the cn = 0 call itself, and two independent
# signals separate the impossible calls from genuine homozygous deletions:
#   depth relative to the sample's diploid baseline — a real homozygous deletion
#   has essentially no reads, an artifact keeps near-diploid coverage; and
#   nbaf, the number of heterozygous SNPs in the segment — impossible inside a
#   homozygous deletion, plentiful in an artifact.
# A segment is only called an artifact when both signals agree, so anything
# either signal calls real is retained. nbaf = NA counts as no het SNPs present.
#
# This runs while the sample's expected-ploidy segments are still in the table:
# the baseline is the median depth across them (cn == ploidy, not a flat 2 — see
# section 3b; a WGD sample's bulk sits at cn == ploidy, not cn == 2).

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && nzchar(a)) a else b
}

artifact_rel_depth_min <- as.numeric(cfg$cnv_filter_params$artifact_rel_depth_min)

# Per-sample coverage baseline at the sample's own expected ploidy. A sample
# without any cn == ploidy segment falls back to its median depth over all
# segments.
.baseline <- cnv_combined %>%
  group_by(sID) %>%
  summarise(
    n_diploid      = sum(cn == ploidy & !is.na(cn)),
    diploid_depth  = median(depth[cn == ploidy & !is.na(cn)], na.rm = TRUE),
    fallback_depth = median(depth, na.rm = TRUE),
    .groups = "drop"
  )

.no_diploid <- .baseline$n_diploid == 0 | is.na(.baseline$diploid_depth)
if (any(.no_diploid)) {
  warning("Samples without any cn == ploidy segments; baseline falls back to all-segment ",
          "median depth: ", paste(.baseline$sID[.no_diploid], collapse = ", "),
          call. = FALSE)
  cat("\nSamples without cn == ploidy segments (baseline from all-segment median depth):",
      paste(.baseline$sID[.no_diploid], collapse = ", "), "\n")
  .baseline$diploid_depth[.no_diploid] <- .baseline$fallback_depth[.no_diploid]
}

cnv_combined <- cnv_combined %>%
  left_join(.baseline %>% select(sID, diploid_depth), by = "sID") %>%
  mutate(rel_depth = depth / diploid_depth)

# A missing cn, rel_depth or nbaf leaves the segment retained.
.artifact_mask <-
  !is.na(cnv_combined$cn)        & cnv_combined$cn == 0 &
  !is.na(cnv_combined$rel_depth) & cnv_combined$rel_depth >= artifact_rel_depth_min &
  !is.na(cnv_combined$nbaf)      & cnv_combined$nbaf > 0

cnv_combined$FILTER[.artifact_mask & cnv_combined$FILTER == "PASS"] <- "artifact_zero_cov"

.zero_cn <- !is.na(cnv_combined$cn) & cnv_combined$cn == 0

cat("\n=== Artifact segments ===\n")
cat(sprintf("Rule: cn == 0 AND depth/diploid-baseline >= %g AND nbaf > 0\n",
            artifact_rel_depth_min))
cat(sprintf("Segments called cn == 0: %d\n", sum(.zero_cn)))
cat(sprintf("Segments flagged artifact_zero_cov: %d in %d samples\n",
            sum(.artifact_mask), n_distinct(cnv_combined$sID[.artifact_mask])))
cat(sprintf("Genomic span flagged: %.2f Gb\n",
            sum(cnv_combined$end[.artifact_mask] - cnv_combined$start[.artifact_mask]) / 1e9))
cat(sprintf("Retained as genuine homozygous deletions: %d\n",
            sum(.zero_cn & !.artifact_mask)))
cat(sprintf("Diploid baseline depth across samples: min %.1fx, median %.1fx, max %.1fx\n",
            min(.baseline$diploid_depth, na.rm = TRUE),
            median(.baseline$diploid_depth, na.rm = TRUE),
            max(.baseline$diploid_depth, na.rm = TRUE)))

if (any(.artifact_mask)) {
  .artifact_coord <- sprintf("%s:%.1f-%.1fMb",
                             cnv_combined$chromosome[.artifact_mask],
                             cnv_combined$start[.artifact_mask] / 1e6,
                             cnv_combined$end[.artifact_mask] / 1e6)
  cat("\nMost recurrent artifact coordinates (samples per coordinate):\n")
  print(head(sort(table(.artifact_coord), decreasing = TRUE), 10))
  rm(.artifact_coord)
}

rm(.baseline, .no_diploid, .artifact_mask, .zero_cn)

# ============================================================================
# 6. Filter segments at the expected normal copy number (no copy number change)
# ============================================================================
# "Normal" is the sample's own resolved ploidy (joined in section 3b), not a flat
# 2 — a whole-genome-duplicated sample (e.g. ploidy 4) sits at cn == 4 across most
# of its genome by construction (cnvkit's --ploidy rescale maps log2 == 0 to
# cn == ploidy), and a flat cn != 2 filter would keep nearly the whole genome as a
# spurious "gain". cn == ploidy is "normal" on autosomes and on chrX in females;
# male chrX/chrY carry half that (hemizygous ancestrally, scaled by any WGD same
# as everything else) and female chrY is always absent regardless of ploidy.
# Filtering a flat cn != 2 on top of that also leaves every male sample's
# chrX/chrY and every female sample's chrY in the table as spurious "loss"
# segments the size of whole chromosome arms.

cnv_combined <- cnv_combined %>%
  mutate(cn_normal = dplyr::case_when(
    !chromosome %in% c("chrX", "chrY") ~ ploidy,
    chromosome == "chrX" & gender == "f" ~ ploidy,
    chromosome == "chrX" & gender == "m" ~ pmax(1L, as.integer(round(ploidy / 2))),
    chromosome == "chrY" & gender == "m" ~ pmax(1L, as.integer(round(ploidy / 2))),
    chromosome == "chrY" & gender == "f" ~ 0L,
    TRUE ~ NA_integer_
  ))

n_unknown_ploidy <- sum(is.na(cnv_combined$cn_normal))
if (n_unknown_ploidy > 0) {
  warning(n_unknown_ploidy, " chrX/chrY segments have an unrecognised gender value; ",
          "retained regardless of cn (expected ploidy unknown).", call. = FALSE)
}

n_before <- nrow(cnv_combined)
n_sexchrom_normal <- sum(cnv_combined$chromosome %in% c("chrX", "chrY") &
                          !is.na(cnv_combined$cn_normal) &
                          cnv_combined$cn == cnv_combined$cn_normal)
cnv_combined <- cnv_combined %>%
  filter(is.na(cn_normal) | cn != cn_normal) %>%
  select(-cn_normal)
n_after <- nrow(cnv_combined)

cat("\n=== Normal-ploidy filter ===\n")
cat("Segments before:", n_before, "\n")
cat("Segments removed (cn == expected ploidy):", n_before - n_after, "\n")
cat("  of which chrX/chrY at the gender-expected ploidy:", n_sexchrom_normal, "\n")
cat("Segments remaining:", n_after, "\n")

cat("\ncn distribution after filter:\n")
print(table(cnv_combined$cn, useNA = "ifany"))

# ============================================================================
# 7. Assess p_ttest and probes in non-diploid segments
# ============================================================================

cat("\n=== Statistical assessment (non-diploid segments) ===\n")
cat("\nprobes distribution:\n")
print(quantile(cnv_combined$probes, probs = c(0, .05, .10, .25, .50, .75, .90, .95, 1),
               na.rm = TRUE))
cat("Segments with probes < 4:", sum(cnv_combined$probes < 4, na.rm = TRUE), "\n")
cat("Segments with probes < 10:", sum(cnv_combined$probes < 10, na.rm = TRUE), "\n")

cat("\np_ttest distribution:\n")
print(quantile(cnv_combined$p_ttest, probs = c(0, .05, .25, .50, .75, .95, .99, 1),
               na.rm = TRUE))
cat("Segments with p_ttest >= 0.05:", sum(cnv_combined$p_ttest >= 0.05, na.rm = TRUE), "\n")
cat("Segments with p_ttest >= 0.01:", sum(cnv_combined$p_ttest >= 0.01, na.rm = TRUE), "\n")

cat("\nweight distribution:\n")
print(quantile(cnv_combined$weight, probs = c(0, .05, .25, .50, .75, .95, .99, 1),
               na.rm = TRUE))

# NOTE on p_ttest and probes:
# Neither column is used as a filter directly here; both are reported so the
# distributions stay visible in the log. p_ttest measures internal homogeneity,
# not whether log2 deviates from 0, and after dropping cn == 2 nearly all segments
# have a very small p_ttest because CBS segmentation already selects internally
# consistent segments — so it separates almost nothing. 'probes' alone is not
# principled either: a low probe count means as little as "spans few on-target
# exons" for a real arm-scale event as it does "no evidence" for an isolated
# noise call — see section 9, which combines it with segment width instead.

# ============================================================================
# 8. Low-probe small-segment filter — CNVkit noise on isolated exons
# ============================================================================
# A low probe count alone doesn't separate real calls from noise: whole-arm
# gains/losses can have just as few *on-target* exonic probes as noise does,
# simply because they span mostly gene-poor sequence — they just carry
# overwhelming independent support elsewhere (hundreds to thousands of
# heterozygous SNPs in nbaf, near-zero p_ttest), which probes doesn't capture.
# Width does separate the two: among PASS segments with probes below the
# threshold, the width distribution has a clean gap (~0.8-1.1 Mb, no segments
# fall in between) between single/double-exon noise and arm-scale calls, so
# adding a width cap catches only the former.

low_probe_probes_max   <- as.numeric(cfg$cnv_filter_params$low_probe_probes_max)
low_probe_width_max_mb <- as.numeric(cfg$cnv_filter_params$low_probe_width_max_mb)

.segment_width_mb <- (cnv_combined$end - cnv_combined$start) / 1e6
.low_probe_small_mask <-
  !is.na(cnv_combined$probes) & cnv_combined$probes < low_probe_probes_max &
  .segment_width_mb < low_probe_width_max_mb

cnv_combined$FILTER[.low_probe_small_mask & cnv_combined$FILTER == "PASS"] <- "low_probe_small_segment"

cat("\n=== Low-probe small-segment filter ===\n")
cat(sprintf("Rule: probes < %g AND width < %g Mb\n", low_probe_probes_max, low_probe_width_max_mb))
cat(sprintf("Segments flagged: %d\n", sum(.low_probe_small_mask)))

rm(.segment_width_mb, .low_probe_small_mask)

# ============================================================================
# 9. FLAGS filter — label segments whose entire gene content is flagged
# ============================================================================
# Unlike SNVs (1 variant = 1 gene), CNV segments span many genes. A segment is
# only labelled if EVERY gene in it is on the flags list, i.e. there are no
# "real" genes present. Segments with gene == "-" (intergenic) stay PASS.

flags_genelist <- if (exists("cfg") && !is.null(cfg$flags_genelist) && nzchar(cfg$flags_genelist)) {
  cfg$flags_genelist
} else {
  "flags_genes.csv"
}
flags <- read.csv(file.path(outdir, flags_genelist), stringsAsFactors = FALSE)
exact_flags  <- flags %>% filter(match_type == "exact")
prefix_flags <- flags %>% filter(match_type == "prefix")

prefix_pattern <- if (nrow(prefix_flags) > 0) {
  paste0("^(", paste(prefix_flags$gene, collapse = "|"), ")[0-9]")
} else {
  NA_character_
}

is_flagged_gene <- function(g) {
  exact_hit  <- g %in% exact_flags$gene
  prefix_hit <- !is.na(prefix_pattern) && grepl(prefix_pattern, g)
  exact_hit | prefix_hit
}

is_all_flagged <- function(gene_string) {
  genes <- trimws(unlist(strsplit(gene_string, ",")))
  genes <- genes[nzchar(genes) & genes != "-"]
  if (length(genes) == 0) return(FALSE)   # intergenic — keep
  all(sapply(genes, is_flagged_gene))
}

flagged_mask <- sapply(cnv_combined$gene, is_all_flagged)

cat("\n=== FLAGS filter ===\n")
cat("Segments flagged (all genes in flags list):", sum(flagged_mask), "\n")
if (sum(flagged_mask) > 0) {
  cat("Flagged segments by gene:\n")
  print(sort(table(cnv_combined$gene[flagged_mask]), decreasing = TRUE))
}

cnv_combined$FILTER[flagged_mask & cnv_combined$FILTER == "PASS"] <- "FLAGS_all_genes"

cat("\n=== FILTER breakdown (non-diploid segments) ===\n")
print(table(cnv_combined$FILTER))
cat("Segments PASS:", sum(cnv_combined$FILTER == "PASS"), "of", nrow(cnv_combined), "\n")

# ============================================================================
# 10. Save outputs
# ============================================================================

qs_save(cnv_combined, file.path(cnv_outdir, "combined_cnvs.qs2"))
cat("\nSaved combined_cnvs.qs2:", nrow(cnv_combined), "rows x",
    ncol(cnv_combined), "columns\n")

cat("\nDone!\n")

w <- warnings()
if (length(w) > 0) {
  cat("\n=== Warnings ===\n")
  print(w)
}

sink(type = "message")
sink()
close(log_con)
