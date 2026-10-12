## filter_SNVs.R
## Reads the WES-snakemake combined SNV table (results/combined/combined_snvs.tsv,
## VEP-annotated, one row per variant per sample, sID as join key), joins project
## metadata, flags artifact-prone genes, applies quality filters, derives the
## per-sample QC flag from the resulting PASS SNV counts, and saves
## combined_snvs.qs2 for downstream analysis.
##
## Inputs:
##   <metadata_dir>/metadata.csv           — patient/sample metadata (Disease, QC, ...)
##   <combined_dir>/combined_snvs.tsv      — WES-snakemake combined SNV table
##   <outdir>/<flags_genelist>              — genes to flag (exact or prefix match);
##                                           config.yaml's flags_genelist names the file
##                                           (default "flags_genes.csv")
##   config.yaml                           — snv_filter_params block sets the germline_ctrl thresholds
##                                           and the strand-bias and germline rescues
##                                           (strandbias_rescue / germline_rescue: true)
##   <outdir>/SNV_consequence_terms.tsv    — VEP Consequence terms, one row per term, with
##                                           pass_consequence/coding_consequence flag
##                                           columns (see snv_consequence_helpers.R); a variant
##                                           stays PASS if ANY of its (possibly compound, "&"-joined)
##                                           Consequence terms has pass_consequence == TRUE
##
## Outputs:
##   <outdir>/analysis-SNVs/combined_snvs.qs2        — full combined table with FILTER column updated
##                                           and a logical `coding` column (any Consequence
##                                           term has coding_consequence == TRUE)
##   <outdir>/analysis-SNVs/combined_snvs_genelisted.xlsx — genelist subset (all variants, FILTER visible;
##                                           DOMAINS PDB chain mappings collapsed to a count)
##   <outdir>/analysis-SNVs/SNV_counts.csv — per tumor sample: ID, sID, n_snv (all variants),
##                                           n_pass_snv (PASS before the consequence filter),
##                                           snv_qc (high_SNV_count if n_pass_snv > 2000, else
##                                           PASS; ignores exclusion)
##   <metadata_dir>/metadata.csv           — rewritten with the QC column: excluded (kept from
##                                           sample_annotations), else high_SNV_count (>2000
##                                           PASS SNVs), else PASS
##
## FILTER labels written:
##   PASS                    — survives all filters (also the pipeline's own Mutect2 PASS)
##   germline                — Mutect2 germline-only call of a tumor-only sample not rescued
##                             (germline_rescue: true; rescued calls enter the filters as PASS)
##   strand_bias             — Mutect2 strand_bias-only call not rescued (strandbias_rescue: true;
##                             rescued calls enter the filters below as PASS)
##   clustered_events, haplotype (and their ";"-joined sets with strand_bias/germline)
##                           — Mutect2 labels kept by the pipeline (keep_filter_labels) and not
##                             rescued (clustered_rescue: true rescues clustered_events/haplotype-only
##                             calls in gene-list and aSHM genes)
##   FLAGS_high_AF           — PASS call in a high_af_flag_genes gene (HLA class I) at
##                             AF_tumor > high_af_flag_min: mapping artefact
##   FLAGS_large_gene        — gene in large-gene blocklist
##   FLAGS_polymorphic       — gene in polymorphic-gene blocklist
##   FLAGS_repetitive        — gene in repetitive-region blocklist
##   low_alt_reads           — AD_ALT_tumor <= 3
##   low_AF                  — AF_tumor < 0.05
##   tumor_only_low_AF       — tumor-only sample, AF_tumor < 0.2
##   germline_ctrl           — matched-normal AF_ctrl above the absolute cut AND
##                             AF_ctrl/AF_tumor above the ratio cut (snv_filter_params)
##   germline_gnomad         — max(gnomADe_AF, gnomADg_AF) > 0.05
##   low_priority_consequence — none of the variant's Consequence terms has pass_consequence == TRUE
##                             (variant is annotated, never dropped; applied after the QC
##                             metadata step so the >2000 PASS-SNV threshold stays calibrated
##                             against the pre-consequence-filter PASS set)

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(qs2)
})

# Quality-filter label set (used in sections 4 and 5)
quality_filter_labels <- c(
  "low_alt_reads", "low_AF", "tumor_only_low_AF", "germline_ctrl", "germline_gnomad"
)

# Sourced by load_WES_results.R, which resolves every path and the config; this script
# has no standalone defaults.
.required <- c("rootdir", "cfg", "combined_dir", "outdir", "metadata_dir")
.missing  <- .required[!vapply(.required, exists, logical(1))]
if (length(.missing) > 0) {
  stop("missing: ", paste(.missing, collapse = ", "),
       " — run via load_WES_results.R --config <config.yaml>", call. = FALSE)
}
rm(.required, .missing)
snv_outdir <- file.path(outdir, "analysis-SNVs")
dir.create(snv_outdir, showWarnings = FALSE, recursive = TRUE)

source(file.path(rootdir, "analysis-SNVs", "snv_consequence_helpers.R"))

log_con <- file(file.path(snv_outdir, "filter_SNVs.log"), open = "wt")
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
# 2. Load combined SNV table
# ============================================================================

snv_path <- file.path(combined_dir, "combined_snvs.tsv")
snv_combined <- as.data.frame(fread(snv_path, sep = "\t", na.strings = c("", "NA")))

cat("\nLoaded", snv_path, "\n")
cat("Combined SNV table:", nrow(snv_combined), "variants from",
    n_distinct(snv_combined$ID), "patients,",
    n_distinct(snv_combined$sID), "samples\n")

# ============================================================================
# 3. Left-join project metadata (Disease, calling_mode, group1, group2 if set, QC)
# ============================================================================

meta_add <- metadata %>% select(sID, Disease, calling_mode, group1, any_of("group2"), QC)
snv_combined <- snv_combined %>%
  left_join(meta_add, by = "sID")

n_unjoined <- sum(is.na(snv_combined$Disease) & !snv_combined$sID %in% metadata$sID)
n_no_meta <- sum(!snv_combined$sID %in% metadata$sID)
if (n_no_meta > 0) {
  cat("Variants from samples not in metadata (dropped):", n_no_meta, "\n")
  cat("  Missing sIDs:", paste(unique(snv_combined$sID[!snv_combined$sID %in% metadata$sID]),
      collapse = ", "), "\n")
  snv_combined <- snv_combined %>% filter(sID %in% metadata$sID)
}

# Reorder: metadata first, then SNV columns in logical order
meta_cols <- intersect(c("group1", "group2", "Disease", "ID", "sID", "sample_type", "gender",
                         "capture_kit", "tumor_fraction", "calling_mode", "QC"),
                       names(snv_combined))
snv_col_order <- c(
  "Position", "SYMBOL", "Variant", "FILTER",
  "AF_tumor", "AF_ctrl", "GT_tumor", "GT_ctrl",
  "DP_tumor", "DP_ctrl",
  "AD_REF_tumor", "AD_REF_ctrl", "AD_ALT_tumor", "AD_ALT_ctrl",
  "normal_cn", "tumor_cn", "tumor_cn1", "tumor_cn2",
  "expected_mutant_copies", "CCF", "Consequence", "coding", "IMPACT",
  "HGVSc", "HGVSp"
)
remaining <- setdiff(names(snv_combined), c(meta_cols, snv_col_order))
snv_combined <- snv_combined[, c(
  meta_cols,
  intersect(snv_col_order, names(snv_combined)),
  remaining
)]

# ============================================================================
# 3a. Hotspots (snv_filter_params$hotspot_list)
# ============================================================================
# on_hotspot: the call sits on a hotspot of hotspot_list (gene, protein_position, ref_aa,
# alt_aa; empty alt_aa = any), matched on VEP SYMBOL, Protein_position and Amino_acids.
# Used by the strand-bias rescue (strandbias_hotspot_exempt) and the germline rescue.

on_hotspot <- rep(FALSE, nrow(snv_combined))
hotspot_path <- cfg$snv_filter_params$hotspot_list %||% ""
if (nzchar(hotspot_path)) {
  need <- c("Protein_position", "Amino_acids")
  if (!all(need %in% names(snv_combined))) {
    stop("hotspot_list needs the VEP columns ", paste(need, collapse = ", "), call. = FALSE)
  }
  hotspots <- read.csv(hotspot_path, stringsAsFactors = FALSE, colClasses = "character")
  hotspots[is.na(hotspots)] <- ""
  ref_aa <- sub("/.*", "", snv_combined$Amino_acids)
  alt_aa <- ifelse(grepl("/", snv_combined$Amino_acids, fixed = TRUE),
                   sub(".*/", "", snv_combined$Amino_acids), ref_aa)
  on_hotspot <- vapply(seq_len(nrow(snv_combined)), function(i) {
    h <- hotspots[hotspots$gene %in% snv_combined$SYMBOL[i], , drop = FALSE]
    isTRUE(any(h$protein_position %in% snv_combined$Protein_position[i] & h$ref_aa %in% ref_aa[i] &
                 (h$alt_aa == "" | h$alt_aa %in% alt_aa[i])))
  }, logical(1))
  rm(need, hotspots, ref_aa, alt_aa)
}

# ============================================================================
# 3b. Strand-bias rescue (snv_filter_params$strandbias_rescue)
# ============================================================================
# The pipeline's final VCF carries calls whose only Mutect2 label is strand_bias (FILTER ==
# "strand_bias"). Each is rescued to PASS by strandbias_rescue_mask() (strandbias_rescue.R), or keeps
# the label. Rescued calls go through the FLAGS and quality filters below like any PASS call.
# A call is recurrent when min_recurring_runs_for_strandbias_filter or more other patients (ID)
# carry it in this table or in strandbias_recurrence_reference (a combined_snvs.qs2 of other
# runs), counting non-excluded samples only. With strandbias_hotspot_exempt, calls on a
# hotspot_list hotspot are never recurrent: a hotspot recurs across patients by definition.

if (isTRUE(as.logical(cfg$snv_filter_params$strandbias_rescue))) {
  source(file.path(rootdir, "analysis-SNVs", "strandbias_rescue.R"))
  sb_p_min    <- as.numeric(cfg$snv_filter_params$strandbias_fisher_p_min)
  sb_min_alt  <- as.integer(cfg$snv_filter_params$strandbias_min_alt_reads)
  sb_af_min   <- as.numeric(cfg$snv_filter_params$strandbias_af_min)
  sb_min_runs <- as.integer(cfg$snv_filter_params$min_recurring_runs_for_strandbias_filter)
  if (is.na(sb_min_runs) || sb_min_runs < 1) {
    stop("min_recurring_runs_for_strandbias_filter must be an integer >= 1", call. = FALSE)
  }
  sb_cols <- c("SB_REF_FWD_tumor", "SB_REF_REV_tumor", "SB_ALT_FWD_tumor", "SB_ALT_REV_tumor")
  if (!all(sb_cols %in% names(snv_combined))) {
    stop("strandbias_rescue needs the ", paste(sb_cols, collapse = ", "),
         " columns; rerun WES-snakemake with keep_strand_bias_calls: true", call. = FALSE)
  }

  # Recurrence pool: non-excluded samples of this table plus the optional reference table
  pool <- snv_combined[!snv_combined$QC %in% "excluded", c("Position", "Variant", "ID")]
  ref_path <- cfg$snv_filter_params$strandbias_recurrence_reference %||% ""
  if (nzchar(ref_path)) {
    ref <- as.data.frame(qs2::qs_read(ref_path))
    ref <- ref[!ref$QC %in% "excluded", c("Position", "Variant", "ID")]
    pool <- rbind(pool, ref)
  }
  var_key <- paste(snv_combined$Position, snv_combined$Variant)
  n_other <- n_other_runs(var_key, snv_combined$ID, paste(pool$Position, pool$Variant), pool$ID)
  sb_hotspot_exempt <- isTRUE(as.logical(cfg$snv_filter_params$strandbias_hotspot_exempt))
  if (sb_hotspot_exempt) {
    if (!nzchar(hotspot_path)) stop("strandbias_hotspot_exempt needs snv_filter_params$hotspot_list", call. = FALSE)
    n_other[on_hotspot] <- 0L
  }

  is_sb <- snv_combined$FILTER == "strand_bias"
  rescued <- strandbias_rescue_mask(snv_combined$FILTER, as.matrix(snv_combined[, sb_cols]),
                                    snv_combined$AD_ALT_tumor, n_other,
                                    sb_p_min, sb_min_alt, sb_min_runs,
                                    af = snv_combined$AF_tumor, af_min = sb_af_min)
  snv_combined$FILTER[rescued] <- "PASS"

  cat("\n=== Strand-bias rescue ===\n")
  cat(sprintf("Fisher p >= %g, AD_ALT >= %d, AF >= %g, fewer than %d other run(s) with the call\n",
              sb_p_min, sb_min_alt, sb_af_min, sb_min_runs))
  cat(sprintf("strand_bias-only calls: %d, rescued to PASS: %d (on a hotspot, recurrence exempt: %d)\n",
              sum(is_sb), sum(rescued), if (sb_hotspot_exempt) sum(rescued & on_hotspot) else 0L))
  rm(pool, ref_path, var_key, n_other, is_sb, rescued, sb_hotspot_exempt)
} else if (any(snv_combined$FILTER == "strand_bias")) {
  # keep_strand_bias_calls was on but the rescue is off: strand_bias calls stay labelled
  cat("\nstrand_bias-only calls present, strandbias_rescue off:",
      sum(snv_combined$FILTER == "strand_bias"), "\n")
}

# ============================================================================
# 3c. Germline rescue (snv_filter_params$germline_rescue)
# ============================================================================
# The pipeline's final VCF carries, for tumor-only samples, calls whose only Mutect2 label is
# germline (not in the panel of normals, gnomAD AF <= tumor_only.germline_rescue_max_gnomad_af).
# Such a call becomes PASS when it is protein-altering (coding_consequence terms of
# SNV_consequence_terms.tsv) in a snv_genelist gene of snv_genelist_diagnosis, or sits on a
# hotspot of snv_filter_params$hotspot_list (gene, protein_position, ref_aa, alt_aa; empty
# alt_aa = any). Rescued calls skip the germline_gnomad filter below: the pipeline's cut already
# bounded their gnomAD AF, and VEP's gnomAD table can differ from the germline resource.

germline_rescued <- rep(FALSE, nrow(snv_combined))
if (isTRUE(as.logical(cfg$snv_filter_params$germline_rescue))) {
  if (!nzchar(hotspot_path) || !nzchar(snv_genelist)) {
    stop("germline_rescue needs snv_filter_params$hotspot_list and snv_genelist", call. = FALSE)
  }
  gl <- read.xlsx(snv_genelist)
  if (nzchar(snv_genelist_diagnosis)) gl <- gl %>% filter(diagnosis == snv_genelist_diagnosis)
  terms_tab <- read_consequence_terms(outdir)
  is_coding <- consequence_has_any(snv_combined$Consequence,
                                   terms_tab$term[terms_tab$coding_consequence])
  in_genelist <- snv_combined$SYMBOL %in% gl$genesymbol

  is_germline <- snv_combined$FILTER == "germline" & snv_combined$calling_mode %in% "tumor-only"
  germline_rescued <- is_germline & (((is_coding %in% TRUE) & in_genelist) | on_hotspot)
  snv_combined$FILTER[germline_rescued] <- "PASS"

  cat("\n=== Germline rescue ===\n")
  cat(sprintf("germline-only calls of tumor-only samples: %d, rescued to PASS: %d (gene list: %d, hotspot: %d)\n",
              sum(is_germline), sum(germline_rescued),
              sum(germline_rescued & is_coding & in_genelist), sum(germline_rescued & on_hotspot)))
  rm(gl, terms_tab, is_coding, in_genelist, is_germline)
} else if (any(snv_combined$FILTER == "germline")) {
  cat("\ngermline-only calls present, germline_rescue off:", sum(snv_combined$FILTER == "germline"), "\n")
}

# ============================================================================
# 3d. Clustered-call rescue (snv_filter_params$clustered_rescue)
# ============================================================================
# The pipeline keeps calls whose FilterMutectCalls labels are a subset of keep_filter_labels.
# Mutect2's clustered_events and haplotype labels also hit real clustered mutations: aSHM of
# PIM1/IGLL5/BCL2 and multi-hit drivers (KLHL14, EZH2). A call whose labels are only
# clustered_events and/or haplotype becomes PASS when its gene is in the snv_genelist or in
# clustered_rescue_genes (aSHM targets), not in clustered_rescue_exclude_genes, with AF_tumor >=
# clustered_rescue_af_min and AD_ALT_tumor >= clustered_rescue_min_alt_reads. Calls also labelled
# strand_bias or germline are not rescued. Rescued calls go through the filters below like any
# PASS call; non-coding aSHM calls end as low_priority_consequence.

if (isTRUE(as.logical(cfg$snv_filter_params$clustered_rescue))) {
  cr_genes_path <- cfg$snv_filter_params$clustered_rescue_genes %||% ""
  if (!nzchar(cr_genes_path) || !nzchar(snv_genelist)) {
    stop("clustered_rescue needs snv_filter_params$clustered_rescue_genes and snv_genelist", call. = FALSE)
  }
  cr_af_min  <- as.numeric(cfg$snv_filter_params$clustered_rescue_af_min)
  cr_min_alt <- as.integer(cfg$snv_filter_params$clustered_rescue_min_alt_reads)
  cr_exclude <- unlist(cfg$snv_filter_params$clustered_rescue_exclude_genes)
  gl <- read.xlsx(snv_genelist)
  if (nzchar(snv_genelist_diagnosis)) gl <- gl %>% filter(diagnosis == snv_genelist_diagnosis)
  ashm <- read.csv(cr_genes_path, stringsAsFactors = FALSE)$gene
  cr_genes <- setdiff(union(gl$genesymbol, ashm), cr_exclude)

  is_clustered <- grepl("^(clustered_events|haplotype)(;(clustered_events|haplotype))*$",
                        snv_combined$FILTER)
  clustered_rescued <- is_clustered & snv_combined$SYMBOL %in% cr_genes &
    !is.na(snv_combined$AF_tumor) & snv_combined$AF_tumor >= cr_af_min &
    !is.na(snv_combined$AD_ALT_tumor) & snv_combined$AD_ALT_tumor >= cr_min_alt
  snv_combined$FILTER[clustered_rescued] <- "PASS"

  cat("\n=== Clustered-call rescue ===\n")
  cat(sprintf("gene list + %d aSHM genes (minus %s), AF >= %g, AD_ALT >= %d\n", length(ashm),
              paste(cr_exclude, collapse = "/"), cr_af_min, cr_min_alt))
  cat(sprintf("clustered_events/haplotype-only calls: %d, rescued to PASS: %d (aSHM genes: %d)\n",
              sum(is_clustered), sum(clustered_rescued),
              sum(clustered_rescued & snv_combined$SYMBOL %in% ashm)))
  rm(cr_genes_path, cr_af_min, cr_min_alt, cr_exclude, gl, ashm, cr_genes, is_clustered,
     clustered_rescued)
}

# ============================================================================
# 4. Apply FLAGS filter
# ============================================================================

flags_genelist <- if (exists("cfg") && !is.null(cfg$flags_genelist) && nzchar(cfg$flags_genelist)) {
  cfg$flags_genelist
} else {
  "flags_genes.csv"
}
flags <- read.csv(file.path(outdir, flags_genelist), stringsAsFactors = FALSE)

exact_flags <- flags %>% filter(match_type == "exact")
prefix_flags <- flags %>% filter(match_type == "prefix")

gene_col <- snv_combined$SYMBOL

# Exact matches
exact_match <- gene_col %in% exact_flags$gene
exact_reasons <- exact_flags$reason[match(gene_col, exact_flags$gene)]

# Prefix matches
prefix_match <- rep(FALSE, length(gene_col))
prefix_reasons <- rep(NA_character_, length(gene_col))
for (i in seq_len(nrow(prefix_flags))) {
  hit <- grepl(paste0("^", prefix_flags$gene[i], "[0-9]"), gene_col)
  # Don't overwrite if already matched by a previous prefix
  new_hit <- hit & !prefix_match
  prefix_match[new_hit] <- TRUE
  prefix_reasons[new_hit] <- prefix_flags$reason[i]
}

# Combine: exact takes priority
is_flagged <- exact_match | prefix_match
reasons <- ifelse(exact_match, exact_reasons, prefix_reasons)

# Update FILTER column
snv_combined$FILTER <- ifelse(
  is_flagged & snv_combined$FILTER == "PASS",
  paste0("FLAGS_", reasons),
  snv_combined$FILTER
)

# High-AF flag: PASS calls in high_af_flag_genes above high_af_flag_min (HLA class I calls at
# AF ~1 are mismapped reads of other HLA alleles)
haf_genes <- unlist(cfg$snv_filter_params$high_af_flag_genes)
if (length(haf_genes)) {
  haf_min <- as.numeric(cfg$snv_filter_params$high_af_flag_min)
  is_high_af <- snv_combined$FILTER == "PASS" & snv_combined$SYMBOL %in% haf_genes &
    !is.na(snv_combined$AF_tumor) & snv_combined$AF_tumor > haf_min
  snv_combined$FILTER[is_high_af] <- "FLAGS_high_AF"
  cat(sprintf("\nFLAGS_high_AF (%s, AF > %g): %d\n", paste(haf_genes, collapse = "/"), haf_min,
              sum(is_high_af)))
  rm(haf_min, is_high_af)
}
rm(haf_genes)

n_pass <- sum(snv_combined$FILTER == "PASS")

cat("\n=== FLAGS filter ===\n")
cat("Variants flagged:", sum(grepl("^FLAGS_", snv_combined$FILTER)), "\n")
cat("Variants remaining PASS:", n_pass, "\n")
cat("FILTER breakdown:\n")
print(table(snv_combined$FILTER))

flagged_genes <- snv_combined %>%
  filter(grepl("^FLAGS_", FILTER)) %>%
  count(SYMBOL, FILTER) %>%
  arrange(desc(n))
cat("\nFlagged genes (top 20):\n")
print(head(flagged_genes, 20))

# ============================================================================
# 5. Apply quality filters
# ============================================================================
# Each filter reclassifies PASS variants to a named label; all five are always applied.
# Non-PASS variants (FLAGS_*) are never touched.
# Priority order: low_alt_reads > low_AF > tumor_only_low_AF > germline_ctrl > germline_gnomad.
# germline_gnomad sits last on purpose: it is the backstop for variants that pass the
# germline_ctrl ratio test but are common in the population.

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && nzchar(a)) a else b
}

# germline_ctrl thresholds (config.yaml snv_filter_params block). When the matched
# "normal" can carry tumor (e.g. PBMC/bone marrow in a leukemic lymphoma), a raised AF_ctrl
# alone does not mean germline — the ratio test separates a germline heterozygote (ratio
# near 1) from tumor-in-normal contamination (ratio well below).
germline_ctrl_af_min    <- as.numeric(cfg$snv_filter_params$germline_ctrl_af_min)
germline_ctrl_ratio_min <- as.numeric(cfg$snv_filter_params$germline_ctrl_ratio_min)

cat("\n=== Quality filters ===\n")
cat("Applied:", paste(quality_filter_labels, collapse = ", "), "\n")
cat(sprintf("germline_ctrl: AF_ctrl > %g and AF_ctrl/AF_tumor > %g\n\n",
            germline_ctrl_af_min, germline_ctrl_ratio_min))

n_pass_before_qf <- sum(snv_combined$FILTER == "PASS")

# VEP can, in principle, emit multiple gnomAD AFs at a multi-allelic site joined by
# "," "|" or "&" depending on annotation options. as.numeric() would silently return NA
# for these, coercing them to 0 and letting germline variants escape the filter.
# Take the max across delimited tokens to correctly classify such sites.
.parse_gnomad_af <- function(x) {
  if (is.numeric(x)) return(ifelse(is.na(x), 0, x))
  vapply(x, function(v) {
    if (is.na(v) || v == "" || v == ".") return(0)
    vals <- suppressWarnings(as.numeric(strsplit(v, "[|,&]")[[1]]))
    max(c(vals[!is.na(vals)], 0))
  }, numeric(1))
}

# AF_ctrl/AF_tumor is undefined when AF_tumor is missing or non-positive; those variants
# are classified on the absolute AF_ctrl test alone.
.ctrl_ratio_germline <- function(af_ctrl, af_tumor) {
  ifelse(is.na(af_tumor) | af_tumor <= 0,
         TRUE,
         af_ctrl / af_tumor > germline_ctrl_ratio_min)
}

snv_combined <- snv_combined %>%
  mutate(FILTER = case_when(
    FILTER != "PASS"                                                    ~ FILTER,
    AD_ALT_tumor <= 3                                                   ~ "low_alt_reads",
    AF_tumor < 0.05                                                     ~ "low_AF",
    is.na(AF_ctrl) & AF_tumor < 0.2                                     ~ "tumor_only_low_AF",
    !is.na(AF_ctrl) & AF_ctrl > germline_ctrl_af_min &
      .ctrl_ratio_germline(AF_ctrl, AF_tumor)                           ~ "germline_ctrl",
    !germline_rescued &
      pmax(.parse_gnomad_af(gnomADe_AF),
           .parse_gnomad_af(gnomADg_AF)) > 0.05                         ~ "germline_gnomad",
    TRUE                                                                ~ "PASS"
  ))

for (lbl in quality_filter_labels) {
  cat(sprintf("  %-22s : %5d variants removed\n", lbl, sum(snv_combined$FILTER == lbl)))
}
n_pass_after_qf <- sum(snv_combined$FILTER == "PASS")
cat(sprintf("\nVariants PASS after quality filters: %d  (removed %d of %d)\n",
            n_pass_after_qf, n_pass_before_qf - n_pass_after_qf, n_pass_before_qf))
rm(.parse_gnomad_af, .ctrl_ratio_germline, germline_rescued)

# ============================================================================
# 6. Update metadata QC column from per-sample PASS SNV count
# ============================================================================
# Samples with an implausibly high PASS SNV count are flagged in metadata.csv. QC is
# "excluded" (set by load_WES_results.R from sample_annotations), else "high_SNV_count",
# else "PASS". The derived values replace the ones joined in section 3, so the saved table
# and the rewritten metadata.csv agree. Per-sample counts and the count-only verdict
# (snv_qc, independent of exclusion) go to analysis-SNVs/SNV_counts.csv.

.snv_qc_threshold <- 2000L

.snv_counts <- metadata %>%
  select(ID, sID) %>%
  left_join(snv_combined %>%
              group_by(sID) %>%
              summarise(n_snv = n(), n_pass_snv = sum(FILTER == "PASS"), .groups = "drop"),
            by = "sID") %>%
  mutate(
    n_snv      = coalesce(n_snv, 0L),
    n_pass_snv = coalesce(n_pass_snv, 0L),
    snv_qc     = ifelse(n_pass_snv > .snv_qc_threshold, "high_SNV_count", "PASS")
  )

write.csv(.snv_counts, file.path(snv_outdir, "SNV_counts.csv"), row.names = FALSE, quote = FALSE)

.meta_path <- file.path(metadata_dir, "metadata.csv")
.meta_qc <- read.csv(.meta_path, stringsAsFactors = FALSE, check.names = FALSE) %>%
  left_join(.snv_counts %>% select(sID, snv_qc), by = "sID") %>%
  mutate(QC = case_when(QC %in% "excluded"           ~ "excluded",
                        snv_qc %in% "high_SNV_count" ~ "high_SNV_count",
                        TRUE                         ~ "PASS")) %>%
  select(-snv_qc)

write.csv(.meta_qc, .meta_path, row.names = FALSE)

snv_combined$QC <- .meta_qc$QC[match(snv_combined$sID, .meta_qc$sID)]
metadata$QC     <- .meta_qc$QC[match(metadata$sID, .meta_qc$sID)]

cat("\n=== Metadata QC ===\n")
cat(sprintf("  QC threshold: >%d PASS SNVs\n", .snv_qc_threshold))
cat(sprintf("  %d PASS, %d high_SNV_count, %d excluded\n", sum(.meta_qc$QC == "PASS"),
            sum(.meta_qc$QC == "high_SNV_count"), sum(.meta_qc$QC == "excluded")))
for (.lbl in c("high_SNV_count", "excluded")) {
  if (any(.meta_qc$QC == .lbl)) {
    cat(sprintf("  %s: %s\n", .lbl, paste(.meta_qc$sID[.meta_qc$QC == .lbl], collapse = ", ")))
  }
}
rm(.snv_qc_threshold, .snv_counts, .meta_path, .meta_qc, .lbl)

# ============================================================================
# 7. Apply consequence allowlist filter
# ============================================================================
# Reclassifies (never drops) currently-PASS variants whose Consequence has no term
# flagged pass_consequence == TRUE in SNV_consequence_terms.tsv. VEP Consequence
# can be compound ("&"-joined SO terms co-occurring on the same transcript, e.g.
# splice_acceptor_variant&coding_sequence_variant&intron_variant), so matching is done
# per-term after splitting on "&", not by a substring/regex test on the whole string —
# a substring exclusion on e.g. "intron_variant" would incorrectly catch splice
# acceptor/donor calls, which are always co-tagged intronic.
# Run after the metadata QC step (section 6) so the >2000 PASS-SNV threshold there stays
# calibrated against the pre-consequence-filter PASS set.

consequence_filter_label <- "low_priority_consequence"

.consequence_term_table <- read_consequence_terms(outdir)
consequence_allowlist <- .consequence_term_table$term[.consequence_term_table$pass_consequence]

.consequence_terms <- strsplit(snv_combined$Consequence, "&", fixed = TRUE)
.has_allowed_term <- vapply(.consequence_terms,
                             function(x) any(x %in% consequence_allowlist),
                             logical(1))

n_pass_before_cf <- sum(snv_combined$FILTER == "PASS")

snv_combined$FILTER <- ifelse(
  snv_combined$FILTER == "PASS" & !.has_allowed_term,
  consequence_filter_label,
  snv_combined$FILTER
)

n_pass_after_cf <- sum(snv_combined$FILTER == "PASS")

cat("\n=== Consequence allowlist filter ===\n")
cat("Allowlist (", file.path(outdir, "SNV_consequence_terms.tsv"), ", pass_consequence == TRUE):\n", sep = "")
cat(paste0("  ", consequence_allowlist), sep = "\n")
cat(sprintf("\nVariants PASS after consequence filter: %d  (reclassified %d of %d)\n",
            n_pass_after_cf, n_pass_before_cf - n_pass_after_cf, n_pass_before_cf))
# coding: TRUE if any Consequence term has coding_consequence == TRUE, independent of FILTER
.coding_terms <- .consequence_term_table$term[.consequence_term_table$coding_consequence]
snv_combined$coding <- consequence_has_any(snv_combined$Consequence, .coding_terms)
snv_combined <- snv_combined %>% relocate(coding, .after = Consequence)
cat(sprintf("coding: %d variants (%d PASS)\n",
            sum(snv_combined$coding), sum(snv_combined$coding & snv_combined$FILTER == "PASS")))

rm(.consequence_terms, .has_allowed_term, .consequence_term_table, .coding_terms)

# ============================================================================
# 8. Save combined table
# ============================================================================

qs_save(snv_combined, file.path(snv_outdir, "combined_snvs.qs2"))
cat("\nSaved combined_snvs.qs2:", nrow(snv_combined), "rows x",
    ncol(snv_combined), "columns\n")
cat(sprintf("  PASS: %d  |  FLAGS: %d  |  Quality-filtered: %d  |  %s: %d\n",
    n_pass_after_cf,
    sum(grepl("^FLAGS_", snv_combined$FILTER)),
    sum(snv_combined$FILTER %in% quality_filter_labels),
    consequence_filter_label,
    sum(snv_combined$FILTER == consequence_filter_label)))

# ============================================================================
# 9. Export genelist subset as Excel (optional — requires snv_genelist in config)
# ============================================================================

if (!exists("snv_genelist")) snv_genelist <- ""
if (!exists("snv_genelist_diagnosis")) snv_genelist_diagnosis <- ""

if (nzchar(snv_genelist) && file.exists(snv_genelist)) {
  genelist <- read.xlsx(snv_genelist)
  if (nzchar(snv_genelist_diagnosis)) {
    genelist <- genelist %>% filter(diagnosis == snv_genelist_diagnosis)
  }
  gl_genes <- unique(genelist$genesymbol)
  cat("\nGene list:", length(gl_genes), "genes",
      if (nzchar(snv_genelist_diagnosis)) paste0("(", snv_genelist_diagnosis, ")") else "",
      "\n")

  # Excel caps a cell at 32,767 characters. VEP DOMAINS lists every PDB chain mapped to
  # the protein (~2,000 for B2M), so the export replaces those entries with their count
  # and truncates anything still over the cap. combined_snvs.qs2 keeps the full column.
  collapse_pdb <- function(domains) {
    vapply(domains, function(d) {
      if (is.na(d)) return(NA_character_)
      parts <- strsplit(d, "&", fixed = TRUE)[[1]]
      is_pdb <- startsWith(parts, "PDB-ENSP_mappings:")
      if (sum(is_pdb) > 1) {
        parts <- c(parts[!is_pdb], paste0("PDB-ENSP_mappings:n=", sum(is_pdb)))
      }
      substr(paste(parts, collapse = "&"), 1, 32000)
    }, character(1), USE.NAMES = FALSE)
  }
  snv_gl <- snv_combined %>%
    filter(SYMBOL %in% gl_genes) %>%
    mutate(DOMAINS = collapse_pdb(DOMAINS))
  write.xlsx(snv_gl, file.path(snv_outdir, "combined_snvs_genelisted.xlsx"), overwrite = TRUE)
  cat("Saved combined_snvs_genelisted.xlsx:", nrow(snv_gl), "rows\n")
} else {
  cat("\nNo SNV genelist configured — skipping genelist export.\n")
}

cat("\nDone!\n")

w <- warnings()
if (length(w) > 0) {
  cat("\n=== Warnings ===\n")
  print(w)
}

sink(type = "message")
sink()
close(log_con)
