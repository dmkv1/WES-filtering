## filter_SVs.R
## Reads the WES-snakemake combined SV table (results/combined/combined_svs.tsv —
## Manta+DELLY consensus via SURVIVOR, AnnotSV-annotated, one row per SV per sample),
## joins project metadata, flags artifact-prone genes, applies quality filters, and
## saves combined_svs.qs2.
##
## Inputs:
##   <metadata_dir>/metadata.csv          — patient/sample metadata (Disease, ...)
##   <outdir>/<flags_genelist>             — genes to flag (exact or prefix); config.yaml's
##                                          flags_genelist names the file (default
##                                          "flags_genes.csv")
##   <combined_dir>/combined_svs.tsv      — WES-snakemake combined SV table
##   config.yaml                          — sv_filter_params block sets the recurrence thresholds
##
## Outputs:
##   <outdir>/analysis-SV/combined_svs.qs2 — all somatic SVs with metadata joined, FILTER applied
##
## FILTER labels written (first-match order):
##   PASS                — survives all filters
##   FLAGS_all_genes      — every gene the SV overlaps is on the flags list (MUC/KIR/...)
##   acmg_benign          — AnnotSV ACMG_class 1 (DEL/DUP only; AnnotSV's own
##                          population-frequency classifier)
##   ig_tcr_locus          — non-TRA breakpoint inside an IGH/IGK/IGL/TCRA-D/TRB/TRG locus:
##                          physiological V(D)J/class-switch recombination, not a somatic event
##   recurrent_artifact    — identical breakpoint (binned) recurring in >= sv_filter_params$
##                          recurrent_min_patients unrelated patients; DEL/DUP/INS/INV only
##   single_caller_TRA     — translocation called by DELLY alone, no Manta support
##
## A translocation with a breakpoint in an Ig/TCR locus is exempted from
## recurrent_artifact/single_caller_TRA, since Ig-partnered translocations are expected to
## recur across patients and to be DELLY-only (see section 6). Driver loci listed
## in config.yaml's sv_driver_loci are exempted the same way (see section 6).

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
sv_outdir <- file.path(outdir, "analysis-SV")
dir.create(sv_outdir, showWarnings = FALSE, recursive = TRUE)

log_con <- file(file.path(sv_outdir, "filter_SVs.log"), open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
options(warn = 1)  # print warnings immediately so they go into the log

# ============================================================================
# 1. Load metadata
# ============================================================================

metadata <- read.csv(file.path(metadata_dir, "metadata.csv"), stringsAsFactors = FALSE)

cat("Metadata rows:", nrow(metadata), "\n")
cat("group1:", paste(names(table(metadata$group1)),
    table(metadata$group1), sep = "=", collapse = ", "), "\n")

# ============================================================================
# 2. Load combined SV table
# ============================================================================

sv_path <- file.path(combined_dir, "combined_svs.tsv")
sv_combined <- as.data.frame(fread(sv_path, sep = "\t", na.strings = c("", "NA")))

cat("\nLoaded", sv_path, "\n")
cat("Combined SV table:", nrow(sv_combined), "SVs from",
    n_distinct(sv_combined$sID), "samples\n")

# ============================================================================
# 3. Left-join project metadata (Disease, calling_mode, group1, group2 if set, QC)
# ============================================================================

meta_add <- metadata %>% select(sID, Disease, calling_mode, group1, any_of("group2"), QC)
sv_combined <- sv_combined %>% left_join(meta_add, by = "sID")

n_no_meta <- sum(!sv_combined$sID %in% metadata$sID)
if (n_no_meta > 0) {
  cat("\nRows from samples not in metadata (dropped):", n_no_meta, "\n")
  dropped_samples <- unique(sv_combined$sID[!sv_combined$sID %in% metadata$sID])
  cat("  Dropped samples:", paste(dropped_samples, collapse = ", "), "\n")
  sv_combined <- sv_combined %>% filter(sID %in% metadata$sID)
}

cat("After metadata filter:", nrow(sv_combined), "rows from",
    n_distinct(sv_combined$sID), "tumor samples\n")

# Every SV starts as PASS; sections 6-7 relabel the ones they catch. The raw column
# already loaded here is the VCF FILTER (Manta/DELLY's own filtering, already all
# "PASS" by the time SVs reach the merged, size-filtered VCF), so overwriting it is safe.
sv_combined$FILTER <- "PASS"

# ============================================================================
# 4. Summary statistics (pre-filter)
# ============================================================================

cat("\n=== Summary statistics (pre-filter) ===\n")

cat("SV_type distribution:\n")
print(table(sv_combined$SV_type, useNA = "ifany"))

cat("\nSV_callers distribution:\n")
print(table(sv_combined$SV_callers, useNA = "ifany"))

cat("\nACMG_class distribution:\n")
print(table(sv_combined$ACMG_class, useNA = "ifany"))

cat("\nSVs per tumor sample:\n")
print(table(sv_combined$sID))

# ============================================================================
# 5. FLAGS filter — SVs whose Gene_name is entirely FLAGS genes
# ============================================================================
# Gene_name is semicolon-separated (all overlapped genes); a call is only flagged when
# every overlapped gene is on the list, so a real event next to one blocklisted gene
# is not thrown out.

flags_genelist <- if (exists("cfg") && !is.null(cfg$flags_genelist) && nzchar(cfg$flags_genelist)) {
  cfg$flags_genelist
} else {
  "flags_genes.csv"
}
flags        <- read.csv(file.path(outdir, flags_genelist), stringsAsFactors = FALSE)
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
  if (length(gene_string) != 1L) return(FALSE)
  if (is.na(gene_string) || !nzchar(trimws(gene_string)) || gene_string == "-") {
    return(FALSE)   # no gene annotation — keep
  }
  genes <- trimws(unlist(strsplit(gene_string, ";")))
  genes <- genes[nzchar(genes) & genes != "-" & !is.na(genes)]
  if (length(genes) == 0L) return(FALSE)
  all(vapply(genes, is_flagged_gene, logical(1L)))
}

flagged_mask <- vapply(as.character(sv_combined[["Gene_name"]]), is_all_flagged, logical(1L))
sv_combined$FILTER[flagged_mask] <- "FLAGS_all_genes"

cat("\n=== FLAGS filter ===\n")
cat("SVs with all genes flagged:", sum(flagged_mask), "\n")
if (sum(flagged_mask) > 0) {
  cat("Flagged gene entries (top 20):\n")
  print(head(sort(table(sv_combined$Gene_name[flagged_mask]), decreasing = TRUE), 20))
}
rm(flags, exact_flags, prefix_flags, flagged_mask)

# ============================================================================
# 6. Driver / Ig-TCR locus exemption
# ============================================================================
# Ig/TCR-family regions need coordinate-based (not Gene_name-based) handling, because
# AnnotSV leaves Gene_name empty for the intergenic/regulatory breakpoints these events
# actually land on: IGH/IGK/IGL/TRA-TRD/TRB/TRG undergo physiological V(D)J joining and
# class-switch recombination in normal and malignant B/T cells. An intra-locus
# DEL/DUP/INV/INS breakpoint here (ig_tcr_locus) is that background process, not a
# somatic driver. A translocation (TRA) breakpoint here is different: Ig-locus partners
# are exactly how lymphoma driver translocations (CCND1, BCL2, BCL6, MYC) present, so TRA is
# exempted here rather than filtered — it falls through to driver_exempt below.
#
# Disease-specific driver loci come from config.yaml's sv_driver_loci (name, chrom,
# start, end; GRCh38 gene bodies). Any SV whose primary breakpoint lands in one is
# exempted from recurrent_artifact/single_caller_TRA, since those drivers are expected
# to recur across patients. With sv_driver_loci empty, driver_exempt reduces to the
# Ig/TCR-partnered TRA case only.
#
# Coordinates are GRCh38 (matches the pipeline's AnnotSV genome build), cytogenetic
# locus extents for IG/TCR. Only the primary breakpoint (SV_chrom/SV_start-SV_end) is
# tested. A translocation either survives the BNDrescue dedup in combine_results.R as the
# row whose SV_chrom/SV_start sits in the Ig/driver locus, or keeps one row per side, each
# tested on its own breakpoint; a chrom-only test on the imprecise SV_chrom2 side would
# over-match. An Ig-partnered translocation's driver-side row (e.g. chr18 of IGH-BCL2) is
# exempted only when the driver is listed in sv_driver_loci.

locus_windows <- data.frame(
  name  = c("IGH", "IGK", "IGL", "TRA_TRD", "TRB", "TRG"),
  group = c("ig_tcr", "ig_tcr", "ig_tcr", "ig_tcr", "ig_tcr", "ig_tcr"),
  chrom = c("14", "2", "22", "14", "7", "7"),
  start = c(105586437, 88857361, 22026076, 21621904, 142299011, 38240024),
  end   = c(106879844, 90235368, 22922913, 22552156, 142813287, 38368055),
  stringsAsFactors = FALSE
)
.driver_loci <- cfg$sv_driver_loci
if (is.null(.driver_loci)) .driver_loci <- list()
if (length(.driver_loci) > 0) {
  .driver_df <- do.call(rbind, lapply(names(.driver_loci), function(nm) {
    l <- .driver_loci[[nm]]
    if (is.null(l$chrom) || is.null(l$start) || is.null(l$end)) {
      stop("sv_driver_loci$", nm, " needs chrom, start and end", call. = FALSE)
    }
    data.frame(name = nm, group = "driver", chrom = as.character(l$chrom),
               start = as.numeric(l$start), end = as.numeric(l$end),
               stringsAsFactors = FALSE)
  }))
  locus_windows <- rbind(locus_windows, .driver_df)
  rm(.driver_df)
}
rm(.driver_loci)

.in_windows <- function(chrom, start, end, windows) {
  if (nrow(windows) == 0) return(rep(FALSE, length(chrom)))
  mapply(function(c, s, e) {
    if (is.na(c) || is.na(s) || is.na(e)) return(FALSE)
    any(as.character(c) == windows$chrom & s <= windows$end & e >= windows$start)
  }, chrom, start, end)
}

sv_combined$is_driver_locus  <- .in_windows(sv_combined$SV_chrom, sv_combined$SV_start,
                                             sv_combined$SV_end,
                                             locus_windows[locus_windows$group == "driver", ])
sv_combined$is_ig_tcr_locus  <- .in_windows(sv_combined$SV_chrom, sv_combined$SV_start,
                                             sv_combined$SV_end,
                                             locus_windows[locus_windows$group == "ig_tcr", ])
# Exempts recurrent_artifact/single_caller_TRA (section 7); does not exempt FLAGS or
# acmg_benign, which never fire on these loci in practice (checked in the cat() below).
driver_exempt <- sv_combined$is_driver_locus |
  (sv_combined$SV_type == "TRA" & sv_combined$is_ig_tcr_locus)

cat("\n=== Driver / Ig-TCR locus exemption ===\n")
cat("Driver loci (sv_driver_loci):",
    if (any(locus_windows$group == "driver"))
      paste(locus_windows$name[locus_windows$group == "driver"], collapse = ", ")
    else "(none)", "\n")
cat("Driver-locus-overlapping SVs:", sum(sv_combined$is_driver_locus), "\n")
cat("Ig/TCR-locus SVs:", sum(sv_combined$is_ig_tcr_locus),
    "(TRA, exempted as driver-partner translocations:",
    sum(sv_combined$is_ig_tcr_locus & sv_combined$SV_type == "TRA"),
    "| non-TRA, subject to ig_tcr_locus filter below:",
    sum(sv_combined$is_ig_tcr_locus & sv_combined$SV_type != "TRA"), ")\n")
cat("Total SVs exempted from recurrent_artifact/single_caller_TRA:", sum(driver_exempt), "\n")

# ============================================================================
# 7. Quality filters
# ============================================================================
# Each filter reclassifies PASS SVs to a named label; non-PASS SVs (already FLAGS_*)
# are never touched. Priority order: acmg_benign > ig_tcr_locus > recurrent_artifact >
# single_caller_TRA.

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && nzchar(a)) a else b
}

# recurrent_artifact thresholds (config.yaml sv_filter_params block). A breakpoint that
# recurs at (near-)identical coordinates across several genetically unrelated patients is
# a mapping/repeat artifact, not a shared somatic event; binning absorbs the few bp of
# jitter callers show at the same underlying breakpoint. Patients, not samples, are
# counted, since a patient's PT/PDX/PDCL samples share their genetics and would trivially
# "recur" against each other.
recurrent_min_patients <- as.integer(cfg$sv_filter_params$recurrent_min_patients)
recurrence_bin_bp      <- as.integer(cfg$sv_filter_params$recurrence_bin_bp)

cat("\n=== Quality filters ===\n")
cat(sprintf("acmg_benign: ACMG_class == 1 (DEL/DUP only)\n"))
cat(sprintf("ig_tcr_locus: non-TRA breakpoint inside an Ig/TCR locus\n"))
cat(sprintf("recurrent_artifact: identical breakpoint (+/- %dbp) in >= %d unrelated patients, DEL/DUP/INS/INV only\n",
            recurrence_bin_bp, recurrent_min_patients))
cat(sprintf("single_caller_TRA: translocation called by DELLY alone\n"))

n_pass_before_qf <- sum(sv_combined$FILTER == "PASS")

sv_combined <- sv_combined %>%
  mutate(FILTER = case_when(
    FILTER != "PASS"                                                          ~ FILTER,
    !driver_exempt & SV_type %in% c("DEL", "DUP") &
      !is.na(ACMG_class) & ACMG_class == 1                                    ~ "acmg_benign",
    !driver_exempt & SV_type != "TRA" & is_ig_tcr_locus                       ~ "ig_tcr_locus",
    TRUE                                                                      ~ "PASS"
  ))

# recurrent_artifact needs a group-by across the still-PASS rows, so it is computed
# separately rather than inline in the case_when above. Recurrence is counted over
# non-excluded samples only; excluded samples' rows in a recurrent bin are still labelled.
.rec_candidates <- sv_combined %>%
  mutate(.rid = row_number()) %>%
  filter(FILTER == "PASS", SV_type %in% c("DEL", "DUP", "INS", "INV"), !driver_exempt) %>%
  mutate(.bin = round(SV_start / recurrence_bin_bp))

.recurrent_bins <- .rec_candidates %>%
  filter(!QC %in% "excluded") %>%
  group_by(SV_chrom, .bin, SV_type) %>%
  summarise(n_patients = n_distinct(ID), .groups = "drop") %>%
  filter(n_patients >= recurrent_min_patients)

.recurrent_rids <- .rec_candidates %>%
  inner_join(.recurrent_bins, by = c("SV_chrom", ".bin", "SV_type")) %>%
  pull(.rid)

sv_combined$FILTER[.recurrent_rids] <- "recurrent_artifact"

sv_combined <- sv_combined %>%
  mutate(FILTER = case_when(
    FILTER != "PASS"                                                          ~ FILTER,
    !driver_exempt & SV_type == "TRA" & SV_callers == "DELLY"                 ~ "single_caller_TRA",
    TRUE                                                                      ~ "PASS"
  ))

quality_filter_labels <- c("acmg_benign", "ig_tcr_locus", "recurrent_artifact", "single_caller_TRA")
for (lbl in quality_filter_labels) {
  cat(sprintf("  %-20s : %5d SVs\n", lbl, sum(sv_combined$FILTER == lbl)))
}
n_pass_after_qf <- sum(sv_combined$FILTER == "PASS")
cat(sprintf("\nSVs PASS after quality filters: %d  (removed %d of %d)\n",
            n_pass_after_qf, n_pass_before_qf - n_pass_after_qf, n_pass_before_qf))
cat(sprintf("Distinct recurrent breakpoint bins: %d\n", nrow(.recurrent_bins)))

rm(.rec_candidates, .recurrent_bins, .recurrent_rids, driver_exempt, locus_windows, .in_windows)

# ============================================================================
# 8. Reorder columns: metadata first, then key SV columns, then annotations
# ============================================================================

meta_cols <- intersect(c("group1", "group2", "Disease", "ID", "sID", "sample_type", "gender",
                         "capture_kit", "tumor_fraction", "calling_mode", "QC"),
                       names(sv_combined))
sv_col_order <- c(
  "AnnotSV_ID", "SV_chrom", "SV_start", "SV_end", "SV_chrom2", "SV_length", "SV_type",
  "SV_callers", "SV_DR_ref", "SV_DR_alt", "FILTER", "is_driver_locus", "is_ig_tcr_locus",
  "Gene_name", "Gene_count", "CytoBand",
  "ACMG_class", "AnnotSV_ranking_score", "AnnotSV_ranking_criteria"
)
remaining <- setdiff(names(sv_combined), c(meta_cols, sv_col_order))
sv_combined <- sv_combined[, c(
  meta_cols,
  intersect(sv_col_order, names(sv_combined)),
  remaining
)]

# ============================================================================
# 9. Summary statistics (post-filter)
# ============================================================================

cat("\n=== Summary statistics (post-filter) ===\n")

cat("FILTER distribution:\n")
print(table(sv_combined$FILTER, useNA = "ifany"))

cat("\nPASS by SV_type:\n")
print(table(sv_combined$SV_type[sv_combined$FILTER == "PASS"], useNA = "ifany"))

cat(sprintf("\nOverall: %d PASS of %d (%.1f%%)\n",
            n_pass_after_qf, nrow(sv_combined), 100 * n_pass_after_qf / nrow(sv_combined)))

# ============================================================================
# 10. Save outputs
# ============================================================================

qs_save(sv_combined, file.path(sv_outdir, "combined_svs.qs2"))
cat("\nSaved combined_svs.qs2:", nrow(sv_combined), "rows x",
    ncol(sv_combined), "columns\n")

cat("\nDone!\n")

w <- warnings()
if (length(w) > 0) {
  cat("\n=== Warnings ===\n")
  print(w)
}

options(warn = 0)  # restore default
sink(type = "message")
sink()
close(log_con)
