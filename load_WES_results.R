## load_WES_results.R
## Entry-point script: reads config.yaml, refreshes metadata.csv from the
## WES-snakemake sample sheet, then sources filter_SNVs.R, filter_CNVs.R
## and filter_SVs.R, which read the pipeline's combined_*.tsv tables directly.
##
## filter_SNVs.R rewrites metadata.csv with the QC column derived from per-sample
## PASS SNV counts, so filter_CNVs.R and filter_SVs.R, sourced afterwards, read the
## same QC values that filter_SNVs.R stored in combined_snvs.qs2.
##
## metadata/metadata.ver holds the name of the experiment directory that last wrote
## metadata.csv.
##
## config.yaml lives with the experiment it configures (results/<date>_WES-filtering/config.yaml),
## not alongside this script, since it pins that experiment's dated paths and tuned thresholds
## — so --config is required rather than defaulting to a path next to this script. Every
## path and project-specific choice (results_dir, metadata_dir, group1/group2) comes
## from the config, and outputs go to the directory holding it; the filter scripts carry none of their own.
## Keys the config leaves out take the defaults of template/config.yaml (config.R).
##
## Usage:
##   Rscript load_WES_results.R --config results/<date>_WES-filtering/config.yaml

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required. Install with: install.packages('yaml')", call. = FALSE)
}

# Auto-detect rootdir from this script's location
.script_path <- if (sys.nframe() > 0 && !is.null(sys.frame(1)$ofile)) {
  sys.frame(1)$ofile
} else {
  sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
}
rootdir <- normalizePath(dirname(.script_path))
rm(.script_path)

source(file.path(rootdir, "config.R"))
cfg <- read_wes_config(config_arg("load_WES_results.R"))
config_path <- cfg$config_path

# Required: results_dir (WES-snakemake results/ directory)
results_dir <- cfg$results_dir
if (is.null(results_dir) || results_dir == "") {
  stop("results_dir must be set in ", config_path, call. = FALSE)
}
results_dir <- normalizePath(results_dir, mustWork = TRUE)
combined_dir <- file.path(results_dir, "combined")
if (!dir.exists(combined_dir)) {
  stop("results_dir does not contain a combined/ subdirectory: ", combined_dir, "\n",
       "  Expected WES-snakemake >=1.0.0 output (combined_snvs.tsv etc).", call. = FALSE)
}

# outdir: the experiment directory holding config.yaml
outdir <- dirname(config_path)
dir.create(file.path(outdir, "analysis-SNVs"), showWarnings = FALSE)
dir.create(file.path(outdir, "analysis-CNVs"), showWarnings = FALSE)
dir.create(file.path(outdir, "analysis-SV"), showWarnings = FALSE)

# Required: metadata_dir (where metadata.csv is written and the hand-edited annotation tables live)
metadata_dir <- cfg$metadata_dir %||% ""
if (!nzchar(metadata_dir)) {
  stop("metadata_dir must be set in ", config_path, call. = FALSE)
}
dir.create(metadata_dir, showWarnings = FALSE, recursive = TRUE)

disease <- cfg$disease %||% NA_character_

# Optional: hand-edited patient-level table (csv keyed by ID), resolved against metadata_dir
# when relative. "" = none.
patient_annotations_path <- cfg$patient_annotations %||% ""
if (nzchar(patient_annotations_path) && !grepl("^/", patient_annotations_path)) {
  patient_annotations_path <- file.path(metadata_dir, patient_annotations_path)
}

# Optional: hand-edited per-sample table (csv keyed by sID), resolved against metadata_dir
# when relative; built by build_sample_annotations.R when missing. A non-empty
# exclude_reason sets QC = "excluded". "" = none.
sample_annotations_path <- cfg$sample_annotations %||% ""
if (nzchar(sample_annotations_path) && !grepl("^/", sample_annotations_path)) {
  sample_annotations_path <- file.path(metadata_dir, sample_annotations_path)
}

# Required: group1 names the column that becomes metadata.csv's group1 — the grouping the
# project compares samples within. Optional: group2, same rule, "" omits the column. Either
# names a samples.tsv-derived column or a patient_annotations column.
group1_col <- cfg$group1 %||% ""
group2_col <- cfg$group2 %||% ""
if (!nzchar(group1_col)) {
  stop("group1 must be set in ", config_path,
       " (name of a metadata column, e.g. sample_type, or a patient_annotations column)",
       call. = FALSE)
}

# Optional: genelist paths (empty string = skip)
snv_genelist          <- cfg$snv_genelist %||% ""
snv_genelist_diagnosis <- cfg$snv_genelist_diagnosis %||% ""
cnv_genelist          <- cfg$cnv_genelist %||% ""

cat("=== load_WES_results.R ===\n")
cat("rootdir:      ", rootdir, "\n")
cat("results_dir:  ", results_dir, "\n")
cat("outdir:       ", outdir, "\n")
cat("metadata_dir: ", metadata_dir, "\n")
cat("disease:      ", disease, "\n")
cat("patient_annotations: ", if (nzchar(patient_annotations_path)) patient_annotations_path else "(none)", "\n")
cat("sample_annotations: ", if (nzchar(sample_annotations_path)) sample_annotations_path else "(none)", "\n")
cat("group1 <- ", group1_col, "\n")
cat("group2 <- ", if (nzchar(group2_col)) group2_col else "(none — column omitted)", "\n")
cat("snv_genelist: ", if (nzchar(snv_genelist)) snv_genelist else "(none — skipping)", "\n")
cat("cnv_genelist: ", if (nzchar(cnv_genelist)) cnv_genelist else "(none — skipping)", "\n\n")

# --- Refresh metadata.csv from WES-snakemake results/metadata/samples.tsv ---------------
# samples.tsv already carries sample_type/gender/capture_kit/tumor_fraction per sID. The
# project-level additions are calling_mode (derived from whether the patient has a CTRL
# sample), every column of patient_annotations but notes (left-joined by ID), Disease (from
# patient_annotations where given, else the config constant), and group1/group2 (copied
# from the columns the config names). An ID absent from samples.tsv is an error. Rows
# sharing an ID are collapsed to one, joining differing values with " | ". A column that
# repeats a samples.tsv column must agree with it per patient and is then dropped.
# sample_annotations is then left-joined by sID under the same rules (sIDs are unique there;
# a duplicate is an error). QC is "excluded" where exclude_reason is set, else NA until
# filter_SNVs.R derives it.
samples_path <- file.path(results_dir, "metadata", "samples.tsv")
if (!file.exists(samples_path)) {
  stop("samples.tsv not found at ", samples_path, call. = FALSE)
}
{
  .ss <- read.delim(samples_path, stringsAsFactors = FALSE)

  .ctrl_ids <- unique(.ss$ID[.ss$sample_type == "CTRL"])
  .ss$calling_mode <- ifelse(.ss$ID %in% .ctrl_ids, "tumor-normal", "tumor-only")

  .meta <- .ss[.ss$sample_type != "CTRL", ]
  names(.meta)[names(.meta) == "sample"] <- "sID"
  .generated_cols <- c("ID", "sID", "sample_type", "gender", "capture_kit",
                       "tumor_fraction", "Disease", "calling_mode")

  .patient_cols <- character(0)
  if (nzchar(patient_annotations_path)) {
    if (!file.exists(patient_annotations_path)) {
      stop("patient_annotations not found: ", patient_annotations_path, call. = FALSE)
    }
    .pat <- read.csv(patient_annotations_path, stringsAsFactors = FALSE,
                     check.names = FALSE, na.strings = c("", "NA"))
    if (!"ID" %in% names(.pat)) {
      stop("patient_annotations (", patient_annotations_path, ") has no ID column", call. = FALSE)
    }
    .unknown <- setdiff(.pat$ID, .ss$ID)
    if (length(.unknown) > 0) {
      stop("patient_annotations: ID(s) not in samples.tsv: ",
           paste(.unknown, collapse = ", "), call. = FALSE)
    }
    .reserved <- intersect(names(.pat), c("sID", "group1", "group2", "QC"))
    if (length(.reserved) > 0) {
      stop("patient_annotations: reserved column name(s): ",
           paste(.reserved, collapse = ", "), call. = FALSE)
    }

    .dup <- unique(.pat$ID[duplicated(.pat$ID)])
    if (length(.dup) > 0) {
      warning("patient_annotations: collapsed duplicated ID(s) ", paste(.dup, collapse = ", "),
              "; differing values joined with \" | \"", call. = FALSE)
      .collapse <- function(x) {
        v <- unique(as.character(x[!is.na(x)]))
        if (length(v) == 0) NA_character_ else paste(v, collapse = " | ")
      }
      .pat <- do.call(rbind, lapply(split(.pat, factor(.pat$ID, levels = unique(.pat$ID))),
                                    function(d) as.data.frame(lapply(d, .collapse),
                                                              check.names = FALSE)))
      .pat[] <- lapply(.pat, type.convert, as.is = TRUE)
      rownames(.pat) <- NULL
      rm(.collapse)
    }

    # notes are for the curator only and stay out of metadata.csv
    .pat$notes <- NULL

    # Columns repeating samples.tsv: each non-NA patient value must equal the single value
    # samples.tsv carries for that patient.
    .shared <- setdiff(intersect(names(.pat), names(.meta)), "ID")
    for (.col in .shared) {
      .bad <- vapply(seq_len(nrow(.pat)), function(i) {
        v <- .pat[[.col]][i]
        if (is.na(v)) return(FALSE)
        g <- unique(as.character(.meta[[.col]][.meta$ID == .pat$ID[i]]))
        !identical(g, as.character(v))
      }, logical(1))
      if (any(.bad)) {
        stop("patient_annotations: column '", .col, "' disagrees with samples.tsv for ",
             paste(.pat$ID[.bad], collapse = ", "), call. = FALSE)
      }
    }
    .pat <- .pat[, setdiff(names(.pat), .shared), drop = FALSE]

    .absent <- setdiff(unique(.meta$ID), .pat$ID)
    if (length(.absent) > 0) {
      warning("patients missing from patient_annotations: ", paste(.absent, collapse = ", "),
              call. = FALSE)
    }
    .meta <- merge(.meta, .pat, by = "ID", all.x = TRUE)
    .patient_cols <- setdiff(names(.pat), c("ID", "Disease"))
    rm(.pat, .dup, .unknown, .reserved, .shared, .absent)
    rm(list = intersect(c(".col", ".bad"), ls(all.names = TRUE)))
  }
  if (!"Disease" %in% names(.meta)) .meta$Disease <- NA_character_
  .meta$Disease[is.na(.meta$Disease)] <- disease

  .sample_cols <- character(0)
  if (nzchar(sample_annotations_path)) {
    if (!file.exists(sample_annotations_path)) {
      source(file.path(rootdir, "build_sample_annotations.R"), local = TRUE)
      build_sample_annotations(samples_path, sample_annotations_path)
    }
    .san <- read.csv(sample_annotations_path, stringsAsFactors = FALSE, check.names = FALSE,
                     na.strings = c("", "NA"))
    if (!"sID" %in% names(.san)) {
      stop("sample_annotations (", sample_annotations_path, ") has no sID column", call. = FALSE)
    }
    .dup <- unique(.san$sID[duplicated(.san$sID)])
    if (length(.dup) > 0) {
      stop("sample_annotations: duplicated sID(s): ", paste(.dup, collapse = ", "), call. = FALSE)
    }
    .unknown <- setdiff(.san$sID, .meta$sID)
    if (length(.unknown) > 0) {
      stop("sample_annotations: sID(s) not among samples.tsv tumor samples: ",
           paste(.unknown, collapse = ", "),
           " (./runall --sample-annotations drops them)", call. = FALSE)
    }
    .reserved <- intersect(names(.san), c("group1", "group2", "QC", .patient_cols))
    if (length(.reserved) > 0) {
      stop("sample_annotations: column name(s) reserved or already in patient_annotations: ",
           paste(.reserved, collapse = ", "), call. = FALSE)
    }
    .absent <- setdiff(.meta$sID, .san$sID)
    if (length(.absent) > 0) {
      warning("samples missing from sample_annotations (./runall --sample-annotations adds them): ",
              paste(.absent, collapse = ", "), call. = FALSE)
    }
    .shared <- setdiff(intersect(names(.san), names(.meta)), "sID")
    for (.col in .shared) {
      .g <- as.character(.meta[[.col]][match(.san$sID, .meta$sID)])
      .v <- as.character(.san[[.col]])
      .bad <- !is.na(.v) & (is.na(.g) | .g != .v)
      if (any(.bad)) {
        stop("sample_annotations: column '", .col, "' disagrees with generated metadata for ",
             paste(.san$sID[.bad], collapse = ", "), call. = FALSE)
      }
    }
    .san <- .san[, setdiff(names(.san), .shared), drop = FALSE]
    if (!"exclude_reason" %in% names(.san)) .san$exclude_reason <- NA_character_
    .meta <- merge(.meta, .san, by = "sID", all.x = TRUE)
    .sample_cols <- setdiff(names(.san), "sID")
    rm(.san, .dup, .unknown, .reserved, .absent, .shared)
    rm(list = intersect(c(".col", ".g", ".v", ".bad"), ls(all.names = TRUE)))
  }
  .meta$QC <- if ("exclude_reason" %in% names(.meta)) {
    ifelse(!is.na(.meta$exclude_reason) & nzchar(trimws(.meta$exclude_reason)),
           "excluded", NA_character_)
  } else {
    NA_character_
  }

  .group_cols <- c(group1_col, if (nzchar(group2_col)) group2_col)
  .missing <- setdiff(.group_cols, names(.meta))
  if (length(.missing) > 0) {
    stop("group column(s) ", paste(.missing, collapse = ", "),
         " in neither samples.tsv nor patient_annotations", call. = FALSE)
  }
  .grp_vals <- lapply(.group_cols, function(col) .meta[[col]])
  .meta$group1 <- .grp_vals[[1]]
  if (nzchar(group2_col)) .meta$group2 <- .grp_vals[[2]]
  .n_no_group <- sum(is.na(.meta$group1))
  if (.n_no_group > 0) {
    warning(.n_no_group, " samples have no group1 (", group1_col, " is NA): ",
            paste(unique(.meta$sID[is.na(.meta$group1)]), collapse = ", "), call. = FALSE)
  }
  rm(.group_cols, .missing, .grp_vals, .n_no_group)

  .meta_path <- file.path(metadata_dir, "metadata.csv")
  .col_order <- c(.generated_cols, "group1", "group2", .patient_cols, .sample_cols, "QC")
  .meta <- .meta[, intersect(.col_order, names(.meta))]
  .meta <- .meta[order(.meta$ID, .meta$sID), ]

  write.csv(.meta, .meta_path, row.names = FALSE)
  # metadata.ver names the experiment directory that wrote metadata.csv
  writeLines(basename(outdir), file.path(metadata_dir, "metadata.ver"))
  cat(sprintf("metadata.csv updated: %d samples (%d excluded)\n\n", nrow(.meta),
              sum(.meta$QC %in% "excluded")))
  rm(.ss, .ctrl_ids, .meta, .meta_path, .col_order, .generated_cols, .patient_cols,
     .sample_cols)
}

# Source filter scripts — they pick up rootdir, config_path, cfg, combined_dir, outdir,
# metadata_dir, disease, and genelist vars from this environment.
cat(">>> Sourcing filter_SNVs.R ...\n")
source(file.path(rootdir, "analysis-SNVs", "filter_SNVs.R"))

cat("\n>>> Sourcing filter_CNVs.R ...\n")
source(file.path(rootdir, "analysis-CNVs", "filter_CNVs.R"))

cat("\n>>> Sourcing filter_SVs.R ...\n")
source(file.path(rootdir, "analysis-SV", "filter_SVs.R"))

cat("\n=== load_WES_results.R complete ===\n")
