## build_sample_annotations.R
## Builds the hand-edited per-sample annotation table (config.yaml's sample_annotations,
## default metadata/sample_annotations.csv): one row per tumor sample in the WES-snakemake
## samples.tsv, keyed by sID, with columns exclude_reason and sample_notes.
##
## Rebuilding never discards edits: rows of sIDs still in samples.tsv keep every value and
## every extra column, new sIDs get empty rows, and rows of sIDs no longer in samples.tsv
## are dropped with a warning naming them.
##
## load_WES_results.R sources this file and builds the table when it does not exist.
##
## Usage (rebuild):
##   Rscript build_sample_annotations.R --config results/<date>_WES-filtering/config.yaml

build_sample_annotations <- function(samples_path, annotations_path) {
  ss <- read.delim(samples_path, stringsAsFactors = FALSE)
  sids <- ss$sample[ss$sample_type != "CTRL"]
  sids <- sids[order(ss$ID[ss$sample_type != "CTRL"], sids)]

  fresh <- data.frame(sID = sids, exclude_reason = NA_character_, sample_notes = NA_character_,
                      stringsAsFactors = FALSE)

  if (file.exists(annotations_path)) {
    old <- read.csv(annotations_path, stringsAsFactors = FALSE, check.names = FALSE,
                    na.strings = c("", "NA"), colClasses = "character")
    stale <- setdiff(old$sID, sids)
    if (length(stale) > 0) {
      warning("sample_annotations: dropped sID(s) no longer in samples.tsv: ",
              paste(stale, collapse = ", "), call. = FALSE)
    }
    new_sids <- setdiff(sids, old$sID)
    for (col in setdiff(names(old), names(fresh))) fresh[[col]] <- NA_character_
    for (col in setdiff(names(fresh), names(old))) old[[col]] <- NA_character_
    kept <- old[old$sID %in% sids, names(fresh), drop = FALSE]
    out  <- rbind(kept, fresh[fresh$sID %in% new_sids, , drop = FALSE])
    out  <- out[match(sids, out$sID), , drop = FALSE]
    cat(sprintf("sample_annotations rebuilt: %d samples (%d added, %d dropped) -> %s\n",
                nrow(out), length(new_sids), length(stale), annotations_path))
  } else {
    out <- fresh
    cat(sprintf("sample_annotations created: %d samples -> %s\n", nrow(out), annotations_path))
  }

  write.csv(out, annotations_path, row.names = FALSE, na = "")
  invisible(out)
}

if (sys.nframe() == 0) {
  .script_path <- sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
  source(file.path(dirname(normalizePath(.script_path)), "config.R"))
  cfg <- read_wes_config(config_arg("build_sample_annotations.R"))
  config_path <- cfg$config_path
  metadata_dir <- cfg$metadata_dir %||% stop("metadata_dir must be set in ", config_path,
                                             call. = FALSE)
  annotations_path <- cfg$sample_annotations %||% stop("sample_annotations must be set in ",
                                                       config_path, call. = FALSE)
  if (!grepl("^/", annotations_path)) {
    annotations_path <- file.path(metadata_dir, annotations_path)
  }
  build_sample_annotations(file.path(cfg$results_dir, "metadata", "samples.tsv"),
                           annotations_path)
}
