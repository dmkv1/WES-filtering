## config.R
## read_wes_config(): an experiment's config.yaml layered over template/config.yaml, which
## holds the defaults for every key. Nested maps (snv_filter_params, ...) merge key by key;
## any other value set in the experiment config replaces the default. sv_driver_loci is
## taken from the experiment config alone, so loci are never inherited.
##
## Path keys are resolved against the experiment directory (the directory holding
## config.yaml) when relative; absolute paths are kept. An empty cytoband selects the
## bundled resources/cytoBand.hg38.txt.gz. patient_annotations and sample_annotations
## stay relative to metadata_dir; flags_genelist stays relative to the experiment directory.

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && nzchar(a)) a else b

.wes_repo_dir <- function() {
  ofile <- if (sys.nframe() > 0) {
    frames <- rev(sys.frames())
    hit <- Filter(function(f) !is.null(f$ofile) && basename(f$ofile) == "config.R", frames)
    if (length(hit)) hit[[1]]$ofile else NULL
  }
  if (is.null(ofile)) stop("config.R must be loaded with source()", call. = FALSE)
  normalizePath(dirname(ofile))
}
.wes_repo <- .wes_repo_dir()

.wes_path_keys     <- c("results_dir", "metadata_dir", "cytoband", "snv_genelist", "cnv_genelist")
.wes_snv_path_keys <- c("hotspot_list", "clustered_rescue_genes", "strandbias_recurrence_reference")

read_wes_config <- function(config_path) {
  if (!file.exists(config_path)) stop("Config file not found: ", config_path, call. = FALSE)
  config_path <- normalizePath(config_path)
  base <- dirname(config_path)

  defaults <- yaml::read_yaml(file.path(.wes_repo, "template", "config.yaml"))
  user <- yaml::read_yaml(config_path)
  if (is.null(user)) user <- list()
  cfg <- utils::modifyList(defaults, user)
  cfg$sv_driver_loci <- user$sv_driver_loci

  resolve <- function(p) {
    if (is.null(p) || !nzchar(p) || grepl("^/", p)) return(p)
    normalizePath(file.path(base, p), mustWork = FALSE)
  }
  if (!nzchar(cfg$cytoband %||% "")) {
    cfg$cytoband <- file.path(.wes_repo, "resources", "cytoBand.hg38.txt.gz")
  }
  for (k in .wes_path_keys) cfg[[k]] <- resolve(cfg[[k]])
  for (k in .wes_snv_path_keys) cfg$snv_filter_params[[k]] <- resolve(cfg$snv_filter_params[[k]])

  cfg$config_path <- config_path
  cfg
}

# Value of --config in commandArgs(); stops with usage when absent.
config_arg <- function(script) {
  args <- commandArgs(trailingOnly = TRUE)
  i <- match("--config", args)
  if (is.na(i) || i == length(args)) {
    stop("--config is required (config.yaml lives in the experiment's results/ directory).\n",
         "  Usage: Rscript ", script, " --config /path/to/config.yaml", call. = FALSE)
  }
  args[i + 1]
}
