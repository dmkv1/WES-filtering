## strandbias_rescue.R
## Strand-bias rescue rule shared by WES-filtering (filter_SNVs.R) and the rescue comparison
## scripts. A call whose only FilterMutectCalls label is strand_bias is rescued (labelled
## PASS) when all hold:
##   - alt reads are split across strands like the ref reads: two-sided Fisher exact test
##     on (ref fwd, ref rev, alt fwd, alt rev) gives p >= fisher_p_min
##   - AD_ALT >= min_alt_reads
##   - AF >= af_min (when af is given)
##   - fewer than min_runs other runs (patients) carry the same call

# Two-sided Fisher p of a n x 4 matrix (ref fwd, ref rev, alt fwd, alt rev); NA for rows with
# missing counts.
sb_fisher_p <- function(sb) {
  sb <- as.matrix(sb)
  vapply(seq_len(nrow(sb)), function(i) {
    x <- sb[i, ]
    if (anyNA(x)) return(NA_real_)
    fisher.test(matrix(x, nrow = 2))$p.value
  }, numeric(1))
}

# Number of runs other than the row's own (ID) that carry the row's variant key in the
# pool (pool_key, pool_id).
n_other_runs <- function(key, id, pool_key, pool_id) {
  pool <- unique(data.frame(key = pool_key, id = pool_id, stringsAsFactors = FALSE))
  n_by_key <- table(pool$key)
  n_all <- as.integer(n_by_key[key]); n_all[is.na(n_all)] <- 0L
  own <- paste(key, id) %in% paste(pool$key, pool$id)
  n_all - own
}

strandbias_rescue_mask <- function(filter, sb, alt_reads, n_other,
                                   fisher_p_min, min_alt_reads, min_runs,
                                   af = NULL, af_min = 0) {
  af_ok <- if (is.null(af)) TRUE else (!is.na(af) & af >= af_min)
  filter == "strand_bias" & af_ok &
    !is.na(sb_fisher_p_safe(filter, sb)) &
    sb_fisher_p_safe(filter, sb) >= fisher_p_min &
    alt_reads >= min_alt_reads &
    n_other < min_runs
}

# Fisher p only for strand_bias rows (the only rows the rule can rescue)
sb_fisher_p_safe <- function(filter, sb) {
  p <- rep(NA_real_, length(filter))
  i <- which(filter == "strand_bias")
  if (length(i)) p[i] <- sb_fisher_p(sb[i, , drop = FALSE])
  p
}
