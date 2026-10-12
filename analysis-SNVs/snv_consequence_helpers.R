## snv_consequence_helpers.R
## Shared VEP-Consequence helpers: filter_SNVs.R's allowlist, plot_snv_filters.R's
## coding split, and any downstream metrics/statistics code. Source this file rather
## than redefining these vectors/functions per script, so the coding definition can't
## drift between them. `outdir` must point at the WES-filtering experiment directory
## holding SNV_consequence_terms.tsv.
##
## VEP's Consequence can be compound ("&"-joined SO terms co-occurring on the same
## transcript). Two different questions need two different matching strategies:
##   - "is this variant coding at all" (inclusion test) -> ANY term matches
##     (same logic filter_SNVs.R already uses for SNV_allowlist)
##   - "which single bucket does this variant go in" (classification, one row one
##     bucket) -> the PRIMARY term only (VEP orders terms most-severe-first), so a
##     stacked count/proportion table doesn't double-count a variant across buckets

## <outdir>/SNV_consequence_terms.tsv is the single source of truth for both VEP
## Consequence term sets used across this pipeline: `pass_consequence`
## (filter_SNVs.R's allowlist — coding + 3' UTR) and `coding_consequence` (the
## stricter coding-only set used here). coding_consequence is the VEP SO-term
## equivalent of the archive's Funcotator coding_types (MISSENSE, NONSENSE,
## FRAME_SHIFT_DEL/INS, SPLICE_SITE, IN_FRAME_INS/DEL, NONSTOP, plus
## protein_altering_variant for VEP calls Funcotator's classifier resolves
## differently) and start_lost. Every coding_consequence
## term is also pass_consequence, but not vice versa (3' UTR terms are
## allowlisted without being coding) — the two columns are not interchangeable.
read_consequence_terms <- function(outdir) {
  path <- file.path(outdir, "SNV_consequence_terms.tsv")
  terms <- read.delim(path, stringsAsFactors = FALSE)
  terms$pass_consequence <- as.logical(terms$pass_consequence)
  terms$coding_consequence <- as.logical(terms$coding_consequence)
  terms
}

.consequence_terms_cache <- new.env(parent = emptyenv())

.coding_consequence_terms <- function() {
  if (is.null(.consequence_terms_cache$terms)) {
    if (!exists("outdir", envir = .GlobalEnv)) {
      stop("snv_consequence_helpers.R: `outdir` must be set before calling is_coding_consequence().",
           call. = FALSE)
    }
    .consequence_terms_cache$terms <- read_consequence_terms(get("outdir", envir = .GlobalEnv))
  }
  with(.consequence_terms_cache$terms, term[coding_consequence])
}

consequence_has_any <- function(consequence, term_set) {
  terms <- strsplit(consequence, "&", fixed = TRUE)
  vapply(terms, function(x) any(x %in% term_set), logical(1))
}

is_coding_consequence <- function(consequence) consequence_has_any(consequence, .coding_consequence_terms())

## is_qc_pass: TRUE for every variant that cleared all quality filters, regardless of
## consequence — i.e. FILTER == "PASS" (QC-clean, allowlisted consequence) OR
## FILTER == "low_priority_consequence" (QC-clean, non-allowlisted consequence: the
## passenger variants filter_SNVs.R reclassifies out of PASS but does not drop).
## This is the "total" burden population; is_coding_consequence() further restricts
## within it to the strict coding subset.
is_qc_pass <- function(FILTER) FILTER %in% c("PASS", "low_priority_consequence")

consequence_primary_term <- function(consequence) {
  vapply(strsplit(consequence, "&", fixed = TRUE), function(x) x[1], character(1))
}

## map_consequence_broad: one broad category per variant, from the primary Consequence
## term. Loosely mirrors the archive's map_varclass() bucket names where a clean VEP
## equivalent exists; Protein_Altering, Splice_Region, Upstream, Downstream and
## Noncoding_Transcript have no Funcotator-classic counterpart and are new buckets.
map_consequence_broad <- function(consequence) {
  primary <- consequence_primary_term(consequence)
  dplyr::case_when(
    primary == "missense_variant"                                     ~ "Missense",
    primary == "stop_gained"                                          ~ "Nonsense",
    primary == "frameshift_variant"                                   ~ "Frameshift",
    primary == "stop_lost"                                            ~ "Nonstop",
    primary == "start_lost"                                           ~ "Start_Codon",
    primary %in% c("inframe_insertion", "inframe_deletion")           ~ "In_Frame",
    primary == "protein_altering_variant"                             ~ "Protein_Altering",
    primary %in% c("splice_acceptor_variant", "splice_donor_variant") ~ "Splice_Site",
    primary %in% c("splice_donor_region_variant",
                   "splice_polypyrimidine_tract_variant",
                   "splice_region_variant")                           ~ "Splice_Region",
    primary == "synonymous_variant"                                   ~ "Silent",
    primary == "3_prime_UTR_variant"                                  ~ "3'UTR",
    primary == "5_prime_UTR_variant"                                  ~ "5'UTR",
    primary == "intron_variant"                                       ~ "Intron",
    primary == "upstream_gene_variant"                                ~ "Upstream",
    primary == "downstream_gene_variant"                              ~ "Downstream",
    primary %in% c("non_coding_transcript_exon_variant",
                   "non_coding_transcript_variant")                   ~ "Noncoding_Transcript",
    TRUE                                                               ~ "Other"
  )
}

## sub6_from_variant: pyrimidine-convention 6-class substitution from the Variant column
## ("REF/ALT", e.g. "C/T") — the archive used Funcotator's refAllele/tumorSeqAllele2 pair,
## which has no equivalent here. Only meaningful for VARIANT_CLASS == "SNV" rows (ref/alt
## are multi-character for indels).
sub6_from_variant <- function(variant_col) {
  parts <- strsplit(variant_col, "/", fixed = TRUE)
  ref <- vapply(parts, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1))
  alt <- vapply(parts, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
  dplyr::case_when(
    ref == "C" & alt == "A" ~ "C>A", ref == "C" & alt == "G" ~ "C>G", ref == "C" & alt == "T" ~ "C>T",
    ref == "T" & alt == "A" ~ "T>A", ref == "T" & alt == "C" ~ "T>C", ref == "T" & alt == "G" ~ "T>G",
    ref == "G" & alt == "T" ~ "C>A", ref == "G" & alt == "C" ~ "C>G", ref == "G" & alt == "A" ~ "C>T",
    ref == "A" & alt == "T" ~ "T>A", ref == "A" & alt == "G" ~ "T>C", ref == "A" & alt == "C" ~ "T>G",
    TRUE ~ NA_character_
  )
}
