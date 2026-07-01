# =============================================================================
# 4_data_prep.R
# Load PCR-collapsed phyloseq objects (unfiltered by rg2 — see design decision
# below), remove non-fungal OTUs via PlutoF SH matching, filter low-depth
# samples, rarefy, fix taxonomy labels, and save analysis-ready datasets.
# Also prepares the full-complexity (per-PCR-replicate) object for models that
# use PCR replicate as a random effect (GLLVM/HMSC), and keeps the rg2-filtered
# objects as a separate, clearly labeled sensitivity check.
#
# Canonical object design (locked 2026-07-01, revised 2026-07-01):
#   alldat / alldat.rfy — PCR replicates collapsed by summing, NO rg2 filter.
#     Primary object for all methods needing independent sample rows
#     (PERMANOVA/dbRDA, iNEXT3D, NRI/NTI, UniFrac).
#   alldat_full — one row per PCR replicate, never collapsed (nopool/pool/
#     pspool). Runs through the exact same non-fungal removal, PlutoF
#     taxonomy upgrade, and taxonomy cleaning as `alldat`; the depth filter
#     uses a very low, effectively pass-through threshold (min_depth_full)
#     since GLLVM/HMSC use PCR replicate as a random effect and don't need
#     fixed-depth filtering the way independent-sample methods do. NOT
#     rarefied — use a per-row read-count offset in the model instead.
#   alldat_rg2_sensitivity — same as alldat but built from the rg2-filtered
#     (OTU must appear in >=2 replicates) objects. Kept only as an optional
#     robustness check; not part of the primary analysis path.
# Collapsing is irreversible (summed counts can't be un-collapsed), so both
# the collapsed and full-complexity objects are persisted here rather than
# derived on demand downstream.
#
# Input:  nopoolps.dada2.rds / poolps.dada2.rds / pspoolps.dada2.rds
#         rg2.nopoolps.tax.rds / rg2.poolps.tax.rds / rg2.pspoolps.tax.rds
#         fullps_nopool.rds / fullps_pool.rds / fullps_pspool.rds
#         PlutoF output CSVs: matches_out_taxonomy_{nopool,pool,pspool}.csv
# Output: eco_analysis.RData              (analysis-ready workspace)
#         phyloseq_comparison_summary.csv (per-step object metrics)
#         nonfungal_reads_summary.csv     (non-fungal removal accounting)
#         nonfungal_kingdom_breakdown.csv (kingdom-level breakdown of removed OTUs)
#
# Reference: 4_data_prep.R @ 65fbffa (Root_fungi_DADA2)
# =============================================================================

library(phyloseq)
library(microViz)
library(fantaxtic)
library(ggplot2)
library(vegan)
library(dplyr)
library(tidyr)

# =============================================================================
# SECTION 1 — LOAD PHYLOSEQ OBJECTS (canonical, sensitivity, full-complexity)
# =============================================================================

setwd("/home/daniel/Ptarmigan/trimmed/mergedPlates/")

# Primary analysis objects: PCR-collapsed (summed), NO rg2 filter,
# taxonomy-annotated (minBoot=80). Built in 2_DADA2_lulu.R Section 5
# (make.phylo, raw/no-rg2 branch); taxonomy in 3_taxonomic_assignment.R
# Section 2. Metadata (Season, Year, TimeBlock, indivID, sex, species,
# sampletype, etc.) is embedded from ITS_metadata_ptarmigan_clean.xlsx.
# nopool is the PRIMARY strategy; pool and pspool are kept for comparison.
nopoolps.dada2  <- readRDS("nopoolps.dada2.rds")
poolps.dada2    <- readRDS("poolps.dada2.rds")
pspoolps.dada2  <- readRDS("pspoolps.dada2.rds")

# rg2-filtered counterparts (OTU must appear in >=2 PCR replicates), kept only
# as an optional sensitivity check — see header note.
rg2.nopoolps.tax  <- readRDS("rg2.nopoolps.tax.rds")
rg2.poolps.tax    <- readRDS("rg2.poolps.tax.rds")
rg2.pspoolps.tax  <- readRDS("rg2.pspoolps.tax.rds")

# Full-complexity counterparts (one row per PCR replicate, never collapsed).
# Built in 2_DADA2_lulu.R Section 5b; taxonomy in 3_taxonomic_assignment.R
# Section 3b. Treated through the same pipeline steps as alldat below, just
# with a pass-through depth filter — see header note.
fullps_nopool  <- readRDS("fullps_nopool.rds")
fullps_pool    <- readRDS("fullps_pool.rds")
fullps_pspool  <- readRDS("fullps_pspool.rds")

stopifnot(all(vapply(
  list(nopoolps.dada2, poolps.dada2, pspoolps.dada2,
       rg2.nopoolps.tax, rg2.poolps.tax, rg2.pspoolps.tax,
       fullps_nopool, fullps_pool, fullps_pspool),
  function(ps) inherits(ps, "phyloseq") && !is.null(tax_table(ps, errorIfNULL=FALSE)),
  logical(1)
)))

cat("Objects loaded:\n")
cat("  nopoolps.dada2:   ", nsamples(nopoolps.dada2),   "samples,",
    ntaxa(nopoolps.dada2),   "taxa\n")
cat("  poolps.dada2:     ", nsamples(poolps.dada2),     "samples,",
    ntaxa(poolps.dada2),     "taxa\n")
cat("  pspoolps.dada2:   ", nsamples(pspoolps.dada2),   "samples,",
    ntaxa(pspoolps.dada2),   "taxa\n")
cat("  rg2.nopoolps.tax: ", nsamples(rg2.nopoolps.tax), "samples,",
    ntaxa(rg2.nopoolps.tax), "taxa\n")
cat("  rg2.poolps.tax:   ", nsamples(rg2.poolps.tax),   "samples,",
    ntaxa(rg2.poolps.tax),   "taxa\n")
cat("  rg2.pspoolps.tax: ", nsamples(rg2.pspoolps.tax), "samples,",
    ntaxa(rg2.pspoolps.tax), "taxa\n")
cat("  fullps_nopool:    ", nsamples(fullps_nopool),    "PCR-replicate samples,",
    ntaxa(fullps_nopool),    "taxa\n")
cat("  fullps_pool:      ", nsamples(fullps_pool),      "PCR-replicate samples,",
    ntaxa(fullps_pool),      "taxa\n")
cat("  fullps_pspool:    ", nsamples(fullps_pspool),    "PCR-replicate samples,",
    ntaxa(fullps_pspool),    "taxa\n")

# Quick metadata check
meta_check <- as.data.frame(sample_data(nopoolps.dada2))
cat("\nSeason distribution (nopool):\n")
print(table(meta_check$Season, useNA="always"))
cat("Year distribution (nopool):\n")
print(table(meta_check$Year, useNA="always"))
cat("TimeBlock distribution (nopool):\n")
print(table(meta_check$TimeBlock, useNA="always"))

# =============================================================================
# SECTION 2 — BUILD ANALYSIS LISTS
# =============================================================================

# alldat: named list; nopool = primary analysis strategy; canonical (no rg2)
alldat <- list(
  nopool = nopoolps.dada2,
  pool   = poolps.dada2,
  pspool = pspoolps.dada2
)

# alldat_rg2_sensitivity: same strategies, rg2-filtered — robustness check only
alldat_rg2_sensitivity <- list(
  nopool = rg2.nopoolps.tax,
  pool   = rg2.poolps.tax,
  pspool = rg2.pspoolps.tax
)

# alldat_full: same strategies, full-complexity (one row per PCR replicate) —
# for GLLVM/HMSC. Treated identically to alldat/alldat_rg2_sensitivity below.
alldat_full <- list(
  nopool = fullps_nopool,
  pool   = fullps_pool,
  pspool = fullps_pspool
)

alldat <- lapply(alldat, function(x) phyloseq_validate(x, remove_undetected=TRUE))
alldat_rg2_sensitivity <- lapply(alldat_rg2_sensitivity, function(x) phyloseq_validate(x, remove_undetected=TRUE))
alldat_full <- lapply(alldat_full, function(x) phyloseq_validate(x, remove_undetected=TRUE))

# Capture initial state for comparison table
summarise_ps <- function(ps, object_label, step) {
  sums <- sample_sums(ps)
  data.frame(
    step         = step,
    object       = object_label,
    n_samples    = nsamples(ps),
    n_taxa       = ntaxa(ps),
    total_reads  = sum(sums),
    min_depth    = min(sums),
    max_depth    = max(sums),
    median_depth = round(median(sums)),
    stringsAsFactors = FALSE
  )
}

step_initial <- dplyr::bind_rows(lapply(names(alldat), function(nm) {
  summarise_ps(alldat[[nm]], nm, "initial")
}))
step_initial_full <- dplyr::bind_rows(lapply(names(alldat_full), function(nm) {
  summarise_ps(alldat_full[[nm]], paste0(nm, "_full"), "initial")
}))

# =============================================================================
# SECTION 3 — REMOVE NON-FUNGAL (HOST/PLANT) OTUs
#
# Taxonomy-unresolved OTUs (no Phylum OR Phylum but no Class) were exported to
# PlutoF SH matching v2.0.0 for BLAST verification. Confirmed non-fungal OTUs
# (mostly Viridiplantae — ingested plant ITS2 in dung) are removed here.
#
# Export was run from nopoolps.dada2 / poolps.dada2 / pspoolps.dada2 — the
# same (no-rg2) objects that are now `alldat` — so the OTU universes match
# directly; results are also applied to alldat_rg2_sensitivity and alldat_full,
# which share the same OTU ID space (intersect() guards against any mismatch).
#
# NOTE: The pool PlutoF run returned 0 matches (run failed on PlutoF server).
#       nopool and pspool results are used. Since nopool is the primary strategy
#       this is acceptable; pool will show no PlutoF-based non-fungal removal
#       (UNITE taxonomy-based filtering still applies via Kingdom assignment).
#
# PlutoF output files (exact names required):
#   nopool -> matches_out_taxonomy_nopool.csv  (270 rows — OK)
#   pool   -> matches_out_taxonomy_pool.csv    (0 rows  — FAILED, skipped)
#   pspool -> matches_out_taxonomy_pspool.csv  (287 rows — OK)
# =============================================================================

# ---- Helper: extract sequence map from refseq slot --------------------------
get_taxa_sequence_map <- function(ps, object_label) {
  seqs <- tryCatch(phyloseq::refseq(ps), error = function(e) NULL)
  if (!is.null(seqs)) {
    return(stats::setNames(as.character(seqs), taxa_names(ps)))
  }
  if (all(grepl("^[ACGTN]+$", taxa_names(ps)))) {
    return(stats::setNames(taxa_names(ps), taxa_names(ps)))
  }
  stop(sprintf(
    "%s has no reference sequences and taxon names are not sequences — cannot export FASTA.",
    object_label
  ))
}

# ---- 3A) Export FASTA of taxonomy-unresolved OTUs for PlutoF ----------------
# Selects: (1) OTUs with Kingdom but no Phylum, and
#          (2) OTUs with Phylum but no Class.
# Both groups may contain plant or other non-fungal sequences misclassified by
# UNITE's naive Bayes at Kingdom level but unresolvable below that.
# NOTE: This export has already been run. PlutoF results are in the working dir.
export_unresolved_fasta <- function(ps, object_label, out_dir = ".") {
  tt      <- as.data.frame(tax_table(ps))
  otu_mat <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)
  reads <- rowSums(otu_mat)

  no_phylum <- is.na(tt$Phylum)
  no_class  <- !is.na(tt$Phylum) & is.na(tt$Class)
  sel       <- no_phylum | no_class

  otu_ids <- taxa_names(ps)[sel]
  seq_map <- get_taxa_sequence_map(ps, object_label)

  export_tbl <- data.frame(
    otu_id       = otu_ids,
    reads_total  = as.numeric(reads[otu_ids]),
    flag         = ifelse(no_phylum[sel], "no_phylum", "phylum_no_class"),
    Kingdom      = tt$Kingdom[sel],
    Phylum       = tt$Phylum[sel],
    Class        = tt$Class[sel],
    stringsAsFactors = FALSE
  )
  export_tbl <- export_tbl[order(export_tbl$reads_total, decreasing = TRUE), ]

  fasta_path <- file.path(out_dir, sprintf("plutof_unresolved_%s.fasta", object_label))
  table_path <- file.path(out_dir, sprintf("plutof_unresolved_%s_table.csv", object_label))

  fasta_lines <- as.vector(rbind(
    paste0(">", export_tbl$otu_id),
    seq_map[export_tbl$otu_id]
  ))
  writeLines(fasta_lines, con = fasta_path)
  write.csv(export_tbl, file = table_path, row.names = FALSE)

  cat(sprintf(
    "  %s: %d unresolved OTUs exported (%d no-phylum | %d phylum-no-class) — %.1f%% of reads\n",
    object_label,
    nrow(export_tbl),
    sum(no_phylum[sel]),
    sum(no_class[sel]),
    100 * sum(export_tbl$reads_total) / sum(reads)
  ))
  cat(sprintf("    -> %s\n", fasta_path))

  invisible(data.frame(
    object     = object_label,
    n_exported = nrow(export_tbl),
    n_no_phylum       = sum(no_phylum[sel]),
    n_phylum_no_class = sum(no_class[sel]),
    reads_exported    = sum(export_tbl$reads_total),
    pct_reads         = round(100 * sum(export_tbl$reads_total) / sum(reads), 1),
    fasta_file = fasta_path,
    table_file = table_path,
    stringsAsFactors = FALSE
  ))
}

# Taxonomy-annotated, no-rg2 objects (already loaded in Section 1) — these are
# the exact objects PlutoF export runs against and that now underlie `alldat`.
taxa_ps_list <- list(nopool = nopoolps.dada2, pool = poolps.dada2, pspool = pspoolps.dada2)

cat("PlutoF export — taxonomy-unresolved OTUs (no Phylum OR Phylum but no Class):\n")
unresolved_export_report <- dplyr::bind_rows(lapply(names(taxa_ps_list), function(nm) {
  export_unresolved_fasta(taxa_ps_list[[nm]], nm)
}))
print(unresolved_export_report)

# Run PlutoF separately for each FASTA and save outputs using EXACT names:
#   nopool -> matches_out_taxonomy_nopool.csv
#   pool   -> matches_out_taxonomy_pool.csv
#   pspool -> matches_out_taxonomy_pspool.csv

# ---- 3B) Read PlutoF CSV per strategy and classify non-fungal OTUs ---------
# Explicit one-to-one mapping: strategy name -> PlutoF CSV file name
plutof_files <- c(
  nopool = "matches_out_taxonomy_nopool.csv",
  pool   = "matches_out_taxonomy_pool.csv",
  pspool = "matches_out_taxonomy_pspool.csv"
)

find_plutof_path <- function(strategy) {
  if (!strategy %in% names(plutof_files)) {
    stop(sprintf("Unknown strategy '%s'. Expected one of: %s",
                 strategy, paste(names(plutof_files), collapse=", ")))
  }

  expected <- unname(plutof_files[[strategy]])
  if (!file.exists(expected)) {
    stop(sprintf(
      "Could not find PlutoF output for '%s'. Expected exact file: %s",
      strategy,
      expected
    ))
  }

  expected
}

read_nonfungal_ids_for_strategy <- function(strategy, valid_ids) {
  plutof_path <- find_plutof_path(strategy)
  plutof_tax <- read.delim(plutof_path, stringsAsFactors=FALSE, check.names=FALSE)
  required_cols <- c("seq_name", "common_taxonomy")
  stopifnot(all(required_cols %in% colnames(plutof_tax)))

  if (nrow(plutof_tax) == 0) {
    cat(sprintf(
      "%s: PlutoF output is empty (run failed) — no non-fungal OTUs removed for this strategy.\n",
      strategy
    ))
    return(character(0))
  }

  class_tbl <- plutof_tax |>
    dplyr::filter(!is.na(seq_name), seq_name != "") |>
    dplyr::mutate(
      is_fungal = grepl("^k__Fungi(;|$)", common_taxonomy),
      is_nonfungal = !is_fungal
    ) |>
    dplyr::group_by(seq_name) |>
    dplyr::summarise(
      any_nonfungal = any(is_nonfungal),
      any_fungal = any(is_fungal),
      n_rows = dplyr::n(),
      .groups = "drop"
    )

  conflict_ids <- class_tbl |>
    dplyr::filter(any_nonfungal & any_fungal) |>
    dplyr::pull(seq_name)
  if (length(conflict_ids) > 0) {
    cat(sprintf("%s: OTUs with mixed fungal/non-fungal PlutoF rows: %d\n", strategy, length(conflict_ids)))
  }

  nonfungal_ids <- class_tbl |>
    dplyr::filter(any_nonfungal) |>
    dplyr::pull(seq_name)

  nonfungal_ids <- intersect(nonfungal_ids, valid_ids)
  cat(sprintf("%s: non-fungal OTUs flagged (present in object): %d\n", strategy, length(nonfungal_ids)))
  nonfungal_ids
}

# Get non-fungal IDs validated against alldat's taxa (== the no-rg2 taxa
# universe used for the PlutoF export itself)
nonfungal_ids_by_strategy <- setNames(
  lapply(names(alldat), function(nm) {
    read_nonfungal_ids_for_strategy(nm, taxa_names(alldat[[nm]]))
  }),
  names(alldat)
)

# ---- Quantify contamination per object before removal -----------------------
audit_nonfungal <- function(ps, object_label, drop_ids) {
  otu_mat <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)
  present_drop <- intersect(drop_ids, rownames(otu_mat))

  reads_total <- sum(otu_mat)
  reads_nonfungal <- if (length(present_drop) > 0) sum(otu_mat[present_drop, , drop=FALSE]) else 0

  data.frame(
    object = object_label,
    taxa_total = nrow(otu_mat),
    nonfungal_taxa_present = length(present_drop),
    reads_total = reads_total,
    reads_nonfungal = reads_nonfungal,
    pct_reads_nonfungal = if (reads_total > 0) round(100 * reads_nonfungal / reads_total, 2) else 0
  )
}

audit_tbl <- dplyr::bind_rows(lapply(names(alldat), function(nm) {
  audit_nonfungal(alldat[[nm]], nm, nonfungal_ids_by_strategy[[nm]])
}))
cat("Non-fungal OTU audit (before removal):\n")
print(audit_tbl)

# ---- Kingdom breakdown: what is being removed? ------------------------------
get_nonfungal_kingdom_breakdown <- function(strategy, nonfungal_ids) {
  if (length(nonfungal_ids) == 0) {
    return(data.frame(strategy=character(0), kingdom=character(0), n_otus=integer(0),
                      stringsAsFactors=FALSE))
  }
  plutof_path <- find_plutof_path(strategy)
  plutof_tax  <- read.delim(plutof_path, stringsAsFactors=FALSE, check.names=FALSE)
  plutof_tax |>
    dplyr::filter(seq_name %in% nonfungal_ids) |>
    dplyr::mutate(kingdom = sub(";.*", "", common_taxonomy)) |>
    dplyr::distinct(seq_name, kingdom) |>
    dplyr::count(kingdom, name="n_otus") |>
    dplyr::mutate(strategy = strategy) |>
    dplyr::select(strategy, kingdom, n_otus)
}

kingdom_breakdown <- dplyr::bind_rows(lapply(names(alldat), function(nm) {
  get_nonfungal_kingdom_breakdown(nm, nonfungal_ids_by_strategy[[nm]])
}))
cat("\nKingdom breakdown of non-fungal OTUs removed:\n")
print(kingdom_breakdown)

# ---- Remove non-fungal OTUs -------------------------------------------------
remove_nonfungal <- function(ps, drop_ids) {
  prune_taxa(!taxa_names(ps) %in% drop_ids, ps)
}

alldat <- setNames(lapply(names(alldat), function(nm) {
  remove_nonfungal(alldat[[nm]], nonfungal_ids_by_strategy[[nm]])
}), names(alldat))

# Same drop lists apply to the rg2 sensitivity and full-complexity objects
# (same OTU ID space — rg2 only zeroes disagreeing-replicate cells and the
# full-complexity object shares the exact OTU numbering, neither changes OTU
# identity).
alldat_rg2_sensitivity <- setNames(lapply(names(alldat_rg2_sensitivity), function(nm) {
  remove_nonfungal(alldat_rg2_sensitivity[[nm]], nonfungal_ids_by_strategy[[nm]])
}), names(alldat_rg2_sensitivity))

alldat_full <- setNames(lapply(names(alldat_full), function(nm) {
  remove_nonfungal(alldat_full[[nm]], nonfungal_ids_by_strategy[[nm]])
}), names(alldat_full))

# ---- 3C) Upgrade taxonomy for PlutoF-confirmed fungal OTUs ------------------
# For OTUs that UNITE left partially unresolved (NA at Phylum or Class) and
# PlutoF confirms as fungal with finer resolution, fill in the NA ranks from
# PlutoF's common_taxonomy. Rules:
#   (1) Only fill ranks that are NA in the UNITE assignment.
#   (2) Any PlutoF rank that conflicts with an existing UNITE value → skip OTU.
#   (3) Tax table column names and prefix format (k__, p__, c__, etc.) are
#       preserved exactly — no format changes.
# If multiple PlutoF rows exist for one OTU, the most resolved hit is used.
update_tax_from_plutof <- function(ps, strategy) {
  rank_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

  plutof_path <- find_plutof_path(strategy)
  plutof_tax  <- read.delim(plutof_path, stringsAsFactors=FALSE, check.names=FALSE)

  if (nrow(plutof_tax) == 0) {
    cat(sprintf("%s: PlutoF output empty — no taxonomy updates applied.\n", strategy))
    return(ps)
  }

  valid_otus <- taxa_names(ps)

  fungal_rows <- plutof_tax |>
    dplyr::filter(seq_name %in% valid_otus,
                  grepl("^k__Fungi(;|$)", common_taxonomy))

  if (nrow(fungal_rows) == 0) {
    cat(sprintf("%s: No fungal OTUs with PlutoF hits — no taxonomy updates.\n", strategy))
    return(ps)
  }

  # For OTUs with multiple PlutoF rows keep the most resolved fungal hit
  best_hits <- fungal_rows |>
    dplyr::group_by(seq_name) |>
    dplyr::slice_max(common_rank, n=1, with_ties=FALSE) |>
    dplyr::ungroup()

  tt <- as.data.frame(tax_table(ps))
  n_updated <- 0

  for (i in seq_len(nrow(best_hits))) {
    otu_id  <- best_hits$seq_name[i]
    tax_str <- best_hits$common_taxonomy[i]

    # Parse semicolon-separated taxonomy; pad to 7 ranks with NA
    parts  <- strsplit(tax_str, ";")[[1]]
    parsed <- rep(NA_character_, 7)
    parsed[seq_along(parts)] <- parts
    parsed[nchar(trimws(parsed)) == 0] <- NA_character_

    current <- as.character(tt[otu_id, rank_cols])

    # Skip if PlutoF disagrees with any existing non-NA UNITE rank (conflict)
    conflict <- any(
      !is.na(current) & !is.na(parsed) & current != parsed,
      na.rm = TRUE
    )
    if (conflict) next

    # Fill only NA ranks
    updated <- ifelse(is.na(current) & !is.na(parsed), parsed, current)

    if (!identical(updated, current)) {
      tt[otu_id, rank_cols] <- updated
      n_updated <- n_updated + 1
    }
  }

  tax_table(ps) <- tax_table(as.matrix(tt))
  cat(sprintf("%s: %d OTUs had taxonomy improved from PlutoF (fungal, finer resolution).\n",
              strategy, n_updated))
  ps
}

alldat <- setNames(lapply(names(alldat), function(nm) {
  update_tax_from_plutof(alldat[[nm]], nm)
}), names(alldat))

alldat_rg2_sensitivity <- setNames(lapply(names(alldat_rg2_sensitivity), function(nm) {
  update_tax_from_plutof(alldat_rg2_sensitivity[[nm]], nm)
}), names(alldat_rg2_sensitivity))

alldat_full <- setNames(lapply(names(alldat_full), function(nm) {
  update_tax_from_plutof(alldat_full[[nm]], nm)
}), names(alldat_full))

# =============================================================================
# SECTION 4 — NON-FUNGAL READS SUMMARY (export CSV)
# =============================================================================

nonfungal_summary <- audit_tbl |>
  dplyr::mutate(
    taxa_fungal        = taxa_total - nonfungal_taxa_present,
    reads_fungal       = reads_total - reads_nonfungal,
    pct_reads_fungal   = round(100 - pct_reads_nonfungal, 2),
    pct_taxa_nonfungal = round(100 * nonfungal_taxa_present / taxa_total, 2)
  ) |>
  dplyr::select(
    object,
    taxa_total, nonfungal_taxa_present, pct_taxa_nonfungal,
    reads_total, reads_nonfungal, pct_reads_nonfungal,
    taxa_fungal, reads_fungal, pct_reads_fungal
  )

write.csv(nonfungal_summary,   "nonfungal_reads_summary.csv",       row.names=FALSE)
write.csv(kingdom_breakdown,   "nonfungal_kingdom_breakdown.csv",   row.names=FALSE)
cat("\nNon-fungal reads summary:\n")
print(nonfungal_summary)

step_postplutof <- dplyr::bind_rows(lapply(names(alldat), function(nm) {
  summarise_ps(alldat[[nm]], nm, "post_PlutoF_filter")
}))
step_postplutof_full <- dplyr::bind_rows(lapply(names(alldat_full), function(nm) {
  summarise_ps(alldat_full[[nm]], paste0(nm, "_full"), "post_PlutoF_filter")
}))

# =============================================================================
# SECTION 5 — REMOVE LOW-DEPTH SAMPLES
# alldat/alldat_rg2_sensitivity use a real threshold (independent-sample
# methods need reliable per-sample depth). alldat_full uses a very low,
# effectively pass-through threshold — GLLVM/HMSC model PCR replicate as a
# random effect and don't need depth-based sample exclusion; this just drops
# empty/near-empty replicates so downstream code doesn't choke on them.
# =============================================================================

min_depth <- 10000
min_depth_full <- 1

cat("\nSample depth range before depth filter:\n")
for (nm in names(alldat)) {
  cat(sprintf("  %-8s min=%d  max=%d\n", nm,
              min(sample_sums(alldat[[nm]])), max(sample_sums(alldat[[nm]]))))
}

depth_filter <- function(ps, threshold) prune_samples(sample_sums(ps) >= threshold, ps)

alldat <- lapply(alldat, depth_filter, threshold = min_depth)
alldat <- lapply(alldat, function(ps) prune_taxa(taxa_sums(ps) > 0, ps))

alldat_rg2_sensitivity <- lapply(alldat_rg2_sensitivity, depth_filter, threshold = min_depth)
alldat_rg2_sensitivity <- lapply(alldat_rg2_sensitivity, function(ps) prune_taxa(taxa_sums(ps) > 0, ps))

alldat_full <- lapply(alldat_full, depth_filter, threshold = min_depth_full)
alldat_full <- lapply(alldat_full, function(ps) prune_taxa(taxa_sums(ps) > 0, ps))

cat("\nSamples remaining after depth filter (>= 10,000 reads):\n")
for (nm in names(alldat)) {
  cat(sprintf("  %-8s %d samples, %d taxa\n", nm, nsamples(alldat[[nm]]), ntaxa(alldat[[nm]])))
}
cat("\nSamples remaining after depth filter, rg2 sensitivity objects:\n")
for (nm in names(alldat_rg2_sensitivity)) {
  cat(sprintf("  %-8s %d samples, %d taxa\n", nm, nsamples(alldat_rg2_sensitivity[[nm]]), ntaxa(alldat_rg2_sensitivity[[nm]])))
}
cat(sprintf("\nSamples remaining after depth filter, full-complexity objects (>= %d reads):\n", min_depth_full))
for (nm in names(alldat_full)) {
  cat(sprintf("  %-8s %d samples, %d taxa\n", nm, nsamples(alldat_full[[nm]]), ntaxa(alldat_full[[nm]])))
}

step_postdepth <- dplyr::bind_rows(lapply(names(alldat), function(nm) {
  summarise_ps(alldat[[nm]], nm, "post_depth_filter")
}))
step_postdepth_full <- dplyr::bind_rows(lapply(names(alldat_full), function(nm) {
  summarise_ps(alldat_full[[nm]], paste0(nm, "_full"), "post_depth_filter")
}))

# =============================================================================
# SECTION 6 — RAREFACTION (even depth per strategy)
# =============================================================================

rarfun <- function(x) {
  rfy <- min(sample_sums(x))
  rarefy_even_depth(x, sample.size=rfy, replace=FALSE, rngseed=1)
}

alldat.rfy <- lapply(alldat, rarfun)
names(alldat.rfy) <- names(alldat)

cat("\nRarefaction depths:\n")
for (nm in names(alldat.rfy)) {
  cat(sprintf("  %-8s rarefied to %d reads/sample (%d samples, %d taxa)\n",
              nm,
              min(sample_sums(alldat.rfy[[nm]])),
              nsamples(alldat.rfy[[nm]]),
              ntaxa(alldat.rfy[[nm]])))
}

step_postrfy <- dplyr::bind_rows(lapply(names(alldat.rfy), function(nm) {
  summarise_ps(alldat.rfy[[nm]], nm, "post_rarefaction")
}))

# =============================================================================
# SECTION 7 — TAXONOMY TABLE CLEANING
# =============================================================================

# Standardize genuinely uninformative tax_table entries (empty string, bare
# rank prefixes like "g__", "unknown"/"NA" text) to real NA. Deliberately NOT
# using microViz::tax_fix() here — it always fills unknown/NA ranks by
# carrying forward the nearest classified parent rank + a suffix (e.g. an
# unassigned Species becomes "g__Thelebolus Genus"), which is misleading for
# this table. Unassigned ranks should just stay NA.
# Matches "", "g__", "unknown", "NA", "NaN", and compound forms like
# "g__unknown" / "g__NA" / "s__Unknown" (all case-insensitive).
tax_placeholder_pattern <- "^([a-zA-Z]__)?(unknown|na|nan)?$"

clean_tax_na <- function(ps) {
  tt <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
  tt[] <- lapply(tt, function(col) {
    col <- trimws(col)
    col[is.na(col) | grepl(tax_placeholder_pattern, col, ignore.case = TRUE)] <- NA_character_
    col
  })
  tax_table(ps) <- tax_table(as.matrix(tt))
  ps
}
label_dup_default <- function(x)
  label_duplicate_taxa(x, "Species", duplicate_label="<tax> <id>")

alldat <- lapply(alldat, clean_tax_na)
alldat.rfy <- lapply(alldat.rfy, clean_tax_na)
alldat_rg2_sensitivity <- lapply(alldat_rg2_sensitivity, clean_tax_na)
alldat_full <- lapply(alldat_full, clean_tax_na)

alldat <- lapply(alldat, label_dup_default)
alldat.rfy <- lapply(alldat.rfy, label_dup_default)
alldat_rg2_sensitivity <- lapply(alldat_rg2_sensitivity, label_dup_default)
alldat_full <- lapply(alldat_full, label_dup_default)

cat("\nTaxonomy cleaned. Example (nopool, first 5 taxa):\n")
print(head(tax_table(alldat$nopool), 5))

step_final <- dplyr::bind_rows(lapply(names(alldat), function(nm) {
  summarise_ps(alldat[[nm]], nm, "final_non_rarefied")
}))
step_final_rfy <- dplyr::bind_rows(lapply(names(alldat.rfy), function(nm) {
  summarise_ps(alldat.rfy[[nm]], nm, "final_rarefied")
}))
step_final_rg2 <- dplyr::bind_rows(lapply(names(alldat_rg2_sensitivity), function(nm) {
  summarise_ps(alldat_rg2_sensitivity[[nm]], paste0(nm, "_rg2"), "final_rg2_sensitivity")
}))
step_final_full <- dplyr::bind_rows(lapply(names(alldat_full), function(nm) {
  summarise_ps(alldat_full[[nm]], paste0(nm, "_full"), "final_full_complexity")
}))

# =============================================================================
# SECTION 8 — EXPORT COMPREHENSIVE COMPARISON TABLE
# Tracks n_samples, n_taxa, total_reads, min/max/median depth at each step.
# =============================================================================

comparison_tbl <- dplyr::bind_rows(
  step_initial,
  step_initial_full,
  step_postplutof,
  step_postplutof_full,
  step_postdepth,
  step_postdepth_full,
  step_final,
  step_final_rfy,
  step_final_rg2,
  step_final_full
)

write.csv(comparison_tbl, "phyloseq_comparison_summary.csv", row.names=FALSE)
cat("\nPhyloseq comparison summary (all steps x all strategies):\n")
print(comparison_tbl)

# =============================================================================
# SECTION 9 — SAVE R WORKSPACE
#
# Objects in eco_analysis.RData:
#   alldat                  — named list (nopool, pool, pspool): PCR replicates
#                             collapsed by summing, NO rg2 filter, non-fungal
#                             removed, depth-filtered, taxonomy-cleaned. NOT
#                             rarefied. CANONICAL object for PERMANOVA/dbRDA
#                             (robust Aitchison), iNEXT3D, NRI/NTI.
#   alldat.rfy              — same list but rarefied to min depth per strategy.
#                             Use for UniFrac (weighted/unweighted) and other
#                             fixed-depth alpha/beta diversity methods.
#   alldat_rg2_sensitivity  — same as alldat but from rg2-filtered (OTU in
#                             >=2 PCR replicates) input. Optional robustness
#                             check only — not the primary analysis path.
#   alldat_full             — named list (nopool, pool, pspool): full-
#                             complexity, one row per PCR replicate. Same
#                             non-fungal removal / PlutoF taxonomy upgrade /
#                             taxonomy cleaning as alldat, but depth-filtered
#                             with a very low pass-through threshold
#                             (min_depth_full) instead of the real 10,000-read
#                             cutoff, and NOT rarefied. For GLLVM/HMSC, which
#                             model PCR replicate (nested in root_id) as a
#                             random effect and use a read-count offset
#                             instead of rarefaction.
#
# Primary analysis strategy: alldat$nopool / alldat.rfy$nopool / alldat_full$nopool
# Pool and pspool kept for sensitivity comparisons.
# =============================================================================

save.image(file="eco_analysis.RData")
cat("\nWorkspace saved to eco_analysis.RData\n")
cat("Load with: load('eco_analysis.RData')\n")
cat("\nNext step: 5_community_composition.R\n")
