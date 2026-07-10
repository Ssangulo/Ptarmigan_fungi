# =============================================================================
# 9_dark_taxa_SH_matching.R
# E1 (dark taxa) exploration: assign UNITE Species Hypothesis (SH) identifiers
# to the retained fungal OTUs of alldat[[1]] via the authoritative TU-NHM
# sh_matching pipeline (the service behind PlutoF / unite.ut.ee), run locally
# in a Singularity container.
#
# Motivation. DADA2 + UNITE naive-Bayes (minBoot=80) leaves a large share of
# OTUs taxonomically unresolved ("dark taxa"; operationally is.na(Genus)). The
# co-author data-check report flags that WINTER samples are dominated by these
# dark taxa -- winter's most persistent OTU, OTU1369, is unclassified below
# Class (Dothideomycetes). This script gives every retained fungal OTU a stable
# SH handle + best-hit lineage + % similarity + a present_in (matches an
# existing SH) vs new-SH (novel) status, so we can state how much of the dark
# winter community is *placeable* to a known SH versus *genuinely novel*.
#
# NB. The pipeline's per-OTU consensus file (matches_out_taxonomy.csv) is the
# same format already produced for non-fungal removal in 4_data_prep.R
# (matches_out_taxonomy_nopool.csv) -- but that step only used common_taxonomy
# and only on the unresolved subset. Here we run SHmatching on the whole
# retained fungal community AND capture the SH codes / new-SH status that step
# discarded.
#
# Requires: Scripts/9_setup_sh_matching.sh run once (builds the container +
#           UNITE SH reference UDBs under /home/daniel/Ptarmigan/sh_matching).
# Run:      conda run -n r_env Rscript Scripts/9_dark_taxa_SH_matching.R
#
# Reference: sh_matching_pub v2.0.4 (github.com/TU-NHM/sh_matching_pub).
# Cite: Abarenkov, Koljalg & Nilsson (2022) BISS 6: e93856.
# =============================================================================

suppressMessages({
  library(phyloseq)
  library(Biostrings)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# ---- Config -----------------------------------------------------------------
sh_base   <- "/home/daniel/Ptarmigan/sh_matching"          # run cwd for the pipeline
sif       <- file.path(sh_base, "sh_matching.sif")
run_id    <- "1"                                            # must be numeric (pipeline asserts)
region    <- "its2"                                         # amplicons are ITS2-only
itsx_step <- "yes"                                          # QC via ITSx; see exclusion check below
force_rerun <- FALSE          # TRUE -> re-run the (slow) container even if outputs exist
primary_thr <- "02"          # primary SH distance threshold: 2.0% (also report 1.5%/3.0%)
report_thrs <- c("015", primary_thr, "03")

data_dir  <- "/home/daniel/Ptarmigan/trimmed/mergedPlates/"
out_dir   <- "/home/daniel/Ptarmigan/models/"
plot_dir  <- "/home/daniel/Ptarmigan/plots/"
dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Headless-server plotting guard (see 4/5 DECISIONS): standing null device so
# no base-graphics call tries to open X11; write PNGs via explicit png().
grDevices::pdf(NULL)
save_png <- function(plot_obj, filename, width, height, dpi = 300, units = "in") {
  png(file.path(plot_dir, filename), width = width, height = height, units = units, res = dpi)
  print(plot_obj); dev.off()
}

# =============================================================================
# SECTION 1 — LOAD alldat[[1]] AND EXPORT REP. SEQUENCES
# =============================================================================
load(file.path(data_dir, "eco_analysis.RData"))
ps <- alldat[[1]]                                # nopool: PCR reps collapsed, fungal-only
stopifnot(!is.null(refseq(ps)))
cat(sprintf("alldat[[1]]: %d samples x %d OTUs\n", nsamples(ps), ntaxa(ps)))

tax <- as.data.frame(as(tax_table(ps), "matrix"), stringsAsFactors = FALSE)
tax$OTU_ID <- rownames(tax)
tax$is_dark <- is.na(tax$Genus)                  # E1 operational definition (is.na Genus)
cat(sprintf("Dark taxa (is.na Genus): %d / %d (%.1f%%)\n",
            sum(tax$is_dark), nrow(tax), 100 * mean(tax$is_dark)))

# Write query fasta: indata/source_<run_id>, headers = OTU IDs
indata_fa <- file.path(sh_base, "indata", paste0("source_", run_id))
writeXStringSet(refseq(ps), filepath = indata_fa, format = "fasta")
cat("Wrote query fasta:", indata_fa, "(", length(refseq(ps)), "seqs )\n")

# =============================================================================
# SECTION 2 — RUN THE SH MATCHING CONTAINER
# =============================================================================
out_zip <- file.path(sh_base, "outdata", paste0("source_", run_id, ".zip"))
if (force_rerun || !file.exists(out_zip)) {
  stopifnot("container not built -- run Scripts/9_setup_sh_matching.sh first" =
              file.exists(sif))
  # run_pipeline.sh positional args: run_id region itsx remove_userdir vsearch usearch05
  cmd <- sprintf(
    "cd %s && singularity run %s /sh_matching/run_pipeline.sh %s %s %s no no no",
    shQuote(sh_base), shQuote(sif), run_id, region, itsx_step)
  cat("Running SH matching (this can take a while)...\n  ", cmd, "\n")
  log_path <- file.path(sh_base, sprintf("run_sh_matching_%s.log", run_id))
  status <- system2("bash", c("-c", shQuote(cmd)),
                    stdout = log_path, stderr = log_path)
  cat("Container exit status:", status, "-- log:", log_path, "\n")
  if (status != 0 || !file.exists(out_zip))
    stop("SH matching run failed; inspect ", log_path)
} else {
  cat("Reusing existing SH matching output:", out_zip,
      "\n(set force_rerun <- TRUE to regenerate)\n")
}

# =============================================================================
# SECTION 3 — PARSE SH MATCHING OUTPUT
# =============================================================================
parse_dir <- file.path(sh_base, "outdata", paste0("parsed_", run_id))
unlink(parse_dir, recursive = TRUE); dir.create(parse_dir, recursive = TRUE)
utils::unzip(out_zip, exdir = parse_dir)

# 3a) Per-OTU consensus taxonomy across thresholds (seq_name = OTU ID)
tax_csv <- file.path(parse_dir, "matches", "matches_out_taxonomy.csv")
sh_tax <- read.delim(tax_csv, stringsAsFactors = FALSE, check.names = FALSE)
# columns: seq_id, seq_name, common_name_selection_status, common_taxonomy,
#          common_rank, matched_sequence, similarity_percentage
names(sh_tax)[names(sh_tax) == "seq_name"] <- "OTU_ID"
sh_tax$similarity_percentage <- suppressWarnings(as.numeric(sh_tax$similarity_percentage))
sh_tax$common_rank <- suppressWarnings(as.numeric(sh_tax$common_rank))

# 3b) SH code + status per threshold. Existing-SH hits live in
# matches_out_<thr>.csv (status "present_in"); new SHs in matches_1_out_<thr>.csv.
# Both share a fixed 9-col tab schema -> parse positionally.
read_threshold <- function(thr) {
  f_hit <- file.path(parse_dir, "matches", sprintf("matches_out_%s.csv", thr))
  f_new <- file.path(parse_dir, "matches", sprintf("matches_1_out_%s.csv", thr))
  cols  <- c("seq_id_tmp","OTU_ID","status","SH_code","SH_taxonomy",
             "compound_cl_code","compound_taxonomy","dup_id","dup_accno")
  grab <- function(f) {
    if (!file.exists(f) || length(readLines(f, n = 2)) < 2) return(NULL)
    # Data rows carry a trailing empty tab-field (one more column than the
    # header). read.delim would silently take col 1 as row.names and shift
    # every column left by one -- so read headerless with padded names for any
    # trailing field(s), then keep the 9 real columns.
    d <- read.table(f, header = FALSE, sep = "\t", quote = "", skip = 1,
                    fill = TRUE, comment.char = "", stringsAsFactors = FALSE,
                    col.names = c(cols, "extra1", "extra2"))
    if (nrow(d) == 0) return(NULL)
    d[, cols]
  }
  out <- dplyr::bind_rows(grab(f_hit), grab(f_new))
  if (is.null(out) || nrow(out) == 0) return(NULL)
  out <- out[, c("OTU_ID","status","SH_code","SH_taxonomy")]
  names(out)[-1] <- paste0(names(out)[-1], "_", thr)
  out
}
sh_thr <- Reduce(function(a, b) full_join(a, b, by = "OTU_ID"),
                 Filter(Negate(is.null), lapply(report_thrs, read_threshold)))

# 3c) ID translation table (OTU_ID <-> internal code) + ITSx-excluded seqs.
# excluded_<run_id>.txt lists INTERNAL codes (i1_Ni), so translate to OTU IDs.
names_path <- file.path(parse_dir, paste0("source_", run_id, "_names"))
id_map   <- read.delim(names_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
code2otu <- setNames(id_map$V1, id_map$V2)             # V1=OTU_ID, V2=internal code
excl_path <- file.path(parse_dir, paste0("excluded_", run_id, ".txt"))
excl_otus <- character(0)
if (file.exists(excl_path)) {
  ex <- read.delim(excl_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  excl_otus <- unique(na.omit(code2otu[ex[[1]]]))
}
cat(sprintf("ITSx-excluded (chimeric/non-ITS) OTUs: %d\n", length(excl_otus)))

# ---- Classify each OTU ------------------------------------------------------
#   placeable_existing_SH : matched a named UNITE SH (status present_in)
#   novel_new_SH          : formed a new SH / singleton within a known compound
#   no_SH_match           : no >=80% match to any UNITE SH ref (deepest dark)
#   excluded_nonITS       : dropped by ITSx as chimeric/broken/non-ITS
status_col <- paste0("status_", primary_thr)
sh <- tax %>%
  left_join(sh_tax[, c("OTU_ID","common_taxonomy","common_rank",
                       "similarity_percentage")], by = "OTU_ID") %>%
  left_join(sh_thr, by = "OTU_ID")
sh$SH_class <- dplyr::case_when(
  sh$OTU_ID %in% excl_otus                               ~ "excluded_nonITS",
  grepl("present", sh[[status_col]], ignore.case = TRUE) ~ "placeable_existing_SH",
  !is.na(sh[[status_col]])                               ~ "novel_new_SH",
  TRUE                                                   ~ "no_SH_match")
cat("\n-- SH classification (all retained OTUs) --\n"); print(table(sh$SH_class))
cat("\n-- SH classification for DARK taxa (is.na Genus) --\n")
print(table(sh$SH_class[sh$is_dark]))

# =============================================================================
# SECTION 4 — SEASONAL OCCUPANCY & READ SHARE (for the winter caveat)
# =============================================================================
otu <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) otu <- t(otu)             # -> samples x taxa
sd  <- as(sample_data(ps), "data.frame")
stopifnot(identical(rownames(otu), rownames(sd)))
rel <- sweep(otu, 1, rowSums(otu), "/")          # per-sample relative abundance

season_stats <- function(season) {
  idx <- which(sd$Season == season)
  data.frame(
    OTU_ID   = colnames(otu),
    occ      = colSums(otu[idx, , drop = FALSE] > 0),          # n samples present
    n_season = length(idx),
    mean_share = colMeans(rel[idx, , drop = FALSE]),           # mean rel. abundance
    stringsAsFactors = FALSE)
}
w <- season_stats("winter"); s <- season_stats("summer")
occ <- data.frame(
  OTU_ID = w$OTU_ID,
  winter_occ = w$occ, winter_n = w$n_season, winter_mean_share = w$mean_share,
  summer_occ = s$occ, summer_n = s$n_season, summer_mean_share = s$mean_share)

sh <- left_join(sh, occ, by = "OTU_ID")

# Sanity check vs the data-check report's OTU1369 figures (~20/26 winter, 11.2%)
o1369 <- sh[sh$OTU_ID == "OTU1369", ]
if (nrow(o1369) == 1)
  cat(sprintf("\nOTU1369 check: winter %d/%d (%.1f%% mean share), summer %d/%d (%.1f%%); SH=%s [%s], %.1f%% id\n",
              o1369$winter_occ, o1369$winter_n, 100*o1369$winter_mean_share,
              o1369$summer_occ, o1369$summer_n, 100*o1369$summer_mean_share,
              o1369[[paste0("SH_code_", primary_thr)]], o1369$SH_class,
              o1369$similarity_percentage))

# =============================================================================
# SECTION 5 — E1 SUMMARIES & THE WINTER-CAVEAT ANSWER
# =============================================================================
# 5a) SH-assignment rates among dark taxa
dark <- sh[sh$is_dark, ]
assign_summary <- dark %>%
  summarise(
    n_dark               = n(),
    placeable_existing_SH = sum(SH_class == "placeable_existing_SH"),
    novel_new_SH          = sum(SH_class == "novel_new_SH"),
    excluded_nonITS       = sum(SH_class == "excluded_nonITS"),
    no_SH_match           = sum(SH_class == "no_SH_match"),
    gained_finer_lineage  = sum(!is.na(common_rank) & common_rank > 5, na.rm = TRUE),
    median_similarity     = median(similarity_percentage, na.rm = TRUE))
cat("\n-- E1: dark-taxa SH assignment summary --\n"); print(assign_summary)
write.csv(assign_summary, file.path(out_dir, "E1_SH_assignment_summary.csv"),
          row.names = FALSE)

# 5b) Winter-dominant dark taxa (ranked by winter mean share), with SH info
winter_dark <- dark %>%
  arrange(desc(winter_mean_share)) %>%
  transmute(OTU_ID, Class, Order,
            SH_code = .data[[paste0("SH_code_", primary_thr)]],
            SH_class, SH_common_taxonomy = common_taxonomy,
            similarity_percentage,
            winter_occ, winter_n, winter_mean_share,
            summer_occ, summer_n, summer_mean_share) %>%
  head(30)
cat("\n-- E1: top winter-dominant dark taxa --\n"); print(head(winter_dark, 12))
write.csv(winter_dark, file.path(out_dir, "E1_winter_dominant_dark_taxa_SH.csv"),
          row.names = FALSE)

# 5c) Full per-OTU SH table (the deliverable)
sh_out <- sh %>%
  select(OTU_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species, is_dark,
         SH_class, common_taxonomy, common_rank, similarity_percentage,
         dplyr::starts_with("SH_code_"), dplyr::starts_with("status_"),
         winter_occ, winter_mean_share, summer_occ, summer_mean_share)
write.csv(sh_out, file.path(out_dir, "dark_taxa_SH_matching.csv"), row.names = FALSE)
saveRDS(sh_out, file.path(data_dir, "dark_taxa_SH.rds"))

# 5d) The winter-caveat answer: per winter sample, what fraction of dark-taxon
#     reads is SH-placeable vs novel vs unresolved.
dark_ids   <- sh$OTU_ID[sh$is_dark]
place_ids  <- sh$OTU_ID[sh$is_dark & sh$SH_class == "placeable_existing_SH"]
winter_ids <- rownames(otu)[sd$Season == "winter"]
frac_tbl <- do.call(rbind, lapply(winter_ids, function(s) {
  tot  <- sum(otu[s, ]); dk <- sum(otu[s, dark_ids])
  plc  <- sum(otu[s, place_ids])
  data.frame(sample = s,
             dark_read_frac      = dk / tot,
             dark_placeable_frac = if (dk > 0) plc / dk else NA_real_)
}))
cat(sprintf("\n-- Winter caveat: dark reads = %.1f%% of winter reads on average;\n   of those, %.1f%% are placeable to an existing UNITE SH --\n",
            100 * mean(frac_tbl$dark_read_frac),
            100 * mean(frac_tbl$dark_placeable_frac, na.rm = TRUE)))
write.csv(frac_tbl, file.path(out_dir, "E1_winter_dark_placeable_fraction.csv"),
          row.names = FALSE)

# =============================================================================
# SECTION 6 — PLOTS
# =============================================================================
# 6a) SH % similarity: dark vs named OTUs
p_sim <- ggplot(sh[!is.na(sh$similarity_percentage), ],
                aes(x = similarity_percentage, fill = is_dark)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.6) +
  scale_fill_manual(values = c(`FALSE` = "grey50", `TRUE` = "firebrick"),
                    labels = c("named (Genus known)", "dark (is.na Genus)"),
                    name = NULL) +
  labs(x = "Best-hit similarity to UNITE SH (%)", y = "OTUs",
       title = "SH-match similarity: dark vs named OTUs") +
  theme_bw()
save_png(p_sim, "E1_SH_similarity_hist.png", width = 7, height = 4.5)

# 6b) Dark-taxa SH class breakdown
p_cls <- ggplot(dark, aes(x = SH_class, fill = SH_class)) +
  geom_bar() + coord_flip() + guides(fill = "none") +
  labs(x = NULL, y = "dark OTUs", title = "E1 dark taxa: placeable vs novel SH") +
  theme_bw()
save_png(p_cls, "E1_dark_taxa_novel_vs_placeable.png", width = 7, height = 4)

cat("\nDone. Tables ->", out_dir, "| plots ->", plot_dir, "\n")
