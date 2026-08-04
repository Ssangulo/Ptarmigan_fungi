# =============================================================================
# 8_functional_guilds.R
# Functional guild assignment (FUNGuild) and seasonal guild read-share analysis
#
# Adapted from the reference pipeline (Angulo Serrano et al., commit 65fbffa).
# The reference version summarised guilds by `habitat` on a root/soil design;
# this study's design is Season (winter/summer) x Year (2022-2024), so all
# habitat/site/elevation groupings have been replaced, and the ericoid-first
# guild hierarchy (2 OTUs here) replaced by the coprophily-first rule shared
# with 6_diversity_analyses.R / 10_hmsc.R.
#
# Workflow:
#   1. Prepare OTU + taxonomy table in R  ->  funguild_input.txt   [gated, see below]
#   2. Run FUNGuild CLI (bash, Python 3.8) ->  *.guilds.txt        [cached]
#   3. Import assignments as an OTU -> guild LOOKUP and build funguild_otu
#   4. Classify OTUs into display guilds and analyse seasonal read shares
#   5. Main figure (3 panels) + tables, staged into Supplementary/
#
# FUNGuild result (cached run): 1144 OTUs in, 700 assigned (61.2%)
# Reference: Nguyen et al. 2016 Methods Ecol & Evol
#
# Requires: alldat, alldat_full from 4_data_prep.R / eco_analysis.RData
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(fungaltraits)
})

setwd("/home/daniel/Ptarmigan/trimmed/mergedPlates/")
loaded_objs <- load("eco_analysis.RData")

# ---- Canonical output dirs (absolute; see CLAUDE.md plotting convention) -----
out_dir  <- "/home/daniel/Ptarmigan/models"     # tables (staging source)
plot_dir <- "/home/daniel/Ptarmigan/plots"
tab_dir  <- "/home/daniel/Ptarmigan/tables"     # CSVs mirrored here
supp_fig <- "/home/daniel/Ptarmigan/Scripts_server/Supplementary/figures"
supp_tab <- "/home/daniel/Ptarmigan/Scripts_server/Supplementary/tables"
for (d in c(out_dir, plot_dir, tab_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Headless plotting: a standing null device removes R's X11 fallback for ANY
# base-graphics call. ggplot figures go through save_png()/save_pdf().
grDevices::pdf(NULL)
save_png <- function(path, plot_obj, width = 8, height = 8, res = 300) {
  grDevices::png(path, width = width, height = height, units = "in", res = res)
  on.exit(grDevices::dev.off())
  print(plot_obj)
}
save_pdf <- function(path, plot_obj, width = 8, height = 8) {
  dev_fun <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  dev_fun(path, width = width, height = height)
  on.exit(grDevices::dev.off())
  print(plot_obj)
}
# write a CSV to both models/ (staging source) and tables/
write_tab <- function(df, name, row.names = FALSE) {
  write.csv(df, file.path(out_dir, name), row.names = row.names)
  write.csv(df, file.path(tab_dir, name), row.names = row.names)
}

# ---- Shared helper (from 6_diversity_analyses.R / 10_hmsc.R) -----------------
otu_mat_of <- function(ps) {
  m <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) m <- t(m)
  m                                              # samples x taxa
}

# =============================================================================
# BASH SETUP (run once outside R)
# =============================================================================
# conda create -n funguild python=3.8
# conda activate funguild
# conda install pandas
# git clone https://github.com/UMNFuN/FUNGuild
# cd FUNGuild/

# =============================================================================
# SECTION 1 -- PREPARE OTU + TAXONOMY TABLE FOR FUNGUILD
# =============================================================================
# Provenance record of how the cached funguild_input.txt / *.guilds.txt were
# produced. GATED OFF by default: the cached .guilds.txt is the source of the
# funguild_otu that Sections 6/7/10 of the appendix already depend on, so a
# casual re-run must not silently regenerate an input that no longer matches it.
# Set REGEN_FUNGUILD_INPUT <- TRUE (and re-run the CLI below) only when
# deliberately refreshing the guild annotation.

REGEN_FUNGUILD_INPUT <- FALSE

ps_fg          <- alldat_full[[1]]                 # the object the cache was built from
otu_mat        <- as(otu_table(ps_fg), "matrix")
otu_table_data <- if (taxa_are_rows(ps_fg)) as.data.frame(otu_mat) else as.data.frame(t(otu_mat))

tax_mat        <- as(tax_table(ps_fg), "matrix")
tax_table_data <- as.data.frame(tax_mat, stringsAsFactors = FALSE)
tax_table_data <- tax_table_data[rownames(otu_table_data), , drop = FALSE]
stopifnot(identical(rownames(otu_table_data), rownames(tax_table_data)))

# ---- Clean taxonomy ---------------------------------------------------------
tax_table_data[is.na(tax_table_data)] <- ""
tax_table_data[] <- lapply(tax_table_data, function(x) {
  # No-op on this study's tax_table (that suffix is a microViz tax_fix artefact,
  # and script 4 deliberately refuses tax_fix) -- retained so this block still
  # reproduces the cached funguild_input.txt exactly as it was generated.
  x <- gsub("\\s+(Kingdom|Phylum|Class|Order|Family|Genus|Species)\\b", "", x)
  # Required: strips script 4's label_duplicate_taxa() suffix, e.g. "s__arescens 1".
  x <- gsub("\\s+\\d+$", "", x)
  trimws(x)
})

# ---- Ensure each rank only contains entries with the correct prefix ----------
rank_prefix <- c(Kingdom = "k__", Phylum = "p__", Class = "c__",
                 Order   = "o__", Family = "f__", Genus = "g__", Species = "s__")
missing <- setdiff(names(rank_prefix), colnames(tax_table_data))
if (length(missing) > 0)
  stop("Taxonomy columns missing: ", paste(missing, collapse = ", "))
for (r in names(rank_prefix)) {
  pref <- rank_prefix[[r]]
  tax_table_data[[r]][tax_table_data[[r]] != "" &
                        !startsWith(tax_table_data[[r]], pref)] <- ""
}

# ---- Remove ranks that repeat the parent rank -------------------------------
ranks <- names(rank_prefix)
for (i in 2:length(ranks)) {
  prev <- ranks[i - 1]; cur <- ranks[i]
  tax_table_data[[cur]][tax_table_data[[cur]] != "" &
                          tax_table_data[[cur]] == tax_table_data[[prev]]] <- ""
}

# ---- Collapse to semicolon-separated taxonomy string ------------------------
tax_table_data$taxonomy_string <- apply(tax_table_data, 1, function(x) {
  paste(x[x != ""], collapse = ";")
})

if (REGEN_FUNGUILD_INPUT) {
  otu_tax_table <- cbind(otu_table_data, taxonomy = tax_table_data$taxonomy_string)
  out_file <- "/home/daniel/Ptarmigan/trimmed/mergedPlates/funguild_input.txt"
  write.table(otu_tax_table, file = out_file, sep = "\t", quote = FALSE, col.names = NA)
  cat("Wrote:", out_file, "-- ", nrow(otu_tax_table), "taxa,",
      ncol(otu_tax_table) - 1, "samples\n")
  cat("NOTE: re-run the FUNGuild CLI (below) before Section 2, or the cached\n",
      "     .guilds.txt will no longer correspond to this input.\n")
} else {
  cat("Section 1: FUNGuild input regeneration SKIPPED (REGEN_FUNGUILD_INPUT = FALSE);\n",
      "           using the cached .guilds.txt as an OTU -> guild lookup.\n")
}

# =============================================================================
# BASH -- RUN FUNGUILD CLI (run outside R, only if REGEN_FUNGUILD_INPUT = TRUE)
# =============================================================================
# conda activate funguild
# cd FUNGuild/
# python Guilds_v1.1.py \
#   -otu /home/daniel/Ptarmigan/trimmed/mergedPlates/funguild_input.txt \
#   -db fungi
#
# Result: .guilds.txt written to same directory as input
# Cached run: 1144 OTUs in, 1086 matching taxonomy records, 700 assigned (61.2%)

# =============================================================================
# SECTION 2 -- IMPORT FUNGUILD RESULTS AS AN OTU -> GUILD LOOKUP
# =============================================================================
# IMPORTANT: only the OTU_ID -> guild columns of the .guilds.txt are used. Its
# own sample columns are NOT read (the cached file carries a superseded 114-
# sample set from before the depth-filter fix); every abundance in this script
# comes from the canonical phyloseq object instead. That makes the cached file's
# sample-set staleness irrelevant while keeping guild calls traceable to the
# exact run behind the funguild_otu used in appendix Sections 6/7/10.

funguild_data <- read.table(
  "/home/daniel/Ptarmigan/trimmed/mergedPlates/funguild_input.guilds.txt",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
)
names(funguild_data) <- gsub(" ", "_", names(funguild_data))
names(funguild_data)[1] <- "OTU_ID"

fg_assigned <- funguild_data[
  !is.na(funguild_data$Confidence_Ranking) &
    funguild_data$Confidence_Ranking != "-", ]
cat(sprintf("FUNGuild lookup: %d of %d OTUs guild-assigned (%.1f%%)\n",
            nrow(fg_assigned), nrow(funguild_data),
            100 * nrow(fg_assigned) / nrow(funguild_data)))

# =============================================================================
# SECTION 2b -- PER-OTU GUILD LOOKUP TABLE (consumed by scripts 6/7/10)
# =============================================================================
funguild_otu_rebuilt <- fg_assigned %>%
  transmute(
    OTU_ID             = OTU_ID,
    taxonomy_string    = taxonomy,
    Taxon              = Taxon,
    Taxon_Level        = Taxon_Level,
    Trophic_Mode       = Trophic_Mode,
    Guild              = Guild,
    Growth_Morphology  = Growth_Morphology,
    Trait              = Trait,
    Confidence_Ranking = Confidence_Ranking,
    Notes              = Notes,
    Citation_Source    = `Citation/Source`
  )
cat("funguild_otu (rebuilt):", nrow(funguild_otu_rebuilt), "guild-assigned OTUs\n")

# =============================================================================
# SECTION 2c -- FUNGALTRAITS GENUS-LEVEL LOOKUP (sensitivity check only)
# =============================================================================
# Kept as a documented cross-check, NOT as an analysis source: fungaltraits_otu
# is a strict subset of funguild_otu (195 vs 700 OTUs) and agrees with it on
# only 37% of Guild strings, because its genus-level "first non-NA" collapsing
# loses the species-level specificity of FUNGuild's per-OTU CLI run.
# Reference: Flores-Moreno et al. 2019 (fungaltraits / "funfun").
#
#   devtools::install_github("ropenscilabs/datastorr")
#   devtools::install_github("traitecoevo/fungaltraits")
# fungal_traits() downloads+caches funtothefun.csv via datastorr on first call.

first_non_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA else x[1]
}

ft <- fungal_traits()
genus_lookup <- ft %>%
  filter(!is.na(Genus), trimws(Genus) != "") %>%
  mutate(Genus_clean = trimws(tolower(Genus))) %>%
  group_by(Genus_clean) %>%
  summarise(
    Trophic_Mode_fg     = first_non_na(trophic_mode_fg),
    Guild_fg            = first_non_na(guild_fg),
    Growth_Form_fg      = first_non_na(growth_form_fg),
    Confidence_fg       = first_non_na(confidence_fg),
    Taxonomic_Level_fg  = first_non_na(taxonomic_level_fg),
    Source_fg           = first_non_na(source_funguild_fg),
    Notes_fg            = first_non_na(notes_fg),
    n_distinct_guild_fg = n_distinct(guild_fg[!is.na(guild_fg)]),
    .groups = "drop"
  )

otu_genus <- data.frame(
  OTU_ID          = rownames(tax_table_data),
  taxonomy_string = tax_table_data$taxonomy_string,
  Genus           = tax_table_data$Genus,
  stringsAsFactors = FALSE
)
otu_genus$Genus_clean <- trimws(tolower(sub("^g__", "", otu_genus$Genus)))

fungaltraits_otu_rebuilt <- otu_genus %>%
  filter(Genus_clean != "") %>%
  inner_join(genus_lookup %>% select(-n_distinct_guild_fg), by = "Genus_clean") %>%
  filter(!is.na(Guild_fg)) %>%
  transmute(
    OTU_ID, taxonomy_string, Genus,
    Trophic_Mode_fg, Guild_fg, Growth_Form_fg,
    Confidence_fg, Taxonomic_Level_fg, Source_fg, Notes_fg
  ) %>%
  as.data.frame()
cat("fungaltraits_otu (rebuilt):", nrow(fungaltraits_otu_rebuilt), "genus-assigned OTUs\n")

# =============================================================================
# SECTION 2d -- NON-DESTRUCTIVE PERSISTENCE INTO THE CANONICAL WORKSPACE
# =============================================================================
# Re-running 4_data_prep.R regenerates eco_analysis.RData from scratch and WIPES
# funguild_otu / fungaltraits_otu (see CLAUDE.md). This block is the restore
# path. It is deliberately conservative:
#   missing    -> rebuild and save (restore)
#   identical  -> skip the save entirely (no needless workspace rewrite)
#   different  -> DO NOT save; warn, dump the rebuild for inspection, and keep
#                 using the workspace copy that Sections 6/7/10 were built on.

check_guild_obj <- function(name, rebuilt) {
  if (!exists(name, inherits = TRUE)) {
    cat(sprintf("  %s: MISSING from workspace -> will restore from rebuild\n", name))
    return(list(action = "save", value = rebuilt))
  }
  current <- get(name, inherits = TRUE)
  if (isTRUE(all.equal(current, rebuilt, check.attributes = FALSE))) {
    cat(sprintf("  %s: unchanged (%d rows) -> no save needed\n", name, nrow(current)))
    return(list(action = "keep", value = current))
  }
  dump_path <- file.path(out_dir, paste0(name, "_rebuilt_DIFFERS.csv"))
  write.csv(rebuilt, dump_path, row.names = FALSE)
  warning(sprintf("%s rebuilt from the cached .guilds.txt DIFFERS from the workspace copy ",
                  name), "(workspace ", nrow(current), " rows vs rebuild ",
          nrow(rebuilt), " rows). NOT overwriting -- appendix Sections 6/7/10 were ",
          "built on the workspace copy. Rebuild dumped to ", dump_path,
          " for inspection.", call. = FALSE)
  cat(sprintf("  %s: DIFFERS -> kept workspace copy, rebuild dumped to %s\n",
              name, dump_path))
  list(action = "keep", value = current)
}

cat("Guild object integrity check:\n")
chk_fg <- check_guild_obj("funguild_otu",     funguild_otu_rebuilt)
chk_ft <- check_guild_obj("fungaltraits_otu", fungaltraits_otu_rebuilt)
funguild_otu     <- chk_fg$value
fungaltraits_otu <- chk_ft$value

if (any(c(chk_fg$action, chk_ft$action) == "save")) {
  save(list = unique(c(loaded_objs, "funguild_otu", "fungaltraits_otu")),
       file = "/home/daniel/Ptarmigan/trimmed/mergedPlates/eco_analysis.RData")
  cat("  -> restored missing guild object(s) into eco_analysis.RData\n")
} else {
  cat("  -> eco_analysis.RData left untouched\n")
}

# =============================================================================
# SECTION 3 -- GUILD CLASSIFICATION (canonical object)
# =============================================================================
# Backbone rule is taken VERBATIM from 10_hmsc.R:193-213 (itself from
# 6_diversity_analyses.R:882-899) so Section 9 and Section 10 of the appendix
# classify OTUs identically. Priority, first match wins:
#   dung_saprotroph > plant_associated > pathotroph > dark_unassigned > other
#
# Two deliberate differences from Section 10, both recorded in the crosswalk:
#   (1) NO prevalence filter -- all 1129 OTUs of alldat$nopool, because this is
#       a descriptive read-share layer, not a model design matrix.
#   (2) The `other` bin is SPLIT for display. In Section 10 `other` mixes
#       FUNGuild-assigned OTUs with 46 that have a genus but no guild call at
#       all; those carry ~10% of winter reads. Acceptable as a model trait,
#       misleading in a read-share figure.

ps_can  <- alldat$nopool                          # 52 samples, PCR reps collapsed
m       <- otu_mat_of(ps_can)                     # samples x taxa, raw counts
otu_ids <- colnames(m)
cat(sprintf("\nCanonical object: %d samples x %d OTUs\n", nrow(m), ncol(m)))

COPRO_GENERA <- c("Sordaria","Podospora","Cercophora","Chaetomium","Schizothecium",
                  "Preussia","Delitschia","Pilobolus","Ascobolus","Saccobolus",
                  "Sporormiella","Coprinopsis","Thelebolus","Coniochaeta")
strip_rank <- function(x) sub("^[a-z]__", "", x)   # "g__Sporormiella" -> "Sporormiella"

guild_str <- funguild_otu$Guild[match(otu_ids, funguild_otu$OTU_ID)]
tax_all   <- data.frame(as(tax_table(ps_can), "matrix"), stringsAsFactors = FALSE)
otu_gen   <- strip_rank(tax_all$Genus[match(otu_ids, rownames(tax_all))])
otu_gen   <- ifelse(otu_gen %in% c("", "NA"), NA, otu_gen)

cat(sprintf("Guild-file coverage of canonical OTUs: %d of %d (%.1f%%) have a FUNGuild call\n",
            sum(!is.na(guild_str)), length(otu_ids),
            100 * mean(!is.na(guild_str))))

is_copro <- (!is.na(guild_str) & grepl("Dung Saprotroph", guild_str)) |
            (!is.na(otu_gen)   & otu_gen %in% COPRO_GENERA)
is_plant <- !is_copro & !is.na(guild_str) &
            grepl("Endophyte|Plant Saprotroph|Plant Pathogen|Epiphyte", guild_str)
is_patho <- !is_copro & !is_plant & !is.na(guild_str) & grepl("Pathogen|Parasite", guild_str)
is_dark  <- !is_copro & !is_plant & !is_patho & is.na(otu_gen) & is.na(guild_str)
guild_class <- ifelse(is_copro, "dung_saprotroph",
               ifelse(is_plant, "plant_associated",
               ifelse(is_patho, "pathotroph",
               ifelse(is_dark,  "dark_unassigned", "other"))))

# ---- Display remap (splits `other` on whether a FUNGuild call exists) --------
GUILD_LEVELS <- c("Dung saprotroph", "Plant-associated", "Pathotroph",
                  "Other assigned", "Unassigned / dark")
guild_display <- ifelse(guild_class == "dung_saprotroph",  "Dung saprotroph",
                 ifelse(guild_class == "plant_associated", "Plant-associated",
                 ifelse(guild_class == "pathotroph",       "Pathotroph",
                 ifelse(guild_class == "other" & !is.na(guild_str), "Other assigned",
                        "Unassigned / dark"))))
guild_display <- factor(guild_display, levels = GUILD_LEVELS)
stopifnot(!anyNA(guild_display))

cat("\nSection 10 guild_class (all canonical OTUs):\n"); print(table(guild_class))
cat("\nSection 9 display guild:\n");                     print(table(guild_display))

# =============================================================================
# SECTION 4 -- SEASONAL READ-SHARE ANALYSIS
# =============================================================================
# Per-SAMPLE read shares (each sample contributes equally regardless of depth),
# then averaged within season -- not pooled reads, which would let the deepest
# samples dominate.

rel <- m / rowSums(m)
stopifnot(all(abs(rowSums(rel) - 1) < 1e-10))

share <- t(rowsum(t(rel), group = as.character(guild_display)))  # samples x guilds
share <- share[, GUILD_LEVELS, drop = FALSE]
stopifnot(all(abs(rowSums(share) - 1) < 1e-10))

md       <- data.frame(as(sample_data(ps_can), "data.frame"), stringsAsFactors = FALSE)
labs_map <- c(winter = "Winter", summer = "Summer")            # appendix convention
Season   <- factor(unname(labs_map[md$Season]), levels = c("Winter", "Summer"))
Year     <- factor(md$Year)
stopifnot(!anyNA(Season), !anyNA(Year), identical(rownames(share), rownames(m)))
cat("\nSeason x Year cell sizes:\n"); print(table(Season, Year))

# ---- Long format for plotting -----------------------------------------------
guild_long <- data.frame(
  SampleID = rep(rownames(share), times = ncol(share)),
  Season   = rep(Season,          times = ncol(share)),
  Year     = rep(Year,            times = ncol(share)),
  Guild    = factor(rep(colnames(share), each = nrow(share)), levels = GUILD_LEVELS),
  share    = as.vector(share),
  stringsAsFactors = FALSE
)

# ---- Winter vs summer per guild ---------------------------------------------
season_stats <- do.call(rbind, lapply(GUILD_LEVELS, function(g) {
  x  <- share[, g]
  xs <- x[Season == "Summer"]; xw <- x[Season == "Winter"]
  wt <- suppressWarnings(wilcox.test(xs, xw))     # ties -> normal approximation
  data.frame(
    Guild            = g,
    n_summer         = length(xs),               n_winter         = length(xw),
    mean_summer      = mean(xs),                 mean_winter      = mean(xw),
    sd_summer        = sd(xs),                   sd_winter        = sd(xw),
    median_summer    = median(xs),               median_winter    = median(xw),
    iqr_summer       = IQR(xs),                  iqr_winter       = IQR(xw),
    prevalence_summer = sum(xs > 0),             prevalence_winter = sum(xw > 0),
    W                = unname(wt$statistic),
    p_value          = wt$p.value,
    stringsAsFactors = FALSE
  )
}))
season_stats$q_value <- p.adjust(season_stats$p_value, method = "BH")

# ---- Robustness: dung saprotroph restricted to COPRO_GENERA genus matches ----
# Drops the OTUs that match only on the FUNGuild guild string, to show the
# headline is not an artefact of the string-matching half of the rule.
copro_genus_only <- !is.na(otu_gen) & otu_gen %in% COPRO_GENERA
copro_share      <- rowSums(rel[, copro_genus_only, drop = FALSE])
wt_copro <- suppressWarnings(wilcox.test(copro_share[Season == "Summer"],
                                         copro_share[Season == "Winter"]))
robust_row <- data.frame(
  Guild             = "Dung saprotroph (COPRO_GENERA genus only, robustness)",
  n_summer          = sum(Season == "Summer"),
  n_winter          = sum(Season == "Winter"),
  mean_summer       = mean(copro_share[Season == "Summer"]),
  mean_winter       = mean(copro_share[Season == "Winter"]),
  sd_summer         = sd(copro_share[Season == "Summer"]),
  sd_winter         = sd(copro_share[Season == "Winter"]),
  median_summer     = median(copro_share[Season == "Summer"]),
  median_winter     = median(copro_share[Season == "Winter"]),
  iqr_summer        = IQR(copro_share[Season == "Summer"]),
  iqr_winter        = IQR(copro_share[Season == "Winter"]),
  prevalence_summer = sum(copro_share[Season == "Summer"] > 0),
  prevalence_winter = sum(copro_share[Season == "Winter"] > 0),
  W                 = unname(wt_copro$statistic),
  p_value           = wt_copro$p.value,
  q_value           = NA_real_,                  # excluded from the BH family
  stringsAsFactors  = FALSE
)
season_stats_out <- rbind(season_stats, robust_row)
cat(sprintf("\nDung saprotroph OTUs: %d total (%d via COPRO_GENERA genus, %d via FUNGuild string only)\n",
            sum(is_copro), sum(copro_genus_only), sum(is_copro & !copro_genus_only)))

cat("\nMean per-sample read share by season:\n")
print(data.frame(Guild   = season_stats$Guild,
                 summer  = round(season_stats$mean_summer, 3),
                 winter  = round(season_stats$mean_winter, 3),
                 q       = signif(season_stats$q_value, 3)), row.names = FALSE)
cat(sprintf("\nReads on a guild-assigned OTU: summer %.0f%%, winter %.0f%%\n",
            100 * mean(rowSums(rel[, !is.na(guild_str), drop = FALSE])[Season == "Summer"]),
            100 * mean(rowSums(rel[, !is.na(guild_str), drop = FALSE])[Season == "Winter"])))

# ---- Season x Year consistency ----------------------------------------------
sy <- aggregate(share ~ Guild + Season + Year, data = guild_long,
                FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
season_year <- data.frame(sy[, c("Guild", "Season", "Year")],
                          mean_share = sy$share[, "mean"],
                          sd_share   = sy$share[, "sd"],
                          n_samples  = sy$share[, "n"])
season_year <- season_year[order(season_year$Guild, season_year$Year, season_year$Season), ]

# Per-year winter-vs-summer test on the dung-saprotroph share: does the guild
# shift repeat in all three years (the functional analogue of H1's null
# Season x Year interaction)?
per_year <- do.call(rbind, lapply(levels(Year), function(y) {
  i  <- Year == y
  xs <- share[i & Season == "Summer", "Dung saprotroph"]
  xw <- share[i & Season == "Winter", "Dung saprotroph"]
  wt <- suppressWarnings(wilcox.test(xs, xw))
  data.frame(Year = y, n_summer = length(xs), n_winter = length(xw),
             mean_summer = mean(xs), mean_winter = mean(xw),
             diff_summer_minus_winter = mean(xs) - mean(xw),
             p_value = wt$p.value, stringsAsFactors = FALSE)
}))
per_year$q_value <- p.adjust(per_year$p_value, method = "BH")
cat("\nDung-saprotroph share, winter vs summer within each year:\n")
print(per_year, row.names = FALSE, digits = 3)
cat(sprintf("Same direction (summer > winter) in %d of %d years\n",
            sum(per_year$diff_summer_minus_winter > 0), nrow(per_year)))

# =============================================================================
# SECTION 5 -- MAIN FIGURE (3 panels)
# =============================================================================
# Palette: Okabe-Ito (colourblind-safe). Dung saprotroph is #E69F00 rather than
# #D55E00 so that #D55E00 keeps reading uniquely as "Summer" in panel B, the
# only panel where Season is colour-mapped.
guild_cols <- c(
  "Dung saprotroph"   = "#E69F00",
  "Plant-associated"  = "#009E73",
  "Pathotroph"        = "#CC79A7",
  "Other assigned"    = "#56B4E9",
  "Unassigned / dark" = "#999999"   # grey reads correctly as "no information"
)
season_cols <- c(Winter = "#0072B2", Summer = "#D55E00")   # repo standard
base_theme  <- theme_classic(base_size = 12) +             # appendix main-figure convention
  theme(strip.background = element_blank(),
        strip.text       = element_text(face = "bold"),
        plot.tag         = element_text(face = "bold", size = 14))
wrap_lab <- c("Dung\nsaprotroph", "Plant-\nassociated", "Pathotroph",
              "Other\nassigned", "Unassigned\n/ dark")
names(wrap_lab) <- GUILD_LEVELS

# ---- Panel A: per-sample stacked composition --------------------------------
# Ordered within season by descending dung-saprotroph share, so the gradient of
# individual-dropping heterogeneity is readable.
ord <- rownames(share)[order(Season, -share[, "Dung saprotroph"])]
pA_dat <- guild_long
pA_dat$SampleID <- factor(pA_dat$SampleID, levels = ord)
n_lab <- setNames(paste0(levels(Season), " (n = ", table(Season), ")"), levels(Season))

pA <- ggplot(pA_dat, aes(SampleID, share, fill = Guild)) +
  geom_col(width = 1, position = position_stack(reverse = TRUE)) +
  facet_grid(~ Season, scales = "free_x", space = "free_x",
             labeller = labeller(Season = n_lab)) +
  scale_fill_manual(values = guild_cols, name = "Functional guild") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = "Dung sample (ordered by dung-saprotroph share)",
       y = "Read share") +
  base_theme +
  theme(axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
        panel.spacing = unit(0.8, "lines"))

# ---- Panel B: guild x season distributions + BH-q marks ---------------------
sig_lab <- function(q) ifelse(is.na(q), "",
                       ifelse(q < 0.001, "***",
                       ifelse(q < 0.01,  "**",
                       ifelse(q < 0.05,  "*", "ns"))))
# Marks sit on one constant row rather than tracking each guild's max, so they
# read as a comparable series instead of floating at five different heights.
ann <- data.frame(
  Guild = factor(GUILD_LEVELS, levels = GUILD_LEVELS),
  y     = 1.04,
  lab   = sig_lab(season_stats$q_value),
  stringsAsFactors = FALSE
)

pB <- ggplot(guild_long, aes(Guild, share, fill = Season)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, alpha = 0.85,
               position = position_dodge(width = 0.75)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
             size = 0.9, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = ann, aes(Guild, y, label = lab), inherit.aes = FALSE,
            size = 3.6, fontface = "bold") +
  scale_fill_manual(values = season_cols, name = "Season") +
  scale_x_discrete(labels = wrap_lab) +
  # Explicit breaks: the 1.09 headroom for the significance row otherwise makes
  # ggplot pick 0/30/60/90%, which clashes with panels A and C.
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.09),
                     breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = "Per-sample read share") +
  base_theme +
  theme(axis.text.x = element_text(size = 9))

# ---- Panel C: Season x Year consistency -------------------------------------
sy_n <- unique(season_year[, c("Season", "Year", "n_samples")])
sy_n <- aggregate(n_samples ~ Season + Year, data = sy_n, FUN = max)

pC <- ggplot(season_year, aes(Year, mean_share, fill = Guild)) +
  geom_col(width = 0.75, position = position_stack(reverse = TRUE)) +
  geom_text(data = sy_n, aes(Year, 1.04, label = paste0("n=", n_samples)),
            inherit.aes = FALSE, size = 2.9, colour = "grey30") +
  facet_wrap(~ Season) +
  scale_fill_manual(values = guild_cols, name = "Functional guild") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0),
                     limits = c(0, 1.08)) +
  labs(x = NULL, y = "Mean read share") +
  base_theme +
  theme(panel.spacing = unit(0.8, "lines"))

# ---- Assemble ---------------------------------------------------------------
fig1 <- (pA / (pB | pC)) +
  plot_layout(heights = c(1, 1.15), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  # legend.box = "vertical": the guild (5 keys) and Season (2 keys) legends do
  # not fit on one 10-inch row -- side by side, "Summer" gets clipped.
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.box.just = "left", legend.margin = margin(t = 0, b = 0),
        legend.key.size = unit(0.5, "cm"), legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))

fig_png <- file.path(plot_dir, "Fig1_guild_season.png")
fig_pdf <- file.path(plot_dir, "Fig1_guild_season.pdf")
save_png(fig_png, fig1, width = 10, height = 9, res = 600)
save_pdf(fig_pdf, fig1, width = 10, height = 9)
cat("\nWrote figure:\n  ", fig_png, "\n  ", fig_pdf, "\n")

# =============================================================================
# SECTION 6 -- TABLES
# =============================================================================

# (a) Guild x season read share + tests
write_tab(season_stats_out, "guild_season_readshare.csv")

# (b) Guild x season x year
write_tab(season_year, "guild_season_year_readshare.csv")

# (c) Per-OTU assignment (all canonical OTUs). NOTE: Trophic_Mode carries
#     leading-whitespace duplicates in the FUNGuild output (12 raw levels
#     collapsing to 8) -- trimws() before any tabulation or e.g. Saprotroph
#     splits into 235 + 6.
troph <- trimws(funguild_otu$Trophic_Mode[match(otu_ids, funguild_otu$OTU_ID)])
conf  <- funguild_otu$Confidence_Ranking[match(otu_ids, funguild_otu$OTU_ID)]
otu_assign <- data.frame(
  OTU_ID             = otu_ids,
  Genus              = otu_gen,
  FUNGuild_Guild     = guild_str,
  Trophic_Mode       = troph,
  Confidence_Ranking = conf,
  guild_class        = guild_class,             # Section 10 (5-level model trait)
  guild_display      = as.character(guild_display),  # Section 9 (5-level display)
  total_reads        = colSums(m),
  prevalence         = colSums(m > 0),
  stringsAsFactors   = FALSE
)
otu_assign <- otu_assign[order(-otu_assign$total_reads), ]
write_tab(otu_assign, "guild_otu_assignment.csv")
cat("\nTrophic modes (trimmed) across canonical OTUs:\n")
print(sort(table(troph, useNA = "no"), decreasing = TRUE))

# (d) Per-sample composition behind panel A
per_sample <- data.frame(SampleID = rownames(share), Season = Season, Year = Year,
                         as.data.frame(share), check.names = FALSE)
write_tab(per_sample, "guild_composition_per_sample.csv")

# (e) Crosswalk so Section 9 and Table S-HMSCa reconcile
crosswalk <- as.data.frame(table(guild_class = guild_class,
                                 guild_display = as.character(guild_display)))
crosswalk <- crosswalk[crosswalk$Freq > 0, ]
names(crosswalk)[3] <- "n_OTUs"
crosswalk$note <- ifelse(
  crosswalk$guild_class == "other" & crosswalk$guild_display == "Unassigned / dark",
  "Section 10 'other' OTUs with a genus but NO FUNGuild call -- split out for display",
  ifelse(crosswalk$guild_class == "other",
         "Section 10 'other' OTUs that DO carry a FUNGuild call", ""))
write_tab(crosswalk, "guild_scheme_crosswalk.csv")
cat("\nSection 10 -> Section 9 crosswalk:\n"); print(crosswalk, row.names = FALSE)
stopifnot(sum(crosswalk$n_OTUs) == length(otu_ids))

# =============================================================================
# SECTION 7 -- STAGE FIGURES AND TABLES INTO Supplementary/
# =============================================================================
figs <- c("Fig1_guild_season.png", "Fig1_guild_season.pdf")
tabs <- c("guild_season_readshare.csv", "guild_season_year_readshare.csv",
          "guild_otu_assignment.csv", "guild_composition_per_sample.csv",
          "guild_scheme_crosswalk.csv")
invisible(file.copy(file.path(plot_dir, figs), supp_fig, overwrite = TRUE))
invisible(file.copy(file.path(out_dir,  tabs), supp_tab, overwrite = TRUE))
cat(sprintf("\nStaged %d figures and %d tables into Supplementary/figures|tables\n",
            length(figs), length(tabs)))
cat("Script 8 complete.\n")
