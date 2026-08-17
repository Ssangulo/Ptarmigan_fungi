# =============================================================================
# 10_hmsc.R
# Hierarchical Modelling of Species Communities (HMSC) -- a joint species
# distribution model (JSDM) for the ptarmigan dung mycobiome.
#
# Purpose: complement the Section-6 GLLVM (H1, per-OTU Season differential
# abundance) with a full JSDM that adds, in one coherent model:
#   (1) variance partitioning across Season/Year and the random structure
#       (biological sample / PCR replicate; the bird level was dropped -- see
#       the note in Section 1);
#   (2) a trait -> response test -- does functional GUILD predict how OTUs
#       respond to Season (the Gamma matrix);
#   (3) a per-OTU Season niche (Beta) cross-checked against the GLLVM;
#   (4) a residual co-occurrence network (Omega); and
#   (5) explanatory vs 2-fold cross-validated predictive R2.
# A secondary phylogeny-augmented variant reports rho (phylogenetic signal in
# seasonal niches) -- expected near 0 for a transient, propagule-sourced dung
# assemblage that is not competitively assembled in situ.
#
# Data: alldat_full[[1]] (== alldat_full$nopool) from 4_data_prep.R -- one row
#       per PCR replicate (109 rows / 55 biological samples / 1143 taxa), NOT
#       PCR-collapsed and NOT depth-filtered (min_depth_full=1000). This is the
#       canonical GLLVM/HMSC object precisely so PCR replicate can be an
#       explicit random level; CLR (compositional) normalisation stands in for
#       the library-size offset that a Gaussian HMSC cannot take. Do not swap in
#       alldat/alldat.rfy here -- those are PCR-collapsed with no per-replicate
#       structure left.
#
# Requires objects from 4_data_prep.R: alldat_full (list nopool/pool/pspool),
#       alldat (for the tax_table), funguild_otu (guild trait; a script-4 re-run
#       WIPES it -- re-run 8_functional_guilds.R to rebuild it if missing).
#
# Reference: no reference-script counterpart -- this is a new JSDM step with no
# equivalent in Root_fungi_DADA2 @ 65fbffa. It reuses this study's own
# conventions from 6_diversity_analyses.R / 7_gllvm.R: the Season(winter ref) x
# Year design, the FUNGuild guild grouping, the indivID->singleton fallback, the
# headless-plotting rule, and the Supplementary staging pattern.
#
# Outputs (canonical dirs):
#   models/  hmsc_clr_fit.rds, hmsc_pa_fit.rds, hmsc_libsize_fit.rds,
#            hmsc_clr_phylo_fit.rds
#   tables/  hmsc_otu_guild.csv, hmsc_convergence_psrf.csv,
#            hmsc_variance_partition.csv, hmsc_gamma_CrI.csv,
#            hmsc_probit_gamma_CrI.csv, hmsc_beta_season.csv,
#            hmsc_probit_beta_season.csv, hmsc_vs_gllvm_season.csv,
#            hmsc_omega_sample.csv, hmsc_predictive_R2.csv,
#            hmsc_predictive_R2_summary.csv, hmsc_rho_phylo.csv
#   plots/   hmsc_variance_partition.png, hmsc_gamma.png, hmsc_beta_season.png,
#            hmsc_vs_gllvm_season.png, hmsc_omega_sample.png,
#            hmsc_R2_explanatory_vs_cv.png, hmsc_trace_*.png
#   staged into Scripts_server/Supplementary/figures|tables (no .qmd yet).
#
# Run: conda run -n r_env Rscript Scripts/10_hmsc.R
#   HMSC_RUN_MODE=pilot   -> fast pipeline / rough-convergence check
#   HMSC_RUN_MODE=production (default) -> long MCMC; run as a background job.
# =============================================================================

suppressMessages({
  library(Hmsc)
  library(coda)
  library(corrplot)
  library(vegan)
  library(phyloseq)
  library(ape)
  library(ggplot2)
})

set.seed(20260717)
setwd("/home/daniel/Ptarmigan/trimmed/mergedPlates/")
load("eco_analysis.RData")
if (!exists("funguild_otu"))
  stop("funguild_otu not found in eco_analysis.RData -- a script-4 re-run wipes it; ",
       "restore it (scratchpad/restore_funguild.R) before running HMSC.")

# ---- Canonical output dirs (absolute) ---------------------------------------
out_dir  <- "/home/daniel/Ptarmigan/models"     # fitted models + (per convention) tables
plot_dir <- "/home/daniel/Ptarmigan/plots"
tab_dir  <- "/home/daniel/Ptarmigan/tables"     # user-requested tables/ dir (CSVs mirrored here)
supp_fig <- "/home/daniel/Ptarmigan/Scripts_server/Supplementary/figures"
supp_tab <- "/home/daniel/Ptarmigan/Scripts_server/Supplementary/tables"
for (d in c(out_dir, plot_dir, tab_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Headless plotting (see CLAUDE.md): a standing null device removes R's X11
# fallback for ANY base-graphics call (corrplot / Hmsc plots included). ggplot
# figures go through save_png(); base-graphics figures use raw png()/dev.off().
grDevices::pdf(NULL)
save_png <- function(path, plot_obj, width = 8, height = 8, res = 300) {
  grDevices::png(path, width = width, height = height, units = "in", res = res)
  on.exit(grDevices::dev.off())
  print(plot_obj)
}
# write a CSV to both models/ (staging) and tables/ (user-requested)
write_tab <- function(df, name, row.names = FALSE) {
  write.csv(df, file.path(out_dir, name), row.names = row.names)
  write.csv(df, file.path(tab_dir, name), row.names = row.names)
}

# ---- Shared helper (from 6_diversity_analyses.R) ----------------------------
otu_mat_of <- function(ps) {
  m <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) m <- t(m)
  m
}

# ---- MCMC intensity ---------------------------------------------------------
RUN_MODE <- Sys.getenv("HMSC_RUN_MODE", "production")
if (RUN_MODE == "pilot") {
  mc <- list(samples = 250, thin = 5, transient = 1250)
} else {
  mc <- list(samples = 1000, thin = 50, transient = 25000)   # 75k iters/chain
}
nChains  <- 4
# CPU politeness. This is a 96-core box; with no BLAS thread cap each parallel
# chain-process fans its matrix ops across ALL cores (nParallel x 96 threads
# oversubscribing 96 cores -> spikes to 100% and thrashing). Cap BLAS to 1
# thread/process and run a small number of parallel chains, so total load is
# ~nParallel cores. NOTE: OpenBLAS reads these at load, so the launch command
# should ALSO export them (OPENBLAS_NUM_THREADS=1 ...); the Sys.setenv here is a
# defensive fallback. Override via HMSC_NPARALLEL / HMSC_BLAS_THREADS.
blas_threads <- Sys.getenv("HMSC_BLAS_THREADS", "1")
Sys.setenv(OMP_NUM_THREADS = blas_threads, OPENBLAS_NUM_THREADS = blas_threads,
           MKL_NUM_THREADS = blas_threads, VECLIB_MAXIMUM_THREADS = blas_threads)
if (requireNamespace("RhpcBLASctl", quietly = TRUE))
  try(RhpcBLASctl::blas_set_num_threads(as.integer(blas_threads)), silent = TRUE)
nParallel <- min(nChains, as.integer(Sys.getenv("HMSC_NPARALLEL", "4")))
cat(sprintf("HMSC run mode = %s : %d samples x thin %d (transient %d), %d chains on %d cores (BLAS threads/proc=%s)\n",
            RUN_MODE, mc$samples, mc$thin, mc$transient, nChains, nParallel, blas_threads))

# =============================================================================
# SECTION 1 -- DATA ASSEMBLY
# =============================================================================
ps   <- alldat_full$nopool                       # PCR-rep level, no depth filter
Ymat <- otu_mat_of(ps)                            # samples (PCR reps) x taxa, counts
md   <- data.frame(as(sample_data(ps), "data.frame"), stringsAsFactors = FALSE)
stopifnot(all(c("Season","Year","Sample_ID_field","indivID","pcr_sample_id") %in% names(md)))

# Prevalence filter: OTUs present in >= 5 PCR-rep rows (matches GLLVM/H3 rule)
MIN_OTU_PREV <- 5
keep_otu <- colSums(Ymat > 0) >= MIN_OTU_PREV
Yk <- Ymat[, keep_otu, drop = FALSE]
libsize <- rowSums(Ymat)                          # library size from the FULL pre-filter table
log_libsize <- as.numeric(scale(log(libsize)))    # z-scored, for the robustness covariate
cat(sprintf("OTUs retained (present in >= %d samples): %d of %d; %d PCR-rep rows / %d biological samples\n",
            MIN_OTU_PREV, ncol(Yk), ncol(Ymat), nrow(Yk), length(unique(md$Sample_ID_field))))

# CLR response (Gaussian). decostand's clr path drops dimnames -- restore.
Yclr <- vegan::decostand(Yk, method = "clr", pseudocount = 1)
dimnames(Yclr) <- dimnames(Yk)
stopifnot(!anyNA(Yclr))

# Presence/absence response (probit robustness)
Ypa <- (Yk > 0) * 1
storage.mode(Ypa) <- "double"

# ---- Fixed-effect design ----------------------------------------------------
XData <- data.frame(
  Season = factor(md$Season, levels = c("winter","summer")),
  Year   = factor(md$Year),
  row.names = rownames(Yk)
)
XData_lib <- cbind(XData, log_libsize = log_libsize[match(rownames(Yk), rownames(Ymat))])
stopifnot(!anyNA(XData$Season), !anyNA(XData$Year))

# ---- Random-effect study design (biological sample / PCR replicate) ---------
# sample: Sample_ID_field (the biological dropping; groups the 2 PCR reps).
# pcr:    pcr_sample_id (unique per row -> the finest, observation-level latent
#         that generates the residual co-occurrence Omega).
# NOTE -- the "individual bird" level was DROPPED. indivID is microsat-confirmed
# for only ~4 repeat birds, so a bird level (with unidentified droppings as
# singletons) was ~1:1 with biological sample (51 vs 55 units); the bird and
# sample latent factors were confounded, which stopped the CLR-Gaussian
# community parameters (Gamma, variance partition) converging in the 3-level
# version (Season Gamma PSRF ~2.6, pilot<->production instability). Sample + PCR
# is the identifiable structure and mirrors the GLLVM's row.eff=~(1|Sample_ID_field).
studyDesign <- data.frame(
  sample = factor(md$Sample_ID_field),
  pcr    = factor(md$pcr_sample_id),
  row.names = rownames(Yk)
)
stopifnot(nlevels(studyDesign$pcr) == nrow(Yk))    # one PCR unit per observation
cat(sprintf("Random levels: %d biological samples, %d PCR replicates (bird level dropped -- see note)\n",
            nlevels(studyDesign$sample), nlevels(studyDesign$pcr)))

rL_sample <- HmscRandomLevel(units = levels(studyDesign$sample))
rL_pcr    <- HmscRandomLevel(units = levels(studyDesign$pcr))
ranLevels <- list(sample = rL_sample, pcr = rL_pcr)

# =============================================================================
# SECTION 2 -- GUILD TRAIT (5-level), reused FUNGuild logic
# =============================================================================
# Priority (first match wins): dung_saprotroph > plant_associated > pathotroph
# > dark_unassigned > other. Coprophily takes precedence (FUNGuild "Dung
# Saprotroph" OR a COPRO_GENERA genus). grepl logic verbatim from
# 6_diversity_analyses.R:882-899; COPRO_GENERA from :711-716.
COPRO_GENERA <- c("Sordaria","Podospora","Cercophora","Chaetomium","Schizothecium",
                  "Preussia","Delitschia","Pilobolus","Ascobolus","Saccobolus",
                  "Sporormiella","Coprinopsis","Thelebolus","Coniochaeta")
strip_rank <- function(x) sub("^[a-z]__", "", x)   # UNITE "g__Sporormiella" -> "Sporormiella"

otu_ids   <- colnames(Yk)
guild_str <- funguild_otu$Guild[match(otu_ids, funguild_otu$OTU_ID)]
tax_all   <- data.frame(as(tax_table(alldat$nopool), "matrix"), stringsAsFactors = FALSE)
otu_gen   <- strip_rank(tax_all$Genus[match(otu_ids, rownames(tax_all))])
otu_gen   <- ifelse(otu_gen %in% c("", "NA"), NA, otu_gen)

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
guild_lvls <- c("dung_saprotroph","plant_associated","pathotroph","dark_unassigned","other")
TrData <- data.frame(guild = factor(guild_class, levels = guild_lvls),
                     row.names = otu_ids)
TrFormula <- ~guild

cat("Guild-trait distribution across retained OTUs:\n"); print(table(TrData$guild))
otu_guild_tbl <- data.frame(OTU_ID = otu_ids, Genus = otu_gen,
                            guild_class = guild_class, FUNGuild = guild_str,
                            stringsAsFactors = FALSE)
write_tab(otu_guild_tbl, "hmsc_otu_guild.csv")

# =============================================================================
# SECTION 3 -- MODEL OBJECTS
# =============================================================================
mk_hmsc <- function(Y, XForm, Xd, distr, phyloTree = NULL) {
  Hmsc(Y = Y, XData = Xd, XFormula = XForm,
       TrData = TrData, TrFormula = TrFormula,
       phyloTree = phyloTree,
       studyDesign = studyDesign, ranLevels = ranLevels,
       distr = distr)
}
m_clr <- mk_hmsc(Yclr, ~Season + Year,               XData,     "normal")
m_pa  <- mk_hmsc(Ypa,  ~Season + Year,               XData,     "probit")
m_lib <- mk_hmsc(Yclr, ~Season + Year + log_libsize, XData_lib, "normal")

# Phylo variant (secondary): prune the tree to the retained OTUs and root it.
tree_obj <- readRDS("/home/daniel/Ptarmigan/trimmed/mergedPlates/tree.rds")
tr <- if (inherits(tree_obj, "phyloseq")) phy_tree(tree_obj) else tree_obj
tr_p <- tryCatch({
  t2 <- ape::keep.tip(tr, intersect(tr$tip.label, otu_ids))
  t2 <- ape::multi2di(t2)                                    # resolve polytomies
  if (!ape::is.rooted(t2)) t2 <- ape::root(t2, outgroup = t2$tip.label[1], resolve.root = TRUE)
  t2$edge.length[t2$edge.length <= 0] <- 1e-8                # vcv needs positive branches
  t2
}, error = function(e) { message("Phylo prune/root failed: ", conditionMessage(e)); NULL })
m_phy <- if (!is.null(tr_p) && setequal(tr_p$tip.label, otu_ids)) {
  mm <- mk_hmsc(Yclr, ~Season + Year, XData, "normal", phyloTree = tr_p)
  # Coarsen the rho grid (default 101 -> 26 points). The phylogenetic Beta
  # update + rho-grid marginal likelihood scale with the rho-grid size and are
  # the entire cost gap vs the main models (~5x slower per iter at 101 points);
  # the pilot posterior for rho occupies ~0.5-0.85, so a 0.04 grid resolves it.
  rv <- seq(0, 1, by = 0.04)
  setPriors(mm, rhopw = cbind(rv, c(0.5, rep(0.5 / (length(rv) - 1), length(rv) - 1))))
} else { message("Phylo variant skipped (tree/OTU tip mismatch)."); NULL }

# =============================================================================
# SECTION 4 -- FIT (MCMC), cached to models/ (delete an .rds to force a refit)
# =============================================================================
fit_or_load <- function(m, path, mcp = mc) {
  if (file.exists(path)) { cat("Loading cached fit:", basename(path), "\n"); return(readRDS(path)) }
  cat("Fitting:", basename(path), "...\n"); t0 <- Sys.time()
  m <- sampleMcmc(m, samples = mcp$samples, thin = mcp$thin, transient = mcp$transient,
                  nChains = nChains, nParallel = nParallel,
                  verbose = max(1, round((mcp$transient + mcp$samples * mcp$thin) / 10)))
  cat(sprintf("  done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  saveRDS(m, path); m
}
m_clr <- fit_or_load(m_clr, file.path(out_dir, "hmsc_clr_fit.rds"))
m_pa  <- fit_or_load(m_pa,  file.path(out_dir, "hmsc_pa_fit.rds"))
m_lib <- fit_or_load(m_lib, file.path(out_dir, "hmsc_libsize_fit.rds"))
# Phylo variant is SECONDARY (its only added deliverable is rho): ~50x slower per
# iter than the main models even with the coarse grid, so it gets a dedicated,
# capped budget independent of RUN_MODE (~3h at 4 chains) -- enough to resolve
# rho, not a full production fit.
mc_phy <- if (RUN_MODE == "pilot") mc else list(samples = 250, thin = 15, transient = 4000)
if (!is.null(m_phy)) m_phy <- fit_or_load(m_phy, file.path(out_dir, "hmsc_clr_phylo_fit.rds"), mc_phy)

# =============================================================================
# SECTION 5 -- CONVERGENCE (Gelman-Rubin PSRF on Beta and Gamma)
# =============================================================================
psrf_summary <- function(m, label) {
  mpost <- convertToCodaObject(m)
  gb <- gelman.diag(mpost$Beta,  multivariate = FALSE)$psrf[, 1]
  gg <- gelman.diag(mpost$Gamma, multivariate = FALSE)$psrf[, 1]
  eb <- effectiveSize(mpost$Beta); eg <- effectiveSize(mpost$Gamma)
  data.frame(model = label,
             beta_psrf_max = round(max(gb, na.rm = TRUE), 3),
             beta_psrf_med = round(median(gb, na.rm = TRUE), 3),
             beta_pct_below_1.1 = round(100 * mean(gb < 1.1, na.rm = TRUE), 1),
             beta_ess_min = round(min(eb, na.rm = TRUE)),
             gamma_psrf_max = round(max(gg, na.rm = TRUE), 3),
             gamma_psrf_med = round(median(gg, na.rm = TRUE), 3),
             gamma_ess_min = round(min(eg, na.rm = TRUE)),
             stringsAsFactors = FALSE)
}
conv <- rbind(psrf_summary(m_clr, "clr"), psrf_summary(m_pa, "pa"))
if (!is.null(m_phy)) conv <- rbind(conv, psrf_summary(m_phy, "clr_phylo"))
write_tab(conv, "hmsc_convergence_psrf.csv")
cat("Convergence (PSRF):\n"); print(conv)

# a few Beta/Gamma traceplots from the headline model
mpost_clr <- convertToCodaObject(m_clr)
png(file.path(plot_dir, "hmsc_trace_beta.png"), width = 9, height = 6, units = "in", res = 200)
plot(mpost_clr$Beta[, 1:min(4, ncol(mpost_clr$Beta[[1]]))]); dev.off()
png(file.path(plot_dir, "hmsc_trace_gamma.png"), width = 9, height = 6, units = "in", res = 200)
plot(mpost_clr$Gamma[, 1:min(4, ncol(mpost_clr$Gamma[[1]]))]); dev.off()

# =============================================================================
# SECTION 6.1 -- VARIANCE PARTITIONING (community-weighted + per-OTU)
# =============================================================================
xcn   <- colnames(m_clr$X)
vp_grp <- ifelse(grepl("^Year", xcn), 2L, 1L)     # intercept + Season -> "Season"; Year dummies -> "Year"
VP <- computeVariancePartitioning(m_clr, group = vp_grp, groupnames = c("Season","Year"))
vp_vals <- VP$vals                                  # components x OTU
vp_tbl <- data.frame(OTU_ID = colnames(vp_vals), t(round(vp_vals, 4)),
                     stringsAsFactors = FALSE, check.names = FALSE)
vp_tbl$Genus <- otu_gen[match(vp_tbl$OTU_ID, otu_ids)]
vp_tbl$guild <- guild_class[match(vp_tbl$OTU_ID, otu_ids)]
write_tab(vp_tbl, "hmsc_variance_partition.csv")
cwm <- rowMeans(vp_vals)                             # community-weighted (mean over OTUs) mean shares
cat("Community-mean variance partition:\n"); print(round(cwm, 3))
png(file.path(plot_dir, "hmsc_variance_partition.png"), width = 11, height = 6, units = "in", res = 800)
plotVariancePartitioning(m_clr, VP, las = 2, cex.names = 0.35,
                         main = "HMSC variance partitioning (CLR abundance)")
dev.off()

# =============================================================================
# SECTION 6.2 -- GAMMA (trait -> Season/Year response) with 95% CrI
# The CLR-Gaussian Gamma converges poorly on this data (see convergence table);
# the presence/absence PROBIT Gamma converges cleanly and is the reliable
# trait->Season inference. Both are written; report the probit.
# =============================================================================
gamma_table <- function(m, mp) {
  ge <- getPostEstimate(m, parName = "Gamma")
  gq <- summary(mp$Gamma)$quantiles
  trcn <- colnames(m$Tr); xc <- colnames(m$X)
  data.frame(
    covariate = rep(xc, times = length(trcn)),
    trait     = rep(trcn, each = length(xc)),
    mean      = round(as.vector(ge$mean),    3),
    support   = round(as.vector(ge$support), 3),
    CrI_2.5   = round(gq[, "2.5%"],  3),
    CrI_97.5  = round(gq[, "97.5%"], 3),
    stringsAsFactors = FALSE
  )
}
gamma_fig <- function(m, path, ttl) {
  ge <- getPostEstimate(m, parName = "Gamma")
  png(path, width = 8, height = 6, units = "in", res = 800)
  plotGamma(m, post = ge, param = "Support", supportLevel = 0.9, main = ttl)
  dev.off()
}
xcn2 <- colnames(m_clr$X)
gamma_tbl <- gamma_table(m_clr, mpost_clr)
write_tab(gamma_tbl, "hmsc_gamma_CrI.csv")
gamma_fig(m_clr, file.path(plot_dir, "hmsc_gamma.png"),
          "Gamma: guild trait -> abundance response (CLR)")
# Probit (occurrence) Gamma -- the converged, reportable trait->Season link
mpost_pa <- convertToCodaObject(m_pa)
gamma_pa <- gamma_table(m_pa, mpost_pa)
write_tab(gamma_pa, "hmsc_probit_gamma_CrI.csv")
gamma_fig(m_pa, file.path(plot_dir, "hmsc_probit_gamma.png"),
          "Gamma: guild trait -> summer occurrence (probit)")
cat("Probit Gamma (trait x covariate) Season rows:\n")
print(gamma_pa[gamma_pa$covariate == "Seasonsummer", ])

# =============================================================================
# SECTION 6.3 -- BETA (per-OTU Season niche) + GLLVM cross-check
# =============================================================================
geB   <- getPostEstimate(m_clr, parName = "Beta")    # covariates x OTU
season_row <- which(xcn2 == "Seasonsummer")
beta_tbl <- data.frame(
  OTU_ID   = colnames(m_clr$Y),
  Genus    = otu_gen,
  guild    = guild_class,
  beta_season_mean = round(geB$mean[season_row, ], 3),
  beta_season_support = round(geB$support[season_row, ], 3),
  stringsAsFactors = FALSE
)
# attach 95% CrI for the Season coefficient, joined by OTU (order-robust).
# Hmsc coda names Beta rows "B[Seasonsummer (C2), OTU1 (S1)]".
bq   <- summary(mpost_clr$Beta)$quantiles
srow <- grep("Seasonsummer", rownames(bq))
otu_of <- sub("^B\\[Seasonsummer \\([^)]*\\), (\\S+) \\(S[0-9]+\\)\\]$", "\\1", rownames(bq)[srow])
cri  <- data.frame(OTU_ID = otu_of,
                   beta_season_CrI2.5  = round(bq[srow, "2.5%"],  3),
                   beta_season_CrI97.5 = round(bq[srow, "97.5%"], 3),
                   stringsAsFactors = FALSE)
beta_tbl <- merge(beta_tbl, cri, by = "OTU_ID", all.x = TRUE)
beta_tbl <- beta_tbl[order(-beta_tbl$beta_season_mean), ]
write_tab(beta_tbl, "hmsc_beta_season.csv")

# Caterpillar of the strongest per-OTU seasonal niches (top+bottom 25 by mean).
bshow <- unique(rbind(head(beta_tbl, 25), tail(beta_tbl, 25)))
bshow$lab <- ifelse(is.na(bshow$Genus), bshow$OTU_ID, paste0(bshow$OTU_ID, " (", bshow$Genus, ")"))
bshow$lab <- factor(bshow$lab, levels = bshow$lab[order(bshow$beta_season_mean)])
p_beta <- ggplot(bshow, aes(beta_season_mean, lab, colour = guild)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_errorbarh(aes(xmin = beta_season_CrI2.5, xmax = beta_season_CrI97.5),
                 height = 0, alpha = 0.6) +
  geom_point() +
  labs(x = "HMSC Beta: Season (summer vs winter, CLR)", y = NULL,
       title = "Per-OTU seasonal niche (top + bottom 25 OTUs)",
       colour = "Guild") +
  theme_bw(base_size = 10)
save_png(file.path(plot_dir, "hmsc_beta_season.png"), p_beta, width = 8, height = 9, res = 800)

# -----------------------------------------------------------------------------
# SECTION 6.3b -- the SAME per-OTU Season niche from the PROBIT model
# Mirrors 6.3 exactly (same columns, same CrI parse), on m_pa instead of m_clr.
# Reason it exists: the CLR Beta is the weaker-converging of the two (see the
# convergence table -- median PSRF ~1.17, 40% of parameters below 1.1) whereas
# the probit Beta converges cleanly (median ~1.01, 99.6% below 1.1). The CLR
# table stays the Section 10.4 headline (it is an ABUNDANCE niche); this one is
# the occurrence-scale counterpart, and is what the main-text figure plots.
# -----------------------------------------------------------------------------
geBp <- getPostEstimate(m_pa, parName = "Beta")
season_row_pa <- which(colnames(m_pa$X) == "Seasonsummer")   # NOT xcn2: that is the CLR design
stopifnot(length(season_row_pa) == 1L)
beta_pa <- data.frame(
  OTU_ID   = colnames(m_pa$Y),
  Genus    = otu_gen,
  guild    = guild_class,
  beta_season_mean    = round(geBp$mean[season_row_pa, ],    3),
  beta_season_support = round(geBp$support[season_row_pa, ], 3),
  stringsAsFactors = FALSE
)
bqp   <- summary(mpost_pa$Beta)$quantiles
srowp <- grep("Seasonsummer", rownames(bqp))
otu_of_pa <- sub("^B\\[Seasonsummer \\([^)]*\\), (\\S+) \\(S[0-9]+\\)\\]$", "\\1", rownames(bqp)[srowp])
# Fail loudly if a coda naming change breaks the parse, rather than silently
# writing a table of NA intervals.
stopifnot(length(srowp) == ncol(m_pa$Y), setequal(otu_of_pa, beta_pa$OTU_ID))
crip <- data.frame(OTU_ID = otu_of_pa,
                   beta_season_CrI2.5  = round(bqp[srowp, "2.5%"],  3),
                   beta_season_CrI97.5 = round(bqp[srowp, "97.5%"], 3),
                   stringsAsFactors = FALSE)
beta_pa <- merge(beta_pa, crip, by = "OTU_ID", all.x = TRUE)
beta_pa <- beta_pa[order(-beta_pa$beta_season_mean), ]
stopifnot(!anyNA(beta_pa$beta_season_CrI2.5), !anyNA(beta_pa$beta_season_CrI97.5))
write_tab(beta_pa, "hmsc_probit_beta_season.csv")
cat(sprintf("Probit Beta(Season): %d OTUs, %d with a 95%% CrI excluding 0, %d at support >= 0.95 or <= 0.05\n",
            nrow(beta_pa),
            sum(beta_pa$beta_season_CrI2.5 > 0 | beta_pa$beta_season_CrI97.5 < 0),
            sum(beta_pa$beta_season_support >= 0.95 | beta_pa$beta_season_support <= 0.05)))

gllvm_path <- file.path(out_dir, "gllvm_perOTU_season_coef.csv")
if (file.exists(gllvm_path)) {
  gl <- read.csv(gllvm_path, stringsAsFactors = FALSE)
  otu_col  <- intersect(c("OTU_ID","OTU","otu","taxon"), names(gl))[1]
  coef_col <- intersect(c("beta_season","Estimate","coef","estimate","season_coef","beta"), names(gl))[1]
  if (!is.na(otu_col) && !is.na(coef_col)) {
    glc <- setNames(gl[, c(otu_col, coef_col)], c("OTU_ID","gllvm_season"))
    # Drop GLLVM near-separation OTUs: their on/off coefs are +/-1000s and would
    # dominate the rank comparison without reflecting graded abundance.
    if ("separation" %in% names(gl)) glc <- glc[!(gl$separation %in% c(TRUE, "TRUE")), , drop = FALSE]
    cmp <- merge(beta_tbl[, c("OTU_ID","beta_season_mean","guild")], glc, by = "OTU_ID")
    write_tab(cmp, "hmsc_vs_gllvm_season.csv")
    rho_s <- suppressWarnings(cor(cmp$beta_season_mean, cmp$gllvm_season, method = "spearman"))
    cat(sprintf("HMSC vs GLLVM Season coef: n=%d overlapping OTUs (near-separation dropped), Spearman rho=%.3f\n",
                nrow(cmp), rho_s))
    p_cmp <- ggplot(cmp, aes(gllvm_season, beta_season_mean, colour = guild)) +
      geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
      geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
      geom_point(alpha = 0.8) +
      labs(x = "GLLVM Season coefficient (Section 6)",
           y = "HMSC Beta Season (CLR)",
           title = "Per-OTU Season niche: HMSC vs GLLVM",
           subtitle = sprintf("Spearman rho = %.3f (n = %d OTUs)", rho_s, nrow(cmp))) +
      theme_bw(base_size = 12)
    save_png(file.path(plot_dir, "hmsc_vs_gllvm_season.png"), p_cmp, width = 7.5, height = 6, res = 800)
  } else cat("GLLVM coef CSV present but expected columns not found -- cross-check skipped.\n")
} else cat("No gllvm_perOTU_season_coef.csv -- run 7_gllvm.R for the Beta cross-check.\n")

# =============================================================================
# SECTION 6.4 -- OMEGA (residual associations) + network plot
# =============================================================================
assoc <- computeAssociations(m_clr)               # one entry per random level (ranLevels order)
s_idx  <- match("sample", names(ranLevels))       # biological-sample level = residual co-occurrence
OmegaS <- assoc[[s_idx]]$mean
suppS  <- assoc[[s_idx]]$support
OmegaS_thr <- OmegaS
OmegaS_thr[suppS < 0.95 & suppS > 0.05] <- 0      # keep only strongly supported associations
write_tab(data.frame(OTU_ID = rownames(OmegaS), round(OmegaS, 3),
                     check.names = FALSE, stringsAsFactors = FALSE),
          "hmsc_omega_sample.csv")
# plot only OTUs with at least one supported association (keeps the matrix legible)
has_assoc <- rowSums(OmegaS_thr != 0) > 1
if (sum(has_assoc) >= 3) {
  M <- OmegaS_thr[has_assoc, has_assoc]
  png(file.path(plot_dir, "hmsc_omega_sample.png"),
      width = 10, height = 10, units = "in", res = 800)
  corrplot(M, method = "color", type = "lower", order = "hclust",
           tl.cex = 0.35, tl.col = "black", diag = FALSE,
           col = colorRampPalette(c("#2166AC","white","#B2182B"))(200),
           title = "HMSC residual associations (biological-sample level, >=0.95 support)",
           mar = c(0,0,2,0))
  dev.off()
} else cat("Too few supported associations to plot Omega network.\n")

# =============================================================================
# SECTION 6.5 -- EXPLANATORY vs 2-FOLD CV PREDICTIVE R2
# =============================================================================
# Explanatory R2 uses the existing posterior (cheap); CV refits the model per
# fold at its stored MCMC settings (expensive), so CV is run only for the
# headline CLR model. Explanatory R2 is reported for CLR and the P/A robustness.
expl_R2 <- function(m, label) {
  mfE <- evaluateModelFit(hM = m, predY = computePredictedValues(m))
  eR2 <- if (!is.null(mfE$R2)) mfE$R2 else mfE$TjurR2
  data.frame(model = label, OTU_ID = colnames(m$Y),
             expl_R2 = round(eR2, 3),
             expl_AUC = if (!is.null(mfE$AUC)) round(mfE$AUC, 3) else NA_real_,
             stringsAsFactors = FALSE)
}
eclr <- expl_R2(m_clr, "clr")
epa  <- expl_R2(m_pa,  "pa")
# 2-fold CV on CLR, folds split by biological sample so PCR reps never leak.
partC <- createPartition(m_clr, nfolds = 2, column = "sample")
mfC   <- evaluateModelFit(hM = m_clr, predY =
           computePredictedValues(m_clr, partition = partC, nParallel = nParallel))
eclr$cv_R2 <- round(if (!is.null(mfC$R2)) mfC$R2 else mfC$TjurR2, 3)
epa$cv_R2  <- NA_real_
r2_per <- rbind(eclr[, c("model","OTU_ID","expl_R2","cv_R2","expl_AUC")],
                epa[,  c("model","OTU_ID","expl_R2","cv_R2","expl_AUC")])
r2_sum <- data.frame(
  model = c("clr","pa"),
  mean_expl_R2 = round(c(mean(eclr$expl_R2, na.rm=TRUE), mean(epa$expl_R2, na.rm=TRUE)), 3),
  mean_cv_R2   = round(c(mean(eclr$cv_R2, na.rm=TRUE), NA_real_), 3),
  stringsAsFactors = FALSE)
write_tab(r2_per, "hmsc_predictive_R2.csv")
write_tab(r2_sum, "hmsc_predictive_R2_summary.csv")
cat("Predictive performance (community mean):\n"); print(r2_sum)
r2c <- eclr
p_r2 <- ggplot(r2c, aes(expl_R2, cv_R2)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey80") +
  geom_point(alpha = 0.7, colour = "#2166AC") +
  labs(x = "Explanatory R2", y = "2-fold CV R2",
       title = "HMSC (CLR) explanatory vs cross-validated fit",
       subtitle = sprintf("mean expl R2 = %.3f, mean CV R2 = %.3f",
                          r2_sum$mean_expl_R2[1], r2_sum$mean_cv_R2[1])) +
  theme_bw(base_size = 12)
save_png(file.path(plot_dir, "hmsc_R2_explanatory_vs_cv.png"), p_r2, width = 6.5, height = 6, res = 800)

# =============================================================================
# SECTION 6.6 -- PHYLO VARIANT: rho posterior + Beta/Gamma-unchanged check
# =============================================================================
if (!is.null(m_phy)) {
  mp <- convertToCodaObject(m_phy)
  rho_draws <- as.matrix(mp$Rho)[, 1]
  rho_tbl <- data.frame(
    param = "rho",
    median = round(median(rho_draws), 3),
    CrI_2.5 = round(quantile(rho_draws, 0.025), 3),
    CrI_97.5 = round(quantile(rho_draws, 0.975), 3),
    P_gt0 = round(mean(rho_draws > 0), 3),
    ess = round(effectiveSize(mp$Rho)),
    psrf = round(tryCatch(gelman.diag(mp$Rho)$psrf[1], error = function(e) NA_real_), 3),
    stringsAsFactors = FALSE
  )
  # does the phylo prior move the headline Season niches / Gamma?
  geB_phy <- getPostEstimate(m_phy, parName = "Beta")$mean[season_row, ]
  geG_phy <- getPostEstimate(m_phy, parName = "Gamma")$mean
  rho_tbl$beta_season_cor_vs_main <- round(cor(geB$mean[season_row, ], geB_phy), 3)
  rho_tbl$gamma_cor_vs_main       <- round(cor(as.vector(geG$mean), as.vector(geG_phy)), 3)
  write_tab(rho_tbl, "hmsc_rho_phylo.csv")
  cat(sprintf("Phylo rho: median=%.3f [%.3f, %.3f], P(rho>0)=%.3f; Beta/Gamma corr vs main = %.3f / %.3f\n",
              rho_tbl$median, rho_tbl$CrI_2.5, rho_tbl$CrI_97.5, rho_tbl$P_gt0,
              rho_tbl$beta_season_cor_vs_main, rho_tbl$gamma_cor_vs_main))
}

# =============================================================================
# SECTION 7 -- STAGE FIGURES/TABLES INTO Supplementary (no .qmd section yet)
# =============================================================================
if (dir.exists(dirname(supp_fig))) {
  dir.create(supp_fig, showWarnings = FALSE, recursive = TRUE)
  dir.create(supp_tab, showWarnings = FALSE, recursive = TRUE)
  figs <- c("hmsc_variance_partition.png","hmsc_gamma.png","hmsc_probit_gamma.png",
            "hmsc_beta_season.png","hmsc_vs_gllvm_season.png","hmsc_omega_sample.png",
            "hmsc_R2_explanatory_vs_cv.png")
  tabs <- c("hmsc_otu_guild.csv","hmsc_convergence_psrf.csv","hmsc_variance_partition.csv",
            "hmsc_gamma_CrI.csv","hmsc_probit_gamma_CrI.csv","hmsc_beta_season.csv",
            "hmsc_probit_beta_season.csv","hmsc_vs_gllvm_season.csv",
            "hmsc_predictive_R2.csv","hmsc_predictive_R2_summary.csv",
            "hmsc_rho_phylo.csv")
  figs <- figs[file.exists(file.path(plot_dir, figs))]
  tabs <- tabs[file.exists(file.path(out_dir,  tabs))]
  invisible(file.copy(file.path(plot_dir, figs), supp_fig, overwrite = TRUE))
  invisible(file.copy(file.path(out_dir,  tabs), supp_tab, overwrite = TRUE))
  cat("Staged HMSC figures/tables into Supplementary/figures|tables\n")
}
cat("10_hmsc.R complete.\n")
