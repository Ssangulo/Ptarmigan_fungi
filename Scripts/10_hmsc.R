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
# seasonal niches).
#
# MODEL STRUCTURE (revised 2026-09-29 -- see hmsc_converge/MEMO_hmsc_convergence.md
# on branch exp/hmsc-converge for the full convergence study):
#   HEADLINE  m_clr  : CLR abundance, AVERAGED over each dropping's PCR
#                      replicates (55 rows), biological-sample latent level
#                      capped at 2 factors, YScale = TRUE.
#   COMPLEMENT m_rep : the same CLR at the PCR-replicate level (109 rows) with
#                      BOTH sample (2 factors) and PCR (1 factor) latent levels.
#   ROBUSTNESS m_pa  : presence/absence probit, PCR-replicate level, unchanged.
# Why: the previous CLR model (7 sample + 5 PCR factors, uncapped) did not
# converge -- its chains disagreed on the Season effects and the variance
# partition, because Season is constant within a dropping and unconstrained
# sample-level factors absorbed it differently in each chain. Capping the
# factors fixes the Season effects; at the PCR-replicate level the sample-vs-PCR
# variance split still has two posterior modes (m_rep), so the headline model
# averages the replicates instead and the PCR-replicate variance is reported
# model-free (Section 6.7). m_clr and m_rep agree per OTU (Pearson ~0.98).
#
# Data: alldat_full[[1]] (== alldat_full$nopool) from 4_data_prep.R -- one row
#       per PCR replicate (109 rows / 55 biological samples / 1143 taxa), NOT
#       PCR-collapsed and NOT depth-filtered (min_depth_full=1000). CLR
#       (compositional) normalisation stands in for the library-size offset
#       that a Gaussian HMSC cannot take; it is computed PER REPLICATE and then
#       averaged for m_clr. Do not swap in alldat/alldat.rfy here -- those are
#       PCR-collapsed by SUMMING counts, which is a different estimator.
#
# UNITS GOTCHA: with YScale = TRUE, Hmsc 3.3-7 stores Beta/Gamma on the per-OTU
# STANDARDISED response scale (not back-transformed). Section 6.3 multiplies
# the Season Beta by each OTU's response SD (YScalePar[2, ]) so the Beta tables
# are in CLR units. Support and the variance partition are scale-invariant;
# Gamma is left on the standardised scale (comparable across OTUs).
#
# Requires objects from 4_data_prep.R: alldat_full (list nopool/pool/pspool),
#       alldat (for the tax_table), funguild_otu (guild trait; a script-4 re-run
#       WIPES it -- re-run 8_functional_guilds.R to rebuild it if missing).
#
# Reference: no reference-script counterpart -- this is a new JSDM step with no
# equivalent in Root_fungi_DADA2 @ 65fbffa. It reuses this study's own
# conventions from 6_diversity_analyses.R / 7_gllvm.R: the Season(winter ref) x
# Year design, the FUNGuild guild grouping, the headless-plotting rule, and the
# Supplementary staging pattern.
#
# Outputs (under HMSC_OUT_ROOT, default /home/daniel/Ptarmigan):
#   models/  hmsc_clr_fit.rds, hmsc_rep_fit.rds, hmsc_pa_fit.rds,
#            hmsc_clr_phylo_fit.rds (only with HMSC_RUN_PHYLO=1)
#   tables/  hmsc_otu_guild.csv, hmsc_convergence_psrf.csv,
#            hmsc_variance_partition.csv, hmsc_rep_variance_partition.csv,
#            hmsc_rep_vp_by_chain.csv, hmsc_gamma_CrI.csv,
#            hmsc_probit_gamma_CrI.csv, hmsc_beta_season.csv,
#            hmsc_rep_beta_season.csv, hmsc_rep_vs_main.csv,
#            hmsc_probit_beta_season.csv, hmsc_vs_gllvm_season.csv,
#            hmsc_omega_sample.csv, hmsc_predictive_R2.csv,
#            hmsc_predictive_R2_summary.csv, hmsc_pcr_replicate_variance.csv,
#            hmsc_pcr_replicate_variance_perOTU.csv, hmsc_rho_phylo.csv
#   plots/   hmsc_variance_partition.png, hmsc_gamma.png, hmsc_probit_gamma.png,
#            hmsc_beta_season.png, hmsc_vs_gllvm_season.png,
#            hmsc_omega_sample.png, hmsc_R2_explanatory_vs_cv.png,
#            hmsc_trace_*.png
#   staged into HMSC_SUPP_DIR/figures|tables.
#
# Run: conda run -n r_env Rscript Scripts/10_hmsc.R
#   HMSC_RUN_MODE=pilot   -> fast pipeline / rough-convergence check
#   HMSC_RUN_MODE=production (default) -> long MCMC; run as a background job.
#   HMSC_OUT_ROOT / HMSC_SUPP_DIR -> redirect ALL outputs (e.g. from a worktree,
#     so an experimental run cannot overwrite main's models/plots/tables).
#   HMSC_RUN_PHYLO=1 -> also fit the phylogeny variant (slow; default off).
#   Launch long runs from a COPY of this file: Rscript reads the script
#   incrementally, so editing it mid-run corrupts the tail of the run.
# =============================================================================

suppressMessages({
  library(Hmsc)
  library(coda)
  library(corrplot)
  library(vegan)
  library(phyloseq)
  library(ape)
  library(ggplot2)
  library(lme4)
})

set.seed(20260717)
setwd("/home/daniel/Ptarmigan/trimmed/mergedPlates/")
load("eco_analysis.RData")
if (!exists("funguild_otu"))
  stop("funguild_otu not found in eco_analysis.RData -- a script-4 re-run wipes it; ",
       "re-run 8_functional_guilds.R (Section 2d rebuilds it) before running HMSC.")

# ---- Output dirs (absolute; overridable so a worktree run stays isolated) ----
OUT_ROOT <- Sys.getenv("HMSC_OUT_ROOT", "/home/daniel/Ptarmigan")
SUPP_DIR <- Sys.getenv("HMSC_SUPP_DIR", "/home/daniel/Ptarmigan/Scripts_server/Supplementary")
CANON_MODELS <- "/home/daniel/Ptarmigan/models"   # read-only inputs (GLLVM table, cached probit)
out_dir  <- file.path(OUT_ROOT, "models")        # fitted models + (per convention) tables
plot_dir <- file.path(OUT_ROOT, "plots")
tab_dir  <- file.path(OUT_ROOT, "tables")        # user-requested tables/ dir (CSVs mirrored here)
supp_fig <- file.path(SUPP_DIR, "figures")
supp_tab <- file.path(SUPP_DIR, "tables")
for (d in c(out_dir, plot_dir, tab_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
cat("Outputs ->", OUT_ROOT, "| staging ->", SUPP_DIR, "\n")

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
cat(sprintf("OTUs retained (present in >= %d samples): %d of %d; %d PCR-rep rows / %d biological samples\n",
            MIN_OTU_PREV, ncol(Yk), ncol(Ymat), nrow(Yk), length(unique(md$Sample_ID_field))))

# CLR response per PCR replicate (Gaussian). decostand's clr path drops dimnames -- restore.
Yclr_rep <- vegan::decostand(Yk, method = "clr", pseudocount = 1)
dimnames(Yclr_rep) <- dimnames(Yk)
stopifnot(!anyNA(Yclr_rep))

# Presence/absence response (probit robustness)
Ypa <- (Yk > 0) * 1
storage.mode(Ypa) <- "double"

# ---- Fixed-effect design (PCR-replicate level) -------------------------------
XData <- data.frame(
  Season = factor(md$Season, levels = c("winter","summer")),
  Year   = factor(md$Year),
  row.names = rownames(Yk)
)
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

# Default (uncapped) random levels -- used ONLY by the probit, whose
# specification is unchanged from the 2026-07-21 fit so its cache stays valid.
rL_sample <- HmscRandomLevel(units = levels(studyDesign$sample))
rL_pcr    <- HmscRandomLevel(units = levels(studyDesign$pcr))
ranLevels <- list(sample = rL_sample, pcr = rL_pcr)

# CAPPED random levels for the CLR models. Uncapped, the sample level took 7
# factors: enough to reproduce the Season contrast (Season is constant within a
# dropping), so each chain split the Season signal differently between Beta and
# the factors. Two sample factors + one PCR factor removes that freedom.
rLc_sample <- setPriors(HmscRandomLevel(units = levels(studyDesign$sample)), nfMax = 2, nfMin = 2)
rLc_pcr    <- setPriors(HmscRandomLevel(units = levels(studyDesign$pcr)),    nfMax = 1, nfMin = 1)
ranLevels_rep <- list(sample = rLc_sample, pcr = rLc_pcr)

# ---- HEADLINE response: CLR averaged over each dropping's PCR replicates ----
# One row per dropping, so the two replicates can never be counted as two
# independent droppings; replicate noise is averaged into the response and is
# reported separately, model-free, in Section 6.7.
drop_id <- factor(md$Sample_ID_field)
stopifnot(all(tapply(md$Season, drop_id, function(z) length(unique(z))) == 1),
          all(tapply(md$Year,   drop_id, function(z) length(unique(z))) == 1))
Yclr <- rowsum(Yclr_rep, drop_id) / as.vector(table(drop_id))   # rows in levels(drop_id) order
first <- match(rownames(Yclr), as.character(drop_id))
XData_drop <- data.frame(Season = XData$Season[first], Year = XData$Year[first],
                         row.names = rownames(Yclr))
studyDesign_drop <- data.frame(sample = factor(rownames(Yclr)), row.names = rownames(Yclr))
ranLevels_drop <- list(sample = setPriors(HmscRandomLevel(units = levels(studyDesign_drop$sample)),
                                          nfMax = 2, nfMin = 2))
cat(sprintf("Headline CLR response: %d droppings (%d with 2 PCR replicates averaged, %d with 1)\n",
            nrow(Yclr), sum(table(drop_id) == 2), sum(table(drop_id) == 1)))

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
mk_hmsc <- function(Y, XForm, Xd, distr, sD = studyDesign, rL = ranLevels,
                    YScale = FALSE, phyloTree = NULL) {
  Hmsc(Y = Y, XData = Xd, XFormula = XForm,
       TrData = TrData, TrFormula = TrFormula,
       phyloTree = phyloTree,
       studyDesign = sD, ranLevels = rL,
       distr = distr, YScale = YScale)
}
m_clr <- mk_hmsc(Yclr,     ~Season + Year, XData_drop, "normal",
                 sD = studyDesign_drop, rL = ranLevels_drop, YScale = TRUE)
m_rep <- mk_hmsc(Yclr_rep, ~Season + Year, XData,      "normal",
                 sD = studyDesign, rL = ranLevels_rep, YScale = TRUE)
m_pa  <- mk_hmsc(Ypa,      ~Season + Year, XData,      "probit")
# (The former library-size robustness fit, m_lib, was never summarised anywhere
# and has been removed.)

# Phylo variant (secondary, off by default): prune the tree to the retained
# OTUs, root it, and attach it to the HEADLINE structure.
RUN_PHYLO <- Sys.getenv("HMSC_RUN_PHYLO", "0") == "1"
m_phy <- NULL
if (RUN_PHYLO) {
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
    mm <- mk_hmsc(Yclr, ~Season + Year, XData_drop, "normal", sD = studyDesign_drop,
                  rL = ranLevels_drop, YScale = TRUE, phyloTree = tr_p)
    # Coarsen the rho grid (default 101 -> 26 points): the rho-grid marginal
    # likelihood scales with the grid size and is the whole cost gap vs the
    # main model.
    rv <- seq(0, 1, by = 0.04)
    setPriors(mm, rhopw = cbind(rv, c(0.5, rep(0.5 / (length(rv) - 1), length(rv) - 1))))
  } else { message("Phylo variant skipped (tree/OTU tip mismatch)."); NULL }
}

# =============================================================================
# SECTION 4 -- FIT (MCMC), cached to models/ (delete an .rds to force a refit)
# =============================================================================
# `reuse`: a read-only cache elsewhere, loaded when `path` is absent. Used for
# the probit, whose specification is unchanged, so a worktree run does not have
# to refit it. The caller checks that the cached data match.
fit_or_load <- function(m, path, mcp = mc, reuse = NULL) {
  if (file.exists(path)) { cat("Loading cached fit:", path, "\n"); return(readRDS(path)) }
  if (!is.null(reuse) && file.exists(reuse)) { cat("Reusing cached fit:", reuse, "\n"); return(readRDS(reuse)) }
  cat("Fitting:", basename(path), "...\n"); t0 <- Sys.time()
  m <- sampleMcmc(m, samples = mcp$samples, thin = mcp$thin, transient = mcp$transient,
                  nChains = nChains, nParallel = nParallel,
                  verbose = max(1, round((mcp$transient + mcp$samples * mcp$thin) / 10)))
  cat(sprintf("  done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  saveRDS(m, path); m
}
m_clr <- fit_or_load(m_clr, file.path(out_dir, "hmsc_clr_fit.rds"))
m_rep <- fit_or_load(m_rep, file.path(out_dir, "hmsc_rep_fit.rds"))
m_pa  <- fit_or_load(m_pa,  file.path(out_dir, "hmsc_pa_fit.rds"),
                     reuse = file.path(CANON_MODELS, "hmsc_pa_fit.rds"))
stopifnot(identical(dim(m_pa$Y), dim(Ypa)), all(m_pa$Y == Ypa),
          identical(colnames(m_pa$Y), otu_ids))          # cached probit = this data
mc_phy <- if (RUN_MODE == "pilot") mc else list(samples = 250, thin = 15, transient = 4000)
if (!is.null(m_phy)) m_phy <- fit_or_load(m_phy, file.path(out_dir, "hmsc_clr_phylo_fit.rds"), mc_phy)

# =============================================================================
# SECTION 5 -- CONVERGENCE (Gelman-Rubin PSRF on Beta and Gamma, plus
# per-chain agreement on the quantities the figures actually report)
# =============================================================================
# PSRF alone missed the old model's failure (OTU1 ESS 2,157 at PSRF 4.2: each
# chain mixed well within its own mode). The per-chain columns compare the
# chains directly: max range across chains of the guild-mean Season Beta
# (standardised units) and of the community-mean variance components.
xcn_of <- function(m) colnames(m$X)
vp_group <- function(m) ifelse(grepl("^Year", xcn_of(m)), 2L, 1L)
chain_checks <- function(m) {
  sr <- which(xcn_of(m) == "Seasonsummer")
  gm <- sapply(m$postList, function(ch) {
    b <- colMeans(t(sapply(ch, function(s) s$Beta[sr, ])))
    tapply(b, TrData$guild[match(colnames(m$Y), otu_ids)], mean)
  })
  vp <- sapply(seq_along(m$postList), function(k) {
    m1 <- m; m1$postList <- m$postList[k]
    rowMeans(computeVariancePartitioning(m1, group = vp_group(m),
                                         groupnames = c("Season", "Year"))$vals)
  })
  c(guild_mean_chain_range_max = max(apply(gm, 1, function(x) diff(range(x)))),
    vp_chain_range_max         = max(apply(vp, 1, function(x) diff(range(x)))))
}
psrf_summary <- function(m, label) {
  mpost <- convertToCodaObject(m)
  gb <- gelman.diag(mpost$Beta,  multivariate = FALSE)$psrf[, 1]
  gg <- gelman.diag(mpost$Gamma, multivariate = FALSE)$psrf[, 1]
  eb <- effectiveSize(mpost$Beta); eg <- effectiveSize(mpost$Gamma)
  gs <- gb[grepl("Seasonsummer", names(gb), fixed = TRUE)]
  cc <- chain_checks(m)
  data.frame(model = label, n_rows = m$ny, n_otu = m$ns,
             beta_psrf_max = round(max(gb, na.rm = TRUE), 3),
             beta_psrf_med = round(median(gb, na.rm = TRUE), 3),
             beta_pct_below_1.1 = round(100 * mean(gb < 1.1, na.rm = TRUE), 1),
             season_psrf_max = round(max(gs), 3),
             season_pct_below_1.1 = round(100 * mean(gs < 1.1), 1),
             beta_ess_min = round(min(eb, na.rm = TRUE)),
             gamma_psrf_max = round(max(gg, na.rm = TRUE), 3),
             gamma_psrf_med = round(median(gg, na.rm = TRUE), 3),
             gamma_ess_min = round(min(eg, na.rm = TRUE)),
             guild_mean_chain_range_max = round(cc[[1]], 3),
             vp_chain_range_max = round(cc[[2]], 3),
             stringsAsFactors = FALSE)
}
conv <- rbind(psrf_summary(m_clr, "clr"), psrf_summary(m_rep, "clr_rep"), psrf_summary(m_pa, "pa"))
if (!is.null(m_phy)) conv <- rbind(conv, psrf_summary(m_phy, "clr_phylo"))
write_tab(conv, "hmsc_convergence_psrf.csv")
cat("Convergence (PSRF + per-chain agreement):\n"); print(conv)

# a few Beta/Gamma traceplots from the headline model
mpost_clr <- convertToCodaObject(m_clr)
png(file.path(plot_dir, "hmsc_trace_beta.png"), width = 9, height = 6, units = "in", res = 200)
plot(mpost_clr$Beta[, 1:min(4, ncol(mpost_clr$Beta[[1]]))]); dev.off()
png(file.path(plot_dir, "hmsc_trace_gamma.png"), width = 9, height = 6, units = "in", res = 200)
plot(mpost_clr$Gamma[, 1:min(4, ncol(mpost_clr$Gamma[[1]]))]); dev.off()

# =============================================================================
# SECTION 6.1 -- VARIANCE PARTITIONING (community-weighted + per-OTU)
# =============================================================================
vp_table <- function(m) {
  VP <- computeVariancePartitioning(m, group = vp_group(m), groupnames = c("Season","Year"))
  vp_tbl <- data.frame(OTU_ID = colnames(VP$vals), t(round(VP$vals, 4)),
                       stringsAsFactors = FALSE, check.names = FALSE)
  vp_tbl$Genus <- otu_gen[match(vp_tbl$OTU_ID, otu_ids)]
  vp_tbl$guild <- guild_class[match(vp_tbl$OTU_ID, otu_ids)]
  list(VP = VP, tbl = vp_tbl)
}
vpc <- vp_table(m_clr); VP <- vpc$VP
write_tab(vpc$tbl, "hmsc_variance_partition.csv")
cat("Community-mean variance partition (headline):\n"); print(round(rowMeans(VP$vals), 3))
png(file.path(plot_dir, "hmsc_variance_partition.png"), width = 11, height = 6, units = "in", res = 800)
plotVariancePartitioning(m_clr, VP, las = 2, cex.names = 0.35,
                         main = "HMSC variance partitioning (CLR abundance, replicate-averaged)")
dev.off()

# Complement: the replicate-level model's partition, pooled AND per chain. Its
# sample-vs-PCR split has two posterior modes (chains 1,4 vs 2,3 in the
# 2026-09-29 fit), so the per-chain table is what should be quoted, as a range.
vpr <- vp_table(m_rep)
write_tab(vpr$tbl, "hmsc_rep_variance_partition.csv")
vp_rep_chain <- t(sapply(seq_along(m_rep$postList), function(k) {
  m1 <- m_rep; m1$postList <- m_rep$postList[k]
  rowMeans(computeVariancePartitioning(m1, group = vp_group(m_rep),
                                       groupnames = c("Season", "Year"))$vals)
}))
vp_rep_chain <- data.frame(chain = c(seq_len(nrow(vp_rep_chain)), "pooled"),
                           round(rbind(vp_rep_chain, rowMeans(vpr$VP$vals)), 4),
                           check.names = FALSE)
write_tab(vp_rep_chain, "hmsc_rep_vp_by_chain.csv")
cat("Replicate-level model, community-mean VP by chain:\n"); print(vp_rep_chain)

# =============================================================================
# SECTION 6.2 -- GAMMA (trait -> Season/Year response) with 95% CrI
# With capped factors the CLR Gamma converges (see convergence table), so it is
# reportable on the abundance scale alongside the probit. CLR Gamma is on the
# STANDARDISED scale (YScale; see header).
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
          "Gamma: guild trait -> abundance response (CLR, standardised)")
cat("CLR Gamma (trait x covariate) Season rows:\n")
print(gamma_tbl[gamma_tbl$covariate == "Seasonsummer", ])
# Probit (occurrence) Gamma -- unchanged model
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
# Season Beta per OTU from the pooled posterior draws, back-transformed to CLR
# units (x each OTU's response SD; YScale gotcha, see header). Support is
# P(Beta > 0), unaffected by the positive rescaling. beta_season_std_mean keeps
# the standardised value (comparable across OTUs).
beta_season_table <- function(m) {
  sr   <- which(colnames(m$X) == "Seasonsummer")
  post <- poolMcmcChains(m$postList)
  Bd   <- t(sapply(post, function(s) s$Beta[sr, ]))            # draws x OTU
  ysd  <- if (!is.null(m$YScalePar)) m$YScalePar[2, ] else rep(1, m$ns)
  Br   <- sweep(Bd, 2, ysd, "*")
  gi   <- match(colnames(m$Y), otu_ids)
  data.frame(OTU_ID = colnames(m$Y), Genus = otu_gen[gi], guild = guild_class[gi],
             beta_season_mean    = round(colMeans(Br), 3),
             beta_season_support = round(colMeans(Bd > 0), 3),
             beta_season_CrI2.5  = round(apply(Br, 2, quantile, 0.025), 3),
             beta_season_CrI97.5 = round(apply(Br, 2, quantile, 0.975), 3),
             beta_season_std_mean = round(colMeans(Bd), 3),
             stringsAsFactors = FALSE)
}
beta_tbl <- beta_season_table(m_clr)
beta_tbl <- beta_tbl[order(-beta_tbl$beta_season_mean), ]
stopifnot(!anyNA(beta_tbl$beta_season_CrI2.5), nrow(beta_tbl) == length(otu_ids))
write_tab(beta_tbl, "hmsc_beta_season.csv")
# Sanity check on the back-transform: for the strongest OTUs the model's Season
# effect must be of the same size as the raw winter/summer difference in mean
# CLR (it is adjusted for Year and shrunk, so not identical).
raw_diff <- colMeans(Yclr[XData_drop$Season == "summer", ]) - colMeans(Yclr[XData_drop$Season == "winter", ])
chk <- head(beta_tbl$OTU_ID, 3)
cat("Back-transform check (model Beta vs raw summer-winter CLR difference):\n")
print(data.frame(OTU = chk, beta = beta_tbl$beta_season_mean[match(chk, beta_tbl$OTU_ID)],
                 raw_diff = round(raw_diff[chk], 3)))
stopifnot(cor(beta_tbl$beta_season_mean, raw_diff[beta_tbl$OTU_ID]) > 0.8)

# Complement: replicate-level model, and its agreement with the headline.
beta_rep <- beta_season_table(m_rep)
beta_rep <- beta_rep[order(-beta_rep$beta_season_mean), ]
write_tab(beta_rep, "hmsc_rep_beta_season.csv")
k <- merge(beta_tbl, beta_rep, by = "OTU_ID", suffixes = c(".main", ".rep"))
res_main <- k$beta_season_CrI2.5.main > 0 | k$beta_season_CrI97.5.main < 0
res_rep  <- k$beta_season_CrI2.5.rep  > 0 | k$beta_season_CrI97.5.rep  < 0
rep_vs_main <- data.frame(
  n_otu = nrow(k),
  pearson  = round(cor(k$beta_season_mean.main, k$beta_season_mean.rep), 3),
  spearman = round(cor(k$beta_season_mean.main, k$beta_season_mean.rep, method = "spearman"), 3),
  median_CrI_width_ratio_rep_over_main = round(median(
    (k$beta_season_CrI97.5.rep - k$beta_season_CrI2.5.rep) /
    (k$beta_season_CrI97.5.main - k$beta_season_CrI2.5.main)), 3),
  resolved95_main = sum(res_main), resolved95_rep = sum(res_rep),
  resolved95_both = sum(res_main & res_rep))
write_tab(rep_vs_main, "hmsc_rep_vs_main.csv")
cat("Replicate-level vs headline Season Beta:\n"); print(rep_vs_main)

# Caterpillar of the strongest per-OTU seasonal niches (top+bottom 25 by mean).
bshow <- unique(rbind(head(beta_tbl, 25), tail(beta_tbl, 25)))
bshow$lab <- ifelse(is.na(bshow$Genus), bshow$OTU_ID, paste0(bshow$OTU_ID, " (", bshow$Genus, ")"))
bshow$lab <- factor(bshow$lab, levels = bshow$lab[order(bshow$beta_season_mean)])
p_beta <- ggplot(bshow, aes(beta_season_mean, lab, colour = guild)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_segment(aes(x = beta_season_CrI2.5, xend = beta_season_CrI97.5, y = lab, yend = lab),
               alpha = 0.6) +
  geom_point() +
  labs(x = "HMSC Beta: Season (summer vs winter, CLR)", y = NULL,
       title = "Per-OTU seasonal niche (top + bottom 25 OTUs)",
       colour = "Guild") +
  theme_bw(base_size = 10)
save_png(file.path(plot_dir, "hmsc_beta_season.png"), p_beta, width = 8, height = 9, res = 800)

# -----------------------------------------------------------------------------
# SECTION 6.3b -- the SAME per-OTU Season niche from the PROBIT model
# Mirrors 6.3 (same columns, same CrI parse), on m_pa instead of m_clr. The
# occurrence-scale counterpart to the abundance niche; supplementary.
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
cat(sprintf("Headline CLR Beta(Season): %d of %d with a 95%% CrI excluding 0; Spearman vs probit %.3f\n",
            sum(beta_tbl$beta_season_CrI2.5 > 0 | beta_tbl$beta_season_CrI97.5 < 0), nrow(beta_tbl),
            cor(beta_tbl$beta_season_mean, beta_pa$beta_season_mean[match(beta_tbl$OTU_ID, beta_pa$OTU_ID)],
                method = "spearman")))

gllvm_path <- file.path(CANON_MODELS, "gllvm_perOTU_season_coef.csv")
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
s_idx  <- match("sample", names(m_clr$ranLevels)) # biological-sample level = residual co-occurrence
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
# fold at its stored MCMC settings, so CV is run only for the headline CLR
# model (55 rows: ~30 min). Explanatory R2 is reported for CLR and the P/A.
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
# 2-fold CV on CLR, folds split by biological sample (one row per dropping here).
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
  geB <- getPostEstimate(m_clr, parName = "Beta"); geG <- getPostEstimate(m_clr, parName = "Gamma")
  season_row <- which(xcn2 == "Seasonsummer")
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
} else cat("Phylo variant not run (HMSC_RUN_PHYLO != 1); hmsc_rho_phylo.csv left as is.\n")

# =============================================================================
# SECTION 6.7 -- PCR-REPLICATE VARIANCE, model-free
# =============================================================================
# The headline model averages the two PCR replicates of each dropping, so it
# has no PCR level to partition. The replicate-to-replicate variance is instead
# measured directly from the 54 droppings with two replicates. Both replicates
# come from the same DNA extract, so their difference is measurement noise.
# Per OTU:
#   share_total  = [mean(d^2)/2] / var(CLR)       d = CLR(rep1) - CLR(rep2);
#                  the /2 because d carries the noise twice.   (Fig 5B column)
#   share_resid  = same numerator / var of the residual after Season + Year
#   share_lme4   = residual / (residual + dropping) variance from
#                  lmer(CLR ~ Season + Year + (1 | dropping)) -- a model-based
#                  check on the same quantity
#   discordance  = among droppings where the OTU is detected in >= 1 replicate,
#                  the fraction where it is detected in only one
#   share_both   = share_total computed only on droppings where BOTH replicates
#                  detect the OTU (abundance noise, no detection dropout)
# Group summaries (all OTUs, by guild) carry a bootstrap 95% CI over droppings.
pairs <- split(seq_len(nrow(Yclr_rep)), studyDesign$sample)
pairs <- pairs[lengths(pairs) == 2]
i1 <- vapply(pairs, `[`, integer(1), 1); i2 <- vapply(pairs, `[`, integer(1), 2)
det <- Yk > 0
share_total_of <- function(a1, a2) vapply(seq_len(ncol(Yclr_rep)), function(j) {
  y1 <- Yclr_rep[a1, j]; y2 <- Yclr_rep[a2, j]
  (mean((y1 - y2)^2) / 2) / var(c(y1, y2))
}, numeric(1))
share_total <- share_total_of(i1, i2)
share_resid <- vapply(seq_len(ncol(Yclr_rep)), function(j) {
  y <- Yclr_rep[c(i1, i2), j]; r <- resid(lm(y ~ Season + Year, data = XData[c(i1, i2), ]))
  (mean((Yclr_rep[i1, j] - Yclr_rep[i2, j])^2) / 2) / var(r)
}, numeric(1))
share_lme4 <- vapply(seq_len(ncol(Yclr_rep)), function(j) {
  d <- data.frame(y = Yclr_rep[c(i1, i2), j], XData[c(i1, i2), ],
                  drop = factor(studyDesign$sample[c(i1, i2)]))
  f <- suppressMessages(suppressWarnings(lmer(y ~ Season + Year + (1 | drop), data = d)))
  vc <- as.data.frame(VarCorr(f)); vc$vcov[vc$grp == "Residual"] / sum(vc$vcov)
}, numeric(1))
pcr_otu <- do.call(rbind, lapply(seq_len(ncol(Yclr_rep)), function(j) {
  p1 <- det[i1, j]; p2 <- det[i2, j]; any1 <- p1 | p2; both <- p1 & p2
  sb <- if (sum(both) >= 5) {
    a <- Yclr_rep[i1[both], j]; b <- Yclr_rep[i2[both], j]; (mean((a - b)^2) / 2) / var(c(a, b))
  } else NA_real_
  data.frame(n_pairs_detected = sum(any1), n_pairs_both = sum(both),
             discordance = if (sum(any1)) sum(any1 & !both) / sum(any1) else NA_real_,
             share_both = sb)
}))
pcr_otu <- data.frame(OTU_ID = otu_ids, Genus = otu_gen, guild = guild_class,
                      share_total = round(share_total, 4), share_resid = round(share_resid, 4),
                      share_lme4 = round(share_lme4, 4), pcr_otu, stringsAsFactors = FALSE)
pcr_otu$discordance <- round(pcr_otu$discordance, 4); pcr_otu$share_both <- round(pcr_otu$share_both, 4)
write_tab(pcr_otu, "hmsc_pcr_replicate_variance_perOTU.csv")

grp_lvls <- c("all", guild_lvls)
grp_of   <- function(g) if (g == "all") rep(TRUE, length(otu_ids)) else guild_class == g
set.seed(20260929)
boot <- replicate(500, {
  kk <- sample(length(pairs), replace = TRUE)
  st <- share_total_of(i1[kk], i2[kk])
  vapply(grp_lvls, function(g) mean(st[grp_of(g)]), numeric(1))
})
pcr_grp <- do.call(rbind, lapply(grp_lvls, function(g) {
  w <- grp_of(g); ci <- quantile(boot[g, ], c(0.025, 0.975))
  data.frame(group = g, n_otu = sum(w),
             share_total_mean = round(mean(share_total[w]), 3),
             share_total_median = round(median(share_total[w]), 3),
             boot_CrI2.5 = round(ci[[1]], 3), boot_CrI97.5 = round(ci[[2]], 3),
             share_resid_mean = round(mean(share_resid[w]), 3),
             share_lme4_mean = round(mean(share_lme4[w]), 3),
             discordance_mean = round(mean(pcr_otu$discordance[w], na.rm = TRUE), 3),
             share_both_median = round(median(pcr_otu$share_both[w], na.rm = TRUE), 3),
             n_otu_share_both = sum(!is.na(pcr_otu$share_both[w])),
             stringsAsFactors = FALSE)
}))
write_tab(pcr_grp, "hmsc_pcr_replicate_variance.csv")
cat(sprintf("PCR-replicate variance (model-free), %d paired droppings:\n", length(pairs)))
print(pcr_grp)
cat(sprintf("share_resid vs share_lme4 per OTU: r = %.3f\n", cor(share_resid, share_lme4)))

# =============================================================================
# SECTION 7 -- STAGE FIGURES/TABLES INTO Supplementary
# =============================================================================
if (dir.exists(SUPP_DIR)) {
  dir.create(supp_fig, showWarnings = FALSE, recursive = TRUE)
  dir.create(supp_tab, showWarnings = FALSE, recursive = TRUE)
  figs <- c("hmsc_variance_partition.png","hmsc_gamma.png","hmsc_probit_gamma.png",
            "hmsc_beta_season.png","hmsc_vs_gllvm_season.png","hmsc_omega_sample.png",
            "hmsc_R2_explanatory_vs_cv.png")
  tabs <- c("hmsc_otu_guild.csv","hmsc_convergence_psrf.csv","hmsc_variance_partition.csv",
            "hmsc_rep_variance_partition.csv","hmsc_rep_vp_by_chain.csv","hmsc_rep_beta_season.csv",
            "hmsc_rep_vs_main.csv","hmsc_gamma_CrI.csv","hmsc_probit_gamma_CrI.csv",
            "hmsc_beta_season.csv","hmsc_probit_beta_season.csv","hmsc_vs_gllvm_season.csv",
            "hmsc_predictive_R2.csv","hmsc_predictive_R2_summary.csv",
            "hmsc_pcr_replicate_variance.csv","hmsc_pcr_replicate_variance_perOTU.csv",
            "hmsc_rho_phylo.csv")
  figs <- figs[file.exists(file.path(plot_dir, figs))]
  tabs <- tabs[file.exists(file.path(out_dir,  tabs))]
  invisible(file.copy(file.path(plot_dir, figs), supp_fig, overwrite = TRUE))
  invisible(file.copy(file.path(out_dir,  tabs), supp_tab, overwrite = TRUE))
  cat("Staged HMSC figures/tables into", SUPP_DIR, "\n")
}
cat("10_hmsc.R complete.\n")
