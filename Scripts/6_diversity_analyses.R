# =============================================================================
# 6_diversity_analyses.R
#
# PART A -- H2 (alpha diversity; preregistered):
#   Hill numbers (q=0,1,2) via coverage-based rarefaction/extrapolation
#   (iNEXT3D) and per-sample values (hillR). Full-sample main effect tested
#   frequentist (mixed model, Season + (1|Year) + (1|indivID)); diet-richness
#   mechanism tested on the plant-matched subset, Bayesian (brms).
#
# PART B -- Phylogenetic diversity & community structure (supplementary; not
#   part of the core H2 preregistration -- kept for completeness/possible
#   later use): ML tree construction, iNEXT3D meanPD, NRI/NTI, beta-NTI, and
#   UniFrac PERMANOVA.
#
# Requires objects from 4_data_prep.R (eco_analysis.RData):
#   alldat, alldat.rfy, alldat_full  (each: list nopool/pool/pspool)
# Requires: /home/daniel/Ptarmigan/plant_ITS/phyloseq_plant_ITS.rds
#   (plant ITS2 phyloseq, samples keyed by Sample_ID_field / RL_#### -- see
#   4_data_prep.R / CLAUDE.md for the join-key gotcha: raw RL_#### string,
#   no suffix stripping)
#
# nopool is the primary strategy throughout (matches scripts 4/5/7).
#
# Reference: 6_diversity_analyses.R @ 65fbffa (Root_fungi_DADA2) -- an Andean
# root study with a forest/subparamo/paramo habitat gradient across 3 sites.
# This study has no habitat/site gradient; the design is Season (winter/
# summer) x Year (2022-2024) x indivID (repeated individuals) -- see H1/H2 in
# project notes. Reference `habitat`/`site` groupings are mapped onto
# Season/Year in Part B below, same mapping used in scripts 5 and 8.
# =============================================================================

library(phyloseq)
library(DECIPHER)
library(phangorn)
library(iNEXT.3D)
library(hillR)
library(lme4)
library(lmerTest)
library(picante)
library(vegan)
library(ape)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(writexl)

setwd("/home/daniel/Ptarmigan/trimmed/mergedPlates/")
load("eco_analysis.RData")

# Absolute output dirs (independent of the data working dir above), matching
# the convention set in 5_community_composition.R.
out_dir  <- "/home/daniel/Ptarmigan/models/"
plot_dir <- "/home/daniel/Ptarmigan/plots/"
dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Plotting strategy: always write PNG via an explicit png() device -------
# Remote/headless server -- see 5_community_composition.R for the full
# rationale (standing null device + explicit png()/dev.off(), never bare
# ggsave() or an interactive graphics window).
grDevices::pdf(NULL)
save_png <- function(plot_obj, filename, width, height, dpi = 800, units = "in") {
  png(file.path(plot_dir, filename), width = width, height = height, units = units, res = dpi)
  print(plot_obj)
  dev.off()
}

# ---- Shared helper: phyloseq otu_table -> samples x taxa matrix -------------
otu_mat_of <- function(ps) {
  m <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) m <- t(m)
  m
}

# ---- Shared helper: indivID with unidentified samples as unique singletons --
# indivID is NA for samples lacking a confirmed microsatellite match; giving
# each a unique placeholder makes it an unshared singleton random intercept
# instead of silently pooling unrelated unidentified individuals together.
# Same convention as indivID_glm in 7_gllvm.R -- reused here for consistency.
add_indivID_glm <- function(md) {
  md$indivID_glm <- ifelse(is.na(md$indivID) | md$indivID == "",
                            paste0("unk_", md$sample), md$indivID)
  md$indivID_glm <- droplevels(factor(md$indivID_glm))
  md
}


# #############################################################################
# PART A -- H2: ALPHA DIVERSITY (preregistered)
# #############################################################################

# =============================================================================
# SECTION 1 -- iNEXT3D: COVERAGE-BASED TAXONOMIC HILL DIVERSITY (q = 0, 1, 2)
# Standardised diversity comparison between timeblocks (winter vs summer).
# =============================================================================

ps <- alldat$nopool

stopifnot(inherits(ps, "phyloseq"))
sd  <- as.data.frame(sample_data(ps))
OTU <- as(otu_table(ps), "matrix");  if (!taxa_are_rows(ps)) OTU <- t(OTU)

# Build incidence matrices (species x sampling-units) per Season
inc_by_season_TD <- lapply(split(rownames(sd), sd$Season), function(samps){
  M <- OTU[, colnames(OTU) %in% samps, drop=FALSE]
  M[M > 0] <- 1L;  storage.mode(M) <- "numeric";  M
})

set.seed(1)
out_TD <- iNEXT3D(data=inc_by_season_TD, diversity="TD", q=c(0,1,2),
                  datatype="incidence_raw", nboot=500)

# Relabel Seasons for plotting
labs_map  <- c(winter="Winter", summer="Summer")
lvl_order <- c("Winter","Summer")

.relabel <- function(df){
  if (!("Assemblage" %in% names(df))) return(df)
  df$Assemblage <- factor(labs_map[df$Assemblage], levels=lvl_order)
  df
}
out_TD$TDInfo   <- .relabel(as.data.frame(out_TD$TDInfo))
out_TD$TDAsyEst <- .relabel(as.data.frame(out_TD$TDAsyEst))
if (!is.null(out_TD$TDiNextEst))
  out_TD$TDiNextEst <- lapply(out_TD$TDiNextEst, function(d) .relabel(as.data.frame(d)))

write.csv(as.data.frame(out_TD$TDAsyEst), file.path(out_dir, "H2_iNEXT3D_TD_asymptotic.csv"),
          row.names=FALSE)

p_TD_1 <- ggiNEXT3D(out_TD, type=1, facet.var="Order.q") +
  labs(x="Sampling units") +
  theme_minimal(base_size=14) +
  theme(plot.title=element_blank(), legend.title=element_blank())

p_TD_3 <- ggiNEXT3D(out_TD, type=3, facet.var="Order.q") +
  labs(x="Sampling units") +
  theme_minimal(base_size=14) +
  theme(plot.title=element_blank(), legend.title=element_blank())

save_png(p_TD_1, "H2_iNEXT_TD_Season_sampleSize.png", width=10, height=6, dpi=800)
save_png(p_TD_3, "H2_iNEXT_TD_by_Season_coverage.png", width=10, height=6, dpi=800)

# ---- Coverage-standardised comparison at shared minimum coverage -----------
# "standardised diversity comparison between timeblocks" (prereg wording):
# rarefy/extrapolate both Seasons to the coverage level of the least-sampled
# Season (Cmin) so TD estimates are compared at equal sample completeness,
# not equal sample size.
Cmin <- min(as.data.frame(out_TD$TDInfo)[["SC(T)"]])
est_TD <- estimate3D(inc_by_season_TD, diversity="TD", q=c(0,1,2),
                     datatype="incidence_raw", base="coverage",
                     level=Cmin, nboot=500)
write.csv(est_TD, file.path(out_dir, "H2_iNEXT3D_TD_coverage_standardised.csv"),
          row.names=FALSE)
cat(sprintf("\nCoverage-standardised TD at shared coverage Cmin=%.4f:\n", Cmin))
print(est_TD)

# =============================================================================
# SECTION 2 -- PER-PCR-REPLICATE HILL DIVERSITY (hillR) + FREQUENTIST MIXED
# MODEL, WITH PCR REPLICATE AS AN EXPLICIT RANDOM EFFECT
# H2 full-sample main effect: diversity ~ Season + (1|Year) + (1|indivID) +
# (1|PCR replicate)
#
# Uses alldat_full[[1]] (== alldat_full$nopool) -- one row per PCR replicate
# (114 rows = 57 dung samples x 2 replicates), NOT the PCR-collapsed
# alldat/alldat.rfy. The two technical replicates of a dung sample are kept
# as separate rows and their pseudoreplication is handled explicitly via a
# (1|pcr_replicate_id) random intercept (pcr_replicate_id groups the 2 reps
# of the same sample via the native `sample` column) -- same convention as
# 7_gllvm.R's row.eff=~(1|indivID_glm)+(1|pcr_replicate_id), just applied to
# a Hill-diversity outcome instead of a GLLVM count model.
#
# "timeblock" in the preregistration = the clean 2-level winter/summer
# `Season` factor (see 7_gllvm.R for the same Season-vs-TimeBlock note --
# `TimeBlock` is a finer Season x Year-block code, S_1/S_2/S_3/W_1/W_2/W_3,
# not the H2 fixed effect).
# =============================================================================

# ---- Depth filter + rarefaction at the PCR-replicate level ------------------
# alldat_full is deliberately NOT depth-filtered/rarefied at the object level
# (4_data_prep.R uses a pass-through threshold, min_depth_full=1, since it's
# built for GLLVM/HMSC-style read-count-offset models). Individual PCR
# replicates range from 4 to ~637,000 reads (verified on current data) --
# hillR needs comparable sampling effort across rows, so computing Hill
# numbers with no floor at all would make a 4-read reaction meaningless and
# would crush the shared rarefaction depth for every other row. We apply the
# same canonical min_depth (10,000 reads) already used for `alldat`, just at
# the replicate level here (stricter than the sample-level bar since it's
# applied before the 2 replicates are combined) -- change PCR_MIN_DEPTH below
# for a looser/tighter bar; ~16/114 replicates currently fall below 10,000.
PCR_MIN_DEPTH <- 10000

ps_pcr <- prune_samples(sample_sums(alldat_full[[1]]) >= PCR_MIN_DEPTH, alldat_full[[1]])
ps_pcr <- prune_taxa(taxa_sums(ps_pcr) > 0, ps_pcr)
cat(sprintf("\nPCR replicates retained at >= %d reads: %d of %d\n",
            PCR_MIN_DEPTH, nsamples(ps_pcr), nsamples(alldat_full[[1]])))

set.seed(1)
rar_depth_pcr <- min(sample_sums(ps_pcr))
ps_pcr_rar <- rarefy_even_depth(ps_pcr, sample.size=rar_depth_pcr,
                                replace=FALSE, rngseed=1, verbose=FALSE)

m_pcr  <- otu_mat_of(ps_pcr_rar)
md_pcr <- as.data.frame(sample_data(ps_pcr_rar))
stopifnot(identical(rownames(m_pcr), rownames(md_pcr)))

md_pcr$Season <- factor(md_pcr$Season, levels=c("winter","summer"))
md_pcr$Year   <- factor(md_pcr$Year)
md_pcr <- add_indivID_glm(md_pcr)
# PCR-replicate random effect: groups the (up to) 2 technical replicates of
# the same dung sample via the native `sample` column.
md_pcr$pcr_replicate_id <- droplevels(factor(md_pcr$sample))

hill_pcr <- data.frame(
  pcr_row     = rownames(m_pcr),
  q0          = hillR::hill_taxa(m_pcr, q=0),
  q1          = hillR::hill_taxa(m_pcr, q=1),
  q2          = hillR::hill_taxa(m_pcr, q=2),
  Season      = md_pcr$Season,
  Year        = md_pcr$Year,
  indivID_glm = md_pcr$indivID_glm,
  pcr_replicate_id = md_pcr$pcr_replicate_id,
  stringsAsFactors = FALSE
)
write.csv(hill_pcr, file.path(out_dir, "H2_hill_per_pcr_replicate.csv"), row.names=FALSE)
cat(sprintf("\nPer-PCR-replicate Hill diversity computed at rarefaction depth = %d reads (%d rows, %d distinct samples)\n",
            rar_depth_pcr, nrow(hill_pcr), nlevels(hill_pcr$pcr_replicate_id)))

hl_pcr <- hill_pcr %>% pivot_longer(c(q0,q1,q2), names_to="q", values_to="value")
p_hill_box <- ggplot(hl_pcr, aes(Season, value)) +
  geom_boxplot(outlier.shape=NA, fill="grey92") +
  geom_jitter(aes(colour=Year), width=0.15, size=1.6, alpha=0.7) +
  facet_wrap(~q, scales="free_y",
             labeller=as_labeller(c(q0="q0 richness", q1="q1 Shannon", q2="q2 Simpson"))) +
  labs(title=sprintf("Per-PCR-replicate Hill diversity (rarefied to %d reads)", rar_depth_pcr),
       y="Hill number", x=NULL) +
  theme_bw(base_size=12)
save_png(p_hill_box, "H2_hill_boxplot_Season_pcr.png", width=10, height=4.5, dpi=800)

# ---- Frequentist mixed model, one per Hill order q = 0, 1, 2 ---------------
# Fallback ladder if the full random-effects structure is singular: drop
# (1|indivID_glm) first (most likely culprit -- many NA-individual singleton
# placeholders), then (1|Year), then plain lm(). (1|pcr_replicate_id) is kept
# as long as possible since it's the effect this section was changed to add.
# Reports whichever model actually converged, flagged in `re_structure`.
fit_h2_lmm <- function(dat, qcol) {
  f_full <- as.formula(sprintf("%s ~ Season + (1|Year) + (1|indivID_glm) + (1|pcr_replicate_id)", qcol))
  fit <- tryCatch(lmerTest::lmer(f_full, data=dat), error=function(e) NULL)
  re_structure <- "Season + (1|Year) + (1|indivID_glm) + (1|pcr_replicate_id)"

  if (is.null(fit) || isSingular(fit)) {
    f2 <- as.formula(sprintf("%s ~ Season + (1|Year) + (1|pcr_replicate_id)", qcol))
    fit <- tryCatch(lmerTest::lmer(f2, data=dat), error=function(e) NULL)
    re_structure <- "Season + (1|Year) + (1|pcr_replicate_id)"
  }
  if (is.null(fit) || isSingular(fit)) {
    f3 <- as.formula(sprintf("%s ~ Season + (1|pcr_replicate_id)", qcol))
    fit <- tryCatch(lmerTest::lmer(f3, data=dat), error=function(e) NULL)
    re_structure <- "Season + (1|pcr_replicate_id)"
  }
  if (is.null(fit) || isSingular(fit)) {
    fit <- lm(as.formula(sprintf("%s ~ Season", qcol)), data=dat)
    re_structure <- "Season (lm, no random effects)"
  }

  co  <- summary(fit)$coefficients
  est <- co["Seasonsummer","Estimate"]
  se  <- co["Seasonsummer","Std. Error"]
  pv  <- if ("Pr(>|t|)" %in% colnames(co)) co["Seasonsummer","Pr(>|t|)"] else co["Seasonsummer","Pr(>|z|)"]

  list(
    summary = data.frame(
      metric = qcol, re_structure = re_structure,
      summer_minus_winter = round(est,3),
      lwr = round(est-1.96*se,3), upr = round(est+1.96*se,3),
      p = round(pv,4)
    ),
    fit = fit
  )
}

h2_lmm_fits <- lapply(c("q0","q1","q2"), function(q) fit_h2_lmm(hill_pcr, q))
names(h2_lmm_fits) <- c("q0","q1","q2")
h2_lmm_tab <- dplyr::bind_rows(lapply(h2_lmm_fits, `[[`, "summary"))

write.csv(h2_lmm_tab, file.path(out_dir, "H2_lmm_season_effect.csv"), row.names=FALSE)
cat("\nH2 full-sample frequentist mixed model (Season effect, summer vs winter ref):\n")
print(h2_lmm_tab)
cat("\nq1 (Shannon) model detail:\n")
print(summary(h2_lmm_fits$q1$fit))

# =============================================================================
# SECTION 3 -- H2 MECHANISM: BAYESIAN DIET-RICHNESS MODEL (matched subset)
# Fungal Hill diversity ~ plant (diet) richness + Season + Year (Year fixed;
# see the note at the model call for why the preregistered REs were dropped).
# H2-mechanism supported if P(slope_diet_richness > 0) > 0.95.
#
# Self-contained from here: plant ITS2 diet data is one row per biological
# dung sample (not per PCR replicate), so this section needs the PCR-
# collapsed fungal object -- it builds its own `ps_rar_sample`/`hill_sample`
# from alldat.rfy$nopool rather than reusing Section 2's now PCR-replicate-
# level `ps_pcr_rar`/`hill_pcr`.
# =============================================================================

# hillR::hill_taxa requires non-negative abundances at even depth for a fair
# richness/evenness comparison across samples -- use the already-rarefied,
# PCR-collapsed canonical object (alldat.rfy), not alldat (raw counts) or
# Section 2's per-replicate object.
ps_rar_sample <- alldat.rfy$nopool
m_rar_sample  <- otu_mat_of(ps_rar_sample)
md_rar_sample <- as.data.frame(sample_data(ps_rar_sample))
stopifnot(identical(rownames(m_rar_sample), rownames(md_rar_sample)))

md_rar_sample$Season <- factor(md_rar_sample$Season, levels=c("winter","summer"))
md_rar_sample$Year   <- factor(md_rar_sample$Year)
md_rar_sample <- add_indivID_glm(md_rar_sample)

hill_sample <- data.frame(
  sample  = rownames(m_rar_sample),
  q0      = hillR::hill_taxa(m_rar_sample, q=0),
  q1      = hillR::hill_taxa(m_rar_sample, q=1),
  q2      = hillR::hill_taxa(m_rar_sample, q=2),
  Season  = md_rar_sample$Season,
  Year    = md_rar_sample$Year,
  indivID_glm = md_rar_sample$indivID_glm,
  stringsAsFactors = FALSE
)
write.csv(hill_sample, file.path(out_dir, "H2_hill_per_sample.csv"), row.names=FALSE)
cat(sprintf("\nPer-sample (PCR-collapsed) Hill diversity computed at rarefaction depth = %d reads (%d samples)\n",
            min(rowSums(m_rar_sample)), nrow(hill_sample)))

# ---- Build the plant-matched subset -----------------------------------------
# Matched = canonical fungal samples (ps_rar_sample, PCR-collapsed + rarefied)
# that also have a plant ITS2 sample, joined on the raw Sample_ID_field key
# (RL_####, no suffix stripping -- see 4_data_prep.R / CLAUDE.md gotcha).
# This is the literal "matched subset" of the preregistration (~31 samples
# per the project design note); it does NOT additionally require
# microsatellite-confirmed species ID.
#
# OPTIONAL STRICTER VARIANT (not used by default): restrict further to
# microsatellite-confirmed L. lagopus only, by uncommenting the line below.
# This trades sample size (~31 -> high-20s) for species-identity certainty;
# an earlier exploratory pass (working/Fable/confirmatory_analysis.R) used
# this stricter variant as its primary analysis set.
plant <- readRDS("/home/daniel/Ptarmigan/plant_ITS/phyloseq_plant_ITS.rds")

md_all <- as.data.frame(sample_data(ps_rar_sample))
matched_ids <- rownames(md_all)[md_all$Sample_ID_field %in% sample_names(plant)]
# matched_ids <- rownames(md_all)[md_all$Sample_ID_field %in% sample_names(plant) &
#                                    !is.na(md_all$NINAmicrosat) & md_all$NINAmicrosat == "success"]

ps_match <- prune_samples(matched_ids, ps_rar_sample)
ps_match <- prune_taxa(taxa_sums(ps_match) > 0, ps_match)

md_match     <- as.data.frame(sample_data(ps_match))
plant_match  <- prune_samples(md_match$Sample_ID_field, plant)
plant_match  <- prune_taxa(taxa_sums(plant_match) > 0, plant_match)

cat(sprintf("\nH2-mechanism matched subset: %d fungal samples with paired plant ITS2 data\n",
            nsamples(ps_match)))

# ---- Per-sample plant (diet) richness ---------------------------------------
# Rarefy plant OTU table to its own minimum depth (independent of the fungal
# rarefaction depth above) for a fair richness comparison across matched
# samples.
plant_otu   <- otu_mat_of(plant_match)
plant_rar   <- rarefy_even_depth(plant_match, sample.size=min(rowSums(plant_otu)),
                                 replace=FALSE, rngseed=1, verbose=FALSE)
plant_rich  <- rowSums(otu_mat_of(plant_rar) > 0)   # names = Sample_ID_field

# ---- Assemble model data frame ----------------------------------------------
# Reuse Section 2's fungal Hill values (same rarefaction) rather than
# recomputing -- PRIMARY_Q selects which Hill order is modelled as the
# mechanism outcome (q1 = Shannon-equivalent, the preregistration's default
# "fungal Hill diversity"; switch to "q0" or "q2" to test richness/Simpson
# instead).
PRIMARY_Q <- "q1"

diet <- data.frame(
  sample          = rownames(md_match),
  fungal_hill     = hill_sample[[PRIMARY_Q]][match(rownames(md_match), hill_sample$sample)],
  plant_richness  = plant_rich[md_match$Sample_ID_field],
  Season          = factor(md_match$Season, levels=c("winter","summer")),
  Year            = factor(md_match$Year),
  indivID         = md_match$indivID,
  stringsAsFactors = FALSE
)
diet <- add_indivID_glm(diet)
diet <- diet[complete.cases(diet[, c("fungal_hill","plant_richness")]), ]

diet$fungal_hill_z    <- as.numeric(scale(diet$fungal_hill))
diet$plant_richness_z <- as.numeric(scale(diet$plant_richness))

write.csv(diet, file.path(out_dir, "H2_diet_richness_matched_data.csv"), row.names=FALSE)
cat(sprintf("Bayesian diet-richness model n = %d (metric: fungal %s)\n", nrow(diet), PRIMARY_Q))

# ---- Bayesian regression (brms) ---------------------------------------------
# TUNABLE PARAMETERS (documented inline so they're easy to revisit):
#   - PRIMARY_Q       : which Hill order is the outcome (set above; "q0"/"q1"/"q2").
#   - priors below     : weakly informative, standardised-slope priors per the
#                         preregistration ("Normal(0, 1) on standardised slopes").
#                         Tighten (e.g. normal(0, 0.5)) for a more skeptical
#                         prior, or widen for a more agnostic one.
#   - family          : gaussian() on z-scored Hill values. If diagnostics
#                         (posterior predictive check / residuals) look poor,
#                         consider modelling fungal_hill directly (unscaled)
#                         with a lognormal() or Gamma(link="log") family
#                         instead of z-scoring + gaussian().
#   - random effects  : NONE. The preregistration named (1|Year) + (1|indivID),
#                         but on this matched subset (n=30) that model does NOT
#                         converge -- 477 divergent transitions, max Rhat 1.12,
#                         RE-SD Bulk_ESS ~50. Cause: 27 of 30 samples are unique
#                         individuals (indivID is unidentifiable, one obs/group)
#                         and Year has only 3 levels (RE-SD poorly estimated).
#                         Year is therefore modelled as a FIXED effect and the
#                         individual RE dropped. Verified (scratch h2_compare.R):
#                         this simplified model converges cleanly (0 divergences,
#                         Rhat 1.001) and LOO is statistically indistinguishable
#                         from the RE models (elpd_diff < 1.6, se_diff ~0.8).
#   - chains/iter/warmup/adapt_delta/cores : standard brms sampler controls.
#   - decision threshold : the preregistration's support criterion is
#                         P(slope_plant_richness > 0) > 0.95 -- change
#                         `support_threshold` below to explore sensitivity.
library(brms)

h2_prior <- c(
  brms::set_prior("normal(0,1)",   class="b"),
  brms::set_prior("normal(0,0.5)", class="Intercept")
)

b_h2 <- brms::brm(
  fungal_hill_z ~ plant_richness_z + Season + Year,
  data    = diet,
  family  = gaussian(),
  prior   = h2_prior,
  chains  = 4, iter = 6000, warmup = 3000, cores = 4,
  control = list(adapt_delta = 0.999, max_treedepth = 15),
  seed    = 1, refresh = 0,
  file    = file.path(out_dir, "H2_brms_diet_richness")
)

draws_slope <- brms::as_draws_df(b_h2)$b_plant_richness_z
support_threshold <- 0.95

h2_mechanism_summary <- data.frame(
  metric   = PRIMARY_Q,
  n        = nrow(diet),
  median   = median(draws_slope),
  lwr95    = unname(quantile(draws_slope, .025)),
  upr95    = unname(quantile(draws_slope, .975)),
  P_gt0    = mean(draws_slope > 0),
  supported = mean(draws_slope > 0) > support_threshold
)
write.csv(h2_mechanism_summary, file.path(out_dir, "H2_mechanism_brms_summary.csv"), row.names=FALSE)

cat("\nH2-mechanism (diet-richness) posterior summary:\n")
print(h2_mechanism_summary)
cat(sprintf("\nH2-mechanism %s (P[slope>0] = %.3f, threshold = %.2f)\n",
            ifelse(h2_mechanism_summary$supported, "SUPPORTED", "NOT supported"),
            h2_mechanism_summary$P_gt0, support_threshold))

# ---- Convergence diagnostics -------------------------------------------------
rh <- brms::rhat(b_h2); es <- brms::neff_ratio(b_h2)
cat(sprintf("Max Rhat = %.4f | Min Neff ratio = %.3f (want Rhat<1.01, Neff ratio not too close to 0)\n",
            max(rh, na.rm=TRUE), min(es, na.rm=TRUE)))

# ---- Posterior predictive check ---------------------------------------------
png(file.path(plot_dir, "H2_brms_diet_richness_ppc.png"), width=900, height=600, res=150)
print(brms::pp_check(b_h2, ndraws=100) + labs(title="H2 mechanism: posterior predictive check"))
dev.off()

# ---- Fitted relationship on the raw (unscaled) diversity/richness scale ----
nd <- data.frame(
  plant_richness_z = seq(min(diet$plant_richness_z), max(diet$plant_richness_z), length=60),
  Season = "winter"
)
fe <- brms::posterior_epred(b_h2, newdata=nd, re_formula=NA)

hill_mean <- mean(diet$fungal_hill); hill_sd <- sd(diet$fungal_hill)
rich_mean <- mean(diet$plant_richness); rich_sd <- sd(diet$plant_richness)
ribbon <- data.frame(
  plant_richness = nd$plant_richness_z * rich_sd + rich_mean,
  fit = apply(fe, 2, median) * hill_sd + hill_mean,
  lwr = apply(fe, 2, quantile, .025) * hill_sd + hill_mean,
  upr = apply(fe, 2, quantile, .975) * hill_sd + hill_mean
)

p_diet <- ggplot(diet, aes(plant_richness, fungal_hill)) +
  geom_ribbon(data=ribbon, aes(plant_richness, ymin=lwr, ymax=upr), inherit.aes=FALSE, alpha=0.2) +
  geom_line(data=ribbon, aes(plant_richness, fit), inherit.aes=FALSE) +
  geom_point(aes(colour=Season), size=3) +
  scale_colour_manual(values=c(winter="#0072B2", summer="#E69F00")) +
  labs(title=sprintf("Fungal Hill %s vs dietary plant richness (matched n=%d)", PRIMARY_Q, nrow(diet)),
       subtitle=sprintf("plant-richness slope P(>0)=%.2f", h2_mechanism_summary$P_gt0),
       x="Plant OTU richness (rarefied)", y=sprintf("Fungal Hill %s", PRIMARY_Q)) +
  theme_bw(base_size=12)
save_png(p_diet, "H2_diet_richness_fit.png", width=7.5, height=5.5, dpi=800)


# #############################################################################
# PART B -- PHYLOGENETIC DIVERSITY & COMMUNITY STRUCTURE (supplementary)
# May not be used in the final analysis; kept for completeness. Adapted only
# to current object names/paths and the Season/Year design -- internal logic
# preserved from the reference script.
# #############################################################################

# =============================================================================
# SECTION 4 -- PHYLOGENETIC TREE CONSTRUCTION
# Build a maximum-likelihood ITS phylogeny from the full-complexity, no-rg2
# object (alldat_full$nopool -- widest taxon universe available). Already run
# once (see tree.rds); re-runs only if the cached tree is missing.
# =============================================================================

tree_path <- "tree.rds"

if (!file.exists(tree_path)) {
  ps <- alldat_full$nopool

  # ---- Align sequences --------------------------------------------------------
  align <- AlignSeqs(DNAStringSet(refseq(ps)), anchor=NA)

  # ---- Build NJ starting tree, then optimise GTR + Gamma + Invariant sites ---
  phang_align <- phyDat(as(align, "matrix"), type="DNA")
  dm     <- dist.ml(phang_align)
  treeNJ <- NJ(dm)
  fit    <- pml(treeNJ, data=phang_align)
  fitGTR <- update(fit, k=4, inv=0.2)
  fit    <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                      optNni=FALSE, optBf=TRUE, optQ=TRUE, optEdge=TRUE,
                      optRooted=FALSE, rearrangement="stochastic",
                      control=pml.control(trace=0))

  tree_data <- merge_phyloseq(ps, fit$tree)
  saveRDS(tree_data, tree_path)
} else {
  cat("Loading precomputed ML tree from", tree_path, "\n")
  tree_data <- readRDS(tree_path)
}

# =============================================================================
# SECTION 5 -- iNEXT3D: PHYLOGENETIC DIVERSITY (meanPD)
# =============================================================================

ps <- tree_data    # phyloseq object with phylogenetic tree (full-complexity)

sd <- as.data.frame(sample_data(ps))
sd$Season <- factor(sd$Season, levels=c("winter","summer"))
sample_data(ps)$Season <- sd$Season

print(table(sample_data(ps)$Season, useNA="ifany"))

OTU <- as(otu_table(ps), "matrix");  if (taxa_are_rows(ps)) OTU <- t(OTU)
sd  <- as.data.frame(sample_data(ps))

inc_by_season_PD <- lapply(split(seq_len(nrow(OTU)), sd$Season), function(i){
  M  <- t(OTU[i,, drop=FALSE]);  M[] <- as.integer(M > 0)
  M[rowSums(M) > 0,, drop=FALSE]
})
inc_by_season_PD <- inc_by_season_PD[vapply(inc_by_season_PD, function(M) nrow(M)>0 & ncol(M)>0, TRUE)]

# Prepare tree: prune, deduplicate labels, root
tr    <- phy_tree(ps)
keep  <- intersect(tr$tip.label, unique(unlist(lapply(inc_by_season_PD, rownames))))
stopifnot(length(keep) >= 2)
tr2   <- keep.tip(tr, keep)
if (any(duplicated(tr2$tip.label))) tr2 <- drop.tip(tr2, which(duplicated(tr2$tip.label)))
tr2$node.label <- NULL
if (!is.rooted(tr2)) tr2 <- midpoint(tr2)

tipset          <- tr2$tip.label
inc_by_season_PD <- lapply(inc_by_season_PD, function(M) M[rownames(M) %in% tipset,, drop=FALSE])
inc_by_season_PD <- inc_by_season_PD[vapply(inc_by_season_PD, function(M) nrow(M)>0 & ncol(M)>0, TRUE)]

set.seed(1)
out_PD <- iNEXT3D(data=inc_by_season_PD, diversity="PD", q=c(0,1,2),
                  datatype="incidence_raw", nboot=500,
                  PDtree=tr2, PDtype="meanPD")

write.csv(as.data.frame(out_PD$PDAsyEst), file.path(out_dir, "PhyloSupp_iNEXT3D_PD_asymptotic.csv"),
          row.names=FALSE)

# =============================================================================
# SECTION 6 -- COMMUNITY PHYLOGENETIC STRUCTURE: NRI / NTI
# Standardised effect sizes for MPD (NRI) and MNTD (NTI) via tip shuffling.
# Canonical object per 4_data_prep.R design: alldat (non-rarefied, PCR-
# collapsed) for NRI/NTI.
# =============================================================================

ps <- alldat$nopool

# Attach pruned tree from tree_data (built in Section 4)
tr_full   <- phy_tree(tree_data)
tx_ps     <- taxa_names(ps)
tr_pruned <- keep.tip(tr_full, intersect(tx_ps, tr_full$tip.label))
if (!is.rooted(tr_pruned)) tr_pruned <- midpoint(tr_pruned)
ps <- prune_taxa(taxa_names(ps) %in% tr_pruned$tip.label, ps)
phy_tree(ps) <- tr_pruned

# Samples x taxa matrix aligned to tree
tr2  <- phy_tree(ps)
comm <- as(otu_table(ps), "matrix");  if (taxa_are_rows(ps)) comm <- t(comm)
keep_taxa <- intersect(colnames(comm), tr2$tip.label)
comm <- comm[, keep_taxa, drop=FALSE];  tr2 <- ape::keep.tip(tr2, keep_taxa)

# Drop ultra-rare taxa to stabilise null distributions
keep_taxa2 <- colSums(comm) >= 200
comm <- comm[, keep_taxa2, drop=FALSE];  tr2 <- ape::keep.tip(tr2, colnames(comm))

md <- as(sample_data(ps), "data.frame") %>%
  mutate(Season=factor(Season, levels=c("winter","summer")), Year=factor(Year))

dist_phy <- cophenetic(tr2)

set.seed(1)
mpd_aw  <- ses.mpd (comm, dist_phy, null.model="taxa.labels", runs=999, abundance.weighted=TRUE)
mntd_aw <- ses.mntd(comm, dist_phy, null.model="taxa.labels", runs=999, abundance.weighted=TRUE)
comm_pa <- 1 * (comm > 0)
mpd_pa  <- ses.mpd (comm_pa, dist_phy, null.model="taxa.labels", runs=999, abundance.weighted=FALSE)
mntd_pa <- ses.mntd(comm_pa, dist_phy, null.model="taxa.labels", runs=999, abundance.weighted=FALSE)

# Compile results
md <- as.data.frame(md)[rownames(comm),, drop=FALSE]
out <- data.frame(
  sample  = rownames(comm),
  NRI_aw  = -mpd_aw$mpd.obs.z,  NTI_aw  = -mntd_aw$mntd.obs.z,
  NRI_pa  = -mpd_pa$mpd.obs.z,  NTI_pa  = -mntd_pa$mntd.obs.z,
  stringsAsFactors=FALSE
)
rownames(out) <- out$sample
out <- cbind(out, md[rownames(out), c("Season","Year"), drop=FALSE])
out <- out[complete.cases(out[, c("NRI_aw","NTI_aw","NRI_pa","NTI_pa")]),]
out$Season <- droplevels(factor(out$Season))
out$Year   <- droplevels(factor(out$Year))

print(table(out$Season, useNA="ifany"))

# Season-level summaries
sum_season <- out %>% group_by(Season) %>%
  summarize(n=n(),
            NRI_aw_med=median(NRI_aw,na.rm=TRUE), NRI_aw_IQR=IQR(NRI_aw,na.rm=TRUE),
            NTI_aw_med=median(NTI_aw,na.rm=TRUE), NTI_aw_IQR=IQR(NTI_aw,na.rm=TRUE),
            NRI_pa_med=median(NRI_pa,na.rm=TRUE),  NRI_pa_IQR=IQR(NRI_pa,na.rm=TRUE),
            NTI_pa_med=median(NTI_pa,na.rm=TRUE),  NTI_pa_IQR=IQR(NTI_pa,na.rm=TRUE),
            .groups="drop")
print(sum_season)
write.csv(sum_season, file.path(out_dir, "PhyloSupp_NRI_NTI_season_summary.csv"), row.names=FALSE)

# Linear models: abundance-weighted NRI and NTI
m_NRI_aw <- lm(NRI_aw ~ Season + Year, data=out)
m_NTI_aw <- lm(NTI_aw ~ Season + Year, data=out)
anova(m_NRI_aw);  anova(m_NTI_aw)

# Residual diagnostics
png(file.path(plot_dir, "PhyloSupp_NTI_aw_residuals.png"), width=1800, height=1600, res=200)
par(mfrow=c(2,2), mar=c(4,4,2,1));  plot(m_NTI_aw);  par(mfrow=c(1,1))
dev.off()

# =============================================================================
# SECTION 7 -- BETA-NTI (pairwise phylogenetic turnover)
# =============================================================================

stopifnot(all(rownames(comm) %in% rownames(md)))

compute_betaNTI <- function(comm_mat, tree, abundance_weighted=TRUE, nperm=999, seed=1) {
  set.seed(seed)
  Dphy <- cophenetic(tree)
  beta_mntd_obs  <- comdistnt(comm_mat, Dphy, abundance.weighted=abundance_weighted)
  null_fun <- function(i) {
    tr_null <- tipShuffle(tree)
    comdistnt(comm_mat, cophenetic(tr_null), abundance.weighted=abundance_weighted)
  }
  null_list <- lapply(seq_len(nperm), null_fun)
  mu  <- Reduce(`+`, null_list) / nperm
  mu2 <- Reduce(`+`, lapply(null_list, function(m) m*m)) / nperm
  sdv <- sqrt(pmax(0, mu2 - mu*mu))
  (beta_mntd_obs - mu) / sdv
}

betaNTI_aw <- compute_betaNTI(comm,    tr2, abundance_weighted=TRUE,  nperm=999, seed=1)
betaNTI_pa <- compute_betaNTI(comm_pa, tr2, abundance_weighted=FALSE, nperm=999, seed=1)

# Tidy long format
# comdistnt() returns a "dist" object, not a plain matrix -- row()/col()/
# upper.tri() below require a matrix, so coerce first (as.matrix.dist keeps
# the Labels as dimnames). md_df's native `sample` column already equals its
# rownames (unlike the reference study), so join on it directly instead of
# rownames_to_column("sample"), which would collide with the existing column.
to_long <- function(mat, md_df) {
  mat <- as.matrix(mat)
  nms <- rownames(mat)
  utx <- upper.tri(mat)
  tibble(
    sample_i = nms[row(mat)[utx]],
    sample_j = nms[col(mat)[utx]],
    betaNTI  = mat[utx]
  ) %>%
    left_join(md_df %>% select(sample, Season, Year),
              by=c("sample_i"="sample")) %>%
    left_join(md_df %>% select(sample, Season, Year),
              by=c("sample_j"="sample"), suffix=c("_i","_j")) %>%
    mutate(
      same_Season = Season_i == Season_j,
      same_Year   = Year_i == Year_j,
      season_pair = paste(pmin(as.character(Season_i), as.character(Season_j)),
                          pmax(as.character(Season_i), as.character(Season_j)), sep="-")
    )
}

betaNTI_aw_long <- to_long(betaNTI_aw, md)

sum_beta <- betaNTI_aw_long %>%
  mutate(group=ifelse(same_Season,"within","between")) %>%
  group_by(group) %>%
  summarize(n_pairs=n(), median=median(betaNTI,na.rm=TRUE), IQR=IQR(betaNTI,na.rm=TRUE),
            prop_gt2=mean(betaNTI>2,na.rm=TRUE), prop_lt2=mean(betaNTI<(-2),na.rm=TRUE),
            prop_abs_lt2=mean(abs(betaNTI)<2,na.rm=TRUE), .groups="drop")
print(sum_beta)
write.csv(sum_beta, file.path(out_dir, "PhyloSupp_betaNTI_summary.csv"), row.names=FALSE)

# =============================================================================
# SECTION 8 -- UNIFRAC PERMANOVA ON PHYLOGENETIC BETA-DIVERSITY
# Canonical object per 4_data_prep.R design: alldat.rfy (rarefied) for
# UniFrac -- rebuilt here with the tree attached rather than reusing
# ps/tr2 from Sections 6-7 (those are deliberately on the non-rarefied
# `alldat` object).
# =============================================================================

ps_rfy <- alldat.rfy$nopool
tr_full_rfy   <- phy_tree(tree_data)
tr_pruned_rfy <- keep.tip(tr_full_rfy, intersect(taxa_names(ps_rfy), tr_full_rfy$tip.label))
if (!is.rooted(tr_pruned_rfy)) tr_pruned_rfy <- midpoint(tr_pruned_rfy)

ps_u <- prune_taxa(taxa_names(ps_rfy) %in% tr_pruned_rfy$tip.label, ps_rfy)
tr_pruned_rfy <- keep.tip(tr_pruned_rfy, taxa_names(ps_u))
phy_tree(ps_u) <- tr_pruned_rfy
stopifnot(!is.null(phy_tree(ps_u))); stopifnot(is.rooted(phy_tree(ps_u)))

UFw  <- UniFrac(ps_u, weighted=TRUE,  normalized=TRUE, parallel=FALSE, fast=TRUE)
UFuw <- UniFrac(ps_u, weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE)
md_u <- as(sample_data(ps_u), "data.frame")
md_u$Season <- factor(md_u$Season, levels=c("winter","summer"))
md_u$Year   <- factor(md_u$Year)

set.seed(1)
per_w  <- adonis2(as.dist(UFw)  ~ Season * Year, by="terms", data=md_u, permutations=9999)
per_uw <- adonis2(as.dist(UFuw) ~ Season * Year, by="terms", data=md_u, permutations=9999)

bd_w_season  <- betadisper(as.dist(UFw),  md_u$Season);  bdw_season_a  <- anova(bd_w_season)
bd_uw_season <- betadisper(as.dist(UFuw), md_u$Season);  bduw_season_a <- anova(bd_uw_season)
bd_w_year    <- betadisper(as.dist(UFw),  md_u$Year);    bdw_year_a    <- anova(bd_w_year)
bd_uw_year   <- betadisper(as.dist(UFuw), md_u$Year);    bduw_year_a   <- anova(bd_uw_year)

per_w;  per_uw
bdw_season_a;  bduw_season_a
bdw_year_a;   bduw_year_a

write.csv(as.data.frame(per_w),  file.path(out_dir, "PhyloSupp_UniFrac_weighted_permanova.csv"))
write.csv(as.data.frame(per_uw), file.path(out_dir, "PhyloSupp_UniFrac_unweighted_permanova.csv"))
