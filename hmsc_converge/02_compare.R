# =============================================================================
# hmsc_converge/02_compare.R  (branch exp/hmsc-converge -- experimental)
# Summarise one new fit in the SAME table schema Figure 5 / appendix Section 9
# read, and compare it with the canonical fits.
#
# Usage: conda run -n r_env Rscript hmsc_converge/02_compare.R <tag>
#   e.g. tag = prod_A_thin50  (reads models/hmsc_v2/<tag>.rds + <tag>_meta.rds)
# Writes models/hmsc_v2/<tag>_tables/{hmsc_beta_season, hmsc_variance_partition,
#   hmsc_convergence_psrf, hmsc_gamma_CrI, hmsc_predictive_R2_summary}.csv
#   + <tag>_compare.txt
#
# UNITS. With YScale = TRUE this Hmsc version stores Beta on the per-OTU
# STANDARDISED response scale (verified 2026-09-29: pilot B0 OTU1 intercept
# -0.86 = (winter-2022 mean 4.27 - 6.07) / 2.01). Season/Year Betas are
# therefore multiplied by each OTU's response SD (YScalePar[2, ]) draw by draw,
# so hmsc_beta_season.csv is in response units (CLR or rCLR) like the canonical
# table. Support and variance partition are scale-invariant; Gamma stays on the
# standardised scale (it is a community-level coefficient across OTUs whose
# responses were put on a common scale, which is the point of YScale).
# =============================================================================

suppressMessages({ library(Hmsc); library(coda) })
tag <- commandArgs(trailingOnly = TRUE)[1]
v2  <- "/home/daniel/Ptarmigan/models/hmsc_v2"
m    <- readRDS(file.path(v2, paste0(tag, ".rds")))
meta <- readRDS(file.path(v2, paste0(tag, "_meta.rds")))
tdir <- file.path(v2, paste0(tag, "_tables")); dir.create(tdir, showWarnings = FALSE)
canon <- "/home/daniel/Ptarmigan/models"
sink(file.path(v2, paste0(tag, "_compare.txt")), split = TRUE)

otu   <- colnames(m$Y)
guild <- unname(meta$guild[otu]); genus <- unname(meta$genus[otu])
xcn   <- colnames(m$X); s_row <- which(xcn == "Seasonsummer")
ysd   <- if (!is.null(m$YScalePar)) m$YScalePar[2, ] else rep(1, m$ns)

# ---- Season Beta draws, back-transformed ------------------------------------
post <- poolMcmcChains(m$postList)
Bd <- t(sapply(post, function(s) s$Beta[s_row, ]))            # draws x OTU (scaled)
Bd_resp <- sweep(Bd, 2, ysd, "*")
beta_tbl <- data.frame(
  OTU_ID = otu, Genus = genus, guild = guild,
  beta_season_mean    = round(colMeans(Bd_resp), 3),
  beta_season_support = round(colMeans(Bd > 0), 3),
  beta_season_CrI2.5  = round(apply(Bd_resp, 2, quantile, 0.025), 3),
  beta_season_CrI97.5 = round(apply(Bd_resp, 2, quantile, 0.975), 3),
  beta_season_std_mean = round(colMeans(Bd), 3),
  stringsAsFactors = FALSE)
beta_tbl <- beta_tbl[order(-beta_tbl$beta_season_mean), ]
write.csv(beta_tbl, file.path(tdir, "hmsc_beta_season.csv"), row.names = FALSE)

# ---- Variance partition (same grouping as 10_hmsc.R 6.1) ---------------------
grp <- ifelse(grepl("^Year", xcn), 2L, 1L)
VP  <- computeVariancePartitioning(m, group = grp, groupnames = c("Season", "Year"))
vp_tbl <- data.frame(OTU_ID = colnames(VP$vals), t(round(VP$vals, 4)),
                     stringsAsFactors = FALSE, check.names = FALSE)
vp_tbl$Genus <- genus[match(vp_tbl$OTU_ID, otu)]
vp_tbl$guild <- guild[match(vp_tbl$OTU_ID, otu)]
write.csv(vp_tbl, file.path(tdir, "hmsc_variance_partition.csv"), row.names = FALSE)

# ---- Convergence (psrf_summary from 10_hmsc.R:286-300) ----------------------
mp <- convertToCodaObject(m)
gb <- gelman.diag(mp$Beta, multivariate = FALSE)$psrf[, 1]
gg <- gelman.diag(mp$Gamma, multivariate = FALSE)$psrf[, 1]
conv <- data.frame(model = "clr",    # label the Fig 5 chunk filters on
                   beta_psrf_max = round(max(gb), 3), beta_psrf_med = round(median(gb), 3),
                   beta_pct_below_1.1 = round(100 * mean(gb < 1.1), 1),
                   beta_ess_min = round(min(effectiveSize(mp$Beta))),
                   gamma_psrf_max = round(max(gg), 3), gamma_psrf_med = round(median(gg), 3),
                   gamma_ess_min = round(min(effectiveSize(mp$Gamma))))
write.csv(conv, file.path(tdir, "hmsc_convergence_psrf.csv"), row.names = FALSE)

# ---- Gamma (standardised scale) ---------------------------------------------
ge <- getPostEstimate(m, parName = "Gamma"); gq <- summary(mp$Gamma)$quantiles
gam <- data.frame(covariate = rep(xcn, times = ncol(m$Tr)),
                  trait = rep(colnames(m$Tr), each = length(xcn)),
                  mean = round(as.vector(ge$mean), 3), support = round(as.vector(ge$support), 3),
                  CrI_2.5 = round(gq[, "2.5%"], 3), CrI_97.5 = round(gq[, "97.5%"], 3))
write.csv(gam, file.path(tdir, "hmsc_gamma_CrI.csv"), row.names = FALSE)

# ---- Explanatory R2 ----------------------------------------------------------
mf <- evaluateModelFit(hM = m, predY = computePredictedValues(m))
r2 <- data.frame(model = tag, mean_expl_R2 = round(mean(mf$R2, na.rm = TRUE), 3))
write.csv(r2, file.path(tdir, "hmsc_predictive_R2_summary.csv"), row.names = FALSE)

# ---- Report + comparisons ----------------------------------------------------
cat(sprintf("=== %s: %d OTUs, %d chains x %d samples, thin %d ===\n",
            tag, m$ns, length(m$postList), m$samples, m$thin))
print(conv, row.names = FALSE)
res95  <- sum(beta_tbl$beta_season_CrI2.5 > 0 | beta_tbl$beta_season_CrI97.5 < 0)
supp95 <- sum(beta_tbl$beta_season_support >= 0.95 | beta_tbl$beta_season_support <= 0.05)
cat(sprintf("\nSeason Beta: %d of %d OTUs with 95%% CrI excluding 0; %d at support >= 0.95 or <= 0.05 (%d positive / %d negative)\n",
            res95, m$ns, supp95, sum(beta_tbl$beta_season_support >= 0.95),
            sum(beta_tbl$beta_season_support <= 0.05)))
cat(sprintf("OTUs pointing summer-ward (mean > 0): %d of %d\n", sum(beta_tbl$beta_season_mean > 0), m$ns))

cat("\nGuild-mean Season Beta (response units; std units in brackets):\n")
gmr <- tapply(beta_tbl$beta_season_mean, beta_tbl$guild, mean)
gms <- tapply(beta_tbl$beta_season_std_mean, beta_tbl$guild, mean)
print(data.frame(n = as.vector(table(beta_tbl$guild)[names(gmr)]),
                 mean_resp = round(gmr, 3), mean_std = round(gms, 3)))

vcols <- setdiff(names(vp_tbl), c("OTU_ID", "Genus", "guild"))
cat("\nVariance partition, community mean (%):\n"); print(round(100 * colMeans(vp_tbl[, vcols]), 1))
cat("By guild (%):\n")
print(round(100 * do.call(rbind, lapply(split(vp_tbl[, vcols], vp_tbl$guild), colMeans)), 1))

cat("\nGamma, Season rows (standardised scale):\n"); print(gam[gam$covariate == "Seasonsummer", ], row.names = FALSE)
cat(sprintf("\nMean explanatory R2: %.3f (canonical CLR: 0.455)\n", r2$mean_expl_R2))

old   <- read.csv(file.path(canon, "hmsc_beta_season.csv"), stringsAsFactors = FALSE)
prob  <- read.csv(file.path(canon, "hmsc_probit_beta_season.csv"), stringsAsFactors = FALSE)
gl    <- read.csv(file.path(canon, "hmsc_vs_gllvm_season.csv"), stringsAsFactors = FALSE)  # separation already dropped
sp <- function(ref, col, lab) {
  r <- setNames(ref[, c("OTU_ID", col)], c("OTU_ID", "ref"))
  k <- merge(beta_tbl[, c("OTU_ID", "beta_season_mean")], r, by = "OTU_ID")
  cat(sprintf("  vs %-30s n = %3d  Spearman rho = %.3f\n", lab, nrow(k),
              cor(k$beta_season_mean, k$ref, method = "spearman")))
}
cat("\nPer-OTU Season Beta rank agreement:\n")
sp(old,  "beta_season_mean", "canonical CLR (non-converged)")
sp(prob, "beta_season_mean", "canonical probit (converged)")
sp(gl,   "gllvm_season",     "GLLVM (separation dropped)")

cat("\nFig 5 OTUs:\n")
f5 <- c("OTU1", "OTU636", "OTU1225", "OTU1320", "OTU92", "OTU1369", "OTU320")
print(merge(data.frame(OTU_ID = f5), beta_tbl[, c("OTU_ID", "Genus", "beta_season_mean", "beta_season_support",
                                                   "beta_season_CrI2.5", "beta_season_CrI97.5")],
            by = "OTU_ID", all.x = TRUE), row.names = FALSE)
sink()
