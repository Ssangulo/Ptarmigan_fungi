# =============================================================================
# hmsc_converge/00_diagnostics.R  (branch exp/hmsc-converge -- experimental)
# One convergence test applied identically to every HMSC fit. Packages the
# per-chain checks that exposed the CLR fit's chain disagreement (2026-09-29).
#
# Usage: conda run -n r_env Rscript hmsc_converge/00_diagnostics.R <fit.rds> [meta.rds]
#   meta.rds (from 01_fit_arms.R) supplies count-based prevalence + guild; for
#   the canonical fits (no meta) prevalence is read off the CLR row-minimum and
#   guild from models/hmsc_otu_guild.csv.
# Writes models/hmsc_v2/<fit>_diag.rds (never next to a canonical fit) (a list of all tables below).
#
# PASS CRITERIA (fixed before any new fit was run):
#   Beta PSRF median <= 1.05, >= 95% < 1.1, max < 1.2
#   ESS >= 400 on the Fig 5 OTUs' Season Beta
#   per-chain guild-mean Season Beta and per-chain community VP within +/-0.05
#   (range across chains <= 0.10); Gamma median PSRF < 1.1 to quote Gamma
# =============================================================================

suppressMessages({ library(Hmsc); library(coda) })

args <- commandArgs(trailingOnly = TRUE)
fit_path  <- args[1]
meta_path <- if (length(args) >= 2) args[2] else NA
m <- readRDS(fit_path)
FIG5_OTUS <- c("OTU1", "OTU636", "OTU1225", "OTU1320", "OTU92", "OTU1369", "OTU320")

if (!is.na(meta_path)) {
  meta  <- readRDS(meta_path)
  prev  <- meta$prev_rows[colnames(m$Y)]
  guild <- meta$guild[colnames(m$Y)]
} else {
  iszero <- sweep(m$Y, 1, apply(m$Y, 1, min), "==")
  prev <- colSums(!iszero)
  gt <- read.csv("/home/daniel/Ptarmigan/models/hmsc_otu_guild.csv", stringsAsFactors = FALSE)
  guild <- setNames(gt$guild_class, gt$OTU_ID)[colnames(m$Y)]
}
stopifnot(!anyNA(guild))
nch <- length(m$postList)
xcn <- colnames(m$X)
s_row <- which(xcn == "Seasonsummer")
grp <- ifelse(grepl("^Year", xcn), 2L, 1L)

# ---- PSRF / ESS (psrf_summary logic from Scripts/10_hmsc.R:286-300) ---------
mp <- convertToCodaObject(m)
gb <- gelman.diag(mp$Beta, multivariate = FALSE)$psrf[, 1]
eb <- effectiveSize(mp$Beta)
gg <- gelman.diag(mp$Gamma, multivariate = FALSE)$psrf[, 1]
nm <- names(gb)
b_cov <- sapply(strsplit(sub("B[", "", nm, fixed = TRUE), " ", fixed = TRUE), function(z) z[1])
b_otu <- sapply(strsplit(nm, ", ", fixed = TRUE), function(z) strsplit(z[2], " ", fixed = TRUE)[[1]][1])

overall <- data.frame(
  beta_psrf_med = median(gb), beta_pct_lt1.1 = 100 * mean(gb < 1.1), beta_psrf_max = max(gb),
  beta_ess_min = min(eb), gamma_psrf_med = median(gg), gamma_psrf_max = max(gg))
by_cov <- do.call(rbind, lapply(split(seq_along(gb), b_cov), function(i)
  data.frame(covariate = b_cov[i[1]], n = length(i), psrf_med = median(gb[i]),
             psrf_p90 = unname(quantile(gb[i], 0.9)), psrf_max = max(gb[i]),
             pct_lt1.1 = 100 * mean(gb[i] < 1.1), ess_min = min(eb[i]))))
si <- which(b_cov == "Seasonsummer")
s_psrf <- setNames(gb[si], b_otu[si]); s_ess <- setNames(eb[si], b_otu[si])
pbin <- cut(prev[names(s_psrf)], c(0, 10, 20, 40, 80, 110))
by_prev <- do.call(rbind, lapply(split(names(s_psrf), pbin, drop = TRUE), function(o)
  data.frame(n = length(o), psrf_med = median(s_psrf[o]), psrf_max = max(s_psrf[o]))))
by_prev$prev_bin <- rownames(by_prev)

# ---- Per-chain Season Beta: Fig 5 OTUs + guild means --------------------------
Bs <- lapply(m$postList, function(ch) t(sapply(ch, function(s) s$Beta[s_row, ])))
f5 <- intersect(FIG5_OTUS, colnames(m$Y))
fig5_chain <- do.call(rbind, lapply(f5, function(o) {
  j <- match(o, colnames(m$Y))
  cm <- sapply(Bs, function(b) mean(b[, j]))
  data.frame(OTU = o, prev = prev[o], t(setNames(round(cm, 2), paste0("ch", seq_len(nch)))),
             range = round(diff(range(cm)), 2), psrf = round(s_psrf[o], 2),
             ess = round(s_ess[o]), check.names = FALSE)
}))
guild_chain <- sapply(Bs, function(b) tapply(colMeans(b), guild, mean))
colnames(guild_chain) <- paste0("ch", seq_len(nch))

# ---- Per-chain variance partition -------------------------------------------
vp_chain <- sapply(seq_len(nch), function(ch) {
  m1 <- m; m1$postList <- m$postList[ch]
  rowMeans(computeVariancePartitioning(m1, group = grp, groupnames = c("Season", "Year"))$vals)
})
colnames(vp_chain) <- paste0("ch", seq_len(nch))
vp_pooled <- rowMeans(computeVariancePartitioning(m, group = grp,
                                                  groupnames = c("Season", "Year"))$vals)

# ---- Season absorption by sample-level latent factors -----------------------
seas_unit <- tapply(as.numeric(m$XData$Season == "summer"), m$studyDesign$sample, mean)
rn <- levels(m$studyDesign$sample)
absorb <- sapply(m$postList, function(ch) {
  idx <- unique(round(seq(1, length(ch), length.out = 20)))
  median(sapply(ch[idx], function(s) {
    E <- s$Eta[[1]]; E <- E[, apply(E, 2, sd) > 0, drop = FALSE]
    if (ncol(E) == 0) return(0)
    max(abs(cor(E, seas_unit[rn])))
  }))
})
nf <- sapply(seq_along(m$ranLevels), function(r) ncol(m$postList[[1]][[1]]$Eta[[r]]))
names(nf) <- names(m$ranLevels)

# ---- Residual variance + posterior predictive check -------------------------
sig <- colMeans(do.call(rbind, lapply(m$postList, function(ch) t(sapply(ch, function(s) s$sigma)))))
pred <- computePredictedValues(m, expected = FALSE, nParallel = 1)   # ny x ns x draws
idx_d <- round(seq(1, dim(pred)[3], length.out = 200))
Yo <- m$Y
ppc <- do.call(rbind, lapply(seq_len(ncol(Yo)), function(j) {
  ok <- !is.na(Yo[, j]); y <- Yo[ok, j]
  rep_sd <- sapply(idx_d, function(d) sd(pred[ok, j, d]))
  rep_sk <- sapply(idx_d, function(d) { z <- pred[ok, j, d]; mean((z - mean(z))^3) / sd(z)^3 })
  obs_sk <- mean((y - mean(y))^3) / sd(y)^3
  data.frame(OTU = colnames(Yo)[j], p_sd = mean(rep_sd >= sd(y)), p_skew = mean(rep_sk >= obs_sk))
}))
ppc_summary <- c(pct_sd_extreme   = 100 * mean(ppc$p_sd   < 0.025 | ppc$p_sd   > 0.975),
                 pct_skew_extreme = 100 * mean(ppc$p_skew < 0.025 | ppc$p_skew > 0.975, na.rm = TRUE))

# ---- Pass / fail --------------------------------------------------------------
crit <- c(
  beta_psrf_med_le1.05   = overall$beta_psrf_med <= 1.05,
  beta_pct_lt1.1_ge95    = overall$beta_pct_lt1.1 >= 95,
  beta_psrf_max_lt1.2    = overall$beta_psrf_max < 1.2,
  fig5_ess_ge400         = all(fig5_chain$ess >= 400),
  guild_chain_range_le.1 = all(apply(guild_chain, 1, function(x) diff(range(x))) <= 0.10),
  vp_chain_range_le.1    = all(apply(vp_chain, 1, function(x) diff(range(x))) <= 0.10),
  gamma_psrf_med_lt1.1   = overall$gamma_psrf_med < 1.1)

# ---- Report -------------------------------------------------------------------
cat(sprintf("\n=== %s  (%d OTUs, %d chains x %d samples, thin %d, transient %d) ===\n",
            basename(fit_path), m$ns, nch, m$samples, m$thin, m$transient))
cat("latent factors per level:", paste(names(nf), nf, sep = "=", collapse = ", "), "\n")
cat("\nOverall:\n"); print(round(overall, 3), row.names = FALSE)
cat("\nBeta PSRF by covariate:\n"); print(by_cov, digits = 3, row.names = FALSE)
cat("\nSeason-Beta PSRF by prevalence (PCR rows present):\n"); print(by_prev, digits = 3, row.names = FALSE)
cat("\nPer-chain Season Beta, Fig 5 OTUs:\n"); print(fig5_chain, row.names = FALSE)
cat("\nPer-chain guild-mean Season Beta:\n"); print(round(guild_chain, 2))
cat("\nPer-chain community VP (pooled in last col):\n"); print(round(cbind(vp_chain, pooled = vp_pooled), 3))
cat("\nmedian max|cor(sample eta, Season)| per chain:", round(absorb, 2), "\n")
cat("residual variance (posterior mean) range:", signif(range(sig), 3), "\n")
cat("PPC: % OTUs with extreme SD p-value:", round(ppc_summary[1], 1),
    "; extreme skew p-value:", round(ppc_summary[2], 1), "\n")
cat("\nPASS CRITERIA:\n"); print(crit)
cat(if (all(crit[-7])) "==> PASSES Beta/VP criteria\n" else "==> FAILS\n")

saveRDS(list(overall = overall, by_cov = by_cov, by_prev = by_prev, fig5_chain = fig5_chain,
             guild_chain = guild_chain, vp_chain = vp_chain, vp_pooled = vp_pooled,
             absorb = absorb, nf = nf, sigma = sig, ppc = ppc, ppc_summary = ppc_summary,
             crit = crit),
        file.path("/home/daniel/Ptarmigan/models/hmsc_v2",
                  sub("\\.rds$", "_diag.rds", basename(fit_path))))
