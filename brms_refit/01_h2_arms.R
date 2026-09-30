# =============================================================================
# brms_refit/01_h2_arms.R  (branch exp/brms-refit -- experimental)
# H2 mechanism (appendix 7.1): alternatives to the current brms model.
# The model frame is script 6's own output (models/H2_diet_richness_matched_
# data.csv); priors and sampler settings are copied VERBATIM from
# Scripts/6_diversity_analyses.R Section 3 / 3b. Only the family / Year
# handling differs by arm:
#
#   G    : gaussian on z-scored q1            (CURRENT primary / drop-Year)
#   SN   : skew_normal on z-scored q1
#   T    : student on z-scored q1             (robust to the two extreme q1)
#   LN   : lognormal on raw q1
#   GA   : Gamma(log) on raw q1               (family of Table S9)
#   YRE  : gaussian, (1 | Year), half-normal(0, 0.5) SD prior (prereg-style)
#   Gm1  : gaussian, S_1_7_P1_8E removed      (influence row, not a candidate)
# each under  ~ plant_richness_z + Season + Year  (SY)  and  + Season  (S).
#
# Run: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
#        conda run -n r_env Rscript brms_refit/01_h2_arms.R
# Outputs ONLY to models/_brms_refit/ -- never touches the canonical H2 caches.
# =============================================================================

suppressMessages({ library(brms); library(posterior); library(phyloseq) })
options(width = 200)

out_dir <- "/home/daniel/Ptarmigan/models/_brms_refit"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

diet <- read.csv("/home/daniel/Ptarmigan/models/H2_diet_richness_matched_data.csv",
                 stringsAsFactors = FALSE)
stopifnot(nrow(diet) == 27)
diet$Season <- factor(diet$Season, levels = c("winter", "summer"))
diet$Year   <- factor(diet$Year)

# ---- 0. Covariate screen ----------------------------------------------------
# Residuals of the current primary structure against every other sample_data
# column plus fungal / plant read depth. Reported only.
load("/home/daniel/Ptarmigan/trimmed/mergedPlates/eco_analysis.RData")
plant <- readRDS("/home/daniel/Ptarmigan/plant_ITS/phyloseq_plant_ITS.rds")
md    <- as(sample_data(alldat$nopool), "data.frame")[diet$sample, ]
pmd   <- as(sample_data(plant), "data.frame")[md$Sample_ID_field, ]
res   <- resid(lm(fungal_hill_z ~ plant_richness_z + Season + Year, data = diet))
cov_df <- cbind(md, plant = pmd[, setdiff(names(pmd), names(md)), drop = FALSE],
                fungal_depth = sample_sums(alldat$nopool)[diet$sample],
                plant_depth  = sample_sums(plant)[md$Sample_ID_field])
screen <- do.call(rbind, lapply(names(cov_df), function(v) {
  x <- cov_df[[v]]; ok <- !is.na(x); nu <- length(unique(x[ok]))
  if (nu < 2) return(NULL)
  xn <- suppressWarnings(as.numeric(as.character(x)))
  if (is.numeric(x) || (sum(!is.na(xn)) == sum(ok) && nu > 8)) {
    ct <- suppressWarnings(cor.test(xn[ok], res[ok], method = "spearman"))
    data.frame(covariate = v, type = "numeric", n = sum(ok), levels = nu,
               stat = unname(ct$estimate), p = ct$p.value)
  } else if (nu <= 8) {
    tb <- table(x[ok]); if (sum(tb >= 2) < 2) return(NULL)
    kt <- kruskal.test(res[ok] ~ factor(x[ok]))
    data.frame(covariate = v, type = "factor", n = sum(ok), levels = nu,
               stat = unname(kt$statistic), p = kt$p.value)
  } else NULL
}))
screen <- screen[order(screen$p), ]
write.csv(screen, file.path(out_dir, "h2_covariate_screen.csv"), row.names = FALSE)
cat("\nH2 covariate screen (residuals of ~ plant_richness_z + Season + Year):\n")
print(screen, digits = 3, row.names = FALSE)

# ---- 1. Arms ----------------------------------------------------------------
prior_z   <- c(set_prior("normal(0,1)",   class = "b"),          # verbatim, script 6
               set_prior("normal(0,0.5)", class = "Intercept"))
prior_raw <- c(set_prior("normal(0,1)",   class = "b"),          # log scale
               set_prior("normal(2.3,1)", class = "Intercept"))  # log(10) ~ median q1

arms <- list(
  G   = list(resp = "fungal_hill_z", family = gaussian(),          prior = prior_z),
  SN  = list(resp = "fungal_hill_z", family = skew_normal(),       prior = prior_z),
  T   = list(resp = "fungal_hill_z", family = student(),           prior = prior_z),
  LN  = list(resp = "fungal_hill",   family = lognormal(),         prior = prior_raw),
  GA  = list(resp = "fungal_hill",   family = Gamma(link = "log"), prior = prior_raw),
  YRE = list(resp = "fungal_hill_z", family = gaussian(),
             prior = c(prior_z, set_prior("normal(0,0.5)", class = "sd")), year_re = TRUE),
  Gm1 = list(resp = "fungal_hill_z", family = gaussian(),          prior = prior_z,
             drop = "S_1_7_P1_8E")
)
structs <- c(SY = "plant_richness_z + Season + Year", S = "plant_richness_z + Season")

sk <- function(x) { m <- mean(x); mean((x - m)^3) / sd(x)^3 }
q1_sd <- sd(diet$fungal_hill)

score <- function(fit, arm, st, dat, raw_scale, comparable) {
  d  <- as_draws_df(fit)$b_plant_richness_z
  np <- nuts_params(fit)
  s  <- summarise_draws(as_draws_df(fit), "rhat", "ess_bulk", "ess_tail")
  s  <- s[!grepl("^(lp__|lprior)", s$variable), ]
  y    <- dat[[if (raw_scale) "fungal_hill" else "fungal_hill_z"]]
  yrep <- posterior_predict(fit, ndraws = 4000)
  ll   <- log_lik(fit); if (!raw_scale) ll <- ll - log(q1_sd)   # Jacobian: z -> raw q1
  lo   <- suppressWarnings(loo::loo(ll, r_eff = loo::relative_eff(exp(ll),
                                   chain_id = rep(1:4, each = nrow(ll) / 4))))
  data.frame(
    arm = arm, structure = st, family = fit$family$family, n = nrow(dat),
    slope_med = median(d), lwr95 = unname(quantile(d, .025)), upr95 = unname(quantile(d, .975)),
    P_gt0 = mean(d > 0), slope_scale = if (raw_scale) "log q1 per SD" else "SD q1 per SD",
    divergences = sum(np$Value[np$Parameter == "divergent__"]),
    max_rhat = max(s$rhat, na.rm = TRUE),
    min_ess_bulk = min(s$ess_bulk, na.rm = TRUE), min_ess_tail = min(s$ess_tail, na.rm = TRUE),
    ppc_p_skew = mean(apply(yrep, 1, sk) >= sk(y)),
    ppc_p_max  = mean(apply(yrep, 1, max) >= max(y)),
    ppc_p_min  = mean(apply(yrep, 1, min) <= min(y)),
    elpd_loo_raw = if (comparable) lo$estimates["elpd_loo", "Estimate"] else NA_real_,
    elpd_se      = if (comparable) lo$estimates["elpd_loo", "SE"] else NA_real_,
    n_pareto_k_bad = sum(lo$diagnostics$pareto_k > 0.7),
    stringsAsFactors = FALSE)
}

rows <- list()
for (a in names(arms)) for (st in names(structs)) {
  A <- arms[[a]]
  if (isTRUE(A$year_re) && st == "S") next          # YRE has no drop-Year counterpart
  rhs <- if (isTRUE(A$year_re)) "plant_richness_z + Season + (1 | Year)" else structs[[st]]
  dat <- if (is.null(A$drop)) diet else diet[diet$sample != A$drop, ]
  cat(sprintf("\n=== arm %s / %s : %s ~ %s ===\n", a, st, A$resp, rhs))
  fit <- brm(as.formula(paste(A$resp, "~", rhs)), data = dat, family = A$family,
             prior = A$prior, chains = 4, iter = 6000, warmup = 3000, cores = 4,
             control = list(adapt_delta = 0.999, max_treedepth = 15),
             seed = 1, refresh = 0,
             file = file.path(out_dir, sprintf("h2_%s_%s", a, st)))
  rows[[paste(a, st)]] <- score(fit, a, st, dat, raw_scale = A$resp == "fungal_hill",
                                comparable = is.null(A$drop))
  print(rows[[paste(a, st)]], digits = 3, row.names = FALSE)
}
sc <- do.call(rbind, rows)
write.csv(sc, file.path(out_dir, "h2_scorecard.csv"), row.names = FALSE)

# ---- 2. Reproduction guard --------------------------------------------------
g <- sc[sc$arm == "G", ]
cat(sprintf("\nReproduction guard: SY median %.4f P %.4f (target 0.2891 / 0.9144); S median %.4f P %.4f (target 0.3837 / 0.9634)\n",
            g$slope_med[g$structure == "SY"], g$P_gt0[g$structure == "SY"],
            g$slope_med[g$structure == "S"],  g$P_gt0[g$structure == "S"]))
stopifnot(abs(g$slope_med[g$structure == "SY"] - 0.2891) < 0.005,
          abs(g$P_gt0[g$structure == "SY"]     - 0.9144) < 0.005,
          abs(g$slope_med[g$structure == "S"]  - 0.3837) < 0.005,
          abs(g$P_gt0[g$structure == "S"]      - 0.9634) < 0.005)
cat("\nH2 scorecard:\n"); print(sc, digits = 3, row.names = FALSE)
