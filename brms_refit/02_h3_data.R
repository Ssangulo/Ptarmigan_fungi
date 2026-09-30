# =============================================================================
# brms_refit/02_h3_data.R  (branch exp/brms-refit -- experimental)
# Shared data build + scoring helpers for the H3 pilot arms. source()d by
# 03_h3_arms.R and 04_h3_compare.R; does not fit or write anything itself.
#
# The matched COUNT object and the long model frame are copied VERBATIM from
# Scripts/6_diversity_analyses.R Section 3c (lines 936-961 and 1004-1017). The
# ONLY change is the plant predictors: script 6 uses decostand("rclr") on six
# genera, which imputes zeros by rank-3 matrix completion (vegan ropt = 3) and
# returns an exactly rank-3 predictor matrix. build_plant_predictors() below
# replaces zeros with a fixed diet share and takes the CLR over the FULL
# composition, which is full rank and independent of plant read depth.
# =============================================================================

suppressMessages({ library(phyloseq); library(brms); library(posterior) })

MIN_OTU_PREV_H3 <- 5      # verbatim, script 6
COPRO_GENERA <- c("Sordaria","Podospora","Cercophora","Chaetomium","Schizothecium",
                  "Preussia","Delitschia","Pilobolus","Ascobolus","Saccobolus",
                  "Sporormiella","Coprinopsis","Thelebolus","Coniochaeta")
PLANT_SETS <- list(
  species = c("Betula_sp", "Vaccinium_myrtillus", "Vaccinium_uliginosum", "Empetrum_nigrum"),
  genus   = c("Betula", "Vaccinium", "Empetrum")
)

otu_mat_of <- function(ps) {
  m <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) m <- t(m)
  m
}

# ---- Matched COUNT object (verbatim, script 6 lines 936-961) ----------------
load("/home/daniel/Ptarmigan/trimmed/mergedPlates/eco_analysis.RData")
stopifnot(exists("funguild_otu"))
plant <- readRDS("/home/daniel/Ptarmigan/plant_ITS/phyloseq_plant_ITS.rds")

ps_cnt  <- alldat_full$nopool                # PCR-rep-level, min_depth_full=1000
md_cnt  <- data.frame(as(sample_data(ps_cnt), "data.frame"), stringsAsFactors=FALSE)
match_r <- rownames(md_cnt)[md_cnt$Sample_ID_field %in% sample_names(plant)]

psf <- prune_samples(match_r, ps_cnt)
psf <- prune_taxa(taxa_sums(psf) > 0, psf)
md_f <- data.frame(as(sample_data(psf), "data.frame"), stringsAsFactors=FALSE)
md_f$Season <- factor(md_f$Season, levels=c("winter","summer"))
md_f$Year   <- factor(md_f$Year)

fom     <- otu_mat_of(psf)                   # PCR-rep rows x fungal OTUs (counts)
libsize <- rowSums(fom)                      # per-PCR-rep library size (all OTUs)
n_bio   <- length(unique(md_f$Sample_ID_field))
pres_bysamp <- rowsum((fom > 0) * 1, group = md_f$Sample_ID_field) > 0   # bio-sample x OTU
prev_f   <- colSums(pres_bysamp)
keep_otu <- names(prev_f)[prev_f >= MIN_OTU_PREV_H3]
fom_k    <- fom[, keep_otu, drop=FALSE]
stopifnot(n_bio == 30, nrow(fom_k) == 60, ncol(fom_k) == 46)

# ---- Plant predictors: fixed-share zero replacement -> CLR over full comp ---
field_ids <- unique(md_f$Sample_ID_field)
plant_m <- prune_samples(field_ids, plant)
plant_m <- prune_taxa(taxa_sums(plant_m) > 0, plant_m)
pom     <- otu_mat_of(plant_m)[field_ids, , drop = FALSE]
ptt     <- data.frame(as(tax_table(plant_m), "matrix"), stringsAsFactors=FALSE)[colnames(pom), ]

.na_lab <- function(x, bad) ifelse(is.na(x) | x %in% bad, NA, x)
g_lab <- .na_lab(ptt$Genus,   c("", "NA", "g__"))
f_lab <- .na_lab(ptt$Family,  c("", "NA", "f__"))
s_lab <- .na_lab(ptt$Species, c("", "NA", "s__"))
lab_genus   <- ifelse(!is.na(g_lab), g_lab,
                      ifelse(!is.na(f_lab), paste0(f_lab, "_fam"), colnames(pom)))
lab_species <- ifelse(!is.na(s_lab), s_lab, lab_genus)

build_plant_predictors <- function(level = c("species", "genus"), delta = 1e-3) {
  level <- match.arg(level)
  lab <- if (level == "species") lab_species else lab_genus
  agg <- t(rowsum(t(pom), group = lab))              # bio samples x plant taxon
  P   <- agg / rowSums(agg)
  Z   <- P                                           # multiplicative zero replacement
  for (i in seq_len(nrow(Z))) {
    z <- Z[i, ] == 0
    Z[i, z]  <- delta
    Z[i, !z] <- Z[i, !z] * (1 - sum(z) * delta)
  }
  clr <- log(Z); clr <- clr - rowMeans(clr)          # CLR over the FULL composition
  keep <- PLANT_SETS[[level]]
  stopifnot(all(keep %in% colnames(clr)))
  X <- scale(clr[, keep, drop = FALSE])
  sv <- svd(X)$d
  stopifnot(!anyNA(X), qr(X)$rank == ncol(X), max(sv) / min(sv) < 10)
  attr(X, "cond") <- max(sv) / min(sv)
  X
}

# ---- Long model frame (verbatim, script 6 lines 1004-1017) ------------------
build_long <- function(pg_z) {
  plant_var <- make.names(colnames(pg_z)); colnames(pg_z) <- plant_var
  field_of <- setNames(md_f$Sample_ID_field, rownames(md_f))
  pg_byrow <- pg_z[field_of[rownames(fom_k)], , drop=FALSE]
  rownames(pg_byrow) <- rownames(fom_k)
  stopifnot(!anyNA(pg_byrow))
  otus <- colnames(fom_k); ns <- nrow(fom_k)
  long <- data.frame(
    row_id          = rep(rownames(fom_k), times=length(otus)),
    OTU             = rep(otus,            each=ns),
    count           = as.vector(fom_k),       # column-major: matches rep() above
    stringsAsFactors = FALSE
  )
  long$OTU             <- factor(long$OTU)
  long$log_libsize     <- log(libsize[long$row_id])
  long$Sample_ID_field <- field_of[long$row_id]
  long$Season          <- md_f$Season[match(long$row_id, rownames(md_f))]
  long$Year            <- md_f$Year[match(long$row_id, rownames(md_f))]
  for (v in plant_var) long[[v]] <- pg_byrow[long$row_id, v]
  list(long = long, plant_var = plant_var)
}

# ---- Guild grouping (verbatim, script 6 lines 1100-1117) --------------------
guild_groups <- function(otu_ids) {
  guild_str <- funguild_otu[match(otu_ids, funguild_otu$OTU_ID), ]$Guild
  is_copro <- !is.na(guild_str) & grepl("Dung Saprotroph", guild_str)
  is_plant <- !is.na(guild_str) & !is_copro &
              grepl("Endophyte|Plant Saprotroph|Plant Pathogen|Epiphyte", guild_str)
  guild_group <- ifelse(is_copro, "coprophilous",
                  ifelse(is_plant, "plant_associated", "other/unassigned"))
  tax_all <- data.frame(as(tax_table(alldat$nopool), "matrix"), stringsAsFactors=FALSE)
  otu_gen <- sub("^[a-z]__", "", tax_all$Genus[match(otu_ids, rownames(tax_all))])
  copro_by_genus <- !is.na(otu_gen) & otu_gen %in% COPRO_GENERA
  group_genus <- ifelse(copro_by_genus, "coprophilous",
                  ifelse(guild_group == "plant_associated", "plant_associated", "other/unassigned"))
  data.frame(OTU_ID = otu_ids, Genus = otu_gen, guild_group = guild_group,
             group_genus = group_genus, stringsAsFactors = FALSE)
}

# ---- One scoring function for every arm (and the cached baseline) -----------
hhi_of <- function(v) { a <- abs(v); s <- sum(a); if (s == 0) NA_real_ else sum((a/s)^2) }

score_h3 <- function(fit, plant_var, tag, ppc_draws = 500) {
  d   <- fit$data
  np  <- nuts_params(fit)
  dr  <- as_draws_df(fit)
  sm  <- summarise_draws(dr, "median", "rhat", "ess_bulk", "ess_tail")
  smx <- sm[!grepl("^(lp__|lprior)$", sm$variable), ]
  med <- setNames(sm$median, sm$variable)
  el  <- tryCatch(rstan::get_elapsed_time(fit$fit), error = function(e) NULL)

  # --- predictive checks ---
  yrep  <- posterior_predict(fit, ndraws = ppc_draws)
  pm    <- rowMeans(yrep); pmax <- apply(yrep, 1, max)
  oz    <- tapply(d$count == 0, d$OTU, mean)
  pz    <- sapply(names(oz), function(o) mean(yrep[, d$OTU == o] == 0))
  maxlib <- max(exp(d$log_libsize))
  lo    <- suppressWarnings(loo(fit))

  # --- per-OTU plant slopes ---
  co     <- coef(fit, summary = FALSE)$OTU
  otu_ids <- dimnames(co)[[2]]
  slopes <- co[, , match(plant_var, dimnames(co)[[3]]), drop = FALSE]   # draws x OTU x plant
  stopifnot(!anyNA(dim(slopes)), dim(slopes)[3] == length(plant_var))
  k      <- length(plant_var)
  s_med  <- apply(slopes, c(2, 3), median)
  s_l95  <- apply(slopes, c(2, 3), quantile, 0.025); s_u95 <- apply(slopes, c(2, 3), quantile, 0.975)
  s_l90  <- apply(slopes, c(2, 3), quantile, 0.05);  s_u90 <- apply(slopes, c(2, 3), quantile, 0.95)
  hhi_pt <- apply(s_med, 1, hhi_of)
  absS   <- abs(slopes)
  hhi_dr <- sapply(seq_along(otu_ids), function(o) { a <- absS[, o, ]; rowSums((a / rowSums(a))^2) })
  # noise-robust index: posterior probability that each plant carries the largest |slope|
  p_top  <- t(sapply(seq_along(otu_ids), function(o)
              tabulate(max.col(absS[, o, ], ties.method = "first"), nbins = k) / dim(absS)[1]))
  dimnames(p_top) <- list(otu_ids, plant_var)
  dom    <- apply(abs(s_med), 1, which.max)
  idx    <- cbind(seq_along(otu_ids), dom)
  res90_dom <- (s_l90[idx] > 0) | (s_u90[idx] < 0)
  res95_any <- rowSums((s_l95 > 0) | (s_u95 < 0)) > 0

  gg <- guild_groups(otu_ids)
  contrast_of <- function(grp, label) {
    ci <- which(grp == "coprophilous"); pi <- which(grp == "plant_associated")
    dd <- rowMeans(hhi_dr[, ci, drop = FALSE]) - rowMeans(hhi_dr[, pi, drop = FALSE])
    data.frame(tag = tag, grouping = label, n_copro = length(ci), n_plant = length(pi),
      copro_hhi_med = median(hhi_pt[ci]), plant_hhi_med = median(hhi_pt[pi]),
      wilcox_p_copro_lower = suppressWarnings(wilcox.test(hhi_pt[ci], hhi_pt[pi], alternative = "less")$p.value),
      P_copro_more_diffuse = mean(dd < 0),
      copro_ptop_med = median(apply(p_top[ci, , drop = FALSE], 1, max)),
      plant_ptop_med = median(apply(p_top[pi, , drop = FALSE], 1, max)),
      stringsAsFactors = FALSE)
  }

  sdv <- function(nm) if (nm %in% names(med)) unname(med[nm]) else NA_real_
  card <- data.frame(
    tag = tag, family = fit$family$family, n_par = nrow(smx),
    divergences = sum(np$Value[np$Parameter == "divergent__"]),
    max_treedepth = max(np$Value[np$Parameter == "treedepth__"]),
    max_rhat = max(smx$rhat, na.rm = TRUE), frac_rhat_gt_1.01 = mean(smx$rhat > 1.01, na.rm = TRUE),
    min_ess_bulk = min(smx$ess_bulk, na.rm = TRUE), min_ess_tail = min(smx$ess_tail, na.rm = TRUE),
    hours_per_chain = if (is.null(el)) NA_real_ else mean(rowSums(el)) / 3600,
    obs_mean = mean(d$count), ppc_mean_lwr = unname(quantile(pm, .025)),
    ppc_mean_med = median(pm), ppc_mean_upr = unname(quantile(pm, .975)),
    ppc_P_max_le_maxlib = mean(pmax <= maxlib), ppc_max_med = median(pmax), obs_max = max(d$count),
    obs_zero = mean(d$count == 0), ppc_zero = mean(yrep == 0),
    otu_zero_cor = cor(oz, pz), otu_zero_maxgap = max(abs(oz - pz)),
    elpd_loo = lo$estimates["elpd_loo", "Estimate"], elpd_se = lo$estimates["elpd_loo", "SE"],
    p_loo = lo$estimates["p_loo", "Estimate"], n_pareto_k_bad = sum(lo$diagnostics$pareto_k > 0.7),
    sd_OTU_Intercept = sdv("sd_OTU__Intercept"), sd_OTU_Season = sdv("sd_OTU__Seasonsummer"),
    sd_SampleOTU = sdv("sd_Sample_ID_field:OTU__Intercept"), shape = sdv("shape"),
    sd_plant_min = min(sapply(plant_var, function(v) sdv(paste0("sd_OTU__", v)))),
    sd_plant_max = max(sapply(plant_var, function(v) sdv(paste0("sd_OTU__", v)))),
    n_otu_dom_resolved90 = sum(res90_dom), n_otu_any_resolved95 = sum(res95_any),
    n_otu_ptop_gt_0.8 = sum(apply(p_top, 1, max) > 0.8),
    stringsAsFactors = FALSE)
  card$pass_convergence <- with(card, divergences == 0 & max_rhat < 1.01 & min_ess_bulk > 400)
  card$pass_ppc <- with(card, obs_mean >= ppc_mean_lwr & obs_mean <= ppc_mean_upr &
                          ppc_P_max_le_maxlib >= 0.95 & otu_zero_cor > 0.9 & otu_zero_maxgap < 0.10)

  per_otu <- data.frame(tag = tag, gg, hhi = hhi_pt, p_top_max = apply(p_top, 1, max),
                        top_plant_by_ptop = plant_var[max.col(p_top, ties.method = "first")],
                        dominant_plant = plant_var[dom], dom_med = s_med[idx],
                        dom_l90 = s_l90[idx], dom_u90 = s_u90[idx],
                        resolved90_dom = res90_dom, resolved95_any = res95_any,
                        stringsAsFactors = FALSE)
  sd_tab <- smx[grepl("^(sd_|shape$|b_)", smx$variable), ]
  list(card = card, per_otu = per_otu, slope_med = s_med, slope_l95 = s_l95, slope_u95 = s_u95,
       p_top = p_top, contrast = rbind(contrast_of(gg$guild_group, "FUNGuild"),
                                       contrast_of(gg$group_genus, "COPRO_GENERA")),
       hyper = as.data.frame(sd_tab), loo = lo)
}
