# =============================================================================
# brms_refit/04_h3_compare.R  (branch exp/brms-refit -- experimental)
# Collect the H3 arm scores written by 03_h3_arms.R, score the CURRENT cached
# fit (models/H3_brms_joint_repRE.rds, read-only) with the same function as the
# baseline row, and write the comparison tables for the memo.
#
# Run: H3_MODE=pilot conda run -n r_env Rscript brms_refit/04_h3_compare.R
# Outputs ONLY to models/_brms_refit/.
# =============================================================================

MODE <- Sys.getenv("H3_MODE", "pilot")
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "02_h3_data.R"))
options(width = 220)
out_dir <- "/home/daniel/Ptarmigan/models/_brms_refit"

# ---- Baseline: the current appendix 7.2 fit, scored identically -------------
base_path <- file.path(out_dir, "h3_baseline_current_score.rds")
if (!file.exists(base_path)) {
  b0 <- readRDS("/home/daniel/Ptarmigan/models/H3_brms_joint_repRE.rds")
  saveRDS(score_h3(b0, c("Betula","Vaccinium","Empetrum","Picea","Calluna","Eriophorum"),
                   "baseline_current"), base_path)
  rm(b0); invisible(gc())
}
files  <- c(base_path, sort(list.files(out_dir, sprintf("^h3_%s_.*_score\\.rds$", MODE), full.names = TRUE)))
scores <- lapply(files, readRDS); names(scores) <- sapply(scores, function(s) s$card$tag)

rbind_fill <- function(l) {                    # arms differ in which columns they carry
  allcol <- unique(unlist(lapply(l, names)))
  do.call(rbind, lapply(l, function(x) { x[setdiff(allcol, names(x))] <- NA; x[allcol] }))
}
card     <- rbind_fill(lapply(scores, `[[`, "card"))
best     <- max(card$elpd_loo[card$tag != "baseline_current"])
card$elpd_diff_vs_best_arm <- card$elpd_loo - best
contrast <- rbind_fill(lapply(scores, `[[`, "contrast"))
dev <- do.call(rbind, lapply(names(scores), function(t) { x <- scores[[t]]
  rbind(if (!is.null(x$dev_mu)) data.frame(tag = t, part = "abundance",  plant = names(x$dev_mu$n_by_plant), n_otu_dev_resolved95 = unname(x$dev_mu$n_by_plant)),
        if (!is.null(x$dev_hu)) data.frame(tag = t, part = "occurrence", plant = names(x$dev_hu$n_by_plant), n_otu_dev_resolved95 = unname(x$dev_hu$n_by_plant))) }))
write.csv(dev, file.path(out_dir, sprintf("h3_%s_dev_resolved.csv", MODE)), row.names = FALSE)
per_otu  <- do.call(rbind, lapply(scores, `[[`, "per_otu"))
write.csv(card,     file.path(out_dir, sprintf("h3_%s_scorecard.csv", MODE)), row.names = FALSE)
write.csv(contrast, file.path(out_dir, sprintf("h3_%s_contrast.csv",  MODE)), row.names = FALSE)
write.csv(per_otu,  file.path(out_dir, sprintf("h3_%s_per_otu.csv",   MODE)), row.names = FALSE)

# ---- Pairwise LOO differences with SE (same response in every arm) ----------
arms_only <- scores[names(scores) != "baseline_current"]
if (length(arms_only) > 1) {
  lc <- loo::loo_compare(lapply(arms_only, `[[`, "loo"))
  write.csv(data.frame(tag = rownames(lc), lc[, c("elpd_diff", "se_diff")]),
            file.path(out_dir, sprintf("h3_%s_loo_compare.csv", MODE)), row.names = FALSE)
  cat("\nLOO comparison (arms only):\n"); print(lc[, c("elpd_diff", "se_diff")], digits = 4)
}

# ---- Do plant slopes move when Year is allowed to vary by OTU? --------------
tg <- names(arms_only)
same <- function(t) sub("_(species|genus)$", "", t)
pg <- do.call(rbind, lapply(grep("_species$", tg, value = TRUE), function(t_sp) {
  t_g <- sub("_species$", "_genus", t_sp); if (!t_g %in% tg) return(NULL)
  a <- scores[[t_sp]]$slope_med; b <- scores[[t_g]]$slope_med[rownames(a), ]
  data.frame(arm = same(t_sp),
             cor_Betula = cor(a[, "Betula_sp"], b[, "Betula"]),
             cor_Empetrum = cor(a[, "Empetrum_nigrum"], b[, "Empetrum"]),
             cor_Vacc_genus_vs_myrtillus = cor(a[, "Vaccinium_myrtillus"], b[, "Vaccinium"]),
             cor_Vacc_genus_vs_uliginosum = cor(a[, "Vaccinium_uliginosum"], b[, "Vaccinium"]),
             stringsAsFactors = FALSE)
}))
if (!is.null(pg)) write.csv(pg, file.path(out_dir, sprintf("h3_%s_species_vs_genus_slopes.csv", MODE)), row.names = FALSE)
sy <- do.call(rbind, lapply(grep("_S_", tg, value = TRUE), function(t_s) {
  t_sy <- sub("_S_", "_SY_", t_s); if (!t_sy %in% tg) return(NULL)
  a <- scores[[t_s]]$slope_med; b <- scores[[t_sy]]$slope_med[rownames(a), colnames(a)]
  data.frame(arm_S = t_s, arm_SY = t_sy, cor_all = cor(as.vector(a), as.vector(b)),
             mean_abs_S = mean(abs(a)), mean_abs_SY = mean(abs(b)),
             max_abs_change = max(abs(a - b)),
             t(setNames(sapply(colnames(a), function(p) cor(a[, p], b[, p])), paste0("cor_", colnames(a)))),
             stringsAsFactors = FALSE)
}))
if (!is.null(sy)) write.csv(sy, file.path(out_dir, sprintf("h3_%s_S_vs_SY_slopes.csv", MODE)), row.names = FALSE)

show <- c("tag","divergences","max_treedepth","max_rhat","min_ess_bulk","hours_per_chain",
          "pass_convergence","obs_mean","ppc_mean_lwr","ppc_mean_upr","ppc_P_max_le_maxlib",
          "otu_zero_cor","otu_zero_maxgap","pass_ppc","elpd_loo","elpd_se","n_pareto_k_bad")
cat("\nH3 scorecard -- adequacy:\n"); print(card[, show], digits = 3, row.names = FALSE)
show2 <- c("tag","sd_OTU_Season","sd_SampleOTU","shape","sd_plant_min","sd_plant_max",
           "n_otu_dom_resolved90","n_otu_any_resolved95","n_otu_ptop_gt_0.8",
           "sd_hu_plant_min","sd_hu_plant_max","n_otu_hu_any_resolved95",
           "n_dev_resolved95_mu","n_dev_resolved95_hu")
cat("\nH3 scorecard -- reported, not used to select:\n"); print(card[, show2], digits = 3, row.names = FALSE)
cat("\nContrast:\n"); print(contrast, digits = 3, row.names = FALSE)
cat("\nSpecies vs genus predictor set, per-OTU slope agreement:\n"); print(pg, digits = 3, row.names = FALSE)
cat("\nSeason vs Season+Year per-OTU structure, plant-slope agreement:\n"); print(sy, digits = 3, row.names = FALSE)
