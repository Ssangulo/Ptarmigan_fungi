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

card     <- do.call(rbind, lapply(scores, `[[`, "card"))
best     <- max(card$elpd_loo[card$tag != "baseline_current"])
card$elpd_diff_vs_best_arm <- card$elpd_loo - best
contrast <- do.call(rbind, lapply(scores, `[[`, "contrast"))
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
           "n_otu_dom_resolved90","n_otu_any_resolved95","n_otu_ptop_gt_0.8")
cat("\nH3 scorecard -- reported, not used to select:\n"); print(card[, show2], digits = 3, row.names = FALSE)
cat("\nContrast:\n"); print(contrast, digits = 3, row.names = FALSE)
cat("\nSeason vs Season+Year per-OTU structure, plant-slope agreement:\n"); print(sy, digits = 3, row.names = FALSE)
