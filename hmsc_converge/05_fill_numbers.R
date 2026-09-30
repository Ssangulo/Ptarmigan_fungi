# =============================================================================
# hmsc_converge/05_fill_numbers.R  (branch exp/hmsc-converge -- one-off)
# Replace the __TOKEN__ placeholders left in Supplementary_Appendix.qmd while
# the revised HMSC was still running. Every value is computed here from the
# STAGED tables (Supplementary/tables), with the same selection rule as the
# fig3-build chunk, so prose, drift guards and figure share one source.
# Tokens inside the fig3-build chunk get plain numbers (guards); tokens in
# prose get formatted text. Refuses to write if any token is left unfilled.
# Run from the worktree root: conda run -n r_env Rscript hmsc_converge/05_fill_numbers.R
# =============================================================================

qmd <- "Supplementary/Supplementary_Appendix.qmd"
td  <- "Supplementary/tables"
b   <- read.csv(file.path(td, "hmsc_beta_season.csv"), stringsAsFactors = FALSE)
vp  <- read.csv(file.path(td, "hmsc_variance_partition.csv"), check.names = FALSE)
pcr <- read.csv(file.path(td, "hmsc_pcr_replicate_variance.csv"), stringsAsFactors = FALSE)
r2  <- read.csv(file.path(td, "hmsc_predictive_R2_summary.csv"), stringsAsFactors = FALSE)
sh  <- read.csv(file.path(td, "dark_taxa_SH_matching.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(b) == 226, nrow(vp) == 226)

gm  <- tapply(b$beta_season_mean, b$guild, mean)
vpa <- 100 * colMeans(vp[, c("Season", "Year", "Random: sample")])
vpg <- 100 * tapply(vp$Season, vp$guild, mean)
n_ci <- sum(b$beta_season_CrI2.5 > 0 | b$beta_season_CrI97.5 < 0)

# Figure 3A selection, exactly as in fig3-build
si <- match(b$OTU_ID, sh$OTU_ID)
b$share <- (sh$winter_mean_share[si] + sh$summer_mean_share[si]) / 2
b$rel <- b$share / sum(b$share)
sel <- b[b$rel >= 0.01, ]
stopifnot(nrow(sel) == 22)
n_res <- sum(sel$beta_season_support >= 0.95 | sel$beta_season_support <= 0.05)
XCLIP <- 6
clip <- sel[sel$beta_season_CrI2.5 < -XCLIP | sel$beta_season_CrI97.5 > XCLIP, ]
fmt_signed <- function(x, d = 2) sub("^-", "−", sprintf(paste0("%+.", d, "f"), x))
clip_txt <- if (nrow(clip) == 0) "No interval extends beyond the ±6 axis." else
  if (nrow(clip) == 1) sprintf("One interval, %s's, extends beyond the axis and is clipped at %s6 with an arrowhead (it reaches %s).",
                               clip$OTU_ID, if (clip$beta_season_CrI2.5 < -XCLIP) "−" else "+",
                               fmt_signed(if (clip$beta_season_CrI2.5 < -XCLIP) clip$beta_season_CrI2.5 else clip$beta_season_CrI97.5)) else
  sprintf("%d intervals extend beyond the axis and are clipped at ±6 with an arrowhead.", nrow(clip))
b$rank <- rank(-b$beta_season_mean)
dung_rank <- round(mean(b$rank[b$guild == "dung_saprotroph"]))
rc <- r2[r2$model == "clr", ]

guard <- c(`__GM_DUNG__` = sprintf("%.3f", gm[["dung_saprotroph"]]),
           `__GM_PLANT__` = sprintf("%.3f", gm[["plant_associated"]]),
           `__GM_PATHO__` = sprintf("%.3f", gm[["pathotroph"]]),
           `__GM_OTHER__` = sprintf("%.3f", gm[["other"]]),
           `__GM_DARK__`  = sprintf("%.3f", gm[["dark_unassigned"]]),
           `__N_CI__`     = as.character(n_ci),
           `__VP_SEASON__` = sprintf("%.2f", vpa[["Season"]]),
           `__VP_YEAR__`   = sprintf("%.2f", vpa[["Year"]]),
           `__VP_SAMPLE__` = sprintf("%.2f", vpa[["Random: sample"]]),
           `__PCR_ALL__`   = sprintf("%.3f", pcr$share_total_mean[pcr$group == "all"]))
prose <- c(`__GM_DUNG__` = fmt_signed(gm[["dung_saprotroph"]]),
           `__GM_DARK__` = fmt_signed(gm[["dark_unassigned"]]),
           `__GM_PATHO__` = fmt_signed(gm[["pathotroph"]]),
           `__GM_OTHER__` = fmt_signed(gm[["other"]]),
           `__N_RES__` = as.character(n_res),
           `__N_CLIP_TEXT__` = clip_txt,
           `__VPG_DUNG__` = sprintf("%.1f", vpg[["dung_saprotroph"]]),
           `__R2_EXPL__` = sprintf("%.2f", rc$mean_expl_R2),
           `__R2_CV__` = sprintf("%.2f", rc$mean_cv_R2),
           `__DUNG_RANK__` = as.character(dung_rank))

x <- readLines(qmd, encoding = "UTF-8")
s <- grep("^```\\{r fig3-build\\}", x); e <- s + which(x[(s + 1):length(x)] == "```")[1]
in_chunk <- seq_along(x) >= s & seq_along(x) <= e
for (k in names(guard)) x[in_chunk] <- gsub(k, guard[[k]], x[in_chunk], fixed = TRUE)
for (k in names(prose)) x[!in_chunk] <- gsub(k, prose[[k]], x[!in_chunk], fixed = TRUE)
left <- grep("__[A-Z0-9_]+__", x, value = TRUE)
left <- left[!grepl("^\\s*#", left)]
if (length(left)) stop("Unfilled tokens remain:\n", paste(left, collapse = "\n"))
writeLines(x, qmd, useBytes = TRUE)

cat("Filled. Guard values:\n"); print(guard)
cat("Prose values:\n"); print(prose)
cat("Guild-mean Season share of variance (%):\n"); print(round(vpg, 1))
