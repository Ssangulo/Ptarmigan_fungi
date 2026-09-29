# =============================================================================
# hmsc_converge/04_v1_disagreement_table.R  (branch exp/hmsc-converge)
# One-off provenance table: how the four chains of the SUPERSEDED CLR HMSC
# (fitted 2026-07-20, 7 sample + 5 PCR latent factors, uncapped) disagreed on
# the quantities Figure 5 reported. Read from the saved diagnostics of that fit
# (00_diagnostics.R on models/hmsc_clr_fit.rds -> models/hmsc_v2/hmsc_clr_fit_diag.rds)
# so the ~200 MB fit need not be reloaded.
#
# Writes Supplementary/tables/hmsc_v1_chain_disagreement.csv (this worktree).
# Run from the worktree root: conda run -n r_env Rscript hmsc_converge/04_v1_disagreement_table.R
# =============================================================================

d <- readRDS("/home/daniel/Ptarmigan/models/hmsc_v2/hmsc_clr_fit_diag.rds")
ch <- paste0("ch", 1:4)

f5 <- d$fig5_chain[d$fig5_chain$OTU %in% c("OTU1", "OTU636", "OTU1369"), ]
lab_otu <- c(OTU1 = "Thelebolus (OTU1) Season Beta", OTU636 = "Sporormiella (OTU636) Season Beta",
             OTU1369 = "Dothideomycetes (OTU1369) Season Beta")
rows_otu <- data.frame(quantity = unname(lab_otu[f5$OTU]), f5[, ch], PSRF = f5$psrf,
                       check.names = FALSE)

g <- d$guild_chain
rows_g <- data.frame(quantity = c("Dung-saprotroph guild mean (CLR units)",
                                  "Unassigned/dark guild mean (CLR units)"),
                     rbind(g["dung_saprotroph", ch], g["dark_unassigned", ch]),
                     PSRF = NA, check.names = FALSE)

v <- d$vp_chain
rows_v <- data.frame(quantity = c("Variance share: Season", "Variance share: dropping",
                                  "Variance share: PCR replicate"),
                     rbind(v["Season", ch], v["Random: sample", ch], v["Random: pcr", ch]),
                     PSRF = NA, check.names = FALSE)

out <- rbind(rows_otu, rows_g, rows_v)
out[, ch] <- round(out[, ch], 2)
write.csv(out, "Supplementary/tables/hmsc_v1_chain_disagreement.csv", row.names = FALSE)
print(out, row.names = FALSE)
cat(sprintf("Superseded model: Beta PSRF median %.2f; Season Beta %% < 1.1 = %.1f\n",
            d$overall$beta_psrf_med, d$by_cov$pct_lt1.1[d$by_cov$covariate == "Seasonsummer"]))
