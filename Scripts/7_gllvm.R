# =============================================================================
# 7_gllvm.R
# Generalised Linear Latent Variable Model (GLLVM) for differential abundance
# of fungal OTUs between winter and summer dung samples (H1)
#
# Model: NB counts ~ Season (winter ref) + offset(log libsize)
#        + random effects: (1|indivID_glm) [individual] + (1|pcr_replicate_id)
#          [PCR replicate pseudoreplication -- pcr_replicate_id groups the 2
#          PCR technical replicates that belong to the same dung sample, via
#          the native `sample` metadata column]
#
# Data: alldat_full[[1]] (== alldat_full$nopool), the full-complexity object
#       from 4_data_prep.R -- one row per PCR replicate (57 samples x 2 reps
#       = 114 rows), NOT PCR-collapsed/rarefied. PCR pseudoreplication is
#       handled explicitly here via the pcr_replicate_id random effect and a
#       read-count offset, per the canonical-object design decision in
#       4_data_prep.R (alldat_full is for GLLVM/HMSC only, for exactly this
#       reason). Do not swap in alldat/alldat.rfy here -- those are already
#       PCR-collapsed and have no per-replicate structure left.
#
# Requires objects from 4_data_prep.R: alldat_full (list: nopool/pool/pspool)
#
# Reference: 7_gllvm.R @ 65fbffa (Root_fungi_DADA2) -- that script modeled
# `habitat` (forest/paramo/subparamo) with site+individual random effects for
# an Andean root study. This study has no habitat/site/elevation gradient;
# the fixed effect of interest is Season (H1: seasonal turnover) and the two
# random effects requested are individual identity and PCR replicate. The
# hardcoded Andean genus list in the original Section 5 (Acephala,
# Hyaloscypha, ...) has no equivalent here and was replaced with a data-
# driven top-N by effect size -- swap in taxa of interest (e.g. Sporormiella/
# Sordariales, the E2 coprophilous screen) once you've inspected genus_table.
#
# SCOPE: this GLLVM is a Season-only model (formula = ~Season). It serves H1 at
# OTU resolution -- which lineages drive the winter/summer turnover, with model
# uncertainty (Section 3b) and a residual co-occurrence ordination (Section 3c).
# It does NOT test H3 (plant-fungus covariation): that needs plant predictors
# and is the deferred diet-matched Bayesian analysis. Do not relabel this H3.
#
# Outputs (canonical dirs): models/gllvm_model_comparison.csv,
# models/gllvm_perOTU_season_coef.csv, models/genus_table_season.{csv,xlsx};
# plots/gllvm_perOTU_season_coef.png, plots/gllvm_residual_ordination.png,
# plots/panel_AB_heatmap_barplot_season.png, plots/gllvm_fit_nb_2_2_residuals.png.
# =============================================================================

library(gllvm)
library(phyloseq)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(writexl)
library(scales)

setwd("/home/daniel/Ptarmigan/trimmed/mergedPlates/")

# --- Data source for the GLLVM section (IMPORTANT) ---------------------------
# The saved fits (models/fit_nb_2*.rds, fitted 2026-07-02) were built on the
# PRE-depth-fix full-complexity object: 57 samples / 114 PCR-replicate rows /
# 327 filtered OTUs. That object now survives only in the pre-fix backup
# workspace `eco_analysis.RData.bak_20260701_partial`. To keep everything in
# this script (the library-size offset, per-OTU abundance weights, and the
# Poisson baseline for the AIC table) internally consistent with those fits, we
# rebuild the model inputs from that backup and later ASSERT the filtered
# response matrix is identical to the saved fits' `$y`. Sections 2-5 of the
# appendix use the corrected current object; Section 6 carries an explicit
# caveat noting this GLLVM predates the depth-filter fix. If you refit the
# GLLVMs on the corrected data, switch this back to load("eco_analysis.RData").
.fit_env <- new.env()
load("eco_analysis.RData.bak_20260701_partial", envir = .fit_env)
alldat_full <- .fit_env$alldat_full

# Canonical output dirs (absolute); saved GLLVM fits live in models/.
out_dir   <- "/home/daniel/Ptarmigan/models"
plot_dir  <- "/home/daniel/Ptarmigan/plots"
model_dir <- "/home/daniel/Ptarmigan/models"
dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Headless plotting convention (see CLAUDE.md): a standing null device removes
# R's X11 fallback for ANY base-graphics call (gllvm's plot()/ordiplot()
# included); every figure is then written through an explicit png() device.
grDevices::pdf(NULL)

save_png <- function(path, plot_obj, width = 8, height = 8, res = 300) {
  grDevices::png(path, width = width, height = height, units = "in", res = res)
  on.exit(grDevices::dev.off())
  print(plot_obj)
}

# =============================================================================
# SECTION 1 — DATA PREPARATION
# =============================================================================

ps <- alldat_full[[1]]

Y  <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) Y <- t(Y)
md <- as(sample_data(ps), "data.frame")

stopifnot(nrow(Y) == nrow(md), all(rownames(Y) == rownames(md)))
stopifnot(all(c("Season","sample","indivID") %in% names(md)))

# ---- Fixed effect: Season (winter ref) ---------------------------------------
# NB: the "TimeBlock" column is a finer season x year-block code (S_1/S_2/S_3/
# W_1/W_2/W_3) -- "Season" is the clean 2-level winter/summer factor that
# matches H1 and is used as the fixed effect below.
season_ref <- "winter"
md$Season  <- relevel(droplevels(factor(md$Season)), ref = season_ref)

# ---- Random effects: individual + PCR replicate (via pcr_replicate_id) ------
# pcr_replicate_id groups the 2 PCR technical replicates of the same dung
# sample, which is what accounts for PCR pseudoreplication (analogous to how
# the reference script nested Unique_ID under site) -- sourced from the
# native `sample` metadata column.
# indivID is NA for ~17 of the 57 samples that lack a confirmed
# microsatellite ID; each gets a unique placeholder so it behaves as an
# unshared singleton random intercept rather than silently pooling unrelated
# unidentified individuals together.
md$pcr_replicate_id <- droplevels(factor(md$sample))
md$indivID_glm  <- ifelse(is.na(md$indivID) | md$indivID == "",
                           paste0("unk_", md$sample), md$indivID)
md$indivID_glm  <- droplevels(factor(md$indivID_glm))

message("PCR replicates per Season:"); print(table(md$Season))
message("Samples (pcr_replicate_id) per Season:")
print(table(unique(md[, c("pcr_replicate_id","Season")])$Season))
message("Individuals (incl. unidentified singletons): ", nlevels(md$indivID_glm))

# ---- Library size offset ----------------------------------------------------
md$libsize     <- rowSums(Y)
md$log_libsize <- log(md$libsize)
print(summary(md$libsize))

# ---- Response matrix: taken directly from the saved fits --------------------
# We deliberately do NOT re-derive the rare-taxa filter here. The saved fits
# carry 327 OTUs, but the current thresholds (below, kept only for provenance)
# yield ~280 from this object -- the thresholds were changed after the Jul-2
# fitting run. Re-filtering would silently desynchronise the modelling matrix
# from the loaded coefficients. Instead we adopt the fits' own stored `$y` as
# the modelling matrix in Section 2, so the response, per-OTU abundance weights,
# and Poisson baseline all match the coefficients exactly. The Season design,
# the random-effect grouping, and the library-size offset are recovered from
# this (backup) object and aligned to the fits' row order once the fits load.
#   Original (superseded) thresholds, for the record:
#     min_total_reads = 1500; presence_threshold = 200; min_prevalence = 2

# =============================================================================
# SECTION 2 — MODEL FITTING (NB + Poisson comparison)
# =============================================================================

# The two negative-binomial fits are expensive and already persisted from the
# original fitting run (each took a long time under VA), so we LOAD rather than
# refit. To reproduce the fits from scratch, run the archived gllvm() calls
# (see git history of this script) with sd.errors=TRUE.
#   fit_nb_2   = NB, individual + PCR-rep RE, 1 latent variable
#   fit_nb_2_2 = NB, individual + PCR-rep RE, 2 latent variables  <-- PRIMARY /
#                reported model: 2 LVs give a 2-D residual (co-occurrence)
#                ordination and are AIC-preferred (see model_comparison below).
fit_nb_2   <- readRDS(file.path(model_dir, "fit_nb_2.rds"))
fit_nb_2_2 <- readRDS(file.path(model_dir, "fit_nb_2_2.rds"))

# Adopt the fits' own stored response as the modelling matrix (exact 327 OTUs x
# 114 rows), and align the backup-derived design + offset to the fits' row
# order. The backup reproduces the samples and per-sample library sizes the fits
# were trained on; matching by rowname recovers the exact offset/grouping.
stopifnot(identical(rownames(fit_nb_2$y), rownames(fit_nb_2_2$y)),
          identical(colnames(fit_nb_2$y), colnames(fit_nb_2_2$y)))
ord <- match(rownames(fit_nb_2_2$y), rownames(md))
if (anyNA(ord))
  stop("Backup object does not cover the saved fits' rows -- cannot recover offset/design.")
md          <- md[ord, , drop = FALSE]
Yf          <- fit_nb_2_2$y
Xseason     <- data.frame(Season = md$Season)
studyDesign <- md[, c("indivID_glm", "pcr_replicate_id")]
stopifnot(identical(rownames(Yf), rownames(md)))
message(sprintf("Adopted saved fits' response: %d rows x %d OTUs; offset + design recovered from backup (libsize %d-%d).",
                nrow(Yf), ncol(Yf), min(md$libsize), max(md$libsize)))

# ---- Poisson baseline (num.lv=0, PCR-replicate RE only) ---------------------
# Not persisted from the original run; fit once (cheapest model) and cache so
# reruns are fast. Justifies the negative-binomial family over Poisson.
pois_path <- file.path(model_dir, "fit_poisson.rds")
if (file.exists(pois_path)) {
  fit_poisson <- readRDS(pois_path)
} else {
  fit_poisson <- gllvm(
    y=Yf, X=Xseason, formula=~Season,
    family=poisson(),
    offset=md$log_libsize,
    row.eff=~(1|pcr_replicate_id),
    studyDesign=studyDesign,
    num.lv=0, sd.errors=FALSE, method="VA"
  )
  saveRDS(fit_poisson, pois_path)
}

# ---- Model comparison (Table S-GLLVMa) --------------------------------------
model_comparison <- data.frame(
  model  = c("Poisson (num.lv=0, PCR-rep RE)",
             "Negative binomial (num.lv=1, indiv + PCR-rep RE)",
             "Negative binomial (num.lv=2, indiv + PCR-rep RE)"),
  object = c("fit_poisson", "fit_nb_2", "fit_nb_2_2"),
  logLik = c(as.numeric(logLik(fit_poisson)),
             as.numeric(logLik(fit_nb_2)),
             as.numeric(logLik(fit_nb_2_2))),
  AIC    = c(AIC(fit_poisson), AIC(fit_nb_2), AIC(fit_nb_2_2)),
  stringsAsFactors = FALSE
)
model_comparison$dAIC <- model_comparison$AIC - min(model_comparison$AIC)
print(model_comparison)
write.csv(model_comparison, file.path(out_dir, "gllvm_model_comparison.csv"),
          row.names = FALSE)

# ---- Residual diagnostics (primary model: fit_nb_2_2) -----------------------
png(file.path(plot_dir, "gllvm_fit_nb_2_2_residuals.png"),
    width=3600, height=2400, res=300)
par(mfrow=c(3,2), mar=c(4,4,2,1), ask=FALSE)
for (i in 1:5) plot(fit_nb_2_2, which=i, var.colors=1)
plot.new()
dev.off()

# =============================================================================
# SECTION 3 — COEFFICIENT EXTRACTION AND TAXONOMY ANNOTATION
# =============================================================================

# Coefficient name gllvm derives from model.matrix(~Season), e.g. "Seasonsummer"
coef_col <- paste0("Season", setdiff(levels(md$Season), season_ref))

coef_nb <- as.data.frame(coef(fit_nb_2_2)$Xcoef) %>% rownames_to_column("taxon")
tax     <- as.data.frame(tax_table(ps))          %>% rownames_to_column("taxon")

# Keep only genus-assigned taxa, exclude Incertae sedis
keep_taxa    <- grepl("^g__", tax$Genus) & !grepl("incertae_sedis", tax$Genus, ignore.case=TRUE)
tax_filt     <- tax[keep_taxa, ]
coef_nb_filt <- coef_nb %>% filter(taxon %in% tax_filt$taxon)

annot_raw <- coef_nb_filt %>%
  left_join(tax_filt %>% select(taxon, Genus, Family, Order, Class, Phylum), by="taxon")

# Clean rank prefixes and trailing rank words
clean_tax <- function(x) {
  x <- gsub("^\\w__","", as.character(x))
  x <- gsub(" (Family|Order|Class|Phylum)$","", x)
  ifelse(is.na(x)|x=="","Unassigned",x)
}
annot_raw <- annot_raw %>% mutate(Genus=clean_tax(Genus), Order=clean_tax(Order))

# Per-OTU abundance weights
otu_abund         <- colSums(Yf)
annot_raw$w_abund <- pmax(1, otu_abund[match(annot_raw$taxon, names(otu_abund))])

# beta_season = log fold-change, summer vs winter (positive = higher in summer)
annot_raw$beta_season <- annot_raw[[coef_col]]

# =============================================================================
# SECTION 3b — PER-OTU SEASON COEFFICIENTS WITH THE MODEL'S OWN UNCERTAINTY
# =============================================================================
# The genus-level summaries below propagate among-OTU spread via a bootstrap.
# Here we instead use each OTU's own Wald standard error from the fitted model
# (sd.errors=TRUE) to give per-OTU 95% CIs and BH-adjusted significance -- the
# statistically strongest use of the GLLVM for "which lineages drive H1". All
# modelled OTUs are included (not only genus-assigned), so dark taxa (e.g. the
# winter-persistent OTU1369) are retained.

Xco <- as.matrix(coef(fit_nb_2_2)$Xcoef)
if (is.null(fit_nb_2_2$sd) || is.null(fit_nb_2_2$sd$Xcoef))
  stop("fit_nb_2_2 has no stored Xcoef standard errors -- refit with sd.errors=TRUE.")
Xse <- as.matrix(fit_nb_2_2$sd$Xcoef)
stopifnot(coef_col %in% colnames(Xco), coef_col %in% colnames(Xse))

otu_coef <- data.frame(
  taxon       = rownames(Xco),
  beta_season = Xco[, coef_col],
  se_season   = Xse[, coef_col],
  row.names   = NULL, stringsAsFactors = FALSE
)
otu_coef$ci_lo <- otu_coef$beta_season - 1.96 * otu_coef$se_season
otu_coef$ci_hi <- otu_coef$beta_season + 1.96 * otu_coef$se_season
otu_coef$z     <- otu_coef$beta_season / otu_coef$se_season
otu_coef$p     <- 2 * pnorm(-abs(otu_coef$z))
otu_coef$p_adj <- p.adjust(otu_coef$p, method = "BH")
otu_coef$sig   <- otu_coef$p_adj < 0.05

# Full taxonomy + abundance weight for every modelled OTU
tax_full <- as.data.frame(tax_table(ps)) %>% rownames_to_column("taxon")
otu_coef <- otu_coef %>%
  left_join(tax_full %>% select(taxon, Phylum, Class, Order, Family, Genus, Species),
            by = "taxon") %>%
  mutate(across(c(Phylum, Class, Order, Family, Genus, Species), clean_tax))
otu_coef$w_abund <- pmax(1, otu_abund[match(otu_coef$taxon, names(otu_abund))])

# Flag near-separated / degenerate estimates: an OTU essentially present in only
# one season pushes the Season coefficient toward +/-Inf with a tiny or
# exploding SE (|beta| in the hundreds, or CI widths of 1e8+). These are real
# "on/off" seasonal switches, but their point estimates/CIs are not
# quantitatively interpretable -- we flag them, keep them in the full table, and
# exclude them from the caterpillar panel (scale) below.
otu_coef$ci_width   <- otu_coef$ci_hi - otu_coef$ci_lo
otu_coef$separation <- !is.finite(otu_coef$z) | otu_coef$se_season <= 0.1 |
                       !is.finite(otu_coef$ci_width) | otu_coef$ci_width > 60 |
                       abs(otu_coef$beta_season) > 50

message(sprintf("Per-OTU Season coefficients: %d OTUs, %d significant at BH q<0.05 (%d up in summer, %d up in winter); %d flagged near-separated",
                nrow(otu_coef), sum(otu_coef$sig, na.rm = TRUE),
                sum(otu_coef$sig & otu_coef$beta_season > 0, na.rm = TRUE),
                sum(otu_coef$sig & otu_coef$beta_season < 0, na.rm = TRUE),
                sum(otu_coef$separation, na.rm = TRUE)))

# Table S-GLLVMb (full per-OTU coefficients, all modelled OTUs)
otu_coef_out <- otu_coef %>%
  select(taxon, Phylum, Class, Order, Family, Genus, Species,
         beta_season, se_season, ci_lo, ci_hi, z, p, p_adj, sig,
         separation, w_abund) %>%
  arrange(desc(beta_season))
write.csv(otu_coef_out, file.path(out_dir, "gllvm_perOTU_season_coef.csv"), row.names = FALSE)

# Figure S-GLLVMa: caterpillar of the 30 most ABUNDANT well-estimated OTUs, each
# with its own model 95% CI. Abundance ranking (rather than |z|, which would
# surface unstable near-separated rare OTUs) keeps the panel to the ecologically
# meaningful players -- including the season-dominant taxa -- on an interpretable
# axis; near-separated OTUs are omitted here but retained in Table S-GLLVMb.
n_show <- 30
cat_df <- otu_coef %>% filter(!separation) %>% slice_max(w_abund, n = n_show)

# Collapse Order to the most abundant few + "Other" so the colour legend stays
# within the Dark2 palette (<= 8 levels).
top_orders <- cat_df %>% count(Order, wt = w_abund, sort = TRUE) %>%
  slice_head(n = 7) %>% pull(Order)
cat_df$Order_col <- ifelse(cat_df$Order %in% top_orders, cat_df$Order, "Other")
cat_df$label <- with(cat_df, ifelse(Genus == "Unassigned",
                                     paste0(taxon, " (", Order, ")"),
                                     paste0(taxon, " (", Genus, ")")))
p_cat <- ggplot(cat_df, aes(x = beta_season, y = reorder(label, beta_season), color = Order_col)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.3, linewidth = 0.6) +
  geom_point(aes(shape = sig), size = 2.4, fill = "white") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                     name = "BH q < 0.05", labels = c(`TRUE` = "yes", `FALSE` = "no")) +
  scale_color_brewer(palette = "Dark2", name = "Order") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(), axis.title.y = element_blank(),
        legend.position = "right") +
  labs(x = expression(paste("Season coefficient ", beta,
                            " (GLLVM log-scale; positive = higher in summer)")),
       title = "Season response of the 30 most abundant OTUs (fit_nb_2_2, model 95% CIs)")
save_png(file.path(plot_dir, "gllvm_perOTU_season_coef.png"), p_cat, width = 9, height = 8)

# =============================================================================
# SECTION 3c — RESIDUAL CO-OCCURRENCE ORDINATION (2 latent variables)
# =============================================================================
# The 2 latent variables of fit_nb_2_2 give a model-based ordination of residual
# OTU co-occurrence AFTER accounting for Season and the random effects -- a
# model-based complement to the unconstrained robust-Aitchison PCoA (Section 4
# of the appendix). Sites coloured by Season; top species loadings shown.

season_pal <- c(winter = "#0072B2", summer = "#D55E00")
site_col   <- season_pal[as.character(md$Season)]

png(file.path(plot_dir, "gllvm_residual_ordination.png"),
    width = 8, height = 7, units = "in", res = 300)
# symbols=TRUE draws sites as coloured points (default is to print each row's
# name as text, which overplots into an unreadable hairball with these long
# sample IDs). ind.spp shows only the few strongest species loadings.
gllvm::ordiplot(fit_nb_2_2, biplot = TRUE, symbols = TRUE, ind.spp = 8,
                s.colors = site_col, s.cex = 1.2,
                main = "GLLVM residual ordination (2 LV, after Season)")
legend("topright", legend = names(season_pal), col = season_pal,
       pch = 1, bty = "n", title = "Season")
dev.off()

# =============================================================================
# SECTION 4 — GENUS- AND ORDER-LEVEL SUMMARIES
# =============================================================================

min_n <- 5

genus_summary <- annot_raw %>%
  group_by(Order, Genus) %>%
  summarise(n_OTUs=n(),
            mean_beta_season  = mean(beta_season, na.rm=TRUE),
            q90_beta_season   = quantile(beta_season, 0.90, na.rm=TRUE),
            max_beta_season   = max(beta_season, na.rm=TRUE),
            wmean_beta_season = weighted.mean(beta_season, w=w_abund, na.rm=TRUE),
            .groups="drop") %>%
  filter(n_OTUs >= min_n) %>%
  arrange(desc(q90_beta_season), desc(wmean_beta_season))

order_summary <- annot_raw %>%
  group_by(Order) %>%
  summarise(n_OTUs=n(),
            median_beta_season = median(beta_season, na.rm=TRUE),
            q90_beta_season    = quantile(beta_season, 0.90, na.rm=TRUE),
            max_beta_season    = max(beta_season, na.rm=TRUE),
            wmean_beta_season  = weighted.mean(beta_season, w=w_abund, na.rm=TRUE),
            .groups="drop") %>%
  filter(n_OTUs >= min_n) %>%
  arrange(desc(q90_beta_season), desc(wmean_beta_season))

head(genus_summary, 15)
head(order_summary, 15)

# ---- Full genus table with abundance-weighted stats and % enrichment --------
min_abs_beta_season <- 2
min_total_w_abund    <- 20000
min_n_OTUs           <- 3

w_median <- function(x, w) {
  x <- as.numeric(x);  w <- as.numeric(w)
  ok <- is.finite(x) & is.finite(w) & !is.na(x) & !is.na(w)
  x <- x[ok];  w <- w[ok]
  if (!length(x)) return(NA_real_)
  x <- x[order(x)];  w <- w[order(x)]
  x[which(cumsum(w)/sum(w) >= 0.5)[1]]
}

genus_table <- annot_raw %>%
  dplyr::group_by(Order, Genus) %>%
  dplyr::summarise(
    n_OTUs=dplyr::n(),
    total_w_abund=sum(w_abund,na.rm=TRUE), mean_w_abund=mean(w_abund,na.rm=TRUE),
    median_beta_season=median(beta_season,na.rm=TRUE),
    mean_beta_season  =mean(beta_season,  na.rm=TRUE),
    q90_beta_season   =quantile(beta_season,0.90,na.rm=TRUE),
    max_beta_season   =max(beta_season,   na.rm=TRUE),
    w_median_beta_season=w_median(beta_season, w_abund),
    w_mean_beta_season  =stats::weighted.mean(beta_season,w_abund,na.rm=TRUE),
    .groups="drop"
  ) %>%
  dplyr::mutate(
    perc_enriched_season   = (exp(median_beta_season)  -1)*100,
    perc_enriched_season_w = (exp(w_median_beta_season)-1)*100
  ) %>%
  dplyr::filter(
    total_w_abund >= min_total_w_abund &
      (abs(median_beta_season) >= min_abs_beta_season | n_OTUs >= min_n_OTUs)
  ) %>%
  dplyr::arrange(dplyr::desc(w_median_beta_season))

head(genus_table, 20)

write_xlsx(genus_table, file.path(out_dir, "genus_table_season.xlsx"))
write.csv( genus_table, file.path(out_dir, "genus_table_season.csv"))

# =============================================================================
# SECTION 5 — HEATMAP: weighted mean β per genus (summer vs winter)
# =============================================================================
# Top/bottom-N by |effect size| in place of the reference script's hardcoded
# Andean genus list (no equivalent a priori list here yet).

top_n <- 15
genera_keep <- genus_table %>%
  slice_max(order_by = abs(w_mean_beta_season), n = top_n) %>%
  pull(Genus)

genus_label_levels <- genus_table %>%
  filter(Genus %in% genera_keep) %>%
  arrange(desc(w_mean_beta_season)) %>%
  pull(Genus)

heatmap_df <- genus_table %>%
  filter(Genus %in% genera_keep) %>%
  select(Order, Genus, total_w_abund, w_mean_beta_season) %>%
  mutate(
    habitat_label = "Summer vs winter",
    Genus = factor(Genus, levels = genus_label_levels)
  )

p_heat <- ggplot(heatmap_df, aes(x=habitat_label, y=Genus, fill=w_mean_beta_season)) +
  geom_tile(color="white", linewidth=0.3) +
  scale_fill_gradient2(low="#3B82F6", mid="white", high="#EF4444", midpoint=0,
                       limits=c(-35,35), oob=squish,
                       name=expression(paste(beta," (weighted mean log fold-change)"))) +
  theme_minimal(base_size=12) +
  theme(panel.grid=element_blank(), axis.title=element_blank(),
        axis.text.y=element_text(size=9), axis.text.x=element_text(size=11,angle=25,hjust=1),
        legend.position="right")

save_png(file.path(plot_dir, "heatmap_genus_season.png"),
         p_heat, width=6, height=10, res=900)

# =============================================================================
# SECTION 6 — BAR PLOT: weighted mean β + bootstrapped 95% CI
# =============================================================================

genus_labels <- genus_table %>%
  filter(Genus %in% genera_keep) %>%
  mutate(Genus_label=sprintf("%s (n=%d, %dk)", Genus, n_OTUs, round(total_w_abund/1000))) %>%
  select(Genus, Genus_label)

boot_wmean <- function(beta, w, R=500) {
  beta <- as.numeric(beta);  w <- as.numeric(w)
  ok   <- is.finite(beta) & is.finite(w) & w>0
  beta <- beta[ok];  w <- w[ok]
  if (length(beta)<2) return(c(lo=NA_real_, hi=NA_real_))
  prob  <- w/sum(w)
  boots <- replicate(R, { idx <- sample(seq_along(beta),replace=TRUE,prob=prob)
  ww <- w[idx]/sum(w[idx]); weighted.mean(beta[idx],ww) })
  quantile(boots, c(0.025,0.975), na.rm=TRUE)
}

genus_ci <- annot_raw %>%
  filter(Genus %in% genera_keep) %>%
  group_by(Genus, Order) %>%
  summarise(w_mean_beta=weighted.mean(beta_season,w_abund,na.rm=TRUE), {
    ci <- boot_wmean(beta_season, w_abund, R=500)
    tibble(ci_lo=ci[1], ci_hi=ci[2])
  }, .groups="drop") %>%
  left_join(genus_labels, by="Genus")

x_lim <- 100
genus_ci <- genus_ci %>%
  mutate(plot_x =pmax(pmin(w_mean_beta,x_lim),-x_lim),
         plot_lo=pmax(pmin(ci_lo,      x_lim),-x_lim),
         plot_hi=pmax(pmin(ci_hi,      x_lim),-x_lim),
         extreme_pos=ci_hi>x_lim, extreme_neg=ci_lo<(-x_lim))

p_bar_w <- ggplot(genus_ci, aes(x=plot_x, y=reorder(Genus_label,-plot_x), color=Order)) +
  geom_errorbarh(aes(xmin=plot_lo,xmax=plot_hi), height=0.3, linewidth=0.8) +
  geom_point(size=3) +
  geom_segment(data=subset(genus_ci,extreme_pos),
               aes(x=x_lim-15,xend=x_lim+5,y=Genus_label,yend=Genus_label),
               arrow=arrow(length=unit(0.25,"cm"),type="closed",ends="last"),
               linewidth=1, color="black", inherit.aes=FALSE) +
  geom_vline(xintercept=0, linetype="dashed", color="grey40") +
  scale_color_brewer(palette="Dark2") +
  coord_cartesian(xlim=c(-x_lim,x_lim+10)) +
  theme_minimal(base_size=12) +
  theme(panel.grid.major.y=element_blank(), axis.title.y=element_blank(),
        legend.position="right", axis.title.x=element_text(margin=margin(t=8))) +
  labs(x=expression(paste("Weighted mean ",beta," (summer vs winter)")), color="Order")

save_png(file.path(plot_dir, "barplot_genus_season.png"),
         p_bar_w, width=8, height=10, res=900)

# =============================================================================
# SECTION 7 — PANEL PLOT (Heatmap + Bar side by side)
# =============================================================================

p_heat_clean <- p_heat +
  labs(title=NULL,subtitle=NULL) +
  theme(plot.title=element_blank(), plot.subtitle=element_blank(),
        legend.title=element_text(size=8), legend.text=element_text(size=7),
        legend.key.height=unit(0.4,"cm"), legend.key.width=unit(0.25,"cm"),
        legend.margin=margin(t=0,r=2,b=0,l=0))

p_bar_clean <- p_bar_w +
  labs(title=NULL,subtitle=NULL) +
  theme(plot.title=element_blank(), plot.subtitle=element_blank(),
        legend.title=element_text(size=8), legend.text=element_text(size=7),
        legend.key.height=unit(0.4,"cm"), legend.key.width=unit(0.4,"cm"),
        legend.margin=margin(t=0,r=2,b=0,l=0),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=10, margin=margin(t=0)))

panel_AB <- p_heat_clean + p_bar_clean +
  plot_layout(widths=c(0.9,1.3)) +
  plot_annotation(tag_levels="A", tag_suffix="")

save_png(file.path(plot_dir, "panel_AB_heatmap_barplot_season.png"),
         panel_AB, width=11, height=6.5, res=900)
