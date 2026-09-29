# =============================================================================
# hmsc_converge/01_fit_arms.R  (branch exp/hmsc-converge -- experimental)
# Fit one HMSC abundance "arm" for the convergence study. Data assembly and the
# guild rule are copied VERBATIM from Scripts/10_hmsc.R Sections 1-2; only the
# response transform, prevalence filter and random-level priors differ by arm.
#
#   A  : CLR (pseudocount 1, over retained OTUs) -- the current Fig 5 estimand;
#        sample level nfMax 2 + PCR level nfMax 1
#   A0 : as A, sample level only
#   B  : rCLR over the FULL table, zeros -> NA (abundance conditional on
#        presence); OTUs in >= 10 PCR rows AND >= 5 droppings;
#        sample nfMax 2 + PCR nfMax 1
#   B0 : as B, sample level only
# YScale = TRUE in every arm.
#
# Run: HMSC_ARM=A HMSC_RUN_MODE=pilot OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
#        conda run -n r_env Rscript hmsc_converge/01_fit_arms.R
# Outputs ONLY to models/hmsc_v2/ -- never touches the canonical hmsc_* caches.
# =============================================================================

suppressMessages({
  library(Hmsc); library(coda); library(vegan); library(phyloseq)
})

ARM      <- Sys.getenv("HMSC_ARM", "A")
RUN_MODE <- Sys.getenv("HMSC_RUN_MODE", "pilot")
THIN     <- as.integer(Sys.getenv("HMSC_THIN", "10"))
stopifnot(ARM %in% c("A", "A0", "B", "B0"))
mc <- if (RUN_MODE == "pilot") {
  list(samples = 250, thin = 10, transient = 2500)
} else {
  list(samples = 1000, thin = THIN, transient = 500 * THIN)
}
nChains <- 4

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE))
  try(RhpcBLASctl::blas_set_num_threads(1L), silent = TRUE)

set.seed(20260929)
setwd("/home/daniel/Ptarmigan/trimmed/mergedPlates/")
load("eco_analysis.RData")
stopifnot(exists("funguild_otu"))
out_dir <- "/home/daniel/Ptarmigan/models/hmsc_v2"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
tag <- if (RUN_MODE == "pilot") sprintf("pilot_%s", ARM) else sprintf("prod_%s_thin%d", ARM, THIN)
fit_path  <- file.path(out_dir, paste0(tag, ".rds"))
meta_path <- file.path(out_dir, paste0(tag, "_meta.rds"))
cat(sprintf("Arm %s, mode %s: %d samples x thin %d (transient %d), %d chains -> %s\n",
            ARM, RUN_MODE, mc$samples, mc$thin, mc$transient, nChains, basename(fit_path)))

otu_mat_of <- function(ps) {
  m <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) m <- t(m)
  m
}

# ---- Data assembly (10_hmsc.R Section 1) ------------------------------------
ps   <- alldat_full$nopool
Ymat <- otu_mat_of(ps)
md   <- data.frame(as(sample_data(ps), "data.frame"), stringsAsFactors = FALSE)

if (ARM %in% c("A", "A0")) {
  keep_otu <- colSums(Ymat > 0) >= 5
  Yk <- Ymat[, keep_otu, drop = FALSE]
  Yresp <- vegan::decostand(Yk, method = "clr", pseudocount = 1)
  dimnames(Yresp) <- dimnames(Yk)
} else {
  # rCLR over the full table: log(x) - mean(log(nonzero x)) per row; zeros -> NA.
  # Equals the nonzero cells of vegan's robust Aitchison transform (Fig 2).
  L <- log(Ymat); L[Ymat == 0] <- NA
  Rclr <- L - rowMeans(L, na.rm = TRUE)
  n_rows  <- colSums(Ymat > 0)
  n_drops <- apply(Ymat > 0, 2, function(z) length(unique(md$Sample_ID_field[z])))
  keep_otu <- n_rows >= 10 & n_drops >= 5
  Yk <- Ymat[, keep_otu, drop = FALSE]
  Yresp <- Rclr[, keep_otu, drop = FALSE]
  stopifnot(all(is.na(Yresp) == (Yk == 0)))
}
prev_rows <- colSums(Yk > 0)
cat(sprintf("OTUs retained: %d; rows %d; %% response cells observed: %.1f\n",
            ncol(Yk), nrow(Yk), 100 * mean(!is.na(Yresp))))

XData <- data.frame(
  Season = factor(md$Season, levels = c("winter", "summer")),
  Year   = factor(md$Year),
  row.names = rownames(Yk)
)
studyDesign <- data.frame(
  sample = factor(md$Sample_ID_field),
  pcr    = factor(md$pcr_sample_id),
  row.names = rownames(Yk)
)
stopifnot(nlevels(studyDesign$pcr) == nrow(Yk))

rL_sample <- setPriors(HmscRandomLevel(units = levels(studyDesign$sample)), nfMax = 2, nfMin = 2)
rL_pcr    <- setPriors(HmscRandomLevel(units = levels(studyDesign$pcr)),    nfMax = 1, nfMin = 1)
if (ARM %in% c("A0", "B0")) {
  ranLevels <- list(sample = rL_sample)
  studyDesign <- studyDesign[, "sample", drop = FALSE]
} else {
  ranLevels <- list(sample = rL_sample, pcr = rL_pcr)
}

# ---- Guild trait (10_hmsc.R Section 2, verbatim rule) -----------------------
COPRO_GENERA <- c("Sordaria","Podospora","Cercophora","Chaetomium","Schizothecium",
                  "Preussia","Delitschia","Pilobolus","Ascobolus","Saccobolus",
                  "Sporormiella","Coprinopsis","Thelebolus","Coniochaeta")
strip_rank <- function(x) sub("^[a-z]__", "", x)
otu_ids   <- colnames(Yk)
guild_str <- funguild_otu$Guild[match(otu_ids, funguild_otu$OTU_ID)]
tax_all   <- data.frame(as(tax_table(alldat$nopool), "matrix"), stringsAsFactors = FALSE)
otu_gen   <- strip_rank(tax_all$Genus[match(otu_ids, rownames(tax_all))])
otu_gen   <- ifelse(otu_gen %in% c("", "NA"), NA, otu_gen)
is_copro <- (!is.na(guild_str) & grepl("Dung Saprotroph", guild_str)) |
            (!is.na(otu_gen)   & otu_gen %in% COPRO_GENERA)
is_plant <- !is_copro & !is.na(guild_str) &
            grepl("Endophyte|Plant Saprotroph|Plant Pathogen|Epiphyte", guild_str)
is_patho <- !is_copro & !is_plant & !is.na(guild_str) & grepl("Pathogen|Parasite", guild_str)
is_dark  <- !is_copro & !is_plant & !is_patho & is.na(otu_gen) & is.na(guild_str)
guild_class <- ifelse(is_copro, "dung_saprotroph",
               ifelse(is_plant, "plant_associated",
               ifelse(is_patho, "pathotroph",
               ifelse(is_dark,  "dark_unassigned", "other"))))
guild_lvls <- c("dung_saprotroph","plant_associated","pathotroph","dark_unassigned","other")
TrData <- data.frame(guild = factor(guild_class, levels = guild_lvls), row.names = otu_ids)
print(table(TrData$guild))

# ---- Model + fit ------------------------------------------------------------
m <- Hmsc(Y = Yresp, XData = XData, XFormula = ~Season + Year,
          TrData = TrData, TrFormula = ~guild,
          studyDesign = studyDesign, ranLevels = ranLevels,
          distr = "normal", YScale = TRUE)

saveRDS(list(arm = ARM, mode = RUN_MODE, mc = mc, prev_rows = prev_rows,
             guild = setNames(guild_class, otu_ids), genus = setNames(otu_gen, otu_ids),
             reads = colSums(Yk)),
        meta_path)

if (file.exists(fit_path)) stop("Fit exists, refusing to overwrite: ", fit_path)
t0 <- Sys.time()
m <- sampleMcmc(m, samples = mc$samples, thin = mc$thin, transient = mc$transient,
                nChains = nChains, nParallel = nChains,
                verbose = max(1, round((mc$transient + mc$samples * mc$thin) / 10)))
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
saveRDS(m, fit_path)
cat("saved", fit_path, "\n")
