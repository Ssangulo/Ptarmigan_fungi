# =============================================================================
# brms_refit/03_h3_arms.R  (branch exp/brms-refit -- experimental)
# Fit and score ONE H3 arm. Data build and scoring come from 02_h3_data.R.
#
#   H3_ARM    L1 : negbinomial, one shape            (current likelihood)
#             L2 : negbinomial, shape ~ 1 + (1 | OTU)
#             L3 : hurdle_negbinomial, hu ~ 1 + Season + (1 + Season || OTU)
#             L4 : zero_inflated_negbinomial, zi ~ 1 + (1 | OTU)
#   H3_STRUCT S  : per-OTU Season + plant slopes
#             SY : per-OTU Season + Year + plant slopes
#   H3_PRED   species (4: Betula, V. myrtillus, V. uliginosum, E. nigrum) | genus (3)
#   H3_MODE   pilot (4 x 1500/750) | prod (4 x 4000/2000)
#
# NOTE on L3: a hurdle model fits mu to the POSITIVE counts only, so an OTU that
# is absent in one season can only express that through hu. hu therefore carries
# Season (global + per OTU); with hu ~ (1 | OTU) alone the 15 season-exclusive
# OTUs would have no way to be seasonal. Plant slopes stay in mu, so under L3
# they are "abundance given presence" slopes -- a different estimand from L1/L2/L4.
#
# Run: H3_ARM=L1 H3_STRUCT=S H3_PRED=species H3_MODE=pilot \
#      OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
#        conda run -n r_env Rscript brms_refit/03_h3_arms.R
# Outputs ONLY to models/_brms_refit/ -- never touches the canonical H3 caches.
# =============================================================================

ARM    <- Sys.getenv("H3_ARM", "L1")
STRUCT <- Sys.getenv("H3_STRUCT", "S")
PRED   <- Sys.getenv("H3_PRED", "species")
MODE   <- Sys.getenv("H3_MODE", "pilot")
ADAPT  <- as.numeric(Sys.getenv("H3_ADAPT_DELTA", "0.95"))
stopifnot(ARM %in% c("L1","L2","L3","L4"), STRUCT %in% c("S","SY"),
          PRED %in% c("species","genus"), MODE %in% c("pilot","prod"))

here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "02_h3_data.R"))
options(width = 200)

out_dir <- "/home/daniel/Ptarmigan/models/_brms_refit"
tag     <- sprintf("h3_%s_%s_%s_%s", MODE, ARM, STRUCT, PRED)
mc      <- if (MODE == "pilot") list(iter = 1500, warmup = 750) else list(iter = 4000, warmup = 2000)

pg_z <- build_plant_predictors(PRED)
bl   <- build_long(pg_z); long <- bl$long; plant_var <- bl$plant_var
cat(sprintf("%s: %d rows, %d OTUs, %d samples; plant predictors (%s, cond %.2f): %s\n", tag,
            nrow(long), nlevels(long$OTU), length(unique(long$Sample_ID_field)), PRED,
            attr(pg_z, "cond"), paste(plant_var, collapse = ", ")))

plant_terms <- paste(plant_var, collapse = " + ")
otu_terms   <- paste(c("1", "Season", if (STRUCT == "SY") "Year", plant_terms), collapse = " + ")
mu_form <- as.formula(sprintf(
  "count ~ 1 + Season + Year + %s + (%s || OTU) + (1 | Sample_ID_field:OTU) + offset(log_libsize)",
  plant_terms, otu_terms))

form <- switch(ARM,
  L1 = bf(mu_form),
  L2 = bf(mu_form, shape ~ 1 + (1 | OTU)),
  L3 = bf(mu_form, hu ~ 1 + Season + (1 + Season || OTU)),
  L4 = bf(mu_form, zi ~ 1 + (1 | OTU)))
fam <- switch(ARM, L1 = negbinomial(), L2 = negbinomial(),
              L3 = hurdle_negbinomial(), L4 = zero_inflated_negbinomial())

# Priors: b / Intercept / plant-slope SDs verbatim from script 6; the OTU
# intercept, Season and Year SDs are widened so season-exclusive OTUs are
# reachable. Distributional parts (shape / hu / zi) keep brms defaults.
wide <- c("Intercept", "Seasonsummer", if (STRUCT == "SY") c("Year2023", "Year2024"))
pri <- c(set_prior("normal(0,1)",        class = "b"),
         set_prior("normal(0,2)",        class = "Intercept"),
         set_prior("student_t(3,0,0.5)", class = "sd"),
         do.call(c, lapply(wide, function(cf)
           set_prior("student_t(3,0,2.5)", class = "sd", group = "OTU", coef = cf))))

fit <- brm(form, data = long, family = fam, prior = pri,
           chains = 4, iter = mc$iter, warmup = mc$warmup, cores = 4,
           control = list(adapt_delta = ADAPT, max_treedepth = 12),
           init_r = 0.5, seed = 1, refresh = 250,
           file = file.path(out_dir, tag))

sc <- score_h3(fit, plant_var, tag)
saveRDS(sc, file.path(out_dir, paste0(tag, "_score.rds")))
cat("\n"); print(t(sc$card)); cat("\n"); print(sc$contrast, digits = 3)
cat("\nhyperparameters:\n"); print(sc$hyper, digits = 3, row.names = FALSE)
