# model_development/: model-checking studies

These are two studies that checked and replaced a fitted model. The pipeline in `Scripts/` does
not run them. They are kept because the appendix cites them as the evidence for its model
choices, and because one script here still produces an appendix table.

Each folder's `README.md` is the memo written when that study was done, dated and in the
working register of the time. It refers to the branch and to the files in place then. The
fitted pilot models are large and are kept outside the repository, on the analysis server:
`models/hmsc_v2/` and `models/_brms_refit/`.

## [`hmsc_convergence/`](hmsc_convergence/) (2026-09-29): HMSC abundance model

**What was wrong.** The four MCMC chains of the first CLR-abundance HMSC fit each settled on a
different answer. The cause was that the uncapped latent factors at the dropping level could
absorb the Season contrast, which is constant within a dropping. Longer chains did not help.

**What was decided.** The study compared several model structures. The model now in
`Scripts/10_hmsc.R` (CLR averaged per dropping, two latent factors) passed every per-chain
check.

**Cited in.** Appendix §9.1 (Table S-HMSCb2) and §10.3 (Figure 3).
`04_v1_disagreement_table.R` is the live producer of Table S-HMSCb2, and the repro harness runs
it from the repository root:

```bash
Rscript model_development/hmsc_convergence/04_v1_disagreement_table.R
```

Scripts: `00_diagnostics.R` (per-chain convergence test), `01_fit_arms.R` (fit one candidate
model), `02_compare.R`, `03_fig5_draft.R` (draft Figure 3 from a candidate fit),
`04_v1_disagreement_table.R`, `05_fill_numbers.R` (one-off fill-in of appendix numbers).

## [`brms_refit/`](brms_refit/) (2026-10-01): H2-mechanism and H3 brms models

**What was wrong.** Neither model had a sampling problem; the problems were structural. In the
H2 diet-richness model, the Gaussian family failed its posterior-predictive check. In the H3
plant–fungus model, the six plant predictors were exactly rank 3, because `decostand("rclr")`
imputes zeros by rank-3 matrix completion. Season was also a single global term, so per-OTU
seasonality leaked into the plant slopes.

**What was decided.** H2 now uses a Gamma (log-link) model. H3 now uses a hurdle
negative-binomial model with full-rank diet-species predictors and per-OTU Season and Year
terms. Both live in `Scripts/6_diversity_analyses.R`.

**Cited in.** Appendix §7.1 (Table S12b) and §7.2.

Scripts: `01_h2_arms.R` (H2 family comparison), `02_h3_data.R` (shared data build),
`03_h3_arms.R` (fit one H3 candidate), `04_h3_compare.R` (H3 comparison tables).
