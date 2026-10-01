# Proposed CLAUDE.md changes (exp/brms-refit) — apply only on Daniel's go-ahead

## Replace in CURRENT FOCUS

Drop the "H3 v2 work is UNCOMMITTED" note (stale) and add:

> H2 mechanism (§7.1) refit as Gamma (log link) on raw q1 — Season + Year primary
> (P = 0.915), Season-only second (P = 0.989). H3 (§7.2) production model chosen = hurdle NB
> with Season, Year and four diet-species slopes in both parts (`brms_refit/`, memo
> `MEMO_brms_refit.md`); port into script 6 §3c + appendix pending.

## New DECISIONS entry

- 2026-10-01: **brms H2/H3 refit** (branch `exp/brms-refit`, memo `brms_refit/MEMO_brms_refit.md`).
  Neither old model had a SAMPLING problem — longer chains change nothing; problems were structural.
  **H2:** Gaussian-on-z failed the PPC skewness check (p = 0); now Gamma(log) on raw q1, same
  family as Table S9. Season + Year P = 0.915 (unchanged in practice), Season-only 0.989
  (CrI excludes 0). Under EVERY family, dropping `S_1_7_P1_8E` gives P ≈ 0.5. Lognormal fits as
  well as Gamma (Δelpd −1.5 ± 0.9) and gives 0.81 — docs-only Table S12b.
  **H3:** old predictors were EXACTLY RANK 3 (see gotcha); Season was global, so the Vaccinium
  slope SD (1.66) was per-OTU seasonality; one NB shape 0.024 predicted means ~5× too high.
  New predictors: zeros → 0.1 % of diet (multiplicative replacement), CLR over the full plant
  composition, Betula / *V. myrtillus* / *V. uliginosum* / *E. nigrum* (3-genus sensitivity).
  Only the HURDLE NB converges and reproduces the data (plain NB: cell RE vs overdispersion not
  separable with 2 reps per cell; ZINB: Rhat ~1.9). Plant signal is in OCCURRENCE only:
  *Coleophoma* OTU45 ↔ *Empetrum*; coprophiles less often detected with Betula. Preregistered
  contrast still null; Thelebolus–Vaccinium gone. Fails the library-size-ceiling PPC (NB tail).
  **Winter diet variation is almost all 2022** (10 of 12 winter 2023/24 samples near-pure birch).

## New GOTCHAs

- **`vegan::decostand(x, "rclr")` imputes zeros by rank-3 matrix completion** (`ropt = 3`) in
  vegan 2.7 — on a sparse matrix the result is EXACTLY rank 3 whatever its width. Never use it
  to build regression predictors; check `qr(X)$rank`.
- **LOO is unusable for the H3 count models** (130–320 of 2,760 points with Pareto k > 0.7);
  compare them by predictive checks instead.
- **Script 6 can now run from a worktree:** `S6_OUT_ROOT`, `S6_SUPP_DIR`,
  `S6_STOP_AFTER_H2MECH=1` (same pattern as 10_hmsc.R). Defaults are unchanged.
- In a shell wait loop, `pgrep -f "<script name>"` matches the loop's own command line and never
  exits — match on something the loop does not contain.
