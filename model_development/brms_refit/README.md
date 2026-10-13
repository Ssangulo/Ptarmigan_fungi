# H2 / H3 brms refit — pilot memo (decision needed)

Branch `exp/brms-refit`, worktree `.worktrees/brms-refit`. Nothing on `main`, in `Scripts/`, in
the appendix or in the canonical `models/H2_*` / `H3_*` files was changed (md5 verified). All
fits and tables are in `/home/daniel/Ptarmigan/models/_brms_refit/`.

**Two picks are needed from Daniel: the H2 primary model and the H3 production model.**

---

## Starting point

Neither current model has a sampling problem (H2: 0 divergences, ESS 6–9k; H3: 0 divergences,
max Rhat 1.004), so longer chains would change nothing. The problems were structural:

- **H2:** the Gaussian family fails its posterior-predictive skewness check (p = 0.000).
- **H3:** the six plant predictors are exactly rank 3 (`decostand("rclr")` imputes zeros by
  rank-3 matrix completion); Season is one global coefficient, so the per-OTU Vaccinium slopes
  were per-OTU seasonality; one shared NB shape (0.024) predicts a mean count ~5× too high.

---

## H2 (§7.1) — `01_h2_arms.R` → `h2_scorecard.csv`, `h2_loo_compare.csv`

Reproduction guard passed exactly: 0.2891 / 0.9144 (Season + Year) and 0.3837 / 0.9634 (Season).
Every arm: 0 divergences, Rhat 1.00, minimum bulk ESS ≥ 2,700.

| Arm | P(slope>0), Season+Year | P(slope>0), Season | PPC skewness p | elpd vs Gamma (± SE) |
|---|---|---|---|---|
| Gaussian on z (current) | 0.914 | 0.963 | 0.000 | −13.1 ± 3.5 |
| **Gamma, log link** | **0.920** | **0.989** | 0.22 | 0 |
| Lognormal | 0.804 | 0.891 | 0.42 | −1.5 ± 0.9 |
| Skew-normal on z | 0.800 | 0.853 | 0.002 | −8.3 ± 3.0 |
| Student-t on z | 0.734 | 0.711 | 0.03 | −10.8 ± 3.4 |
| Gaussian, `(1 | Year)` | 0.939 | — | 0.000 | −13.4 ± 4.0 |

Gamma slope (log q1 per SD of plant richness): +0.215 [−0.09, 0.52] with Season + Year;
+0.335 [0.05, 0.63] with Season only.

**Covariate screen** (`h2_covariate_screen.csv`): no recorded variable explains the residuals —
sample mass p = 0.19, fungal depth p = 0.31; plate, extraction batch, sex, date, plant depth
p > 0.44.

**Two things that do not change with the family:**

1. **One sample carries the slope.** With `S_1_7_P1_8E` removed, P(slope > 0) is 0.50–0.59
   (Gaussian), 0.44–0.49 (lognormal), 0.52–0.57 (Gamma).
2. **Gamma and lognormal fit equally well** (difference 1.5 ± 0.9) and give 0.92 vs 0.80. Gamma
   models the arithmetic mean, so the one high sample counts for more; lognormal models the
   geometric mean.

**Recommendation: make Gamma (log link), Season + Year, the primary.** It is the only arm that
meets both priorities — support at least as high as now, and a model that fits (it passes all
three predictive checks and beats the Gaussian by 13.1 ± 3.5 elpd). It is also the family
Table S9 already uses, which was chosen on 2026-09-29 for unrelated reasons. Report lognormal
and the current Gaussian beside it, and keep the leave-one-out row: the honest sentence is
still "positive, carried by one sample, 0.95 cleared only without Year".

---

## H3 (§7.2) — `02_h3_data.R`, `03_h3_arms.R`, `04_h3_compare.R`

### Predictors (as decided)

Zeros → 0.1 % of diet, CLR over the full composition, z-scored. Four species (Betula,
*V. myrtillus*, *V. uliginosum*, *E. nigrum*): full rank, condition number 4.1. Three genera:
3.6. Predictors correlate ≥ 0.98 across replacement values 0.01–0.5 %. The count frame is
identical to script 6's (30 samples, 60 PCR replicates, 46 OTUs).

### Arms (pilot: 4 chains × 1500 / 750)

| Arm | Likelihood | Divergences | Max Rhat | Predicted mean, 95 % (obs 2,114) | Per-OTU zero fraction r | Verdict |
|---|---|---|---|---|---|---|
| current fit | NB, one shape | 0 | 1.00 | 1,730 – 3.3 M | 0.66 | fits badly |
| L1 | NB, one shape | 604 / 0 | 1.55 / 1.03 | 26 k – 49 M | 0.90 | fails |
| L2 | NB, per-OTU shape | 2 / 0 | 1.03 / 1.01 | 47 k – 5 G | 0.99 | fails |
| L4 | zero-inflated NB | 961 / 568 | 1.97 / 1.93 | — | — | did not converge |
| L3 | hurdle NB | 0 / 0 | 1.04 / 1.01 | 1,730 – 10.9 k | 0.99 | fits |
| L5 | hurdle NB, plants also in occurrence part | 0 / 0 | 1.02 / 1.01 | 1,770 – 10.0 k | 0.99 | fits |
| **L6** | L5 + Year in occurrence part | 0 | 1.02 | 1,734 – 10.8 k | 0.99 | fits |

(Pairs are per-OTU Season / per-OTU Season + Year.)

- **Plain NB fails for a structural reason:** with two replicates per sample × OTU cell, the
  cell random effect and NB overdispersion are not separable. The cell SD ran to ~6.4 and
  predicted means are hundreds of times too high.
- **Hurdle NB is the only family that converges and reproduces the data.** NB shape 1.25, cell
  SD 1.9 (the PCR-replicate control is now doing something; it was 0.16).
- **Against the pre-set criteria, no arm passes everything.** The hurdle arms pass the mean and
  zero-fraction checks and fail the library-size ceiling (predicted maximum ≤ largest library in
  13–18 % of draws, criterion ≥ 95 %; median predicted maximum 1.2–1.5 M against 569 k). At
  pilot length minimum bulk ESS is 200–340 (criterion 400); production length should clear it.
- **LOO could not be used as planned:** every arm has 128–316 unreliable points of 2,760.

**Deviations from the approved plan:** (1) the hurdle arm's zero part carries Season — without
it season-exclusive OTUs could not be seasonal; (2) L5 and L6 were added after the first round,
because a plain hurdle model only lets plants act on abundance given presence.

### What the hurdle arms say (pilot length — provisional)

- **Abundance given presence: no plant specificity.** Per-OTU slope SDs 0.14–0.28; no OTU
  deviates resolvably from the community slope in any arm.
- **Occurrence is where plants matter.** Per-OTU SD of the occurrence slope: Betula 0.90,
  *E. nigrum* 0.45, *V. myrtillus* 0.26, *V. uliginosum* 0.13. OTU-specific deviations with a
  95 % CrI excluding 0 (L6): 9 for Betula, 2 for *Empetrum*, 0 for either *Vaccinium*.
- ***Empetrum*:** OTU45 *Coleophoma* (plant-associated) and OTU2 (Myriangiales) are detected
  more often when *Empetrum* is in the diet (deviation −0.76 [−1.48, −0.14] and −0.62
  [−1.38, −0.04] on the P(zero) logit scale). *Coleophoma empetri* is an *Empetrum* leaf fungus:
  this is the ingested-passenger pattern in a plant-associated taxon. Same two OTUs at genus level.
- **Betula:** *Sporormiella* OTU636, *Coniochaeta* OTU1225 and four other OTUs are detected
  less often in birch-heavy diets; dark OTU1369 and OTU691 more often. This reads as a
  birch-diet ↔ winter-type community syndrome that persists within season.
- ***V. myrtillus*:** one community-wide effect (+0.58: nearly all OTUs detected less often),
  no OTU-specific signal. A richness effect, not specificity.
- **Robust to Year:** occurrence slopes correlate 0.98–0.99 between L5 (no Year) and L6.
- **The preregistered contrast is still null.** P(coprophilous more diffuse) is 0.48–0.50 on
  abundance and 0.36–0.39 on occurrence (L6); Wilcoxon p 0.72 / 0.85. H3 as preregistered
  remains not supported.
- **The Thelebolus–Vaccinium association is gone** in every hurdle arm.
- **Species vs genus:** Betula and *Empetrum* slopes agree (r 0.97–0.98); genus Vaccinium
  tracks *V. myrtillus* (0.81), not *V. uliginosum* (0.45).

### Recommendation

**Production model = L6, four species:** hurdle NB; global and per-OTU Season, Year and plant
slopes in both the abundance and the occurrence part; `(1 | Sample:OTU)`; library-size offset.
4 chains × 4000 / 2000, about 0.6 h per chain. Three-genus version as the sensitivity fit.
L5 (no Year in the occurrence part) is the simpler alternative and gives the same slopes.

If this goes ahead, §7.2 changes in substance: the result becomes "no abundance specificity;
occurrence-level associations with Betula and *Empetrum*, one of them a plausible passenger;
preregistered contrast null", and the "resolved" flag should be defined on OTU-specific
deviations, since a slope shared by all OTUs is not specificity.

---

## Still open / not done

- §7.3 (GLLVM) uses the same rank-3 matrix — decision deferred to after this memo.
- Nothing is ported into `6_diversity_analyses.R` or the appendix; that needs a separate go-ahead.
- Text defects seen, not fixed: §7.2 calls Vaccinium "winter browse" (it is the summer diet
  here; winter is near-pure Betula); §7.1 says "27 of the 27 samples are distinct individuals"
  (24 are); script 6 comments still say n = 30 for H2; §7.1's "primarily at the seasonal scale"
  (already flagged 2026-09-10).

## Files (`models/_brms_refit/`)

`h2_scorecard.csv` · `h2_loo_compare.csv` · `h2_covariate_screen.csv` · `h3_pilot_scorecard.csv` ·
`h3_pilot_contrast.csv` · `h3_pilot_per_otu.csv` · `h3_pilot_dev_resolved.csv` ·
`h3_pilot_S_vs_SY_slopes.csv` · `h3_pilot_species_vs_genus_slopes.csv` · `h3_pilot_loo_compare.csv` ·
per-arm `*_score.rds` (slope matrices) and fits `h2_*.rds`, `h3_pilot_*.rds`.

---

# Production results (2026-10-01, after Daniel's picks)

**Picks:** H2 primary = Gamma (log link), Season + Year, with Gamma Season-only as the second
model, replacing the Gaussian models in the appendix; a small family-choice table kept
(docs-only Table S12b). H3 = L6 (hurdle NB, Year in both parts), four species, plus the
three-genus sensitivity.

## H2 — ported (script 6 Section 3/3b + appendix §7.1, commits f41b608, 71d32ec)

Run from the worktree via the new `S6_OUT_ROOT` / `S6_SUPP_DIR` / `S6_STOP_AFTER_H2MECH`
overrides. Every upstream Part A output (Hill tables, Tables S8/S9 and their draws, the matched
frame) came out **byte-identical** to the committed appendix copies.

| Model (Gamma, log link) | Slope (log, per SD) | 95 % CrI | ×q1 per SD | P(slope>0) |
|---|---|---|---|---|
| Season + Year (primary) | +0.216 | [−0.10, 0.51] | 1.24 | **0.915** |
| Season only | +0.336 | [0.05, 0.63] | 1.40 | **0.989** |

**Note on 0.915 vs the pilot's 0.920:** script 6 centres the intercept prior on
log(median outcome) for every Hill order (needed for the q0/q2 sensitivity fits), where the
pilot used a fixed Normal(2.3, 1). That shift moves the Season + Year P by 0.005, about the size
of its Monte Carlo error. The Season + Year support is therefore the **same** as the old
Gaussian (0.914); what Gamma buys is a model that fits (skewness PPC p = 0.23 vs 0.000;
+12.9 ± 3.5 elpd) and a Season-only model whose 95 % CrI now excludes zero.

Hill-order sensitivity (Table S13, Gamma): q0 0.78 / 0.92, q1 0.91 / 0.99, q2 0.95 / 0.99
(Season + Year / Season only). Leave-one-out (Table S14, Gamma GLM): `S_1_7_P1_8E` is still the
most influential sample; all 27 n−1 refits non-significant under Season + Year.

## H3 — production fits done, NOT ported

`h3_prod_scorecard.csv`, `h3_prod_contrast.csv`, `h3_prod_dev_resolved.csv`, `h3_prod_per_otu.csv`.

| | Four species | Three genera |
|---|---|---|
| Divergences / max Rhat / min bulk ESS | 0 / 1.00 / 1,034 | 0 / 1.01 / 879 |
| Convergence criterion | **pass** | **pass** |
| Predicted mean 95 % (obs 2,114) | 1,784 – 8,555 | 1,687 – 9,931 |
| Per-OTU zero fraction r / max gap | 0.988 / 0.037 | 0.990 / 0.034 |
| Predicted max ≤ largest library | 13 % of draws (**fails**, criterion 95 %) | 15 % (**fails**) |
| Hours per chain | 0.51 | 0.46 |

Results match the pilot:

- **Abundance given presence:** no OTU deviates from the community plant slope (0 resolved).
- **Occurrence, Betula** (OTU SD 0.91): 8 OTUs deviate. Coprophiles *Sporormiella* OTU636,
  *Coniochaeta* OTU1225 and OTU1305, plus OTU483, OTU554 *Exobasidium*, OTU758, OTU97, are
  detected less often in birch-heavy diets; dark OTU1369 more often.
- **Occurrence, *E. nigrum*** (OTU SD 0.46): OTU45 *Coleophoma* −0.76 [−1.50, −0.17] and OTU2
  Myriangiales −0.63 [−1.36, −0.02], detected more often with *Empetrum* in the diet. Same two
  at genus level.
- ***V. myrtillus*:** community-wide occurrence slope +0.58 (fewer OTUs detected), OTU SD 0.25,
  no OTU-specific signal. *V. uliginosum*: nothing.
- **Preregistered contrast:** P(coprophilous more diffuse) 0.48 (abundance) / 0.37
  (occurrence); Wilcoxon p 0.45 / 0.72. **H3 as preregistered: not supported.**

## H3 ported + §7.3 rebuilt (2026-10-01, commits 78da398 and the appendix commit after it)

**Script 6 §3c** now fits the hurdle model (four diet species, plus the three-genus sensitivity)
and writes `H3_perOTU_plant_coef` (total slope AND OTU deviation, both parts), `H3_specificity_index`,
`H3_specificity_contrast` (per part), `H3_model_summary` (community slopes + deviation SDs),
`H3_model_checks`, `H3_species_vs_genus`. Re-running it from the worktree reproduced the
production fit exactly (same min ESS 1,034, same hyperparameters). Run control renamed:
`S6_STOP_AFTER=H2MECH|H3`.

**Script 7 §8** uses the same predictors; caches renamed `gllvm_H3sp_*` (the old `gllvm_H3_*`
caches can no longer reload by mistake); `S7_OUT_ROOT` / `S7_SUPP_DIR` added; new agreement
tables against script 6 (`gllvm_H3_vs_brms_by_plant`, `gllvm_H3_vs_brms_resolved`). The H1
sections reproduced all three canonical H1 tables byte for byte.

**Correction to the pilot memo:** I expected the full-rank predictors to tame the GLLVM's
exploding unpooled slopes. They did the opposite — median |slope| 3.7 / 90th pct 14.0 / max
28.8 (was 1.9 / 7.0 / 15.9). The cause is per-OTU fixed Season + Year + four plants on 30
samples with 80 % zeros; the rank-3 matrix had been constraining them. AIC now prefers the plant
model over the null by ≈ 220 (it preferred the null before). Contrast still null (P 0.39).
Agreement with the Bayesian model: occurrence Spearman 0.49 / 0.50 / 0.71 (Betula / E. nigrum /
V. uliginosum), ~0 for V. myrtillus and for abundance; of the 10 Bayesian-resolved associations
8 have the same GLLVM sign and 6 a GLLVM 95 % CI excluding 0 (incl. Coleophoma–Empetrum).

**Appendix:** §7.2 and §7.3 rewritten — methods, captions, tables and factual Result paragraphs.
Both Interpretation paragraphs are now labelled slots for Daniel; the pre-refit interpretations
are kept in the documentation build only. §7.2 Limitation kept, counts updated.

## Still needs a go-ahead

1. Interpretation prose for §7.1 (the flagged "seasonal scale" sentences), §7.2 and §7.3 — Daniel.
2. Merge `exp/brms-refit` → main, and the CLAUDE.md update (`CLAUDE_md_proposal.md`).
