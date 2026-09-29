# HMSC abundance model — convergence study (2026-09-29)

Branch `exp/hmsc-converge`. Nothing on `main`, in `Supplementary/` or in the canonical
`models/hmsc_*` caches was touched (md5s verified). Fits and tables are in `models/hmsc_v2/`,
drafts in `plots/hmsc_v2/`.

## 1. What was wrong with the published CLR model

All four chains converged, but to **different answers**. The pooled posterior in Figure 5 and
Table S-HMSCc is therefore a mixture of disagreeing chains (`diag_canonical_clr.log`):

| | ch1 | ch2 | ch3 | ch4 |
|---|---|---|---|---|
| OTU1 *Thelebolus* Season β | 0.36 | 3.03 | 0.69 | −0.03 |
| Dung-saprotroph guild mean | 0.16 | 1.10 | 0.28 | 0.02 |
| VP Season | 0.11 | 0.28 | 0.06 | 0.05 |
| VP dropping / PCR | .57/.27 | .18/.42 | .45/.42 | .78/.14 |

- Only **10%** of Season β have PSRF < 1.1. The "40%" in Table S-HMSCb is diluted by the Year terms.
- The effective sample size **looked fine** (OTU1 ESS 2,157 at PSRF 4.2). A chain stuck in its
  own mode looks well mixed from the inside, which is why the problem went unnoticed.
- **Cause.** The model ran 7 dropping-level + 5 PCR-level latent factors. Season is constant
  within a dropping, so those factors can absorb the Season contrast: they correlate with Season
  at |r| ≈ 0.57 in 2 of 4 chains. Each chain split the Season signal differently between β and
  the factors. Confirmed by arm A5: restoring 5 factors brings the leak back (|r| ≈ 0.45), and
  Season VP drops to 10%.
- **Longer chains would not have fixed it.** The problem is structural.

## 2. What was tried

All arms use the same data and guild trait as `10_hmsc.R`. Pilots ran 5k iterations,
production 75k (thin 50), 4 chains each.

| Arm | Change | Result |
|---|---|---|
| **A** | Same CLR model; factors capped at dropping 2 + PCR 1; `YScale` | Season β converge (median PSRF 1.03, Fig 5 OTUs ESS > 3,000, guild means stable). **But the dropping-vs-PCR split has two modes** (chains 1,4: .52/.14; chains 2,3: .43/.27), identical at 5k and 75k iterations, so thin 100 won't help. Fails 2 of 7 criteria (max PSRF 1.41, VP). |
| A0 | A without the PCR level | Near-pass at pilot length |
| A5 | A with 5 dropping factors | Season leaks back into the factors, β shrink. Confirms the cause. |
| **C** | CLR averaged over each dropping's 2 PCR replicates; dropping level only (nf 2) | **Passes every criterion.** Beta PSRF max 1.005, 100% < 1.1; chains agree within ±0.02 on every Fig 5 OTU, guild mean and VP component; Gamma converges (max 1.002). |
| B | Hurdle abundance: rCLR, zeros → NA (abundance where present), 85 OTUs | Fits the data's shape (see §4), but still doesn't fully converge (max PSRF 2.26). Mainly useful as a finding, see §4. |

**A and C agree OTU by OTU:** Season β Pearson 0.98, Spearman 0.92. A's CrIs are about 13%
narrower, which is expected: within-dropping residual correlation between replicates ≈ 0.6
remains, so replicates are partly counted twice. Resolved sets overlap on 64 of 72–74 OTUs.

## 3. What changes in the results (C shown; A in brackets)

| Quantity | Published (non-converged) | New |
|---|---|---|
| OTUs with 95% CrI excl. 0 | 4 / 226 | **74** (72) |
| Fig 5A OTUs at 0.95 support | 9 / 22, all winter-end | **16** / 22 (18), coprophiles now resolved (support 1.00) |
| *Thelebolus* / *Sporormiella* / *Coniochaeta* β (CLR units) | 1.01 / 1.10 / 0.97 | **4.18 / 3.52 / 3.43** (3.66 / 3.34 / 2.78) |
| OTU1369 β | −1.68 | **−4.97** (−4.91) |
| Guild means: dung sap / other / plant / dark / patho | +0.39 / +0.04 / −0.01 / −0.15 / −0.17 | **+1.32 / +0.22 / −0.10 / −0.31 / −0.33** |
| VP community: Season / Year / dropping / PCR | 11 / 8 / 49 / 31 % | **31.5 / 18.5 / 50.0 / —** (19 / 12 / 47 / 21, bimodal) |
| Guild → Season Gamma (abundance) | not quotable (PSRF 1.46) | **quotable.** Dung sap +0.31 (support 0.99); dark, pathotroph, plant significantly below dung sap |
| Spearman vs probit / GLLVM | 0.82 / 0.70 | 0.76 / 0.75 (0.82 / 0.72) |

**Claims that survive:** dung saprotrophs are the most summer-shifted guild; dark taxa and
pathotrophs lean winter; OTU1369 is the most winter-shifted OTU; the ranking agrees with the
probit and the GLLVM.

**Claims that must change:**
- "Dung saprotroph is the ONLY positive guild" is false; "other" is +0.22.
- "Every resolved OTU is at the winter end" and "a large but individually uncertain coprophile
  effect" are false: the coprophiles are now the best-resolved OTUs.
- "Season 11% ≈ PERMANOVA R² 0.11" must go. The 11% was an artefact of the absorption, and HMSC
  VP partitions explained variance, not total variance, so it was never comparable to R².
- The panel-A axis needs roughly ±6, not ±2.5.

## 4. Caveats to carry, whichever model is chosen

- **Distribution misfit.** The Gaussian can't reproduce the skew of zero-heavy CLR data.
  Posterior-predictive skew is extreme for about 90% of OTUs, in the published model too. Mean
  effects are robust to this; the CrIs are approximate.
- **The Season signal is mostly about presence, not abundance (arm B).** Given an OTU is
  present, only 7 of 85 OTUs differ by season, and those estimates barely track the CLR or probit
  (ρ 0.37 / 0.10). *Thelebolus* (+1.9) and Sporormiaceae OTU1320 (+2.6) do bloom more in
  summer where present. *Sporormiella*'s summer signal is purely about how often it turns up.
  This is worth one sentence in §9.4. It fits the existing "coprophiles are what summer dung
  is made of" framing, but sharpens it.
- **Year effects.** Year takes 37% of variance in the where-present abundance model, which
  could be a sequencing-run effect. Worth a look before the Year terms are interpreted.
- **Units gotcha.** With `YScale = TRUE`, this Hmsc version stores β on the per-OTU
  **standardised** scale. `02_compare.R` multiplies by `YScalePar[2, ]` to return to CLR units.

## 5. Proposed edits to main (only after approval)

- `Scripts/10_hmsc.R`: cap factors with `setPriors(rL, nfMax=2, nfMin=2)` (plus PCR `nfMax=1`
  if A is chosen); `YScale = TRUE`; if C, average CLR over replicates before `Hmsc()`.
  Back-transform β in §6.3 (the `YScalePar` gotcha); add the per-chain checks from
  `00_diagnostics.R` to §5.
- **The new model needs its own fresh fits.** The old caches must be renamed (not deleted) or
  `fit_or_load` silently reloads them. The CLR CV (§6.5, multi-hour) and the phylo variant were
  built on the old structure and need re-running, or dropping, if they are to be quoted.
- **Appendix:**
  - §9.1: model description, Table S-HMSCb, and the false "well behaved" sentence.
  - §9.2: Table S-HMSCc, and drop the PERMANOVA-match claim.
  - §9.3: the abundance Gamma can now be quoted next to the probit.
  - §9.4: Beta tables, the Figure S-HMSCd ρ, and the occurrence-vs-abundance paragraph.
  - Figure 5 caption and `fig5-build` guards: XCLIP, VP_EXP, GM_EXP, n_ci, the psrf guard, and
    the PCR bar if C.
- The new §6 sentence on main ("OTU-level seasonal pattern is carried by HMSC in Section 9")
  becomes true with either A or C.
