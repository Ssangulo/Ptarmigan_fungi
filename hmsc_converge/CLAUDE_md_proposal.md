# Proposed CLAUDE.md changes (apply only on approval, at merge time)

## 1. New DECISIONS entry (top of the list)

- 2026-09-29: **HMSC abundance model REPLACED — the published CLR fit had not converged**
  (branch `exp/hmsc-converge`; study + memo in `Scripts_server/hmsc_converge/`).
  The old CLR HMSC (replicate level, uncapped latent factors: 7 sample + 5 PCR) had four chains
  converged to four different answers — OTU1 Season β 0.36/3.03/0.69/−0.03, Season VP 5–28%
  — so Figure 5 and Table S-HMSCc quoted a pooled mixture (dung-sap mean +0.39, Season 11%).
  PSRF hid it (OTU1 ESS 2,157 at PSRF 4.2: each chain mixed well inside its own mode).
  **Cause:** Season is constant within a dropping, so enough sample-level latent factors
  reproduce the Season contrast; each chain split the signal differently between β and ΛH.
  Longer chains do not fix it (identical modes at 5k and 75k iterations).
  **New structure (`Scripts/10_hmsc.R`):** HEADLINE `m_clr` = CLR per replicate, then
  AVERAGED per dropping (55 rows), sample level only, `nfMax = 2`, `YScale = TRUE`;
  COMPLEMENT `m_rep` = same CLR at replicate level (109 rows), sample `nfMax 2` + PCR `nfMax 1`;
  probit unchanged (reused from cache with a data-identity guard). Unused `m_lib` removed.
  **Why averaged, not replicate-level, is the headline:** the replicate-level dropping/PCR split
  is weakly identified (54 pairs) — one of two independent fits (same spec, different seed) had
  chains in two modes (.52/.14 vs .43/.27) — and it counts replicates as partly independent
  (CrIs 0.84× the headline's). Headline passes every per-chain check (ranges ≤ 0.003).
  They agree per OTU: Pearson 0.97, Spearman 0.89.
  **Results that changed:** 74/226 OTUs resolved at 95% CrI (was 4); coprophiles now the best
  resolved (*Thelebolus* +4.17, *Sporormiella* +3.52, *Coniochaeta* +3.44 CLR units); OTU1369
  −5.00; VP Season 31.5 / Year 18.5 / dropping 49.9%; CLR Gamma now converges and is quotable
  (dung-sap Season +0.31 std units, support 0.994). "Only positive guild" is FALSE ("other"
  weakly +) and "Season 11% ≈ PERMANOVA R² 0.11" is DROPPED (artefact; VP is a share of
  explained variance, never comparable with R²). Spearman vs probit 0.76, vs GLLVM 0.75.
  **PCR-replicate variance is reported model-free** (§6.7 of script 10, Table S-HMSCc2, Fig 5B
  italic column): replicate variance = mean(d²)/2 over the 54 pairs ÷ total CLR variance →
  32% [boot 28–40]; lme4 check 34% (r = 0.998 per OTU); dung sap 17%, pathotroph 50%. It is
  almost all DETECTION dropout (46% of detections are in one replicate only); when both
  replicates detect an OTU, noise is a median 3% of variance.
  **Still pending:** phylo variant (§9.6, old structure — flagged in the appendix) and the
  probit's looser per-chain agreement (guild-mean chain range 0.24 on the probit scale).
  **GOTCHAs:** (1) with `YScale = TRUE`, Hmsc 3.3-7 stores Beta/Gamma on the per-OTU
  STANDARDISED scale — script 10 multiplies Season β by `YScalePar[2, ]` for CLR units.
  (2) `Rscript` reads a script INCREMENTALLY — editing a script while it runs corrupts the tail;
  launch long runs from a frozen copy. (3) `pgrep -f`/`pkill -f` with a pattern that appears in
  your own shell command matches the shell itself (exit 144) — match on `/proc/<pid>/environ` or
  a unique path instead. (4) Script 10 honours `HMSC_OUT_ROOT` / `HMSC_SUPP_DIR` so a worktree
  run cannot overwrite main's `models/`, `plots/`, `tables/`, `Supplementary/`.
  **Superseded:** the 2026-07-21 entry's "variance partition Season 11%/Year 8%/sample 49%/PCR
  31% (Season ~ PERMANOVA R2=0.11)" and "the CLR-Gaussian abundance Gamma does NOT converge";
  the 2026-09-11b Figure 5 entry's numbers (4 of 226, 9 of 22 all winter, guild means, ±2.5 clip,
  PSRF 1.17) — its selection rule and layout decisions still stand.

## 2. Edits to existing text

- CURRENT FOCUS, §10.5 sentence: "Figure 5 = HMSC synthesis (… 2 panels, both now the
  ABUNDANCE/CLR model …)" → add "(revised 2026-09-29: replicate-averaged converged model +
  model-free PCR column; see DECISIONS)".
- 2026-07-21 entry: append "SUPERSEDED in part 2026-09-29 — see that entry" after its Results.
- 2026-09-11b entry: append the same pointer.
