# Reproducibility audit — findings (2026-10-02)

Branch `exp/repro-audit` (worktree `.worktrees/repro-audit`). Nothing on `main` was changed.
`main`'s Table S4 / §7.1 commits from today were merged into the branch before this final run.
The machine-generated `summary.md` in this folder is the full record; this file is the
human-curated reading of it and of the two earlier runs (`2026-10-02_1032_cached`,
`2026-10-02_1147_cached`).

## Bottom line

From the committed code plus the saved model fits:

- **Analysis objects.** A fresh `4_data_prep.R` rebuilds `alldat`, `alldat.rfy`, `alldat_full` and
  `alldat_rg2_sensitivity` exactly.
- **Tables.** All **70 of 70** CSV tables the appendix reads come back byte-identical.
- **Main-text figures.** All five are byte-identical as PNG and JPG; the PDFs render pixel-identically (only the embedded date differs).
- **Supplementary figures.** **31 of 35** supplementary PNGs are byte-identical. The other 4 are the stale §3.1 barplots (decision 1 below).
- **Quoted numbers.** **422 of 422** registered numbers in the prose and captions match their tables.
- **Producers.** Every file under `Supplementary/tables|figures` (116) is regenerated during the run by the producer `repro/manifest.csv` names.
- **Pointer chunks.** Each script-backed section has one, and each renders in both builds.
- **Runtime.** All scripts 4–10, every build chunk and both renders run without error. Cached run: 52 min, 45 of them script 6's phylogenetic Part B.

## Decisions for Daniel

1. **§3.1 Figures S2–S5 are stale.** They were made on 2026-07-03 from an older `alldat_full`.
   - Script 5 now builds a **55**-sample object, not the 57 the appendix text states (§3 and §10 documentation notes).
   - The regenerated barplots differ visibly: taxonomy labels and sample composition both changed. Compare with `git diff` after a run.
   - Options:
     - (a) accept the regenerated figures and change "57" to "55";
     - (b) switch §3.1 to the canonical 52-sample object, as Fig. 1A already uses;
     - (c) keep them as they are, documented as stale.
   - Until you decide, the committed figures are left untouched and every run reports them.
2. **Overnight cold refit** (`./repro/run_repro.sh cold`, ~15–20 h).
   - It checks that the fits themselves reproduce, not just everything downstream of them.
   - MCMC refits will not be bit-identical, so they are judged at a 0.05 relative tolerance and then by eye against Monte-Carlo error.
3. **Numbers quoted but stored in no table** (section C). Do you want them written to tables, so the registry can cover them too?
4. **Merge** `exp/repro-audit` (`git merge --no-ff`), and apply `repro/CLAUDE_md_proposal.md`.

## A. Fixed on the branch (all verified by the final run)

| Where | Problem | Fix |
|---|---|---|
| `4_data_prep.R` | trailing bare `funguild_otu` / `fungaltraits_otu` lines errored under `Rscript` after the save (non-zero exit) | commented out |
| scripts 4/5/8/9 | wrote only to the canonical paths, so they could not be run in isolation; script 8 could rewrite `eco_analysis.RData` | `S4_OUT_DIR`, `S5/S8/S9_OUT_ROOT` + `_SUPP_DIR` (defaults unchanged) |
| `10_hmsc.R` | no way to force a cold refit | `HMSC_NO_REUSE=1` |
| scripts 5, 7 | §3.1 barplots and all §6 outputs were copied into `Supplementary/` by hand | the scripts now stage them |
| §5 | Tables S8/S9 + Figures S7/S8 had two producers (`sec5-build` and script 6, byte-identical) | script 6 is the single producer; `sec5-build` keeps Tables S10/S11/S3 and reads script 6's per-sample Hill table; proven 13/13 byte-identical |
| `8_functional_guilds.R` | crosswalk note shown in Table S-GUILDa said "Section 10 'other'" (HMSC is §9 since 2026-08-18); all its section refs used the old numbering | updated |
| Figure S-GUILDa | unseeded jitter (figure changed on every run) | `seed = 1` |
| Figure S-GLLVMc | unseeded bootstrap CIs | `set.seed(1)` |
| GLLVM residual diagnostics | randomised quantile residuals unseeded (not shown) | `set.seed(1)` |
| Figure S-GLLVMb | committed PNG drawn by an older plotting version (overlapping text labels); same coordinates | regenerated: points, as the caption describes |
| `hill_size_helper.R`, Fig. 4 provenance | claimed to be the "single source" for §5.2 and called Fig. 4 "Figure 3" | reworded: fig4-build sources it, script 6 carries the copy, the build guard checks they agree |
| `build_main_figure.R` | "Section 11" | Section 10 |
| §8 | no "Documentation notes — Section 8" bin (every other section has one) | added (empty) |

Correction to the plan: the suspected wrong filename "`2_lulu.R`" was a false alarm. The
appendix cites `2_DADA2_lulu.R` correctly.

## B. Informational

- **Script 6 Part B** (phylogenetic, not in the appendix) runs cleanly, but its five
  `PhyloSupp_*` tables differ from the copies in canonical `models/`. Those copies predate
  the July depth fix (CLAUDE.md already notes Part B was never re-run). Re-run Part B, or
  delete it if the phylogenetic analyses are definitively out.
- **Table S5b vs S7c.** The rank-3 row of Table S7c repeats the Table S5b test with
  different permutation p-values (0.775/0.009/0.116 vs 0.789/0.005/0.127). The §4.3 prose
  already explains this (separately seeded streams).
- **Stale documentation-only prose** (yours to rewrite, not touched):
  - the closing status callout (end of §10) still uses the pre-2026-08-18 numbering ("9 (functional guilds), 10 (HMSC) and 11 (main-text figures)");
  - the overview paragraph still lists "phylogenetic community-structure null models" among the remote computations.

## C. Numbers quoted but stored in no table (triage of `unverified_numbers.csv`)

124 numbers in the prose are tied to no table. Most are legitimate: settings, seeds,
literature values, and history ("ρ ≈ 0.72 reported from the earlier model"). These are
**computed by the pipeline but only printed**, so a future change could silently
invalidate them:

| Section | Numbers | Computed in |
|---|---|---|
| 3.2 | medians 47.5% / 0.02%, Wilcoxon p = 0.0012 | `sec5-build` (printed) |
| 3.4 | 34.6% placeable winter dark burden | script 9 writes `E1_winter_dark_placeable_fraction.csv` to `models/` but does not stage it |
| 3.4 | 21 of 79 to species, 11 of 19 gain a genus; 554 vs 553 `present_in`; 248 PlutoF chimera flags; median ~97% | script 9 logs |
| 4.1 / 10.2 | PCoA axis shares 58.0 / 31.3 / 10.7; full-rank 13.7 / 9.5 / 23.2; Season per axis 15.3 / 4.0 / 10.3 → 8.9 / 1.2 / 1.1 | `sec4-build` (printed) |
| 4.4 | 171 / 265 / 693 shared / winter-only / summer-only; 9 / 51 / 11 at ≥ 5 detections; Sørensen R² 0.113 | `sec4-build` |
| 5 | coverage ≥ 0.99999 in 52, exactly 1 in 46, 43 without singletons; library medians 551,236 vs 209,602 (2.6×); Pielou vs size-std evenness r = 0.995 | script 6 / `sec5-build` |
| 6.3 | winter LV2 variance 3.4×, Fligner p = 0.013; mean distance to centroid 24.8 vs 13.6 | script 7 / `sec4-build` |
| 8.3 | within-year q = 0.027 / 0.053 / 0.74; FungalTraits 195 OTUs, 17%, 37%, 49% | script 8 (printed) |
| 9.2 / 9.4 | r = 0.998; depth ratio 1.31; ρ = 0.21, p = 0.13; probit → Gamma within 0.04 | script 10 |
| 10.3 | 22 OTUs, 73% of modelled reads, 16 resolved | `fig3-build` |
| 10.5 | Wilcoxon p = 1.9 × 10⁻⁵ | script 9 (printed) |

Writing these to small CSVs would let the registry cover them. It needs a few lines in
each producer and touches no result.
