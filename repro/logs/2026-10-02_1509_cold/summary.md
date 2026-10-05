# Reproducibility run `2026-10-02_1509_cold`

**Verdict: FINDINGS -- see below**

- mode: `cold` (every model refitted from scratch)
- git commit  : ef5c00d84acb79b491573fb1eadde6753e3ddfbd   
- run date    : 2026-10-02 15:09:55 CEST   
- sandbox: `/home/daniel/Ptarmigan/repro_sandbox/2026-10-02_1509_cold`

Scripts 1-3 are not re-run; their outputs are pinned in `inputs.sha256`.

## 1. Steps

| step | status | exit | minutes | log |
|---|---|---|---|---|
| 00_env | PASS | 0 | 0 | steps/00_env.log |
| 00_checksums | PASS | 0 | 0 | steps/00_checksums.log |
| 01_script4 | PASS | 0 | 0 | steps/01_script4.log |
| 01_compare_workspace | PASS | 0 | 0 | steps/01_compare_workspace.log |
| 05_script5 | PASS | 0 | 0 | steps/05_script5.log |
| 08_script8 | PASS | 0 | 0 | steps/08_script8.log |
| 09_script9 | PASS | 0 | 0 | steps/09_script9.log |
| 06_script6 | PASS | 0 | 111 | steps/06_script6.log |
| 07_script7 | PASS | 0 | 0 | steps/07_script7.log |
| 10_script10 | PASS | 0 | 689 | steps/10_script10.log |
| 20_sec2-build | PASS | 0 | 1 | steps/20_sec2-build.log |
| 20_sec3-build | PASS | 0 | 0 | steps/20_sec3-build.log |
| 20_sec4-build | PASS | 0 | 0 | steps/20_sec4-build.log |
| 20_sec5-build | PASS | 0 | 0 | steps/20_sec5-build.log |
| 20_fig1-build | PASS | 0 | 0 | steps/20_fig1-build.log |
| 20_fig2-build | PASS | 0 | 0 | steps/20_fig2-build.log |
| 20_fig3-build | PASS | 0 | 0 | steps/20_fig3-build.log |
| 20_fig4-build | PASS | 0 | 0 | steps/20_fig4-build.log |
| 20_fig5-build | PASS | 0 | 0 | steps/20_fig5-build.log |
| 21_hmsc_v1_table | PASS | 0 | 0 | steps/21_hmsc_v1_table.log |
| 30_render_submission | PASS | 0 | 0 | steps/30_render_submission.log |
| 31_render_docs | PASS | 0 | 0 | steps/31_render_docs.log |
| 32_fence_check | PASS | 0 | 0 | steps/32_fence_check.log |
| 33_stub_visibility | PASS | 0 | 0 | steps/33_stub_visibility.log |

## 2. Script 4 workspace gate

PASS -- a fresh `4_data_prep.R` rebuilds `alldat`, `alldat.rfy`, `alldat_full` and `alldat_rg2_sensitivity` exactly.

| object | role | status | detail |
|---|---|---|---|
| alldat | GATE | PASS | identical |
| alldat_full | GATE | PASS | identical |
| alldat_rg2_sensitivity | GATE | PASS | identical |
| alldat.rfy | GATE | PASS | identical |
| fungaltraits_otu | added_by_script8 | EXPECTED | only in canonical |
| funguild_otu | added_by_script8 | EXPECTED | only in canonical |
| .Random.seed | other | INFO | only in regenerated |
| unresolved_export_report | other | INFO | Component “fasta_file”: 3 string mismatches / Component “table_file”: 3 string mismatches |

## 3. Outputs vs committed / canonical

Counts by comparison set and status:

| set | status | Freq |
|---|---|---|
| appendix | DRIFT | 5 |
| sandbox_vs_canonical | DRIFT | 11 |
| appendix | identical | 98 |
| sandbox_vs_canonical | identical | 56 |
| script_vs_appendix_copy | identical | 56 |
| appendix | same_pixels | 6 |
| appendix | VISUAL | 6 |
| appendix | within_tol | 1 |
| sandbox_vs_canonical | within_tol | 1 |

Appendix files (and script copies of them) not identical or equal in value / pixels:

| set | file | regenerated | status | detail |
|---|---|---|---|---|
| appendix | Supplementary/figures/Order_barplot_ordered_season.png | TRUE | VISUAL | 14.9754% of pixels differ, max channel diff 0.953 |
| appendix | Supplementary/figures/class_barplot_ordered_season.png | TRUE | VISUAL | 32.8344% of pixels differ, max channel diff 1.000 |
| appendix | Supplementary/figures/genus_barplot_ordered_season.png | TRUE | VISUAL | 57.0972% of pixels differ, max channel diff 1.000 |
| appendix | Supplementary/figures/hmsc_R2_explanatory_vs_cv.png | TRUE | VISUAL | 7.5899% of pixels differ, max channel diff 1.000 |
| appendix | Supplementary/figures/hmsc_probit_gamma.png | TRUE | VISUAL | 11.9427% of pixels differ, max channel diff 0.020 |
| appendix | Supplementary/figures/phylum_barplot_ordered_season.png | TRUE | VISUAL | 13.7528% of pixels differ, max channel diff 0.933 |
| appendix | Supplementary/tables/hmsc_convergence_psrf.csv | TRUE | DRIFT | max rel diff 0.727 (col vp_chain_range_max); 0 text cells differ |
| appendix | Supplementary/tables/hmsc_predictive_R2.csv | TRUE | DRIFT | max rel diff 2 (col cv_R2); 0 text cells differ |
| appendix | Supplementary/tables/hmsc_predictive_R2_summary.csv | TRUE | DRIFT | max rel diff 0.15 (col mean_cv_R2); 0 text cells differ |
| appendix | Supplementary/tables/hmsc_probit_beta_season.csv | TRUE | DRIFT | max rel diff 2 (col beta_season_mean); 516 text cells differ (e.g. guild) |
| appendix | Supplementary/tables/hmsc_probit_gamma_CrI.csv | TRUE | DRIFT | max rel diff 1.6 (col CrI_97.5); 0 text cells differ |

For information -- script outputs NOT shown in the appendix that differ from the canonical models/ copy (does not affect the verdict):

| set | file | regenerated | status | detail |
|---|---|---|---|---|
| sandbox_vs_canonical | guild_scheme_crosswalk.csv | TRUE | DRIFT | max rel diff 0 (col -); 2 text cells differ (e.g. note) |
| sandbox_vs_canonical | hmsc_convergence_psrf.csv | TRUE | DRIFT | max rel diff 0.727 (col vp_chain_range_max); 0 text cells differ |
| sandbox_vs_canonical | hmsc_predictive_R2_summary.csv | TRUE | DRIFT | max rel diff 0.15 (col mean_cv_R2); 0 text cells differ |
| sandbox_vs_canonical | hmsc_predictive_R2.csv | TRUE | DRIFT | max rel diff 2 (col cv_R2); 0 text cells differ |
| sandbox_vs_canonical | hmsc_probit_beta_season.csv | TRUE | DRIFT | max rel diff 2 (col beta_season_mean); 516 text cells differ (e.g. guild) |
| sandbox_vs_canonical | hmsc_probit_gamma_CrI.csv | TRUE | DRIFT | max rel diff 1.6 (col CrI_97.5); 0 text cells differ |
| sandbox_vs_canonical | PhyloSupp_betaNTI_summary.csv | TRUE | DRIFT | max rel diff 0.868 (col median); 0 text cells differ |
| sandbox_vs_canonical | PhyloSupp_iNEXT3D_PD_asymptotic.csv | TRUE | DRIFT | max rel diff 0.259 (col s.e.); 0 text cells differ |
| sandbox_vs_canonical | PhyloSupp_NRI_NTI_season_summary.csv | TRUE | DRIFT | max rel diff 0.147 (col NTI_aw_IQR); 0 text cells differ |
| sandbox_vs_canonical | PhyloSupp_UniFrac_unweighted_permanova.csv | TRUE | DRIFT | max rel diff 0.98 (col Pr(>F)); 0 text cells differ |
| sandbox_vs_canonical | PhyloSupp_UniFrac_weighted_permanova.csv | TRUE | DRIFT | max rel diff 0.98 (col Pr(>F)); 0 text cells differ |

## 4. Manifest

| check | status | Freq |
|---|---|---|
| pointer_stub_coverage | PASS | 8 |
| pointer_stub_outputs | PASS | 11 |
| producer_names_file | PASS | 116 |
| regenerated | PASS | 116 |

Failures and exceptions:

_none_

## 5. Quoted numbers

Registry: 406 of 422 claims pass.

| id | section | file | quoted | value | value_ok | snippet_ok | status | detail |
|---|---|---|---|---|---|---|---|---|
| S9_pa_psrf_med | 9.1 | tables/hmsc_convergence_psrf.csv | 1.01 | 1.016 | FALSE | TRUE | FAIL | table gives 1.02 at 2 dp |
| S9_pa_psrf_pct | 9.1 | tables/hmsc_convergence_psrf.csv | 99.6 | 85.2 | FALSE | TRUE | FAIL | table gives 85.2 at 1 dp |
| S9_pa_chain_guild | 9.1 | tables/hmsc_convergence_psrf.csv | 0.24 | 0.6 | FALSE | TRUE | FAIL | table gives 0.60 at 2 dp |
| S9_pa_gamma_dung | 9.3 | tables/hmsc_probit_gamma_CrI.csv | 0.96 | 0.969 | FALSE | TRUE | FAIL | table gives 0.97 at 2 dp |
| S9_pa_gamma_patho | 9.3 | tables/hmsc_probit_gamma_CrI.csv | 0.85 | 0.892 | FALSE | TRUE | FAIL | table gives 0.89 at 2 dp |
| S9_pa_gamma_dung_s8 | 8.3 | tables/hmsc_probit_gamma_CrI.csv | 0.96 | 0.969 | FALSE | TRUE | FAIL | table gives 0.97 at 2 dp |
| S9_n_resolved_pa | 9.4 | tables/hmsc_probit_beta_season.csv | 55 | 58 | FALSE | TRUE | FAIL | table gives 58 at 0 dp |
| S9_pa_vs_clr_rho | 9.4 | tables/hmsc_beta_season.csv | 0.76 | 0.7892553 | FALSE | TRUE | FAIL | table gives 0.79 at 2 dp |
| S9_pa_n_pos | 9.4 | tables/hmsc_probit_beta_season.csv | 157 | 160 | FALSE | TRUE | FAIL | table gives 160 at 0 dp |
| S9_pa_med | 9.4 | tables/hmsc_probit_beta_season.csv | 0.70 | 0.8035 | FALSE | TRUE | FAIL | table gives 0.80 at 2 dp |
| S9_pa_dung_rank | 9.4 | tables/hmsc_probit_beta_season.csv | 125 | 120.5357 | FALSE | TRUE | FAIL | table gives 121 at 0 dp |
| S9_pa_otu1369 | 9.4 | tables/hmsc_probit_beta_season.csv | -2.25 | -2.439 | FALSE | TRUE | FAIL | table gives -2.44 at 2 dp |
| S9_pa_otu1369_lo | 9.4 | tables/hmsc_probit_beta_season.csv | -3.53 | -3.803 | FALSE | TRUE | FAIL | table gives -3.80 at 2 dp |
| S9_pa_otu1369_hi | 9.4 | tables/hmsc_probit_beta_season.csv | -1.08 | -1.148 | FALSE | TRUE | FAIL | table gives -1.15 at 2 dp |
| S9_r2_cv | 9.7 | tables/hmsc_predictive_R2_summary.csv | 0.03 | 0.04 | FALSE | TRUE | FAIL | table gives 0.04 at 2 dp |
| S10_3_r2cv | 10.3 | tables/hmsc_predictive_R2_summary.csv | 0.03 | 0.04 | FALSE | TRUE | FAIL | table gives 0.04 at 2 dp |

Sweep of all other numbers in prose and captions (788): convention 77, design_constant 118, registry 469, UNVERIFIED 124. Only the registry verifies a number. UNVERIFIED ones are listed in `unverified_numbers.csv`: quoted in the prose but tied to no table -- each is either worth adding to the registry (and, if no table stores it yet, to a table) or legitimately table-free (literature values, settings, dates).

