# Reproducibility run `2026-10-02_1250_cached`

**Verdict: FINDINGS -- see below**

- mode: `cached` (fitted models reused as copies; everything downstream of a fit regenerated)
- git commit  : 55ce64f336be5ad3993bd3ebfb07814b3489fb50   
- run date    : 2026-10-02 12:50:51 CEST   
- sandbox: `/home/daniel/Ptarmigan/repro_sandbox/2026-10-02_1250_cached`

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
| 06_script6 | PASS | 0 | 45 | steps/06_script6.log |
| 07_script7 | PASS | 0 | 0 | steps/07_script7.log |
| 10_script10 | PASS | 0 | 6 | steps/10_script10.log |
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
| sandbox_vs_canonical | DRIFT | 6 |
| appendix | identical | 106 |
| sandbox_vs_canonical | identical | 62 |
| script_vs_appendix_copy | identical | 56 |
| appendix | same_pixels | 6 |
| appendix | VISUAL | 4 |

Appendix files (and script copies of them) not identical or equal in value / pixels:

| set | file | regenerated | status | detail |
|---|---|---|---|---|
| appendix | Supplementary/figures/Order_barplot_ordered_season.png | TRUE | VISUAL | 14.9754% of pixels differ, max channel diff 0.953 |
| appendix | Supplementary/figures/class_barplot_ordered_season.png | TRUE | VISUAL | 32.8344% of pixels differ, max channel diff 1.000 |
| appendix | Supplementary/figures/genus_barplot_ordered_season.png | TRUE | VISUAL | 57.0972% of pixels differ, max channel diff 1.000 |
| appendix | Supplementary/figures/phylum_barplot_ordered_season.png | TRUE | VISUAL | 13.7528% of pixels differ, max channel diff 0.933 |

For information -- script outputs NOT shown in the appendix that differ from the canonical models/ copy (does not affect the verdict):

| set | file | regenerated | status | detail |
|---|---|---|---|---|
| sandbox_vs_canonical | guild_scheme_crosswalk.csv | TRUE | DRIFT | max rel diff 0 (col -); 2 text cells differ (e.g. note) |
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

Registry: 422 of 422 claims pass.

_none_

Sweep of all other numbers in prose and captions (788): convention 77, design_constant 118, registry 469, UNVERIFIED 124. Only the registry verifies a number. UNVERIFIED ones are listed in `unverified_numbers.csv`: quoted in the prose but tied to no table -- each is either worth adding to the registry (and, if no table stores it yet, to a table) or legitimately table-free (literature values, settings, dates).

