#!/usr/bin/env bash
# =============================================================================
# repro/run_repro.sh -- regenerate every number in the Supplementary Appendix
# and check it against what is committed.
#
#   ./repro/run_repro.sh cached          # scripts 4-10 + appendix chunks, fitted
#                                        # models reused (copies) -- a few hours
#   ./repro/run_repro.sh cold            # same, every model refitted from
#                                        # scratch -- ~15-20 h, run in background
#   ./repro/run_repro.sh checks <RUN_ID> # re-run only the checks for a run
#
# Where things go:
#   * Every numbered script writes to a SANDBOX outside the repo
#     ($REPRO_SANDBOX/<RUN_ID>/{models,plots,tables}); the canonical
#     ~/Ptarmigan/models|plots|tables and eco_analysis.RData are only READ.
#   * Everything the appendix reads is staged into THIS checkout's
#     Supplementary/, so `git diff` there is the byte-level drift detector.
#     For that reason the run refuses to start unless Supplementary/ and
#     Scripts/ are clean.
#   * Logs go to repro/logs/<RUN_ID>/ (small text, committed).
#
# Scripts 1-3 (demultiplexing, DADA2/LULU, taxonomy) are not re-run: their
# outputs are pinned by checksum in inputs.sha256.
# =============================================================================
set -uo pipefail

MODE="${1:-cached}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RSCRIPT=(conda run --no-capture-output -n r_env Rscript)
CANON=/home/daniel/Ptarmigan
DATA="$CANON/trimmed/mergedPlates"
SANDBOX_ROOT="${REPRO_SANDBOX:-$CANON/repro_sandbox}"

# -- checks-only mode ---------------------------------------------------------
if [[ "$MODE" == "checks" ]]; then
  RUN_ID="${2:?usage: run_repro.sh checks <RUN_ID>}"
  LOG="$REPO/repro/logs/$RUN_ID"; SANDBOX="$SANDBOX_ROOT/$RUN_ID"
  [[ -d "$LOG" ]] || { echo "no such run: $LOG"; exit 2; }
  "${RSCRIPT[@]}" "$REPO/repro/run_checks.R" "$REPO" "$LOG" "$SANDBOX"
  exit $?
fi
[[ "$MODE" == "cached" || "$MODE" == "cold" ]] || { echo "mode must be cached|cold|checks"; exit 2; }

# -- refuse to start on a dirty tree (git diff must mean "this run changed it")
dirty=$(git -C "$REPO" status --porcelain -- Supplementary Scripts hmsc_converge)
if [[ -n "$dirty" && "${REPRO_ALLOW_DIRTY:-0}" != "1" ]]; then
  echo "Supplementary/ or Scripts/ has uncommitted changes -- commit or stash first:"
  echo "$dirty"; exit 2
fi

RUN_ID="$(date +%Y-%m-%d_%H%M)_$MODE"
LOG="$REPO/repro/logs/$RUN_ID"
SANDBOX="$SANDBOX_ROOT/$RUN_ID"
mkdir -p "$LOG/steps" "$SANDBOX"/{models,plots,tables,frozen,workspace}
STEPS="$LOG/steps.tsv"
printf "step\tstatus\texit\tminutes\tlog\n" > "$STEPS"
date +%s > "$LOG/.run_start_epoch"
echo "$MODE" > "$LOG/.mode"
echo "RUN $RUN_ID | repo $REPO | sandbox $SANDBOX"

# Threads: HMSC oversubscribes BLAS on this 96-core box unless capped.
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 HMSC_BLAS_THREADS=4 HMSC_NPARALLEL=4

# step <name> <cmd...> : run, time, log; record PASS/FAIL. Never aborts the run
# (except where a gate below says so) so one failure doesn't hide the rest.
step() {
  local name="$1"; shift
  local f="$LOG/steps/$name.log" t0=$(date +%s)
  echo "[$(date +%H:%M)] $name ..."
  { echo "# $name"; echo "# cmd: $*"; echo "# started $(date -Is)"; } > "$f"
  "$@" >> "$f" 2>&1
  local rc=$?
  local mins=$(( ($(date +%s) - t0) / 60 ))
  local st=PASS; [[ $rc -ne 0 ]] && st=FAIL
  printf "%s\t%s\t%d\t%d\tsteps/%s.log\n" "$name" "$st" "$rc" "$mins" "$name" >> "$STEPS"
  echo "[$(date +%H:%M)] $name -> $st ($mins min)"
  return $rc
}

# -- 0. environment + upstream checksums --------------------------------------
step 00_env "${RSCRIPT[@]}" "$REPO/repro/env_info.R" "$REPO" "$LOG/env.txt"
step 00_checksums bash -c "cd / && sha256sum \
  '$DATA/merged_ITS.RData' '$DATA/eco_analysis.RData' '$DATA/tree.rds' \
  '$DATA/nopoolps.dada2.rds' '$DATA/fullps_nopool.rds' \
  '$DATA/matches_out_taxonomy_nopool.csv' '$DATA/funguild_input.guilds.txt' \
  '$CANON/plant_ITS/phyloseq_plant_ITS.rds' \
  '$CANON/sh_matching/outdata/source_1.zip' '$CANON/models/source_67536.zip' \
  '$CANON/models/fit_nb_2.rds' '$CANON/models/fit_nb_2_2.rds' \
  '$CANON/models/hmsc_v2/hmsc_clr_fit_diag.rds' \
  '$CANON/reference_db/sh_general_release_dynamic_s_19.02.2025.fasta' \
  > '$LOG/inputs.sha256'"

# Freeze the scripts: Rscript reads a file incrementally, so an edit made while
# a long run is going would corrupt its tail. Run from copies.
cp "$REPO"/Scripts/*.R "$SANDBOX/frozen/"
S="$SANDBOX/frozen"

# -- cached mode: seed the sandbox with COPIES of the fitted models -----------
# (cp, never hardlinks: saveRDS truncates in place and would write through.)
if [[ "$MODE" == "cached" ]]; then
  while read -r f; do
    [[ -z "$f" || "$f" == \#* ]] && continue
    cp -p "$CANON/models/$f" "$SANDBOX/models/" || echo "WARN: cached fit missing: $f" >> "$LOG/steps/seed.log"
  done < "$REPO/repro/cached_fits.txt"
fi

# -- 1. script 4 -> sandbox workspace, then GATE on object identity -----------
step 01_script4 env S4_OUT_DIR="$SANDBOX/workspace" "${RSCRIPT[@]}" "$S/4_data_prep.R"
if ! step 01_compare_workspace "${RSCRIPT[@]}" "$REPO/repro/compare_workspace.R" \
       "$DATA/eco_analysis.RData" "$SANDBOX/workspace/eco_analysis.RData" "$LOG/workspace_compare.csv"; then
  echo "GATE: regenerated workspace differs from the canonical one -- downstream scripts"
  echo "      read the canonical workspace, so their results would not be a valid check. Stopping."
  "${RSCRIPT[@]}" "$REPO/repro/run_checks.R" "$REPO" "$LOG" "$SANDBOX"
  exit 1
fi

# -- 2-4. analysis scripts ------------------------------------------------------
SUPP="$REPO/Supplementary"
step 05_script5 env S5_OUT_ROOT="$SANDBOX" S5_SUPP_DIR="$SUPP" "${RSCRIPT[@]}" "$S/5_community_composition.R"
step 08_script8 env S8_OUT_ROOT="$SANDBOX" S8_SUPP_DIR="$SUPP" "${RSCRIPT[@]}" "$S/8_functional_guilds.R"
step 09_script9 env S9_OUT_ROOT="$SANDBOX" S9_SUPP_DIR="$SUPP" "${RSCRIPT[@]}" "$S/9_dark_taxa_SH_matching.R"
step 06_script6 env S6_OUT_ROOT="$SANDBOX" S6_SUPP_DIR="$SUPP" "${RSCRIPT[@]}" "$S/6_diversity_analyses.R"
step 07_script7 env S7_OUT_ROOT="$SANDBOX" S7_SUPP_DIR="$SUPP" "${RSCRIPT[@]}" "$S/7_gllvm.R"
if [[ "$MODE" == "cold" ]]; then
  step 10_script10 env HMSC_OUT_ROOT="$SANDBOX" HMSC_SUPP_DIR="$SUPP" HMSC_RUN_PHYLO=1 \
       HMSC_RUN_CV=1 HMSC_NO_REUSE=1 "${RSCRIPT[@]}" "$S/10_hmsc.R"
else
  step 10_script10 env HMSC_OUT_ROOT="$SANDBOX" HMSC_SUPP_DIR="$SUPP" HMSC_RUN_PHYLO=1 \
       HMSC_RUN_CV=0 "${RSCRIPT[@]}" "$S/10_hmsc.R"
fi

# -- 5. appendix build chunks (source of truth for Sections 2/4/5 + main figs) --
for lab in sec2-build sec3-build sec4-build sec5-build fig1-build fig2-build fig3-build fig4-build fig5-build; do
  step "20_${lab}" bash -c "cd '$SUPP' && ${RSCRIPT[*]} build_main_figure.R $lab"
done
step 21_hmsc_v1_table bash -c "cd '$REPO' && ${RSCRIPT[*]} hmsc_converge/04_v1_disagreement_table.R"

# -- 6. render both profiles -------------------------------------------------------
step 30_render_submission bash -c "cd '$SUPP' && conda run --no-capture-output -n r_env quarto render Supplementary_Appendix.qmd --to html"
step 31_render_docs bash -c "cd '$SUPP' && conda run --no-capture-output -n r_env quarto render Supplementary_Appendix.qmd --to html --profile docs --output Appendix_documentation.html"
step 32_fence_check bash -c "cd '$SUPP' && ! grep -o '[A-Za-z0-9_.]*:::[^ <\"]\{0,20\}' Appendix_HTML.html Appendix_documentation.html | grep -v 'vegan:::\.calc_rclr'"

step 33_stub_visibility python3 "$REPO/repro/stub_visibility.py"

# -- 7. checks + summary ---------------------------------------------------------------
"${RSCRIPT[@]}" "$REPO/repro/run_checks.R" "$REPO" "$LOG" "$SANDBOX"
rc=$?
echo "Summary: $LOG/summary.md"
exit $rc
