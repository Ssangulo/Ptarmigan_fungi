# repro/run_checks.R -- run the three checks for one reproducibility run and
# write <log>/summary.md.  Called at the end of run_repro.sh, or on its own via
# `run_repro.sh checks <RUN_ID>`.
# Usage: Rscript run_checks.R <repo> <log_dir> <sandbox_dir>
args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(args[1]); log <- args[2]; sandbox <- args[3]
here <- file.path(repo, "repro")
for (f in c("compare_outputs.R", "check_manifest.R", "check_claims.R")) source(file.path(here, f))

mode      <- readLines(file.path(log, ".mode"))
run_start <- as.POSIXct(as.numeric(readLines(file.path(log, ".run_start_epoch"))), origin = "1970-01-01")
tol       <- if (mode == "cold") 0.05 else 1e-8     # cold: flag, then judge against MC error by hand

steps <- utils::read.delim(file.path(log, "steps.tsv"), stringsAsFactors = FALSE)
cmp   <- compare_outputs(repo, sandbox, run_start, tol)
utils::write.csv(cmp, file.path(log, "diff_outputs.csv"), row.names = FALSE)
man   <- check_manifest(repo, cmp)
utils::write.csv(man, file.path(log, "manifest_check.csv"), row.names = FALSE)
reg   <- check_registry(repo)
utils::write.csv(reg, file.path(log, "claims_result.csv"), row.names = FALSE)
sw    <- sweep_numbers(repo, utils::read.csv(file.path(here, "claims.csv"), stringsAsFactors = FALSE))
utils::write.csv(sw[sw$category == "UNMATCHED", ], file.path(log, "unmatched_numbers.csv"), row.names = FALSE)
utils::write.csv(sw, file.path(log, "sweep_all_numbers.csv"), row.names = FALSE)
wsc <- file.path(log, "workspace_compare.csv")
ws  <- if (file.exists(wsc)) utils::read.csv(wsc, stringsAsFactors = FALSE) else NULL

# ---- summary.md ---------------------------------------------------------------
md <- character(0); w <- function(...) md <<- c(md, paste0(...))
tbl <- function(d) {
  if (is.null(d) || !nrow(d)) return(w("_none_\n"))
  d[] <- lapply(d, function(x) gsub("\\|", "/", as.character(x)))
  w("| ", paste(names(d), collapse = " | "), " |")
  w("|", paste(rep("---", ncol(d)), collapse = "|"), "|")
  for (i in seq_len(nrow(d))) w("| ", paste(d[i, ], collapse = " | "), " |")
  w("")
}
env <- readLines(file.path(log, "env.txt"), warn = FALSE)
ok_out  <- c("identical", "same_values", "same_pixels", "same_content", "within_tol", "NO_CANONICAL")
bad_out <- cmp[!cmp$status %in% ok_out, ]
app_not_regen <- man[man$check == "regenerated" & man$status == "FAIL", ]
fail_steps <- steps[steps$status != "PASS", ]
fail_man   <- man[man$status == "FAIL", ]
fail_reg   <- reg[reg$status == "FAIL", ]
gate_ok    <- !is.null(ws) && all(ws$status[ws$role == "GATE"] == "PASS")

verdict <- if (!nrow(fail_steps) && !nrow(bad_out) && !nrow(fail_man) && !nrow(fail_reg) && gate_ok)
  "PASS" else "FINDINGS -- see below"
w("# Reproducibility run `", basename(log), "`\n")
w("**Verdict: ", verdict, "**\n")
w("- mode: `", mode, "` (", if (mode == "cached") "fitted models reused as copies; everything downstream of a fit regenerated" else
  "every model refitted from scratch", ")")
w("- ", grep("^git commit", env, value = TRUE), "  ")
w("- ", grep("^run date", env, value = TRUE), "  ")
w("- sandbox: `", sandbox, "`\n")
w("Scripts 1-3 are not re-run; their outputs are pinned in `inputs.sha256`.\n")

w("## 1. Steps\n"); tbl(steps)
w("## 2. Script 4 workspace gate\n")
if (is.null(ws)) w("_not run_\n") else {
  w(if (gate_ok) "PASS -- a fresh `4_data_prep.R` rebuilds `alldat`, `alldat.rfy`, `alldat_full` and `alldat_rg2_sensitivity` exactly.\n"
    else "**FAIL** -- see `workspace_compare.csv`.\n")
  tbl(ws[ws$role != "other" | !ws$status %in% c("PASS"), ][, c("object", "role", "status", "detail")])
}
w("## 3. Outputs vs committed / canonical\n")
w("Counts by comparison set and status:\n")
tbl(as.data.frame(table(set = cmp$set, status = cmp$status), stringsAsFactors = FALSE)[
  as.data.frame(table(cmp$set, cmp$status))$Freq > 0, ])
w("Everything not identical (or equal in value / pixels):\n"); tbl(bad_out)
w("## 4. Manifest\n")
tbl(as.data.frame(table(check = man$check, status = man$status), stringsAsFactors = FALSE)[
  as.data.frame(table(man$check, man$status))$Freq > 0, ])
w("Failures and exceptions:\n"); tbl(man[man$status %in% c("FAIL", "EXCEPTION"), ])
w("## 5. Quoted numbers\n")
w(sprintf("Registry: %d of %d claims pass.\n", sum(reg$status == "PASS"), nrow(reg)))
tbl(fail_reg)
w(sprintf("Sweep of all other numbers in prose and captions (%d): ", nrow(sw)),
  paste(sprintf("%s %d", names(table(sw$category)), table(sw$category)), collapse = ", "),
  ". Unmatched ones are listed in `unmatched_numbers.csv` for triage (sample sizes, thresholds and literature values legitimately have no table).\n")
writeLines(md, file.path(log, "summary.md"))
cat(md, sep = "\n")
quit(status = if (verdict == "PASS") 0 else 1)
