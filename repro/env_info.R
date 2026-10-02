# repro/env_info.R -- record what a reproducibility run ran on.
# Usage: Rscript env_info.R <repo> <out_file>
args <- commandArgs(trailingOnly = TRUE)
repo <- args[1]; out <- args[2]
git <- function(...) tryCatch(system2("git", c("-C", repo, ...), stdout = TRUE, stderr = TRUE),
                             error = function(e) "git unavailable")
pk <- c("phyloseq", "vegan", "iNEXT.3D", "brms", "cmdstanr", "rstan", "Hmsc", "gllvm",
        "lme4", "ggplot2", "patchwork", "dplyr", "knitr", "dada2", "Biostrings")
ver <- vapply(pk, function(p) if (requireNamespace(p, quietly = TRUE))
  as.character(utils::packageVersion(p)) else "not installed", character(1))
sink(out)
cat("run date    :", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("host        :", Sys.info()[["nodename"]], "| cores", parallel::detectCores(), "\n")
cat("git commit  :", git("rev-parse", "HEAD"), "\n")
cat("git branch  :", git("rev-parse", "--abbrev-ref", "HEAD"), "\n")
st <- git("status", "--porcelain")
cat("git dirty   :", if (length(st)) paste(length(st), "paths") else "clean", "\n")
cat("quarto      :", tryCatch(system2("quarto", "--version", stdout = TRUE), error = function(e) NA), "\n")
cat("BLAS caps   : OPENBLAS", Sys.getenv("OPENBLAS_NUM_THREADS"), "| OMP", Sys.getenv("OMP_NUM_THREADS"), "\n\n")
cat("key packages:\n"); print(data.frame(version = ver), right = FALSE)
cat("\n"); print(sessionInfo())
sink()
