# =============================================================================
# build_main_figure.R -- rebuild a main-text figure from the appendix itself.
#
# Main-text figures have no script of their own: the code chunk in
# Supplementary_Appendix.qmd IS the source of truth (Section 11). Those chunks
# are `eval: false`, so rendering never re-runs them; this runner pulls one out
# by its chunk label and executes it.
#
# Usage, FROM THIS DIRECTORY (so the chunk's relative paths resolve exactly as
# they do at render time):
#   conda run -n r_env Rscript build_main_figure.R fig1-build
#
# knitr::purl() cannot be used for this: it comments out `eval: false` chunks.
# =============================================================================

label <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(label)) stop("usage: Rscript build_main_figure.R <chunk-label>")

qmd <- "Supplementary_Appendix.qmd"
if (!file.exists(qmd)) stop("run this from the Supplementary/ directory (no ", qmd, " here)")

L    <- readLines(qmd, warn = FALSE)
open <- grep(sprintf("^```\\{r %s\\}$", label), L)
if (length(open) != 1L)
  stop(sprintf("expected exactly 1 chunk labelled '%s', found %d", label, length(open)))
close <- grep("^```$", L)
close <- close[close > open][1]
if (is.na(close)) stop("unterminated chunk: ", label)

code <- L[(open + 1):(close - 1)]
code <- code[!grepl("^#\\|", code)]              # drop the chunk options

cat(sprintf("Running chunk '%s' (%d lines) from %s\n\n", label, length(code), qmd))
eval(parse(text = code), envir = new.env(parent = globalenv()))
cat("\nDone.\n")
