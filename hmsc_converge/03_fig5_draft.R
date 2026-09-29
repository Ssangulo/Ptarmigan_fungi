# =============================================================================
# hmsc_converge/03_fig5_draft.R  (branch exp/hmsc-converge -- experimental)
# Draft Figure 5 from a new fit WITHOUT copying or editing the figure code: the
# fig5-build chunk is read from main's Supplementary_Appendix.qmd, its input
# (tab_dir) and output (out_dir) are repointed at models|plots/hmsc_v2, and its
# drift guards are downgraded to messages -- so every guard that trips is a
# printed list of exactly which quoted numbers the new model moves.
#
# Usage (run from anywhere):
#   conda run -n r_env Rscript hmsc_converge/03_fig5_draft.R <tag>
# Needs 02_compare.R to have written models/hmsc_v2/<tag>_tables/ first.
# Writes plots/hmsc_v2/<tag>_Fig5_draft.{png,jpg,pdf} + <tag>_fig5_console.txt
# =============================================================================

tag  <- commandArgs(trailingOnly = TRUE)[1]
qmd  <- "/home/daniel/Ptarmigan/Scripts_server/Supplementary/Supplementary_Appendix.qmd"
tdir <- file.path("/home/daniel/Ptarmigan/models/hmsc_v2", paste0(tag, "_tables"))
odir <- "/home/daniel/Ptarmigan/plots/hmsc_v2"
stopifnot(dir.exists(tdir))

# The chunk also reads dark_taxa_SH_matching.csv (lineage + read share), which
# the new fit does not change -- link it into the v2 table dir.
supp_tab <- "/home/daniel/Ptarmigan/Scripts_server/Supplementary/tables"
file.copy(file.path(supp_tab, "dark_taxa_SH_matching.csv"), tdir, overwrite = TRUE)

lines <- readLines(qmd)
s <- grep("^```\\{r fig5-build\\}", lines)
e <- s + which(lines[(s + 1):length(lines)] == "```")[1]
code <- lines[(s + 1):(e - 1)]
code <- code[!grepl("^#\\|", code)]

sub1 <- function(pat, rep) {
  hit <- grep(pat, code)
  if (length(hit) != 1) stop("expected exactly one match for: ", pat)
  code[hit] <<- rep
}
sub1('^tab_dir <- if \\(dir.exists\\("tables"\\)\\)', sprintf('tab_dir <- "%s"; if (FALSE)', tdir))
sub1('^out_dir <- "figures/main"', sprintf('out_dir <- "%s"', odir))
sub1('^save_fig\\(file.path\\(out_dir, "Fig5_hmsc_synthesis"\\)',
     sprintf('save_fig(file.path(out_dir, "%s_Fig5_draft"), fig5, 11.5, 6.4)', tag))

env <- new.env()
env$stopifnot <- function(...) {
  cl <- deparse(sys.call(), width.cutoff = 500L)
  tryCatch(base::stopifnot(...),
           error = function(err) message("GUARD MOVED: ", conditionMessage(err), "\n    in: ",
                                         substr(paste(cl, collapse = " "), 1, 160)))
}
sink(file.path(odir, paste0(tag, "_fig5_console.txt")), split = TRUE)
eval(parse(text = code), envir = env)
sink()
cat("Draft written:", file.path(odir, paste0(tag, "_Fig5_draft.png")), "\n")
