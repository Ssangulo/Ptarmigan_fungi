# repro/compare_workspace.R -- does a fresh 4_data_prep.R run rebuild the SAME
# analysis objects as the canonical eco_analysis.RData?
#
# Usage: Rscript compare_workspace.R <canonical.RData> <regenerated.RData> <out.csv>
#
# Exit 0 only if every GATE object is identical (all.equal on the phyloseq
# components). Objects that later scripts add to the canonical workspace
# (funguild_otu, fungaltraits_otu from script 8) are reported, not failed.
# Other differing objects (intermediates, helper functions) are reported as
# INFO -- downstream scripts only use the gate objects.
suppressMessages(library(phyloseq))
args <- commandArgs(trailingOnly = TRUE)
canon <- new.env(); regen <- new.env()
load(args[1], envir = canon); load(args[2], envir = regen)

GATE <- c("alldat", "alldat.rfy", "alldat_full", "alldat_rg2_sensitivity")
ADDED_LATER <- c("funguild_otu", "fungaltraits_otu")

ps_parts <- function(ps) list(
  otu = as(otu_table(ps), "matrix"), taxa_are_rows = taxa_are_rows(ps),
  tax = if (!is.null(tax_table(ps, FALSE))) as(tax_table(ps), "matrix"),
  sam = if (!is.null(sample_data(ps, FALSE))) data.frame(sample_data(ps)),
  refseq = if (!is.null(refseq(ps, FALSE))) as.character(refseq(ps)),
  tree = !is.null(phy_tree(ps, FALSE)))
norm <- function(x) {
  if (inherits(x, "phyloseq")) return(ps_parts(x))
  if (is.list(x) && !is.data.frame(x)) return(lapply(x, norm))
  if (is.function(x)) return(deparse(x))
  x
}
cmp <- function(a, b) {
  r <- all.equal(norm(a), norm(b), check.attributes = TRUE)
  if (isTRUE(r)) "identical" else paste(head(r, 3), collapse = " | ")
}

nm <- union(ls(canon, all.names = TRUE), ls(regen, all.names = TRUE))
res <- do.call(rbind, lapply(nm, function(n) {
  inC <- exists(n, canon, inherits = FALSE); inR <- exists(n, regen, inherits = FALSE)
  detail <- if (inC && inR) cmp(get(n, canon), get(n, regen)) else
            if (inC) "only in canonical" else "only in regenerated"
  role <- if (n %in% GATE) "GATE" else if (n %in% ADDED_LATER) "added_by_script8" else "other"
  status <- if (detail == "identical") "PASS" else
            if (role == "GATE") "FAIL" else
            if (role == "added_by_script8" && detail == "only in canonical") "EXPECTED" else "INFO"
  data.frame(object = n, role = role, status = status, detail = substr(detail, 1, 300))
}))
res <- res[order(match(res$role, c("GATE", "added_by_script8", "other")), res$object), ]
write.csv(res, args[3], row.names = FALSE)
print(res[res$role != "other" | res$status != "PASS", ], row.names = FALSE)
cat(sprintf("\n%d objects compared: %d identical, %d differ/unpaired\n",
            nrow(res), sum(res$status == "PASS"), sum(res$status != "PASS")))
miss <- setdiff(GATE, res$object[res$status == "PASS"])
if (length(miss)) { cat("GATE FAILED for:", paste(miss, collapse = ", "), "\n"); quit(status = 1) }
cat("GATE PASSED: all canonical analysis objects reproduce exactly.\n")
