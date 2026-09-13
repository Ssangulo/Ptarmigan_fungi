# Shared helper: size-standardised Hill numbers via iNEXT.3D.
#
# SINGLE SOURCE for the appendix. Sourced by the `sec5-build` and `fig3-build`
# chunks of Supplementary_Appendix.qmd so Section 5 and main-text Figure 3
# cannot drift apart. `Scripts/6_diversity_analyses.R` carries a verbatim copy
# of the same two functions (it must stay standalone, like every numbered
# pipeline script) -- if you change one, change the other.
#
# Requires: phyloseq (for the otu_mat_of input convention) and iNEXT.3D.

otu_mat_of <- function(ps) {
  m <- as(phyloseq::otu_table(ps), "matrix")
  if (phyloseq::taxa_are_rows(ps)) m <- t(m)
  m
}

# ---- Shared helper: size-standardised Hill numbers via iNEXT.3D -------------
# Replaces the previous `rarefy_even_depth() %>% hillR::hill_taxa()` pattern
# used in Sections 2 and 3. Both compute the SAME estimand -- the Hill number
# of a sample of `level` reads drawn without replacement -- but they differ in
# how they get there:
#   hillR on a rarefied object = ONE random draw from that distribution, so the
#     answer moves with the rarefaction seed. Verified against 200 independent
#     draws: the committed draw sat a median 0.68-0.82 draw-SD from the
#     expectation (up to 3.3 SD; up to 12.7 OTUs at q0).
#   estimate3D(base="size") = the ANALYTIC EXPECTATION, computed from the full
#     raw counts. Same check put it 0.05-0.07 draw-SD from the expectation.
# So this is the same quantity estimated ~10x more precisely, deterministically,
# and without discarding ~97% of the reads. Input must be RAW counts.
#
# `level` defaults to the minimum library size, i.e. exactly the depth the old
# rarefaction used, so the two are directly comparable.
#
# GOTCHA: iNEXT.3D refuses any assemblage with fewer than 5 observed species
# ("the number of observed species should be at least five"). A handful of
# low-richness rows here trip that. Rather than drop them (which would silently
# change n), those are computed from the exact closed-form rarefaction
# expectations of Chao et al. 2014 (Ecol Monogr 84:45-67) -- Hurlbert's formula
# for q0, the expected plug-in Shannon entropy for q1, and the exact expected
# inverse Simpson for q2. Verified against estimate3D on rows where both run:
# max absolute difference 8.5e-10, i.e. the identical estimand.
hill_rarefy_exact <- function(x, m) {
  x <- x[x > 0]; n <- sum(x)
  stopifnot(m <= n)
  if (m == n) { p <- x/n
    return(c(q0 = length(x), q1 = exp(-sum(p*log(p))), q2 = 1/sum(p^2))) }
  q0 <- sum(1 - exp(lchoose(n - x, m) - lchoose(n, m)))
  lcnm <- lchoose(n, m); H <- 0
  for (xi in x) {                       # X_i ~ Hypergeometric(n, x_i, m)
    k  <- seq.int(max(1L, m - (n - xi)), min(xi, m))
    kk <- k/m
    H  <- H + sum(exp(lchoose(xi, k) + lchoose(n - xi, m - k) - lcnm) * (-kk*log(kk)))
  }
  q2 <- 1/(1/m + (1 - 1/m) * sum(x*(x-1))/(n*(n-1)))
  c(q0 = q0, q1 = exp(H), q2 = q2)
}

hill_inext_size <- function(ps, level = NULL, label = "") {
  m <- otu_mat_of(ps)                                   # rows = assemblages
  if (is.null(level)) level <- min(rowSums(m))
  stopifnot(level >= 1, all(rowSums(m) >= level))
  S   <- rowSums(m > 0)
  big <- rownames(m)[S >= 5]; small <- rownames(m)[S < 5]

  out <- data.frame(row = rownames(m), q0 = NA_real_, q1 = NA_real_, q2 = NA_real_,
                    stringsAsFactors = FALSE)
  if (length(big)) {
    d <- as.data.frame(t(m[big, , drop = FALSE]))
    d <- d[rowSums(d) > 0, , drop = FALSE]
    e <- iNEXT.3D::estimate3D(d, diversity = "TD", q = c(0,1,2), datatype = "abundance",
                              base = "size", level = level, nboot = 0)
    w <- reshape(e[, c("Assemblage","Order.q","qTD")], idvar = "Assemblage",
                 timevar = "Order.q", direction = "wide")
    names(w) <- c("row","q0","q1","q2")
    i <- match(w$row, out$row); out[i, c("q0","q1","q2")] <- w[, c("q0","q1","q2")]
  }
  for (r in small)
    out[match(r, out$row), c("q0","q1","q2")] <- as.list(hill_rarefy_exact(m[r, ], level))

  stopifnot(!anyNA(out[, c("q0","q1","q2")]), identical(out$row, rownames(m)))
  cat(sprintf("%sSize-standardised Hill numbers at m = %d reads: %d assemblages (%d via estimate3D, %d via exact formula, S.obs < 5)\n",
              if (nzchar(label)) paste0(label, ": ") else "", level, nrow(out), length(big), length(small)))
  out
}
