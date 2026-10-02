# repro/check_claims.R -- do the numbers QUOTED in the appendix prose still
# match the tables they come from?  Sourced by run_checks.R.
#
# (1) Registry. claims.csv, one row per headline number:
#   id        short unique id (S4_permanova_season_R2)
#   section   appendix section (4.2)
#   file      table it comes from, relative to Supplementary/
#   expr      R expression returning ONE number; the table is bound to `d`
#             (e.g. d$R2[d$term == "Season"], or 100 * mean(d$x[d$Season == "winter"]));
#             tab("tables/other.csv") reads a second table when a claim joins two
#   digits    decimal places the prose rounds to
#   quoted    the number exactly as the prose writes it (no thousands commas)
#   snippet   a verbatim piece of the qmd containing that number
# A claim PASSES when round(value, digits) equals `quoted` AND the snippet is
# still in the qmd -- so it fails whether the data moved or the prose was edited.
#
# (2) Sweep. Every other number in the prose and captions is matched against
# the tables at the precision it is quoted to; the leftovers go to
# unmatched_numbers.csv for a human to triage (most are sample sizes,
# thresholds or literature values that need no table).

check_registry <- function(repo) {
  qmd <- paste(readLines(file.path(repo, "Supplementary/Supplementary_Appendix.qmd"), warn = FALSE),
               collapse = "\n")
  cl <- utils::read.csv(file.path(repo, "repro/claims.csv"), stringsAsFactors = FALSE,
                        na.strings = "", colClasses = "character")
  cache <- list()
  get_tab <- function(f) {
    if (is.null(cache[[f]])) cache[[f]] <<- utils::read.csv(file.path(repo, "Supplementary", f),
                                                           check.names = FALSE, stringsAsFactors = FALSE)
    cache[[f]]
  }
  out <- lapply(seq_len(nrow(cl)), function(i) {
    r <- cl[i, ]
    val <- tryCatch({
      d <- get_tab(r$file)
      x <- eval(parse(text = r$expr), list(d = d, tab = get_tab), globalenv())
      x <- suppressWarnings(as.numeric(gsub("[,%]", "", as.character(x))))
      if (length(x) != 1) stop(sprintf("expr returned %d values", length(x)))
      x
    }, error = function(e) structure(NA_real_, err = conditionMessage(e)))
    err <- attr(val, "err")
    dg  <- as.integer(r$digits)
    q   <- as.numeric(gsub("[,−]", "", sub("^−", "-", r$quoted)))
    # inclusive half-unit: a value exactly on a rounding boundary (often a CSV
    # that itself stores a rounded number) cannot contradict either rounding
    value_ok <- !is.na(val) && isTRUE(abs(val - q) <= 10^(-dg) / 2 + 1e-9)
    tie <- value_ok && isTRUE(abs(abs(val - q) - 10^(-dg) / 2) < 1e-9)
    snip_ok  <- !is.na(r$snippet) && grepl(r$snippet, qmd, fixed = TRUE)
    data.frame(id = r$id, section = r$section, file = r$file, quoted = r$quoted,
               value = if (is.na(val)) NA_character_ else format(signif(val, 8)),
               value_ok = value_ok, snippet_ok = snip_ok,
               status = if (value_ok && snip_ok) "PASS" else "FAIL",
               detail = if (!is.null(err)) err else
                        if (tie && snip_ok) "rounding tie (table precision cannot resolve the last digit)" else
                        if (!value_ok) sprintf("table gives %s at %d dp", format(round(val, dg), nsmall = dg), dg) else
                        if (!snip_ok) "snippet no longer in the qmd (prose edited?)" else "")
  })
  do.call(rbind, out)
}

# ---- sweep -------------------------------------------------------------------
prose_units <- function(lines) {
  # Returns data.frame(line, section, text): prose lines outside chunks, plus
  # long string literals inside chunks (table captions).
  in_chunk <- FALSE; in_yaml <- FALSE; in_comment <- FALSE; sec <- "front"
  out <- list()
  for (i in seq_along(lines)) {
    l <- lines[i]
    if (i == 1 && l == "---") { in_yaml <- TRUE; next }
    if (in_yaml) { if (l == "---") in_yaml <- FALSE; next }
    if (grepl("^```", l)) { in_chunk <- !in_chunk && grepl("^```\\{", l); next }
    if (in_chunk) {
      s <- regmatches(l, gregexpr('"[^"]{40,}"', l))[[1]]
      if (length(s)) out[[length(out) + 1]] <- data.frame(line = i, section = sec, kind = "caption",
                                                          text = paste(s, collapse = " "))
      next
    }
    if (grepl("<!--", l)) in_comment <- TRUE
    if (in_comment) { if (grepl("-->", l)) in_comment <- FALSE; next }
    h <- regmatches(l, regexpr("^#{2,3} +(Section +)?[0-9]+(\\.[0-9]+)?", l))
    if (length(h)) sec <- sub("^#{2,3} +(Section +)?", "", h)
    if (grepl("^#|^!\\[|^:::", l) || !nzchar(trimws(l))) next
    out[[length(out) + 1]] <- data.frame(line = i, section = sec, kind = "prose", text = l)
  }
  do.call(rbind, out)
}

extract_numbers <- function(txt) {
  t <- gsub("`[^`]*`", " ", txt)                                  # inline code / filenames
  t <- gsub("\\]\\([^)]*\\)", "]", t)                              # link targets
  t <- gsub("(Sections?|§|Tables?|Figures?|Figs?\\.|Fig|Data|OTU|PC|PCoA|q|S|H|E|R|SH|lv|v)\\s?[0-9][0-9A-Za-z.–-]*", " ", t, perl = TRUE)
  t <- gsub("\\b(19|20)[0-9]{2}(-[0-9]{2}-[0-9]{2})?\\b", " ", t, perl = TRUE)  # years / dates
  t <- gsub("[a-zA-Z][0-9]+", " ", t, perl = TRUE)                 # ids like rclr0, log10, R2
  m <- regmatches(t, gregexpr("(?<![0-9.])[−-]?[0-9]{1,3}(,[0-9]{3})+(\\.[0-9]+)?|(?<![0-9.])[−-]?[0-9]+(\\.[0-9]+)?", t, perl = TRUE))[[1]]
  m
}

table_values <- function(repo, files) {
  v <- unlist(lapply(files, function(f) {
    p <- file.path(repo, "Supplementary", f)
    if (!file.exists(p) || !grepl("\\.csv$", f)) return(NULL)
    d <- tryCatch(utils::read.csv(p, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(d)) return(NULL)
    unlist(lapply(d, function(x) suppressWarnings(as.numeric(gsub("[,%]", "", as.character(x))))))
  }))
  v <- v[is.finite(v)]
  unique(c(v, 100 * v))
}

sweep_numbers <- function(repo, registry) {
  lines <- readLines(file.path(repo, "Supplementary/Supplementary_Appendix.qmd"), warn = FALSE)
  pu <- prose_units(lines)
  all_files <- unique(unlist(regmatches(lines, gregexpr("tables/[A-Za-z0-9_./-]+\\.csv", lines))))
  # files referenced inside each top-level section
  top <- sub("\\..*$", "", pu$section)
  sec_files <- tapply(seq_len(nrow(pu)), top, function(ix) {
    rng <- range(pu$line[ix]); seg <- lines[rng[1]:rng[2]]
    unique(unlist(regmatches(seg, gregexpr("tables/[A-Za-z0-9_./-]+\\.csv", seg))))
  }, simplify = FALSE)
  pool_all <- table_values(repo, all_files)
  pool_sec <- lapply(sec_files, function(f) table_values(repo, f))
  reg_snips <- registry$snippet
  convention <- c("95", "0.95", "0.05", "2.5", "97.5", "50", "100", "0.5", "1.1", "1.01")
  rows <- list()
  for (k in seq_len(nrow(pu))) {
    nums <- extract_numbers(pu$text[k])
    for (s in nums) {
      n <- as.numeric(gsub("[,]", "", sub("^−", "-", s)))
      if (is.na(n)) next
      dp <- if (grepl("\\.", s)) nchar(sub("^.*\\.", "", s)) else 0
      if (dp == 0 && abs(n) < 10) next                                  # trivial small integers
      covered <- any(vapply(reg_snips, function(sn) grepl(s, sn, fixed = TRUE) &&
                              grepl(sn, pu$text[k], fixed = TRUE), logical(1)))
      hit <- function(pool) any(abs(round(pool, dp) - n) < 10^(-dp) / 2 |
                                abs(round(abs(pool), dp) - abs(n)) < 10^(-dp) / 2)
      cat <- if (covered) "registry" else
             if (s %in% convention) "convention" else
             if (!is.null(pool_sec[[top[k]]]) && hit(pool_sec[[top[k]]])) "section_table" else
             if (hit(pool_all)) "other_table" else "UNMATCHED"
      rows[[length(rows) + 1]] <- data.frame(line = pu$line[k], section = pu$section[k],
        kind = pu$kind[k], number = s, category = cat,
        context = substr(gsub("\\s+", " ", pu$text[k]), max(1, regexpr(s, pu$text[k], fixed = TRUE) - 60),
                         regexpr(s, pu$text[k], fixed = TRUE) + 60))
    }
  }
  do.call(rbind, rows)
}
