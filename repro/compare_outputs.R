# repro/compare_outputs.R -- did the run regenerate the committed appendix
# files, and do they match?  Sourced by run_checks.R.
#
# Three comparisons:
#   (a) every tracked file under Supplementary/tables|figures vs git HEAD
#       (the run staged fresh copies over them in this checkout);
#   (b) every CSV the scripts wrote to the sandbox models/ vs the canonical
#       ~/Ptarmigan/models copy of the same name (covers script outputs the
#       appendix does not show directly);
#   (c) same-named CSVs written BOTH to the sandbox by a script and to
#       Supplementary/tables by an appendix chunk (double producers).

read_csv_safe <- function(path) tryCatch(
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
  error = function(e) NULL)

# Value-level CSV comparison. Returns status + detail.
#   identical      : bytes equal
#   same_values    : bytes differ, every value equal (formatting / row order of names)
#   within_tol     : max relative numeric difference <= tol
#   DRIFT          : larger numeric differences or changed text cells
#   STRUCTURE      : different columns / row count
compare_csv <- function(old, new, tol) {
  if (identical(unname(tools::md5sum(old)), unname(tools::md5sum(new))))
    return(list(status = "identical", detail = ""))
  a <- read_csv_safe(old); b <- read_csv_safe(new)
  if (is.null(a) || is.null(b)) return(list(status = "DRIFT", detail = "unreadable CSV"))
  if (!identical(names(a), names(b)) || nrow(a) != nrow(b))
    return(list(status = "STRUCTURE", detail = sprintf(
      "cols %d->%d, rows %d->%d; cols only old: %s; only new: %s", ncol(a), ncol(b),
      nrow(a), nrow(b), paste(head(setdiff(names(a), names(b)), 5), collapse = ","),
      paste(head(setdiff(names(b), names(a)), 5), collapse = ","))))
  max_rel <- 0; worst <- ""; txt_bad <- 0; txt_col <- ""
  for (j in names(a)) {
    x <- a[[j]]; y <- b[[j]]
    if (is.numeric(x) && is.numeric(y)) {
      na_mismatch <- sum(is.na(x) != is.na(y))
      ok <- !is.na(x) & !is.na(y)
      d  <- abs(x[ok] - y[ok]) / pmax(abs(x[ok]), abs(y[ok]), 1e-12)
      d[abs(x[ok] - y[ok]) < 1e-12] <- 0
      m <- if (length(d)) max(d) else 0
      if (na_mismatch > 0) m <- Inf
      if (m > max_rel) { max_rel <- m; worst <- j }
    } else {
      bad <- sum(as.character(x) != as.character(y), na.rm = TRUE) + sum(is.na(x) != is.na(y))
      if (bad > 0) { txt_bad <- txt_bad + bad; txt_col <- j }
    }
  }
  detail <- sprintf("max rel diff %.3g (col %s); %d text cells differ%s", max_rel,
                    if (nzchar(worst)) worst else "-", txt_bad,
                    if (txt_bad) paste0(" (e.g. ", txt_col, ")") else "")
  status <- if (max_rel == 0 && txt_bad == 0) "same_values" else
            if (txt_bad == 0 && max_rel <= tol) "within_tol" else "DRIFT"
  list(status = status, detail = detail)
}

compare_png <- function(old, new) {
  if (identical(unname(tools::md5sum(old)), unname(tools::md5sum(new))))
    return(list(status = "identical", detail = ""))
  a <- tryCatch(png::readPNG(old), error = function(e) NULL)
  b <- tryCatch(png::readPNG(new), error = function(e) NULL)
  if (is.null(a) || is.null(b)) return(list(status = "VISUAL", detail = "unreadable PNG"))
  if (!identical(dim(a), dim(b)))
    return(list(status = "VISUAL", detail = sprintf("dims %s -> %s",
      paste(dim(a), collapse = "x"), paste(dim(b), collapse = "x"))))
  dd <- abs(a - b)
  m  <- if (length(dim(dd)) == 3) do.call(pmax, lapply(seq_len(dim(dd)[3]), function(k) dd[, , k])) else dd
  frac <- mean(m > 1/255)
  list(status = if (frac == 0) "same_pixels" else "VISUAL",
       detail = sprintf("%.4f%% of pixels differ, max channel diff %.3f", 100 * frac, max(m)))
}

# PDFs embed creation dates / IDs; compare with those lines removed.
compare_pdf <- function(old, new) {
  if (identical(unname(tools::md5sum(old)), unname(tools::md5sum(new))))
    return(list(status = "identical", detail = ""))
  strip <- function(f) {
    r <- readBin(f, "raw", file.info(f)$size)
    txt <- rawToChar(r[r != as.raw(0)])
    gsub("/(CreationDate|ModDate) *\\([^)]*\\)|/ID *\\[[^]]*\\]", "", txt, useBytes = TRUE)
  }
  if (identical(strip(old), strip(new)))
    return(list(status = "same_content", detail = "differs only in creation date / ID"))
  # Compressed streams (cairo) can carry a timestamp the text strip can't see:
  # rasterise every page of both at 150 dpi and compare pixels instead.
  if (nzchar(Sys.which("pdftoppm"))) {
    d <- tempfile(); dir.create(d)
    system2("pdftoppm", c("-r", "150", "-png", shQuote(old), file.path(d, "old")))
    system2("pdftoppm", c("-r", "150", "-png", shQuote(new), file.path(d, "new")))
    po <- sort(list.files(d, "^old", full.names = TRUE)); pn <- sort(list.files(d, "^new", full.names = TRUE))
    if (length(po) && length(po) == length(pn)) {
      rs <- mapply(compare_png, po, pn, SIMPLIFY = FALSE)
      if (all(vapply(rs, function(r) r$status %in% c("identical", "same_pixels"), logical(1))))
        return(list(status = "same_pixels", detail = "PDF bytes differ; rendered pages identical at 150 dpi"))
      return(list(status = "VISUAL", detail = paste("rendered PDF differs:", rs[[1]]$detail)))
    }
  }
  list(status = "VISUAL", detail = "PDF content differs")
}

compare_jpg <- function(old, new) {
  if (identical(unname(tools::md5sum(old)), unname(tools::md5sum(new))))
    return(list(status = "identical", detail = ""))
  list(status = "VISUAL", detail = "JPEG bytes differ (check the matching PNG)")
}

# xlsx files embed timestamps: compare every sheet's values instead.
compare_xlsx <- function(old, new) {
  if (identical(unname(tools::md5sum(old)), unname(tools::md5sum(new))))
    return(list(status = "identical", detail = ""))
  if (!requireNamespace("readxl", quietly = TRUE))
    return(list(status = "VISUAL", detail = "xlsx differs (readxl unavailable)"))
  so <- readxl::excel_sheets(old); sn <- readxl::excel_sheets(new)
  if (!identical(so, sn)) return(list(status = "STRUCTURE", detail = "different sheets"))
  same <- all(vapply(so, function(s) isTRUE(all.equal(
    as.data.frame(readxl::read_excel(old, s)), as.data.frame(readxl::read_excel(new, s)))), logical(1)))
  if (same) list(status = "same_values", detail = "xlsx bytes differ (timestamps); every sheet equal")
  else list(status = "DRIFT", detail = "xlsx cell values differ")
}

compare_file <- function(old, new, tol) {
  ext <- tolower(tools::file_ext(new))
  switch(ext,
    csv = compare_csv(old, new, tol),
    png = compare_png(old, new),
    pdf = compare_pdf(old, new),
    jpg = , jpeg = compare_jpg(old, new),
    xlsx = compare_xlsx(old, new),
    if (identical(unname(tools::md5sum(old)), unname(tools::md5sum(new))))
      list(status = "identical", detail = "") else list(status = "VISUAL", detail = "binary differs"))
}

compare_outputs <- function(repo, sandbox, run_start, tol) {
  git <- function(...) system2("git", c("-C", repo, ...), stdout = TRUE, stderr = FALSE)
  tracked <- git("ls-files", "Supplementary/tables", "Supplementary/figures")
  tmp <- tempfile(); dir.create(tmp)
  rows <- lapply(tracked, function(rel) {
    new <- file.path(repo, rel)
    if (!file.exists(new))
      return(data.frame(set = "appendix", file = rel, regenerated = FALSE,
                        status = "DELETED", detail = "file missing after run"))
    regen <- file.info(new)$mtime >= run_start
    old <- file.path(tmp, gsub("/", "__", rel))
    system2("git", c("-C", repo, "show", paste0("HEAD:", rel)), stdout = old)
    r <- compare_file(old, new, tol)
    data.frame(set = "appendix", file = rel, regenerated = regen,
               status = r$status, detail = r$detail)
  })
  # untracked new files the run produced in Supplementary/tables|figures
  untracked <- git("ls-files", "--others", "--exclude-standard",
                   "Supplementary/tables", "Supplementary/figures")
  if (length(untracked)) rows <- c(rows, list(data.frame(
    set = "appendix", file = untracked, regenerated = TRUE,
    status = "NEW_UNTRACKED", detail = "written by the run, not committed")))

  # (b) sandbox script outputs vs canonical models/
  sb_csv <- list.files(file.path(sandbox, "models"), pattern = "\\.csv$", full.names = TRUE)
  canon_dir <- "/home/daniel/Ptarmigan/models"
  rows_b <- lapply(sb_csv, function(f) {
    cf <- file.path(canon_dir, basename(f))
    if (!file.exists(cf))
      return(data.frame(set = "sandbox_vs_canonical", file = basename(f), regenerated = TRUE,
                        status = "NO_CANONICAL", detail = "no same-named file in canonical models/"))
    r <- compare_csv(cf, f, tol)
    data.frame(set = "sandbox_vs_canonical", file = basename(f), regenerated = TRUE,
               status = r$status, detail = r$detail)
  })

  # (c) double producers: script copy (sandbox) vs appendix-chunk copy (Supplementary)
  supp_csv <- file.path(repo, "Supplementary/tables", basename(sb_csv))
  dbl <- file.exists(supp_csv)
  rows_c <- lapply(which(dbl), function(i) {
    r <- compare_csv(supp_csv[i], sb_csv[i], tol)
    data.frame(set = "script_vs_appendix_copy", file = basename(sb_csv[i]), regenerated = TRUE,
               status = r$status, detail = r$detail)
  })
  out <- do.call(rbind, c(rows, rows_b, rows_c))
  rownames(out) <- NULL
  out
}
