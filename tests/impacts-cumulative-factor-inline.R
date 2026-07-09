# Regression tests for cumulative impacts() with factor() terms written
# directly in formulas.
#
# These tests intentionally avoid pre-created factor columns in data. The
# current supported workflow is to write factors explicitly in the formula,
# e.g. factor(CP), factor(HIGHINC), or INC * factor(CP).
#
# Conservative goal: protect the current behavior that cumulative impacts with
# dydx = "exact" agree with dydx = "numeric" for inline factor terms. We do not
# test last cumulative = from.region for factors because these objects currently
# use different names/structures.

if (!requireNamespace("spdep", quietly = TRUE)) {
  cat("Package 'spdep' is not available; skipping cumulative factor-inline impacts tests.\n")
  q("no")
}

suppressPackageStartupMessages({
  library(spldv)
  library(spdep)
})

flatten_impacts <- function(x) {
  if (is.null(x)) {
    return(stats::setNames(numeric(0), character(0)))
  }
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (is.matrix(x)) {
    y <- as.numeric(x)
    rn <- rownames(x)
    cn <- colnames(x)
    if (is.null(rn)) rn <- as.character(seq_len(nrow(x)))
    if (is.null(cn)) cn <- as.character(seq_len(ncol(x)))
    names(y) <- as.vector(outer(rn, cn, paste, sep = "::"))
    return(y)
  }
  if (is.list(x)) {
    return(unlist(x, recursive = TRUE, use.names = TRUE))
  }
  y <- as.numeric(x)
  names(y) <- names(x)
  y
}

last_cumulative <- function(x) {
  if (is.null(x)) {
    return(stats::setNames(numeric(0), character(0)))
  }
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (is.matrix(x)) {
    y <- as.numeric(x[, ncol(x), drop = TRUE])
    names(y) <- rownames(x)
    if (is.null(names(y))) names(y) <- as.character(seq_along(y))
    return(y)
  }
  if (is.list(x)) {
    pieces <- lapply(names(x), function(nm) {
      y <- last_cumulative(x[[nm]])
      if (length(y)) names(y) <- paste(nm, names(y), sep = ".")
      y
    })
    return(unlist(pieces, recursive = FALSE, use.names = TRUE))
  }
  y <- as.numeric(x)
  names(y) <- names(x)
  y
}

assert_close_impacts <- function(exact, numeric, label, tol = 1e-6) {
  e <- flatten_impacts(exact)
  n <- flatten_impacts(numeric)

  if (!length(e) || !length(n)) {
    stop(label, ": empty impacts object", call. = FALSE)
  }
  if (any(!is.finite(e)) || any(!is.finite(n))) {
    stop(label, ": non-finite impact value", call. = FALSE)
  }
  if (!setequal(names(e), names(n))) {
    stop(label, ": exact and numeric cumulative impacts have different names", call. = FALSE)
  }

  common <- intersect(names(e), names(n))
  max_diff <- max(abs(e[common] - n[common]), na.rm = TRUE)
  if (!is.finite(max_diff) || max_diff > tol) {
    stop(label, ": exact and numeric cumulative impacts differ by ",
         format(max_diff, digits = 8), " > ", tol, call. = FALSE)
  }

  le <- last_cumulative(exact)
  ln <- last_cumulative(numeric)
  if (!setequal(names(le), names(ln))) {
    stop(label, ": exact and numeric last cumulative columns have different names", call. = FALSE)
  }
  common_last <- intersect(names(le), names(ln))
  max_diff_last <- max(abs(le[common_last] - ln[common_last]), na.rm = TRUE)
  if (!is.finite(max_diff_last) || max_diff_last > tol) {
    stop(label, ": exact and numeric last cumulative impacts differ by ",
         format(max_diff_last, digits = 8), " > ", tol, call. = FALSE)
  }

  invisible(max(max_diff, max_diff_last))
}

fit_gmm <- function(formula) {
  suppressWarnings(
    sbinaryGMM(
      formula   = formula,
      data      = dat,
      listw     = lw,
      link      = "probit",
      type      = "onestep",
      winitial  = "identity",
      s.matrix  = "iid",
      nins      = 1,
      fastmom   = TRUE,
      verbose   = FALSE
    )
  )
}

raw_cumulative <- function(fit, variable, dydx, Q = 8L, from.unit = 20L) {
  spldv:::dydx.bingmm(
    theta = coef(fit), obj = fit, data = fit$data,
    result = "cumulative", variable = variable, from.unit = from.unit,
    dydx = dydx, Q = Q, het = TRUE, verbose = FALSE
  )
}

run_case <- function(formula, variables, label, Q = 8L) {
  fit <- fit_gmm(formula)
  for (variable in variables) {
    exact <- raw_cumulative(fit, variable = variable, dydx = "exact", Q = Q)
    numeric <- raw_cumulative(fit, variable = variable, dydx = "numeric", Q = Q)
    assert_close_impacts(exact, numeric, paste(label, variable, "Q", Q))
  }
  invisible(TRUE)
}

data(oldcol, package = "spdep")
dat <- COL.OLD
dat$CRIMED <- as.numeric(dat$CRIME > 35)
dat$HIGHINC <- dat$INC > stats::median(dat$INC, na.rm = TRUE)
lw <- spdep::nb2listw(COL.nb, style = "W", zero.policy = TRUE)

run_case(CRIMED ~ INC + HOVAL + factor(CP),
         variables = "CP",
         label = "factor(CP)")

run_case(CRIMED ~ INC * factor(CP) + HOVAL,
         variables = c("CP", "INC"),
         label = "INC * factor(CP)")

run_case(CRIMED ~ INC + HOVAL + factor(HIGHINC),
         variables = "HIGHINC",
         label = "factor(HIGHINC)")

cat("cumulative factor-inline impacts regression tests passed\n")
