# Regression tests for impacts() with factor() terms written directly in formulas.
# These tests intentionally avoid pre-created factor columns in data.  The
# current supported workflow is to write factors explicitly in the formula,
# e.g. factor(CP) or INC * factor(CP).

if (!requireNamespace("spdep", quietly = TRUE)) {
  cat("Package 'spdep' is not available; skipping factor-inline impacts tests.\n")
  q("no")
}

suppressPackageStartupMessages({
  library(spldv)
  library(spdep)
})

flatten_impacts <- function(x) {
  if (is.list(x)) {
    out <- unlist(x, recursive = TRUE, use.names = TRUE)
  } else {
    out <- as.numeric(x)
    names(out) <- names(x)
  }
  out
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

  common <- intersect(names(e), names(n))
  if (!length(common)) {
    stop(label, ": exact and numeric impacts have no common names", call. = FALSE)
  }

  if (!setequal(names(e), names(n))) {
    stop(label, ": exact and numeric impacts have different names", call. = FALSE)
  }

  max_diff <- max(abs(e[common] - n[common]), na.rm = TRUE)
  if (!is.finite(max_diff) || max_diff > tol) {
    stop(label, ": exact and numeric impacts differ by ",
         format(max_diff, digits = 8), " > ", tol, call. = FALSE)
  }

  invisible(max_diff)
}

run_case <- function(formula, variable, label, from.unit = 20L) {
  fit <- suppressWarnings(
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

  theta <- coef(fit)

  s_exact <- spldv:::dydx.bingmm(
    theta = theta, obj = fit, data = fit$data,
    result = "summary", dydx = "exact", het = TRUE, verbose = FALSE
  )
  s_numeric <- spldv:::dydx.bingmm(
    theta = theta, obj = fit, data = fit$data,
    result = "summary", dydx = "numeric", het = TRUE, verbose = FALSE
  )
  assert_close_impacts(s_exact, s_numeric, paste(label, "summary"))

  f_exact <- spldv:::dydx.bingmm(
    theta = theta, obj = fit, data = fit$data,
    result = "from.region", variable = variable, from.unit = from.unit,
    dydx = "exact", het = TRUE, verbose = FALSE
  )
  f_numeric <- spldv:::dydx.bingmm(
    theta = theta, obj = fit, data = fit$data,
    result = "from.region", variable = variable, from.unit = from.unit,
    dydx = "numeric", het = TRUE, verbose = FALSE
  )
  assert_close_impacts(f_exact, f_numeric, paste(label, "from.region"))

  invisible(TRUE)
}

data(oldcol, package = "spdep")
dat <- COL.OLD
dat$CRIMED <- as.numeric(dat$CRIME > 35)
lw <- spdep::nb2listw(COL.nb, style = "W", zero.policy = TRUE)

run_case(CRIMED ~ INC + HOVAL + factor(CP),
         variable = "CP",
         label = "factor(CP)")

run_case(CRIMED ~ INC * factor(CP) + HOVAL,
         variable = "CP",
         label = "INC * factor(CP)")

cat("factor-inline impacts regression tests passed\n")
