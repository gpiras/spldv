##### Functions for sbinaryLGMM ####

#' Estimation of SAR for binary models using LGMM
#' 
#' @name sbinaryLGMM
#' @import Matrix car stats spatialreg methods
#' @export 
sbinaryLGMM <- function(formula, data, subset, na.action, 
                        listw = NULL, 
                        nins  = 2, # Number of instruments
                        link  = c("logit", "probit"),
                        ...){
  # Based on McMillen's code 
  link  <- match.arg(link)
  
  # Spatial weight matrix
  if (!inherits(listw, "listw")) stop("No neighbourhood list")
  listw       <- as(listw, "CsparseMatrix")
  
  # Model frame
  callT    <- match.call(expand.dots = TRUE)
  callF    <- match.call(expand.dots = FALSE)
  mf       <- callT
  m        <- match(c("formula", "data"), names(mf), 0L)
  mf       <- mf[c(1L, m)]
  f1       <- Formula(formula)
  #Check if there exists WX variables 
  if (length(f1)[2L] == 2L) mixed <- TRUE else mixed <- FALSE 
  mf$formula <- f1 
  mf[[1L]] <- as.name("model.frame")
  mf       <- eval(mf, parent.frame())
  nframe   <- length(sys.calls())
  
  # Get variables and checks
  y  <- model.response(mf)
  if (any(is.na(y))) stop("NAs in dependent variable")
  X  <- model.matrix(f1, data = mf, rhs = 1)
  if (mixed){
    x.for.w <- model.matrix(f1, data = mf, rhs = 2)
    name.wx <- colnames(x.for.w)
    WX      <- crossprod(t(listw), x.for.w)
    name.wx <- name.wx[which(name.wx != "(Intercept)")]
    WX      <- WX[, name.wx, drop = FALSE]
    #colnames(WX) <- paste("W", name.wx, sep = "_")
    colnames(WX) <- paste0("W.", name.wx)
    if (any(is.na(WX))) stop("NAs in WX variable")
    X <- cbind(X, WX)
  }
  if (any(is.na(X))) stop("NAs in dependent variables")
  N  <- nrow(X)
  K  <- ncol(X)
  sn <- nrow(listw)
  if (N != sn) stop("Number of spatial units in W is different to the number of data")
  
  ## Generate initial instruments
  Z <- make.instruments(listw, x = X, q = nins)
  H <- cbind(X, Z)
  # Let linearly independent columns
  H <- H[, qr(H)$pivot[seq_len(qr(H)$rank)]]
  P <- ncol(H)
  if (P < K) stop("Underspecified model")
  
  ## Starting values
  sbinary <- glm(y ~ as.matrix(X) - 1, family = binomial(link = link), data = mf)
  b_init  <- sbinary$coef # Initial values of beta from a standard binary model
  
  ## First step is standard binary of y on X and compute the generalized residuals 
  pfun <- switch(link,
                 "probit" = pnorm,
                 "logit"  = plogis)
  dfun <- switch(link,
                 "probit" = dnorm,
                 "logit"  = dlogis)
  ddfun <- switch(link,
                  "logit"  = function(x) (1 - 2 * pfun(x)) * pfun(x) * (1 - pfun(x)),
                  "probit" = function(x) -x * dnorm(x))  
  ai   <- as.vector(crossprod(t(X), b_init))
  q    <- 2*y - 1
  #fa   <- dfun(q*ai)
  fa   <- pmax(dfun(q*ai), .Machine$double.eps)
  Fa   <- pmax(pfun(q*ai), .Machine$double.eps)
  ffa  <- ddfun(q*ai) 
  u0   <- q * (fa/Fa)
  grad <- -1 * as.vector(q^2 * ((ffa * Fa - fa^2) / (Fa^2))) # Common vector of the derivative
  #grad  <- as.vector(u0 * (u0 + ai))
  Gb    <- grad * X
  Grho  <- grad * crossprod(t(listw), ai)
  #Grho  <- grad * (W %*% ai)
  G     <- cbind(Gb, Grho)
  
  # Second step
  Ghat <- H %*% solve(crossprod(H)) %*% crossprod(H, G) #Predicted values of G
  epsilon <- u0 + as.vector(crossprod(t(Gb), b_init))
  #epsilon <- u0 + Gb %*% b_init
  fit <- lm(epsilon ~ as.matrix(Ghat) + 0)
  bhat <- coef(fit)
  names(bhat) <- c(colnames(X), "rho") 
  
  ## Saving results
  out <- structure(
    list(
      coefficients = bhat, 
      call         = callF,
      X            = X, 
      H            = H, 
      y            = y, 
      listw        = listw,
      link         = link, 
      fit          = fit
    ), 
    class = "binlgmm"
  )
  out
}

#### S3 methods
#' @rdname sbinaryLGMM
#' @method vcov binlgmm
#' @import stats
#' @export 
vcov.binlgmm <- function(object, ...){
  V <- hccm(object$fit)
  colnames(V) <- rownames(V) <- names(coef(object$fit))
  return(V)
}

#' Get Model Summaries for use with "mtable" for objects of class binlgmm
#' 
#' A generic function to collect coefficients and summary statistics from a \code{binlgmm} object. It is used in \code{mtable}
#' 
#' @param obj a \code{binlgmm} object,
#' @param alpha level of the confidence intervals,
#' @param ... further arguments,
#' 
#' @details For more details see package \pkg{memisc}.
#' @import stats
#' @importFrom memisc getSummary
#' @export 
getSummary.binlgmm <- function(obj, alpha = 0.05, ...){
  smry <- summary(obj)
  coef <- smry$Coef
  lower <- coef[, 1] - coef[, 2] * qnorm(alpha/2)
  upper <- coef[, 1] + coef[, 2] * qnorm(alpha/2)
  coef <- cbind(coef, lower, upper)
  colnames(coef) <- c("est", "se", "stat", "p", "lwr", "upr")
  N <-  nrow(obj$X)
  sumstat <- c(logLik = NA, deviance = NA, AIC = NA, BIC = NA, N = N, 
               LR = NA, df = NA, p = NA, Aldrich.Nelson = NA, McFadden = NA, Cox.Snell = NA,
               Nagelkerke = NA)
  list(coef = coef, sumstat = sumstat, contrasts = obj$contrasts,
       xlevels = NULL, call = obj$call)
}

#' @rdname sbinaryLGMM
#' @method print binlgmm
#' @import stats
#' @export 
print.binlgmm <- function(x, 
                          digits = max(3, getOption("digits") - 3),
                          ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print.default(format(drop(x$coefficients), digits = digits), print.gap = 2,
                quote = FALSE)
  cat("\n")
  invisible(x)
}

#' @rdname sbinaryLGMM
#' @method summary binlgmm
#' @import stats
#' @export
summary.binlgmm <- function(object, ...){
  b <- coef(object$fit)
  names(b) <- names(object$coefficients)
  std.err <- sqrt(diag(vcov(object)))
  z <- b / std.err
  p <- 2 * (1 - pnorm(abs(z)))
  CoefTable <- cbind(b, std.err, z, p)
  colnames(CoefTable) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  object$CoefTable    <- CoefTable
  class(object)       <- c("summary.binlgmm", "binlgmm")
  return(object)
}

#' @rdname sbinaryLGMM
#' @method print summary.binlgmm
#' @import stats
#' @export
print.summary.binlgmm <- function(x,
                                  digits = max(3, getOption("digits") - 2),
                                  width = getOption("width"),
                                  ...)
{
  cat("        ------------------------------------------------------------\n")
  cat("                      SLM Binary Model by Linearized GMM \n")
  cat("        ------------------------------------------------------------\n")
  
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  
  cat("\nCoefficients:\n")
  printCoefmat(x$CoefTable, digits = digits, P.values = TRUE, has.Pvalue = TRUE)
  
  cat(paste("\nSample size:", signif(nrow(x$X), digits)))
  invisible(x)
}

