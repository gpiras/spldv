### Marginal Effects (Impacts) ----

#' @title Estimation of the average marginal effects for SARB models.
#' 
#' @aliases impacts
#' @import spatialreg
#' @export impacts
#' @title Estimation of the average marginal effects for SARB model estimated using GMM procedures.
#' 
#' @description Obtain the average marginal effects from \code{bingmm} or \code{binlgmm} class model.
#' 
#' @param obj an object of class \code{bingmm}, \code{binlgmm}.
#' @param object an object of class \code{impacts.bingmm} for \code{summary} methods.
#' @param x an object of class \code{impacts.bingmm} for \code{print} methods. 
#' @param vcov an estimate of the asymptotic variance-covariance matrix of the parameters for a \code{bingmm} or \code{binlgmm} object.
#' @param vce string indicating what kind of variance-covariance matrix of the estimate should be computed when using \code{effect.bingmm}. For the one-step GMM estimator, the options are \code{"robust"} and \code{"ml"}. For the two-step GMM estimator, the options are \code{"robust"}, \code{"efficient"} and \code{"ml"}. The option \code{"vce = ml"} is an exploratory method that evaluates the VC of the RIS estimator using the GMM estimates.
#' @param het logical. If \code{TRUE} (the default), then the heteroskedasticity is taken into account when computing the average marginal effects. 
#' @param atmeans logical. If \code{FALSE} (the default), then the average marginal effects are computed at the unit level. 
#' @param type string indicating which method is used to compute the standard errors of the average marginal effects. If \code{"mc"}, then the Monte Carlo approximation is used. If \code{"delta"}, then the Delta Method is used.
#' @param R numerical. Indicates the number of draws used in the Monte Carlo approximation if \code{type = "mc"}.
#' @param approximation logical. If \code{TRUE} then \eqn{(I - \lambda W)^{-1}} is approximated as \eqn{I + \lambda W + \lambda^2 W^2 + \lambda^3 W^3 + ... +\lambda^q W^q}. The default is \code{FALSE}.
#' @param pw numeric. The power used for the approximation \eqn{I + \lambda W + \lambda^2 W^2 + \lambda^3 W^3 + ... +\lambda^q W^q}. The default is 5.
#' @param tol Argument passed to \code{mvrnorm}: tolerance (relative to largest variance) for numerical lack of positive-definiteness in the coefficient covariance matrix.
#' @param empirical logical. Argument passed to \code{mvrnorm} (default \code{FALSE}): if \code{TRUE}, the coefficients and their covariance matrix specify the empirical not population mean and covariance matrix
#' @param digits the number of digits.
#' @param ... further arguments. Ignored.
#' @param x an object of class \code{impacts.bingmm}.
#' 
#' @details 
#' 
#' Let the model be:
#' 
#' \deqn{
#' y^*= X\beta + WX\gamma + \lambda W y^* + \epsilon = Z\delta + \lambda Wy^{*} + \epsilon
#' }
#' 
#' where  \eqn{y = 1} if \eqn{y^*>0} and 0 otherwise; \eqn{\epsilon \sim N(0, 1)} if \code{link = "probit"} or \eqn{\epsilon \sim L(0, \pi^2/3)} if \code{link = "logit"}.
#' 
#' The marginal effects respect to variable \eqn{x_r} can be computed as
#' 
#' \deqn{
#' diag(f(a))D^{-1}_{\lambda}A^{-1}_{\lambda}\left(I_n\beta_r + W\gamma_r\right) = C_r(\theta)
#' }
#' 
#' where \eqn{f()} is the pdf, which depends on the assumption of the error terms; \eqn{diag} is the operator that creates a \eqn{n \times n} diagonal matrix; \eqn{A_{\lambda}= (I -\lambda W)}; and \eqn{D_{\lambda}} is a diagonal matrix whose elements represent the square root of the diagonal elements of the variance-covariance matrix of \eqn{u = A_{\lambda}^{-1}\epsilon}. 
#' 
#' We implement these three summary measures: (1) The average total effects, \eqn{ATE_r  = n^{-1}i_n'C_{r}i_n}, (2) The average direct effects, \eqn{ADE_r  = n^{-1}tr(C_{r})}, and (3) the average indirect effects, \eqn{ATE_r - ADE_r}. 
#' 
#' The standard errors of the average total, direct and indirect effects can be estimated using either Monte Carlo (MC) approximation, which takes into account the sampling distribution of \eqn{\theta}, or Delta Method. 
#' 
#' @examples
#' \donttest{
#' # Data set
#' data(oldcol, package = "spdep")
#' 
#' # Create dependent (dummy) variable
#' COL.OLD$CRIMED <- as.numeric(COL.OLD$CRIME > 35)
#' 
#' # Two-step (Probit) GMM estimator
#' ts <- sbinaryGMM(CRIMED ~ INC + HOVAL| HOVAL,
#'                 link = "probit", 
#'                 listw = spdep::nb2listw(COL.nb, style = "W"), 
#'                 data = COL.OLD, 
#'                 type = "twostep")
#'                 
#' # Marginal effects using Delta Method
#' summary(impacts(ts, type = "delta"))
#' 
#' # Marginal effects using MC with 100 draws
#' summary(impacts(ts, type = "mc", R = 100))
#' 
#' # Marginal effects using efficient VC matrix
#' summary(impacts(ts, type = "delta", vce = "efficient"))
#' 
#' # Marginal effects using efficient VC matrix and ignoring the heteroskedasticity
#' summary(impacts(ts, type = "delta", vce = "efficient", het = FALSE))
#'}
#' @return An object of class \code{impacts.bingmm}. 
#' @seealso \code{\link[spldv]{sbinaryGMM}}, \code{\link[spldv]{sbinaryLGMM}}.
#' @author Mauricio Sarrias and Gianfranco Piras. 
#' @method impacts bingmm
#' @importFrom numDeriv jacobian
#' @importFrom MASS mvrnorm
#' @export
impacts.bingmm <- function(obj,
                          vcov = NULL,
                          vce = c("robust", "efficient", "ml"), 
                          het = TRUE,
                          atmeans = FALSE, 
                          type = c("mc", "delta"), 
                          R = 100,
                          approximation = FALSE,
                          pw  = 5, 
                          tol = 1e-06, 
                          empirical = FALSE,
                          ...){
  # Type of standard errors
  type  <- match.arg(type)
  vce   <- match.arg(vce)
  
  # Variance covariance matrix
  if (is.null(vcov)){
    V <- vcov(obj, vce = vce)
  } else {
    V <- vcov
    n.param <- length(coef(obj))
    if (dim(V)[1L] != n.param | dim(V)[2L] != n.param)  stop("dim of vcov are not the same as the estimated parameters")
  }
  mu <- coef(obj)
  
  # Make effects
  if (type == "delta"){
    me <- dydx.bingmm(coeff = mu, object = obj, het = het, atmeans = atmeans, approximation = approximation, pw = pw)
    # Make Jacobian (use numerical jacobian)
    jac <- numDeriv::jacobian(dydx.bingmm, mu, object = obj, het = het, atmeans = atmeans, approximation = approximation, pw = pw)
    se <- sqrt(diag(jac %*% V %*% t(jac))) 
  } else {
    W            <- obj$listw
    sym          <- all(W == t(W))
    omega        <- eigen(W, only.values = TRUE, symmetric = sym)
    interval     <- if (is.complex(omega$values)) 1 / range(Re(omega$values)) else 1 / range(omega$values)
    samples      <- MASS::mvrnorm(n = R, mu = mu, Sigma = V, tol = tol, empirical = empirical)
    check        <- ((samples[, length(mu)] > interval[1]) & (samples[, length(mu)] < interval[2]))
    if (any(!check)) samples <- samples[check, ]
    sres         <- apply(samples, 1, dydx.bingmm, object = obj, het = het, atmeans = atmeans, approximation = approximation, pw = pw)
    me           <- apply(sres, 1, mean)
    se           <- apply(sres, 1, sd)
  }
  
  # Save results
  z  <-  me / se
  p  <- 2 * pnorm(-abs(z))
  results        <- cbind(`dydx` = me, `Std. error` = se, `z value` = z, `Pr(> z)` = p)
  #object$margins <- results
  class(results)  <- c("impacts.bingmm")
  return(results)
}

#' @rdname impacts.bingmm
#' @method impacts binlgmm
#' @importFrom numDeriv jacobian
#' @importFrom MASS mvrnorm
#' @export
impacts.binlgmm <- function(obj,
                           vcov = NULL,
                           het = TRUE,
                           atmeans = FALSE, 
                           type = c("mc", "delta"), 
                           R = 100,
                           approximation = FALSE,
                           pw  = 5, 
                           tol = 1e-06, 
                           empirical = FALSE,
                           ...){
  # Type of standard errors
  type  <- match.arg(type)

  # Variance covariance matrix
  if (is.null(vcov)){
    V <- vcov(obj)
  } else {
    V <- vcov
    n.param <- length(coef(obj))
    if (dim(V)[1L] != n.param | dim(V)[2L] != n.param)  stop("dim of vcov are not the same as the estimated parameters")
  }
  mu <- coef(obj)
  
  # Make effects
  if (type == "delta"){
    me <- dydx.bingmm(coeff = mu, object = obj, het = het, atmeans = atmeans, approximation = approximation, pw = pw)
    # Make Jacobian (use numerical jacobian)
    jac <- numDeriv::jacobian(dydx.bingmm, mu, object = obj, het = het, atmeans = atmeans, approximation = approximation, pw = pw)
    se <- sqrt(diag(jac %*% V %*% t(jac))) 
  } else {
    W            <- obj$listw
    sym          <- all(W == t(W))
    omega        <- eigen(W, only.values = TRUE, symmetric = sym)
    interval     <- if (is.complex(omega$values)) 1 / range(Re(omega$values)) else 1 / range(omega$values)
    samples      <- MASS::mvrnorm(n = R, mu = mu, Sigma = V, tol = tol, empirical = empirical)
    check        <- ((samples[, length(mu)] > interval[1]) & (samples[, length(mu)] < interval[2]))
    if (any(!check)) samples <- samples[check, ]
    sres         <- apply(samples, 1, dydx.bingmm, object = obj, het = het, atmeans = atmeans, approximation = approximation, pw = pw)
    me           <- apply(sres, 1, mean)
    se           <- apply(sres, 1, sd)
  }
  
  # Save results
  z  <-  me / se
  p  <- 2 * pnorm(-abs(z))
  results        <- cbind(`dydx` = me, `Std. error` = se, `z value` = z, `Pr(> z)` = p)
  #object$margins <- results
  class(results)  <- c("impacts.bingmm")
  return(results)
}

#' @rdname impacts.bingmm
#' @method print impacts.bingmm
#' @export
print.impacts.bingmm <- function(x, ... ){
  cat("The average total effects are:\n")
  k <- nrow(x) / 3
  cat("Estimate(s):", x[1:k, 1], "\n")
}

#' @rdname impacts.bingmm
#' @method summary impacts.bingmm
#' @export
summary.impacts.bingmm <- function(object, ...){
  CoefTable      <- object
  summary        <- list(CoefTable = CoefTable)
  class(summary) <- "summary.impacts.bingmm"
  summary
}

#' @rdname impacts.bingmm
#' @method print summary.impacts.bingmm
#' @export
print.summary.impacts.bingmm <- function(x, digits = max(3, getOption("digits") - 3), ...){
  k <- nrow(x$CoefTable) / 3
  cat("------------------------------------------------------", fill = TRUE)
  cat("(a) Total effects :\n")
  cat("------------------------------------------------------",fill = TRUE)
  
  printCoefmat(x$CoefTable[1:k, , drop = FALSE], digits = digits, signif.legend = FALSE)
  
  cat("\n------------------------------------------------------", fill = TRUE)
  cat("(b) Direct effects :\n")
  cat("------------------------------------------------------",fill = TRUE)
  
  printCoefmat(x$CoefTable[(k + 1):(k*2), , drop = FALSE], digits = digits, signif.legend = FALSE)
  
  cat("\n------------------------------------------------------", fill = TRUE)
  cat("(c) Indirect effects :\n")
  cat("------------------------------------------------------",fill = TRUE)
  
  printCoefmat(x$CoefTable[(2*k + 1):(k*3), , drop = FALSE], digits = digits)
}

dydx.bingmm <- function(coeff, object, het = TRUE, atmeans = FALSE, approximation = TRUE, pw = 5, ...){
  ### FIXME: Check if the matrix is row-standardized
  theta.hat <- coeff
  lambda    <- theta.hat["lambda"]
  beta      <- theta.hat[1:length(theta.hat) - 1]
  cons      <- grep("(Intercept)", names(beta))
  has.cons  <- length(cons) > 0
  if (has.cons){
    betaNi <- beta[which(names(beta) != "(Intercept)")]
  } else {
    betaNi <- beta
  }
  dfun <- switch(object$link,
                 "probit" = dnorm,
                 "logit"  = dlogis)
  # to Make "a"
  X        <- object$X
  N        <- nrow(X)
  I        <- sparseMatrix(i = 1:N, j = 1:N, x = 1)
  W        <- object$listw
  A        <- I - lambda * W
  B        <- if (approximation) app_W(W, lambda, pw)  else solve(A)
  #trs      <- makeTraces(W, o = o, iiter = 50)
  #lambdavec <- lambda^(0:(o-1))
  if (atmeans) X <- matrix(rep(colMeans(X), nrow(X)), ncol = ncol(X), byrow = TRUE) # X is a n * k matrix of regressor means
  if (het){
    Sigma_u <- tcrossprod(B)
    Di      <- sparseMatrix(i = 1:N, j = 1:N, x = 1 / sqrt(diag(Sigma_u)))
    a       <- Di %*% B %*% X %*% beta
    fa      <- dfun(as.numeric(a))
    dfa     <- sparseMatrix(i = 1:N, j = 1:N, x = fa)
    #dfa    <- dfa %*% Di
    dM      <- dfa %*% Di %*% B
    #trs    <- makeTraces(dfa %*% Di %*% W, o = o, iiter = 50)
  } else {
    a        <- B %*% X %*% beta
    fa       <- dfun(as.numeric(a))
    dfa      <- sparseMatrix(i = 1:N, j = 1:N, x = fa)
    #trs       <- makeTraces(dfa %*% W, o = o, iiter = 50)
    dM       <- dfa %*% B
  }
  
  # Make effects
  b.names  <- names(betaNi)
  x.vars   <- all.vars(formula(object$formula, lhs = 0, rhs = 1))
  ones     <- rep(1, N)
  TE <- DE <- c()
  for (k in 1:length(x.vars)){
    name.wx  <- paste0("lag_", x.vars[k]) 
    has.wlag <-  any(b.names == name.wx)
    if (has.wlag) coef <- I * betaNi[x.vars[k]] + betaNi[name.wx]*W  else coef <- I*betaNi[x.vars[k]]
    ME       <- dM  %*% coef
    TE       <- c(TE, sum(ME) / N)
    #TE       <- c(TE, n^(-1) * t(ones)%*% (diag(dfa)) * betaNi[x.vars[k]] / (1- lambda))
    DE       <- c(DE, sum(diag(ME)) / N)
    #DE       <- c(DE, (t(diag(dfa)) %*% trs %*% lambdavec) * coef / N)
  }
  IE <- TE - DE
  names(TE) <- names(DE) <- names(IE) <- x.vars
  #out    <- list(TE = TE, DE = DE, IE = IE)
  out <- c(TE, DE, IE)
  return(out)
}


### From spatialprobit
# makeTraces <- function(W, o = 100, iiter = 50) 
# {
#   n <- nrow(W)
#   trW_i <- matrix(data = 0, nrow = n, ncol = o)
#   u <- matrix(rnorm(n * iiter), nrow = n, ncol = iiter)
#   xx <- u
#   trW_i[, 1] <- apply(u * as.matrix(xx), 1, sum)
#   for (i in 2:o) {
#     xx <- W %*% xx
#     trW_i[, i] <- apply(u * as.matrix(xx), 1, sum)
#   }
#   trW_i <- trW_i/iiter
#   return(trW_i)
# }
