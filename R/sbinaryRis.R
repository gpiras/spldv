##### Functions for sbinaryRis ####

#' @title Estimation of spatial probit model for binary outcomes using RIS (GHK) simulator
#' 
#' @description Estimation of spatial probit model  using RIS-normal (a.k.a GHK) simulator. The models can be the SAR or SEM probit model model. The SAR probit model has the following structure: 
#' 
#' \deqn{
#' y^*= X\beta + WX\gamma + \lambda W y^* + \epsilon = Z\delta + \lambda Wy^{*} + \epsilon = A_{\lambda}^{-1}Z\delta + u, 
#' }
#' where  \eqn{y = 1} if \eqn{y^*>0} and 0 otherwise, \eqn{Z = (X, WX)}, \eqn{\delta = (\beta', \gamma')'}, \eqn{u = A_{\lambda}^{-1}\epsilon} with \eqn{A_{\lambda} = (I - \lambda W)}, 
#' and \eqn{\epsilon \sim N(0, I)}. The SEM probit model has the following structure:
#' 
#' \deqn{
#' y^* = X\beta + WX\gamma + u = Z\delta  + u
#' }
#'where \eqn{y = 1} if \eqn{y^*>0} and 0 otherwise, \eqn{Z = (X, WX)}, \eqn{\delta = (\beta', \gamma')'}, 
#'\eqn{u = \rho W u + \epsilon} such that \eqn{u = A_{\rho}^{-1}\epsilon}, and \eqn{\epsilon \sim N(0, I)}, 
#'
#' @name sbinaryRis
#' @param formula a symbolic description of the model of the form \code{y ~ x | wx} where \code{y} is the binary dependent variable, \code{x} are the independent variables. The variables after \code{|} are those variables that enter spatially lagged: \eqn{WX}. The variables in the second part of \code{formula} must also appear in the first part. This rules out situations in which one of the regressors can be specified only in lagged form.
#' @param data the data of class \code{data.frame}.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param na.action a function which indicates what should happen when the data contains \code{NA}s. 
#' @param listw object. An object of class \code{listw}, \code{matrix}, or \code{Matrix}.  
#' @param R numerical. The number of draws used in RIS (GHK) simulator. 
#' @param model string. A string indicating which model to estimate. It can be \code{"SAR"} for the spatial autoregressive spatial model or \code{"SEM"} for the spatial error model. 
#' @param varcov string. A string indicating over which variance-covariance matrix to apply the Chokesly factorization. 
#' @param start if not \code{NULL}, the user must provide a vector of initial parameters for the optimization procedure. When \code{start = NULL}, \code{sbinaryRis} uses the traditional Probit estimates as initial values for the parameters, and the correlation between \eqn{y} and \eqn{Wy} as initial value for \eqn{\lambda} or \eqn{\rho}.
#' @param approximation logical. If \code{TRUE} (the default) then \eqn{(I - \lambda W)^{-1}} or \eqn{(I - \rho W)^{-1}} is approximated as \eqn{I + \lambda W + \lambda^2 W^2 + \lambda^3 W^3 + ... +\lambda^q W^q}. If  \code{FALSE}, then the inverse is computed without approximations.
#' @param pw numeric. The power used for the approximation \eqn{I + \lambda W + \lambda^2 W^2 + \lambda^3 W^3 + ... +\lambda^q W^q}. The default is 5.
#' @param print.init logical. If \code{TRUE} the initial parameters used in the optimization of the first step are printed. 
#' @param Qneg logical. Whether to construct the diagonal elements of \eqn{Q}. If \code{Qneg = FALSE}, then \eqn{q_{ii} = 2y_{i} - 1}. If  \code{Qneg = TRUE}, then \eqn{q_{ii} = 1- 2y_{i}}.
#' @param ... additional arguments passed to \code{maxLik}.
#' @param x,object,  an object of class \code{bingmm}.
#' @param eigentol the standard errors are only calculated if the ratio of the smallest and largest eigenvalue of the Hessian matrix is less than \code{eigentol}.  Otherwise the Hessian is treated as singular. 
#' @param digits the number of digits
#' @details 
#' 
#' The models are estimated by simulating the probabilities using the RIS-normal (GHK) simulator. The aim is to evaluate 
#' the multivariate density function \eqn{P(\upsilon = Q u < s)}, where \eqn{Q} is a diagonal matrix with entries \eqn{q_{ii} = 2y_i - 1} and the \eqn{n\times 1} vector
#' \eqn{s} depends on whether the model is SAR or SEM. If \code{model = "SAR"}, then \eqn{s = QA_{\lambda}^{-1}Z\delta} where \eqn{A_{\lambda} = (I - \lambda W)}; 
#' if \code{model = "SEM"}, then \eqn{s = QZ\delta}. 
#' 
#' Let \eqn{\Sigma_{\upsilon} = QVar(u)Q'} be the variance-covariance model of the transformed model. If \code{model = "SAR"}
#'  \eqn{\Sigma_{\upsilon} = Q\Sigma_{\lambda}Q'}, where \eqn{\Sigma_{\lambda} = (A_{\lambda}'A_{\lambda})^{-1}}. 
#'  If \code{model = "SEM"}, then \eqn{\Sigma_{\upsilon} = Q\Sigma_{\rho}Q'}, where \eqn{\Sigma_{\rho} = (A_{\rho}'A_{\rho})^{-1}}.
#'  
#'  Since \eqn{\Sigma_{\upsilon}} is positive definite, there exists a Cholesky decomposition such that
#'  \eqn{C'C = \Sigma_{\upsilon}^{-1}}, where \eqn{C} is the upper triangular Cholesky matrix and \eqn{\Sigma_{\upsilon}^{-1}} is the precision matrix. 
#'  Let \eqn{B = C^{-1}}. Then the random vector \eqn{Qu} can be replaced by \eqn{Qu = C^{-1}\eta = B\eta}, where \eqn{\eta} is
#'  a vector of standard normal variables. Then, the upper limit of the integration becomes \eqn{Qu = B \eta < s}, which 
#'  can be written as \eqn{\eta < B^{-1}s = \nu}. 
#'  
#'  The RIS simulator is implemented by drawing a large number \eqn{R} of random vector \eqn{\eta} and computing \eqn{\eta_{rj}} recursively
#'  for \eqn{j = 1, \ldots, n}. The parameters are estimated using the simulated maximum likelihood (SML) function:
#'  
#'  \deqn{\ln \tilde{L}(\theta) = \ln \left(\frac{1}{R}\sum_{r = 1}^R \tilde{p}_r\right), }
#'   
#'  where:
#'  
#'  \deqn{\tilde{p}_r = \prod_{j = 1}^n \Phi(\hat{\eta}_{jr}),}
#'  
#'  and \eqn{\Phi(\cdot)} is the univariate CDF of the standard normal density. 
#'  
#'  By default, \code{sbinaryRis} compute the SML using the Cholesky transformation on \eqn{\Sigma_{\upsilon}^{-1}}, \code{varcov = "invsigma"}. 
#'  The transformation can also be applied to \eqn{\Sigma_{\upsilon}} using \code{varcov = "sigma"}, which is slower than the previous option. 
#'  
#'  This estimator can take several minutes for large datasets. Thus, by default the inverse matrices \eqn{A_{\lambda}^{-1}} and \eqn{A_{\rho}^{-1}}
#'   are approximated using the Leontief expansion.
#'   
#' @author Mauricio Sarrias and Gianfranco Piras. 
#' @return An object of class "\code{binris}", a list with elements:
#' \item{varcov}{the matrix over which the Chokesly factorization is applied,}
#' \item{model}{type of model that was fitted,}
#' \item{Qneg}{matrix Q used,}
#' \item{approximation}{a logical value indicating whether approximation was used to compute the inverse matrix,}
#' \item{pw}{the powers for the approximation,}
#' \item{call}{the matched call,}
#' \item{X}{the X matrix, which contains also WX if the second part of the \code{formula} is used, }
#' \item{y}{the dependent variable,}
#' \item{listw}{the spatial weight matrix,}
#' \item{formula}{the formula,}
#' \item{R}{number of draws,}
#' \item{mf}{model frame.}
#' @examples
#' \donttest{
# Data set
#'data(oldcol, package = "spdep")
#'
#'# Create dependent (dummy) variable
#'COL.OLD$CRIMED <- as.numeric(COL.OLD$CRIME > 35)
#'
#'# Estimate SAR probit model using RIS simulator using Sigma_v^{-1}
#'ris_sar <- sbinaryRis(CRIMED ~ INC + HOVAL, 
#'                      data = COL.OLD,
#'                      R = 50, 
#'                      listw = spdep::nb2listw(COL.nb, style = "W"), 
#'                      print.level = 2)
#'summary(ris_sar)
#'
#'# Estimate SAR probit model using RIS simulator using Sigma_v
#'ris_sar_2  <- sbinaryRis(CRIMED ~ INC + HOVAL, 
#'                         data = COL.OLD,
#'                         R = 50, 
#'                         listw = spdep::nb2listw(COL.nb, style = "W"), 
#'                         varcov = "sigma", 
#'                         print.level = 2)
#'summary(ris_sar_2)
#'
#'# Estimate SDM probit model using RIS simulator
#'ris_sdm <- sbinaryRis(CRIMED ~ INC + HOVAL | INC + HOVAL, 
#'                      data = COL.OLD,
#'                      R = 50, 
#'                      listw = spdep::nb2listw(COL.nb, style = "W"), 
#'                      print.level = 2)
#'summary(ris_sdm)
#'
#'# Estimate SEM probit model using RIS simulator
#'ris_sem <- sbinaryRis(CRIMED ~ INC + HOVAL | INC, 
#'                      data = COL.OLD,
#'                      R = 50, 
#'                      listw = spdep::nb2listw(COL.nb, style = "W"), 
#'                      model = "SEM")
#'summary(ris_sem)
#'
#'}
#' @references 
#' 
#' Beron, K. J., & Vijverberg, W. P. (2004). Probit in a spatial context: a Monte Carlo analysis. In Advances in spatial econometrics: methodology, tools and applications (pp. 169-195). Berlin, Heidelberg: Springer Berlin Heidelberg.
#' 
#' Fleming, M. M. (2004). Techniques for estimating spatially dependent discrete choice models. In Advances in spatial econometrics (pp. 145-168). Springer, Berlin, Heidelberg.
#' 
#' Pace, R. K., & LeSage, J. P. (2016). Fast simulated maximum likelihood estimation of the spatial probit model capable of handling large samples. In Spatial Econometrics: Qualitative and Limited Dependent Variables (pp. 3-34). Emerald Group Publishing Limited.
#' 
#' Piras, G., & Sarrias, M. (2023). One or Two-Step? Evaluating GMM Efficiency for Spatial Binary Probit Models. Journal of choice modelling, 48, 100432. 
#' @seealso \code{\link[spldv]{sbinaryGMM}}, \code{\link[spldv]{impacts.binris}}.
#' @keywords models
#' @rawNamespace import(Matrix,  except = c(cov2cor, toeplitz, update)) 
#' @import stats methods Formula maxLik
#' @importFrom sphet listw2dgCMatrix
#' @export 
sbinaryRis <- function(formula, data, subset, na.action, 
                       listw = NULL, 
                       R = 20,                  
                       model = c("SAR", "SEM"),
                       varcov = c("invsigma", "sigma"), 
                       approximation = TRUE,
                       pw = 5,
                       start = NULL,
                       Qneg = FALSE,
                       print.init = FALSE, 
                       ...){
  # ============================-
  # 1. Obtain arguments ----
  # ============================-
  W      <- listw
  R      <- R
  model  <- match.arg(model)
  varcov <- match.arg(varcov) 
  
  # ============================-
  # 2. Spatial weight matrix (W): as CsparseMatrix ----
  # ============================-
  if(!inherits(listw,c("listw", "Matrix", "matrix"))) stop("Neighbourhood list or listw format unknown")
  if(inherits(listw,"listw"))   W    <- sphet::listw2dgCMatrix(listw)	
  if(inherits(listw,"matrix"))  W    <- Matrix(listw)	
  if(inherits(listw,"Matrix"))  W    <- listw	
  
  # ============================-
  # 3. Model frame ----
  # ============================-
  callT    <- match.call(expand.dots = TRUE)
  callF    <- match.call(expand.dots = FALSE)
  mf       <- callT
  m        <- match(c("formula", "data"), names(mf), 0L)
  mf       <- mf[c(1L, m)]
  f1       <- Formula::Formula(formula)
  if (length(f1)[2L] == 2L) Durbin <- TRUE else Durbin <- FALSE 
  mf$formula <- f1 
  mf[[1L]] <- as.name("model.frame")
  mf       <- eval(mf, parent.frame())
  nframe   <- length(sys.calls())
  
  # ============================-
  # 4. Optimization defaults ----
  # ============================-
  if (is.null(callT$method)) callT$method  <- 'bhhh'
  if (is.null(callT$iterlim)) callT$iterlim <- 1000
  
  # ============================-
  # 5. Obtain variables ----
  # ============================-
  y  <- model.response(mf)
  if (any(is.na(y))) stop("NAs in dependent variable")
  if (!all(y %in% c(0, 1, TRUE, FALSE))) stop("All dependent variables must be either 0, 1, TRUE or FALSE")
  if (!is.numeric(y)) y <- as.numeric(y)
  X  <- model.matrix(f1, data = mf, rhs = 1)
  if (Durbin){
    x.for.w  <- model.matrix(f1, data = mf, rhs = 2)
    name.wx  <- colnames(x.for.w)
    ## Check variables are in the first part of formula
    check.wx <- !(name.wx %in% colnames(X))
    if (sum(check.wx) > 0) stop("Some variables in WX do not appear in X. Check the formula")
    WX           <- crossprod(t(W), x.for.w)
    name.wx      <- name.wx[which(name.wx != "(Intercept)")]
    WX           <- WX[, name.wx, drop = FALSE] # Delete the constant from the WX
    colnames(WX) <- paste0("lag_", name.wx)
    if (any(is.na(WX))) stop("NAs in WX variable")
    X <- cbind(X, WX)
  }
  if (any(is.na(X))) stop("NAs in independent variables")
  N  <- nrow(X)
  K  <- ncol(X)
  sn <- nrow(W)
  if (N != sn) stop("Number of spatial units in W is different to the number of data")
  
  # ============================-
  # 5. Starting values ----
  # ============================-
  if (is.null(start)){
    sbinary      <- glm.fit(as.matrix(X), y, family = binomial(link = "probit"))
    b_init       <- sbinary$coef
    Wy           <- as.numeric(crossprod(t(W), y))
    lambda.init  <- cor(y, Wy)
    theta        <- c(b_init, lambda.init)
    names(theta) <- switch(model, 
                           "SAR" = c(colnames(X), "lambda"), 
                           "SEM" = c(colnames(X), "rho"))
  } else {
    theta <- start
    if (length(start) != length(c(colnames(X), "lambda"))) stop("Incorrect number of intial parameters")
    names(theta)  <- switch(model, 
                            "SAR" = c(colnames(X), "lambda"), 
                            "SEM" = c(colnames(X), "rho"))
  }
  if (print.init) {
    cat("\nStarting Values:\n")
    print(theta)
  } 
  
  # ============================-
  # 6. Optimization ----
  # ============================-
  opt <- callT
  opt$start <- theta
  m <- match(c('method', 'print.level', 'iterlim',
               'start', 'tol', 'ftol', 'steptol', 'fixed', 'constraints', 
               'control', 'finalHessian', 'reltol'),
             names(opt), 0L)
  opt <- opt[c(1L, m)]
  opt[[1]]     <- as.name('maxLik')
  if (model == "SAR") opt$logLik   <- as.name('lls_sar') else opt$logLik <- as.name('lls_sem')
  opt[c('y', 'X')]  <- list(as.name('y'), as.name('X'))
  opt$W             <- as.name('W')
  opt$R             <- as.name('R')
  opt$varcov        <- as.name('varcov')
  opt$Qneg          <- as.name('Qneg')
  opt$approximation <- as.name('approximation')
  opt$pw            <- as.name('pw')
  out <- eval(opt, sys.frame(which = nframe))
  
  # ============================-
  # 7. Save results ----
  # ============================-
  out$varcov        <- varcov
  out$model         <- model
  out$Qneg          <- Qneg
  out$approximation <- approximation
  out$pw            <- pw
  out$call          <- callT
  out$X             <- X
  out$y             <- y
  out$listw         <- W
  out$formula       <- f1
  out$R             <- R
  out$mf            <- mf
  class(out)        <- c("binris", class(out))        
  return(out)
}

############################---
# S3 method for binris class ----
#############################---

#' @rdname sbinaryRis
#' @method terms binris
#' @export
terms.binris <- function(x, ...){
  formula(x$formula)
}

#' @rdname sbinaryRis
#' @method estfun binris
#' @importFrom sandwich estfun
#' @export estfun.binris
estfun.binris <- function(x, ...){
  class(x) <- c("maxLik", "maxim")
  estfun(x, ...)
}

#' @rdname sbinaryRis
#' @method bread binris
#' @importFrom sandwich bread
#' @export bread.binris
bread.binris <- function(x, ...){
  class(x) <- c("maxLik", "maxim")
  bread(x, ...)
}

#' @rdname sbinaryRis
#' @import stats
#' @export
df.residual.binris <- function(object, ...){
  return(nrow(object$gradientObs) - length(coef(object)))
}

#' @rdname sbinaryRis
#' @method vcov binris
#' @import stats
#' @export 
vcov.binris <- function(object, ...){
  class(object) <- c("maxLik", "maxim")
  vcov(object, ...)
}

#' @rdname sbinaryRis
#' @export
coef.binris <- function(object, ...){
  class(object) <- c("maxLik", "maxim")
  coef(object, ...)
}

#' @rdname sbinaryRis
#' @export 
logLik.binris <- function(object, ...){
  structure(object$maximum, df = length(coef(object)), nobs = nrow(object$gradientObs), class = "logLik")
}


#' @rdname sbinaryRis
#' @method print binris
#' @import stats
#' @export 
print.binris <- function(x, ...){
  cat("Simulated Maximum Likelihood Estimation\n")
  cat(maximType(x), ", ", nIter(x), " iterations\n", sep = "")
  cat("Return code ", returnCode(x), ": ", returnMessage(x), 
      "\n", sep = "")
  if (!is.null(x$estimate)) {
    cat("Log-Likelihood:", x$maximum)
    cat(" (", sum(activePar(x)), " free parameter(s))\n", 
        sep = "")
    cat("Estimate(s):", x$estimate, "\n")
  }
}

#' @rdname sbinaryRis
#' @method summary binris
#' @import stats
#' @importFrom miscTools stdEr
#' @export
summary.binris <- function(object, eigentol = 1e-12, ...){
  result    <- object$maxim
  nParam    <- length(coef(object))
  activePar <- activePar(object)
  if ((object$code < 100) & !is.null(coef(object))) {
    t <- coef(object)/stdEr(object, eigentol = eigentol)
    p <- 2 * pnorm(-abs(t))
    t[!activePar(object)] <- NA
    p[!activePar(object)] <- NA
    results <- cbind(Estimate = coef(object), `Std. error` = stdEr(object, 
                                                                   eigentol = eigentol), 
                     `z value` = t, `Pr(> z)` = p)
  }
  else {
    results <- NULL
  }
  summary <- list(maximType = object$type, iterations = object$iterations, 
                  returnCode = object$code, returnMessage = object$message, 
                  loglik = object$maximum, estimate = results, fixed = !activePar, 
                  NActivePar = sum(activePar), constraints = object$constraints, 
                  model = object$model, varcov = object$varcov, R = object$R, N = nrow(object$X), approximation = object$approximation)
  class(summary) <- "summary.binris"
  summary
}

#' @rdname sbinaryRis
#' @method print summary.binris
#' @import stats
#' @export
print.summary.binris <- function(x, digits = max(3, getOption("digits") - 2),
                                 ...){
  cat("--------------------------------------------\n")
  cat("Simulated ML of Spatial Probit Model using RIS \n")
  cat("Model: ", x$model, ", using R: ", x$R, "\n", sep = "")
  cat("Cholesky applied on ", x$varcov, " matrix\n", sep = "")
  cat("Approximation used?: ", x$approximation, "\n", sep = "")
  cat("Sample size: ", x$N, "\n", sep = "")
  cat(maximType(x), ", ", nIter(x), " iterations\n", sep = "")
  cat("Return code ", returnCode(x), ": ", returnMessage(x), 
      "\n", sep = "")
  if (!is.null(x$estimate)) {
    cat("Log-Likelihood:", x$loglik, "\n")
    cat(x$NActivePar, " free parameters\n")
    cat("Estimates:\n")
    printCoefmat(x$estimate, digits = digits)
  }
  if (!is.null(x$constraints)) {
    cat("\nWarning: constrained likelihood estimation.", 
        "Inference is probably wrong\n")
    cat("Constrained optimization based on", x$constraints$type, 
        "\n")
    if (!is.null(x$constraints$code)) 
      cat("Return code:", x$constraints$code, "\n")
    if (!is.null(x$constraints$message)) 
      cat(x$constraints$message, "\n")
    cat(x$constraints$outer.iterations, " outer iterations, barrier value", 
        x$constraints$barrier.value, "\n")
  }
  cat("--------------------------------------------\n")
}


#' Get Model Summaries for use with "mtable" for objects of class binris
#' 
#' A generic function to collect coefficients and summary statistics from a \code{binris} object. It is used in \code{mtable}
#' 
#' @param obj a \code{binris} object,
#' @param alpha level of the confidence intervals,
#' @param ... further arguments,
#' 
#' @details For more details see package \pkg{memisc}.
#' @return A list with an array with coefficient estimates and a vector containing the model summary statistics. 
#' @importFrom memisc getSummary
#' @method getSummary binris
#' @export 
getSummary.binris <- function(obj, alpha = 0.05, ...){
  smry <- summary(obj)$estimate
  coef <- smry
  lower <- coef[, 1] - coef[, 2] * qnorm(alpha/2)
  upper <- coef[, 1] + coef[, 2] * qnorm(alpha/2)
  coef <- cbind(coef, lower, upper)
  colnames(coef) <- c("est", "se", "stat", "p", "lwr", "upr")
  N <-  nrow(obj$gradientObs)
  sumstat <- c(logLik = logLik(obj), deviance = NA, AIC = AIC(obj), BIC = BIC(obj), N = N, 
               LR = NA, df = NA, p = NA, Aldrich.Nelson = NA, McFadden = NA, Cox.Snell = NA,
               Nagelkerke = NA)
  list(coef = coef, sumstat = sumstat, contrasts = obj$contrasts,
       xlevels = NULL, call = obj$call)
}

############################---
# Additional functions ----
#############################---

# GHK algorithm when Sigma = Sigma_v
ghk_l <- function(s, Sigma, R){
  #R:number of simulations
  n   <- ncol(Sigma)
  eta <- matrix(NA, nrow = n, ncol = R)
  p   <- matrix(NA, nrow = n, ncol = R)
  C   <- t(chol(Sigma))
  # First observation
  nu_1     <-  s[1] / C[1, 1]
  zeta_1   <- m.draws(R)
  eta[1, ] <- qnorm(zeta_1 * pnorm(nu_1))
  p[1, ]   <- pnorm(eta[1, ])
  for(j in 2:n){
    nu       <- (s[j] - C[j, 1:(j - 1)] %*% eta[1:(j - 1), ]) / C[j, j]
    zeta_j   <- m.draws(R) 
    eta[j, ] <- qnorm(zeta_j * pnorm(nu))
    p[j, ]   <- pnorm(nu)
  }
  return(p)
}

# GHK algorithm when Sigma = Sigma_v^{-1}
ghk_u <- function(s, Sigma, R){
  #R:number of simulations
  n        <- ncol(Sigma)
  eta      <- matrix(NA, nrow = n, ncol = R)
  p        <- matrix(NA, nrow = n, ncol = R)
  C        <- chol(Sigma) #Upper triangular 
  #B        <- solve(C)
  B        <- backsolve(C, diag(ncol(C)), upper.tri =  TRUE)
  # Last observation
  nu_n     <-  s[n] / B[n, n]
  zeta_n   <- m.draws(R)
  eta[n, ] <- qnorm(zeta_n * pnorm(nu_n))
  p[n, ]   <- pnorm(rep(nu_n, R))
  # Loop for j 
  for(j in (n - 1):1){
    nu       <- (s[j] - crossprod(B[j, (j + 1):n], eta[(j + 1):n, ])) / B[j, j]
    p[j, ]   <- pnorm(nu)
    if (j > 1) {
      zeta_j   <- m.draws(R)
      eta[j, ] <- qnorm(zeta_j * pnorm(nu))
    }
  }
  return(p)
}


# Antithetical sampling
m.draws <- function(R){
  u <- halton(length = (R/2))
  u <- c(u, 1 - u)
  return(u)
}

# Halton pseudo draws
halton <- function(prime = 3, length = 100, drop = 10){
  halt <- 0
  t <- 0
  while (length(halt) < length + drop) {
    t <- t + 1
    halt <- c(halt, rep(halt, prime - 1) + rep(seq(1, prime - 1, 1) / prime ^ t, each = length(halt)))
  }
  halt[(drop + 1):(length + drop)]
}

# Log-likelihood function for SEM model 
lls_sem <- function(theta, y, X, W, R,
                    varcov = c("sigma", "invsigma"), 
                    Qneg =  TRUE, 
                    approximation = FALSE, 
                    pw = 5){
  varcov   <- match.arg(varcov)
  K        <- ncol(X)
  N        <- nrow(X)
  beta     <- theta[1:K]
  rho      <- theta[K + 1]
  I        <- sparseMatrix(i = 1:N, j = 1:N, x = 1)
  A        <- I - rho * W
  Xb       <- crossprod(t(X), beta)
  if (Qneg){
    Q      <- sparseMatrix(i = 1:N, j = 1:N, x = as.vector(1 - 2*y))
    s      <- - crossprod(t(Q), Xb)
  } else {
    Q      <- sparseMatrix(i = 1:N, j = 1:N, x = as.vector(2*y - 1))
    s      <- crossprod(t(Q), Xb)
  }
  if (varcov == "sigma"){
    A_i      <- if(approximation) app_W(W, rho, pw) else solve(A)
    Sigma_u  <- tcrossprod(A_i)
    Sigma_v  <- crossprod(t(Q), tcrossprod(Sigma_u, Q))
    pir      <- ghk_l(s = s, Sigma = Sigma_v, R = R)
  } else {
    Q_i       <- sparseMatrix(i = 1:N, j = 1:N, x =  1 / diag(Q))
    Sigma_vi  <- tcrossprod(crossprod(Q_i, tcrossprod(A)), t(Q_i))
    pir       <- ghk_u(s = s, Sigma = Sigma_vi, R = R)
  }
  pi_s     <- apply(pmax(pir, .Machine$double.eps), 1, mean)
  ll       <- log(pi_s)
  return(ll)
}

# Log-likelihood function for SAR model 
lls_sar <- function(theta, y, X, W, R,
                    varcov = c("sigma", "invsigma"), 
                    Qneg =  TRUE, 
                    approximation = FALSE, 
                    pw = 5){
  varcov   <- match.arg(varcov)
  K        <- ncol(X)
  N        <- nrow(X)
  beta     <- theta[1:K]
  lambda   <- theta[K + 1]
  I        <- sparseMatrix(i = 1:N, j = 1:N, x = 1)
  A        <- I - lambda * W
  Xb       <- crossprod(t(X), beta)
  A_i      <- if(approximation) app_W(W, lambda, pw) else solve(A)
  if (Qneg){
    Q      <- sparseMatrix(i = 1:N, j = 1:N, x = as.vector(1 - 2*y))
    QA_i   <- crossprod(t(Q), A_i)
    s      <- - crossprod(t(QA_i), Xb)
  } else {
    Q      <- sparseMatrix(i = 1:N, j = 1:N, x = as.vector(2*y - 1))
    QA_i   <- crossprod(t(Q), A_i)
    s      <- crossprod(t(QA_i), Xb)
  }
  if (varcov == "sigma"){
    Sigma_u  <- tcrossprod(A_i)
    Sigma_v  <- crossprod(t(Q), tcrossprod(Sigma_u, Q))
    pir      <- ghk_l(s = s, Sigma = Sigma_v, R = R)
  } else {
    Q_i       <- sparseMatrix(i = 1:N, j = 1:N, x =  1 / diag(Q))
    Sigma_vi  <- tcrossprod(crossprod(Q_i, tcrossprod(A)), t(Q_i))
    pir       <- ghk_u(s = s, Sigma = Sigma_vi, R = R)
  }
  pi_s     <- apply(pmax(pir, .Machine$double.eps), 1, mean)
  ll       <- log(pi_s)
  return(ll)
}



