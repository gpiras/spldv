##### Functions for sbinaryGMM ####


#' @title Estimation of SAR for binary dependent models using GMM
#' 
#' @description Estimation of SAR model for binary dependent variables (either Probit or Logit), using one- or two-step GMM estimator. The type of model supported has the following structure: 
#' 
#' \deqn{
#' y^*= X\beta + WX\gamma + \lambda W y^* + \epsilon = Z\delta + \lambda Wy^{*} + \epsilon
#' }
#' where  \eqn{y = 1} if \eqn{y^*>0} and 0 otherwise; \eqn{\epsilon \sim N(0, 1)} if \code{link = "probit"} or \eqn{\epsilon \sim L(0, \pi^2/3)} if \code{link = "logit"}.
#' 
#' @name sbinaryGMM
#' @param formula a symbolic description of the model of the form \code{y ~ x | wx} where \code{y} is the binary dependent variable, \code{x} are the independent variables. The variables after \code{|} are those variables that enter spatially lagged: \eqn{WX}. The variables in the second part of \code{formula} must also appear in the first part. This rules out situations in which one of the regressors can be specified only in lagged form.
#' @param data the data of class \code{data.frame}.
#' @param listw object. An object of class \code{listw}, \code{matrix}, or \code{Matrix}.  
#' @param nins numerical. Order of instrumental-variable approximation; as default \code{nins = 2}, such that \eqn{H = (Z, WZ, W^2Z)} are used as instruments.  
#' @param link string. The assumption of the distribution of the error term; it can be either \code{link = "probit"} (the default) or \code{link = "logit"}. 
#' @param winitial string. A string indicating the initial moment-weighting matrix \eqn{\Psi}; it can be either \code{winitial = "optimal"} (the default) or \code{winitial = "identity"}. 
#' @param s.matrix string. Only valid of \code{type = "twostep"} is used. This is a string indicating the type of variance-covariance matrix \eqn{\hat{S}} to be used in the second-step procedure; it can be \code{s.matrix = "robust"} (the default) or \code{s.matrix = "iid"}.
#' @param type string. A string indicating whether the one-step (\code{type = "onestep"}), or two-step GMM (\code{type = "twostep"}) should be computed.
#' @param gradient logical. Only for testing procedures. Should the analytic gradient be used in the GMM optimization procedure? \code{TRUE} as default. If \code{FALSE}, then the numerical gradient is used. 
#' @param start if not \code{NULL}, the user must provide a vector of initial parameters for the optimization procedure. When \code{start = NULL}, \code{sbinaryGMM} uses the traditional Probit or Logit estimates as initial values for the parameters, and the correlation between \eqn{y} and \eqn{Wy} as initial value for \eqn{\lambda}.
#' @param cons.opt logical. Should a constrained optimization procedure for \eqn{\lambda} be used? \code{FALSE} as default.  
#' @param approximation logical. If \code{TRUE} then \eqn{(I - \lambda W)^{-1}} is approximated as \eqn{I + \lambda W + \lambda^2 W^2 + \lambda^3 W^3 + ... +\lambda^q W^q}. The default is \code{FALSE}.
#' @param pw numeric. The power used for the approximation \eqn{I + \lambda W + \lambda^2 W^2 + \lambda^3 W^3 + ... +\lambda^q W^q}. The default is 5.
#' @param tol.solve Tolerance for \code{solve()}.
#' @param verbose logical. If \code{TRUE}, the code reports messages and some values during optimization. 
#' @param print.init logical. If \code{TRUE} the initial parameters used in the optimization of the first step are printed. 
#' @param ... additional arguments passed to \code{maxLik}.
#' @param x,object,  an object of class \code{bingmm}
#' @param vce string. A string indicating what kind of standard errors should be computed when using \code{summary}. For the one-step GMM estimator, the options are \code{"robust"} and \code{"ml"}. For the two-step GMM estimator, the options are \code{"robust"}, \code{"efficient"} and \code{"ml"}. The option \code{"vce = ml"} is an exploratory method that evaluates the VC of the RIS estimator using the GMM estimates.
#' @param R numeric. Only valid if \code{vce = "ml"}. It indicates the number of draws used to compute the simulated probability in the RIS estimator.  
#' @param method string. Only valid if \code{vce = "ml"}. It indicates the algorithm used to compute the Hessian matrix of the RIS estimator. The defult is \code{"bhhh"}.  
#' @param digits the number of digits
#' @details 
#' 
#' The data generating process is:
#' 
#' \deqn{
#' y^*= X\beta + WX\gamma + \lambda W y^* + \epsilon = Z\delta + \lambda Wy^{*} + \epsilon
#' }
#' where  \eqn{y = 1} if \eqn{y^*>0} and 0 otherwise; \eqn{\epsilon \sim N(0, 1)} if \code{link = "probit"} or \eqn{\epsilon \sim L(0, \pi^2/3)} if \code{link = "logit"}.. The general GMM
#'   estimator minimizes 
#' \deqn{
#'  J(\theta) = g'(\theta)\hat{\Psi} g(\theta)
#' }
#' where \eqn{\theta = (\beta, \gamma, \lambda)} and 
#' \deqn{
#' g = n^{-1}H'v
#' }
#' where \eqn{v} is the generalized residuals. Let \eqn{Z = (X, WX)}, then the instrument matrix \eqn{H} contains the linearly independent
#' columns of \eqn{H = (Z, WZ, ..., W^qZ)}. The one-step GMM estimator minimizes \eqn{J(\theta)} setting either 
#' \eqn{\hat{\Psi} = I_p} if \code{winitial = "identity"} or \eqn{\hat{\Psi} = (H'H/n)^{-1}} if \code{winitial = "optimal"}. The two-step GMM estimator
#' uses an additional step to achieve higher efficiency by computing the variance-covariance matrix of the moments \eqn{\hat{S}} to weight the sample moments.
#' This matrix is computed using the residuals or generalized residuals from the first-step, which are consistent. This matrix is computed as
#'  \eqn{\hat{S} = n^{-1}\sum_{i = 1}^n h_i(f^2/(F(1 - F)))h_i'} if \code{s.matrix = "robust"} or 
#'   \eqn{\hat{S} = n^{-1}\sum_{i = 1}^n \hat{v}_ih_ih_i'}, where \eqn{\hat{v}} are the first-step generalized residuals. 
#'   
#' @author Mauricio Sarrias and Gianfranco Piras. 
#' @return An object of class ``\code{bingmm}'', a list with elements:
#' \item{coefficients}{the estimated coefficients,}
#' \item{call}{the matched call,}
#' \item{callF}{the full matched call,}  
#' \item{X}{the X matrix, which contains also WX if the second part of the \code{formula} is used, }
#' \item{H}{the H matrix of instruments used,}
#' \item{y}{the dependent variable,}
#' \item{listw}{the spatial weight matrix,}
#' \item{link}{the string indicating the distribution of the error term,}
#' \item{Psi}{the moment-weighting matrix used in the last round,}
#' \item{type}{type of model that was fitted,}
#' \item{s.matrix}{the type of S matrix used in the second round,}
#' \item{winitial}{the moment-weighting matrix used for the first step procedure}
#' \item{opt}{object of class \code{maxLik},}
#' \item{approximation}{a logical value indicating whether approximation was used to compute the inverse matrix,}
#' \item{pw}{the powers for the approximation,}
#' \item{formula}{the formula.}
#' @examples
#' \donttest{
#' # Data set
#' data(oldcol, package = "spdep")
#' 
#' # Create dependent (dummy) variable
#' COL.OLD$CRIMED <- as.numeric(COL.OLD$CRIME > 35)
#' 
#' # Two-step (Probit) GMM estimator
#' ts <- sbinaryGMM(CRIMED ~ INC + HOVAL,
#'                 link = "probit", 
#'                 listw = spdep::nb2listw(COL.nb, style = "W"), 
#'                 data = COL.OLD, 
#'                 type = "twostep",
#'                 verbose = TRUE)
#'                 
#'# Robust standard errors
#'summary(ts)
#'# Efficient standard errors
#'summary(ts, vce = "efficient")
#'
#'# One-step (Probit) GMM estimator 
#'os <- sbinaryGMM(CRIMED ~ INC + HOVAL,
#'                 link = "probit", 
#'                 listw = spdep::nb2listw(COL.nb, style = "W"), 
#'                 data = COL.OLD, 
#'                 type = "onestep",
#'                 verbose = TRUE)
#'summary(os)
#'
#'# One-step (Logit) GMM estimator with identity matrix as initial weight matrix
#'os_l <- sbinaryGMM(CRIMED ~ INC + HOVAL,
#'                   link = "logit", 
#'                   listw = spdep::nb2listw(COL.nb, style = "W"), 
#'                   data = COL.OLD, 
#'                   type = "onestep",
#'                   winitial = "identity", 
#'                   verbose = TRUE)
#'summary(os_l)
#'
#'# Two-step (Probit) GMM estimator with WX
#'ts_wx <- sbinaryGMM(CRIMED ~ INC + HOVAL| INC + HOVAL,
#'                    link = "probit", 
#'                    listw = spdep::nb2listw(COL.nb, style = "W"), 
#'                    data = COL.OLD, 
#'                    type = "twostep",
#'                    verbose = FALSE)
#'summary(ts_wx)
#'
#'# Constrained two-step (Probit) GMM estimator 
#'ts_c <- sbinaryGMM(CRIMED ~ INC + HOVAL,
#'                   link = "probit", 
#'                   listw = spdep::nb2listw(COL.nb, style = "W"), 
#'                   data = COL.OLD, 
#'                   type = "twostep",
#'                   verbose = TRUE, 
#'                   cons.opt = TRUE)
#'summary(ts_c)
#'}
#' @references 
#' 
#' Pinkse, J., & Slade, M. E. (1998). Contracting in space: An application of spatial statistics to discrete-choice models. Journal of Econometrics, 85(1), 125-154.
#' 
#' Fleming, M. M. (2004). Techniques for estimating spatially dependent discrete choice models. In Advances in spatial econometrics (pp. 145-168). Springer, Berlin, Heidelberg.
#' 
#' Klier, T., & McMillen, D. P. (2008). Clustering of auto supplier plants in the United States: generalized method of moments spatial logit for large samples. Journal of Business & Economic Statistics, 26(4), 460-471.
#' 
#' LeSage, J. P., Kelley Pace, R., Lam, N., Campanella, R., & Liu, X. (2011). New Orleans business recovery in the aftermath of Hurricane Katrina. Journal of the Royal Statistical Society: Series A (Statistics in Society), 174(4), 1007-1027.
#' 
#' Piras, G., & Sarrias, M. (2023). One or Two-Step? Evaluating GMM Efficiency for Spatial Binary Probit Models. Journal of choice modelling, 48, 100432. 
#' 
#' Piras, G,. & Sarrias, M. (2023). GMM Estimators for Binary Spatial Models in R. Journal of Statistical Software, 107(8), 1-33. 
#' @seealso \code{\link[spldv]{sbinaryLGMM}}, \code{\link[spldv]{impacts.bingmm}}.
#' @keywords models
#' @rawNamespace import(Matrix,  except = c(cov2cor, toeplitz, update)) 
#' @import stats methods Formula maxLik
#' @importFrom sphet listw2dgCMatrix
#' @export 
sbinaryGMM <- function(formula, 
                       data, 
                       listw = NULL, 
                       nins = 2,                                        # number of instruments
                       link = c("probit", "logit"),                     # probit or logit?
                       winitial = c("optimal", "identity"),             # initial moment-weighing matrix
                       s.matrix = c("robust", "iid"),                   # estimate for second round moment-weighing matrix
                       type = c("onestep", "twostep"),                  # one- and two-step procedure
                       gradient = TRUE,                                 # use the gradient of J for the minimization?
                       start    = NULL,                                 # vector of starting values 
                       cons.opt = FALSE,                                # use constrained optimization for lambda?
                       approximation = FALSE,                           # use inverse approximation?
                       verbose = TRUE,                                  # print messages?
                       print.init = FALSE,                              # print initial values?
                       pw = 5,                                          # default powers for the approximation 
                       tol.solve = .Machine$double.eps,                 # tolarance for solve
                       ...){
  # winitial: Weight matrix of moment for first step. 
  #      1. If optimal, then  Psi=((1/n)Z'Z)^{-1} Eq(13)
  #      2. If identity, then Psi = I Eq(16)
  # s.matrix is the estimated variance matrix for the moments
  #     1. If iid then S = n^1sum_i u_i^2h_ih_i'
  #     2. If robust then as in the paper
  # In this version, we use the notation in paper. 
  
  # Obtain arguments 
  winitial    <- match.arg(winitial)
  type        <- match.arg(type)
  s.matrix    <- match.arg(s.matrix)
  link        <- match.arg(link)
  
  # Spatial weight matrix (W): as CsparseMatrix
  if(!inherits(listw,c("listw", "Matrix", "matrix"))) stop("Neighbourhood list or listw format unknown")
  if(inherits(listw,"listw"))   W    <- sphet::listw2dgCMatrix(listw)	
  if(inherits(listw,"matrix"))  W    <- Matrix(listw)	
  if(inherits(listw,"Matrix"))  W    <- listw	
  
  # Model frame
  callT    <- match.call(expand.dots = TRUE)
  callF    <- match.call(expand.dots = FALSE)
  mf       <- callT
  m        <- match(c("formula", "data"), names(mf), 0L)
  mf       <- mf[c(1L, m)]
  f1       <- Formula(formula)
  if (length(f1)[2L] == 2L) Durbin <- TRUE else Durbin <- FALSE 
  mf$formula <- f1 
  mf[[1L]] <- as.name("model.frame")
  mf       <- eval(mf, parent.frame())
  nframe   <- length(sys.calls())
  
  # Optimization default controls if not added
  if (is.null(callT$method)) callT$method  <- 'bfgs'
  if (is.null(callT$iterlim)) callT$iterlim <- 10000
  # if (is.null(callT$reltol)) callT$reltol  <- 1e-7
  if (is.null(callT$reltol)) callT$reltol  <- 1e-6
  callT$finalHessian <- FALSE # We do not require the Hessian. This speeds the optimization procedure. 
  
  # Get variables and run some checks
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
  
  ## Generate initial instruments 
  Z <- make.instruments(W, x = X, q = nins)
  H <- cbind(X, Z)
  # Just linearly independent columns
  H <- H[, qr(H)$pivot[seq_len(qr(H)$rank)]]
  P <- ncol(H)
  if (P < K) stop("Underspecified model")
  if (any(is.na(H))) stop("NAs in the instruments")
  
  ## Starting values for optimization of GMM
  if (is.null(start)){
    # Initial values from probit model for beta
    #sbinary      <- glm(y ~ as.matrix(X) - 1, family = binomial(link = link), data = mf)
    sbinary      <- glm.fit(as.matrix(X), y, family = binomial(link = link))
    b_init       <- sbinary$coef
    Wy           <- as.numeric(crossprod(t(W), y))
    lambda.init  <- cor(y, Wy)
    theta        <- c(b_init, lambda.init)
    names(theta) <- c(colnames(X), "lambda")
  } else {
    theta <- start
    if (length(start) != length(c(colnames(X), "lambda"))) stop("Incorrect number of intial parameters")
    names(theta) <- c(colnames(X), "lambda")
  }
  if (print.init) {
    cat("\nStarting Values:\n")
    print(theta)
  } 
  
  # lambda_space and constrained optimization if requested. 
  if (cons.opt){
    sym          <- all(W == t(W))
    omega        <- eigen(W, only.values = TRUE, symmetric = sym)
    lambda_space <- if (is.complex(omega$values)) 1 / range(Re(omega$values)) else 1 / range(omega$values)
    # A %*% theta + B >= 0 
    A <- rbind(c(rep(0, K), 1), 
               c(rep(0, K), -1))  # Matrix of restrictions. See help(maxLik)
    B <- cbind(c(-1* (lambda_space[1] + .Machine$double.eps), lambda_space[2] - .Machine$double.eps))
    callT$constraints <- list(ineqA = A, ineqB = B)
  }
  
  ## Initial moment-weighing Psi matrix: it can be (N^{-1}H'H)^{-1} or identity matrix I_P
  #Psi <- if (winitial == "optimal") solve(crossprod(H) / N) else Diagonal(ncol(H))
  Psi <- if (winitial == "optimal") chol2inv(chol(crossprod(H) / N)) else Diagonal(ncol(H))
  # Optimization for nonlinear one-step GMM estimator 
  if (verbose) cat("\nFirst-step GMM optimization based on", winitial, "initial weight matrix \n")
  opt <- callT
  opt$start <- theta
  m <- match(c('method', 'print.level', 'iterlim',
               'start','tol', 'ftol', 'steptol', 'fixed', 'constraints', 
               'control', 'finalHessian', 'reltol'),
             names(opt), 0L)
  opt <- opt[c(1L, m)]
  opt[[1]]      <- as.name('maxLik')
  opt$logLik    <- as.name('J_minML')
  opt$gradient  <- as.name('gradient')
  opt$link      <- as.name('link')
  opt$tol.solve <- as.name('tol.solve')
  opt$listw     <- as.name('W')
  opt$approximation    <- as.name('approximation')
  opt$pw <- as.name('pw')
  opt[c('y', 'X', 'H', 'Psi')] <- list(as.name('y'), 
                                       as.name('X'), 
                                       as.name('H'), 
                                       as.name('Psi'))
  x     <- eval(opt, sys.frame(which = nframe))
  b_hat <- coef(x)

  # Optimization for nonlinear two-step GMM estimator if requested
  if (type == "twostep") {
    # Make S matrix
    S   <- makeS(b_hat, y, X, H, W, link, wmatrix = s.matrix, approximation, pw)
    Psi <- solve(S, tol = tol.solve) 
    
    if (verbose) cat("\nSecond-step GMM optimization using S moment-weighing matrix\n")
    opt$start     <- b_hat
    opt[c('Psi')] <- list(as.name('Psi'))
    x             <- eval(opt, sys.frame(which = nframe))
    b_hat         <- coef(x)
  }
  
  ## Saving results
  out <- structure(
    list(
      coefficients  = b_hat, 
      call          = callT,
      callF         = callF,
      X             = X, 
      H             = H, 
      y             = y, 
      listw         = W,
      link          = link, 
      Psi           = Psi, 
      type          = type, 
      s.matrix      = s.matrix, 
      winitial      = winitial, 
      opt           = x, 
      approximation = approximation,
      pw            = pw, 
      formula       = f1, 
      tol.solve     = tol.solve
    ), 
    class = "bingmm"
  )
  out
}

##### Other functions ----
# Approximate (I-\rho W)^-1
app_W <- function(listw, lambda, pw){
  n           <- dim(listw)[1]
  wi          <- listw
  lambda_loop <- lambda 
  Dn          <- Matrix::Diagonal(n) 
  pW          <- Matrix::Matrix(0, n, n)
  i <- 1
  while (i < pw) {
    wi          <- listw %*% wi 
    lambda_loop <- lambda * lambda_loop
    pW          <- pW + (lambda_loop * wi)
    i           <- i + 1
  }
  pW <- Dn + lambda*listw +  pW
}

# Moment function
momB_slm <- function(start, y, X, H, listw, link, approximation, pw, tol.solve = .Machine$double.eps){
  # This function generates: 
  #  (1) generalized residuals
  #  (2) the moment conditions 
  #  (3) the gradient dv/dtheta for SLM binary models
  K       <- ncol(X)
  N       <- nrow(X)
  beta    <- start[1:K]
  lambda  <- start[K + 1]
  I       <- sparseMatrix(i = 1:N, j = 1:N, x = 1)
  W       <- listw
  A       <- I - lambda * W
  ##FIXME: approximation can be computed in a better way?   
  B       <- if (approximation) app_W(W, lambda, pw)  else solve(A, tol = tol.solve)
  #B      <- qr.solve(A)   
  
  Sigma_u <- tcrossprod(B) #Equation (3)
  sigma   <- sqrt(diag(Sigma_u))
  ai      <- as.vector(crossprod(t(B), crossprod(t(X), beta)) /  sigma)  #Equation (4) (only a_i) and Equation(6)
  
  # F, f and f' for probit or logit
  pfun <- switch(link,
                 "probit" = pnorm,
                 "logit"  = plogis)
  dfun <- switch(link,
                 "probit" = dnorm,
                 "logit"  = dlogis)
  ddfun <- switch(link,
                  "logit"  = function(x) (1 - 2 * pfun(x)) * pfun(x) * (1 - pfun(x)),
                  "probit" = function(x) -x * dnorm(x))  
  q    <- 2*y - 1
  ##FIXME: how to avoid overflow more efficiently?
  #fa   <- dfun(q*ai)
  fa   <- pmax(dfun(q*ai), .Machine$double.eps)
  Fa   <- pmax(pfun(q*ai), .Machine$double.eps)
  #Fa   <- pmin(Fa, 0.99999)# This is to avoid 0 in generalized residuals
  ffa  <- ddfun(q*ai) 
  v    <- q * (fa/Fa)                                                            # Generalized residuals
  g    <- crossprod(H, v) / N                                                    # Moment conditions P * 1
  
  ## FIXME el calculo del gradiente depende dell'argument 
  der  <- as.vector(q^2 * ((ffa * Fa - fa^2) / (Fa^2)))                          # Common vector of the derivative
  Gb   <- crossprod(t(B), X) *  matrix(der/sigma, nrow = N, ncol = K)            # Gradient of beta: N * K Equation (B2) in Appendix
  Drho    <- (1 / (2 * sigma)) * diag(crossprod(t(crossprod(t(Sigma_u), A + t(A))), crossprod(t(W), Sigma_u))) # Equation (B4)
  #Drho    <- (1 / (2 * sigma)) * diag(Sigma_u %*% (A + t(A)) %*% W %*% Sigma_u)
  grho <- der * (crossprod(t(B), crossprod(t(W), ai)) - (1 / sigma) * ai * Drho) # Gradient of rho Equation (B3) (combines with Equation(B4))
  G    <- as.matrix(cbind(Gb, grho))                                             # N* (K + 1)  Jacobian matrix
  
  out <- list(g = g, G = G, v = v, u = y - pfun(as.vector(ai)), vu = (fa^2)/(Fa*(1 - Fa)), Sigma_u = Sigma_u)
  return(out)
}

# Objective function to be minimized
J_minML <- function(start, y, X, H, listw, link, Psi, gradient, approximation, pw, tol.solve){
  getR <- momB_slm(start, y, X, H, listw, link, approximation, pw, tol.solve)
  g <-  getR$g
  J <-  -1 * crossprod(g, crossprod(t(Psi), g)) # g'Psi g 
  
  if (gradient) {
    N  <- nrow(X)
    G  <- getR$G # N x (K + 1) matrix
    Gr <- -2 * t(crossprod(H, G / N)) %*% crossprod(t(Psi), g)
    attr(J, "gradient") <- as.vector(Gr)
  }
  return(J)
}

# Make S matrix: var-cov of moments
makeS <- function(b_hat, y, X, H, listw, link, wmatrix, approximation, pw){
  N <- nrow(X)
  evm <- momB_slm(start = b_hat, y, X, H, listw, link, approximation, pw)
  if (wmatrix == "iid") {
    u_hat <- evm$v # predicted generalized residuals
    Shat <- 0
    for (i in 1:N) {
      Shat <- Shat + (u_hat[i] ^ 2 * tcrossprod(H[i,])) #klier and McMillen (2008)
    }
    Shat <- Shat / N
  } 
  if (wmatrix == "robust"){
    vu       <- evm$vu # This is f^2 / ((1 - F)*F)
    Shat <- 0
    for (i in 1:N) {
      Shat <- Shat + (H[i, ] %*% t(H[i, ] * vu[i]))    ## Equation(15)
    }
    #Shat <- Shat / (N^2)
    Shat <- Shat / N
  }
  return(Shat)
}

### Overtest
#over.test <- function(object, ...){
#  H <- object$H
#  K <- length(coef(object))
#  P <- ncol(H)
#  N <- nrow(H)
#  J <- -N * object$opt$maximum 
#  out <- list(statistic =  J, df = P - K, p.value =  pchisq(J, P - K, lower.tail =  FALSE))
#  return(out)
#}

make.instruments <- function(listw, x, q = 3){
  # This function creates the instruments (WX, ...,W^qX) as in K&P
  W <- listw
  # Drop constant (if any)
  names.x <- colnames(x)
  if (names.x[1] == "(Intercept)") x <- matrix(x[,-1], dim(x)[1], dim(x)[2] - 1)
  names.x <- names.x[which(names.x != "(Intercept)")]
  sq1 <- seq(1, ncol(x) * q, ncol(x))
  sq2 <- seq(ncol(x), ncol(x) * q, ncol(x))
  Zmat <- matrix(NA, nrow = nrow(x), ncol = ncol(x) * q)
  names.ins <- c()
  for (i in 1:q) {
    Zmat[, sq1[i]:sq2[i]] <- as.matrix(W %*% x)
    x <- Zmat[, sq1[i]:sq2[i]]
    names.ins <- c(names.ins, paste(paste(replicate(i, "W"), collapse = ""), names.x, sep = "*"))
  }
  colnames(Zmat) <- names.ins
  return(Zmat)
}


### S3 methods ----

#' @rdname sbinaryGMM
#' @export
coef.bingmm <- function(object, ...){
  object$coefficients
}


#' @rdname sbinaryGMM
#' @method vcov bingmm
#' @export 
vcov.bingmm <- function(object, vce = c("robust", "efficient", "ml"), method = "bhhh", R = 1000, tol.solve = .Machine$double.eps, ...){
  # vce: indicates the estimator for the variance-covariance matrix when using two-step estimator
  #       1- if efficient, then the lowest bound of the vc is used. 
  vce           <- match.arg(vce)
  R             <- R
  method        <- method
  link          <- object$link
  listw         <- object$listw
  type          <- object$type 
  theta         <- object$coefficients
  s.matrix      <- object$s.matrix
  winitial      <- object$winitial
  approximation <- object$approximation
  pw            <- object$pw
  K             <- length(theta)
  y             <- object$y
  H             <- object$H
  X             <- object$X
  Psi           <- object$Psi
  N             <- nrow(X)
  # FIXME:  Si traes el gradiente desde el optimizer no necessitas calcularlo otra vez aca...
  G        <- momB_slm(start = theta, y, X, H, listw, link, approximation, pw)$G # Gradient evaluated at the optimum
  if (vce == "efficient" && type == "onestep") stop("Efficient VCE for one-step estimator not yet implemented")
  if (vce == "efficient"){
    G_bar  <- crossprod(H, G)
    # Similar to Stata, we use the the Psi that was used to compute the final-round estimate 
    V      <- N * solve(crossprod(G_bar, crossprod(t(Psi), G_bar)), tol = tol.solve)
  }
  if (vce == "robust"){
    # Based on Stata: Psi is the weight matrix requested with wmatrix() and it is calculated based on the residuals obtained after the first estimation step. 
    S     <- makeS(theta, y, X, H, listw, link, s.matrix, approximation, pw)
    G_bar <- crossprod(H, G)
    pan   <- solve(crossprod(G_bar, crossprod(t(Psi), G_bar)), tol = tol.solve)
    queso <- crossprod(t(crossprod(G_bar, Psi)), crossprod(t(crossprod(t(S), Psi)), G_bar))
    V     <- N * pan %*% queso %*% pan
  }
  if (vce == "ml"){
    opt <- object$call
    opt$start <- theta
    m <- match(c('start'),
               names(opt), 0L)
    opt$iterlim <- 0
    opt <- opt[c(1L, m)]
    opt[[1]]     <- as.name('maxLik')
    opt$logLik   <- as.name('ll_ris')
    W <- listw
    opt$method <- as.name("method")
    opt$W    <- as.name('W')
    opt$R    <- as.name('R')
    opt[c('y', 'X')] <- list(as.name('y'), as.name('X'))
    x <- eval(opt)
    return(vcov(x))
  }
  colnames(V) <- rownames(V) <- names(theta)
  return(V)
}


#' @rdname sbinaryGMM
#' @method print bingmm
#' @export 
print.bingmm <- function(x, 
                         digits = max(3, getOption("digits") - 3),
                         ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print.default(format(drop(coef(x)), digits = digits), print.gap = 2,
                quote = FALSE)
  cat("\n")
  invisible(x)
}


#' @rdname sbinaryGMM
#' @method summary bingmm
#' @export
summary.bingmm <- function(object, vce = c("robust", "efficient", "ml"), method = "bhhh", R = 1000, tol.solve = .Machine$double.eps, ...){
  vce                 <- match.arg(vce)
  b                   <- object$coefficients
  std.err             <- sqrt(diag(vcov(object, vce = vce, method = method, R = R, tol.solve = tol.solve)))
  z                   <- b / std.err
  p                   <- 2 * (1 - pnorm(abs(z)))
  CoefTable           <- cbind(b, std.err, z, p)
  colnames(CoefTable) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  object$CoefTable    <- CoefTable
  class(object)       <- c("summary.bingmm", "bingmm")
  return(object)
}


#' @rdname sbinaryGMM
#' @method print summary.bingmm
#' @export
print.summary.bingmm <- function(x,
                                 digits = max(5, getOption("digits") - 3),
                                 ...)
{
  cat("        ------------------------------------------------------------\n")
  cat("                      SLM Binary Model by GMM \n")
  cat("        ------------------------------------------------------------\n")
  
  
  cat("\nCall:\n")
  cat(paste(deparse(x$callF), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  
  cat("\nCoefficients:\n")
  printCoefmat(x$CoefTable, digits = digits, P.values = TRUE, has.Pvalue = TRUE)
  
  cat(paste("\nSample size:", signif(nrow(x$X), digits)), "\n")
  #cat(paste("\nJ at the optimum:", signif(-1 * x$opt$maximum, digits)))
  #cat(paste("\nInstruments:", paste(colnames(x$H))))
  invisible(x)
}

#' Get Model Summaries for use with "mtable" for objects of class bingmm
#' 
#' A generic function to collect coefficients and summary statistics from a \code{bingmm} object. It is used in \code{mtable}
#' 
#' @param obj a \code{bingmm} object,
#' @param alpha level of the confidence intervals,
#' @param ... further arguments,
#' 
#' @details For more details see package \pkg{memisc}.
#' @return A list with an array with coefficient estimates and a vector containing the model summary statistics. 
#' @importFrom memisc getSummary
#' @method getSummary bingmm
#' @export 
getSummary.bingmm <- function(obj, alpha = 0.05, ...){
  if (inherits(obj, c("summary.bingmm"))){
    coef <- obj$CoefTable
  } else {
    smry <- summary(obj)
    coef <- smry$Coef
  }
  lower <- coef[, 1] - coef[, 2] * qnorm(alpha/2)
  upper <- coef[, 1] + coef[, 2] * qnorm(alpha/2)
  coef <- cbind(coef, lower, upper)
  colnames(coef) <- c("est", "se", "stat", "p", "lwr", "upr")
  N <-  nrow(obj$X)
  #sumstat <- c(logLik = NA, deviance = NA, AIC = NA, BIC = NA, N = N, 
  #             LR = NA, df = NA, p = NA, Aldrich.Nelson = NA, McFadden = NA, Cox.Snell = NA,
  #             Nagelkerke = NA)
  sumstat <- c(N = N)
  list(coef = coef, sumstat = sumstat, contrasts = obj$contrasts,
       xlevels = NULL, call = obj$call)
}


