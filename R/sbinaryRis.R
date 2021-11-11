### Ris simulator  ----

#' Estimation of SAR for binary models using RIS
#' 
#' @name sbinaryRis
#' @import Matrix maxLik stats spatialreg
#' @export 
sbinaryRis <- function(formula, data, subset, na.action, 
                        listw = NULL, 
                        R = 1000, # Number of draws
                        method = "bhhh", 
                        rho.init = NULL, 
                        print.init = TRUE, 
                        ll = c("ll_ris", "ll_ris2"),
                        eigen = FALSE, 
                        ...){
  # Obtain arguments 
  W  <- listw
  R  <- R
  ll <- match.arg(ll)
  
  # Model Frame
  callT    <- match.call(expand.dots = TRUE)
  mf       <- callT
  m        <- match(c("formula", "data"), names(mf), 0L)
  mf       <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf       <- eval(mf, parent.frame())
  nframe   <- length(sys.calls()) # This is for evaluation
  
  # Get variables and globals
  y <- model.response(mf)
  X <- model.matrix(formula, mf)
  N <- nrow(X)
  K <- ncol(X)
  sn <- nrow(W)
  if (N != sn) stop("Number of spatial units in W is different to the number of data")
  
  ## Starting values
  sbinary <- glm(formula, family = binomial(link = "probit"), data = mf)
  b_init  <- sbinary$coef # Initial values of beta from a standard binary model
  Wy <- as.numeric(crossprod(t(W), y))
  if (is.null(rho.init)) rho.init <- cor(y, Wy)
  theta <- c(b_init, rho.init)
  names(theta) <- c(colnames(X), "rho")
  # if (print.init) {
  #  cat("\nStarting Values:\n")
  # print(theta)
  #} 
  
  # Optimization 
  opt <- callT
  opt$start <- theta
  m <- match(c('method', 'print.level', 'iterlim',
               'start','tol', 'ftol', 'steptol', 'fixed', 'constraints', 
               'control', 'finalHessian'),
             names(opt), 0L)
  opt <- opt[c(1L, m)]
  opt[[1]]     <- as.name('maxLik')
  if (ll == "ll_ris") opt$logLik   <- as.name('ll_ris') else opt$logLik <- as.name('ll_ris2') 
  opt$W    <- as.name('W')
  opt$R    <- as.name('R')
  opt$eigen <- as.name('eigen')
  opt[c('y', 'X')] <- list(as.name('y'), as.name('X'))
  
  x <- eval(opt, sys.frame(which = nframe))
  return(x)
}

ll_ris <- function(start, y, X, W, R, eigen = FALSE){
  K       <- ncol(X)
  N       <- nrow(X)
  beta    <- start[1:K]
  rho     <- start[K + 1]
  I       <- sparseMatrix(i = 1:N, j = 1:N, x = 1)
  W       <- as(W, "CsparseMatrix")
  # A <- I + rhoW + rho^2W^2 + rho^3W^3
  A       <- I - rho * W
  AI      <- solve(A)
  Sigma_u <- tcrossprod(AI)
  Z       <- diag(as.vector(1 - 2*y))
  # Generate upsilon
  V       <- - Z %*% AI %*% X %*% beta # Z(I -rho * W)^{-1}Xb
  # Generate Omega of upsilon
  Omega_V <- Z %*% Sigma_u %*% t(Z)
  # Generate upper-triangular matrix C = A^{-1}, where A = cholesky
  #check <- t(chol(as.matrix(solve(Omega_V)))) %*% chol(as.matrix(solve(Omega_V))) # check
  C        <- solve(chol(as.matrix(solve(Omega_V)))) 
  #C        <- solve(chol(solve(Omega_V))) 
  Eta <- p <- matrix(NA, nrow = N, ncol = R)
  a        <- rep(1/C[1, 1] * V[1], R)
  #Eta[1, ] <- qnorm(draws(R) * pnorm(b[N,]))
  
  #p <- rep(1, R)
  for(j in 1:N){
    if (j > 1) a <- 1/C[j, j] * (V[j] - t(C[j, 1:(j - 1)]) %*% Eta[1:(j - 1), ]) # vector or R * 1
    a <- pmax(pnorm(a), .Machine$double.eps)
    p[j, ] <- a
    if (j < N) Eta[j, ] <- qnorm(draws(R) * a)
  }
  p_sim <- apply(p, 1, mean)
  ll <- log(p_sim)
  return(ll)
}

ll_ris2 <- function(start, y, X, W, R, eigen = FALSE){
  K       <- ncol(X)
  N       <- nrow(X)
  beta    <- start[1:K]
  rho     <- start[K + 1]
  I       <- sparseMatrix(i = 1:N, j = 1:N, x = 1)
  W       <- as(W, "CsparseMatrix")
  A       <- I - rho * W
  AI      <- solve(A)
  Sigma_u <- tcrossprod(AI)
  Z       <- diag(as.vector(1 - 2*y))
  # Generate upsilon
  V       <- - Z %*% AI %*% X %*% beta
  # Generate Omega of upsilon
  Omega_V <- Z %*% Sigma_u %*% t(Z)
  # Generate upper-triangular matrix C = A^{-1}, where A = cholesky
  #check <- t(chol(as.matrix(solve(Omega_V)))) %*% chol(as.matrix(solve(Omega_V))) # check
  C           <- solve(chol(as.matrix(solve(Omega_V)))) 
  Eta <- eta0 <-  matrix(NA, nrow = N, ncol = R)
  eta0[N, ]   <- rep(1/C[N, N] * V[N], R)
  Eta[N, ]    <- qnorm(draws(R) * pnorm(eta0[N, ]))
  
  
  for(j in (N-1):1){
    eta0[j, ] <- 1/C[j, j] * (V[j] - t(C[j, (j + 1):N]) %*% Eta[(j + 1):N, ]) # vector or R * 1
    if (j > 1) Eta[j, ] <- qnorm(draws(R) * pnorm(eta0[j, ]))
  }
  pi  <- pmax(pnorm(eta0), .Machine$double.eps)
  pi_s <- apply(pi, 1, mean)
  ll <- log(pi_s)
  return(ll)
}


# Generate random draw from importance density function
draws <- function(R){
  #q <- runif(R/2)
  u <- halton(length = (R/2))
  u <- c(u, 1-u)
  #qnorm(q * pnorm(upper.bound))
  return(u)
}


halton <- function(prime = 3, length = 100, drop = 10){
  halt <- 0
  t <- 0
  while (length(halt) < length + drop) {
    t <- t + 1
    halt <- c(halt, rep(halt, prime - 1) + rep(seq(1, prime - 1, 1) / prime ^ t, each = length(halt)))
  }
  halt[(drop + 1):(length + drop)]
}

#' Get Model Summaries for use with "mtable" for objects of class maxLik
#' 
#' A generic function to collect coefficients and summary statistics from a \code{maxLik} object. It is used in \code{mtable}
#' 
#' @param obj a \code{maxLik} object,
#' @param alpha level of the confidence intervals,
#' @param ... further arguments,
#' 
#' @details For more details see package \pkg{memisc}.
#' @import stats
#' @importFrom memisc getSummary
#' @export 
getSummary.maxLik <- function(obj, alpha = 0.05, ...){
  smry <- summary(obj)$estimate
  coef <- smry
  lower <- coef[, 1] - coef[, 2] * qnorm(alpha/2)
  upper <- coef[, 1] + coef[, 2] * qnorm(alpha/2)
  coef <- cbind(coef, lower, upper)
  colnames(coef) <- c("est", "se", "stat", "p", "lwr", "upr")
  N <-  nrow(obj$gradientObs)
  sumstat <- c(logLik = NA, deviance = NA, AIC = NA, BIC = NA, N = N, 
               LR = NA, df = NA, p = NA, Aldrich.Nelson = NA, McFadden = NA, Cox.Snell = NA,
               Nagelkerke = NA)
  list(coef = coef, sumstat = sumstat, contrasts = obj$contrasts,
       xlevels = NULL, call = obj$call)
}
