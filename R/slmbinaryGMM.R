## This function estimates SLM binary models of the form:
## y^* = (I - rho W)^{-1}XB + I - rho W)^{-1} * epsilon
## using standard GMM (not linerized) :()

slmbinaryGMM <- function(formula, data, subset, na.action, 
                         listw = NULL, 
                         instruments = 2, 
                         link = c("logit", "probit"),
                         winitial = c("optimal", "identity"),
                         wmatrix  = c("robust", "iid"), 
                         type = c("onestep", "twostep"), 
                         vce = c("robust", "unadjusted"), 
                         gradient = TRUE, 
                         print.init = TRUE,
                         rho.init = 0, 
                         ...){
  # Obtain arguments 
  require("spdep")
  winitial <- match.arg(winitial)
  type     <- match.arg(type)
  vce      <- match.arg(vce)
  wmatrix  <- match.arg(wmatrix)
  link     <- match.arg(link)
  
  # Model Frame
  callT <- match.call(expand.dots = TRUE)
  mf <- callT
  m  <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  nframe  <- length(sys.calls()) # This is for eval
  
  # Optimization default
  if (is.null(callT$method)) callT$method <- 'bfgs'
  if (is.null(callT$iterlm)) callT$iterlim <- 100000
  callT$finalHessian <- FALSE
  
  # Get variables and Globals
  y <- model.response(mf)
  X <- model.matrix(formula, mf)
  N <- nrow(X)
  K <- ncol(X)
  sn <- length(listw$neighbours)
  if (n != sn) stop("number of spatial units in W is different to the number of data")
  
  ## Generate instruments Z
  Z <- make.inst(X, listw, p =  instruments)
  L <- ncol(Z)
  if (L < K) stop("Underspecified model")
  
  ## Starting values for optimization of GMM
  b_init <- glm(formula, family = binomial(link = link), data = mf)$coef
  theta <- c(b_init, rho.init)
  names(theta) <- c(colnames(X), "rho")
  if (print.init) {
    cat("\nStarting Values:\n")
    print(theta)
  } 
  
  ## Initial W matrix
  if (winitial == "optimal") {
    W <- solve(crossprod(Z) / N)
  }  else {
    W <- diag(1, nrow = ncol(Z), ncol = ncol(Z))
  }
  
  
  # Optimization for nonlinear GMM estimation (Note that BFGS can take more iterations)
  cat("\nFirst step GMM optimization based on", winitial, "initial weight matrix \n")
  require("maxLik")
  opt <- callT
  opt$start <- theta
  m <- match(c('method', 'print.level', 'iterlim',
               'start','tol', 'ftol', 'steptol', 'fixed', 'constraints', 
               'control', 'finalHesian'),
             names(opt), 0L)
  opt <- opt[c(1L, m)]
  opt[[1]] <- as.name('maxLik')
  opt$logLik <- as.name('J_minML')
  opt$gradient <- as.name('gradient')
  opt$link <- as.name('link')
  opt$listw <- as.name('listw')
  opt[c('y', 'X', 'Z', 'W')] <- list(as.name('y'), 
                                     as.name('X'), 
                                     as.name('Z'), 
                                     as.name('W'))
  
  x <- eval(opt, sys.frame(which = nframe))
  b_hat <- coef(x)


  # Proceed with second step if requested 
  if (type == "twostep") {
    S <- makeS(b_hat, y, X, Z, listw, link, wmatrix)
    W <- solve(S) 
    cat("\nSecond step GMM optimization based on", wmatrix, "weight matrix \n")
    opt$start <- b_hat
    opt[c('W')] <- list(as.name('W'))
    x <- eval(opt, sys.frame(which = nframe))
    b_hat <- coef(x)
  }
  
  
  ## Saving results
  out <- structure(
                    list(
                          coefficients = b_hat, 
                          call         = callT, 
                          X            = X, 
                          Z            = Z, 
                          y            = y, 
                          listw        = listw,
                          link         = link, 
                          W            = W, 
                          type         = type, 
                          wmatrix      = wmatrix, 
                          winitial     = winitial, 
                          opt          = x, 
                          vce          = vce
                    ), 
                    class = "bingmm"
                  )
  out
}


make.inst <- function(X, listw, p){
  # if X has a constant, drop it to construct the instruments
  if (colnames(X)[1] == "(Intercept)") {
    Xnc <- matrix(X[, -1], dim(X)[1], dim(X)[2] - 1)
    colnames(Xnc) <- colnames(X)[-1]
  } else {
    Xnc <- X
  }
  WX  <- matrix(NA, dim(X)[1], dim(Xnc)[2]) 
  WWX <- matrix(NA, dim(X)[1], dim(Xnc)[2])
  name.wx  <- c()
  name.wwx <- c()
  if (p == 1) {
    for (i in 1:dim(WX)[2])	{
      name.wx <- c(name.wx, paste("w.", colnames(Xnc)[i], sep = ""))
      WX[, i] <- lag.listw(listw, Xnc[, i])
    }
    Z <- cbind(WX)
    colnames(Z) <- name.wx
  }	
  else{ 
    for (i in 1:dim(WX)[2]) {
      name.wx <- c(name.wx, paste("w.", colnames(Xnc)[i], sep = ""))
      WX[, i] <- lag.listw(listw, Xnc[, i])
    }
    for (i in 1:dim(WWX)[2]) { 
      name.wwx <- c(name.wwx, paste("ww.", colnames(Xnc)[i], sep = ""))
      WWX[, i] <- lag.listw(listw, WX[, i])
    }
    Z  <- cbind(WX, WWX)
    colnames(Z) <- c(name.wx, name.wwx)
  }
  Z <- cbind(X, Z)
  return(Z)
}

vcov.bingmm <- function(obj, ...){
  vce   <- obj$vce
  type  <- obj$type
  link  <- obj$link
  listw <- obj$listw
  theta <- obj$coefficients
  wmatrix <- obj$wmatrix
  y     <- obj$y
  Z     <- obj$Z
  X     <- obj$X
  W     <- obj$W
  N     <- nrow(X)
  D     <- momB_slm(start = theta, y, X, Z, listw, link)$D
  S     <- makeS(theta, y, X, Z, listw, link, wmatrix)
  
  if (vce == "robust") {
    pan <- solve(t(D) %*% Z %*% W %*% t(Z) %*% D)
    queso <- t(D) %*% Z %*% W %*% S %*% W %*% t(Z) %*% D 
    V <- N * pan %*% queso %*% pan
  } else {
    if (type == "onestep") {
      sigma <- crossprod(u_hat) / N
      V <- (1/N) * drop(sigma) * solve((t(D) %*% Z / N) %*% solve(S) %*% (t(Z) %*% D / N))
    } else {
      V <- N * solve(t(D) %*% Z %*% solve(S) %*% t(Z) %*% D)
    }
  }
  colnames(V) <- rownames(V) <- names(theta)
  return(V)
}

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

summary.bingmm <- function(object, ...){
  b <- object$coefficients
  std.err <- sqrt(diag(vcov(object)))
  z <- b / std.err
  p <- 2 * (1 - pnorm(abs(z)))
  CoefTable <- cbind(b, std.err, z, p)
  colnames(CoefTable) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  object$CoefTable    <- CoefTable
  class(object)       <- c("summary.bingmm", "bingmm")
  return(object)
}



print.summary.bingmm <- function(x,
                                 digits = max(3, getOption("digits") - 2),
                                 width = getOption("width"),
                                 ...)
{
  cat("        ------------------------------------------------------------\n")
  cat("                      SLM Binary Model by GMM \n")
  cat("        ------------------------------------------------------------\n")
  
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  
  cat("\nCoefficients:\n")
  printCoefmat(x$CoefTable, digit = digits, P.value = TRUE, has.Pvalue = TRUE)
  
  cat(paste("\nSample size:", signif(nrow(x$X), digits)))
  cat(paste("\nJ at the optimum:", signif(-1 * x$opt$maximum, digits)))
  invisible(x)
}



momB_slm <- function(start, 
                     y, 
                     X, 
                     Z,
                     listw, 
                     link){
  # This function generates: 
  #  (1) generalized residuals
  #  (2) the moment conditions 
  #  (3) the gradient dv/dtheta for SLM binary models
  K       <- ncol(X)
  N       <- nrow(X)
  beta    <- start[1:K]
  rho     <- start[K + 1]
  W       <- listw2mat(listw)
  A       <- diag(N) - rho * W
  B       <- solve(A)
  Omega_u <- solve(crossprod(A))
  sigma   <- sqrt(diag(Omega_u))
  ai      <- crossprod(t(B), crossprod(t(X), beta)) /  sigma
  Drho    <- (1 / (2 * sigma)) * diag(Omega_u %*% (A + t(A)) %*% W %*% Omega_u)
  
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
    fa   <- dfun(q*ai)
    Fa   <- pfun(q*ai)
    ffa  <- ddfun(q*ai) 
    v    <- q * (fa/Fa)                                                            # Generalized residuals
    g    <- crossprod(Z, v) / N                                                    # Moment conditions
    der  <- as.vector(q^2 * ((ffa * Fa - fa^2) / (Fa^2)))                          # Common vector of the derivative
    Gb   <- crossprod(t(B), X) *  matrix(der/sigma, nrow = N, ncol = K)            # Gradient of beta: N * K
    grho <- der * (crossprod(t(B), crossprod(t(W), ai)) - (1 / sigma) * ai * Drho) # Gradient of rho
    D    <- cbind(Gb, grho)                                                        # N * K Jacobian matrix
  
  out <- list(g = g, D = D, v = v, u = y - pfun(ai), h = fa/(Fa*(1 - Fa)))
  return(out)
}

# Objective function
J_minML <- function(start, y, X, Z, listw, link, W, gradient){
  getR <- momB_slm(start, y, X, Z, listw, link)
  g <-  getR$g
  J <-  -1 * crossprod(g, crossprod(t(W), g)) # g'Wg 
  
  if (gradient) {
    N <- nrow(X)
    D <- getR$D # N x K matrix
    G <- -2 * t(crossprod(Z, D / N)) %*% crossprod(t(W), g)
    attr(J, "gradient") <- t(G)
  }
  return(J)
}

# Make S matrix: var-cov of moments
makeS <- function(b_hat, y, X, Z, listw, link, wmatrix){
  N <- nrow(X)
  evm <- momB_slm(start = b_hat, y, X, Z, listw, link)
  if (wmatrix == "iid") {
    u_hat <- evm$v
    Shat <- 0
    for (i in 1:N) {
      Shat <- Shat + (u_hat[i] ^ 2 * tcrossprod(Z[i,]))
    }
    Shat <- Shat / N
  } else {
    h       <- evm$h
    u       <- evm$u
    Shat <- 0
    for (i in 1:N) {
      Shat <- Shat + (h[i]^2 * u[i] ^ 2 * tcrossprod(Z[i,]))
    }
    Shat <- Shat / N
  }
  return(Shat)
}

### Overtest
over.test <- function(object, ...){
  Z <- object$Z
  K <- length(coef(object))
  L <- ncol(Z)
  N <- nrow(Z)
  J <- -N * object$opt$maximum 
  out <- list(statistic =  J, df = L - K, p.value =  pchisq(J, L - K, lower.tail =  FALSE))
  return(out)
}
