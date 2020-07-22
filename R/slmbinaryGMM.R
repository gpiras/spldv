## This function estimates SLM binary models of the form:
## y^* = (I - rho W)^{-1}XB + I - rho W)^{-1} * epsilon
## using standard GMM (not linerized)

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
  
  # Get variables and Globals
  y <- model.response(mf)
  X <- model.matrix(formula, mf)
  N <- nrow(X)
  K <- ncol(X)
  sn <- length(listw$neighbours)
  if (n != sn) stop("number of spatial units in W is different to the number of data")
  
  ## Generate instruments Z
  
  # if X has a constant, drop it to construct the instruments
  if (colnames(X)[1] == "(Intercept)") {
    Xnc <- matrix(X[, -1], dim(X)[1], dim(X)[2] - 1)
  } else {
    Xnc <- X
  }
  WX  <- matrix(NA, dim(X)[1], dim(Xnc)[2]) 
  WWX <- matrix(NA, dim(X)[1], dim(Xnc)[2])
  if (instruments == 1) {
    for (i in 1:dim(WX)[2])	WX[, i] <- lag.listw(listw, Xnc[, i])
    Z <- cbind(WX)
  }	
  else{ 
    for (i in 1:dim(WX)[2])	   WX[, i] <- lag.listw(listw, Xnc[, i])
    for (i in 1:dim(WWX)[2])	WWX[, i] <- lag.listw(listw, WX[, i])
    Z  <- cbind(WX, WWX)
  }
  Z <- cbind(X, Z)
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
               'control'),
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
    u_hat <- momB_slm(start = b_hat, y, X, Z, listw, link)$v
    S <- makeS(Z = Z, ehat = u_hat)
    if (wmatrix == "robust") W <- solve(S) else W <- solve(drop(crossprod(u_hat) / N) * crossprod(Z) / N)
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

vcov.bingmm <- function(obj, ...){
  vce   <- obj$vce
  type  <- obj$type
  link  <- obj$link
  listw <- obj$listw
  theta <- obj$coefficients
  y     <- obj$y
  Z     <- obj$Z
  X     <- obj$X
  W     <- obj$W
  N     <- nrow(X)
  u_hat <- momB_slm(start = theta, y, X, Z, listw, link)$v
  S     <- makeS(Z = Z, ehat = u_hat)
  D     <- momB_slm(start = theta, y, X, Z, listw, link)$D
  
  if (vce == "robust") {
    pan <- solve((t(D) %*% Z / N) %*% W %*% (t(Z) %*% D / N))
    queso <- (t(D) %*% Z / N) %*% W %*% S %*% W %*% (t(Z) %*% D / N)
    V <- (1 / N) * pan %*% queso %*% pan
  } else {
    if (type == "onestep") {
      sigma <- crossprod(u_hat) / N
      V <- (1/N) * drop(sigma) * solve((t(D) %*% Z / N) %*% W %*% (t(Z) %*% D / N))
    } else {
      V <- (1/N) * solve((t(D) %*% Z / N) %*% W %*% (t(Z) %*% D / N))
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
  cat(paste("\nJ at the optimum:", signif(x$opt$maximum, digits)))
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
  ai      <- (B %*% crossprod(t(X), beta)) /  sigma
  Drho    <- Omega_u %*% (W + t(W) - 2*rho * crossprod(W)) %*% Omega_u
  
  if (link == "logit") {
    F    <- plogis
    f    <- dlogis
    v    <- y - F(ai)
    g    <- (t(Z) %*% v) / N
    Gb   <- -1 * drop(f(ai)) * (B %*% X) / matrix(sigma, nrow = N, ncol = K)
    #Grho <- -1 * drop(f(ai)) * (B %*% W %*% ai - (B %*% crossprod(t(X), beta) * diag(Drho) /  (2*sigma^2)))
    Grho <- -1 * drop(f(ai)) * (B %*% W %*% ai - (B %*% crossprod(t(X), beta) * diag(Drho) / (2 * sigma^3)))
    D    <- cbind(Gb, Grho)
  }
  if (link == "probit") {
    F    <- pnorm
    f    <- dnorm
    ff   <- function(x) -x * dnorm(x)
    q    <- 2*y - 1
    fa   <- f(q*ai)
    Fa   <- F(q*ai)
    ffa  <- ff(q*ai) 
    v    <- q * sigma * (fa/Fa)
    g    <- (t(Z) %*% v) / N
    bs   <- as.vector((ffa * Fa - fa*fa) / (Fa^2))
    Gb   <- as.vector(q^2 * sigma * bs) * ((B %*% X) / matrix(sigma, nrow = N, ncol = K))
    grho1 <- (1 / (2 * sigma)) *  diag(Drho) * (fa/Fa)
    grho2 <- as.vector(sigma * bs * q) * (B %*% W %*% ai - (B %*% crossprod(t(X), beta) * diag(Drho) / (2 * sigma^3)))
    grho  <- q*(grho1 + grho2)
    D     <- cbind(Gb, grho)
  }
  
  out <- list(g = g, D = D, v = v)
  return(out)
}

# Objective function
J_minML <- function(start, y, X, Z, listw, link, W, gradient){
  getR <- momB_slm(start, y, X, Z, listw, link)
  g <-  getR$g
  J <-  -1*t(g) %*% W %*% g
  
  if (gradient) {
    N <- nrow(X)
    D <- getR$D # N x K matrix
    G <- -2 * t(t(Z) %*% D/N) %*% W %*% g
    attr(J, "gradient") <- t(G)
  }
  return(J)
}

# Make S matrix: var-cov of moments
makeS <- function(Z, ehat){
  # Create S hat
  N <- nrow(Z)
  Shat <- 0
  for (i in 1:N) {
    Shat <- Shat + (ehat[i] ^ 2 * tcrossprod(Z[i,]))
  }
  return(Shat/N)
}


#### I. Check SLM probit ----
library("McSpatial")
set.seed(2222)
cmap <- readShapePoly(system.file("maps/CookCensusTracts.shp",
                                  package="McSpatial"))
cmap <- cmap[cmap$CHICAGO==1&cmap$CAREA!="O'Hare",]
lmat <- coordinates(cmap)
dnorth <- geodistance(lmat[,1],lmat[,2], -87.627800, 
                      41.881998, dcoor=TRUE)$dnorth
cmap <- cmap[dnorth>0,]
wmat <- makew(cmap)$wmat
n = nrow(wmat)
alpha  <- -1.5
beta   <- 0.5
rho <- 0.6
x <- rnorm(n, 4, 2)
u <- rnorm(n)
A <- solve(diag(n) - rho*wmat)
ystar <- A %*% (alpha + beta * x) + A %*% u
y <- as.numeric(ystar > 0)
data <- as.data.frame(cbind(y, x))

## Check two-step robust VCE with gradient
check1 <- slmbinaryGMM(y ~ x, 
                       listw = mat2listw(wmat), 
                       data = data, 
                       instruments = 1, 
                       link = "probit", 
                       winitial = "optimal", 
                       type = "twostep",
                       wmatrix = "robust",
                       gradient = TRUE, 
                       print.level = 2)
summary(check1)


## Check two-step robust VCE without gradient
check2 <- slmbinaryGMM(y ~ x, 
                       listw = mat2listw(wmat), 
                       data = data, 
                       instruments = 1, 
                       link = "probit", 
                       winitial = "optimal", 
                       type = "twostep",
                       wmatrix = "robust",
                       gradient = FALSE, 
                       print.level = 2)
summary(check2)


## Check McMillien function
check3 <- gmmprobit(y ~ x, 
                       wmat = wmat)
check3

check4 <- spprobit(y ~ x, 
                    wmat = wmat)
check4

## Check with more instruments
check5 <- slmbinaryGMM(y ~ x, 
                       listw = mat2listw(wmat), 
                       data = data, 
                       instruments = 2, 
                       link = "probit", 
                       winitial = "optimal", 
                       type = "twostep",
                       wmatrix = "robust",
                       gradient = TRUE, 
                       print.level = 2)
summary(check5)


#### I. Check SLM logit ----
library("McSpatial")
set.seed(1111)
cmap <- readShapePoly(system.file("maps/CookCensusTracts.shp",
                                  package="McSpatial"))
cmap <- cmap[cmap$CHICAGO==1&cmap$CAREA!="O'Hare",]
lmat <- coordinates(cmap)
dnorth <- geodistance(lmat[,1],lmat[,2], -87.627800, 
                      41.881998, dcoor=TRUE)$dnorth
cmap <- cmap[dnorth>0,]
wmat <- makew(cmap)$wmat
n = nrow(wmat)
alpha  <- -1.5
beta   <- 0.5
rho    <- 0.6
x <- rnorm(n, 4, 2)
u <- rlogis(n)
A <- solve(diag(n) - rho*wmat)
ystar <- A %*% (alpha + beta * x) + A %*% u
y <- as.numeric(ystar > 0)
data <- as.data.frame(cbind(y, x))

## Check two-step robust VCE with gradient
check1 <- slmbinaryGMM(y ~ x, 
                       listw = mat2listw(wmat), 
                       data = data, 
                       instruments = 1, 
                       link = "logit", 
                       winitial = "optimal", 
                       type = "twostep",
                       wmatrix = "robust",
                       gradient = TRUE, 
                       print.level = 2)
summary(check1)


## Check two-step robust VCE without gradient
check2 <- slmbinaryGMM(y ~ x, 
                       listw = mat2listw(wmat), 
                       data = data, 
                       instruments = 1, 
                       link = "logit", 
                       winitial = "optimal", 
                       type = "twostep",
                       wmatrix = "robust",
                       gradient = FALSE, 
                       print.level = 2)
summary(check2)


## Check McMillien function
check3 <- gmmlogit(y ~ x, 
                    wmat = wmat)
check3

check4 <- splogit(y ~ x, 
                   wmat = wmat)
check4

## Check with more instruments
check5 <- slmbinaryGMM(y ~ x, 
                       listw = mat2listw(wmat), 
                       data = data, 
                       instruments = 2, 
                       link = "logit", 
                       winitial = "optimal", 
                       type = "twostep",
                       wmatrix = "robust",
                       gradient = TRUE, 
                       print.level = 2)
summary(check5)



