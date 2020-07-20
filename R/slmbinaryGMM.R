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
                         print.init = TRUE){
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
  
  #Starting values for optimization of GMM
  b_init <- glm(formula, family = binomial(link = link), data = mf)$coef
  start <- c(b_init, 0)
  names(start) <- c(colnames(X), "rho")
  #if (print.init) cat(paste("\nInitial values for first step:\n", start))
  
  # Initial W matrix
  if (winitial == "optimal") {
    W <- solve(crossprod(Z) / N)
  }  else {
    W <- diag(1, nrow = ncol(Z), ncol = ncol(Z))
  }
  
  gr <- if (gradient) G_min else NULL
  # Optimization for nonlinear GMM estimation (Note that BFGS can take more iterations)
  # Todo: improve arguments on optimization
  
  cat("\nFirst step GMM optimization based on", winitial, "initial weight matrix \n")
  opt <- optim(fn = J_min, 
               gr = gr,   
               par = start,
               y = y, 
               X = X, 
               Z = Z, 
               W = W,
               link = link,
               listw = listw,
               method = "BFGS", 
               control = list(trace = 1, maxit = 10000))
  # Get results under one step
  b_hat <- opt$par
  u_hat <- momB_slm(start = b_hat, y, X, Z, listw, link)$v
  S <- makeS(Z = Z, ehat = u_hat)
  #print(opt)
  
  # Proceed with second step if requested 
  if (type == "twostep") {
    if (wmatrix == "robust") {
      W <- solve(S)
    }
    if (wmatrix == "iid") {
      W <- solve(drop(crossprod(u_hat) / N) * crossprod(Z) / N)
    }
    cat("\nSecond step GMM optimization based on", wmatrix, "weight matrix \n")
    opt <- optim(fn = J_min, 
                 gr = gr,   
                 par = b_hat,
                 y = y, 
                 X = X, 
                 Z = Z, 
                 W = W,
                 link = link,
                 listw = listw,
                 method = "BFGS", 
                 control = list(trace = 1, maxit = 10000))
    b_hat <- opt$par
    u_hat <- momB_slm(start = b_hat, y, X, Z, listw, link)$v
    S     <- makeS(Z = Z, ehat = u_hat)
  }
  
  # Generate VCOV
  D <- momB_slm(start = b_hat, y, X, Z, listw, link)$D
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
  
  se <- sqrt(diag(V))
  t       <- b_hat / se
  pv      <- pnorm(abs(t), lower.tail = FALSE) * 2
  table   <- cbind(b_hat, se, t, pv)
  
  cat(paste("\nFinal value of J:", opt$value, "\n"))
  cat(paste("\nType of GMM:", type, "\n"))
  cat(paste("\nInitial W matrix:", winitial, "\n"))
  cat(paste("\nVCE:", vce, "\n\n"))
  cat(paste("\nEstimates from GMM binary model \n\n"))
  colnames(table) <- c("Estimate", "Std. Error", "t-value", "Pr(>|t|)")
  printCoefmat(table)
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
  ai      <- B %*% crossprod(t(X), beta) /  sigma
  Drho    <- Omega_u %*% (W + t(W) - 2*rho * crossprod(W)) %*% Omega_u
  
  if (link == "logit") {
    F    <- plogis
    f    <- dlogis
    v    <- y - F(ai)
    g    <- (t(Z) %*% v) / N
    Gb   <- -1 * drop(f(ai)) * X / sigma
    #Grho <- -1 * drop(f(ai)) * (solve(A) %*% W %*% ai - (solve(A) %*% crossprod(t(X), beta) * diag(dlamb) /  sigma^2))
    Grho <- -1 * drop(f(ai)) * (B %*% W %*% ai - (B %*% crossprod(t(X), beta) * diag(Drho) / (2 * sigma^3)))
    D    <- cbind(Gb, Grho)
  }
  out <- list(g = g, D = D, v = v)
  return(out)
}

# Objective function
J_min <- function(start, y, X, Z, listw, link, W){
  getR <- momB_slm(start, y, X, Z, listw, link)
  g <-  getR$g
  J <-  t(g) %*% W %*% g
  return(J)
}

# Gradient of objective function
G_min <- function(start, y, X, Z, listw, link, W){
  getR <- momB_slm(start, y, X, Z, listw, link)
  N <- nrow(X)
  g <- getR$g
  D <- getR$D # N x K matrix
  G <- 2 * t(t(Z) %*% D/N) %*% W %*% g
  return(G)
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

## Check
library("McSpatial")
set.seed(9947)
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
lambda <- 0.6
x <- rnorm(n, 4, 2)
u <- rlogis(n, location = 0, scale = 1)
A <- solve(diag(n) - lambda*wmat)
ystar <- A %*% (alpha + beta * x) + A %*% u
y <- as.numeric(ystar > 0)

data <- as.data.frame(cbind(y, x))

## Check two-step robust VCE
slmbinaryGMM(y ~ x, 
             listw = mat2listw(wmat), 
             data = data, 
             instruments = 1, 
             link = "logit", 
             winitial = "optimal", 
             type = "twostep",
             wmatrix = "robust",
             gradient = TRUE)

slmbinaryGMM(y ~ x, 
             listw = mat2listw(wmat), 
             data = data, 
             instruments = 1, 
             link = "logit", 
             winitial = "optimal", 
             type = "twostep",
             wmatrix = "robust",
             gradient = FALSE)

## Check two-step unajusted VCE
slmbinaryGMM(y ~ x, 
             listw = mat2listw(wmat), 
             data = data, 
             instruments = 1, 
             link = "logit", 
             winitial = "optimal", 
             type = "twostep",
             wmatrix = "robust",
             vce = "unadjusted", 
             gradient = TRUE)

## Check with McMillen functions
fit <- splogit(y ~ x,  wmat = wmat)

# As Lagunes: onestep uses W = I
slmbinaryGMM(y ~ x, 
             listw = mat2listw(wmat), 
             data = data, 
             instruments = 1, 
             link = "logit", 
             winitial = "identity", 
             type = "twostep",
             wmatrix = "robust",
             gradient = TRUE)

