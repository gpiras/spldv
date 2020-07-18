#### Replicate estimation of spatial probit models #####

# Binary GMM: this function estimates a GMM-binary model
binary_gmm <- function(formula, data, subset, na.action, 
                       link = c("probit", "logit"),
                       winitial = c("optimal", "identity"),
                       wmatrix  = c("robust", "iid"), 
                       type = c("onestep", "twostep"), 
                       vce = c("robust", "unadjusted")){
  require("Formula")
  winitial <- match.arg(winitial)
  type     <- match.arg(type)
  vce      <- match.arg(vce)
  link     <- match.arg(link)
  wmatrix  <- match.arg(wmatrix)
  
  #Process the call
  callT <- match.call(expand.dots = TRUE)
  callF <- match.call(expand.dots = FALSE) #... elements of the matched call are remove
  #callF: keep only the arguments which should go into the model frame
  
  #Set up the model frame using Formula method
  formula <- callF$formula <- Formula(formula)
  nframe <- length(sys.calls())
  
  #Model frame
  mf <- callT
  m <- match(c("formula", "data", "subset", "na.action", "weights"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$formula <- formula
  mf[[1L]] <- as.name("model.frame")
  
  #Get model frame after change formula
  mf <- eval(mf, sys.frame(which = nframe))
  
  #Extract variables and parameters
  y <- model.response(mf)
  X <- model.matrix(formula, data = mf, rhs = 1) # Variables
  Z <- model.matrix(formula, data = mf, rhs = 2) # Instruments
  L <- ncol(Z)
  K <- ncol(X)
  N <- nrow(X)
  
  #Warning for underspecified model
  if (L < K) stop("Underspecified model")
  
  #Starting values for optimization of GMM
  start <- rep(0, K)
  names(start) <- colnames(X)
  
  # Initial W matrix
  if (winitial == "optimal") {
    W <- solve(crossprod(Z)/N)
  }  else {
    W <- diag(1, nrow = ncol(Z), ncol = ncol(Z))
  }

  # Optimization for nonlinear GMM estimation (Note tha BFGS can take more iterations)
  # Todo: improve arguments on optimization
  opt <- optim(fn = J, 
               gr = Grad,  
               par = start,
                y = y, 
                X = X, 
                Z = Z, 
                W = W, 
                link = link, 
               method = "BFGS", 
               control = list(trace = 1, maxit = 10000))
  
  # Get results under one step
  b_hat <- opt$par
  if (link == "logit") F <- plogis else F <- pnorm
  u_hat <- y - F(crossprod(t(X), b_hat))
  S     <- makeS(Z = Z, ehat = u_hat)
  #print(opt)
  
  # Proceed with second step if requested 
  if (type == "twostep") {
    W <- if (wmatrix == "robust") solve(S) else  solve(drop(crossprod(u_hat) / N) * crossprod(Z) / N)
    opt <- optim(fn = J, 
                 gr = Grad,  
                 par = b_hat,
                 y = y, 
                 X = X, 
                 Z = Z, 
                 W = W, 
                 link = link, 
                 method = "BFGS", 
                 control = list(trace = 1, maxit = 10000))
    b_hat <- opt$par
    u_hat <- y - F(crossprod(t(X), b_hat))
    S     <- makeS(Z = Z, ehat = u_hat)
  }

  # Generate VCOV
  D <- D_binary(b_hat = b_hat, X = X, link = link)
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


# Return the objective function to be minimized for GMM under binary model
J <- function(start, y, X, Z, W, link){
  # This function returns  J =  (1/n u'Z)W(1/nZ'u) for binary models: a scalar
  if (link == "logit") F <- plogis else F <- pnorm
  xb <- crossprod(t(X), start)
  u <- y - F(xb)
  n <- nrow(X)
  g <- (t(Z) %*% u) / n
  #g <- (t(Z) %*% u)
  J <-  t(g) %*% W %*% g
  #print(J)
  return(J)
}

# This function retrieves the Jacobian: d J / d theta
Grad <- function(start, y, X, Z, W, link){
  # Retrieve Jacobian: K x 1
  if (link == "logit") F <- plogis else F <- pnorm
  xb <- crossprod(t(X), start)
  u <- y - F(xb)
  n <- nrow(X)
  g <- (t(Z) %*% u) / n
  #g <- (t(Z) %*% u)
  D <- D_binary(b_hat = start, X = X, link = link)
  #Grad <- 2 * n *  t(t(Z) %*% D/n) %*% W %*% g
  Grad <- 2 * t(t(Z) %*% D/n) %*% W %*% g
  return(Grad) 
}

# This function retrieves d u/ dtheta of dimension n x K
D_binary <- function(b_hat, X, link = link){
  # Returns a N x K matrix G = d u/ dtheta for binary models
  if (link == "logit") f <- dlogis else f <- dnorm
  xb <- crossprod(t(X), b_hat)
  D <- -1 * drop(f(xb)) * X
  return(D)
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



# Check 1: Check with stata
library("haven")
auto <- read_dta("~/Dropbox/work/Spatial/auto.dta")

# Check onestep procedure, W = optimal, VCE = unajusted
gmm_check1 <- binary_gmm(foreign ~ gear_ratio + length + headroom | gear_ratio + length + headroom, 
                        data = auto, 
                        link = "probit", 
                        winitial = "optimal", 
                        vce = "unadjusted")

# Check onestep procedure, W = optimal, VCE = Robust
gmm_check2 <- binary_gmm(foreign ~ gear_ratio + length + headroom | gear_ratio + length + headroom, 
                        data = auto, 
                        link = "probit", 
                        winitial = "optimal", 
                        vce = "robust")

# Check twostep procedure, W = optimal, VCE = unadjusted, Wmatrix = robust
gmm_check3 <- binary_gmm(foreign ~ gear_ratio + length + headroom | gear_ratio + length + headroom, 
                        data = auto, 
                        link = "probit", 
                        winitial = "optimal", 
                        vce = "unadjusted",
                        type = "twostep")

# Check twostep procedure, W = optimal, VCE = robust, Wmatrix = robust
gmm_check4 <- binary_gmm(foreign ~ gear_ratio + length + headroom | gear_ratio + length + headroom, 
                         data = auto, 
                         link = "probit", 
                         winitial = "optimal", 
                         vce = "robust",
                         type = "twostep")

# Check twostep procedure, W = optimal, VCE = unadjusted, Wmatrix = iid
gmm_check5 <- binary_gmm(foreign ~ gear_ratio + length + headroom | gear_ratio + length + headroom, 
                         data = auto, 
                         link = "probit", 
                         winitial = "optimal", 
                         wmatrix = "iid", 
                         vce = "unadjusted",
                         type = "twostep")

# Check twostep procedure, W = optimal, VCE = robust, Wmatrix = iid
gmm_check6 <- binary_gmm(foreign ~ gear_ratio + length + headroom | gear_ratio + length + headroom, 
                         data = auto, 
                         link = "probit", 
                         winitial = "optimal", 
                         wmatrix = "iid", 
                         vce = "robust",
                         type = "twostep")

# Check 2
N <- 10000
set.seed(1234)
x1 <- rnorm(N)
x2 <- rnorm(N)
y_start <- 1 - 2*x1 + 1 * x2 + rnorm(N)
y <- as.numeric(y_start > 0)
data <- as.data.frame(cbind(y, x1, x2))


gmm_check <- binary_gmm(y ~ x1 + x2 | x1 + x2, 
                        data = data, 
                        link = "probit", 
                        winitial = "optimal", 
                        vce = "unadjusted")