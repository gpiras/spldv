rm(list=ls())
set.seed("19771225")
library(spdep)
library(Matrix)
library(sphet)
mc.count <- 10000
N <- c(49, 100, 225, 400, 900)
lmb <- c(0, 0.2, 0.4, 0.6, 0.8)

MCR1 <- matrix(,nrow = length(lmb), ncol = length(N))
MCR2 <- matrix(,nrow = length(lmb), ncol = length(N))

for(n in N){

  sqn <- sqrt(n)
  Wlist <- nb2listw(cell2nb(sqn,sqn, type = "queen"))
  Ws <- listw2dgCMatrix(Wlist)
  WpW <- crossprod(Ws)
  WW <- Ws %*% Ws
  S <- WpW + WW
  tr <- sum(diag(S))
  sqtr <- sqrt(tr)

  I <- Diagonal(n)
  x <- cbind(1, runif(n,2,6))
  b <- c(1,1)
  xb <- x %*% b
  xpx <- crossprod(x)
  xpxi <- solve(xpx)
  xpxixp <- xpxi %*% t(x)
  P <- x %*% xpxixp
  M <- I - P

  c1 <- -n/((n-ncol(x))*sqtr)
  c2 <- n^2/((n-ncol(x))*(n-ncol(x)+2)*tr)
  EI <- c1 *sum(diag(P %*% Ws))
  TrMWMWp <- sum(diag(M %*% Ws %*% M %*% t(Ws)))
  TrMWMW <- sum(diag(M %*% Ws %*% M %*% Ws))
  TrMW2 <- (sum(diag(M %*% Ws)))^2
  TrTerms <- TrMWMWp + TrMWMW + TrMW2
  VI <- c2 * TrTerms - EI^2

  Dec1 <- vector("list", length = length(lmb))
  Dec2 <- vector("list", length = length(lmb))

  for(lm in lmb ){

    Ilw <- I - lm*Ws
    Ilwi <- solve(Ilw)
    count = 1

    dec1 <- numeric()
    dec2 <- numeric()
    while(count <= mc.count){
      eps <- rnorm(n)
      # The data are generated from an error model
      y <-xb +  Ilwi %*% eps
      # print(var(as.numeric(xb))/var(as.numeric(y)))
      ols <- lm(as.matrix(y) ~ x-1)
      uh <- residuals(ols)
      wuh <- Ws %*% uh
      num <- as.numeric(crossprod(uh, wuh))
      den <- 1/n * crossprod(uh) * sqtr
      MI <- 	as.numeric(num/den)
      SMI <- (MI - EI)/ sqrt(VI)
      dec1 <- c(dec1, ifelse(abs(MI) > 1.96, 1,0)	)
      dec2 <- c(dec2, ifelse(abs(SMI) > 1.96, 1,0))

      count <- count + 1
    }

    Dec1[[which(lmb ==lm )]] <- dec1
    Dec2[[which(lmb ==lm )]] <- dec2

  }
  MCR1[,which(N == n )] <- unlist(lapply(Dec1,mean))
  MCR2[,which(N == n )] <- unlist(lapply(Dec2,mean))
}

rownames(MCR1) <- c("\\lambda = 0", "\\lambda = 0.2", "\\lambda = 0.4", "\\lambda = 0.6", "\\lambda = 0.8")

rownames(MCR2) <- c("\\lambda = 0", "\\lambda = 0.2", "\\lambda = 0.4", "\\lambda = 0.6", "\\lambda = 0.8")

MCR1
MCR2
write.csv(MCR1, "MI_mc1.csv")
write.csv(MCR2, "MI_mc2.csv")




# # # # # # # #
rm(list=ls())
set.seed("19771225")
library(spdep)
mc.count <- 10000
N <- c(49, 100, 225, 400, 900)
lmb <- c(0, 0.2, 0.4, 0.6, 0.8)
rho <- 0.4


MCR <- matrix(,nrow = length(lmb), ncol = length(N))

for(n in N){

  sqn <- sqrt(n)
  Wlist <- nb2listw(cell2nb(sqn,sqn, type = "queen"))
  Ws <- as(Wlist, "CsparseMatrix")
  WpW <- crossprod(Ws)
  WW <- Ws %*% Ws
  S <- WpW + WW
  tr <- sum(diag(S))
  sqtr <- sqrt(tr)

  I <- Diagonal(n)
  x <- cbind(1, runif(n,2,6))
  b <- c(1,1)
  xb <- x %*% b

  Dec <- vector("list", length = length(lmb))

  IrWi <- solve(I - 0.4 * Ws)

  for(lm in lmb ){

    Ilw <- I - lm*Ws
    Ilwi <- solve(Ilw)
    count = 1

    dec <- numeric()

    while(count <= mc.count){
      eps <- rnorm(n)
      y <-  IrWi %*% (xb +  Ilwi%*% eps)
      # print(var(as.numeric(xb))/var(as.numeric(y)))
      ols <- lm(as.matrix(y) ~ x-1)
      uh <- residuals(ols)
      wuh <- Ws %*% uh
      num <- as.numeric(crossprod(uh, wuh))
      den <- 1/n * crossprod(uh) * sqtr
      MI <- 	as.numeric(num/den)

      dec <- c(dec, ifelse(abs(MI) > 1.96, 1,0)	)

      count <- count + 1
    }

    Dec[[which(lmb ==lm )]] <- dec


  }
  MCR[,which(N == n )] <- unlist(lapply(Dec,mean))
}

rownames(MCR) <- c("\\lambda = 0", "\\lambda = 0.2", "\\lambda = 0.4", "\\lambda = 0.6", "\\lambda = 0.8")

MCR

write.csv(MCR, "MI_mc.csv")




# # # Censored regression: TOBIT
install.packages("censReg")
install.packages("AER")
library( "censReg" )
library("AER")
data( "Affairs", package = "AER" )
estResult <- censReg( affairs ~ age + yearsmarried + religiousness + + occupation + rating, data = Affairs )

betas <- coefficients(estResult)[-length(coefficients(estResult))]
sigma<- exp(coefficients(estResult)[length(coefficients(estResult))])
xbhat <- cbind(1, Affairs$age, Affairs$yearsmarried, Affairs$religiousness, Affairs$occupation, Affairs$rating) %*% betas
xbhats <- xbhat / sigma
pdfxb <- pnorm(xbhats )
dfxb <- dnorm(xbhats)
pdxb <- dfxb / pdfxb
fv <- sigma *dfxb *(xbhats + pdxb)
res <- Affairs$affairs - fv



