rm(list=ls())
library(haven)
library(boot)
library(Matrix)

MROZ <- read_dta("~/Downloads/MROZ.dta")

y <- MROZ$inlf
x <- cbind(1,MROZ$nwifeinc,MROZ$educ,MROZ$exper,MROZ$expersq,  MROZ$age, MROZ$kidslt6,MROZ$kidsge6)

###################
###LOGIT MODEL#####
###################

# n0 = 1
# z <- log(( y + 0.5)/(n0 - y + 0.5))
# beta0 <- coefficients(lm(z ~x-1))
# betan <- rep(0.1,ncol(x))
#i = 1
#max.it = 1
#IWLS
#while(sum(abs(beta0-betan)) >= 1e-8 || max.it == 500){

  if(i==1) eta <- x %*% beta0
  else eta <- x %*% betan
  beta0 <- betan
muh <- inv.logit(as.numeric(eta) )
zs <- eta + (y - muh)/(muh*(n0-muh)) *n0
wi <- muh *(n0 - muh)/n0
W <- Diagonal(nrow(x))
diag(W) <- wi
betan <- solve(crossprod(x,W%*%x)) %*% crossprod(x, W%*%zs)
ll = -sum(-log(1+exp(x%*%betan))) - sum(y* x%*%betan)
#ll = (sum(y*  log(exp(x %*% betan)/( 1+ exp(x %*% betan))) + sum((1-y)*log(1-(exp(x %*% betan)/( 1+ exp(x %*% betan)))))))
print(ll)
i = i +1
}

# logit.17.1 <- glm(inlf ~ nwifeinc + educ + exper + expersq + age +
#                     kidslt6 + kidsge6,
#                   data = MROZ,
#                   family = binomial(link = "logit"))

#log.likelihood
loglik.logit <- function(par, y = NULL,x){
  ll = -sum(-log(1+exp(x%*%par))) - sum(y* x%*%par)
  #print(ll)
}
opt <- nlminb(coefficients(lm(y~as.matrix(x)-1)),loglik.logit, y = y, x = as.matrix(x))
opt$objective
opt$par
opt

############################
#### PROBIT MODEL######
############################



# probit.17.1 <- glm(inlf ~ nwifeinc + educ + exper + expersq + age +
#                     kidslt6 + kidsge6,
#                   data = MROZ,
#                   family = binomial(link = "probit"))

loglik.probit <- function(par, y = NULL,x){
  ll = -sum(log(pnorm((2*y-1)*x%*%par)))
  print(ll)
}
opt <- nlminb(coefficients(lm(y~X-1)),loglik.probit, y = y, x = X)
opt$objective
opt$par
