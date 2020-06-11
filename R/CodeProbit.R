#######################################################
# Lab: Binary Choice Model
# By: Mauricio Sarrias
# Crated on: Sept 15, 2016
# Last modification: June 11, 2020
# Warning: Not completely check
#######################################################



# I. Load the data base ====

library("foreign")
setwd("/home/gpiras/R dev lib/Logit Family/")
binlfp2 <- read.dta("binlfp2.dta", convert.factors = FALSE)
names(binlfp2)


# II. Create model frame using Formula package ====
library("Formula")
F1  <- Formula(lfp ~ k5 + k618 + age + wc + hc)
mf1 <- model.frame(F1, data = binlfp2)
X <- model.matrix(F1, data = binlfp2, rhs = 1)
y <- model.response(mf1)

# III. Function for probit ====
Probit.ml <- function(par, y, X, type = c("logit", "probit")) {
    #Compute LL
    if (type == "logit") {
      F <- plogis
      f <- dlogis
      ff <- function(x) plogis(x) * (1 - plogis(x)) * (1 - 2 * plogis(x))
    } else {
      F <- pnorm
      f <- dnorm
      ff <- function(x) -x * dnorm(x)
    }

    par <- as.matrix(par)
    qi <- 2 * y - 1
    xb <- crossprod(t(X), par)
    ll <- log(F(qi * xb))

    #Compute gradient
    lambdai <- qi * (f(qi * xb) / F(qi * xb))
    attr(ll, "gradient") <- drop(lambdai) * X

    #Hessian
    hi <- drop(((ff(qi * xb) / F(qi * xb)) - ((f(qi * xb) / F(qi * xb))) ^ 2)) * qi ^ 2 * X
    H <- crossprod(hi, X)
    attr(ll, 'hessian') <- H
    ll
  }


# IV. Function for probit ====
par <- lm(y ~ X - 1)$coeff  #starting values
library("maxLik")

## BHHH Method
t1 <- proc.time()
Est1 <- maxLik(logLik = Probit.ml, start = par, y = y, X = X, method = "BHHH")
print(round((proc.time() - t1)[3]/60,9))

#"Newton-Raphson" Method
t1 <- proc.time()
Est2 <- maxLik(logLik = Probit.ml, start = par, y = y, X = X, method = "Newton-Raphson", type = "logit")
print(round((proc.time() - t1)[3]/60,9))

#"BFGS" Method
t1 <- proc.time()
Est3 <- maxLik(logLik = Probit.ml, start = par,y = y, X = X, method = "BFGS")
print(round((proc.time() - t1)[3]/60,9))


attributes(Est1)
Est1$estimate
Est2$estimate
Est3$estimate

  glm(y ~X-1, family = binomial(link = "probit"))
## Get standard errors

G <- Est1$gradientObs
V <- 0
for (i in 1:nrow(X)) {
  V <- V + crossprod(t(G[i,]))
}
Var1 <- solve(V)
sd1  <- sqrt(diag(Var1))

H <- Est2$hessian
Var2 <- -solve(H)
sd2 <- sqrt(diag(Var2))

H2   <- Est3$hessian
Var3 <- -solve(H2)
sd3  <- sqrt(diag(Var3))

Table <- cbind(Est1$estimate, sd1, sd2, sd3)
round(Table,7)
