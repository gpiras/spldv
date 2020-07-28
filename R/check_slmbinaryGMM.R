### Checks for slmbinaryGMM.R

setwd("~/Dropbox/work/Spatial/spldv/R")
source("slmbinaryGMM.R")


#### I. Check SLM probit ----
library("McSpatial")
library("spdep")
set.seed(26051986)
n  <- 225
W  <- cell2nb(sqrt(n), sqrt(n))
rho <- 0.6
IA <- invIrM(W, rho)
alpha   <- -1
beta1   <-  1
beta2   <-  1

x1 <- runif(n)
x2 <- runif(n)
u <- rnorm(n)

ystar <- IA %*% (alpha + beta1 * x1 + beta2 * x2 + u)
y <- as.numeric(ystar > 0)
data <- as.data.frame(cbind(y, x1, x2))


check1 <- slmbinaryGMM(y ~ x1 + x2, 
                       listw = nb2listw(W), 
                       data = data, 
                       instruments = 1, 
                       link = "probit", 
                       winitial = "optimal", 
                       type = "twostep",
                       wmatrix = "robust",
                       gradient = TRUE, 
                       print.level = 2)
summary(check1)


check2 <- slmbinaryGMM(y ~ x1 + x2, 
                       listw = nb2listw(W), 
                       data = data, 
                       instruments = 1, 
                       link = "probit", 
                       winitial = "optimal", 
                       type = "twostep",
                       wmatrix = "iid",
                       gradient = TRUE, 
                       print.level = 2)
summary(check2)


######### Unadjusted matrix
check3 <- slmbinaryGMM(y ~ x1 + x2, 
                       listw = nb2listw(W), 
                       data = data, 
                       instruments = 1, 
                       link = "probit", 
                       winitial = "optimal", 
                       type = "twostep",
                       wmatrix = "robust",
                       gradient = TRUE, 
                       print.level = 2, 
                       vce = "unadjusted")
summary(check3)

check4 <- slmbinaryGMM(y ~ x1 + x2, 
                       listw = nb2listw(W), 
                       data = data, 
                       instruments = 1, 
                       link = "probit", 
                       winitial = "optimal", 
                       type = "twostep",
                       wmatrix = "iid",
                       gradient = TRUE, 
                       print.level = 2, 
                       vce = "unadjusted")
summary(check4)

check5 <- slmbinaryGMM(y ~ x1 + x2, 
                       listw = nb2listw(W), 
                       data = data, 
                       instruments = 1, 
                       link = "probit", 
                       winitial = "identity", 
                       type = "twostep",
                       wmatrix = "robust",
                       gradient = TRUE, 
                       print.level = 2, 
                       vce = "unadjusted")
summary(check5)

check5 <- slmbinaryGMM(y ~ x1 + x2, 
                       listw = nb2listw(W), 
                       data = data, 
                       instruments = 1, 
                       link = "probit", 
                       winitial = "identity", 
                       type = "twostep",
                       wmatrix = "robust",
                       gradient = TRUE, 
                       print.level = 2, 
                       vce = "robust")
summary(check5)


check6 <- slmbinaryGMM(y ~ x1 + x2, 
                       listw = nb2listw(W), 
                       data = data, 
                       instruments = 1, 
                       link = "probit", 
                       winitial = "optimal", 
                       type = "onestep",
                       wmatrix = "robust",
                       gradient = TRUE, 
                       print.level = 2, 
                       vce = "robust")
summary(check6)




## Check McMillien function
check7 <- gmmprobit(y ~ x1 + x2,
                    data = data, 
                    wmat = nb2mat(W))
check7

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

## Bayesian
library("spatialprobit")
X <- cbind(1, x)
check6 <- sarprobit(y ~ X-1, 
                    W = mat2listw(wmat),
                    ndraw = 1000, 
                    burn.in = 200, 
                    thinning = 1, 
                    m = 10)
summary(check6)


#### I. Check SLM logit ----
source("slmbinaryGMM.R")
library("McSpatial")
set.seed(0)
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



