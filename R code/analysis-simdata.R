################################################################################
##
## Analysis of simdata
##
## Author: Brandon M. Greenwell
##
################################################################################

## Load packages ---------------------------------------------------------------
library(lattice)
library(nlme)
library(lme4)
library(investr)
library(boot)

## Simulated data --------------------------------------------------------------
simRCD <- function(n = 10, m = 20, fixed = c(0, 1), vars = c(0, 0, 0.1)) {
  subject <- rep(1:m, each = n)
  x <- rep(seq(from = 0, to = 1, length = n), times = m) 
  B0 <- rnorm(m, mean = 0, sd = sqrt(vars[1]))
  B1 <- rnorm(m, mean = 0, sd = sqrt(vars[2]))
  y <- rnorm(m*n, mean = fixed[1]+B0[subject] + (fixed[2]+B1[subject])*x, 
             sd = sqrt(vars[3]))
  data.frame(x = x, y = y, subject = factor(subject))
}
set.seed(1234)
simdata <- simRCD(m = 15, n = 30, fixed = c(0, 1), vars = c(0.01, 0, 0.001))

## Scatterplot
xyplot(y ~ x, data = simdata, groups = subject, type = "l")

## Models ----------------------------------------------------------------------
fit.ols <- lm(y ~ x, data = simdata)
fit.lme <- lme(y ~ x, random = ~ 1|subject, data = simdata)
fit.gls <- gls(y ~ x, data = simdata, 
               correlation = corCompSymm(form = ~ 1|subject))

## Compare coefficients
coefs <- rbind(coef(fit.ols), fixef(fit.lme), coef(fit.gls))
rownames(coefs) <- c("OLS", "LMM", "GLS")
colnames(coefs) <- c("beta0", "beta1")
coefs

## Compare variance-covariance matrices (fixed effects)
vcov(fit.ols)
vcov(fit.lme)
vcov(fit.gls)

## Calibration -----------------------------------------------------------------

## Unknowns
y0 <- 0.75
x0.est <- calibrate(fit.ols, y0 = y0)$estimate

## Inversion and Wald intervals for OLS model
calibrate(fit.ols, y0 = y0)
calibrate(fit.ols, y0 = y0, interval = "Wald")

## Inversion and Wald intervals for LMM model
source("/home/w108bmg/Desktop/Dissertation-knitr/R code/invest-lme.R")
invest.lme(fit.lme, y0 = y0, lower = 0, upper = 1.5, tol = 1e-10)
invest.lme(fit.lme, y0 = y0, interval = "Wald")
x0.est + c(qnorm(0.025), qnorm(0.975))*fit.gls$sigma/beta[2]

## Wald interval for GLS model
beta <- as.numeric(coef(fit.gls))
covmat <- diag(3)
covmat[1:2, 1:2] <- vcov(fit.gls)
covmat[3, 3] <- fit.gls$sigma^2
params <- c(beta0 = beta[1], beta1 = beta[2], y0 = y0)
gstring <- "(y0 - beta0) / beta1"
dm <- car:::deltaMethod(params, g = gstring, vcov. = covmat)
rownames(dm) <- ""
dm

## Distribution-free interval for LMM model
resq <- as.numeric(quantile(resid(fit.ols), c(0.025, 0.975)))
c((y0-resq[2] - beta[1])/beta[2], (y0-resq[1] - beta[1])/beta[2])

## Parametric bootstrap interval for LMM model
fit.lmer <- lmer(y ~ x + (1|subject), data = simdata)
bootFun <- function(.) {
  beta <- as.numeric(fixef(.))
  Y0 <- rnorm(1, mean = y0, sd = fit.gls$sigma)
  #Y0 <- y0 + sample(resid(fit.gls), size = 1)
  (Y0 - beta[1]) / beta[2]
}
bootFun0 <- function(.) {
  beta <- as.numeric(fixef(.))
  (y0 - beta[1]) / beta[2]
}
set.seed(1010)
x0.pb <- bootMer2(fit.lmer, FUN = bootFun, FUN0 = bootFun0, nsim = 9999)
boot.ci(x0.pb, type = c("norm", "basic", "perc"))
hist(x0.pb$t, breaks = 50, freq = FALSE, col = "grey", border = "white",
     main = "")
abline(v = x0.est, lwd = 2)
save(x0.pb, file = "/home/w108bmg/Desktop/Dissertation-knitr/Data/simdata-parboot.RData")

## Matrix of different intervals
int1 <- unlist(calibrate(fit.ols, y0 = y0)[c("lower", "upper")])
int2 <- unlist(calibrate(fit.ols, y0 = y0, interval = "Wald")[c("lower", "upper")])
int3 <- invest.lme(fit.lme, y0 = y0, lower = 0, upper = 1.5, tol = 1e-10)[2:3]
int4 <- invest.lme(fit.lme, y0 = y0, interval = "Wald")[2:3]
int5 <- boot.ci(res, type = c("norm", "basic", "perc"))$normal[2:3]
int6 <- boot.ci(res, type = c("norm", "basic", "perc"))$basic[4:5],
int7 <- boot.ci(res, type = c("norm", "basic", "perc"))$percent[4:5]
int8 <- c((y0-resq[2] - beta[1])/beta[2], (y0-resq[1] - beta[1])/beta[2])
cimat <- rbind(c(int1, int1[2] - int1[1]),
               c(int2, int2[2] - int2[1]),
               c(int3, int3[2] - int3[1]), 
               c(int4, int4[2] - int4[1]),
               c(int5, int5[2] - int5[1]),
               c(int6, int6[2] - int6[1]),
               c(int7, int7[2] - int7[1]),
               c(int8, int8[2] - int8[1]))
colnames(cimat) <- c("Lower", "Upper", "Length")
rownames(cimat) <- c("OLS-inversion", "OLS-Wald", "LMM-inversion", "LMM-Wald",
                     "Parboot-norm", "Parboot-basic", "Parboot-perc", 
                     "Dist-free")
cimat
