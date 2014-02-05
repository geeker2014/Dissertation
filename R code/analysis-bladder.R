################################################################################
##
## Analysis of bladder data
##
## Author: Brandon M. Greenwell
##
################################################################################

## Load packages
library(lattice)
library(nlme)
library(lme4)
#library(RLRsim)
library(investr)
library(sandwich)
# library(rjags)

## Bladder data ----------------------------------------------------------------
load("/home/w108bmg/Desktop/Dissertation/Data/bladder.RData")
head(bladder)

## Scatterplot of original data
xyplot(HD ~ volume, groups = subject, data = bladder)
xyplot(HD ~ volume | subject, data = bladder)

## Scatterplot of log transformed data (linear)
bladder2 <- data.frame(logHD = log(bladder$HD), logVolume = log(bladder$volume), 
                       subject = bladder$subject)
xyplot(logHD ~ logVolume, groups = subject, data = bladder2, type = "b",
       cex = 0.5, pch = 16)
xyplot(logHD ~ logVolume | subject, data = bladder2, type = "b",
       cex = 0.5, pch = 16)

## Create grouped data objects
bladder <- na.omit(bladder)
Bladder <- groupedData(HD ~ volume | subject, data = bladder) 
Bladder2 <- groupedData(logHD ~ logVolume | subject, data = bladder2) 

## Analyses of original data ---------------------------------------------------

## Plot of individual regression coefficients
bladder.lmList <- lmList(HD ~ poly(volume, degree = 2) | subject, data = bladder)
bladder.int <- intervals(bladder.lmList)
attributes(bladder.int)$dimnames[[3]] <- c("beta0", "beta1", "beta2")
plot(bladder.int) # suggests beta0 and beta1 are random

## Linear model
bladder.lm <- lm(HD ~ volume, data = bladder)
bladder.lm2 <- lm(HD ~ volume + I(volume^2), data = bladder)
calibrate(bladder.lm, y0 = 75) # inversion limits for linear model
calibrate(bladder.lm, y0 = 75, interval = "Wald") # Wald limits for linear model
invest(bladder.lm2, y0 = 75) # inversion limits for quadratic model
invest(bladder.lm2, y0 = 75, interval = "Wald") # Wald limits for quadratic model

cov2cor(vcov(bladder.lm2))

## Fit mixed-effects models
bladder.lme <- lme(HD ~ volume + I(volume^2), random = ~ volume|subject, data = Bladder)
bladder.lme2 <- update(bladder.lme, random = pdDiag(~volume))
anova(bladder.lme, bladder.lme2) # bladder.lme2 is better
plot(augPred(bladder.lme, length.out = 101, level = 0:1))
plot(augPred(bladder.lme2, length.out = 101, level = 0:1))

## Fixed-effects variance-covariance matrices
round(vcov(bladder.lm2), 5)
round(sandwich(bladder.lm2), 5)
round(vcov(bladder.lme2), 5)
covb <- vcov(bladder.lme2)

## Point estimate --------------------------------------------------------------
y0 <- 70
x0.est <- uniroot(function(x) {
  predict(bladder.lme2, list(volume = x), level = 0) - y0  
}, interval = range(bladder$volume), tol = 1e-10, maxiter = 1000)$root


## Distribution-free calibration interval --------------------------------------
resq <- as.numeric(quantile(resid(bladder.lm2), c(0.025, 0.975)))
c(uniroot(function(x) {
  predict(bladder.lme2, list(volume = x), level = 0) + resq[2] - y0
}, interval = range(bladder$volume), tol = 1e-10, maxiter = 1000)$root,
  uniroot(function(x) {
    predict(bladder.lme2, list(volume = x), level = 0) + resq[1] - y0
  }, interval = range(bladder$volume), tol = 1e-10, maxiter = 1000)$root)


## Parametric bootstrap interval for LMM model ---------------------------------
bladder.lmer <- lmer(HD ~ volume + I(volume^2) + (0+1|subject) + 
                       (0+volume|subject), data = bladder)
sd.y0 <- sqrt(VarCorr(bladder.lmer)[[1]][1] + 
                x0.est^2*VarCorr(bladder.lmer)[[2]][1] + sigma(bladder.lmer)^2)

bootFun <- function(.) {
  covb <- as.matrix(vcov(.))        # (X' V^-1 X)^-1
  beta.boot <- as.numeric(fixef(.)) # fixed effects
  v1.boot <- VarCorr(.)[[1]][1]     # intercept
  v2.boot <- VarCorr(.)[[2]][1]     # linear
  v3.boot <- sigma(.)^2             # residual
  y0.boot <- rnorm(1, mean = 70, sd = sd.y0)
  x0.boot <- (-beta.boot[2]+sqrt(beta.boot[2]^2-4*beta.boot[3]*
                                   (beta.boot[1]-y0.boot)))/(2*beta.boot[3])
  mu0.boot <-  as.numeric(crossprod(beta.boot, c(1, x0.est, x0.est^2)))
  var.y0 <- v1.boot + x0.est^2*v2.boot + v3.boot
  var.mu0 <- t(c(1, x0.est, x0.est^2))%*%covb%*%c(1, x0.est, x0.est^2)
  Z.boot <- (y0.boot-mu0.boot)/sqrt(var.y0+var.mu0)
  c(x0.boot, Z.boot, beta.boot)
}
bootFun0 <- function(.) {
  covb <- as.matrix(vcov(.))        # (X' V^-1 X)^-1
  beta.boot <- as.numeric(fixef(.)) # fixed effects
  v1.boot <- VarCorr(.)[[1]][1]     # intercept
  v2.boot <- VarCorr(.)[[2]][1]     # linear
  v3.boot <- sigma(.)^2             # residual
  y0.boot <- 70
  x0.boot <- (-beta.boot[2]+sqrt(beta.boot[2]^2-4*beta.boot[3]*
                                   (beta.boot[1]-y0.boot)))/(2*beta.boot[3])
  mu0.boot <-  as.numeric(crossprod(beta.boot, c(1, x0.est, x0.est^2)))
  var.y0 <- v1.boot + x0.est^2*v2.boot + v3.boot
  var.mu0 <- t(c(1, x0.est, x0.est^2))%*%covb%*%c(1, x0.est, x0.est^2)
  W.boot <- (y0.boot-mu0.boot)/sqrt(var.y0+var.mu0)
  c(x0.boot, W.boot, beta.boot)
}
set.seed(0101)
bladder.pb <- bootMer2(bladder.lmer, FUN = bootFun, FUN0 = bootFun0, nsim = 9999)
# head(as.data.frame(bladder.pb))
# e1071:::skewness(bladder.pb$t[,1])
beta.boot <- bladder.pb$t[,-1]
covb.boot <- var(beta.boot)

bladder.t <- as.data.frame(bladder.pb)
e1071:::skewness(bladder.pb$t[,1])
apply(bladder.pb$t[,1:2], 2, quantile, prob=c(0.025, 0.975))

Z.boot <- as.numeric(bladder.pb$t)
hist(Z.boot, freq=F, br=30)
curve(dnorm(x), lwd=3, add=T)
quantile(Z.boot, c(0.025, 0.975))
w.025 <- quantile(Z.boot, 0.025)
w.975 <- quantile(Z.boot, 0.975)


## Bootstrap adjusted inversion interval ---------------------------------------
predFun <- function(x) {
  z <- list(x)
  names(z) <- "volume"
  fit <- predict(bladder.lme2, newdata = z, level = 0)
  se.fit <- sqrt(diag(cbind(1, unlist(z), unlist(z)^2) %*% 
                        bladder.lme2$varFix %*% 
                        t(cbind(1, unlist(z), unlist(z)^2))))
  list(fit = fit, se.fit = se.fit)
}
invFun.upper <- function(x) { 
  z <- list(x)
  names(z) <- "volume"
  pred <- predFun(x)
  (y0 - pred$fit)/sqrt((var.y0 + pred$se.fit^2)) - w.975
}
invFun.lower <- function(x) { 
  z <- list(x)
  names(z) <- "volume"
  pred <- predFun(x)
  (y0 - pred$fit)/sqrt((var.y0 + pred$se.fit^2)) - w.025
}
c(uniroot(invFun.upper, interval = c(min(bladder$volume), x0.est), 
          tol = 1e-10, maxiter = 1000)$root, 
  uniroot(invFun.lower, interval = c(x0.est, max(bladder$volume)), 
          tol = 1e-10, maxiter = 1000)$root)
# (57.9428, 138.1559)

boot.ci(v0.pb, type = c("norm", "basic", "perc"))
hist(v0.pb$t, breaks = 50, freq = FALSE, col = "grey", border = "white",
     main = "")
abline(v = x0.est, lwd = 2)
save(bladder.pb, file = "/home/w108bmg/Desktop/Dissertation-knitr/Data/bladder-parboot.RData")


## Wald calibration interval ---------------------------------------------------
beta <- as.numeric(fixef(bladder.lme2))
var.y0 <- getVarCov(bladder.lme2)[1, 1] + 
  x0.est^2*getVarCov(bladder.lme2)[2, 2] + summary(bladder.lme2)$sigma^2
covmat <- diag(4)
covmat[1:3, 1:3] <- covb.boot # covb.boot
covmat[4, 4] <- var.y0
params <- c(beta0 = beta[1], beta1 = beta[2], beta2 = beta[3], Y0 = y0)
gstring <- "(-beta1 + sqrt(beta1^2 - 4*beta2*(beta0-Y0))) / (2*beta2)"
dm <- car:::deltaMethod(params, g = gstring, vcov. = covmat)
x0.est + qnorm(c(0.025, 0.975))*dm$SE


## Inversion interval ----------------------------------------------------------
predFun <- function(x) {
  z <- list(x)
  names(z) <- "volume"
  fit <- predict(bladder.lme2, newdata = z, level = 0)
  se.fit <- sqrt(diag(cbind(1, unlist(z), unlist(z)^2) %*% 
                        bladder.lme2$varFix %*% 
                        t(cbind(1, unlist(z), unlist(z)^2))))
  list(fit = fit, se.fit = se.fit)
}
invFun.bounds <- function(x) { 
  z <- list(x)
  names(z) <- "volume"
  pred <- predFun(x)
  (y0 - pred$fit)^2/(var.y0 + pred$se.fit^2) - qnorm(0.975)^2
}
c(uniroot(invFun.bounds, interval = c(min(bladder$volume), x0.est), 
          tol = 1e-10, maxiter = 1000)$root, 
  uniroot(invFun.bounds, interval = c(x0.est, max(bladder$volume)), 
          tol = 1e-10, maxiter = 1000)$root)


# ## Parametric bootstrap for transformed data -----------------------------------
# bladder.lmer <- lmer(logHD ~ logVolume + (0+1|subject) + (0+logVolume|subject), 
#                      data = bladder2)
# beta <- as.numeric(fixef(bladder.lmer))
# x0.est <- (log(70) - beta[1]) / beta[2]
# sd.y0 <- sqrt(VarCorr(bladder.lmer)[[1]][1] + 
#                 x0.est^2*VarCorr(bladder.lmer)[[2]][1] + sigma(bladder.lmer)^2)
# bootFun <- function(.) {
#   beta.boot <- as.numeric(fixef(.))
#   y0.boot <- rnorm(1, mean = log(70), sd = sd.y0)
#   (y0.boot - beta.boot[1]) / beta.boot[2]
# }
# bootFun0 <- function(.) {
#   beta.boot <- as.numeric(fixef(.))
#   (log(70) - beta.boot[1]) / beta.boot[2]
# }
# set.seed(0202)
# bladder.log.pb <- bootMer2(bladder.lmer, FUN = bootFun, FUN0 = bootFun0, 
#                            nsim = 9999)
# save(bladder.pb, bladder.log.pb, 
#      file = "/home/w108bmg/Desktop/Dissertation-knitr/Data/bladder-parboot.RData")


## Bayesian credible intervals -------------------------------------------------

## Model file for JAGS 
# model.string <- "
# model {
# 
#   ## Likelihood
#   for (i in 1:N) {
#     y[i] ~ dnorm(mu[i], tau.error) 
#     mu[i] <- a0[subject[i]] + a1[subject[i]]*x[i] + mu.a3*x[i]^2
#   }
# 
#   ## Calibration
#   y0 ~ dnorm(a0.hat+a00 + (a1.hat+a01)*x0 + a2.hat*x0^2, tau.e) 
# 
#   ## Priors for random intercepts
#   for (j in 1:m) {
#     a0[j] ~ dnorm(mu.a0, tau.a0)
#     a1[j] ~ dnorm(mu.a1, tau.a1)
#   }
# 
#   ## Priors for fixed effects
#   mu.a0 ~ dnorm(0, 1.0E-6)
#   mu.a1 ~ dnorm(0, 1.0E-6)
#   mu.a2 ~ dnorm(0, 1.0E-6)
# 
#   ## Priors for calibration
#   a00 ~ dnorm(mu.a0, tau.a0)
#   a01 ~ dnorm(mu.a1, tau.a1)
#   x0 ~ dunif(0, 200)
# 
#   ## Prior for variance components
#   tau.a0 ~ dgamma(1.0E-3, 1.0E-3)
#   tau.a1 ~ dgamma(1.0E-3, 1.0E-3)
#   tau.e ~ dgamma(1.0E-3, 1.0E-3)
# 
#   ## Additional calculation
#   sigma.e <- 1/sqrt(tau.e)
#   sigma.a0 <- 1/sqrt(tau.a0)
#   sigma.a1 <- 1/sqrt(tau.a1)
# 
# }
# "

model.string <- "
model {

  ## Likelihood
  for (i in 1:N) {
    y[i] ~ dnorm(a0[subject[i]] + a1[subject[i]]*x[i], tau.e) 
  }
  y0 ~ dnorm(b0 + b1*x0, tau.e) 

  ## Priors for random coefficients
  for (j in 1:m) {
    a0[j] ~ dnorm(mu.a0, tau.a0)
    a1[j] ~ dnorm(mu.a1, tau.a1)
  }

  ## Priors for fixed effects
  mu.a0 ~ dnorm(0, 1.0E-6)
  mu.a1 ~ dnorm(0, 1.0E-6)

  ## Priors for calibration
  b0 ~ dnorm(mu.a0, tau.a0)
  b1 ~ dnorm(mu.a1, tau.a1)
  x0 ~ dgamma(1.0E-3, 1.0E-3)

  ## Prior for variance components
  tau.a0 ~ dgamma(1.0E-3, 1.0E-3)
  tau.a1 ~ dgamma(1.0E-3, 1.0E-3)
  tau.e ~ dgamma(1.0E-3, 1.0E-3)

  ## Additional calculation
  sigma.e <- 1/sqrt(tau.e)
  sigma.a0 <- 1/sqrt(tau.a0)
  sigma.a1 <- 1/sqrt(tau.a1)
  exp.x0 <- exp(x0)

}
"
writeLines(model.string, con = "/home/w108bmg/Desktop/bladder-jags/model.txt")

# Specify data, as a list
data.list = list(N = nrow(bladder), m = length(unique(bladder$subject)), 
                 subject = unclass(bladder$subject), x = log(bladder$volume), 
                 y = log(bladder$HD), y0 = log(70))

## Initialize chains
fit <- lme(log.HD ~ log.volume, data = Bladder.log, random = pdDiag(~log.volume))
beta <- as.numeric(fixef(fit))
x0.est <- calibrate(lm(log.HD ~ log.volume, data = bladder.log), y0 = log(70))
initsFun <- function() { 
  list(mu.a0 = beta[1], mu.a1 = beta[2], x0 = 4.5794, 
       tau.e = 1/summary(fit)$sigma^2, 
       tau.a0 = 1/getVarCov(fit)[1, 1], tau.a1 = 1/getVarCov(fit)[2, 2], 
       a0 = ranef(fit)[, 1] + beta[1], a1 = ranef(fit)[, 2] + beta[2], 
       b0 = beta[1], b1 = beta[1]) 
}
adapt.steps <- 1000    # number of steps to "tune" the samplers
burnin.steps <- 10000  # number of steps to "burn-in" the samplers
n.chains <- 1          # number of chains to run
n.saved.steps <- 50000 # total number of steps in chains to save
thin.steps <- 1        # number of steps to "thin" (1 = keep every step)
n.per.chain <- ceiling((n.saved.steps*thin.steps) / n.chains) # steps per chain

# Run JAGS model
jags.mod <- jags.model("/home/w108bmg/Desktop/bladder-jags/model.txt", 
                       data = data.list, inits = initsFun, n.chains = n.chains, 
                       n.adapt = adapt.steps)
update(jags.mod, n.iter = 10000) # burn-in

params <- c("mu.a0", "mu.a1", "x0", "sigma.e", "sigma.a0", "sigma.a1")
params.coda <- coda.samples(jags.mod, params, n.iter = 10000)
summary(params.coda, quantile = c(0.025, 0.975))

# x0.coda <- coda.samples(jags.mod, "x0", n.iter = 10000)
# x0.samp <- as.matrix(x0.coda)
# 
# beta0.coda <- coda.samples(jags.mod, "beta0", n.iter = 10000)
# beta0.samp <- as.matrix(beta0.coda)
# 
# beta1.coda <- coda.samples(jags.mod, "beta1", n.iter = 10000)
# beta1.samp <- as.matrix(beta1.coda)
# 
# beta2.coda <- coda.samples(jags.mod, "beta2", n.iter = 10000)
# beta2.samp <- as.matrix(beta2.coda)
# 
# params <- cbind(x0.samp, beta0.samp, beta1.samp, beta2.samp)
# apply(params, 2, mean)
# apply(params, 2, var)

params <- c("x0", "beta0", "beta1", "beta2")
params.coda <- coda.samples(jags.mod, params, n.iter = 10000)
plot(params.coda, dens = FALSE)
autocorr.plot(params.coda)

## Analyses of transformed data ------------------------------------------------

bladder2 <- 