################################################################################
##
## Simulation study for linear calibration with grouped data.
##
## Author: Brandon M. Greenwell
## Date: 2/4/2014
##
## TODO:
##   (1) Should Y0 be simulated ahead of time; for example,
##       Y0 <- rnorm(nsim, 0, var.y0)? Is this even possible?
##
################################################################################

## Load packages ---------------------------------------------------------------
library(lattice) # for plotting grouped data
library(plyr)    # for working with the data
library(lme4)    # for fitting mixed models and bootstrap inference
library(nlme)    # for fitting mixed models
library(boot)    # for calculating bootstrap confidence intervals
source("/home/w108bmg/Desktop/Dissertation/R code/bootMer2.R")
source("/home/w108bmg/Desktop/Dissertation/R code/bootMer2_parallel.R")

## Simulation setup ------------------------------------------------------------

## Parameters
params <- list(
  nsim = 1000,                   # simulation size
  m = 30,                        # number of subjects
  n = 20,                        # sample size per subject
  beta = c(0, 2),                # fixed effecs
  theta = c(0.001, 0.05, 0.001), # variance components
  y0 = c(0, 0.5, 1, 1.5, 2)[4]   # true observed response: 0, 0.5, 1.0, 2.0
)
params$x0 <- (params$y0 - params$beta[1]) / params$beta[2]
params$var.y0 <- params$theta[1] + params$theta[2]*params$x0^2 + 
  params$theta[3]

## Function to simulate data from a random intercept and slope model
simData <- function(n = params$n, m = params$m, beta = params$beta, 
                    theta = params$theta) {
  subject <- rep(1:m, each = n)
  x <- rep(seq(from = 0, to = 1, length = n), times = m) 
  B0 <- rnorm(m, mean = 0, sd = sqrt(theta[1]))
  B1 <- rnorm(m, mean = 0, sd = sqrt(theta[2]))
  y <- rnorm(m*n, mean = beta[1]+B0[subject] + (beta[2]+B1[subject])*x, 
             sd = sqrt(theta[3]))
  data.frame(x = x, y = y, subject = factor(subject))
}

## Inverse estimate function
x0Fun <- function(object, y0 = params$y0) {
  beta <- as.numeric(fixef(object))
  (y0 - beta[1]) / beta[2]
}

## Function for coverting C.I. list into a two-column matrix
list2Matrix <- function(object) {
  matrix(unlist(object), nrow = length(object), ncol = length(object[[1]]), 
         byrow = TRUE)
}

## Function to summarize confidence interval
.summarize <- function(x) {
  .coverage <- if (x[1] <= params$x0 && params$x0 <= x[2]) 1 else 0
  .length <- x[2] - x[1]
  c(.coverage, .length)
}

.summarizeBoot <- function(x) {
  .fun <- function(z) {
    c(.summarize(z[1, ]), .summarize(z[2, ]), .summarize(z[3, ]))
  }
  (d <- ldply(x, .fun))
  PB <- apply(d[, 1:2], 2, mean)
  PB.inv1 <- apply(d[, 3:4], 2, mean)
  PB.inv2 <- apply(d[, 5:6], 2, mean)
  res <- rbind(PB, PB.inv1, PB.inv2)
  colnames(res) <- c("Coverage", "Length")
  res
}

## Function to calculate the Wald-based C.I.
waldCI <- function(.data) {
  ## FIXME: Should this be calculated based on the original model?
  Y0 <- rnorm(1, mean = params$y0, sd = sqrt(params$var.y0))
  ## Fit model using lme4 package and estimate x0 and Var(Y0)
  mod <- lmer(y ~ x + (0+1|subject) + (0+x|subject), data = .data)
  x0.est <- x0Fun(mod, y0 = Y0)
  var.y0 <- VarCorr(mod)[[1]][1] + x0.est^2*VarCorr(mod)[[2]][1] + 
    sigma(mod)^2
  ## Delta method
  beta <- as.numeric(fixef(mod))
  covmat <- diag(3)
  covmat[1:2, 1:2] <- as.matrix(vcov(mod))
  covmat[3, 3] <- var.y0
  params <- c(beta0 = beta[1], beta1 = beta[2], y0 = Y0)
  gstring <- "(y0 - beta0) / beta1"
  dm <- car:::deltaMethod(params, g = gstring, vcov. = covmat)
  rownames(dm) <- ""
  dm$Estimate + qnorm(c(0.025, 0.975))*dm$SE
}

## Function to calculate the inversion interval
invCI <- function(.data, q1 = qnorm(0.025), q2 = qnorm(0.975), Y0) {
  ## FIXME: Should this be calculated based on the original model?
  if (missing(Y0)) Y0 <- rnorm(1, mean = params$y0, sd = sqrt(params$var.y0))
  ## Fit model using nlme package and estimate x0 and Var(Y0)
  mod <- lme(y ~ x, random = list(subject = pdDiag(~x)), data = .data)
  x0.est <- x0Fun(mod, y0 = Y0)
  var.y0 <- getVarCov(mod)[1, 1] + getVarCov(mod)[2, 2]*x0.est^2 + 
    summary(mod)$sigma^2
  ## Prediction function that also returns standard error
  predFun <- function(x) {
    z <- list("x" = x)
    fit <- predict(mod, newdata = z, level = 0)
    se.fit <- sqrt(diag(cbind(1, unlist(z)) %*% mod$varFix %*% 
                          t(cbind(1, unlist(z)))))
    list(fit = fit, se.fit = se.fit)
  }
  ## Inverse functions for calculating confidence limits
  invFun1 <- function(x) { 
    z <- list(x)
    names(z) <- "volume"
    pred <- predFun(x)
    (Y0 - pred$fit)/sqrt((var.y0 + pred$se.fit^2)) - q2
  }
  invFun2 <- function(x) { 
    z <- list(x)
    names(z) <- "volume"
    pred <- predFun(x)
    (Y0 - pred$fit)/sqrt((var.y0 + pred$se.fit^2)) - q1
  }
  ## Find roots of inverse function
  lower <- uniroot(invFun1, interval = c(-1, x0.est), tol = 1e-10, 
                   maxiter = 1000)$root
  upper <- uniroot(invFun2, interval = c(x0.est, 2), tol = 1e-10, 
                   maxiter = 1000)$root
  c(lower, upper)
}

pbootCI <- function(.data, R = 999, .parallel = TRUE) {
  
  ## FIXME: Should this be calculated based on the original model?
  Y0 <- rnorm(1, mean = params$y0, sd = sqrt(params$var.y0))
  
  ## Fit model using lme4 package and estimate x0 and Var(Y0)
  mod <- lmer(y ~ x + (0+1|subject) + (0+x|subject), data = .data)
  x0.est <- x0Fun(mod, y0 = Y0)
  var.y0 <- VarCorr(mod)[[1]][1] + VarCorr(mod)[[2]][1]*x0.est^2 + 
    sigma(mod)^2
  
  ## Function to calculate bootstrap estimate
  bootFun <- function(.) {
    
    ## Extract model components
    covb <- as.matrix(vcov(.))        # (X' V^-1 X)^-1
    beta.boot <- as.numeric(fixef(.)) 
    
    ## Calculate estimates
    y0.boot <- rnorm(1, mean = Y0, sd = sqrt(var.y0))
    x0.boot <- x0Fun(., y0 = y0.boot)
    mu0.boot <-  as.numeric(crossprod(beta.boot, c(1, x0.est)))
    
    ## FIXME: Should these variances be calculated at x0.est or x0.boot?
    var1.y0 <- VarCorr(.)[[1]][1] + VarCorr(.)[[2]][1]*x0.boot^2 + sigma(.)^2
    var1.mu0 <- t(c(1, x0.boot)) %*% covb %*% c(1, x0.boot)
    var2.y0 <- VarCorr(.)[[1]][1] + VarCorr(.)[[2]][1]*x0.est^2 + sigma(.)^2
    var2.mu0 <- t(c(1, x0.est)) %*% covb %*% c(1, x0.est)
    
    Q1.boot <- (y0.boot - mu0.boot)/sqrt(var1.y0+var1.mu0)
    Q2.boot <- (y0.boot - mu0.boot)/sqrt(var2.y0+var2.mu0)
    
    c(x0.boot, Q1.boot, Q2.boot)
    
  } 
  
  ## Function that returns original estimate (i.e., no random y0)
  bootFun0 <- function(.) {
    
    ## Extract model components
    covb <- as.matrix(vcov(.))        # (X' V^-1 X)^-1
    beta.boot <- as.numeric(fixef(.)) 
    
    ## Calculate estimates
    y0.boot <- Y0
    x0.boot <- x0Fun(., y0 = y0.boot)
    mu0.boot <-  as.numeric(crossprod(beta.boot, c(1, x0.est)))

    ## FIXME: Should these variances be calculated at x0.est or x0.boot?
    var1.y0 <- VarCorr(.)[[1]][1] + VarCorr(.)[[2]][1]*x0.boot^2 + sigma(.)^2
    var1.mu0 <- t(c(1, x0.boot)) %*% covb %*% c(1, x0.boot)
    var2.y0 <- VarCorr(.)[[1]][1] + VarCorr(.)[[2]][1]*x0.est^2 + sigma(.)^2
    var2.mu0 <- t(c(1, x0.est)) %*% covb %*% c(1, x0.est)
    
    Q1.boot <- (y0.boot - mu0.boot)/sqrt(var1.y0+var1.mu0)
    Q2.boot <- (y0.boot - mu0.boot)/sqrt(var2.y0+var2.mu0)
    
    c(x0.boot, Q1.boot, Q2.boot)
    
  }
  
  ## Calculate quantiles of bootstrap sample
  x0.pb <- if (.parallel) {
    bootMer2_parallel(mod, FUN = bootFun, FUN0 = bootFun0, nsim = R,
                      parallel = "multicore", ncpus = 4)
  } else {
    bootMer2(mod, FUN = bootFun, FUN0 = bootFun0, nsim = R)
  }
#   as.numeric(quantile(x0.pb$t, c(0.025, 0.975)))
#   x0.pb
  quants1 <- as.numeric(quantile(x0.pb$t[, 2], c(0.025, 0.975)))
  quants2 <- as.numeric(quantile(x0.pb$t[, 3], c(0.025, 0.975)))
  rbind(as.numeric(quantile(x0.pb$t[, 1], c(0.025, 0.975))),
        invCI(.data, q1 = quants1[1], q2 = quants1[2], Y0 = Y0),
        invCI(.data, q1 = quants2[1], q2 = quants2[2], Y0 = Y0))
  
}
# 
# getBootCI <- function(x, .which = 1) {
#   quantile(x$t[, .which], c(0.025, 0.975))
# }
# 
# res <- llply(dfs, pbootCI, .progress = "text")
# res2 <- llply(res, getBootCI, .which = 2)
# apply(apply(list2Matrix(res2), 1, .summarize), 1, mean)

## Simulation ------------------------------------------------------------------

## Simulate data frames
set.seed(5746)
dfs <- rlply(params$nsim, simData)

## Fit models to sample data
set.seed(8306)
simdata <- simData()
xyplot(y ~ x, groups = subject, data = simdata, type = "b",
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...) 
         panel.abline(h = params$y0, v = params$x0)
})
plot(intervals(lmList(y ~ I(x-mean(x)) | subject, data = simdata)))
mod.lme4 <- lmer(y ~ x + (0+1|subject) + (0+x|subject), data = simdata)
mod.nlme <- lme(y ~ x, random = list(subject = pdDiag(~x)), data = simdata)
x0Fun(mod.lme4)
x0Fun(mod.nlme)

## Simulation for the Wald-based interval
wald.cis <- list2Matrix(llply(dfs, waldCI, .progress = "text"))
apply(apply(wald.cis, 1, .summarize), 1, mean)

## Simulation for the inversion interval
inv.cis <- list2Matrix(llply(dfs, invCI, .progress = "text"))
apply(apply(inv.cis, 1, .summarize), 1, mean)

## Simulation for the PB percentile interval -----------------------------------
pboot.cis <- llply(dfs, pbootCI, .progress = "text")
.summarizeBoot(pboot.cis)
## Two particular occasions produced:
## [1] 0.9500000 0.3340199
## [1] 0.9200000 0.3395194

## Test parallel version
# system.time(res1 <- pbootCI(simdata))
# system.time(res2 <- pbootCI(simdata, .parallel = TRUE))
# rbind(res1, res2)

## Save results
cal_lmm_sim1 <- list(  
  params = params,
  wald.cis = wald.cis,
  inv.cis = inv.cis,
  pboot.cis = pboot.cis
)
save(cal_lmm_sim1, file = "/home/w108bmg/Desktop/Dissertation/Data/cal_lmm_sim1.RData")
