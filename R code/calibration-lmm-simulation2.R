## Simulation study for linear calibration with grouped data.
##
## Author: Brandon M. Greenwell
## Date: 2/4/2014
##
## TODO:
##   (1) Should Y0 be simulated ahead of time; for example,
##       Y0 <- rnorm(nsim, 0, var.y0)? Is this even possible?

## Load packages ---------------------------------------------------------------
library(lattice)
library(plyr)
library(lme4)
library(nlme)
library(boot)
source("/home/w108bmg/Desktop/Dissertation/R code/bootMer2.R")
source("/home/w108bmg/Desktop/Dissertation/R code/bootMer2_parallel.R")

## Simulating data -------------------------------------------------------------

## Function to simulate data from a random intercept and slope model with a 
## quadratic term
simRCD <- function(n = 10, m = 20, fixed = c(0, 1, -0.5), 
                   vars = c(0.0001, 0.0075, 0.001)) {
  subject <- rep(1:m, each = n)
  x <- rep(seq(from = 0, to = 1, length = n), times = m) 
  B0 <- rnorm(m, mean = 0, sd = sqrt(vars[1]))
  B1 <- rnorm(m, mean = 0, sd = sqrt(vars[2]))
  y <- rnorm(m*n, mean = fixed[1]+B0[subject] + (fixed[2]+B1[subject])*x + 
               fixed[3]*x^2, 
             sd = sqrt(vars[3]))
  data.frame(x = x, y = y, subject = factor(subject))
}

## True population parameters
beta0 <- 0                       # fixed intercept
beta1 <- 1                       # fixed slope
beta2 <- -0.5
var0 <- 0.0001                   # variance of random intercepts
var1 <- 0.0075                   # variance of random slopes
var.error <- 0.001               # error variance
y0.true <- 0.4                   # observed response
x0.true <- (-beta1 + sqrt(beta1^2 - 4*beta2*(beta0-y0.true))) / (2*beta2)
var.y0.true <- var0 + var1*x0.true^2 + var.error

## Simulate data frames
nsim <- 100 # change to 10 or 100 for testing
set.seed(5746)
dfs <- rlply(nsim, simRCD(n = 10, m = 20, fixed = c(0, 1, -0.5), 
                          vars = c(0.0001, 0.0075, 0.001)))

## Fit models to sample data
set.seed(8306)
simdata <- simRCD(n = 10, m = 20, fixed = c(0, 1, -0.5), 
                  vars = c(0.0001, 0.0075, 0.001))
xyplot(y ~ x, groups = subject, data = simdata, type = "b",
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...) 
         panel.abline(h = y0.true, v = x0.true)
})
plot(intervals(lmList(y ~ poly(x, degree = 2) | subject, data = simdata)))
mod.lme4 <- lmer(y ~ x + I(x^2) + (0+1|subject) + (0+x|subject), data = simdata)
mod.nlme <- lme(y ~ x + I(x^2), random = list(subject = pdDiag(~x)), data = simdata)

## Inverse estimate function
x0Fun <- function(object, y0 = y0.true) {
  beta <- as.numeric(fixef(object))
  (-beta[2] + sqrt(beta[2]^2 - 4*beta[3]*(beta[1]-y0))) / (2*beta[3])
}
x0.est <- x0Fun(mod.lme4)
x0Fun(mod.lme4)
x0Fun(mod.nlme)

## Function for coverting C.I. list into a two-column matrix
list2Matrix <- function(object) {
  m <- length(object)
  n <- length(object[[1]])
  matrix(unlist(object), nrow = m, ncol = n, byrow = TRUE)
}

## Function to summarize confidence interval
.summarize <- function(x, x0 = x0.true) {
  .coverage <- if (x[1] <= x0 && x0 <= x[2]) 1 else 0
  .length <- x[2] - x[1]
  c(.coverage, .length)
}
  
## Simulate coverage probability and length for the Wald-based method ----------

## Function to calculate the Wald-based C.I.
waldCI <- function(.data, y0 = y0.true) {
  
  ## FIXME: Should this be calculated based on the original model?
  Y0 <- rnorm(1, mean = y0, sd = sqrt(var.y0.true))
  
  ## Fit model using lme4 package and estimate x0 and Var(Y0)
  mod <- lmer(y ~ x + I(x^2) + (0+1|subject) + (0+x|subject), data = .data)
  x0.est <- x0Fun(mod, y0 = Y0)
  var.y0 <- VarCorr(mod)[[1]][1] + x0.est^2*VarCorr(mod)[[2]][1] + 
    sigma(mod)^2
  
  ## Delta method
  beta <- as.numeric(fixef(mod))
  covmat <- diag(4)
  covmat[1:3, 1:3] <- as.matrix(vcov(mod))
  covmat[4, 4] <- var.y0
  params <- c(beta0 = beta[1], beta1 = beta[2], beta2 = beta[3], y0 = Y0)
  gstring <- "(-beta1 + sqrt(beta1^2 - 4*beta2*(beta0-Y0))) / (2*beta2)"
  dm <- car:::deltaMethod(params, g = gstring, vcov. = covmat)
  rownames(dm) <- ""
  dm$Estimate + qnorm(c(0.025, 0.975))*dm$SE
}

## Apply to each data frame
wald.cis <- list2Matrix(llply(dfs, waldCI, .progress = "text"))
apply(apply(wald.cis, 1, .summarize), 1, mean)

## Simulate coverage probability and length for the inversion method -----------

## Function to calculate the inversion interval
invCI <- function(.data, y0 = y0.true) {

  ## FIXME: Should this be calculated based on the original model?
  Y0 <- rnorm(1, mean = y0, sd = sqrt(var.y0.true))
  
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
  
  ## Inverse function for calculating confidence limits
  invFun.bounds <- function(x) { 
    z <- list("x" = x)
    pred <- predFun(x)
    (Y0 - pred$fit)^2/(var.y0 + pred$se.fit^2) - qnorm(0.975)^2
  }
  
  ## Find roots of inverse function
  c(uniroot(invFun.bounds, interval = c(0, x0.est), tol = 1e-10, 
            maxiter = 1000)$root, 
    uniroot(invFun.bounds, interval = c(x0.est, 2), tol = 1e-10, 
            maxiter = 1000)$root)
  
}

## Apply to each data frame
inv.cis <- list2Matrix(llply(dfs, invCI, .progress = "text"))
apply(apply(inv.cis, 1, .summarize), 1, mean)

## Simulate coverage probability and length for the PB approach ----------------

pbootCI <- function(.data, y0 = y0.true, R = 999, .parallel = TRUE) {
  
  ## FIXME: Should this be calculated based on the original model?
  Y0 <- rnorm(1, mean = y0, sd = sqrt(var.y0.true))
  
  ## Fit model using lme4 package and estimate x0 and Var(Y0)
  mod <- lmer(y ~ x + (0+1|subject) + (0+x|subject), data = .data)
  x0.est <- x0Fun(mod, y0 = Y0)
  var.y0 <- VarCorr(mod)[[1]][1] + VarCorr(mod)[[2]][1]*x0.est^2 + 
    sigma(mod)^2
  
  ## Function to calculate bootstrap estimate
  bootFun <- function(.) {
    y0.boot <- rnorm(1, mean = Y0, sd = sqrt(var.y0))
    x0Fun(., y0 = y0.boot)
  }
  
  ## Function to return original estimate
  bootFun0 <- function(.) {
    x0Fun(., y0 = Y0)
  }
  
  ## Calculate quantiles of bootstrap sample
  x0.pb <- if (.parallel) {
    bootMer2_parallel(mod, FUN = bootFun, FUN0 = bootFun0, nsim = R,
             parallel = "multicore", ncpus = 4)
  } else {
    bootMer2(mod, FUN = bootFun, FUN0 = bootFun0, nsim = R)
  }
  as.numeric(quantile(x0.pb$t, c(0.025, 0.975)))
  
}

## Test parallel version
system.time(res1 <- pbootCI(simdata))
system.time(res2 <- pbootCI(simdata, .parallel = TRUE))
rbind(res1, res2)

## Apply to each data frame
pboot.cis <- list2Matrix(llply(dfs, pbootCI, .progress = "text"))
apply(apply(pboot.cis, 1, .summarize), 1, mean)
## One particular occasion produced:
## [1] 0.9500000 0.3340199