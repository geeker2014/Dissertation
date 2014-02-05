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
library(boot)

## Simulating data -------------------------------------------------------------
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

simdata <- simRCD(m = 15, n = 30, fixed = c(0, 2), vars = c(0.001, 0.05, 0.001))
xyplot(y ~ x, groups = subject, data = simdata)
fit <- lmer(y ~ x + (0+1|subject) + (0+x|subject), data = simdata)

## True population parameters
beta0 <- 0                       # fixed intercept
beta1 <- 2                       # fixed slope
var0 <- 0.001                    # variance of random intercepts
var1 <- 0.05                     # variance of random slopes
var.error <- 0.001               # error variance
y0.true <- 1.5                   # observed response
x0.true <- (y0 - beta0) / beta1  # true unknown
var.y0.true <- var0 + var1*x0.true^2 + var.error

## Simulate nsim data frames
nsim <- 100
set.seed(5746)
dfs <- rlply(nsim, simRCD(m = 15, n = 30, fixed = c(0, 2), 
                          vars = c(0.001, 0.05, 0.001)))

## Inverse estimate function
x0Fun <- function(object, y0 = y0.true) {
  beta <- as.numeric(fixef(fit))
  (y0 - beta[1]) / beta[2]
}
x0.est <- x0Fun(fit)

## Function for coverting C.I. list into a two-column matrix
list2Matrix <- function(object) {
  m <- length(object)
  n <- length(object[[1]])
  matrix(unlist(object), nrow = m, ncol = n, byrow = TRUE)
}

## Function to check coverage
checkCoverage <- function(x, x0 = x0.true) {
  if (x[1] <= x0 && x0 <= x[2]) 1 else 0
}
  
## Simulate coverage probability and length for the Wald-based method ----------

## Function to calculate the Wald-based C.I.
waldCI <- function(.data, y0 = y0.true) {
  fit <- lmer(y ~ x + (0+1|subject) + (0+x|subject), data = .data)
  beta <- as.numeric(fixef(fit))
  covmat <- diag(3)
  covmat[1:2, 1:2] <- as.matrix(vcov(fit))
  var.y0 <- VarCorr(fit)[[1]][1] + x0.est^2*VarCorr(fit)[[2]][1] + 
                  sigma(fit)^2
  covmat[3, 3] <- var.y0
  
  ## FIXME: Should this be calculated based on the original model?
  Y0 <- rnorm(1, mean = y0, sd = sqrt(var.y0.true))
  
  params <- c(beta0 = beta[1], beta1 = beta[2], y0 = Y0)
  gstring <- "(y0 - beta0) / beta1"
  dm <- car:::deltaMethod(params, g = gstring, vcov. = covmat)
  rownames(dm) <- ""
  dm$Estimate + qnorm(c(0.025, 0.975))*dm$SE
}

## Apply to each data frame
wald.cis <- list2Matrix(llply(dfs, waldCI, .progress = "text"))
mean(apply(wald.cis, 1, checkCoverage))

## Simulate coverage probability and length for the inversion method -----------

## Function to calculate the inversion interval
invCI <- function(.data) {
  
  ## Fit model using nlme package and estimate Var(Y0)
  mod <- nlme(y ~ x, random = list(subject = pdDiag(~x)), data = .data)
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
    (y0 - pred$fit)^2/(var.y0 + pred$se.fit^2) - qnorm(0.975)^2
  }
  
  ## Find roots of inverse function
  c(uniroot(invFun.bounds, interval = c(0, x0.est), tol = 1e-10, 
            maxiter = 1000)$root, 
    uniroot(invFun.bounds, interval = c(x0.est, 1), tol = 1e-10, 
            maxiter = 1000)$root)
  
}