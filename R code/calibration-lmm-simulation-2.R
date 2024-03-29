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
library(rootSolve)
library(Rcpp)
sourceCpp("/home/w108bmg/Desktop/simfuns.cpp")
source("/home/w108bmg/Desktop/Dissertation/R code/bootMer2.R")
source("/home/w108bmg/Desktop/Dissertation/R code/bootMer2_parallel.R")

## Simulation setup ------------------------------------------------------------

## Parameters
params <- list(
  nsim = 1000,                      # simulation size
  m = 30,                           # number of subjects
  n = 20,                           # sample size per subject
  beta = c(0, 3, -1),               # fixed effecs
  theta = c(0.0001, 0.05, 0.001),   # variance components
  y0 = c(0, 0.5, 1, 1.5, 2)[5] # true observed response
)
params$x0 <- xinv_cpp(params$beta[3], params$beta[2], params$beta[1]-params$y0, 
                      lower = 0, upper = 1)
params$var.y0 <- params$theta[1] + params$theta[2]*params$x0^2 + params$theta[3]

# (-params$beta[2] + sqrt(params$beta[2]^2 - 4*params$beta[3]*(params$beta[1]-params$y0))) / (2*params$beta[3])

f1 <- function(x) 3*x - x^2         # quadratic
f2 <- function(x) x^3 - 3*x^2 + 3*x # cubic
par(mfrow = c(2, 1))
curve(f1, from = 0, to = 1)
curve(f2, from = 0, to = 1)

## Function to simulate data from a random intercept and slope model
simData <- function(n = params$n, m = params$m, beta = params$beta, 
                    theta = params$theta) {
  subject <- rep(1:m, each = n)
  x <- rep(seq(from = 0, to = 1, length = n), times = m) 
  B0 <- rnorm(m, mean = 0, sd = sqrt(theta[1]))
  B1 <- rnorm(m, mean = 0, sd = sqrt(theta[2]))
  y <- rnorm(m*n, mean = beta[1]+B0[subject] + (beta[2]+B1[subject])*x + 
               beta[3]*x^2, 
             sd = sqrt(theta[3]))
  data.frame(x = x, y = y, subject = factor(subject))
}

## Inverse estimate function (just use xinv_cpp)
xinv <- function(object, y0 = params$y0, lower = 0, upper = 1.5) {
  fixed <- as.numeric(fixef(object))
  xinv_cpp(fixed[3], fixed[2], fixed[1]-y0, lower, upper)
}

# solveable <- function(object, y0, lower = 0, upper = 1) {
#   coefs <- as.numeric(fixef(object))
#   roots <- try(xinv_cpp(coefs[3], coefs[2], coefs[1]-y0, lower, upper), 
#                silent = TRUE)
#   if (inherits(roots, "try-error")) FALSE else TRUE
# }
# 
# ## Function to summarize confidence interval
# simulationSummary <- function(x, boot = FALSE) {
#   
#   ## Function for coverting list of C.I.'s into a two-column matrix
#   list2Matrix <- function(object) {
#     matrix(unlist(object), nrow = length(object), ncol = length(object[[1]]), 
#            byrow = TRUE)
#   }
#   
#   ## Function that returns coverage and length of a single C.I.
#   coverage.and.length <- function(x) {
#     .coverage <- if (x[1] <= params$x0 && params$x0 <= x[2]) 1 else 0
#     .length <- x[2] - x[1]
#     c(.coverage, .length)
#   }
#   
#   if (!boot) {
#     
#     ci.mat <- list2Matrix(x)
#     res <- apply(apply(ci.mat, 1, coverage.and.length), 1, mean)
#     names(res) <- c("Coverage", "Length")
#     res
#     
#   } else {
#     
#     ## Data frame of coverge and length estimates
#     d <- ldply(x, function(x) {
#       c(coverage.and.length(x[1, ]), # normal
#         coverage.and.length(x[2, ]), # basic
#         coverage.and.length(x[3, ]), # percentile
#         coverage.and.length(x[4, ])) # adjusted inversion
#     })
#     
#     ## Calculate mean coverage and mean length
#     boot.norm <- apply(d[, 1:2], 2, mean)
#     boot.basic <- apply(d[, 3:4], 2, mean)
#     boot.perc <- apply(d[, 5:6], 2, mean)
#     boot.inv <- apply(d[, 7:8], 2, mean)
#     res <- rbind(boot.norm, boot.basic, boot.perc, boot.inv)
#     colnames(res) <- c("Coverage", "Length")
#     res
#     
#   }
#   
# }
# 
# ## Function to calculate the Wald-based C.I.
# waldCI <- function(.data) {
#   
#   mod <- lme(y ~ x + I(x^2), random = list(subject = pdDiag(~x)), data = .data)
#   repeat {
#     Y0 <- rnorm(1, mean = params$y0, sd = sqrt(params$var.y0))
#     if (solveable(mod, y0 = Y0)) break
#   }
#   x0.est <- x0Fun(mod, y0 = Y0)
#   var.y0 <- getVarCov(mod)[1, 1] + getVarCov(mod)[2, 2]*x0.est^2 + 
#     summary(mod)$sigma^2
#   
#   ## Delta method
#   beta <- as.numeric(fixef(mod))
#   covmat <- diag(4)
#   covmat[1:3, 1:3] <- as.matrix(vcov(mod))
#   covmat[4, 4] <- var.y0
#   params <- c(beta0 = beta[1], beta1 = beta[2], beta2 = beta[3], y0 = Y0)
#   gstring <- "(-beta1 + sqrt(beta1^2 - 4*beta2*(beta0-y0))) / (2*beta2)"
#   dm <- car:::deltaMethod(params, g = gstring, vcov. = covmat)
#   rownames(dm) <- ""
#   dm$Estimate + qnorm(c(0.025, 0.975))*dm$SE
#   
# }
# 
# ## Function to calculate the inversion interval
# invCI <- function(.data, q1 = qnorm(0.025), q2 = qnorm(0.975), Y0) {
#   
#   mod <- lme(y ~ x + I(x^2), random = list(subject = pdDiag(~x)), data = .data)
#   if (missing(Y0)) {
#     repeat {
#       Y0 <- rnorm(1, mean = params$y0, sd = sqrt(params$var.y0))
#       if (solveable(mod, y0 = Y0)) break
#     }
#   }
#   x0.est <- x0Fun(mod, y0 = Y0)
#   var.y0 <- getVarCov(mod)[1, 1] + getVarCov(mod)[2, 2]*x0.est^2 + 
#     summary(mod)$sigma^2
#   
#   ## Prediction function that also returns standard error
#   predFun <- function(x) {
#     z <- list("x" = x)
#     fit <- predict(mod, newdata = z, level = 0)
#     se.fit <- sqrt(diag(cbind(1, unlist(z), unlist(z)^2) %*% mod$varFix %*% 
#                           t(cbind(1, unlist(z), unlist(z)^2))))
#     list(fit = fit, se.fit = se.fit)
#   }
#   
# #   ## Inverse function for calculating confidence limits
# #   invFun.bounds <- function(x) { 
# #     z <- list("x" = x)
# #     pred <- predFun(x)
# #     (Y0 - pred$fit)^2/(var.y0 + pred$se.fit^2) - qnorm(0.975)^2
# #   }
# #   ## Find roots of inverse function
# #   c(uniroot(invFun.bounds, interval = c(-1, x0.est), tol = 1e-10, 
# #             maxiter = 1000)$root, 
# #     uniroot(invFun.bounds, interval = c(x0.est, 2), tol = 1e-10, 
# #             maxiter = 1000)$root)
#   
#   ## Find roots of predictive pivot
#   invFun1 <- function(x) { 
#     z <- list(x)
#     names(z) <- "volume"
#     pred <- predFun(x)
#     (Y0 - pred$fit)/sqrt((var.y0 + pred$se.fit^2)) - q2
#   }
#   invFun2 <- function(x) { 
#     z <- list(x)
#     names(z) <- "volume"
#     pred <- predFun(x)
#     (Y0 - pred$fit)/sqrt((var.y0 + pred$se.fit^2)) - q1
#   }
#   roots <- c(uniroot.all(invFun1, interval = c(-1, 4), tol = 1e-10, 
#                          maxiter = 1000),
#              uniroot.all(invFun2, interval = c(-1, 4), tol = 1e-10, 
#                          maxiter = 1000))
#   sort(roots)[1:2]
# 
# }
# 
# pbootCI <- function(.data, R = 999, .parallel = TRUE) {
# 
#   ## Fit model and calculate estimates of x0 and Var(Y0)
#   mod <- lmer(y ~ x + I(x^2) + (0+1|subject) + (0+x|subject), data = .data)
#   repeat {
#     Y0 <- rnorm(1, mean = params$y0, sd = sqrt(params$var.y0))
#     if (solveable(mod, y0 = Y0)) break
#   }
#   x0.est <- x0Fun(mod, y0 = Y0)
#   var.y0 <- VarCorr(mod)[[1]][1] + VarCorr(mod)[[2]][1]*x0.est^2 + 
#     sigma(mod)^2
# 
#   ## Function to calculate bootstrap estimate
#   bootFun <- function(.) {
#     
#     ## Extract model components
#     covb <- as.matrix(vcov(.))        # (X' V^-1 X)^-1
#     beta.boot <- as.numeric(fixef(.)) 
#     
#     ## Calculate estimates
#     y0.boot <- rnorm(1, mean = Y0, sd = sqrt(var.y0))
#     x0.boot <- x0Fun(., y0 = y0.boot)
#     mu0.boot <-  as.numeric(crossprod(beta.boot, c(1, x0.est, x0.est^2)))
#     
#     ## FIXME: Should these variances be calculated at x0.est or x0.boot?
#     var1.y0 <- VarCorr(.)[[1]][1] + VarCorr(.)[[2]][1]*x0.boot^2 + sigma(.)^2
#     var1.mu0 <- t(c(1, x0.boot, x0.boot^2)) %*% covb %*% c(1, x0.boot, x0.boot^2)
#     var2.y0 <- VarCorr(.)[[1]][1] + VarCorr(.)[[2]][1]*x0.est^2 + sigma(.)^2
#     var2.mu0 <- t(c(1, x0.est, x0.est^2)) %*% covb %*% c(1, x0.est, x0.est^2)
#     
#     Q1.boot <- (y0.boot - mu0.boot)/sqrt(var1.y0+var1.mu0)
#     Q2.boot <- (y0.boot - mu0.boot)/sqrt(var2.y0+var2.mu0)
#     
#     c(x0.boot, Q1.boot, Q2.boot)
#     
#   } 
#   
#   ## Function that returns original estimate (i.e., no random y0)
#   bootFun0 <- function(.) {
#     
#     ## Extract model components
#     covb <- as.matrix(vcov(.))        # (X' V^-1 X)^-1
#     beta.boot <- as.numeric(fixef(.)) 
#     
#     ## Calculate estimates
#     y0.boot <- Y0
#     x0.boot <- x0Fun(., y0 = y0.boot)
#     mu0.boot <-  as.numeric(crossprod(beta.boot, c(1, x0.est, x0.est^2)))
#     
#     ## FIXME: Should these variances be calculated at x0.est or x0.boot?
#     var1.y0 <- VarCorr(.)[[1]][1] + VarCorr(.)[[2]][1]*x0.boot^2 + sigma(.)^2
#     var1.mu0 <- t(c(1, x0.boot, x0.boot^2)) %*% covb %*% c(1, x0.boot, x0.boot^2)
#     var2.y0 <- VarCorr(.)[[1]][1] + VarCorr(.)[[2]][1]*x0.est^2 + sigma(.)^2
#     var2.mu0 <- t(c(1, x0.est, x0.est^2)) %*% covb %*% c(1, x0.est, x0.est^2)
#     
#     Q1.boot <- (y0.boot - mu0.boot)/sqrt(var1.y0+var1.mu0)
#     Q2.boot <- (y0.boot - mu0.boot)/sqrt(var2.y0+var2.mu0)
#     
#     c(x0.boot, Q1.boot, Q2.boot)
#     
#   }
#   
#   ## Return bootstrap samples
#   x0.pb <- if (.parallel) {
#              bootMer2_parallel(mod, FUN = bootFun, FUN0 = bootFun0, nsim = R,
#                                parallel = "multicore", ncpus = 4)
#            } else {
#              bootMer2(mod, FUN = bootFun, FUN0 = bootFun0, nsim = R)
#            }
# #   quants1 <- as.numeric(quantile(x0.pb$t[, 2], c(0.025, 0.975), na.rm = TRUE))
# #   quants2 <- as.numeric(quantile(x0.pb$t[, 3], c(0.025, 0.975), na.rm = TRUE))
# #   rbind(as.numeric(quantile(x0.pb$t[, 1], c(0.025, 0.975), na.rm = TRUE)),
# #         invCI(.data, q1 = quants1[1], q2 = quants1[2], Y0 = Y0),
# #         invCI(.data, q1 = quants2[1], q2 = quants2[2], Y0 = Y0))
#   x0.pb
#   
# }

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
# plot(intervals(lmList(y ~ poly(x, degree = 2) | subject, data = simdata)))
mod.lme4 <- lmer(y ~ x + I(x^2) + (0+1|subject) + (0+x|subject), data = simdata)
mod.nlme <- lme(y ~ x + I(x^2), random = list(subject = pdDiag(~x)), data = simdata)
xinv(mod.lme4)
xinv(mod.nlme)

## Simulation for the Wald-based interval
wald.cis <- list2Matrix(llply(dfs, waldCI, .progress = "text"))
apply(apply(wald.cis, 1, .summarize), 1, mean)

## Simulation for the inversion interval
inv.cis <- list2Matrix(llply(dfs, invCI, .progress = "text"))
apply(apply(inv.cis, 1, .summarize), 1, mean)

## Simulation for the PB percentile interval -----------------------------------
pboot.cis <- llply(dfs, pbootCI, .progress = "text")
save(pboot.cis, file = "/home/w108bmg/Desktop/Dissertation/Data/pboot_quad_20.RData")

summarizeBoot <- function(z) {
  fun.1 <- function(x) x[1, ]
  fun.2 <- function(x) x[2, ]
  fun.3 <- function(x) x[3, ]
  cis.1 <- apply(apply(ldply(z, fun.1), 1, .summarize), 1, mean)
  cis.2 <- apply(apply(ldply(z, fun.2), 1, .summarize), 1, mean)
  cis.3 <- apply(apply(ldply(z, fun.3), 1, .summarize), 1, mean)
  names(cis.1) <- names(cis.2) <- names(cis.3) <- c("Coverage", "Length")
  rbind(PB = cis.1, PB.inv1 = cis.2, PB.inv2 = cis.3)
}


anyMissing <- function(.) {
  any(is.na(.))
}

is.missing <- sapply(pboot.cis, anyMissing)
not.missing <- !is.missing
pboot.cis.not.missing <- (pboot.cis[!.missing])
summarizeBoot(pboot.cis.not.missing)

num.not.missing <- length(pboot.cis.not.missing)
i <- num.not.missing + 1
while (i <= 1000) {
  cis <- pbootCI(simData())
  if (!any(is.na(cis))) {
    pboot.cis.not.missing <- append(pboot.cis.not.missing, list(cis))
    i <- i + 1
  }
}
summarizeBoot(pboot.cis.not.missing)
