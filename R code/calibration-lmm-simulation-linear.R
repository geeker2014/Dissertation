################################################################################
##
## Simulation study for linear calibration with grouped data.
##
## Author: Brandon M. Greenwell
## Date: 2/4/2014
##
################################################################################

## Load packages ---------------------------------------------------------------
library(lattice) # for plotting grouped data
library(plyr)    # for working with the data
library(lme4)    # for fitting mixed models and bootstrap inference
library(nlme)    # for fitting mixed models
library(boot)    # for calculating bootstrap confidence intervals
source("/home/w108bmg/Desktop/Dissertation/R code/Bootstrap functions.R")
source("/home/w108bmg/Desktop/Dissertation/R code/uniroot2.R")

## Functions -------------------------------------------------------------------

## Function to simulate data from a random intercept and slope model
genData <- function(m = 15, n = 10, beta = c(0, 1), 
                    theta = c(0.01, 0.05, 0.001)) {
  subject <- rep(1:m, each = n)
  x <- rep(seq(from = 0, to = 1, length = n), times = m) 
  B0 <- rnorm(m, mean = 0, sd = sqrt(theta[1]))
  B1 <- rnorm(m, mean = 0, sd = sqrt(theta[2]))
  y <- rnorm(m*n, mean = beta[1]+B0[subject] + (beta[2]+B1[subject])*x, 
             sd = sqrt(theta[3]))
  data.frame(x = x, y = y, subject = factor(subject))
}

## Function to calculate inverse estimate given a mixed model object
xest <- function(object, y0) {
  fixed <- as.numeric(fixef(object))
  (y0 - fixed[1]) / fixed[2]
}

## Function for coverting list of C.I.'s into a two-column matrix
list2Matrix <- function(object) {
  matrix(unlist(object), nrow = length(object), ncol = length(object[[1]]), 
         byrow = TRUE)
}

## Function that returns coverage and length of a single C.I.
covlen <- function(x, x0) {
  .coverage <- if (x[1] <= x0 && x0 <= x[2]) 1 else 0
  .length <- x[2] - x[1]
  c(.coverage, .length)
}

## Function to summarize confidence intervals
simulationSummary <- function(x, params, boot = FALSE) {
  if (!boot) { # summarize nonbootstrap intervals
    ci.mat <- list2Matrix(x)
    res <- apply(apply(ci.mat, 1, covlen, x0 = params$x0), 1, mean)
    names(res) <- c("Coverage", "Length")
    res
  } else { # summarize bootstrap intervals
    d <- ldply(x, function(x) {
      c(covlen(x[1, ], x0 = params$x0), # normal
        covlen(x[2, ], x0 = params$x0), # basic
        covlen(x[3, ], x0 = params$x0), # percentile
        covlen(x[4, ], x0 = params$x0)) # adjusted inversion
    })
    boot.norm <- apply(d[, 1:2], 2, mean)
    boot.basic <- apply(d[, 3:4], 2, mean)
    boot.perc <- apply(d[, 5:6], 2, mean)
    boot.inv <- apply(d[, 7:8], 2, mean)
    res <- rbind(boot.norm, boot.basic, boot.perc, boot.inv)
    colnames(res) <- c("Coverage", "Length")
    res
  }
}

## Function to calculate the Wald-based C.I.
waldCI <- function(.data, params, Y0) {
  ## FIXME: Should this be calculated based on the original model?
  if (missing(Y0)) Y0 <- rnorm(1, mean = params$y0, sd = sqrt(params$var.y0))
  ## Fit model using lme4 package and estimate x0 and Var(Y0)
  mod <- lmer(y ~ x + (0+1|subject) + (0+x|subject), data = .data)
  x0.est <- xest(mod, y0 = Y0)
  var.y0 <- VarCorr(mod)[[1]][1] + x0.est^2*VarCorr(mod)[[2]][1] + 
    sigma(mod)^2
  ## Delta method
  beta <- as.numeric(fixef(mod))
  covmat <- diag(3)
  covmat[1:2, 1:2] <- as.matrix(vcov(mod))
  covmat[3, 3] <- var.y0
  dm <- car:::deltaMethod(c(beta0 = beta[1], beta1 = beta[2], y0 = Y0), 
                          g = "(y0 - beta0) / beta1", vcov. = covmat)
  rownames(dm) <- ""
  dm$Estimate + qnorm(c(0.025, 0.975))*dm$SE
}

## Function to calculate the inversion interval
invCI <- function(.data, params, q1 = qnorm(0.025), q2 = qnorm(0.975), Y0) {
  ## FIXME: Should this be calculated based on the original model?
  if (missing(Y0)) Y0 <- rnorm(1, mean = params$y0, sd = sqrt(params$var.y0))
  ## Fit model using nlme package and estimate x0 and Var(Y0)
  mod <- lme(y ~ x, random = list(subject = pdDiag(~x)), data = .data)
  x0.est <- xest(mod, y0 = Y0)
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
  fun1 <- function(x) { 
    pred <- predFun(x)
    (Y0 - pred$fit)/sqrt((var.y0 + pred$se.fit^2)) - q2
  }
  fun2 <- function(x) { 
    pred <- predFun(x)
    (Y0 - pred$fit)/sqrt((var.y0 + pred$se.fit^2)) - q1
  }
  ## Find roots of inverse function
  lower <- uniroot2(fun1, interval = c(-3, x0.est), tol = 1e-10, 
                    maxiter = 1000, extend = "left")$root
  upper <- uniroot2(fun2, interval = c(x0.est, 4), tol = 1e-10, 
                    maxiter = 1000, extend = "right")$root
  c(lower, upper)
}

## Function to calculate the inversion interval without LMM
inv2CI <- function(.data, params, Y0) {
  if (missing(Y0)) Y0 <- rnorm(1, mean = params$y0, sd = sqrt(params$var.y0))
  cal <- calibrate(.data[, 1:2], y0 = Y0, interval = "inversion")
  c(cal$lower, cal$upper)
}


## Function to calculate bootstrap intervals
pbootCI <- function(.data, params, R = 999, .parallel = TRUE, Y0) {
  ## FIXME: Should this be calculated based on the original model?
  if (missing(Y0)) Y0 <- rnorm(1, mean = params$y0, sd = sqrt(params$var.y0))
  ## Fit model using lme4 package and estimate x0 and Var(Y0)
  mod <- lmer(y ~ x + (0+1|subject) + (0+x|subject), data = .data)
  x0.est <- xest(mod, y0 = Y0)
  var.y0 <- VarCorr(mod)[[1]][1] + VarCorr(mod)[[2]][1]*x0.est^2 + 
    sigma(mod)^2
  ## Function to calculate bootstrap estimate
  bootFun <- function(.) {
    ## Extract model components
    covb <- as.matrix(vcov(.))        # (X' V^-1 X)^-1
    beta.boot <- as.numeric(fixef(.)) 
    ## Calculate estimates
    y0.boot <- rnorm(1, mean = Y0, sd = sqrt(var.y0))
    x0.boot <- xest(., y0 = y0.boot)
    mu0.boot <-  as.numeric(crossprod(beta.boot, c(1, x0.est)))
    ## FIXME: Should these variances be calculated at x0.est or x0.boot?
    var.y0 <- VarCorr(.)[[1]][1] + VarCorr(.)[[2]][1]*x0.est^2 + sigma(.)^2
    var.mu0 <- t(c(1, x0.est)) %*% covb %*% c(1, x0.est)
    Q.boot <- (y0.boot - mu0.boot)/sqrt(var.y0+var.mu0)
    c(x0.boot, Q.boot)
  } 
  ## Function that returns original estimate (i.e., no random y0)
  bootFun0 <- function(.) {
    ## Extract model components
    covb <- as.matrix(vcov(.))        # (X' V^-1 X)^-1
    beta.boot <- as.numeric(fixef(.)) 
    ## Calculate estimates
    y0.boot <- Y0
    x0.boot <- xest(., y0 = y0.boot)
    mu0.boot <-  as.numeric(crossprod(beta.boot, c(1, x0.est)))
    ## FIXME: Should these variances be calculated at x0.est or x0.boot?
    var.y0 <- VarCorr(.)[[1]][1] + VarCorr(.)[[2]][1]*x0.est^2 + sigma(.)^2
    var.mu0 <- t(c(1, x0.est)) %*% covb %*% c(1, x0.est)
    Q.boot <- (y0.boot - mu0.boot)/sqrt(var.y0+var.mu0)
    c(x0.boot, Q.boot)
  }
  ## Generate bootstrap samples
  x0.pb <- if (.parallel) {
    bootMer2_parallel(mod, FUN = bootFun, FUN0 = bootFun0, nsim = R,
                      parallel = "multicore", ncpus = 4)
  } else {
    bootMer2(mod, FUN = bootFun, FUN0 = bootFun0, nsim = R)
  }
  ## Calculate bootstrap CIs
  x0.pb.ci <- suppressWarnings(boot.ci(x0.pb, 
                                       type = c("norm", "basic", "perc")))
  q.quant <- as.numeric(quantile(x0.pb$t[, 2], c(0.025, 0.975)))
  rbind(x0.pb.ci$normal[2:3],
        x0.pb.ci$basic[4:5],
        x0.pb.ci$perc[4:5],
        invCI(.data, q1 = q.quant[1], q2 = q.quant[2], Y0 = Y0))  
}

## Extra -----------------------------------------------------------------------

## Plot and fit models to sample data
# set.seed(8306)
# simdata <- genData()
# p <- xyplot(y ~ x, , data = simdata, type = "n")
# xyplot(y ~ x, groups = subject, data = simdata, type = "b", alpha = 0.75,
#        panel = function(x, y, ...) {
#          panel.xyplot(x, y, ...) 
#          panel.abline(params$beta, lwd = 3)
#          panel.segments(p$x.limits[1], params$y0, params$x0, params$y0)
#          panel.arrows(params$x0, params$y0, params$x0, p$y.limits[1]) 
# })
# # plot(intervals(lmList(y ~ I(x-mean(x)) | subject, data = simdata)))
# mod.lme4 <- lmer(y ~ x + (0+1|subject) + (0+x|subject), data = simdata)
# mod.nlme <- lme(y ~ x, random = list(subject = pdDiag(~x)), data = simdata)
# xest(mod.lme4)
# xest(mod.nlme)

## Test functions on sample data
# res <- rbind(waldCI(simdata, Y0 = 2),
#              invCI(simdata, Y0 = 2),
#              pbootCI(simdata, Y0 = 2))
# colnames(res) <- c("lower", "upper")
# rownames(res) <- c("wald", "inversion", "norm", "basic", "perc", "adjusted")
# cbind(res, length = apply(res, 1, diff))

## Parameter lists for each simulation -----------------------------------------

## Parameters 1
params1 <- list(
  beta = c(0, 1),               # fixed effecs
  theta = c(0.01, 0.05, 0.001), # variance components
  x0 = 0                        # true unknown
)
params1$y0 <- params1$beta[1] + params1$beta[2]*params1$x0 
params1$var.y0 <- params1$theta[1] + params1$theta[2]*params1$x0^2 + 
  params1$theta[3]

## Parameters 2
params2 <- list(
  beta = c(0, 1),               # fixed effecs
  theta = c(0.01, 0.05, 0.001), # variance components
  x0 = 0.25                     # true unknown
)
params2$y0 <- params2$beta[1] + params2$beta[2]*params2$x0 
params2$var.y0 <- params2$theta[1] + params2$theta[2]*params2$x0^2 + 
  params2$theta[3]

## Parameters 3
params3 <- list(
  beta = c(0, 1),               # fixed effecs
  theta = c(0.01, 0.05, 0.001), # variance components
  x0 = 0.5                      # true unknown
)
params3$y0 <- params3$beta[1] + params3$beta[2]*params3$x0 
params3$var.y0 <- params3$theta[1] + params3$theta[2]*params3$x0^2 + 
  params3$theta[3]

## Parameters 4
params4 <- list(
  beta = c(0, 1),               # fixed effecs
  theta = c(0.01, 0.05, 0.001), # variance components
  x0 = 0.75                     # true unknown
)
params4$y0 <- params4$beta[1] + params4$beta[2]*params4$x0 
params4$var.y0 <- params4$theta[1] + params4$theta[2]*params4$x0^2 + 
  params4$theta[3]

## Parameters 5
params5 <- list(
  beta = c(0, 1),               # fixed effecs
  theta = c(0.01, 0.05, 0.001), # variance components
  x0 = 1                        # true unknown
)
params5$y0 <- params5$beta[1] + params5$beta[2]*params5$x0 
params5$var.y0 <- params5$theta[1] + params5$theta[2]*params5$x0^2 + 
  params5$theta[3]

## Simulations -----------------------------------------------------------------

## Simulate data frames
set.seed(5746)
dfs <- rlply(1000, genData)

## Simulations for the Wald-based interval -------------------------------------
# wald.cis1 <- llply(dfs, waldCI, params = params1, .progress = "text")
# wald.cis2 <- llply(dfs, waldCI, params = params2, .progress = "text")
# wald.cis3 <- llply(dfs, waldCI, params = params3, .progress = "text")
# wald.cis4 <- llply(dfs, waldCI, params = params4, .progress = "text")
# wald.cis5 <- llply(dfs, waldCI, params = params5, .progress = "text")
# simulationSummary(wald.cis1, params1)
# simulationSummary(wald.cis2, params2)
# simulationSummary(wald.cis3, params3)
# simulationSummary(wald.cis4, params4)
# simulationSummary(wald.cis5, params5)

## Simulations for the inversion intervals -------------------------------------
# inv.cis1 <- llply(dfs, invCI, params = params1, .progress = "text")
# inv.cis2 <- llply(dfs, invCI, params = params2, .progress = "text")
# inv.cis3 <- llply(dfs, invCI, params = params3, .progress = "text")
# inv.cis4 <- llply(dfs, invCI, params = params4, .progress = "text")
# inv.cis5 <- llply(dfs, invCI, params = params5, .progress = "text")
# simulationSummary(inv.cis1, params1)
# simulationSummary(inv.cis2, params2)
# simulationSummary(inv.cis3, params3)
# simulationSummary(inv.cis4, params4)
# simulationSummary(inv.cis5, params5)

# inv2.cis1 <- llply(dfs, inv2CI, params = params1, .progress = "text")
# inv2.cis2 <- llply(dfs, inv2CI, params = params2, .progress = "text")
# inv2.cis3 <- llply(dfs, inv2CI, params = params3, .progress = "text")
# inv2.cis4 <- llply(dfs, inv2CI, params = params4, .progress = "text")
# inv2.cis5 <- llply(dfs, inv2CI, params = params5, .progress = "text")
# rbind(simulationSummary(inv2.cis1, params1),
#       simulationSummary(inv2.cis2, params2),
#       simulationSummary(inv2.cis3, params3),
#       simulationSummary(inv2.cis4, params4),
#       simulationSummary(inv2.cis5, params5))

## Simulation for the PB percentile intervals (~ 8 hrs) ------------------------
# pb.cis1 <- llply(dfs, pbootCI, params = params1, .progress = "text")
# pb.cis2 <- llply(dfs, pbootCI, params = params2, .progress = "text")
# pb.cis3 <- llply(dfs, pbootCI, params = params3, .progress = "text")
pb.cis4 <- llply(dfs, pbootCI, params = params4, .progress = "text")
# pb.cis5 <- llply(dfs, pbootCI, params = params5, .progress = "text")
# simulationSummary(pb.cis1, params1, boot = TRUE)
# simulationSummary(pb.cis2, params2, boot = TRUE)
# simulationSummary(pb.cis3, params3, boot = TRUE)
simulationSummary(pb.cis4, params4, boot = TRUE)
# simulationSummary(pb.cis5, params5, boot = TRUE)

save(wald.cis1, wald.cis2, wald.cis3, wald.cis4, wald.cis5, inv.cis1, inv.cis2,
     inv.cis3, inv.cis4, inv.cis5, pb.cis1, pb.cis5,
     file = "/home/w108bmg/Desktop/sim_results_linear4.RData")
