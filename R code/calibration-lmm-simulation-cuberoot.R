## 
rm(list = ls())

## Load packages and source code -----------------------------------------------
library(lattice) # for plotting grouped data
library(plyr)    # for working with the data
library(lme4)    # for fitting mixed models and bootstrap inference
library(nlme)    # for fitting mixed models
library(boot)    # for calculating bootstrap confidence intervals
source("/home/w108bmg/Desktop/Dissertation/R code/uniroot2.R")
source("/home/w108bmg/Desktop/Dissertation/R code/Bootstrap functions.R")
trellis.device(color = FALSE)

## Functions -------------------------------------------------------------------

## Cube root function
curt <- function(x) sign(x) * abs(x)^(1/3)

## Function to generate data
genData <- function(n = params$n, m = params$m, beta = params$beta, 
                    theta = params$theta) {
  subject <- rep(1:m, each = n)
  x <- rep(seq(from = 0, to = 1, length = n), times = m) 
  B0 <- rnorm(m, mean = 0, sd = sqrt(theta[1]))
  B1 <- rnorm(m, mean = 0, sd = sqrt(theta[2]))
  y <- rnorm(m*n, mean = beta[1]+B0[subject] + 
               (beta[2]+B1[subject])*x + beta[3]*curt(x), 
             sd = sqrt(theta[3]))
  data.frame(x = x, y = y, subject = factor(subject))
}

## Function to calculate inverse estimate given a mixed model object
xest <- function(object, y0 = params$y0, lower, upper, tol = 1e-10, 
                 maxiter = 1000, extend = "both", ...) {
  fixed <- as.numeric(fixef(object))
  f <- function(x) fixed[1] + fixed[2]*x + fixed[3]*curt(x) - y0
  uniroot2(f, lower = lower, upper = upper, tol = tol, 
           maxiter = maxiter, extend = extend, ...)$root
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

## Function to summarize confidence interval
mc.summary <- function(x, boot = FALSE) {
  
  ## Function for coverting list of C.I.'s into a two-column matrix
  list2Matrix <- function(object) {
    matrix(unlist(object), nrow = length(object), ncol = length(object[[1]]), 
           byrow = TRUE)
  }
  
  ## Function that returns coverage and length of a single C.I.
  coverage.and.length <- function(x) {
    .cov <- if (x[1] <= params$x0 && params$x0 <= x[2]) 1 else 0
    .len <- x[2] - x[1]
    c(.cov, .len)
  }
  
  ## Summarize non bootstrap intervals
  if (!boot) {
    ci.mat <- list2Matrix(x)
    res <- apply(apply(ci.mat, 1, coverage.and.length), 1, mean)
    names(res) <- c("Coverage", "Length")
    res
  ## Summarize bootstrap intervals
  } else {
    d <- ldply(x, function(x) {
      c(coverage.and.length(x[1, ]), # normal
        coverage.and.length(x[2, ]), # basic
        coverage.and.length(x[3, ]), # percentile
        coverage.and.length(x[4, ])) # adjusted inversion
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
wald <- function(.data, Y0) {
  
  ## FIXME: Should this be calculated based on the original model?
  if (missing(Y0)) Y0 <- rnorm(1, mean = params$y0, sd = sqrt(params$var.y0))
  
  ## Fit model using lme4 package and estimate x0 and Var(Y0)
  mod <- lme(y ~ x + curt(x), random = list(subject = pdDiag(~x)), 
             data = .data)  
  x0.est <- xest(mod, y0 = Y0, lower = params$x0-2, upper = params$x0+2)
  var.y0 <- getVarCov(mod)[1, 1] + getVarCov(mod)[2, 2]*x0.est^2 + 
    summary(mod)$sigma^2
  
  ## Function of parameters whose gradient is required
  dmFun <- function(pars) {
    fun <- function(x) (t(c(1, x, curt(x))) %*% pars[-length(pars)]) - 
      pars[length(pars)]
    uniroot2(fun, lower = params$x0-2, upper = params$x0+2, tol = 1e-10, 
            maxiter = 1000, extend = "both")$root
  }

  ## Assign parameter names, calculate gradient, and return standard error
  covmat <- diag(4)
  covmat[1:3, 1:3] <- as.matrix(vcov(mod))
  covmat[4, 4] <- var.y0
  pars <- c(fixef(mod), Y0)
  gv <- attr(numericDeriv(quote(dmFun(pars)), "pars"), "gradient")
  se <- as.numeric(sqrt(gv %*% covmat %*% t(gv)))
  x0.est + qnorm(c(0.025, 0.975))*se
  
}

## Function to calculate the inversion interval
inversion <- function(.data, q1 = qnorm(0.025), q2 = qnorm(0.975), Y0, 
                      lower = params$x0-1, upper = params$x0+1, 
                      ...) {
  
  ## Fit model using nlme package and estimate x0 and Var(Y0)
  mod <- lme(y ~ x + curt(x), random = list(subject = pdDiag(~x)), 
             data = .data)   
  
  ## FIXME: Should this be calculated based on the original model?
  if (missing(Y0)) Y0 <- rnorm(1, mean = params$y0, sd = sqrt(params$var.y0))
  x0.est <- xest(mod, y0 = Y0, lower = lower, upper = upper)
  var.y0 <- getVarCov(mod)[1, 1] + getVarCov(mod)[2, 2]*x0.est^2 + 
    summary(mod)$sigma^2
  
  ## FIXME: can the following code be fixed to avoid errors in the bootstrap
  ##        simulation? The problem seems to occur in finding the roots of fun1 
  ##        and fun2.
  
  ## Prediction function that also returns standard error
  predFun <- function(x) {
    z <- list("x" = x)
    fit <- predict(mod, newdata = z, level = 0)
    se.fit <- sqrt(diag(cbind(1, unlist(z), curt(unlist(z))) %*% mod$varFix %*% 
                          t(cbind(1, unlist(z), curt(unlist(z))))))
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
  
  ## Check curves
  curve(fun1, from = x0.est-10, to = x0.est+10, lwd = 2, ylab = "Bounds", 
        ylim = c(-10, 10))
  curve(fun2, lwd = 2, col = "red", add = TRUE)
  abline(h = 0, v = x0.est, lty = 2)
  
  ## Find roots based on newer bounds
  .lower <- try(uniroot2(fun1, interval = c(lower, x0.est), tol = 1e-10, 
                         maxiter = 1000, extend = "left", ...)$root,
                silent = TRUE)
  .upper <- try(uniroot2(fun2, interval = c(x0.est, upper), tol = 1e-10, 
                         maxiter = 1000, extend = "right", ...)$root, 
                silent = TRUE)
  if (inherits(.lower, "try-error")) .lower <- -Inf
  if (inherits(.upper, "try-error")) .upper <- Inf
  
  ## Find roots of inverse function
  c(.lower, .upper)
  
}

## Function to calculate bootstrap intervals
parboot <- function(.data, R = 999, .parallel = TRUE, Y0) {
  
  ## Fit model using lme4 package and estimate x0 and Var(Y0)
  mod <- lmer(y ~ x + curt(x) + (0+1|subject) + (0+x|subject), 
              data = .data)
  
  ## FIXME: Should this be calculated based on the original model?
  if (missing(Y0)) Y0 <- rnorm(1, mean = params$y0, sd = sqrt(params$var.y0))
  x0.est <- xest(mod, y0 = Y0, lower = params$x0-2, upper =  params$x0+2)
  var.y0 <- VarCorr(mod)[[1]][1] + VarCorr(mod)[[2]][1]*x0.est^2 + 
    sigma(mod)^2
  
  ## Function to calculate bootstrap estimates
  bootFun <- function(.) {

    ## Calculate bootstrap inverse estimate
    if (all(getME(., "y") == .data$y)) {
      y0.boot <- Y0 
    } else {
      y0.boot <- rnorm(1, mean = Y0, sd = sqrt(var.y0))
    }
    x0.boot <- xest(., y0 = y0.boot, lower = params$x0-2, upper =  params$x0+2)

    ## Calculate bootstrap predictive pivot
    ## FIXME: Should variances be calculated at x0.est or x0.boot?
    covb <- as.matrix(vcov(.))        # (X' V^-1 X)^-1
    beta.boot <- as.numeric(fixef(.))
    mu0.boot <-  as.numeric(crossprod(beta.boot, c(1, x0.est, curt(x0.est))))
    var.y0 <- VarCorr(.)[[1]][1] + VarCorr(.)[[2]][1]*x0.est^2 + sigma(.)^2
    var.mu0 <- t(c(1, x0.est, curt(x0.est))) %*% covb %*% 
      c(1, x0.est, curt(x0.est))
    Q.boot <- (y0.boot - mu0.boot)/sqrt(var.y0+var.mu0)
        
    ## Return bootstrap estimates
    c(x0.boot, Q.boot)
    
  }
  
  ## Generate bootstrap samples
  x0.pb <- if (.parallel) {
    bootMer_parallel(mod, FUN = bootFun, nsim = R, parallel = "multicore", 
                     ncpus = 4)
  } else {
    bootMer(mod, FUN = bootFun, nsim = R)
  }

  ## Calculate bootstrap CIs
  x0.pb.ci <- boot.ci(x0.pb, type = c("norm", "basic", "perc"))
  q.quant <- as.numeric(quantile(x0.pb$t[, 2], c(0.025, 0.975), na.rm = TRUE))
  rbind(x0.pb.ci$normal[2:3],
        x0.pb.ci$basic[4:5],
        x0.pb.ci$perc[4:5],
        inversion(.data, q1 = q.quant[1], q2 = q.quant[2], Y0 = Y0))  
  
}

## Simulation parameters -------------------------------------------------------

## Parameters
params <- list(
  nsim = 1000,                  # simulation size
  m = 15,                       # number of subjects
  n = 10,                       # sample size per subject
  beta = c(0, 0.1, 0.9),        # fixed effecs
  theta = c(0.01, 0.05, 0.001), # variance components
  x0 = 0.25                     # true unknown
)
params$y0 <- params$beta[1] + params$beta[2]*params$x0 + 
  params$beta[3]*curt(params$x0)
params$var.y0 <- params$theta[1] + params$theta[2]*params$x0^2 + 
  params$theta[3]

## Sample data and fits
set.seed(101)
simdata <- genData(n = 10, m = 15)
p <- xyplot(y ~ x, , data = simdata, type = "n")
xyplot(y ~ x, groups = subject, data = simdata, type = "l", alpha = 1, 
       lty = 2, panel = function(x, y, ...) {
         panel.xyplot(x, y, ...) 
         panel.curve(params$beta[1] + params$beta[2]*x + params$beta[3]*curt(x), 
                     lwd = 3, from = 0, to = 1)
         panel.segments(p$x.limits[1], params$y0, params$x0, params$y0)
         panel.arrows(params$x0, params$y0, params$x0, p$y.limits[1]) 
       })
fit.nlme <- lme(y ~ x + curt(x), random = list(subject = pdDiag(~x)), 
               data = simdata)
fit.lme4 <- lmer(y ~ x + curt(x) + (0+1|subject) + (0+x|subject), 
                 data = simdata)
est1 <- xest(fit.nlme, lower = -0.5, upper = 0.5)
est2 <- xest(fit.lme4, lower = -0.5, upper = 0.5)
c(est1, est2)
all.equal(est1, est2)
wald(simdata, Y0 = params$y0)
inversion(simdata, Y0 = params$y0)

## Simulations -----------------------------------------------------------------

## Path for savinf results
path <- "/home/w108bmg/Desktop/Simulation results/0.25"

## Simulate data frames
set.seed(5746)
dfs <- rlply(params$nsim, genData)

## Simulation for the Wald-based interval --------------------------------------
mc.wald <- llply(dfs, wald, .progress = "text")
mc.summary(mc.wald)
save(mc.wald, file = paste(path, "mc.wald.RData", sep = "/"))

## Simulation for the inversion interval ---------------------------------------
mc.inversion <- llply(dfs, inversion, .progress = "text")
mc.summary(mc.inversion)
save(mc.inversion, file = paste(path, "mc.inversion.RData", sep = "/"))

## Simulation for the PB intervals (~ 8 hrs) -----------------------------------
# mc.parboot <- llply(dfs, parboot, .progress = "text")
# round(mc.summary(mc.parboot, boot = TRUE), 4)
# save(mc.parboot, file = paste(path, "mc.parboot.RData", sep = "/"))

## Splitting it up is safer!
mc.parboot1 <- llply(dfs[1:100], parboot, .progress = "text")
mc.summary(mc.parboot1, boot = TRUE)
save(mc.parboot1, file = paste(path, "mc.parboot1.RData", sep = "/"))
  
# > mc.summary(mc.parboot1, boot = TRUE)
#            Coverage   Length
# boot.norm      0.87 6.232612
# boot.basic     0.83 3.793693
# boot.perc      0.94 3.793693

mc.parboot2 <- llply(dfs[101:200], parboot, .progress = "text")
mc.summary(mc.parboot2, boot = TRUE)
save(mc.parboot2, file = paste(path, "mc.parboot2.RData", sep = "/"))
     
mc.parboot3 <- llply(dfs[201:300], parboot, .progress = "text") 
mc.summary(mc.parboot3, boot = TRUE)
save(mc.parboot3, file = paste(path, "mc.parboot3.RData", sep = "/"))
     
mc.parboot4 <- llply(dfs[301:400], parboot, .progress = "text") 
mc.summary(mc.parboot4, boot = TRUE)
save(mc.parboot4, file = paste(path, "mc.parboot4.RData", sep = "/"))
     
mc.parboot5 <- llply(dfs[401:500], parboot, .progress = "text")
mc.summary(mc.parboot5, boot = TRUE)
save(mc.parboot5, file = paste(path, "mc.parboot5.RData", sep = "/"))
     
mc.parboot6 <- llply(dfs[501:600], parboot, .progress = "text") 
mc.summary(mc.parboot6, boot = TRUE)
save(mc.parboot6, file = paste(path, "mc.parboot6.RData", sep = "/"))
     
mc.parboot7 <- llply(dfs[601:700], parboot, .progress = "text")
mc.summary(mc.parboot7, boot = TRUE)
save(mc.parboot7, file = paste(path, "mc.parboot7.RData", sep = "/"))
     
mc.parboot8 <- llply(dfs[701:800], parboot, .progress = "text")
mc.summary(mc.parboot8, boot = TRUE)
save(mc.parboot8, file = paste(path, "mc.parboot8.RData", sep = "/"))
     
mc.parboot9 <- llply(dfs[801:900], parboot, .progress = "text")
mc.summary(mc.parboot9, boot = TRUE)
save(mc.parboot9, file = paste(path, "mc.parboot9.RData", sep = "/"))
     
mc.parboot10 <- llply(dfs[901:1000], parboot, .progress = "text")
mc.summary(mc.parboot10, boot = TRUE)
save(mc.parboot10, file = paste(path, "mc.parboot10.RData", sep = "/"))
     
mc.parboot <- c(mc.parboot1, mc.parboot2, mc.parboot3, mc.parboot4, mc.parboot5, 
                mc.parboot6, mc.parboot7, mc.parboot8, mc.parboot9, mc.parboot10)
mc.summary(mc.parboot, boot = TRUE)
save(mc.parboot, file = paste(path, "mc.parboot.RData", sep = "/"))
