## n-th root mean response curve -----------------------------------------------
f <- function(x) 0.1*x + 0.9*curt(x)
curve(f, xlab = "x", ylab = "E(Y|x)")
set.seed(101)
y0 <- rnorm(100000, mean = params$y0, sd = sqrt(params$var.y0))
x0 <- numeric(length(y0))
for (i in 1:length(y0)) {
  x0[i] <- uniroot(function(x) { f(x) - y0[i] }, lower = -1, 
                   upper = 2, tol = 1e-10, maxiter = 1000)$root
}
par(mfrow = c(1, 2))
hist(y0, br = 50, freq = FALSE, col = "gray", border = "white")
hist(x0, br = 50, freq = FALSE, col = "gray", border = "white")
# curve(dnorm(x, mean = mean(x0), sd = sd(x0)), lwd = 2, add = TRUE)
# lines(density(x0), lwd = 2, lty = 2, col = "red2")
# qqnorm(x0)
# qqline(x0)





## Load packages and source code -----------------------------------------------
library(lattice) # for plotting grouped data
library(plyr)    # for working with the data
library(lme4)    # for fitting mixed models and bootstrap inference
library(nlme)    # for fitting mixed models
library(boot)    # for calculating bootstrap confidence intervals
source("/home/w108bmg/Desktop/Scripts/R files/compustat/roots.R")
source("/home/w108bmg/Desktop/Dissertation/R code/Parametric bootstrap functions.R")
trellis.device(color = FALSE)
## Simulation setup ------------------------------------------------------------

## Cube root function
curt <- function(x) sign(x) * abs(x)^(1/3)

## Parameters
params <- list(
  nsim = 1000,                    # simulation size
  m = 15,                         # number of subjects
  n = 10,                         # sample size per subject
  beta = c(0, 0.1, 0.9),          # fixed effecs
  theta = c(0.01, 0.05, 0.001), # variance components
  x0 = 0.4                        # true unknown
)
params$y0 <- params$beta[1] + params$beta[2]*params$x0 + 
  params$beta[3]*curt(params$x0)
params$var.y0 <- params$theta[1] + params$theta[2]*params$x0^2 + 
  params$theta[3]

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
                 maxiter = 1000, safe = FALSE, frac = 0.1) {
  .range <- c(lower, upper)
  fixed <- as.numeric(fixef(object))
  f <- function(x) fixed[1] + fixed[2]*x + fixed[3]*curt(x) - y0
  if (safe) {
    if (f(.range[1]) * f(.range[2]) > 0) {
      repeat {
        .range <- extendrange(.range, f = frac)
        if (f(.range[1]) * f(.range[2]) <= 0) break
      }
    }
    uniroot(f, lower = .range[1], upper = .range[2], tol = tol, 
            maxiter = maxiter)$root
  } else {
    uniroot(f, lower = lower, upper = upper, tol = tol, maxiter = maxiter)$root
  }
}

## Function to calculate inverse estimate given a mixed model object based on
## Newton's method
xest_newton <- function(object, y0 = params$y0, start, tol = 1e-10, 
                        maxiter = 1000) {
  fixed <- as.numeric(fixef(object))
  f <- function(x) fixed[1] + fixed[2]*x + fixed[3]*curt(x) - y0
  g <- function(x) fixed[2] + fixed[3]/(3*curt(x^2))
  newton(f, g, start = start, tol = tol, maxiter = maxiter)$root
}

## Function to calculate inverse estimate given a mixed model object based on
## the secant method
xest_secant <- function(object, y0 = params$y0, start, tol = 1e-20, 
                        maxiter = 1000) {
  fixed <- as.numeric(fixef(object))
  f <- function(x) fixed[1] + fixed[2]*x + fixed[3]*curt(x) - y0
  secant(f, start = start, tol = tol, maxiter = maxiter)$root
}

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
est1 <- xest(fit.nlme, lower = 0, upper = 1)
est2 <- xest(fit.lme4, lower = 0, upper = 1)
est3 <- xest_secant(fit.nlme, start = c(0, 2))
est4 <- xest_secant(fit.lme4, start = c(0, 2))
c(est1, est2, est3, est4)
all.equal(est1, est2, est3, est4)

# fit.list <- lmList(y ~ x + curt(x) | subject, data = simdata)
# plot(intervals(fit.list))

## Simulated data frames
set.seed(5746)
dfs <- rlply(params$nsim, genData)

## Function to summarize confidence interval
simulationSummary <- function(x, boot = FALSE) {
  
  ## Function for coverting list of C.I.'s into a two-column matrix
  list2Matrix <- function(object) {
    matrix(unlist(object), nrow = length(object), ncol = length(object[[1]]), 
           byrow = TRUE)
  }
  
  ## Function that returns coverage and length of a single C.I.
  coverage.and.length <- function(x) {
    .coverage <- if (x[1] <= params$x0 && params$x0 <= x[2]) 1 else 0
    .length <- x[2] - x[1]
    c(.coverage, .length)
  }
  
  ## Summarize non bootstrap intervals
  if (!boot) {
    
    ci.mat <- list2Matrix(x)
    res <- apply(apply(ci.mat, 1, coverage.and.length), 1, mean)
    names(res) <- c("Coverage", "Length")
    res
    
  ## Summarize bootstrap intervals
  } else {
    
    ## Data frame of coverge and length estimates
    d <- ldply(x, function(x) {
      c(coverage.and.length(x[1, ]), # normal
        coverage.and.length(x[2, ]), # basic
        coverage.and.length(x[3, ]), # percentile
        coverage.and.length(x[4, ])) # adjusted inversion
    })
    
    ## Calculate mean coverage and mean length
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
waldCI <- function(.data, Y0) {
  
  ## FIXME: Should this be calculated based on the original model?
  if (missing(Y0)) Y0 <- rnorm(1, mean = params$y0, sd = sqrt(params$var.y0))
  
  ## Fit model using lme4 package and estimate x0 and Var(Y0)
  mod <- lme(y ~ x + curt(x), random = list(subject = pdDiag(~x)), 
             data = simdata)  
  x0.est <- xest(mod, y0 = Y0, lower = -2, upper = 3)
  var.y0 <- getVarCov(mod)[1, 1] + getVarCov(mod)[2, 2]*x0.est^2 + 
    summary(mod)$sigma^2
  
  ## Function of parameters whose gradient is required
  dmFun <- function(pars) {
    fun <- function(x) (t(c(1, x, curt(x))) %*% pars[-length(pars)]) - 
      pars[length(pars)]
    uniroot(fun, lower = -2, upper = 3, tol = 1e-10, maxiter = 1000)$root
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
invCI <- function(.data, q1 = qnorm(0.025), q2 = qnorm(0.975), Y0, safe = FALSE,
                  lower = -2, upper = 3) {
  
  ## FIXME: Should this be calculated based on the original model?
  if (missing(Y0)) Y0 <- rnorm(1, mean = params$y0, sd = sqrt(params$var.y0))
  
  ## Fit model using nlme package and estimate x0 and Var(Y0)
  mod <- lme(y ~ x + curt(x), random = list(subject = pdDiag(~x)), 
             data = simdata)   
  x0.est <- xest(mod, y0 = Y0, lower = -2, upper = 3, safe = TRUE)
  var.y0 <- getVarCov(mod)[1, 1] + getVarCov(mod)[2, 2]*x0.est^2 + 
    summary(mod)$sigma^2
  
  ## FIXME: can the following code be fixed to avoid errors in the bootstrap
  ## simulation? The problem seems to occur in finding the roots of fun1 and 
  ## fun2.
  
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
  curve(fun1, from = x0.est-1, to = x0.est+2, lwd = 2)
  curve(fun2, lwd = 2, col = "red", add = TRUE)
  abline(h = 0, v = x0.est, lty = 2)
  
  ## Safe version
  if (safe) {
    
    ## Original lower and upper bounds
    .lwr <- lower
    .upr <- upper
    
    ## Extend lower bound until root is contained
    if (fun1(.lwr) * fun1(x0.est) > 0) {
      repeat {
        .lwr <- .lwr - 0.1
        if (fun1(.lwr) * fun1(x0.est) <= 0) break
      }
    }
    
    ## Extend upper bound until root is contained
    if (fun2(x0.est) * fun2(.upr) > 0) {
      repeat {
        .upr <- .upr + 0.1
        if (fun2(x0.est) * fun2(.upr) <= 0) break
      }
    }
    
    ## Find roots based on newer bounds
    .lower <- uniroot(fun1, interval = c(.lwr, x0.est), tol = 1e-10, 
                      maxiter = 1000)$root
    .upper <- uniroot(fun2, interval = c(x0.est, .upr), tol = 1e-10, 
                      maxiter = 1000)$root
    
  ## Unsafe version
  } else {
    
    ## Find roots based on fixed bounds
    .lower <- uniroot(fun1, interval = c(lower, x0.est), tol = 1e-10, 
                      maxiter = 1000)$root
    .upper <- uniroot(fun2, interval = c(x0.est, upper), tol = 1e-10, 
                      maxiter = 1000)$root
    
  }
  
  ## Find roots of inverse function
  c(.lower, .upper)
  
}

## Function to calculate bootstrap intervals
pbootCI <- function(.data, R = 999, .parallel = TRUE, Y0) {
  
  ## FIXME: Should this be calculated based on the original model?
  if (missing(Y0)) Y0 <- rnorm(1, mean = params$y0, sd = sqrt(params$var.y0))
  
  ## Fit model using lme4 package and estimate x0 and Var(Y0)
  mod <- lmer(y ~ x + curt(x) + (0+1|subject) + (0+x|subject), 
              data = .data)
  x0.est <- xest(mod, y0 = Y0, lower = -2, upper = 3, safe = TRUE)
#   x0.est <- xest_secant(mod, y0 = Y0, start = c(0, 2)) # use secant method
  var.y0 <- VarCorr(mod)[[1]][1] + VarCorr(mod)[[2]][1]*x0.est^2 + 
    sigma(mod)^2
  
  ## Function to calculate bootstrap estimate
  bootFun <- function(.) {
    
    ## Extract model components
    covb <- as.matrix(vcov(.))        # (X' V^-1 X)^-1
    beta.boot <- as.numeric(fixef(.)) 
    
    ## Calculate estimates
    y0.boot <- rnorm(1, mean = Y0, sd = sqrt(var.y0))
    x0.boot <- xest(., y0 = y0.boot, lower = -2, upper = 3, safe = TRUE)
#     x0.boot <- xest_secant(., y0 = y0.boot, start = c(0, 2))
    mu0.boot <-  as.numeric(crossprod(beta.boot, c(1, x0.est, curt(x0.est))))
    
    ## FIXME: Should variances be calculated at x0.est or x0.boot?
    var.y0 <- VarCorr(.)[[1]][1] + VarCorr(.)[[2]][1]*x0.est^2 + sigma(.)^2
    var.mu0 <- t(c(1, x0.est, curt(x0.est))) %*% covb %*% 
      c(1, x0.est, curt(x0.est))
    Q.boot <- (y0.boot - mu0.boot)/sqrt(var.y0+var.mu0)
    
    ## Return bootstrap estimates
    c(x0.boot, Q.boot)
    
  } 
  
  ## Function that returns original estimate (i.e., no random y0)
  bootFun0 <- function(.) {
    
    ## Extract model components
    covb <- as.matrix(vcov(.))        # (X' V^-1 X)^-1
    beta.boot <- as.numeric(fixef(.)) 
    
    ## Calculate estimates
    y0.boot <- Y0
    x0.boot <- xest(., y0 = y0.boot, lower = -2, upper = 3, safe = TRUE)
#     x0.boot <- xest_secant(., y0 = y0.boot, start = c(0, 2))
    mu0.boot <-  as.numeric(crossprod(beta.boot, c(1, x0.est, curt(x0.est))))
    
    ## FIXME: Should variances be calculated at x0.est or x0.boot?
    var.y0 <- VarCorr(.)[[1]][1] + VarCorr(.)[[2]][1]*x0.est^2 + sigma(.)^2
    var.mu0 <- t(c(1, x0.est, curt(x0.est))) %*% covb %*% 
      c(1, x0.est, curt(x0.est))
    Q.boot <- (y0.boot - mu0.boot)/sqrt(var.y0+var.mu0)
    
    ## Return bootstrap estimates
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
  q.quant <- as.numeric(quantile(x0.pb$t[, 2], c(0.025, 0.975), na.rm = TRUE))
  rbind(x0.pb.ci$normal[2:3],
        x0.pb.ci$basic[4:5],
        x0.pb.ci$perc[4:5],
        invCI(.data, q1 = q.quant[1], q2 = q.quant[2], Y0 = Y0, safe = TRUE))  
  
}

## Test CI functions -----------------------------------------------------------
res <- rbind(waldCI(simdata, Y0 = params$y0),
             invCI(simdata, Y0 = params$y0),
             pbootCI(simdata, Y0 = params$y0))
colnames(res) <- c("Coverage", "Length")
rownames(res) <- c("Wald", "Inversion", "Bootstrap - norm", 
                   "Bootstrap - basic", "Bootstrap - perc", "Bootstrap - inv")
res

# > res
# Coverage    Length
# Wald               0.003440265 0.8587347
# Inversion          0.126357431 1.0422462
# Bootstrap - norm  -0.104777797 0.8793157
# Bootstrap - basic -0.179456688 0.7445643
# Bootstrap - perc   0.117610714 1.0416317
# Bootstrap - inv    0.119380259 1.0948971

## Simulation for the Wald-based interval --------------------------------------
wald.cis <- llply(dfs, waldCI, .progress = "text")
round(simulationSummary(wald.cis), 4)

## Simulation for the inversion interval ---------------------------------------
inv.cis <- llply(dfs, invCI, .progress = "text")
round(simulationSummary(inv.cis), 4)

## Simulation for the PB intervals (~ 8 hrs) -----------------------------------
pb.cis <- llply(dfs, pbootCI, .progress = "text")
round(simulationSummary(pb.cis, boot = TRUE), 4)
system.time(pbootCI(simdata, Y0 = params$y0))

## Try splitting it up
pb.cis1 <- llply(dfs[1:100], pbootCI, .progress = "text")
pb.cis2 <- llply(dfs[101:200], pbootCI, .progress = "text")
pb.cis3 <- llply(dfs[201:300], pbootCI, .progress = "text")
pb.cis4 <- llply(dfs[301:400], pbootCI, .progress = "text")
pb.cis5 <- llply(dfs[401:500], pbootCI, .progress = "text")
pb.cis6 <- llply(dfs[501:600], pbootCI, .progress = "text")
pb.cis7 <- llply(dfs[601:700], pbootCI, .progress = "text")
pb.cis8 <- llply(dfs[701:800], pbootCI, .progress = "text")
pb.cis9 <- llply(dfs[801:900], pbootCI, .progress = "text")
pb.cis10 <- llply(dfs[901:1000], pbootCI, .progress = "text")

# pb.cis <- rbind(pb.cis1, ..., pb.cis10)
round(simulationSummary(pb.cis1, boot = TRUE), 4)
