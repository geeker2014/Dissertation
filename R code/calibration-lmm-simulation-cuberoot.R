## n-th root mean response curve -----------------------------------------------
root <- 3
f <- function(x) x^(1/root)
xinv <- function(y0) y0^root
curve(f, xlab = "x", ylab = "E(Y|x)")
set.seed(101)
y0 <- rnorm(100000, mean = 0.7368063, sd = 0.1224745)
x0 <- xinv(y0)
hist(y0, br = 50, freq = FALSE, col = "gray", border = "white")
hist(x0, br = 50, freq = FALSE, col = "gray", border = "white", xlim = c(-0.5, 2))
curve(dnorm(x, mean = mean(x0), sd = sd(x0)), lwd = 2, add = TRUE)
lines(density(x0), lwd = 2, lty = 2, col = "red2")
qqnorm(x0)
qqline(x0)

genData <- function(n = 10, m = 15, beta = c(0, 1), 
                    theta = c(0.01, 0.025, 0.001)) {
  subject <- rep(1:m, each = n)
  x <- rep(seq(from = 0, to = 1, length = n), times = m) 
  B0 <- rnorm(m, mean = 0, sd = sqrt(theta[1]))
  B1 <- rnorm(m, mean = 0, sd = sqrt(theta[2]))
  y <- rnorm(m*n, mean = beta[1]+B0[subject] + 
               (beta[2]+B1[subject])*(x^(1/root)), 
             sd = sqrt(theta[3]))
  data.frame(x = x, y = y, subject = factor(subject))
}

xyplot(y ~ x, groups = subject, data = genData(n = 10, m = 1000), type = "l",
       alpha = 0.5)

plot(NULL, xlim = c(0, 1.25), ylim = c(0, 1.5), type = "n")
rect((0.4^(1/3)-3*sqrt(params$var.y0))^3, 0.4^(1/3)-3*sqrt(params$var.y0), 
     (0.4^(1/3)+3*sqrt(params$var.y0))^3, 0.4^(1/3)+3*sqrt(params$var.y0),
     col = "lightgray", border = "gray")
abline(v = 0.4)
abline(h = 0.4^(1/3))
curve(x^(1/3), lwd = 2, add = TRUE)
x0.dens <- density(x0)
lines(x0.dens$x, 0.1*x0.dens$y+1.1)



## Load packages and source code -----------------------------------------------
library(lattice) # for plotting grouped data
library(plyr)    # for working with the data
library(lme4)    # for fitting mixed models and bootstrap inference
library(nlme)    # for fitting mixed models
library(boot)    # for calculating bootstrap confidence intervals
source("/home/w108bmg/Desktop/Dissertation/R code/Parametric bootstrap functions.R")
trellis.device(color = FALSE)
## Simulation setup ------------------------------------------------------------

## Parameters
params <- list(
  nsim = 1000,                    # simulation size
  m = 15,                         # number of subjects
  n = 10,                         # sample size per subject
  beta = c(0, 1),                 # fixed effecs
  theta = c(0.01, 0.025, 0.001), # variance components
  x0 = 0.4                        # true unknown
)
params$y0 <- params$beta[1] + params$beta[2]*params$x0^(1/3)
params$var.y0 <- params$theta[1] + params$theta[2]*params$x0^2 + 
  params$theta[3]

## Function to generate data
genData <- function(n = 10, m = 15, beta = c(0, 1), 
                    theta = c(0.001, 0.015, 0.001)) {
  subject <- rep(1:m, each = n)
  x <- rep(seq(from = 0, to = 1, length = n), times = m) 
  B0 <- rnorm(m, mean = 0, sd = sqrt(theta[1]))
  B1 <- rnorm(m, mean = 0, sd = sqrt(theta[2]))
  y <- rnorm(m*n, mean = beta[1]+B0[subject] + 
               (beta[2]+B1[subject])*(x^(1/root)), 
             sd = sqrt(theta[3]))
  data.frame(x = x, y = y, subject = factor(subject))
}

## Function to calculate inverse estimate given a mixed model object
xest <- function(object, y0 = params$y0) {
  fixed <- as.numeric(fixef(object))
  ((y0 - fixed[1]) / fixed[2])^3
}

## Sample data and fits
set.seed(101)
simdata <- genData(n = 10, m = 15)
p <- xyplot(y ~ x, , data = simdata, type = "n")
xyplot(y ~ x, groups = subject, data = simdata, type = "l", alpha = 1, lty = 2,
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...) 
         panel.curve(x^(1/3), lwd = 3, from = 0, to = 1)
         panel.segments(p$x.limits[1], params$y0, params$x0, params$y0)
         panel.arrows(params$x0, params$y0, params$x0, p$y.limits[1]) 
       })
simdata$x2 <- simdata$x^(1/3)
fit.list <- lmList(y ~ I(x^(1/3) - mean(x^(1/3))) | subject, data = simdata)
fit.nlme <- lme(y ~ I(x^(1/3)), random = list(subject = pdDiag(~I(x^(1/3)))), 
               data = simdata)
fit.lme4 <- lmer(y ~ I(x^(1/3)) + (0+1|subject) + (0+I(x^(1/3))|subject), 
                 data = simdata)
xest(fit.nlme, y0 = 0.75)
plot(intervals(fit.list))


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
  
  if (!boot) {
    
    ci.mat <- list2Matrix(x)
    res <- apply(apply(ci.mat, 1, coverage.and.length), 1, mean)
    names(res) <- c("Coverage", "Length")
    res
    
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
  mod <- lmer(y ~ I(x^(1/3)) + (0+1|subject) + (0+I(x^(1/3))|subject), 
              data = .data)
  x0.est <- xest(mod, y0 = Y0)
  var.y0 <- VarCorr(mod)[[1]][1] + x0.est^2*VarCorr(mod)[[2]][1] + 
    sigma(mod)^2
  
  ## Delta method
  beta <- as.numeric(fixef(mod))
  covmat <- diag(3)
  covmat[1:2, 1:2] <- as.matrix(vcov(mod))
  covmat[3, 3] <- var.y0
  params <- c(beta0 = beta[1], beta1 = beta[2], y0 = Y0)
  gstring <- "((y0 - beta0) / beta1)^3"
  dm <- car:::deltaMethod(params, g = gstring, vcov. = covmat)
  rownames(dm) <- ""
  dm$Estimate + qnorm(c(0.025, 0.975))*dm$SE
  
}

## Function to calculate the inversion interval
invCI <- function(.data, q1 = qnorm(0.025), q2 = qnorm(0.975), Y0) {
  
  ## FIXME: Should this be calculated based on the original model?
  if (missing(Y0)) Y0 <- rnorm(1, mean = params$y0, sd = sqrt(params$var.y0))
  
  ## Fit model using nlme package and estimate x0 and Var(Y0)
  mod <- lme(y ~ I(x^(1/3)), random = list(subject = pdDiag(~I(x^(1/3)))), 
             data = .data)
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
  lower <- uniroot(invFun1, interval = c(0, x0.est), tol = 1e-10, 
                   maxiter = 1000)$root
  upper <- uniroot(invFun2, interval = c(x0.est, 4), tol = 1e-10, 
                   maxiter = 1000)$root
  c(lower, upper)
  
}

## Function to calculate bootstrap intervals
pbootCI <- function(.data, R = 999, .parallel = TRUE, Y0) {
  
  ## FIXME: Should this be calculated based on the original model?
  if (missing(Y0)) Y0 <- rnorm(1, mean = params$y0, sd = sqrt(params$var.y0))
  
  ## Fit model using lme4 package and estimate x0 and Var(Y0)
  mod <- lmer(y ~ I(x^(1/3)) + (0+1|subject) + (0+I(x^(1/3))|subject), 
              data = .data)
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
  x0.pb
  ## Calculate bootstrap CIs
#   x0.pb.ci <- suppressWarnings(boot.ci(x0.pb, 
#                                        type = c("norm", "basic", "perc")))
#   q.quant <- as.numeric(quantile(x0.pb$t[, 2], c(0.025, 0.975)))
#   rbind(x0.pb.ci$normal[2:3],
#         x0.pb.ci$basic[4:5],
#         x0.pb.ci$perc[4:5],
#         #         as.numeric(quantile(x0.pb$t[, 1], c(0.025, 0.975))),
#         invCI(.data, q1 = q.quant[1], q2 = q.quant[2], Y0 = Y0))  
  
}

