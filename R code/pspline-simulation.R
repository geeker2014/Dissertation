################################################################################
##
## Monte Carlo Study for P-spline Calibration
##
################################################################################

## Load libraries --------------------------------------------------------------
library(lattice)
library(investr)
library(plyr)
library(parallel)
source("/home/w108bmg/Desktop/Dissertation-knitr/R code/pspline.R")

## Simulation setup ------------------------------------------------------------

## Simulation parameters
params <- list(
  nsim <- 100,
  f = function(x) sin(pi*x - pi/2)/2 + 0.5,
  m = 3,
  n = 10,
  sd = 0.05,
  x0 = 0.75
)

## Function to simulate scatterplot data
simData <- function(m = params$m, n = params$n, sd = params$sd) {
  x <- rep(seq(from = 0, to = 1, length = params$n), each = params$m)
  y <- params$f(x) + rnorm(params$m*params$n, sd = params$sd)
  data.frame(x = x, y = y)
}

## Function to simulate coverage probability for pspline
simulationSummary <- function(nsim = 10000, degree = degree, 
                              mean.response = FALSE, adjust = TRUE) {
  
  ## Simulate data frames
  dataframes <- rlply(nsim, simData)
  
  ## Function to calculate confidence interval
  fun <- function(df) { 
    fit <- pspline(df, degree = degree)
    Y0 <- f(params$x0) + rnorm(1, sd = params$sd)
    res <- try(invest(fit, y0 = Y0, mean.response = mean.response, 
                      adjust = adjust), silent = TRUE)
    if (inherits(res, "try-error") || any(is.na(res))) {
      repeat {
        newdf <- simData()
        newfit <- pspline(newdf, degree = degree)
        newY0 <- f(params$x0) + rnorm(1, sd = params$sd)
        res <- try(invest(newfit, y0 = newY0, mean.response = mean.response, 
                          adjust = adjust), silent = TRUE)
        if (!inherits(res, "try-error") && !any(is.na(res))) break
      }
    }
    c(if (res["lower"] < params$x0 && res["upper"] > params$x0) 1 else 0,
      abs(res["upper"] - res["lower"]))
  }
  
  z <- matrix(unlist(mclapply(dataframes, fun, mc.cores = 4)), 
              nrow = nsim, ncol = 2, byrow = TRUE)
  colnames(z) <- c("coverage", "length")
  apply(z, 2, mean)
  
}

## Function to simulate coverage probability for linear model
mcCoverage2 <- function(f, m, n, sd = 0.05, nsim = 10000, x0 = 0.75) {
  dataframes <- rlply(nsim, simData(f = f, m = m, n = n, sd = sd))
  fun <- function(df) { 
    Y0 <- f(x0) + rnorm(1, sd = sd)
    res <- calibrate(cbind(f(df$x), df$y), y0 = Y0)
    c(if (res[["lower"]] < f(x0) && res[["upper"]] > f(x0)) 1 else 0,
      abs(res[["upper"]] - res[["lower"]]))
  }
  z <- matrix(unlist(lapply(dataframes, fun)), nrow = nsim, ncol = 2, byrow = TRUE)
  colnames(z) <- c("coverage", "length")
  apply(z, 2, mean)
}

## Function to simulate coverage probability for linear model
mcCoverage4 <- function(f, m, n, sd = 0.05, nsim = 10000, x0 = 0.75) {
  dataframes <- rlply(nsim, simData(f = f, m = m, n = n, sd = sd))
  fun <- function(df) { 
    Y0 <- f(x0)
    res <- calibrate(cbind(f(df$x), df$y), y0 = Y0, mean.response = TRUE)
    c(if (res[["lower"]] < f(x0) && res[["upper"]] > f(x0)) 1 else 0,
      abs(res[["upper"]] - res[["lower"]]))
  }
  z <- matrix(unlist(lapply(dataframes, fun)), nrow = nsim, ncol = 2, byrow = TRUE)
  colnames(z) <- c("coverage", "length")
  apply(z, 2, mean)
}

## Calibration simulation (P-spline) -------------------------------------------
sim1 <- NULL
for (i in 1:3) {
  for (j in c(10, 30, 50)) {
    for (k in 1:3) {
      sim1 <- rbind(sim1, c(mcCoverage(f, m = k, n = j, degree = i), i, j, k))
    }
  }
}

## Covert to data frame and add extra columns
sim1 <- as.data.frame(sim1)
names(sim1) <- c("coverage", "length", "degree", "n", "m")
xyplot(coverage ~ n | as.factor(degree), groups = m, 
       type = "b", auto.key = list(space = "right", title = expression(m)), 
       xlab = "Sample size", ylab = "Coverage probability", data = sim1, 
       main = "P-spline", pch = 19, abline = list(h = 0.95))
#save(sim1, file = "/home/w108bmg/Desktop/Dissertation-knitr/Data/sim1.RData")


## Calibration simulation (LM) -------------------------------------------------
sim2 <- NULL
for (i in c(10, 30, 50)) {
  for (j in 1:3) {
    sim2 <- rbind(sim2, c(mcCoverage2(f, m = j, n = i), i, j))
  }
}

## Covert to data frame and add extra columns
sim2 <- as.data.frame(sim2)
names(sim2) <- c("coverage", "length", "n", "m")
xyplot(coverage ~ n, groups = m, 
       type = "b", auto.key = list(space = "right", title = expression(m)), 
       xlab = "Sample size", ylab = "Coverage probability", data = sim2, 
       main = "Fieller", pch = 19, abline = list(h = 0.95))
#save(sim2, file = "/home/w108bmg/Desktop/Dissertation-knitr/Data/sim2.RData")

## Regulation simulation (P-spline) --------------------------------------------
sim3 <- NULL
for (i in 1:3) {
  for (j in c(10, 30, 50)) {
    for (k in 1:3) {
      sim3 <- rbind(sim3, c(mcCoverage3(f, m = k, n = j, degree = i), i, j, k))
    }
  }
}

## Covert to data frame and add extra columns
sim3 <- as.data.frame(sim3)
names(sim3) <- c("coverage", "length", "degree", "n", "m")
xyplot(coverage ~ n | as.factor(degree), groups = m, 
       type = "b", auto.key = list(space = "right", title = expression(m)), 
       xlab = "Sample size", ylab = "Coverage probability", data = sim3, 
       main = "P-spline", pch = 19, abline = list(h = 0.95))
#save(sim3, file = "/home/w108bmg/Desktop/Dissertation-knitr/Data/sim3.RData")

mcCoverage3(f, m = 1, n = 10, degree = 3)
# > sim3[which(is.na(sim3[,1])),]
#                   co         le d  n m  
# [1,]      0.91970000 0.10777020 3 10 1
# [2,]      0.94430000 0.06065590 3 10 3
# [3,]      0.94050000 0.05985015 3 30 1
# [4,]      0.94620000 0.02597487 3 50 3
# sim3[19,1:2] <- c(0.91970000, 0.10777020)
# sim3[21,1:2] <- c(0.94430000, 0.06065590)
# sim3[22,1:2] <- c(0.94050000, 0.05985015)
# sim3[27,1:2] <- c(0.94620000, 0.02597487)

## Regulation simulation (LM) --------------------------------------------------
sim4 <- NULL
for (i in c(10, 30, 50)) {
  for (j in 1:3) {
    sim4 <- rbind(sim4, c(mcCoverage4(f, m = j, n = i), i, j))
  }
}

## Covert to data frame and add extra columns
sim4 <- as.data.frame(sim4)
names(sim4) <- c("coverage", "length", "n", "m")
xyplot(coverage ~ n, groups = m, 
       type = "b", auto.key = list(space = "right", title = expression(m)), 
       xlab = "Sample size", ylab = "Coverage probability", data = sim4, 
       main = "Fieller", pch = 19, abline = list(h = 0.95))
#save(sim4, file = "/home/w108bmg/Desktop/Dissertation-knitr/Data/sim4.RData")