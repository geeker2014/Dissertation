## Simulation study for nonparametric calibration using penalized spline
## regression

## Load packages
library(lattice)
library(investr)
library(plyr)
library(parallel)

##==============================================================================
##
## Monte Carlo study (using parallel backend)
##
##==============================================================================

f <- function(x) sin(pi*x - pi/2)/2 + 0.5

## Logistic function
f1 <- function(x, theta = c(1, 1, 7, 0.5)) {
  theta[1] / (theta[2] + exp(-theta[3]*(x - theta[4])))
}
# curve(f1, las = 1, xlab = expression(x), ylab = expression(f(x)))

## Sine wave
f2 <- function(x) {
  sin(pi*x - pi/2)/2 + 0.5
}
# curve(f2, las = 1, xlab = expression(x), ylab = expression(f(x)))

f3 <- function(x, Asym = 1, R0 = 0, lrc = 1.5) {
  SSasymp(x, Asym = Asym, R0 = R0, lrc = lrc)
}
# curve(f3, las = 1, xlab = expression(x), ylab = expression(f(x)))

f4 <- function(x, theta = c(0.1, -0.5)) {
  1 - (x - 1)^2
}
# curve(f4, las = 1, xlab = expression(x), ylab = expression(f(x)))

## Function to simulate scatterplot data
simData <- function(f, m = 3, n = 10, sd = 0.05, from = 0, to = 1) {
  x <- rep(seq(from = from, to = to, length = n), each = m)
  y <- f(x) + rnorm(m*n, sd = sd)
  data.frame(x = x, y = y)
}

## Simulation function

#
# FIXME: x0 and y0 should not both be 0.5!
#

mcCoverage <- function(f, m, n, sd = 0.05, num.sim = 5000, x0 = 0.5, degree, ...) 
{
  dataframes <- rlply(num.sim, simData(f = f, m = m, n = n, sd = sd))
  fun <- function(df) {
    fit <- with(df, lmmSmooth(x, y, degree = degree))
    Y0 <- f(x0) + rnorm(1, sd = sd)
    res <- tryCatch(invest(fit, y0 = Y0, ...), 
                    error = function(e) NULL)
    if (is.null(res)) {
      repeat {
        newdf <- simData(f = f, m = m, n = n, sd = sd)
        newfit <- with(newdf, lmmSmooth(x, y, degree = degree))
        newY0 <- f(x0) + rnorm(1, sd = sd)
        res <- tryCatch(invest(newfit, y0 = newY0, ...), 
                        error = function(e) NULL)
        if (!is.null(res)) break
      }
    }
    c(if (res["lower"] < x0 && res["upper"] > x0) 1 else 0,
      abs(res["upper"] - res["lower"]))
  }
  z <- matrix(unlist(mclapply(dataframes, fun, mc.cores = 4)), 
              nrow = num.sim, ncol = 2, byrow = TRUE)
  colnames(z) <- c("coverage", "length")
  apply(z, 2, mean)
}

## For mean.response = TRUE
mcCoverage2 <- function(f, m, n, sd = 0.05, num.sim = 5000, x0 = 0.5, degree, ...) 
{
  dataframes <- rlply(num.sim, simData(f = f, m = m, n = n, sd = sd))
  fun <- function(df) {
    fit <- with(df, lmmSmooth(x, y, degree = degree))
    res <- tryCatch(invest(fit, y0 = f(x0), mean.response = TRUE, ...), 
                    error = function(e) NULL)
    if (is.null(res)) {
      repeat {
        newdf <- simData(f = f, m = m, n = n, sd = sd)
        newfit <- with(newdf, lmmSmooth(x, y, degree = degree))
        res <- tryCatch(invest(newfit, y0 = f(x0), mean.response = TRUE, ...), 
                        error = function(e) NULL)
        if (!is.null(res)) break
      }
    }
    c(if (res["lower"] < x0 && res["upper"] > x0) 1 else 0,
      abs(res["upper"] - res["lower"]))
  }
  z <- matrix(unlist(mclapply(dataframes, fun, mc.cores = 4)), 
              nrow = num.sim, ncol = 2, byrow = TRUE)
  colnames(z) <- c("coverage", "length")
  apply(z, 2, mean)
}

## Simulation for logistic function --------------------------------------------
plot(simData(f1), las = 1, xlab = expression(x), ylab = expression(f(x)))
curve(f1, lwd = 2, add = T)

sim.f1 <- matrix(0, nrow = 36, ncol = 2)
sim.f1[1, ]  <- mcCoverage(f1, m = 1, n = 10, degree = 1, corrected = T)
sim.f1[2, ]  <- mcCoverage(f1, m = 2, n = 10, degree = 1, corrected = T)
sim.f1[3, ]  <- mcCoverage(f1, m = 3, n = 10, degree = 1, corrected = T)
sim.f1[4, ]  <- mcCoverage(f1, m = 1, n = 30, degree = 1, corrected = T)
sim.f1[5, ]  <- mcCoverage(f1, m = 2, n = 30, degree = 1, corrected = T)
sim.f1[6, ]  <- mcCoverage(f1, m = 3, n = 30, degree = 1, corrected = T)
sim.f1[7, ]  <- mcCoverage(f1, m = 1, n = 50, degree = 1, corrected = T)
sim.f1[8, ]  <- mcCoverage(f1, m = 2, n = 50, degree = 1, corrected = T)
sim.f1[9, ]  <- mcCoverage(f1, m = 3, n = 50, degree = 1, corrected = T)
sim.f1[10, ] <- mcCoverage(f1, m = 1, n = 10, degree = 1, corrected = F)
sim.f1[11, ] <- mcCoverage(f1, m = 2, n = 10, degree = 1, corrected = F)
sim.f1[12, ] <- mcCoverage(f1, m = 3, n = 10, degree = 1, corrected = F)
sim.f1[13, ] <- mcCoverage(f1, m = 1, n = 30, degree = 1, corrected = F)
sim.f1[14, ] <- mcCoverage(f1, m = 2, n = 30, degree = 1, corrected = F)
sim.f1[15, ] <- mcCoverage(f1, m = 3, n = 30, degree = 1, corrected = F)
sim.f1[16, ] <- mcCoverage(f1, m = 1, n = 50, degree = 1, corrected = F)
sim.f1[17, ] <- mcCoverage(f1, m = 2, n = 50, degree = 1, corrected = F)
sim.f1[18, ] <- mcCoverage(f1, m = 3, n = 50, degree = 1, corrected = F)

sim.f1[19, ] <- mcCoverage(f1, m = 1, n = 10, degree = 2, corrected = T)
sim.f1[20, ] <- mcCoverage(f1, m = 2, n = 10, degree = 2, corrected = T)
sim.f1[21, ] <- mcCoverage(f1, m = 3, n = 10, degree = 2, corrected = T)
sim.f1[22, ] <- mcCoverage(f1, m = 1, n = 30, degree = 2, corrected = T)
sim.f1[23, ] <- mcCoverage(f1, m = 2, n = 30, degree = 2, corrected = T)
sim.f1[24, ] <- mcCoverage(f1, m = 3, n = 30, degree = 2, corrected = T)
sim.f1[25, ] <- mcCoverage(f1, m = 1, n = 50, degree = 2, corrected = T)
sim.f1[26, ] <- mcCoverage(f1, m = 2, n = 50, degree = 2, corrected = T)
sim.f1[27, ] <- mcCoverage(f1, m = 3, n = 50, degree = 2, corrected = T)
sim.f1[28, ] <- mcCoverage(f1, m = 1, n = 10, degree = 2, corrected = F)
sim.f1[29, ] <- mcCoverage(f1, m = 2, n = 10, degree = 2, corrected = F)
sim.f1[30, ] <- mcCoverage(f1, m = 3, n = 10, degree = 2, corrected = F)
sim.f1[31, ] <- mcCoverage(f1, m = 1, n = 30, degree = 2, corrected = F)
sim.f1[32, ] <- mcCoverage(f1, m = 2, n = 30, degree = 2, corrected = F)
sim.f1[33, ] <- mcCoverage(f1, m = 3, n = 30, degree = 2, corrected = F)
sim.f1[34, ] <- mcCoverage(f1, m = 1, n = 50, degree = 2, corrected = F)
sim.f1[35, ] <- mcCoverage(f1, m = 2, n = 50, degree = 2, corrected = F)
sim.f1[36, ] <- mcCoverage(f1, m = 3, n = 50, degree = 2, corrected = F)

## Covert to data frame and add extra columns
sim.f1 <- as.data.frame(sim.f1)
names(sim.f1) <- c("coverage", "length")
sim.f1$m <- rep(1:3, times = 12)
sim.f1$n <- rep(rep(c(10, 30, 50), each = 3), times = 4)
sim.f1$degree <- rep(1:2, each = 18)
sim.f1$corrected <- rep(rep(c("yes", "no"), each = 9), times = 2)
xyplot(coverage ~ n | as.factor(degree)*as.factor(corrected), groups = m, 
       type = "b", auto.key = list(space = "right", title = expression(m)), 
       xlab = "Coverage probability", ylab = "Sample size", data = sim.f1, 
       pch = 19, abline = list(h = 0.95))
write.csv(sim.f1, "/home/w108bmg/Desktop/Dissertation-knitr/R code/Simulation data/sim.f1.csv", row.names = F)

## Simulation for sine wave ----------------------------------------------------
plot(simData(f2), las = 1, xlab = expression(x), ylab = expression(f(x)))
curve(f2, lwd = 2, add = T)

sim.f2 <- matrix(0, nrow = 36, ncol = 2)
sim.f2[1, ]  <- mcCoverage(f2, m = 1, n = 10, degree = 1, corrected = T)
sim.f2[2, ]  <- mcCoverage(f2, m = 2, n = 10, degree = 1, corrected = T)
sim.f2[3, ]  <- mcCoverage(f2, m = 3, n = 10, degree = 1, corrected = T)
sim.f2[4, ]  <- mcCoverage(f2, m = 1, n = 30, degree = 1, corrected = T)
sim.f2[5, ]  <- mcCoverage(f2, m = 2, n = 30, degree = 1, corrected = T)
sim.f2[6, ]  <- mcCoverage(f2, m = 3, n = 30, degree = 1, corrected = T)
sim.f2[7, ]  <- mcCoverage(f2, m = 1, n = 50, degree = 1, corrected = T)
sim.f2[8, ]  <- mcCoverage(f2, m = 2, n = 50, degree = 1, corrected = T)
sim.f2[9, ]  <- mcCoverage(f2, m = 3, n = 50, degree = 1, corrected = T)
sim.f2[10, ] <- mcCoverage(f2, m = 1, n = 10, degree = 1, corrected = F)
sim.f2[11, ] <- mcCoverage(f2, m = 2, n = 10, degree = 1, corrected = F)
sim.f2[12, ] <- mcCoverage(f2, m = 3, n = 10, degree = 1, corrected = F)
sim.f2[13, ] <- mcCoverage(f2, m = 1, n = 30, degree = 1, corrected = F)
sim.f2[14, ] <- mcCoverage(f2, m = 2, n = 30, degree = 1, corrected = F)
sim.f2[15, ] <- mcCoverage(f2, m = 3, n = 30, degree = 1, corrected = F)
sim.f2[16, ] <- mcCoverage(f2, m = 1, n = 50, degree = 1, corrected = F)
sim.f2[17, ] <- mcCoverage(f2, m = 2, n = 50, degree = 1, corrected = F)
sim.f2[18, ] <- mcCoverage(f2, m = 3, n = 50, degree = 1, corrected = F)

sim.f2[19, ] <- mcCoverage(f2, m = 1, n = 10, degree = 2, corrected = T)
sim.f2[20, ] <- mcCoverage(f2, m = 2, n = 10, degree = 2, corrected = T)
sim.f2[21, ] <- mcCoverage(f2, m = 3, n = 10, degree = 2, corrected = T)
sim.f2[22, ] <- mcCoverage(f2, m = 1, n = 30, degree = 2, corrected = T)
sim.f2[23, ] <- mcCoverage(f2, m = 2, n = 30, degree = 2, corrected = T)
sim.f2[24, ] <- mcCoverage(f2, m = 3, n = 30, degree = 2, corrected = T)
sim.f2[25, ] <- mcCoverage(f2, m = 1, n = 50, degree = 2, corrected = T)
sim.f2[26, ] <- mcCoverage(f2, m = 2, n = 50, degree = 2, corrected = T)
sim.f2[27, ] <- mcCoverage(f2, m = 3, n = 50, degree = 2, corrected = T)
sim.f2[28, ] <- mcCoverage(f2, m = 1, n = 10, degree = 2, corrected = F)
sim.f2[29, ] <- mcCoverage(f2, m = 2, n = 10, degree = 2, corrected = F)
sim.f2[30, ] <- mcCoverage(f2, m = 3, n = 10, degree = 2, corrected = F)
sim.f2[31, ] <- mcCoverage(f2, m = 1, n = 30, degree = 2, corrected = F)
sim.f2[32, ] <- mcCoverage(f2, m = 2, n = 30, degree = 2, corrected = F)
sim.f2[33, ] <- mcCoverage(f2, m = 3, n = 30, degree = 2, corrected = F)
sim.f2[34, ] <- mcCoverage(f2, m = 1, n = 50, degree = 2, corrected = F)
sim.f2[35, ] <- mcCoverage(f2, m = 2, n = 50, degree = 2, corrected = F)
sim.f2[36, ] <- mcCoverage(f2, m = 3, n = 50, degree = 2, corrected = F)

## Covert to data frame and add extra columns
sim.f2 <- as.data.frame(sim.f2)
names(sim.f2) <- c("coverage", "length")
sim.f2$m <- rep(1:3, times = 12)
sim.f2$n <- rep(rep(c(10, 30, 50), each = 3), times = 4)
sim.f2$degree <- rep(1:2, each = 18)
sim.f2$corrected <- rep(rep(c("yes", "no"), each = 9), times = 2)
xyplot(coverage ~ n | as.factor(degree)*as.factor(corrected), groups = m, 
       type = "b", auto.key = list(space = "right", title = expression(m)), 
       xlab = "Coverage probability", ylab = "Sample size", data = sim.f2, 
       pch = 19, abline = list(h = 0.95))
write.csv(sim.f2, "/home/w108bmg/Desktop/Dissertation-knitr/R code/Simulation data/sim.f2.csv", row.names = F)

## Asymptotic ------------------------------------------------------------------
plot(simData(f3), las = 1, xlab = expression(x), ylab = expression(f(x)))
curve(f3, lwd = 2, add = T)

sim.f3 <- matrix(0, nrow = 36, ncol = 2)
sim.f3[1, ]  <- mcCoverage(f3, m = 1, n = 10, degree = 1, corrected = T)
sim.f3[2, ]  <- mcCoverage(f3, m = 2, n = 10, degree = 1, corrected = T)
sim.f3[3, ]  <- mcCoverage(f3, m = 3, n = 10, degree = 1, corrected = T)
sim.f3[4, ]  <- mcCoverage(f3, m = 1, n = 30, degree = 1, corrected = T)
sim.f3[5, ]  <- mcCoverage(f3, m = 2, n = 30, degree = 1, corrected = T)
sim.f3[6, ]  <- mcCoverage(f3, m = 3, n = 30, degree = 1, corrected = T)
sim.f3[7, ]  <- mcCoverage(f3, m = 1, n = 50, degree = 1, corrected = T)
sim.f3[8, ]  <- mcCoverage(f3, m = 2, n = 50, degree = 1, corrected = T)
sim.f3[9, ]  <- mcCoverage(f3, m = 3, n = 50, degree = 1, corrected = T)
sim.f3[10, ] <- mcCoverage(f3, m = 1, n = 10, degree = 1, corrected = F)
sim.f3[11, ] <- mcCoverage(f3, m = 2, n = 10, degree = 1, corrected = F)
sim.f3[12, ] <- mcCoverage(f3, m = 3, n = 10, degree = 1, corrected = F)
sim.f3[13, ] <- mcCoverage(f3, m = 1, n = 30, degree = 1, corrected = F)
sim.f3[14, ] <- mcCoverage(f3, m = 2, n = 30, degree = 1, corrected = F)
sim.f3[15, ] <- mcCoverage(f3, m = 3, n = 30, degree = 1, corrected = F)
sim.f3[16, ] <- mcCoverage(f3, m = 1, n = 50, degree = 1, corrected = F)
sim.f3[17, ] <- mcCoverage(f3, m = 2, n = 50, degree = 1, corrected = F)
sim.f3[18, ] <- mcCoverage(f3, m = 3, n = 50, degree = 1, corrected = F)

sim.f3[19, ] <- mcCoverage(f3, m = 1, n = 10, degree = 2, corrected = T)
sim.f3[20, ] <- mcCoverage(f3, m = 2, n = 10, degree = 2, corrected = T)
sim.f3[21, ] <- mcCoverage(f3, m = 3, n = 10, degree = 2, corrected = T)
sim.f3[22, ] <- mcCoverage(f3, m = 1, n = 30, degree = 2, corrected = T)
sim.f3[23, ] <- mcCoverage(f3, m = 2, n = 30, degree = 2, corrected = T)
sim.f3[24, ] <- mcCoverage(f3, m = 3, n = 30, degree = 2, corrected = T)
sim.f3[25, ] <- mcCoverage(f3, m = 1, n = 50, degree = 2, corrected = T)
sim.f3[26, ] <- mcCoverage(f3, m = 2, n = 50, degree = 2, corrected = T)
sim.f3[27, ] <- mcCoverage(f3, m = 3, n = 50, degree = 2, corrected = T)
sim.f3[28, ] <- mcCoverage(f3, m = 1, n = 10, degree = 2, corrected = F)
sim.f3[29, ] <- mcCoverage(f3, m = 2, n = 10, degree = 2, corrected = F)
sim.f3[30, ] <- mcCoverage(f3, m = 3, n = 10, degree = 2, corrected = F)
sim.f3[31, ] <- mcCoverage(f3, m = 1, n = 30, degree = 2, corrected = F)
sim.f3[32, ] <- mcCoverage(f3, m = 2, n = 30, degree = 2, corrected = F)
sim.f3[33, ] <- mcCoverage(f3, m = 3, n = 30, degree = 2, corrected = F)
sim.f3[34, ] <- mcCoverage(f3, m = 1, n = 50, degree = 2, corrected = F)
sim.f3[35, ] <- mcCoverage(f3, m = 2, n = 50, degree = 2, corrected = F)
sim.f3[36, ] <- mcCoverage(f3, m = 3, n = 50, degree = 2, corrected = F)

## Covert to data frame and add extra columns
sim.f3 <- as.data.frame(sim.f3)
names(sim.f3) <- c("coverage", "length")
sim.f3$m <- rep(1:3, times = 12)
sim.f3$n <- rep(rep(c(10, 30, 50), each = 3), times = 4)
sim.f3$degree <- rep(1:2, each = 18)
sim.f3$corrected <- rep(rep(c("yes", "no"), each = 9), times = 2)
xyplot(coverage ~ n | as.factor(degree)*as.factor(corrected), groups = m, 
       type = "b", auto.key = list(space = "right", title = expression(m)), 
       xlab = "Coverage probability", ylab = "Sample size", data = sim.f3, 
       pch = 19, abline = list(h = 0.95))
write.csv(sim.f3, "/home/w108bmg/Desktop/Dissertation-knitr/R code/Simulation data/sim.f3.csv", row.names = F)

## Simulation for quadratic ----------------------------------------------------
plot(simData(f4), las = 1, xlab = expression(x), ylab = expression(f(x)))
curve(f4, lwd = 2, add = T)

sim.f4 <- matrix(0, nrow = 36, ncol = 2)
sim.f4[1, ]  <- mcCoverage(f4, m = 1, n = 10, degree = 1, corrected = T)
sim.f4[2, ]  <- mcCoverage(f4, m = 2, n = 10, degree = 1, corrected = T)
sim.f4[3, ]  <- mcCoverage(f4, m = 3, n = 10, degree = 1, corrected = T)
sim.f4[4, ]  <- mcCoverage(f4, m = 1, n = 30, degree = 1, corrected = T)
sim.f4[5, ]  <- mcCoverage(f4, m = 2, n = 30, degree = 1, corrected = T)
sim.f4[6, ]  <- mcCoverage(f4, m = 3, n = 30, degree = 1, corrected = T)
sim.f4[7, ]  <- mcCoverage(f4, m = 1, n = 50, degree = 1, corrected = T)
sim.f4[8, ]  <- mcCoverage(f4, m = 2, n = 50, degree = 1, corrected = T)
sim.f4[9, ]  <- mcCoverage(f4, m = 3, n = 50, degree = 1, corrected = T)
sim.f4[10, ] <- mcCoverage(f4, m = 1, n = 10, degree = 1, corrected = F)
sim.f4[11, ] <- mcCoverage(f4, m = 2, n = 10, degree = 1, corrected = F)
sim.f4[12, ] <- mcCoverage(f4, m = 3, n = 10, degree = 1, corrected = F)
sim.f4[13, ] <- mcCoverage(f4, m = 1, n = 30, degree = 1, corrected = F)
sim.f4[14, ] <- mcCoverage(f4, m = 2, n = 30, degree = 1, corrected = F)
sim.f4[15, ] <- mcCoverage(f4, m = 3, n = 30, degree = 1, corrected = F)
sim.f4[16, ] <- mcCoverage(f4, m = 1, n = 50, degree = 1, corrected = F)
sim.f4[17, ] <- mcCoverage(f4, m = 2, n = 50, degree = 1, corrected = F)
sim.f4[18, ] <- mcCoverage(f4, m = 3, n = 50, degree = 1, corrected = F)

sim.f4[19, ] <- mcCoverage(f4, m = 1, n = 10, degree = 2, corrected = T)
sim.f4[20, ] <- mcCoverage(f4, m = 2, n = 10, degree = 2, corrected = T)
sim.f4[21, ] <- mcCoverage(f4, m = 3, n = 10, degree = 2, corrected = T)
sim.f4[22, ] <- mcCoverage(f4, m = 1, n = 30, degree = 2, corrected = T)
sim.f4[23, ] <- mcCoverage(f4, m = 2, n = 30, degree = 2, corrected = T)
sim.f4[24, ] <- mcCoverage(f4, m = 3, n = 30, degree = 2, corrected = T)
sim.f4[25, ] <- mcCoverage(f4, m = 1, n = 50, degree = 2, corrected = T)
sim.f4[26, ] <- mcCoverage(f4, m = 2, n = 50, degree = 2, corrected = T)
sim.f4[27, ] <- mcCoverage(f4, m = 3, n = 50, degree = 2, corrected = T)
sim.f4[28, ] <- mcCoverage(f4, m = 1, n = 10, degree = 2, corrected = F)
sim.f4[29, ] <- mcCoverage(f4, m = 2, n = 10, degree = 2, corrected = F)
sim.f4[30, ] <- mcCoverage(f4, m = 3, n = 10, degree = 2, corrected = F)
sim.f4[31, ] <- mcCoverage(f4, m = 1, n = 30, degree = 2, corrected = F)
sim.f4[32, ] <- mcCoverage(f4, m = 2, n = 30, degree = 2, corrected = F)
sim.f4[33, ] <- mcCoverage(f4, m = 3, n = 30, degree = 2, corrected = F)
sim.f4[34, ] <- mcCoverage(f4, m = 1, n = 50, degree = 2, corrected = F)
sim.f4[35, ] <- mcCoverage(f4, m = 2, n = 50, degree = 2, corrected = F)
sim.f4[36, ] <- mcCoverage(f4, m = 3, n = 50, degree = 2, corrected = F)

## Covert to data frame and add extra columns
sim.f4 <- as.data.frame(sim.f4)
names(sim.f4) <- c("coverage", "length")
sim.f4$m <- rep(1:3, times = 12)
sim.f4$n <- rep(rep(c(10, 30, 50), each = 3), times = 4)
sim.f4$degree <- rep(1:2, each = 18)
sim.f4$corrected <- rep(rep(c("yes", "no"), each = 9), times = 2)
xyplot(coverage ~ n | as.factor(degree)*as.factor(corrected), groups = m, 
       type = "b", auto.key = list(space = "right", title = expression(m)), 
       xlab = "Coverage probability", ylab = "Sample size", data = sim.f4, 
       pch = 19, abline = list(h = 0.95))
write.csv(sim.f4, "/home/w108bmg/Desktop/Dissertation-knitr/R code/Simulation data/sim.f4.csv", row.names = F)

################################################################################
##
## New Monte Carlo Study
##
################################################################################

## Load libraries --------------------------------------------------------------
library(lattice)
library(investr)
library(plyr)
library(parallel)
source("/home/w108bmg/Desktop/Dissertation-knitr/R code/lmmSmooth.R")

## Simulation functions --------------------------------------------------------

## True f(x)
f <- function(x) sin(pi*x - pi/2)/2 + 0.5

## Function to simulate scatterplot data
simData <- function(f, m = 3, n = 10, sd = 0.05, from = 0, to = 1) {
  x <- rep(seq(from = from, to = to, length = n), each = m)
  y <- f(x) + rnorm(m*n, sd = sd)
  data.frame(x = x, y = y)
}

## Function to simulate coverage probability for pspline
mcCoverage <- function(f, m, n, sd = 0.05, num.sim = 1000, x0 = 0.75, degree, ...) {
  dataframes <- rlply(num.sim, simData(f = f, m = m, n = n, sd = sd))
  fun <- function(df) { 
    fit <- with(df, lmmSmooth(x, y, degree = degree))
    Y0 <- f(x0) + rnorm(1, sd = sd)
    res <- tryCatch(invest(fit, y0 = Y0, ...), 
                    error = function(e) NULL)
    if (is.null(res)) {
      repeat {
        newdf <- simData(f = f, m = m, n = n, sd = sd)
        newfit <- with(newdf, lmmSmooth(x, y, degree = degree))
        newY0 <- f(x0) + rnorm(1, sd = sd)
        res <- tryCatch(invest(newfit, y0 = newY0, ...), 
                        error = function(e) NULL)
        if (!is.null(res)) break
      }
    }
    c(if (res["lower"] < x0 && res["upper"] > x0) 1 else 0,
      abs(res["upper"] - res["lower"]))
  }
  z <- matrix(unlist(mclapply(dataframes, fun, mc.cores = 4)), 
              nrow = num.sim, ncol = 2, byrow = TRUE)
  colnames(z) <- c("coverage", "length")
  apply(z, 2, mean)
}

## Function to simulate coverage probability for linear model
mcCoverage <- function(f, m, n, sd = 0.05, num.sim = 1000, x0 = 0.75, degree, ...) {
  dataframes <- rlply(num.sim, simData(f = f, m = m, n = n, sd = sd))
  fun <- function(df) { 
    fit <- lm(y~f(x), data = df)
    Y0 <- f(x0) + rnorm(1, sd = sd)
    res <- tryCatch(invest(fit, y0 = Y0, ...), error = function(e) NULL)
    if (is.null(res)) {
      repeat {
        newdf <- simData(f = f, m = m, n = n, sd = sd)
        newfit <- lm(y~f(x), data = newdf)
        newY0 <- f(x0) + rnorm(1, sd = sd)
        res <- tryCatch(invest(newfit, y0 = newY0, ...), error = function(e) NULL)
        if (!is.null(res)) break
      }
    }
    c(if (res["lower"] < x0 && res["upper"] > x0) 1 else 0,
      abs(res["upper"] - res["lower"]))
  }
  z <- matrix(unlist(mclapply(dataframes, fun, mc.cores = 4)), 
              nrow = num.sim, ncol = 2, byrow = TRUE)
  colnames(z) <- c("coverage", "length")
  apply(z, 2, mean)
}

## Simulation (calibration) ----------------------------------------------------
sim1 <- matrix(0, nrow = 54, ncol = 2)
sim1[1, ]  <- mcCoverage(f, m = 1, n = 10, degree = 1, corrected = T, mean.response = FALSE)
sim1[2, ]  <- mcCoverage(f, m = 2, n = 10, degree = 1, corrected = T, mean.response = FALSE)
sim1[3, ]  <- mcCoverage(f, m = 3, n = 10, degree = 1, corrected = T, mean.response = FALSE)
sim1[4, ]  <- mcCoverage(f, m = 1, n = 30, degree = 1, corrected = T, mean.response = FALSE)
sim1[5, ]  <- mcCoverage(f, m = 2, n = 30, degree = 1, corrected = T, mean.response = FALSE)
sim1[6, ]  <- mcCoverage(f, m = 3, n = 30, degree = 1, corrected = T, mean.response = FALSE)
sim1[7, ]  <- mcCoverage(f, m = 1, n = 50, degree = 1, corrected = T, mean.response = FALSE)
sim1[8, ]  <- mcCoverage(f, m = 2, n = 50, degree = 1, corrected = T, mean.response = FALSE)
sim1[9, ]  <- mcCoverage(f, m = 3, n = 50, degree = 1, corrected = T, mean.response = FALSE)
sim1[10, ] <- mcCoverage(f, m = 1, n = 10, degree = 1, corrected = F, mean.response = FALSE)
sim1[11, ] <- mcCoverage(f, m = 2, n = 10, degree = 1, corrected = F, mean.response = FALSE)
sim1[12, ] <- mcCoverage(f, m = 3, n = 10, degree = 1, corrected = F, mean.response = FALSE)
sim1[13, ] <- mcCoverage(f, m = 1, n = 30, degree = 1, corrected = F, mean.response = FALSE)
sim1[14, ] <- mcCoverage(f, m = 2, n = 30, degree = 1, corrected = F, mean.response = FALSE)
sim1[15, ] <- mcCoverage(f, m = 3, n = 30, degree = 1, corrected = F, mean.response = FALSE)
sim1[16, ] <- mcCoverage(f, m = 1, n = 50, degree = 1, corrected = F, mean.response = FALSE)
sim1[17, ] <- mcCoverage(f, m = 2, n = 50, degree = 1, corrected = F, mean.response = FALSE)
sim1[18, ] <- mcCoverage(f, m = 3, n = 50, degree = 1, corrected = F, mean.response = FALSE)

sim1[19, ] <- mcCoverage(f, m = 1, n = 10, degree = 2, corrected = T, mean.response = FALSE)
sim1[20, ] <- mcCoverage(f, m = 2, n = 10, degree = 2, corrected = T, mean.response = FALSE)
sim1[21, ] <- mcCoverage(f, m = 3, n = 10, degree = 2, corrected = T, mean.response = FALSE)
sim1[22, ] <- mcCoverage(f, m = 1, n = 30, degree = 2, corrected = T, mean.response = FALSE)
sim1[23, ] <- mcCoverage(f, m = 2, n = 30, degree = 2, corrected = T, mean.response = FALSE)
sim1[24, ] <- mcCoverage(f, m = 3, n = 30, degree = 2, corrected = T, mean.response = FALSE)
sim1[25, ] <- mcCoverage(f, m = 1, n = 50, degree = 2, corrected = T, mean.response = FALSE)
sim1[26, ] <- mcCoverage(f, m = 2, n = 50, degree = 2, corrected = T, mean.response = FALSE)
sim1[27, ] <- mcCoverage(f, m = 3, n = 50, degree = 2, corrected = T, mean.response = FALSE)
sim1[28, ] <- mcCoverage(f, m = 1, n = 10, degree = 2, corrected = F, mean.response = FALSE)
sim1[29, ] <- mcCoverage(f, m = 2, n = 10, degree = 2, corrected = F, mean.response = FALSE)
sim1[30, ] <- mcCoverage(f, m = 3, n = 10, degree = 2, corrected = F, mean.response = FALSE)
sim1[31, ] <- mcCoverage(f, m = 1, n = 30, degree = 2, corrected = F, mean.response = FALSE)
sim1[32, ] <- mcCoverage(f, m = 2, n = 30, degree = 2, corrected = F, mean.response = FALSE)
sim1[33, ] <- mcCoverage(f, m = 3, n = 30, degree = 2, corrected = F, mean.response = FALSE)
sim1[34, ] <- mcCoverage(f, m = 1, n = 50, degree = 2, corrected = F, mean.response = FALSE)
sim1[35, ] <- mcCoverage(f, m = 2, n = 50, degree = 2, corrected = F, mean.response = FALSE)
sim1[36, ] <- mcCoverage(f, m = 3, n = 50, degree = 2, corrected = F, mean.response = FALSE)

sim1[37, ] <- mcCoverage(f, m = 1, n = 10, degree = 3, corrected = T, mean.response = FALSE)
sim1[38, ] <- mcCoverage(f, m = 2, n = 10, degree = 3, corrected = T, mean.response = FALSE)
sim1[39, ] <- mcCoverage(f, m = 3, n = 10, degree = 3, corrected = T, mean.response = FALSE)
sim1[40, ] <- mcCoverage(f, m = 1, n = 30, degree = 3, corrected = T, mean.response = FALSE)
sim1[41, ] <- mcCoverage(f, m = 2, n = 30, degree = 3, corrected = T, mean.response = FALSE)
sim1[42, ] <- mcCoverage(f, m = 3, n = 30, degree = 3, corrected = T, mean.response = FALSE)
sim1[43, ] <- mcCoverage(f, m = 1, n = 50, degree = 3, corrected = T, mean.response = FALSE)
sim1[44, ] <- mcCoverage(f, m = 2, n = 50, degree = 3, corrected = T, mean.response = FALSE)
sim1[45, ] <- mcCoverage(f, m = 3, n = 50, degree = 3, corrected = T, mean.response = FALSE)
sim1[46, ] <- mcCoverage(f, m = 1, n = 10, degree = 3, corrected = F, mean.response = FALSE)
sim1[47, ] <- mcCoverage(f, m = 2, n = 10, degree = 3, corrected = F, mean.response = FALSE)
sim1[48, ] <- mcCoverage(f, m = 3, n = 10, degree = 3, corrected = F, mean.response = FALSE)
sim1[49, ] <- mcCoverage(f, m = 1, n = 30, degree = 3, corrected = F, mean.response = FALSE)
sim1[50, ] <- mcCoverage(f, m = 2, n = 30, degree = 3, corrected = F, mean.response = FALSE)
sim1[51, ] <- mcCoverage(f, m = 3, n = 30, degree = 3, corrected = F, mean.response = FALSE)
sim1[52, ] <- mcCoverage(f, m = 1, n = 50, degree = 3, corrected = F, mean.response = FALSE)
sim1[53, ] <- mcCoverage(f, m = 2, n = 50, degree = 3, corrected = F, mean.response = FALSE)
sim1[54, ] <- mcCoverage(f, m = 3, n = 50, degree = 3, corrected = F, mean.response = FALSE)

## Covert to data frame and add extra columns
sim1 <- as.data.frame(sim1)
names(sim1) <- c("coverage", "length")
sim1$coverage <- round(sim1$coverage, 2)
sim1$m <- rep(1:3, times = 18)
sim1$n <- rep(rep(c(10, 30, 50), each = 3), times = 6)
sim1$degree <- rep(1:3, each = 18)
sim1$corrected <- rep(rep(c("yes", "no"), each = 9), times = 3)
xyplot(coverage ~ n | as.factor(degree)*as.factor(corrected), groups = m, 
       type = "b", auto.key = list(space = "right", title = expression(m)), 
       xlab = "Sample size", ylab = "Coverage probability", data = sim1, 
       main = "Calibration", pch = 19, abline = list(h = 0.95))
save(sim1, file = "/home/w108bmg/Desktop/Dissertation-knitr/Data/sim1.RData")

## Simulation (regulation) -----------------------------------------------------
sim2 <- matrix(0, nrow = 54, ncol = 2)
sim2[1, ]  <- mcCoverage(f, m = 1, n = 10, degree = 1, corrected = T, mean.response = TRUE)
sim2[2, ]  <- mcCoverage(f, m = 2, n = 10, degree = 1, corrected = T, mean.response = TRUE)
sim2[3, ]  <- mcCoverage(f, m = 3, n = 10, degree = 1, corrected = T, mean.response = TRUE)
sim2[4, ]  <- mcCoverage(f, m = 1, n = 30, degree = 1, corrected = T, mean.response = TRUE)
sim2[5, ]  <- mcCoverage(f, m = 2, n = 30, degree = 1, corrected = T, mean.response = TRUE)
sim2[6, ]  <- mcCoverage(f, m = 3, n = 30, degree = 1, corrected = T, mean.response = TRUE)
sim2[7, ]  <- mcCoverage(f, m = 1, n = 50, degree = 1, corrected = T, mean.response = TRUE)
sim2[8, ]  <- mcCoverage(f, m = 2, n = 50, degree = 1, corrected = T, mean.response = TRUE)
sim2[9, ]  <- mcCoverage(f, m = 3, n = 50, degree = 1, corrected = T, mean.response = TRUE)
sim2[10, ] <- mcCoverage(f, m = 1, n = 10, degree = 1, corrected = F, mean.response = TRUE)
sim2[11, ] <- mcCoverage(f, m = 2, n = 10, degree = 1, corrected = F, mean.response = TRUE)
sim2[12, ] <- mcCoverage(f, m = 3, n = 10, degree = 1, corrected = F, mean.response = TRUE)
sim2[13, ] <- mcCoverage(f, m = 1, n = 30, degree = 1, corrected = F, mean.response = TRUE)
sim2[14, ] <- mcCoverage(f, m = 2, n = 30, degree = 1, corrected = F, mean.response = TRUE)
sim2[15, ] <- mcCoverage(f, m = 3, n = 30, degree = 1, corrected = F, mean.response = TRUE)
sim2[16, ] <- mcCoverage(f, m = 1, n = 50, degree = 1, corrected = F, mean.response = TRUE)
sim2[17, ] <- mcCoverage(f, m = 2, n = 50, degree = 1, corrected = F, mean.response = TRUE)
sim2[18, ] <- mcCoverage(f, m = 3, n = 50, degree = 1, corrected = F, mean.response = TRUE)

sim2[19, ] <- mcCoverage(f, m = 1, n = 10, degree = 2, corrected = T, mean.response = TRUE)
sim2[20, ] <- mcCoverage(f, m = 2, n = 10, degree = 2, corrected = T, mean.response = TRUE)
sim2[21, ] <- mcCoverage(f, m = 3, n = 10, degree = 2, corrected = T, mean.response = TRUE)
sim2[22, ] <- mcCoverage(f, m = 1, n = 30, degree = 2, corrected = T, mean.response = TRUE)
sim2[23, ] <- mcCoverage(f, m = 2, n = 30, degree = 2, corrected = T, mean.response = TRUE)
sim2[24, ] <- mcCoverage(f, m = 3, n = 30, degree = 2, corrected = T, mean.response = TRUE)
sim2[25, ] <- mcCoverage(f, m = 1, n = 50, degree = 2, corrected = T, mean.response = TRUE)
sim2[26, ] <- mcCoverage(f, m = 2, n = 50, degree = 2, corrected = T, mean.response = TRUE)
sim2[27, ] <- mcCoverage(f, m = 3, n = 50, degree = 2, corrected = T, mean.response = TRUE)
sim2[28, ] <- mcCoverage(f, m = 1, n = 10, degree = 2, corrected = F, mean.response = TRUE)
sim2[29, ] <- mcCoverage(f, m = 2, n = 10, degree = 2, corrected = F, mean.response = TRUE)
sim2[30, ] <- mcCoverage(f, m = 3, n = 10, degree = 2, corrected = F, mean.response = TRUE)
sim2[31, ] <- mcCoverage(f, m = 1, n = 30, degree = 2, corrected = F, mean.response = TRUE)
sim2[32, ] <- mcCoverage(f, m = 2, n = 30, degree = 2, corrected = F, mean.response = TRUE)
sim2[33, ] <- mcCoverage(f, m = 3, n = 30, degree = 2, corrected = F, mean.response = TRUE)
sim2[34, ] <- mcCoverage(f, m = 1, n = 50, degree = 2, corrected = F, mean.response = TRUE)
sim2[35, ] <- mcCoverage(f, m = 2, n = 50, degree = 2, corrected = F, mean.response = TRUE)
sim2[36, ] <- mcCoverage(f, m = 3, n = 50, degree = 2, corrected = F, mean.response = TRUE)

sim2[37, ] <- mcCoverage(f, m = 1, n = 10, degree = 3, corrected = T, mean.response = TRUE)
sim2[38, ] <- mcCoverage(f, m = 2, n = 10, degree = 3, corrected = T, mean.response = TRUE)
sim2[39, ] <- mcCoverage(f, m = 3, n = 10, degree = 3, corrected = T, mean.response = TRUE)
sim2[40, ] <- mcCoverage(f, m = 1, n = 30, degree = 3, corrected = T, mean.response = TRUE)
sim2[41, ] <- mcCoverage(f, m = 2, n = 30, degree = 3, corrected = T, mean.response = TRUE)
sim2[42, ] <- mcCoverage(f, m = 3, n = 30, degree = 3, corrected = T, mean.response = TRUE)
sim2[43, ] <- mcCoverage(f, m = 1, n = 50, degree = 3, corrected = T, mean.response = TRUE)
sim2[44, ] <- mcCoverage(f, m = 2, n = 50, degree = 3, corrected = T, mean.response = TRUE)
sim2[45, ] <- mcCoverage(f, m = 3, n = 50, degree = 3, corrected = T, mean.response = TRUE)
sim2[46, ] <- mcCoverage(f, m = 1, n = 10, degree = 3, corrected = F, mean.response = TRUE)
sim2[47, ] <- mcCoverage(f, m = 2, n = 10, degree = 3, corrected = F, mean.response = TRUE)
sim2[48, ] <- mcCoverage(f, m = 3, n = 10, degree = 3, corrected = F, mean.response = TRUE)
sim2[49, ] <- mcCoverage(f, m = 1, n = 30, degree = 3, corrected = F, mean.response = TRUE)
sim2[50, ] <- mcCoverage(f, m = 2, n = 30, degree = 3, corrected = F, mean.response = TRUE)
sim2[51, ] <- mcCoverage(f, m = 3, n = 30, degree = 3, corrected = F, mean.response = TRUE)
sim2[52, ] <- mcCoverage(f, m = 1, n = 50, degree = 3, corrected = F, mean.response = TRUE)
sim2[53, ] <- mcCoverage(f, m = 2, n = 50, degree = 3, corrected = F, mean.response = TRUE)
sim2[54, ] <- mcCoverage(f, m = 3, n = 50, degree = 3, corrected = F, mean.response = TRUE)

## Covert to data frame and add extra columns
sim2 <- as.data.frame(sim2)
names(sim2) <- c("coverage", "length")
sim2$coverage <- round(sim1$coverage, 2)
sim2$m <- rep(1:3, times = 18)
sim2$n <- rep(rep(c(10, 30, 50), each = 3), times = 6)
sim2$degree <- rep(1:3, each = 18)
sim2$corrected <- rep(rep(c("yes", "no"), each = 9), times = 3)
xyplot(coverage ~ n | as.factor(degree)*as.factor(corrected), groups = m, 
       type = "b", auto.key = list(space = "right", title = expression(m)), 
       xlab = "Sample size", ylab = "Coverage probability", data = sim2, 
       main = "Regulation", pch = 19, abline = list(h = 0.95))
save(sim2, file = "/home/w108bmg/Desktop/Dissertation-knitr/Data/sim2.RData")
