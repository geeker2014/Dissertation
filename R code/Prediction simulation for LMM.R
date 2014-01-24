## Load packages ---------------------------------------------------------------
library(lattice)
library(nlme)
library(rjags)
library(car)

## Function to simulate data ---------------------------------------------------
simRCD <- function(n = 10, m = 20, fixed = c(0, 1), vars = c(0, 0, 0.1)) {
  subject <- rep(1:m, each = n)
  x <- rep(seq(from = 0, to = 1, length = n), times = m) 
  B0 <- rnorm(m, mean = 0, sd = sqrt(vars[1]))
  B1 <- rnorm(m, mean = 0, sd = sqrt(vars[2]))
  y <- rnorm(m*n, mean = fixed[1]+B0[subject] + (fixed[2]+B1[subject])*x, 
             sd = sqrt(vars[3]))
  data.frame(x = x, y = y, subject = subject)
}

set.seed(1234)
simdata <- simRCD(m = 15, n = 30, fixed = c(0, 1), vars = c(0.01, 0, 0.001))
xyplot(y ~ x, data = simdata, groups = subject, type = "l",
       scales = list(tck = c(1, 0)))

## Calibration parameters ------------------------------------------------------
x0 <- 0.5
y0 <- 0.5

## Delta method ----------------------------------------------------------------
fit <- lme(y ~ x, random = ~ 1|subject, data = simdata)
covmat <- diag(3)
covmat[1:2, 1:2] <- vcov(fit)
covmat[3,3] <- summary(fit)$sigma^2 + getVarCov(fit)[1]

## Calibration
dm.cal <- deltaMethod(c(b0 = -0.03644237, b1 = 1.00586080, Y0 = 0.5), 
                     "(Y0 - b0)/b1", vcov. = covmat)
c(dm.cal$Estimate - qnorm(0.975)*dm.cal$SE,
  dm.cal$Estimate + qnorm(0.975)*dm.cal$SE)

## Regulation
dm.reg <- deltaMethod(fit, "(0.5-b0)/b1", parameterNames = c("b0", "b1"))
c(dm.reg$Estimate - qnorm(0.975)*dm.reg$SE,
  dm.reg$Estimate + qnorm(0.975)*dm.reg$SE)

## Inverting asymptotic prediction interval ------------------------------------
var.error <- summary(fit)$sigma^2
var.alpha <- getVarCov(fit)[1]
covb <- vcov(fit)
x.vec <- c(1, 0.5)
x.vec %*% (covb + var.alpha) %*% t(x.vec) + var.error


## Inverting distribution free prediction interval -----------------------------
res <- lm(y ~ x, data = simdata)$residuals
a <- quantile(res, 0.025) 
b <- quantile(res, 0.975)
newx <- seq(from = 0, to = 1, length = 500)
lower <- fixef(fit)[1]+fixef(fit)[2]*newx + a
upper <- fixef(fit)[1]+fixef(fit)[2]*newx + b
plot(y ~ x, data = simdata, col = "grey")
lines(newx, lower, col = "red2", lwd = 2)
lines(newx, upper, col = "red2", lwd = 2)
abline(h = 0.5, lty = 2, col = "red2")

df.lower <- ((y0-b) - fixef(fit)[1])/fixef(fit)[2]
df.upper <- ((y0-a) - fixef(fit)[1])/fixef(fit)[2]

## Distribution free inversion interval
abline(v = c(df.lower, df.upper), lty = 2, col = "red2")

## Wald interval
abline(v = c(dm.cal$Estimate - qnorm(0.975)*dm.cal$SE,
             dm.cal$Estimate + qnorm(0.975)*dm.cal$SE), lty = 2, col = "blue")

distFree <- function(data, y0, level = 0.95) {
  mod <- lme(y ~ x, random = ~ 1|subject, data = data)
  res <- lm(y ~ x, data = data)$residuals
  a <- quantile(res, 1-(1+level)/2) 
  b <- quantile(res, (1+level)/2)
  c(lower = as.numeric(((y0-b) - fixef(fit)[1])/fixef(fit)[2]), 
    upper = as.numeric(((y0-a) - fixef(fit)[1])/fixef(fit)[2]))
}

library(plyr)
library(parallel)
mcCoverage <- function(num.sim = 1000) {
  dataframes <- rlply(num.sim, simRCD(m = 50, n = 50, fixed = c(0, 1), 
                                      vars = c(0.01, 0, 0.001)))
  fun <- function(df) {
    Y0 <- rnorm(1, mean = x0, sd = sqrt(0.01 + 0.001))
    res <- distFree(df, y0 = Y0)
    if (res["lower"] < x0 && res["upper"] > x0) 1 else 0
  }
  mean(unlist(mclapply(dataframes, fun, mc.cores = 4)))
}
mcCoverage()

## Bayesian prediction interval ------------------------------------------------

## Model file for JAGS 
model.string <- "
model {

  ## Likelihood
  for (i in 1:N) {
      y[i] ~ dnorm(mu[i], tau.error) 
      mu[i] <- beta0 + beta1*x[i] + alpha[subject[i]]
  }

  ## Priors for fixed effects
  beta0 ~ dnorm(0, 1.0E-6)
  beta1 ~ dnorm(0, 1.0E-6)

  ## Priors for random intercepts
  for (j in 1:m) {
    alpha[j] ~ dnorm(0, tau.alpha)
  }

  ## Calibration
  y0 ~ dnorm(beta0 + beta1*x0, 1/total.var) 
  #x0 ~ dunif(0, 1) 
  x0 ~ dnorm(0, 1.0E-6)

  ## Prior for variance components
  tau.alpha ~ dgamma(1.0E-6, 1.0E-6)
  tau.error ~ dgamma(1.0E-6, 1.0E-6)

  ## Additional calculation
  sigma.error <- 1/sqrt(tau.error)
  sigma.alpha <- 1/sqrt(tau.alpha)
  total.var <- pow(sigma.error, 2) + pow(sigma.alpha, 2)

  ## Predicting new observations
#   for (i in 1:N) {
#     epsilon.star[i] ~ dnorm(0, tau.error)
#     y.star[i] <- mu[i] + epsilon.star[i]
#   }
  for (i in 1:100) {
    epsilon.star[i] ~ dnorm(0, 1/total.var)
    y.star[i] <- beta0 + beta1*newx[i] + epsilon.star[i]
  }

}
"
writeLines(model.string, con = "/home/w108bmg/Desktop/simdata-jags/model.txt")

# Specify data, as a list
data.list = list(N = 15*30, m = 15, y0 = 0.5, subject = simdata$subject, 
                 x = simdata$x, y = simdata$y, 
                 newx = seq(from = 0, to = 1, length = 100))

## Initialize chains
initsFun <- function() { 
  list(beta0 = 0, beta1 = 1, tau.error = 10, tau.alpha = 1000,
       alpha = rnorm(data.list$m, sd = sqrt(0.01))) 
}

adapt.steps <- 1000    # number of steps to "tune" the samplers
burnin.steps <- 10000  # number of steps to "burn-in" the samplers
n.chains <- 1          # number of chains to run
n.saved.steps <- 50000 # total number of steps in chains to save
thin.steps <- 1        # number of steps to "thin" (1 = keep every step)
n.per.chain <- ceiling((n.saved.steps*thin.steps) / n.chains) # steps per chain

# Run JAGS model
jags.mod <- jags.model("/home/w108bmg/Desktop/simdata-jags/model.txt", 
                       data = data.list, inits = initsFun, n.chains = n.chains, 
                       n.adapt = adapt.steps)
update(jags.mod, n.iter = 10000) # burn-in

x0.coda <- coda.samples(jags.mod, "x0", n.iter = 10000)
x0.samp <- as.matrix(x0.coda)

## Generate posterior samples from density of mu.star
mu.coda <- coda.samples(jags.mod, variable.names = "mu.star", n.iter = 10000)
mu.samp <- as.matrix(mu.coda)

## Generate posterior samples from density of y.star
pred.coda <- coda.samples(jags.mod, variable.names = "y.star", n.iter = 10000)
pred.samp <- as.matrix(pred.coda)

newx <- seq(from = 0, to = 1, length = 100)
mu <- apply(mu.post, 2, mean)
pred <- apply(pred.samp, 2, quantile, prob = c(0.025, 0.975))
plot(y ~ x, data = simdata, type = "n")
polygon(c(newx, rev(newx)), c(pred[1,], rev(pred[2,])),
        col = "skyblue", border = "skyblue")
points(y ~ x, data = simdata, cex = 0.5)
