## Load packages ---------------------------------------------------------------
library(investr)
library(rjags)
library(coda)
library(RColorBrewer)
set1 <- brewer.pal(9, "Set1")
dark2 <- brewer.pal(8, "Dark2")

## Import ELISA data frame -----------------------------------------------------
elisa <- read.csv("/home/w108bmg/Desktop/Dissertation/Data/elisa.csv", header = T)

## Source pspline code
source("/home/w108bmg/Desktop/Dissertation/R code/pspline.R")

## Fit models ------------------------------------------------------------------
mod1 <- nls(resp ~ b1 + (b2 - b1)/(1 + exp(b4*(log(conc) - b3))), 
            start = list(b1 = 25, b2 = 1, b3 = 1, b4 = 1), data = elisa)
mod2 <- with(elisa, pspline(conc, resp, degree = 2))
mod3 <- with(elisa, pspline(conc, resp, degree = 3))

## Calibration for the point y0 = 20
cal1 <- invest(mod1, y0 = 20)
cal2 <- invest(mod2, y0 = 20)
cal3 <- invest(mod3, y0 = 20)

## Plot fits with 95% prediction bands

# should fix plotFit so that ylim can be specified!

par(mfrow = c(1, 2))
plotFit(mod1, interval = "prediction", lty.pred = 1, extend.range = F,
        main = "Four-parameter logistic function")
abline(h = 20, col = "purple", lwd = 2)
arrows(cal1$lower, 20, cal1$lower, -0.5, col = "purple", lwd = 2)
arrows(cal1$upper, 20, cal1$upper, -0.5, col = "purple", lwd = 2)
plot(mod2, interval = "prediction", main = "Quadratic P-spline", ylim = c(0.5, 29.5))
abline(h = 20, col = "purple", lwd = 2)
arrows(cal2[2], 20, cal2[2], -0.5, col = "purple", lwd = 2)
arrows(cal2[3], 20, cal2[3], -0.5, col = "purple", lwd = 2)

## Bayesian FPLM ---------------------------------------------------------------

## Model file for JAGS 
model1 <- "
model {

## Likelihood
for (i in 1:length(y)) {
  y[i] ~ dnorm(beta[1]+(beta[2]-beta[1])/(1+exp(beta[4]*(log(x[i])-beta[3]))), tau.error) 
}

## Priors for fixed effects
for (i in 1:4) {
  beta[i] ~ dnorm(0, 1.0E-6)
}

## For calibration
y0 ~ dnorm(beta[1]+(beta[2]-beta[1])/(1+exp(beta[4]*(log(x0)-beta[3]))), tau.error) 
x0 ~ dunif(0, 50) 

## Priors for precision parameters
tau.error ~ dgamma(1.0E-6, 1.0E-6)

## Other calculations
sigma.error <- 1/sqrt(tau.error)

## Predict new observations
for (i in 1:250) {
  epsilon.star[i] ~ dnorm(0, tau.error)
  y.star[i] <- mu.star[i] + epsilon.star[i]
  mu.star[i] <- beta[1]+(beta[2]-beta[1])/(1+exp(beta[4]*(log(newx[i])-beta[3])))
}

}
"
model1.file <- "/home/w108bmg/Desktop/Dissertation/R code/JAGS models/model1.txt"
writeLines(model1, con = model1.file)

## Inputs for JAGS
data.list <- list(y = elisa$resp, x = elisa$conc, y0 = 20, 
                  newx = seq(from = 0, to = 50, length = 250))

initsFun <- function() { 
  list(beta = coef(elisa.nls), tau.error = 0.1, x0 = 9.1837) 
}

adapt.steps <- 1000    # number of steps to "tune" the samplers
burnin.steps <- 10000  # number of steps to "burn-in" the samplers
n.chains <- 1          # number of chains to run
n.saved.steps <- 50000 # total number of steps in chains to save
thin.steps <- 1        # number of steps to "thin" (1 = keep every step)
n.per.chain <- ceiling((n.saved.steps*thin.steps) / n.chains) # steps per chain

# Run JAGS model
sim <- jags.model(model.file, data = data.list, inits = initsFun, 
                  n.chains = n.chains, n.adapt = adapt.steps)
update(sim, n.iter = 10000) # burn-in

## Generate posterior samples from density of x0
x0.coda <- coda.samples(sim, variable.names = "x0", n.iter = 50000)
x0.post <- as.numeric(as.matrix(x0.coda))
dens.x0 <- density(x0.post)
dens.x0$x[which.max(dens.x0$y)]
quantile(x0.post, c(0.025, 0.975))

## Generate posterior samples from density of mu.star
mu.coda <- coda.samples(sim, variable.names = "mu.star", n.iter = 50000)
mu.post <- as.matrix(mu.coda)

## Generate posterior samples from density of y.star
pred.coda <- coda.samples(sim, variable.names = "y.star", n.iter = 50000)
pred.post <- as.matrix(pred.coda)

## Generate posterior samples from density of beta (fixed effects)
beta.coda <- coda.samples(sim, variable.names = "beta", n.iter = 50000)
beta.post <- as.matrix(beta.coda)

## Write posterior samples to external file
save(x0.coda, x0.post, mu.coda, mu.post, pred.coda, pred.post, 
     beta.coda, beta.post, 
     file = "/home/w108bmg/Desktop/Dissertation-knitr/Data/elisa-jags-nls.RData")

## Bayesian Quadratic P-spline -------------------------------------------------

## Model file for JAGS 
model2 <- "
model {
  
  ## Likelihood
  for (i in 1:length(y)) {
    y[i] ~ dnorm(inprod(beta[], X[i,]) + inprod(alpha[], Z[i,]), tau.error) 
  }
  
  ## Priors for fixed effects
  for (j in 1:degree+1) {
    beta[j] ~ dnorm(0, 1.0E-6)
  }
  
  ## Priors for random effects
  for (k in 1:num.knots) {
    alpha[k] ~ dnorm(0, tau.alpha)
  }
  
  ## Calibration
  for (j in 1:degree+1) {
    X0[j] <- pow(x0, j-1)
  }
  for (k in 1:num.knots) {
    u0[k] <- (x0-knots[k])*step(x0-knots[k])
    Z0[k] <- pow(u0[k], degree)
  }
  y0 ~ dnorm(inprod(beta[], X0[]) + inprod(alpha[], Z0[]), tau.error) 
  x0 ~ dunif(0, 50) 
  #x0 ~ dnorm(0, 1.0E-6) 
  #x0 ~ dgamma(1.0E-3, 1.0E-3)
   
  ## Priors for precision parameters
  tau.error ~ dgamma(1.0E-3, 1.0E-3)
  tau.alpha ~ dgamma(1.0E-3, 1.0E-3)
  
  ## Other calculations
  sigma.error <- 1/sqrt(tau.error)
  sigma.alpha <- 1/sqrt(tau.alpha)
  lambda <- pow(sigma.alpha, 2)/pow(sigma.error, 2)
  
  ## Predict new observations
  for (i in 1:250) {
    epsilon.star[i] ~ dnorm(0, tau.error)
    y.star[i] <- mu.star[i] + epsilon.star[i]
    mu.star[i] <- inprod(beta[], newX[i,]) + inprod(alpha[], newZ[i,])
  }
  
}
"
model2.file <- "/home/w108bmg/Desktop/Dissertation/R code/JAGS models/model2.txt"
writeLines(model2, model2.file)

## Calculate design matrices
knots <- mod2$knots
X <- cbind(1, poly(elisa$conc, degree = 2, raw = TRUE))
Z <- outer(elisa$conc, knots, "-")
Z <- (Z * (Z > 0)) ^ 2
newx <- seq(from = 0, to = 50, length = 250)
newX <- cbind(1, poly(newx, degree = 2, raw = TRUE))
newZ <- outer(newx, knots, "-")
newZ <- (newZ * (newZ > 0)) ^ 2

## Inputs for JAGS
data.list <- list(
  y = elisa$resp,
  y0 = 20,
  X = X,
  Z = Z,
  newX = newX,
  newZ = newZ,
  knots = c(0.1683333, 0.6166667, 1.8500000, 5.5000000, 16.7500000),
  num.knots = 5, 
  degree = 2 # degree of polynomial
)
# initsFun <- function() { 
#   list(beta = rep(0, 3), alpha = rep(0, 5), tau.error = 0.1, tau.alpha = 0.01,
#        x0 = 10) 
# }
# initsFun <- function() { 
#   list(beta = rnorm(3, mean = mod2$beta.hat), 
#        alpha = rnorm(5, sd = mod2$var.components["u"]), 
#        tau.error = 0.1, tau.alpha = 0.01,
#        x0 = runif(1, min = 0, max = 20)) 
# }

adapt.steps <- 10000   # number of steps to "tune" the samplers
burnin.steps <- 20000  # number of steps to "burn-in" the samplers
n.chains <- 3          # number of chains to run
n.saved.steps <- 30000 # total number of steps in chains to save
thin.steps <- 5        # number of steps to "thin" (1 = keep every step)
n.per.chain <- ceiling((n.saved.steps*thin.steps) / n.chains) # steps per chain

# Run JAGS model
sim <- jags.model(model2.file, data = data.list, #inits = initsFun, 
                  n.chains = n.chains, n.adapt = adapt.steps)
update(sim, n.iter = burnin.steps) # burn-in

## Generate posterior samples from density of x0
x0.psp2.coda <- coda.samples(sim, variable.names = "x0", n.iter = n.per.chain, 
                            thin = thin.steps)
x0.psp2.post <- as.numeric(as.matrix(x0.psp2.coda))
x0.psp2.dens <- density(x0.psp2.post)
x0.psp2.mode <- x0.psp2.post[which.max(x0.psp2.dens$y)]

par(mfrow = c(1, 2))
traceplot(x0.psp2.coda, col = set1)
abline()
cumuplot(x0.psp2.coda, auto.layout = FALSE)
autocorr.plot(x0.psp2.coda, auto.layout = FALSE)
heidel.diag(x0.psp2.coda)

## Generate posterior samples from density of mu.star
mu.psp2.coda <- coda.samples(sim, variable.names = "mu.star", n.iter = n.saved.steps)
mu.psp2.post <- as.matrix(mu.psp2.coda)

## Generate posterior samples from density of y.star
pred.psp2.coda <- coda.samples(sim, variable.names = "y.star", n.iter = n.saved.steps)
pred.psp2.post <- as.matrix(pred.psp2.coda)

## Generate posterior samples from density of beta (fixed effects)
beta.psp2.coda <- coda.samples(sim, variable.names = "beta", n.iter = n.saved.steps)
beta.psp2.post <- as.matrix(beta.psp2.coda)

## Generate posterior samples from density of alpha (random effects)
alpha.psp2.coda <- coda.samples(sim, variable.names = "alpha", n.iter = n.saved.steps)
alpha.psp2.post <- as.matrix(alpha.psp2.coda)

## Generate posterior samples from density of lambda (smoothing parameter)
lambda.psp2.coda <- coda.samples(sim, variable.names = "lambda", n.iter = n.saved.steps)
lambda.psp2.post <- as.numeric(as.matrix(lambda.psp2.coda))

## Write posterior samples to external file
save(x0.coda, x0.post, mu.coda, mu.post, pred.coda, pred.post, 
     beta.coda, beta.post, alpha.coda, alpha.post, lambda.coda, lambda.post, 
     file = "/home/w108bmg/Desktop/elisa-jags/elisa-jags.RData")

## Plot smooths with 95% prediction intervals
par(las = 1, cex.axis = 0.9, mfrow = c(1, 2))

mu <- apply(mu.post, 2, mean)
pred <- apply(pred.post, 2, quantile, prob = c(0.025, 0.975))
plot(elisa, type = "n", main = "Bayesian smooth", ylim = c(0, 30),
     xlab = "Concentration", ylab = "Response")
polygon(c(newx, rev(newx)), c(pred[1,], rev(pred[2,])),
        col = "skyblue", border = "skyblue")
points(elisa)
lines(newx, mu, lwd = 2)
dens.x0 <- density(x0.post)
polygon(dens.x0$x, (10/max(dens.x0$y))*dens.x0$y*dens.x0$y, col = "plum2", 
        border = "black", lwd = 2)
rug(HDIofMCMC(x0.post), lwd = 3, lend = 2)
abline(v = dens.x0$x[which.max(dens.x0$y)], lty = 2)

pspline.lwr <- predict(mod2, newx, interval = "prediction")$lwr
pspline.upr <- predict(mod2, newx, interval = "prediction")$upr
plot(elisa, type = "n", main = "Mixed model smooth",
     xlab = "Concentration", ylab = "Response")
abline(v = mod2$knots, col = "darkgrey", lty = 3)
polygon(c(newx, rev(newx)), c(pspline.lwr, rev(pspline.upr)),
        col = "skyblue", border = "skyblue")
points(elisa)
lines(newx, predict(mod2, newx), lwd = 2)
abline(h = 20, col = "purple", lwd = 2)
segments(cal2[2], 20, cal2[2], -5, col = "purple", lwd = 2)
segments(cal2[3], 20, cal2[3], -5, col = "purple", lwd = 2)

## Bayesian Cubic P-spline -----------------------------------------------------

## Model file for JAGS 
model3 <- "
model {

## Likelihood
for (i in 1:length(y)) {
y[i] ~ dnorm(inprod(beta[], X[i,]) + inprod(alpha[], Z[i,]), tau.error) 
}

## Priors for fixed effects
for (j in 1:degree+1) {
beta[j] ~ dnorm(0, 1.0E-6)
}

## Priors for random effects
for (k in 1:num.knots) {
alpha[k] ~ dnorm(0, tau.alpha)
}

## Calibration
for (j in 1:degree+1) {
X0[j] <- pow(x0, j-1)
}
for (k in 1:num.knots) {
u0[k] <- (x0-knots[k])*step(x0-knots[k])
Z0[k] <- pow(u0[k], degree)
}
y0 ~ dnorm(inprod(beta[], X0[]) + inprod(alpha[], Z0[]), tau.error) 
x0 ~ dunif(0, 50) 
#x0 ~ dnorm(0, 1.0E-6) 
#x0 ~ dgamma(1.0E-3, 1.0E-3)

## Priors for precision parameters
tau.error ~ dgamma(1.0E-3, 1.0E-3)
tau.alpha ~ dgamma(1.0E-3, 1.0E-3)

## Other calculations
sigma.error <- 1/sqrt(tau.error)
sigma.alpha <- 1/sqrt(tau.alpha)
lambda <- pow(sigma.alpha, 2)/pow(sigma.error, 2)

## Predict new observations
for (i in 1:250) {
epsilon.star[i] ~ dnorm(0, tau.error)
y.star[i] <- mu.star[i] + epsilon.star[i]
mu.star[i] <- inprod(beta[], newX[i,]) + inprod(alpha[], newZ[i,])
}

}
"
model3.file <- "/home/w108bmg/Desktop/Dissertation/R code/JAGS models/model3.txt"
writeLines(model3, model3.file)

## Calculate design matrices
knots <- mod3$knots
X <- cbind(1, poly(elisa$conc, degree = 3, raw = TRUE))
Z <- outer(elisa$conc, knots, "-")
Z <- (Z * (Z > 0)) ^ 3
newx <- seq(from = 0, to = 50, length = 250)
newX <- cbind(1, poly(newx, degree = 3, raw = TRUE))
newZ <- outer(newx, knots, "-")
newZ <- (newZ * (newZ > 0)) ^ 3

## Inputs for JAGS
data.list <- list(
  y = elisa$resp,
  y0 = 20,
  X = X,
  Z = Z,
  newX = newX,
  newZ = newZ,
  knots = mod3$knots,
  num.knots = length(mod3$knots), 
  degree = 3 # degree of polynomial
)

adapt.steps <- 10000   # number of steps to "tune" the samplers
burnin.steps <- 20000  # number of steps to "burn-in" the samplers
n.chains <- 3          # number of chains to run
n.saved.steps <- 30000 # total number of steps in chains to save
thin.steps <- 5        # number of steps to "thin" (1 = keep every step)
n.per.chain <- ceiling((n.saved.steps*thin.steps) / n.chains) # steps per chain

# Run JAGS model
sim <- jags.model(model3.file, data = data.list, #inits = initsFun, 
                  n.chains = n.chains, n.adapt = adapt.steps)
update(sim, n.iter = burnin.steps) # burn-in

## Generate posterior samples from density of x0
x0.psp3.coda <- coda.samples(sim, variable.names = "x0", n.iter = n.per.chain, 
                             thin = thin.steps)
x0.psp3.post <- as.numeric(as.matrix(x0.psp3.coda))
x0.psp3.dens <- density(x0.psp3.post)
x0.psp3.mode <- x0.psp3.post[which.max(x0.psp3.dens$y)]

par(mfrow = c(1, 2))
traceplot(x0.psp3.coda, col = set1)
abline()
cumuplot(x0.psp3.coda, auto.layout = FALSE)
autocorr.plot(x0.psp3.coda, auto.layout = FALSE)
heidel.diag(x0.psp3.coda)

## Bootstrap cubic smoothing spline --------------------------------------------
elisa.smooth <- with(elisa, smooth.spline(conc, resp))
elisa.small <- with(elisa, smooth.spline(conc, resp, 
                                         spar = elisa.smooth$spar/2))
elisa.big <- with(elisa, smooth.spline(conc, resp, spar = elisa.smooth$spar*2))

plot(elisa)
lines(elisa.smooth, col = "black", lty = 1, lwd = 1)
lines(elisa.small, col = "red", lty = 2, lwd = 1)
lines(elisa.big, col = "blue", lty = 4, lwd = 1)

res <- (elisa$resp - elisa.small$y)#/sqrt(1 - elisa.small$lev)
elisa.mle <- data.frame(bigfit = elisa.big$y, res = res)
xpoints <- c(10, 20, 25, 30, 35, 45)
elisaFun <- function(data, x) {
  y.smooth <- smooth.spline(data$conc, data$resp)
  predict(y.smooth, x)$y
  
#   Y0 <- y0 + 
  
}
elisaGen <- function(data, mle) {
  d <- data
  i <- sample(c(1:nrow(data)), replace = TRUE)
  d$resp <- mle$bigfit + mle$res[i]
  d
}
elisa.boot <- boot(elisa, elisaFun, R = 9999, sim = "parametric",
                   ran.gen = elisaGen, mle = elisa.mle, x = xpoints)
mu.big <- predict(elisa.big, xpoints)$y
mu <- predict(elisa.smooth, xpoints)$y
ylims <- apply(elisa.boot$t, 2, quantile, c(0.025, 0.975))
lower <- mu - (ylims[2,]-mu.big)
upper <- mu - (ylims[1,]-mu.big)
points(c(xpoints, xpoints), c(lower, upper), pch = 15, col = "red")