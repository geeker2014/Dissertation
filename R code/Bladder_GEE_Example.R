# Load libraries
library(gee)
library(geepack)
library(nlme)

# Set up data
subject <- rep(1:23, times = 8)
volume <- rep(c(10, 25, 50, 75, 100, 125, 150, 175), each = 23)
HD <- c( 13.2,  11.1,  10.3,    NA,   4.8,   7.7,    NA,   5.9,   1.9,   6.5,  19.8, 14.6,
    NA,    NA,   9.7,  17.2,  10.6,  19.3,   8.5,   6.9,   8.1,  14.8,  13.7,  27.4,
  27.5,  15.0,  10.0,  18.6,  12.6,  24.0,  28.4,  12.5,  16.7,  29.6,  27.1,  14.0,
  18.7,  20.3,  35.8,  23.6,  37.4,  31.3,  23.7,  22.0,  34.3,  28.5,  41.6,  58.1,
  34.2,  28.8,  29.9,  31.4,  46.9,  44.4,  26.8,  30.6,  51.7,  49.8,  19.1,  35.8,
  38.9,  41.4,  49.9,  58.6,  54.8,  44.0,  39.1,  58.5,  41.5,  60.1,  78.8,  49.4,
  46.4,  39.4,  45.3,  50.4,  70.7,  54.4,  41.8,  72.2,  67.5,  39.2,  49.6,  65.1,
  69.7,  67.7,  73.7,  78.3,  65.7,  44.7,  72.1,  59.8,  73.9,  91.5,  71.3,  54.8,
    NA,  48.0,  67.8,  89.4,  63.1,  49.6,  81.9,  79.1,  48.7,  65.6,  65.1,  81.9,
  87.7,  79.4,  93.0,  80.3,  68.9,  90.9,  77.5,  85.5,  98.3,  81.3,  69.4,    NA,
  66.6,  81.0, 105.8,  83.5,  60.8,  95.1,  95.1,  67.0,  85.3,  86.9,  96.6,  89.3,
 102.6,    NA,  93.6,  93.3, 105.0,  92.9,  95.6, 111.4,  94.0,  73.9,    NA,    NA,
  91.2, 113.5, 114.5,  80.1, 115.4, 109.8,  72.7,  90.4,  98.6, 115.0, 108.0, 110.9,
    NA,  99.2, 102.4, 117.5,  99.4, 107.4, 121.0, 104.3,    NA,    NA,    NA,  99.8,
 127.3, 124.0,  87.1,    NA,    NA,    NA,    NA, 107.2, 117.0, 114.8, 122.4,    NA,
 112.2, 104.7, 124.2, 113.0)
bladder <- data.frame(subject = subject, HD = HD, volume = volume)
bladder2 <- na.omit(bladder)

# Model formula
mf <- formula(HD ~ volume)

# Fit LM
bladder.lm <- lm(mf, data = bladder2)
(covb.lm <- vcov(bladder.lm))

# Use GEE
bladder.gee <- gee(mf, data = bladder2, id = subject, corstr = "exchangeable", tol = 1e-5, maxiter = 100)
(covb.gee <- bladder.gee$robust.variance)

# Use GLS
#bladder.gls <- gls(mf, data = bladder2, correlation = corSymm(form = ~ 1 | subject))
bladder.gee <- geeglm(HD ~ volume, data = bladder2, id = subject, corstr = "exchangeable")
bladder.gee <- geeglm(HD ~ volume, data = bladder2, id = subject, corstr = "ar1")

# Inverse estimation for mu = 60 (x0 = 77.60075)
x0 <- (60 - coef(bladder.lm)[1]) / coef(bladder.lm)[2]
x0.nse <- sqrt((1/(coef(bladder.lm)[2]^2)) * (covb.lm[1, 1]  + 2*x0*covb.lm[1, 2]  + x0^2*covb.lm[2, 2]))
x0.rse <- sqrt((1/(coef(bladder.lm)[2]^2)) * (covb.gee[1, 1] + 2*x0*covb.gee[1, 2] + x0^2*covb.gee[2, 2]))
x0.nint <- x0 + qt(c(0.025, 0.975), 166-2)*x0.nse
x0.rint <- x0 + qt(c(0.025, 0.975), 166-2)*x0.rse

# Plot data and results
xx <- seq(from = 0, to = 200, length = 1000)
out.conf <- as.data.frame(predict(bladder.lm, list(volume = xx), interval = "confidence"))
out.pred <- as.data.frame(predict(bladder.lm, list(volume = xx), interval = "prediction"))
with(bladder2, plot(volume, HD,
  pch = 1,
  cex = 1,
  col = "black",
  las = 1,
  main = "GEE vs SLR",
  xlab = "True volume (ml)",
  ylab = "Height (ml) times depth (ml)",
  panel.first = {
    polygon(c(xx, rev(xx)), c(out.pred$lwr, rev(out.pred$upr)), col = grey(0.9), border = grey(0.9))
    polygon(c(xx, rev(xx)), c(out.conf$lwr, rev(out.conf$upr)), col = grey(0.7), border = grey(0.7))
    lines(xx, out.conf$fit, lwd = 2, col = "black")
    abline(h = 60, lty = 1, col = "black")
  }))
segments(x0, -10, x0, 60)
rug(x0.nint, col = "blue2", lwd = 2)
rug(x0.rint, col = "red2", lwd = 2)






timeorder <- rep(1:5, 6)
tvar <- timeorder + rnorm(length(timeorder))
idvar <- rep(1:6, each=5)
uuu <- rep(rnorm(6), each=5)
yvar <- 1 + 2*tvar + uuu + rnorm(length(tvar))
simdat <- data.frame(idvar, timeorder, tvar, yvar)
head(simdat,12)
mod1 <- geeglm(yvar~tvar, id=idvar, data=simdat, corstr="ar1")
mod1
