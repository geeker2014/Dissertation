## Load packages and data
library(SemiPar)
data(fossil)
attach(fossil)

## Using pspline.R -------------------------------------------------------------
fit.spm <- spm(strontium.ratio ~ f(age, basis = "trunc.poly", degree = 2))
knotloc <- fit.spm$info$pen$knots[[1]]
fit.lme <- pspline(age, strontium.ratio, degree = 2, knots = knotloc)

## Plot
#plot(fit.spm)
plot(fit.lme, interval = "confidence", adjusted = F)
xseq <- seq(from = min(age), to = max(age), length = 500)
pred <- predict(fit.lme, newdata = xseq, interval = "confidence")
lines(xseq, pred$fit, col = "red", lty = 2, lwd = 2)
lines(xseq, pred$lwr, col = "red", lty = 2, lwd = 2)
lines(xseq, pred$upr, col = "red", lty = 2, lwd = 2)
pred2 <- predict(fit.lme, newdata = xseq, interval = "confidence", adjusted = FALSE)
lines(xseq, pred2$lwr, lty = 2)
lines(xseq, pred2$upr, lty = 2)

## Usin lmmSmooth.R ------------------------------------------------------------
fit.spm <- spm(strontium.ratio ~ f(age, basis = "trunc.poly", degree = 2))
knotloc <- fit.spm$info$pen$knots[[1]]
fit.lme <- lmmSmooth(age, strontium.ratio, degree = 2, knots = knotloc)

## Plot
plot(fit.spm)
xseq <- seq(from = min(age), to = max(age), length = 500)
pred <- predict(fit.lme, newdata = xseq, interval = "prediction")
lines(xseq, pred$fit, col = "red", lty = 2, lwd = 2)
lines(xseq, pred$lwr, col = "red", lty = 2, lwd = 2)
lines(xseq, pred$upr, col = "red", lty = 2, lwd = 2)
pred2 <- predict(fit.lme, newdata = xseq, interval = "prediction", corrected = F)
lines(xseq, pred2$lwr, lty = 2)
lines(xseq, pred2$upr, lty = 2)

