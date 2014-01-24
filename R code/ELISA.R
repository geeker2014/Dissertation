x <- rep(c( 0.0000,  0.0750,  0.1025,  0.1350,  0.1850,  0.2500,
            0.4000,  0.5500,  0.7500,  1.0000,  1.3750,  1.8500,
            2.5000,  3.2500,  4.5000,  6.0000,  8.2500, 11.2500,
           15.0000, 20.2500, 27.5000, 37.0000, 50.0000), each = 4)
y <- c( 1.700,  1.660,  1.950,  2.070, 
        1.910,  2.270,  2.110,  2.390,
        2.220,  2.250,  3.260,  2.920,
        2.800,  2.940,  2.380,  2.700,
        2.780,  2.640,  2.710,  2.850,
        3.540,  2.860,  3.150,  3.320,
       
        3.910,  3.830,  4.880,  4.210,
        4.540,  4.470,  4.790,  5.680,
        6.060,  5.070,  5.000,  5.980,
        5.840,  5.790,  6.100,  7.810, 
        7.310,  7.080,  7.060,  6.870,
        9.880, 10.120,  9.220,  9.960,
       
       11.040, 10.460, 10.880, 11.650,
       13.510, 15.470, 14.210, 13.920,
       16.070, 14.670, 14.780, 15.210,
       17.340, 16.850, 16.740, 16.870,
       18.980, 19.850, 18.750, 18.510,
       21.666, 21.218, 19.790, 22.669,
        
       23.206, 22.239, 22.436, 22.597,
       23.922, 24.871, 23.815, 24.871,
       25.478, 25.874, 24.907, 24.871,
       24.441, 25.874, 25.748, 27.270,
       29.580, 26.698, 26.536, 27.181)
elisa <- data.frame(x = x, y = y)
fit.nls <- nls(y ~ b1 + (b2 - b1)/(1 + exp(b4*(log(x) - b3))), 
               start = list(b1 = 25, b2 = 1, b3 = 1, b4 = 1), data = elisa)

library(investr)
source("/home/w108bmg/Desktop/Dissertation-knitr/R code/lmmSmooth.R")
invest(fit.nls, y0 = 20)

fit.lmm <- lmmSmooth(elisa$x, elisa$y, degree = 2)
plot(fit.lmm, interval = "prediction")
invest(fit.nls, y0 = 20)
invest(fit.lmm, y0 = 20, corrected = T)

## Bootstrap cubic smoothing spline (css)
library(boot)
fit.css <- smooth.spline(elisa$x, elisa$y)
x0.fun <- function(obj, y0 = 20) {
  uniroot(function(x) {predict(obj, x = x)$y - y0}, interval = c(0.01, 250), 
          tol = 1e-10, maxiter = 1000)$root
}
x0.orig <- x0.fun(fit.css)
res <- elisa$y - fit.css$y
res <- res - mean(res)
boot.data <- data.frame(elisa, res = res, fit = fit.css$y)
boot.fun <- function(data, i) {
  d <- data
  boot.fit <- smooth.spline(d$x, d$fit+d$res[i])
  Y0 <- 20 + sample(d$res, size = 1)
  x0 <- if (all(i == 1:nrow(d))) {
          x0.orig 
        } else {
          tryCatch(x0.fun(boot.fit, y0 = Y0), error = function(e) NULL)
        }
  while (is.null(x0)) {
    newfit <- smooth.spline(d$x, d$fit+d$res[i])
    Y0 <- 20 + sample(d$res, size = 1)
    x0 <- tryCatch(x0.fun(newfit, y0 = Y0), 
                   error = function(e) NULL)
  }
  x0
}
sim <- boot(boot.data, boot.fun, R = 99)
boot.ci(sim)