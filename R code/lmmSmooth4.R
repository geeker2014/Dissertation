## Smooth based on mixed model -------------------------------------------------
lmmSmooth<- function(x, y, degree = 1, knots, num.knots, 
                     method = c("REML", "ML")) {
  ## Load required packages
  require(splines)
  require(nlme)
  
  ## Catch errors
  if (!(degree %in% 1:3)) {
    stop("degree must be either 1, 2, or 3")
  }
  
  ## Spline basis matrix
  if (missing(knots)) {
    if (missing(num.knots)) {
      num.knots <- max(5, min(floor(length(unique(x))/4), 35))
    }
    knotseq <- seq(from = 0, to = 1, length = num.knots + 2)
    knots <- as.numeric(quantile(unique(x), knotseq)[-c(1, (num.knots + 2))])
  } else {
    num.knots <- length(knots)
  }
  Z <- outer(x, knots, "-")
  Z <- (Z * (Z > 0)) ^ degree
  colnames(Z) <- paste("k", 1:ncol(Z), sep = "")
  
  ## Fit LME model
  d <- data.frame(y = y, x = x, Z, subject = factor(rep(1, length(x))))
  
  rmod <- as.formula(paste("~-1+", paste(colnames(Z), collapse = "+")))
  fmod <- as.formula(paste("y ~ poly(x, degree =", degree, ", raw = TRUE)"))
  fit.lme <- lme(fixed.model, data = d, method = match.arg(method),
                 random = list(subject = pdIdent(random.model)))  
  
  beta.hat <- as.numeric(fit.lme$coef$fixed)
  u.hat <- as.numeric(unlist(fit.lme$coef$random))
  
  ## Variance components and smoothing parameter
  var.e.hat <- fit.lme$sigma^2
  var.u.hat <- as.numeric(var.e.hat * 
                            exp(2 * unlist(fit.lme$modelStruct)))
  var.rat <- var.e.hat / var.u.hat
  spar <- var.rat ^ (1 / (2 * degree))
  
  ## Compute fitted values, key matrices, and residual df
  fit <- as.numeric(fitted(fit.lme))
  X <- cbind(1, poly(x, degree = degree, raw = TRUE))
  C.mat <- cbind(X, Z) # C matrix
  D.mat <- diag(c(rep(0, degree + 1), rep(1, num.knots))) # D matrix
  
  ## System can become computationally singular
  A.mat <- qr.solve(crossprod(C.mat) + var.rat*D.mat, tol = 1e-10)
  
  S.mat <- C.mat %*% tcrossprod(A.mat, C.mat) # smoother matrix
  df.res <- nrow(X) - 2*sum(diag(S.mat)) + sum(diag(tcrossprod(S.mat)))
  
  ## Return list of results
  res <- list(x = x, fit = fit, y = y, knots = knots, degree = degree, 
              beta.hat = beta.hat, u.hat = u.hat, spar = spar, 
              var.rat = var.rat, C.mat = C.mat, D.mat = D.mat, S.mat = S.mat,
              var.components = c(error = var.e.hat, u = var.u.hat), 
              df.res = df.res)
  class(res) <- "lmmSmooth"
  res
}

## Pedict method for objects of class "lmmSmooth" ------------------------------
predict.lmmSmooth <- function(object, newdata,  
                              interval = c("none", "confidence", "prediction"),
                              level = 0.95, corrected = TRUE) {
  
  ## Extract needed components from model
  degree <- object$degree
  beta.hat <- object$beta.hat
  u.hat <- object$u.hat
  
  ## Create design matrices and compute fitted values
  newx <- if (missing(newdata)) object$x else newdata
  Z <- outer(newx, object$knots, "-")
  Z <- (Z * (Z > 0)) ^ degree
  X <- cbind(1, poly(newx, degree = degree, raw = TRUE))
  f.hat <- X %*% beta.hat + Z %*% u.hat
  
  ## Confidence/prediction intervals
  interval <- match.arg(interval)
  if (interval != "none") {
    C.mat <- cbind(X, Z)
    A.mat <- qr.solve(crossprod(object$C.mat) + object$var.rat*object$D.mat, 
                      tol = 1e-10)
    S.mat <- C.mat %*% A.mat %*% t(C.mat)
    df.res <- object$df.res
    se <- if (interval == "confidence") {
      if (corrected) {
        sqrt(object$var.components)[1] * sqrt(diag(S.mat))
      } else {
        sqrt(object$var.components)[1] * sqrt(diag(S.mat %*% t(S.mat)))
      }
    } else {
      if (corrected) {
        sqrt(object$var.components)[1]*sqrt(1+diag(S.mat))
      } else {
        sqrt(object$var.components)[1]*sqrt(1+diag(S.mat%*%t(S.mat)))
      }
    }
    res <- list(fit = f.hat)
    res$lwr <- f.hat - qt((1 + level)/2, df = df.res) * se
    res$upr <- f.hat + qt((1 + level)/2, df = df.res) * se
  } else {
    res <- f.hat
  }
  
  ## Return list of results
  res
  
}

## Plot method for objects of class "lmmSmooth" --------------------------------
plot.lmmSmooth <- function(object, ..., n = 500, lwd.smooth = 2, lty.smooth = 1, 
                           col.smooth = "black", knots = TRUE, 
                           col.knots = "grey", lwd.knots = 1, lty.knots = 3,
                           interval = c("none", "confidence", "prediction")) {
  
  ## Extract needed components from lmmSmooth object
  x <- object$x
  y <- object$y
  newx <- seq(from = min(x), to = max(x), length = n)
  degree <- object$degree
  beta.hat <- object$beta.hat
  u.hat <- object$u.hat
  
  ## Create design matrices and compute fitted values
  Z <- outer(newx, object$knots, "-")
  Z <- (Z * (Z > 0)) ^ degree
  X <- cbind(1, poly(newx, degree = degree, raw = TRUE))
  fit <- X %*% beta.hat + Z %*% u.hat
  
  ## Plot data with smoothed line and confidence/prediction bands (if requested)
  plot(x, y, ...)
  lines(newx, fit, col = col.smooth, lwd = lwd.smooth, 
        lty = lty.smooth)
  if (knots) {
    abline(v = object$knots, col = col.knots, lwd = lwd.knots, lty = lty.knots)
  }
  interval <- match.arg(interval)
  if (interval != "none") {
    conf <- predict.lmmSmooth(object, newdata = newx, interval = interval)
    lines(newx, conf$lwr)
    lines(newx, conf$upr)
  }
  
}

## lmmSmooth method for invest function ----------------------------------------
invest.lmmSmooth <- function(object, y0, level = 0.95, mean.response = FALSE, 
                             corrected = TRUE, tol = 1e-10, maxiter = 1000) {
  
  ## TODO: 
  ##  * allow y0 to be a vector
  
  if (length(y0) > 1) stop("only a single value for y0 is allowed")
  interval <- if(mean.response) "confidence" else "prediction"
  
  ## Functions for uniroot
  invEstFun <- function(x) {
    predict.lmmSmooth(object, newdata = x) - y0
  }
  invLwrFun <- function(x) {
    predict.lmmSmooth(object, newdata = x, interval = interval,
                      level = level, corrected = corrected)$lwr - y0
  }
  invUprFun <- function(x) {
    predict.lmmSmooth(object, newdata = x, interval = interval,
                      level = level, corrected = corrected)$upr - y0
  } 
  
  ## Find roots
  xrange <- range(object$x)
  est <- uniroot(invEstFun, interval = xrange, tol = tol, 
                 maxiter = maxiter)$root
  lwr <- uniroot(invLwrFun, interval = xrange, tol = tol, 
                 maxiter = maxiter)$root
  upr <- uniroot(invUprFun, interval = xrange, tol = tol, 
                 maxiter = maxiter)$root
  
  ## Return vector of results
  c(estimate = est, lower = min(c(lwr, upr)), upper = max(c(lwr, upr)))
  
}
