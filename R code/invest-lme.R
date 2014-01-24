## Warning: only for the simple linear random intercept model 
invest.lme <- function(object, y0, interval = c("inversion", "Wald"), 
                       level = 0.95, mean.respone = FALSE, lower, upper, 
                       tol = .Machine$double.eps^0.25, maxiter = 1000, ...) {
  
  ## Data frame and variables
  d <- eval(object$call$data)
  yvar <- all.vars(formula(object)[[2]])
  xvar <- intersect(all.vars(formula(object)[[3]]), colnames(d))
  if (missing(lower)) {
    lower <- min(d[, xvar])
  }
  if (missing(upper)) {
    upper <- max(d[, xvar])
  }
  
  ## Variance of new observations
  var.y0 <- summary(object)$sigma^2 + getVarCov(object)[1]
  
#   yname <- all.vars(formula(simdata.lme))[1]
#   xname <- all.vars(formula(simdata.lme))[2]
#   yvar <- d[[yname]]
#   xvar <- d[[xname]]
  
  ## Compute estimate and confidence bounds
  interval <- match.arg(interval)
  if (interval == "inversion") { # Inversion interval
    
    ## Alternate prediction function
    predFun <- function(x) {
      z <- list(x)
      names(z) <- xvar
      fit <- predict(object, newdata = z, level = 0)
      se.fit <- sqrt(diag(cbind(1, unlist(z)) %*% object$varFix %*% 
                         t(cbind(1, unlist(z)))))
      list(fit = fit, se.fit = se.fit)
    }
    
    ## Calculate ML estimate
    invFun.est <- function(x) {
      z <- list(x)
      names(z) <- xvar
      predict(object, newdata = z, level = 0) - y0
    }
    x0.est <- uniroot(invFun.est, interval = c(lower, upper), tol = tol, 
                      maxiter = maxiter)$root
    
    ## Calculate confidence bounds
    invFun.bounds <- function(x) { 
      z <- list(x)
      names(z) <- xvar
      pred <- predFun(x)
      (y0 - pred$fit)^2/(var.y0 + pred$se.fit^2) - qnorm(0.975)^2
    }
    lwr <- uniroot(invFun.bounds, interval = c(lower, x0.est), tol = tol, 
                   maxiter = maxiter)$root
    upr <- uniroot(invFun.bounds, interval = c(x0.est, upper), tol = tol, 
                   maxiter = maxiter)$root
    
    ## Print results
    c(estimate = x0.est, lower = lwr, upper = upr)
    
  } else { # Wald interval
    
    ## Calculate covariance matrix
    covmat <- diag(3)
    covmat[1:2, 1:2] <- vcov(object)
    covmat[3,3] <- var.y0
    b <- as.numeric(fixef(object))
    
    ## Delta method
    object.dm <- car:::deltaMethod(c(b0 = b[1], b1 = b[2], Y0 = y0), 
                                   "(Y0 - b0)/b1", vcov. = covmat)
    
    ## Print results
    c(estimate = object.dm$Estimate, 
      lower = object.dm$Estimate - qnorm(0.975)*object.dm$SE,
      upper = object.dm$Estimate + qnorm(0.975)*object.dm$SE,
      se = object.dm$SE)
  }
  
}