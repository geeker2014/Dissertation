#' Fit a Penalized Regression Spline
#' 
#' Fits a penalized regression spline to the supplied data using a linear mixed 
#' model representation.
#' 
#' @param x A vector giving the values of the explanatory variable.
#' @param y A vector giving the values of the response variable.
#' @param degree Degree of piecewise polynomial, default is 1 for piecewise
#'               linear splines.
#' @param knots The internal breakpoints that define the spline.
#' @param num.knots The number of knots.
#' @param basis This option is currently ignored.
#' @param random This option is currently ignored.
#' @param method Method used for estimating parameters in the linear mixed
#'               effects model, default is \code{"REML"} for restricted
#'               maximum likelihood.
#' @param ... Additional optional arguments to be passed on to the function
#'            \code{lme}.
#' @return An object of class \code{pspline} containing the following 
#'         components:
#'   \itemize{
#'     \item{x}{...}
#'     \item{fit}{...}
#'     \item{y}{...}
#'     \item{knots}{...}
#'     \item{degree}{...}
#'     \item{beta.hat}{...}
#'     \item{u.hat}{...}
#'     \item{spar}{...}
#'     \item{var.rat}{...}
#'     \item{C.mat}{...}
#'     \item{D.mat}{...}
#'     \item{S.mat}{...}
#'     \item{var.components}{...}
#'     \item{df.res}{Residual degrees of freedom.}
#'     \item{fit.lme}{...}
#'     \item{call}{...}
#'   }
#' @references 
#' Ruppert, D., Wand, M. & Carroll, R. (2003). Semiparametric regression. 
#' Cambridge New York: Cambridge University Press.
#' @rdname pspline
#' @aliases print.pspline
#' @export
pspline <- function(x, ...) {
  UseMethod("pspline")
}

#' @rdname pspline
#' @export
#' @method pspline data.frame
pspline.data.frame <- function(x, ...) {
  res <- pspline.default(x = x[[1]], y = x[[2]], ...)
  cl <- match.call()
  cl[[1L]] <- as.name("pspline")
  res$call <- cl
  res
}

#' @rdname pspline
#' @export
#' @method pspline matrix
pspline.matrix <- function(x, ...) {
  res <- pspline.default(x = x[, 1], y = x[, 2], ...)
  cl <- match.call()
  cl[[1L]] <- as.name("pspline")
  res$call <- cl
  res
}

#' @rdname pspline
#' @export
#' @method pspline default
pspline.default <- function(x, y, degree = 1, knots, num.knots, 
                   basis = c("trunc.poly", "radial"), 
                   random, method = c("REML", "ML"), ...) {
  
  ## TODO:
  ##   (1) Fix unadjusted intervals!
  ##   (2) Omit load packages section before adding to package.
  
  ## Load packages -------------------------------------------------------------
  require(splines)
  require(nlme)
  
  ## Catch errors --------------------------------------------------------------
  if (!(degree %in% 1:3)) {
    stop("degree must be either 1, 2, or 3")
  }
  
  ## Spline basis matrix -------------------------------------------------------
  if (missing(knots)) {
    if (missing(num.knots)) {
      num.knots <- max(5, min(floor(length(unique(x))/4), 35))
    }
    knotseq <- seq(from = 0, to = 1, length = num.knots + 2)
    knots <- as.numeric(quantile(unique(x), knotseq)[-c(1, (num.knots + 2))])
  } else {
    num.knots <- length(knots)
  }
  basis <- match.arg(basis)
  if (basis == "trunc.poly") {
    Z <- outer(x, knots, "-")
    Z <- (Z * (Z > 0)) ^ degree
    colnames(Z) <- paste("k", 1:ncol(Z), sep = "")
  } else {
    stop(paste(basis, "basis functions are not yet supported"))
    Z0 <- (abs(outer(x, knots, "-"))) ^ (2*degree - 1)
    svd.omega <- svd((abs(outer(knots, knots, "-"))) ^ (2*degree - 1))
    sqrt.omega <- t(svd.omega$v %*% (t(svd.omega$u)*sqrt(svd.omega$d)))
    Z <- t(solve(sqrt.omega, t(Z0)))
  }
  
  ## Fit LMM -------------------------------------------------------------------
  d <- data.frame(y = y, x = x, Z, subject = factor(rep(1, length(x))))
  random.model <- as.formula(paste("~-1+", paste(colnames(Z), collapse = "+")))
  fixed.model <- as.formula(paste("y ~ poly(x, degree =", degree, 
                                  ", raw = TRUE)"))
  fit.lme <- lme(fixed.model, data = d, method = match.arg(method),
                 random = list(subject = pdIdent(random.model)), ...)                          
  beta.hat <- as.numeric(fit.lme$coef$fixed)
  u.hat <- as.numeric(unlist(fit.lme$coef$random))
  
  ## Extract variance components and calculate smoothing parameter -------------
  var.e.hat <- fit.lme$sigma^2
  var.u.hat <- as.numeric(var.e.hat* 
                            exp(2*unlist(fit.lme$modelStruct)))
  var.rat <- var.e.hat / var.u.hat
  spar <- var.rat ^ (1 / (2 * degree))
  
  ## Calculate fitted values, key matrices, and residual df --------------------
  fit <- as.numeric(fitted(fit.lme))
  X <- cbind(1, poly(x, degree = degree, raw = TRUE))
  C.mat <- cbind(X, Z) # C matrix
  D.mat <- diag(c(rep(0, degree + 1), rep(1, num.knots))) 
  A.mat <- qr.solve(crossprod(C.mat) + var.rat*D.mat, tol = 1e-10) #############
  S.mat <- C.mat %*% tcrossprod(A.mat, C.mat) 
  df.res <- nrow(X) - 2*sum(diag(S.mat)) + sum(diag(tcrossprod(S.mat)))
  
  ## Return list of results ----------------------------------------------------
  cl <- match.call()
  cl[[1L]] <- as.name("pspline")
  res <- list(x = x, fit = fit, y = y, knots = knots, degree = degree, 
              beta.hat = beta.hat, u.hat = u.hat, spar = spar, A.mat = A.mat,
              C.mat = C.mat, D.mat = D.mat, S.mat = S.mat, var.rat = var.rat,
              var.components = c(error = var.e.hat, u = var.u.hat), 
              df.res = df.res, fit.lme = fit.lme, call = cl)
  class(res) <- "pspline"
  res
}

## Experimental function to compute spline basis matrix
splineBM <- function(x, basis, knots, degree) {
  if (basis == "trunc.poly") {
    Z <- outer(x, knots, "-")
    Z <- (Z * (Z > 0)) ^ degree
    colnames(Z) <- paste("k", 1:ncol(Z), sep = "")
  } else {
    Z0 <- (abs(outer(x, knots, "-"))) ^ (2*degree - 1)
    svd.omega <- svd((abs(outer(knots, knots, "-"))) ^ (2*degree - 1))
    sqrt.omega <- t(svd.omega$v %*% (t(svd.omega$u)*sqrt(svd.omega$d)))
    Z <- t(solve(sqrt.omega, t(Z0)))
  }
  Z
}

## Pedict method for objects of class "pspline" --------------------------------
predict.pspline <- function(object, newdata, se.fit = FALSE, 
                            interval = c("none", "confidence", "prediction"),
                            level = 0.95) {
  
  ## Extract needed components from model
  degree <- object$degree
  beta.hat <- object$beta.hat
  u.hat <- object$u.hat
  
  ## Create design matrices and compute fitted values
  newx <- if (missing(newdata)) object$x else newdata
  Z <- outer(newx, object$knots, "-")
  Z <- (Z * (Z > 0)) ^ degree
  X <- cbind(1, poly(newx, degree = degree, raw = TRUE))
  f.hat <- as.numeric(X %*% beta.hat + Z %*% u.hat)
  
  ## Standard error of fitted value
  if (se.fit || interval != "none") {
    C.mat <- cbind(X, Z)
    A.mat <- object$A.mat
    S.mat <- C.mat %*% tcrossprod(A.mat, C.mat)
    
    ## Calculate standard error of fit
    sigma.f <- sqrt(object$var.components[1] * diag(S.mat%*%t(S.mat)))
  }
  
  ## Confidence/prediction intervals
  interval <- match.arg(interval)
  if (interval != "none") {
    df.res <- object$df.res
    se <- if (interval == "confidence") {
            sigma.f
          } else {
            sqrt(object$var.components[1] + sigma.f^2)
          }
    
    ## Store results in a list
    res <- list(fit = f.hat)
    res$lwr <- f.hat - qt((1 + level)/2, df = df.res) * se
    res$upr <- f.hat + qt((1 + level)/2, df = df.res) * se
    if (se.fit) res$se.fit <- sigma.f
  } else {
    res <- if (se.fit) list(fit = f.hat, se.fit = sigma.f) else f.hat
  }
  
  ## Return list of results
  res
  
}

## Returns the residual degrees-of-freedom extracted from a P-spline object
df.residual.pspline <- function(object) {
  object$df.res
}

## Plot method for objects of class "pspline" ----------------------------------
plot.pspline <- function(object, ..., n = 500, lwd.smooth = 2, lty.smooth = 1, 
                         col.smooth = "black", knots = TRUE, adjusted = TRUE,
                         col.knots = "grey", lwd.knots = 1, lty.knots = 3,
                         interval = c("none", "confidence", "prediction")) {
  
  ## Extract needed components from pspline object
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
  plot(x, y, panel.first = {
    if (knots) {
      abline(v = object$knots, col = col.knots, lwd = lwd.knots, lty = lty.knots)
    }
  }, ...) 
  lines(newx, fit, col = col.smooth, lwd = lwd.smooth, 
        lty = lty.smooth)
  interval <- match.arg(interval)
  if (interval != "none") {
    conf <- predict.pspline(object, newdata = newx, interval = interval, 
                            adjusted = adjusted)
    lines(newx, conf$lwr)
    lines(newx, conf$upr)
  }
  
}

## pspline method for invest function ------------------------------------------
invest.pspline <- function(object, y0, level = 0.95, mean.response = FALSE, 
                           tol = 1e-10, maxiter = 1000) {
  
  if (length(y0) > 1) stop("only a single value for y0 is allowed")
  interval <- if(mean.response) "confidence" else "prediction"
                
  ## Functions for uniroot
  invEstFun <- function(x) {
    predict.pspline(object, newdata = x) - y0
  }
  invLwrFun <- function(x) {
    predict.pspline(object, newdata = x, interval = interval, 
                    level = level)$lwr - y0
  }
  invUprFun <- function(x) {
    predict.pspline(object, newdata = x, interval = interval,
                    level = level)$upr - y0
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

## Simulate method for 'pspline' objects ---------------------------------------
simulate.pspline <- function(object, nsim = 1, seed = NULL, ...) {
  
  ## Extract fitted LMM
  object <- object$fit.lme
  
  ## Random seed
  if(!is.null(seed)) {
    set.seed(seed)
  } 
  if(!exists(".Random.seed", envir = .GlobalEnv)) {
    runif(1) # initialize the RNG if necessary
  }
  
  ## Groups
  groups <- object$groups[[1]]
  unique.groups <- unique(groups)
  individuals <- as.matrix(unique.groups)
  if (is.numeric(individuals)) {
    individuals <- unique.groups[individuals]
  }
  
  ## Code from the function simulateY in the nlmeU package
  v.list <- nlme::getVarCov(object, individuals, type = "marginal")
  fit <- fitted(object, level = 0)
  ch.v.list <- lapply(v.list, chol)
  nx <- nsim * length(fit)
  noise.mtx <- matrix(rnorm(nx), nrow = length(fit), ncol = nsim)
  dim.n   <-   sapply(ch.v.list, ncol)  # Number of records per subject
  cdim.n1 <- cumsum(dim.n)
  cdim.n0 <- cumsum(dim.n) - dim.n + 1
  cdim.n  <- cbind(cdim.n0, cdim.n1)
  t.list <- vector(length(dim.n), mode = "list")
  t.list[] <- 1:length(dim.n)
  res <- lapply(t.list, FUN = function(el) {
    ch.v <- ch.v.list[[el]]
    ix <- cdim.n[el, ]
    i0 <- ix[1]
    i1 <- ix[2]
    noise.x <- noise.mtx[i0:i1, ]
    t(ch.v) %*% noise.x
  })
  res.all <- NULL 
  for (i in 1:length(dim.n)) {
    res.all <- rbind(res.all, res[[i]])
  }
  res.all + fit
  
}


