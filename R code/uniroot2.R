uniroot2 <- function (f, interval, ..., lower = min(interval), 
                      upper = max(interval), f.lower = f(lower, ...), 
                      f.upper = f(upper, ...), 
                      tol = .Machine$double.eps^0.25, maxiter = 1000,
                      extend = c("none", "left", "right", "both"), frac = 0.1,
                      maxk = 1000) 
{
  extend <- match.arg(extend)
  if (!missing(interval) && length(interval) != 2L) {
    stop("'interval' must be a vector of length 2")
  }
  if (!is.numeric(lower) || !is.numeric(upper) || lower >= upper) {
    stop("lower < upper  is not fulfilled")
  }
  if (is.na(f.lower)) {
    stop("f.lower = f(lower) is NA")
  }
  if (is.na(f.upper)) {
    stop("f.upper = f(upper) is NA")
  }
  if (f.lower * f.upper > 0) {
    if (extend != "none") {
      k <- 1
      if (extend == "left") { # extend lower bound
        repeat {
          lower <- lower - frac
          f.lower = f(lower, ...)
          k <- k + 1
          if (f.lower * f.upper <= 0 || k >= maxk) break
        }
      } else if (extend == "right") { # extend upper bound
        repeat {
          upper <- upper + frac
          f.upper = f(upper, ...)
          k <- k + 1
          if (f.lower * f.upper <= 0 || k >= maxk) break
        }
      } else { # extend lower bound and upper bound
        repeat {
          lower <- lower - frac
          upper <- upper + frac
          f.lower = f(lower, ...)
          f.upper = f(upper, ...)
          k <- k + 1
          if (f.lower * f.upper <= 0 || k >= maxk) break
        }
      }
    } else {
      stop("f() values at end points not of opposite sign")
    }
  }
  uniroot(f, lower = lower, upper = upper, tol = tol, maxiter = maxiter)
}

