uniroot2 <- function (f, interval, ..., lower = min(interval), 
                      upper = max(interval), f.lower = f(lower, ...), 
                      f.upper = f(upper, ...), 
                      tol = .Machine$double.eps^0.25, maxiter = 1000,
                      extend = c("none", "left", "right", "both")) 
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
      if (extend == "left") { # extend lower bound
        repeat {
          lower <- lower - frac
          f.lower = f(lower, ...)
          if (f.lower * f.upper <= 0) break
        }
      } else if (extend == "right") { # extend upper bound
        repeat {
          upper <- upper + frac
          f.upper = f(upper, ...)
          if (f.lower * f.upper <= 0) break
        }
      } else { # extend lower bound and upper bound
        repeat {
          lower <- lower - frac
          upper <- upper + frac
          f.lower = f(lower, ...)
          f.upper = f(upper, ...)
          if (f.lower * f.upper <= 0) break
        }
      }
    } else {
      stop("f() values at end points not of opposite sign")
    }
  }
           
  val <- .External2(C_zeroin2, function(arg) f(arg, ...), lower, 
                    upper, f.lower, f.upper, tol, as.integer(maxiter))
  iter <- as.integer(val[2L])
  if (iter < 0) {
    warning(sprintf(ngettext(maxiter, "_NOT_ converged in %d iteration", 
                             "_NOT_ converged in %d iterations"), maxiter), domain = NA)
    iter <- maxiter
  }
  list(root = val[1L], f.root = f(val[1L], ...), iter = iter, 
       estim.prec = val[3L])
}