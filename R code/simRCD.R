#' Function to simulate data from a (balanced) random coefficient model
#'
#' @param n The number of observations per subject.
#' @param m The number of subjects.
#' @param fixed A vector of length 2 representing the fixed effects.
#' @param vars A vector of length 3 representing the variance components.
#' @return A data frame with m*n rows and 3 columns: x, y, and subject.
simRCD <- function(n = 10, m = 5, fixed = c(0, 1), vars = c(0.5, 0, 0.25)) {
  subject <- rep(1:m, each = n)
  x <- rep(seq(from = 0, to = 1, length = n), times = m) 
  B0 <- rnorm(m, mean = 0, sd = sqrt(vars[1]))
  B1 <- rnorm(m, mean = 0, sd = sqrt(vars[2]))
  y <- rnorm(m * n, fixed[1] + B0[subject] + fixed[2]*x + B1[subject]*x, 
             sd = sqrt(vars[3]))
  data.frame(x = x, y = y, subject = subject)
}

library(lattice)
d1 <- simRCD(vars = c(0.2, 0, 0.001))
d2 <- simRCD(vars = c(0, 0.1, 0.001))
d3 <- simRCD(vars = c(0.2, 0.1, 0.001))
p1 <- xyplot(y ~ x, groups = subject, data = d1, type = c("p", "r"))
p2 <- xyplot(y ~ x, groups = subject, data = d2, type =  c("p", "r"))
p3 <- xyplot(y ~ x, groups = subject, data = d3, type =  c("p", "r"))
print(p1, pos = c(0, 0, 1/3, 1), more = TRUE) 
print(p2, pos = c(1/3, 0, 2/3, 1), more = TRUE) 
print(p3, pos = c(2/3, 0, 1, 1), more = FALSE) 