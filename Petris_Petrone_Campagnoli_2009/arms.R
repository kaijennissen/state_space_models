# Petris, G. & Petrone, S. & Campagnoli, P. (2009).
# Dynamic Linear Models with R.
library(dlm)

# Adaptive Rejection Metropolis Sampling ----------------------------------
support <- function(x) {
  return(as.numeric(-1 < x[2] && x[2] < 1 &&
    -2 < x[1] &&
    (x[1] < 1 || crossprod(x - c(1, 0)) < 1)))
}
Min.log <- log(.Machine$double.xmin) + 10
logf <- function(x) {
  if (x[1] < 0) {
    return(log(1 / 4))
  } else
  if (crossprod(x - c(1, 0)) < 1) {
    return(log(1 / pi))
  }
  return(Min.log)
}
x <- as.matrix(expand.grid(seq(-2.2, 2.2, length = 40), seq(-1.1, 1.1, length = 40)))
y <- sapply(1:nrow(x), function(i) support(x[i, ]))
plot(x, type = "n", asp = 1)
points(x[y == 1, ], pch = 1, cex = 1, col = "green")
z <- arms(c(0, 0), logf, support, 5000)
points(z, pch = 20, cex = 0.5, col = "blue")
polygon(c(-2, 0, 0, -2), c(-1, -1, 1, 1))
curve(-sqrt(1 - (x - 1)^2), 0, 2, add = TRUE)
curve(sqrt(1 - (x - 1)^2), 0, 2, add = TRUE)
sum(z[, 1] < 0) # sampled points in the square
sum(apply(t(z) - c(1, 0), 2, crossprod) < 1) # sampled points in the circle
