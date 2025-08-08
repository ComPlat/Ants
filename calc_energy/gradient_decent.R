gradient_descent <- function(f, x0,
                             lr = 0.01, tol = 1e-6,
                             max_iter = 1000,
                             h = 1e-5, ...) {
  backtrack <- function(f, x, grad, fx, lr, beta = 0.5) {
    while (TRUE) {
      x_new <- x - lr * grad
      if (f(x_new) < fx) break
      lr <- lr * beta
      if (lr < 1e-12) break
    }
    list(x_new = x_new, lr = lr)
  }
  x <- x0
  for (i in 1:max_iter) {
    fx <- f(x, ...)
    grad <- numeric(length(x))
    for (j in seq_along(x)) {
      xh <- x
      xh[j] <- xh[j] + h
      grad[j] <- (f(xh, ...) - fx) / h
    }

    temp <- backtrack(f, x, grad, fx, lr)
    x_new <- temp$x_new
    lr <- temp$lr
    # x_new <- x - lr * grad

    if (sqrt(sum((x_new - x)^2)) < tol) {
      message("Converged in ", i, " iterations")
      return(x_new)
    }

    x <- x_new
  }

  stop("Did not converge")
}
