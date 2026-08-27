.positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 || x != as.integer(x)) {
    stop(name, " must be a positive integer.", call. = FALSE)
  }
  as.integer(x)
}

.nonnegative_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 0 || x != as.integer(x)) {
    stop(name, " must be a non-negative integer.", call. = FALSE)
  }
  as.integer(x)
}

.check_matrix_dim <- function(x, nr, nc, name) {
  if (!is.numeric(x) || !identical(dim(x), c(nr, nc)) || any(!is.finite(x))) {
    stop(name, " must be a finite numeric ", nr, " by ", nc, " matrix.",
         call. = FALSE)
  }
  invisible(x)
}

.safe_solve <- function(x) {
  tryCatch(solve(x), error = function(e) MASS::ginv(x))
}

.stabilize_covariance <- function(x, epsilon = 1e-8) {
  x <- (x + t(x)) / 2
  decomposition <- eigen(x, symmetric = TRUE)
  if (min(decomposition$values) < epsilon) {
    x <- decomposition$vectors %*%
      (pmax(decomposition$values, epsilon) * t(decomposition$vectors))
    x <- (x + t(x)) / 2
  }
  x
}
