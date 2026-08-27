# SpatMix

`SpatMix` implements spatial mixture models for complete or incomplete data,
including spatial Gaussian mixtures and mixtures of spatial factor analyzers.

The package has one main fitting function:

- `spatmix()` fits a spatial mixture model. Missing values
  automatically activate the built-in matrix partial EM updates.

## Install

```r
install.packages("SpatMix")
```

## Examples

The two simulations below use the same sigmoid covariance as `spatmix()`.

```r
library(SpatMix)

sigmoid <- function(d, beta, shift = 3) {
  z0 <- plogis(-shift)
  (plogis(beta * d - shift) - z0) /
    (plogis(2 * beta - shift) - z0)
}

spatial_cov <- function(coords, beta, alpha) {
  if (!is.list(coords)) coords <- list(coords)
  dims <- lengths(coords)
  decay <- lapply(seq_along(coords), function(j) {
    d <- as.matrix(dist(coords[[j]]))
    d <- 2 * d / max(d)
    matrices <- lapply(seq_along(coords), function(k) {
      if (j == k) sigmoid(d, beta) else matrix(1, dims[k], dims[k])
    })
    Reduce(kronecker, rev(matrices))
  })
  p <- prod(dims)
  J <- matrix(1, p, p) - diag(p)
  alpha[1] * J + alpha[2] * Reduce("+", decay) + alpha[3] * diag(p)
}
```

### Spatial-only data

```r
set.seed(1)
n <- 150
coords <- 1:10
beta <- c(3, 7)
alpha <- list(c(1, -0.30, 1.2), c(1, -0.45, 1.4))

Xi <- list()
for (g in 1:2) Xi[[g]] <- spatial_cov(coords, beta[g], alpha[[g]])

X <- rbind(
  MASS::mvrnorm(n, rep(0, 10), Xi[[1]]),
  MASS::mvrnorm(n, rep(2.5, 10), Xi[[2]])
)
truth <- rep(1:2, each = n)

fit <- spatmix(
  X, G = 2, coords = coords, spatial_decay = "sigmoid",
  sigmoid_ctrl = list(init = 4, lower = 0.1, upper = 12),
  max_iter = 50, tol = 0.01, verbose = FALSE
)

ord <- order(colMeans(fit$means[, 1, ]))
table(truth, fitted = match(fit$cluster, ord))
round(rbind(truth = beta, fitted = fit$sigmoid[1, ord]), 2)
round(cbind(
  truth.1 = alpha[[1]], fitted.1 = fit$alpha[, ord[1]],
  truth.2 = alpha[[2]], fitted.2 = fit$alpha[, ord[2]]
), 2)
```

### Spatial and non-spatial data

Here `Xi` is the spatial covariance and
`Omega = Lambda %*% t(Lambda) + diag(Psi)` is the non-spatial covariance.

```r
set.seed(16)
n <- 60
coords <- list(x = c(0, 0.5, 2), y = c(0, 0.5, 2))
p <- prod(lengths(coords))
q <- 2
beta <- c(2, 4)
alpha <- list(c(1, -0.30, 1.2), c(1, -0.45, 1.4))

Xi <- list()
for (g in 1:2) Xi[[g]] <- spatial_cov(coords, beta[g], alpha[[g]])

Lambda <- list(
  matrix(c(1, 0.6), ncol = 1),
  matrix(c(-0.5, 0.9), ncol = 1)
)
Psi <- list(c(0.5, 0.4), c(0.4, 0.7))
Omega <- list()
for (g in 1:2) Omega[[g]] <- tcrossprod(Lambda[[g]]) + diag(Psi[[g]])

M <- list(matrix(0, p, q), matrix(3, p, q))
z <- rbind(
  MASS::mvrnorm(n, as.vector(M[[1]]), kronecker(Omega[[1]], Xi[[1]])),
  MASS::mvrnorm(n, as.vector(M[[2]]), kronecker(Omega[[2]], Xi[[2]]))
)
X <- array(t(z), dim = c(p, q, 2 * n))

fit <- spatmix(
  X, G = 2, r = 1, coords = coords, nknots = 1, degree = 2,
  spatial_decay = "ispline", decay_group = 1,
  mean_structure = "constrained", max_iter = 15, tol = 0.01,
  spatial_max_iter = 20, verbose = FALSE
)

ord <- order(sapply(1:2, function(g) mean(fit$means[, , g])))
table(truth = rep(1:2, each = n), fitted = match(fit$cluster, ord))

d <- seq(0, 2, length.out = 200)
basis <- splines2::iSpline(
  d, knots = fit$knots[[1]], degree = 2, intercept = TRUE,
  Boundary.knots = range(d)
)
fitted_decay <- basis %*% fit$coordinate_beta[[1]][, ord]
plot(d, sigmoid(d, beta[1]), type = "l", lwd = 2,
     col = "firebrick", ylim = c(0, 1),
     xlab = "Normalized distance", ylab = "Decay")
lines(d, fitted_decay[, 1], col = "firebrick", lwd = 2, lty = 2)
lines(d, sigmoid(d, beta[2]), col = "navy", lwd = 2)
lines(d, fitted_decay[, 2], col = "navy", lwd = 2, lty = 2)
legend("bottomright",
       c("Component 1: sigmoid", "Component 1: I-spline",
         "Component 2: sigmoid", "Component 2: I-spline"),
       col = c("firebrick", "firebrick", "navy", "navy"),
       lty = c(1, 2, 1, 2), lwd = 2, bty = "n")
```

Point coordinates supplied as a numeric vector or matrix use one Euclidean
distance matrix. A list of coordinate vectors defines a complete Cartesian
grid with one spatial term per axis. For grid coordinates, use `decay_group`
to share decay parameters across axes. The default `1` places every axis in
one group, `decay_group = NULL` keeps the axes separate, and
`decay_group = c(1, 1, 2)` groups the first two axes separately from the third.
Axes in the same group share the decay parameters and the corresponding
spatial covariance coefficient.

Use `mean_structure = "constrained"` for a mean that is constant across spatial
locations. The default `mean_structure = "unconstrained"` estimates a separate
mean at each spatial location.

Use `init = "covariance"` when components may have similar means but different
covariance structures. It applies k-means to squared centered observations
before fitting the mixture.

The sigmoid uses distances normalized to `[0, 2]` and one positive decay
parameter per decay group. Bounds and the fixed shift can be changed with,
for example, `sigmoid_ctrl = list(init = 1, lower = 0.001, upper = 20,
shift = 3)`.
