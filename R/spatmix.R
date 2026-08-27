#' Fit a spatial mixture model
#'
#' Fits a spatial Gaussian mixture or a mixture with a non-spatial
#' factor-analyzer covariance structure to complete or incomplete data. The row
#' covariance is a linear spatial covariance model with decay represented by
#' monotone I-splines or a normalized sigmoid function. When `X` contains
#' missing values, conditional means and covariance blocks are updated by the
#' matrix partial EM (MPEM) algorithm.
#'
#' `coords` supports two spatial layouts. A numeric vector or matrix gives the
#' point coordinates of the `p` spatial locations and produces one Euclidean
#' distance matrix. A list of coordinate vectors, such as
#' `list(x = 1:8, y = 1:8)`, defines a complete Cartesian grid and produces one
#' axis-wise spatial term per list element. Grid locations must follow the row
#' order returned by `do.call(expand.grid, coords)`, with the first axis varying
#' fastest.
#'
#' @param X A numeric `p` by `q` by `n` array, whose dimensions correspond to
#'   spatial locations, non-spatial variables, and observations. Spatial-only
#'   data may instead be supplied as an `n` by `p` matrix, with observations in
#'   rows and spatial locations in columns. `NA` values invoke the MPEM updates.
#' @param G Number of mixture components.
#' @param r Number of column-covariance factors per component. Required for
#'   array input with a non-spatial dimension; use `NULL` or `0` for
#'   spatial-only data.
#' @param coords A numeric vector or `p`-row numeric matrix of point
#'   coordinates, or a list of numeric vectors defining the axes of a complete
#'   Cartesian grid. For grid coordinates, `prod(lengths(coords))` must equal
#'   `p`.
#' @param nknots Number of internal I-spline knots, given as a scalar or one
#'   value per grid axis. Internal knots are placed at empirical quantiles of
#'   the positive pairwise distances. Point coordinates use one set of knots.
#'   Ignored for sigmoid decay.
#' @param degree I-spline degree, given as a scalar or one value per grid axis.
#'   Point coordinates use one degree. Ignored for sigmoid decay.
#' @param common_noise If `TRUE`, use one row-noise parameter within each
#'   component; otherwise estimate a location-specific diagonal.
#' @param mean_structure Mean structure for each mixture component.
#'   `"unconstrained"` estimates a separate mean at each spatial location;
#'   `"constrained"` estimates a mean that is constant across spatial locations
#'   for each non-spatial variable.
#' @param max_iter Maximum number of outer EM iterations.
#' @param tol Absolute convergence tolerance for successive log-likelihood
#'   values.
#' @param seed Seed used for initialization.
#' @param init Initialization method. `"kmeans"` applies k-means to the
#'   observations, `"covariance"` applies k-means to squared centered
#'   observations for covariance-driven clusters, and `"random"` uses random
#'   initial responsibilities.
#' @param mpem_sweeps Number of conditional-moment coordinate sweeps per MPEM
#'   update. Used only when `X` contains missing values.
#' @param cov_warm_start Weight in `[0, 1]` assigned to the previous conditional
#'   covariance block at the start of each MPEM update; the remaining weight is
#'   assigned to a diagonal approximation. Used only when `X` contains missing
#'   values.
#' @param spatial_max_iter Maximum projected-gradient iterations for each
#'   I-spline spatial update. Ignored for sigmoid decay.
#' @param spatial_tol Tolerance for the I-spline or sigmoid decay update.
#' @param verbose If `TRUE`, print the log-likelihood at each iteration.
#' @param spatial_decay Spatial decay representation: `"ispline"` for monotone
#'   I-splines or `"sigmoid"` for a normalized one-parameter sigmoid function.
#' @param decay_group Labels specifying which grid-coordinate axes share decay
#'   parameters. The default `1` places all axes in one group; `NULL` gives each
#'   axis a distinct group; and, for example, `c(1, 1, 2)` groups the first two
#'   axes separately from the third. Within each mixture component, axes in the
#'   same group share both the decay parameters and the corresponding spatial
#'   covariance coefficient. Point coordinates produce one joint distance
#'   matrix, so `decay_group` has no nontrivial effect for that layout. Axes that
#'   share an I-spline decay must have the same number of basis functions.
#' @param sigmoid_ctrl Optional list with elements `init`, `lower`, `upper`, and
#'   `shift`. The first three may be scalars, values supplied per grid axis, or
#'   values supplied per decay group; `shift` must be a scalar. Sigmoid
#'   distances are normalized to `[0, 2]` separately for each axis. By default,
#'   the upper bound is determined from the nearest positive normalized
#'   distance.
#' @return An object of class `spatmixfit`. It is a list containing:
#'
#' * `call`: the matched function call.
#' * `converged`, `iterations`, and `log_likelihood`: convergence information
#'   and the log-likelihood sequence.
#' * `BIC` and `n_parameters`: the criterion
#'   `2 * log-likelihood - log(n) * n_parameters` and its parameter count. Larger
#'   `BIC` values indicate better models.
#' * `proportions`, `responsibility`, and `cluster`: fitted mixing proportions,
#'   an `n` by `G` responsibility matrix, and hard cluster assignments.
#' * `means`: a `p` by `q` by `G` array of component means.
#' * `row_covariance`, `row_precision`, `col_covariance`, and `col_precision`:
#'   fitted row and column covariance matrices and their precisions.
#' * `alpha`: spatial covariance coefficients, with one column per component.
#' * `coordinate_beta`: a list of fitted decay parameters, with one element per
#'   decay group.
#' * `sigmoid`: fitted sigmoid parameters when `spatial_decay = "sigmoid"`, and
#'   `NULL` otherwise. It is a matrix for one decay group and a list of matrices
#'   for multiple groups.
#' * `loading` and `uniqueness`: factor-analyzer loadings and uniquenesses for
#'   the non-spatial covariance.
#' * `imputation` and `has_missing`: the completed data in the same layout as
#'   the input and an indicator of whether `X` contained missing values.
#' * `knots` and `coordinate_dimensions`: the I-spline knots and spatial layout
#'   dimensions.
#' * `settings`: the model settings used for the fit.
#' @examples
#' sigmoid <- function(d, beta, shift = 3) {
#'   z0 <- plogis(-shift)
#'   (plogis(beta * d - shift) - z0) /
#'     (plogis(2 * beta - shift) - z0)
#' }
#'
#' spatial_cov <- function(coords, beta, alpha) {
#'   if (!is.list(coords)) coords <- list(coords)
#'   dims <- lengths(coords)
#'   decay <- lapply(seq_along(coords), function(j) {
#'     d <- as.matrix(dist(coords[[j]]))
#'     d <- 2 * d / max(d)
#'     matrices <- lapply(seq_along(coords), function(k) {
#'       if (j == k) sigmoid(d, beta) else matrix(1, dims[k], dims[k])
#'     })
#'     Reduce(kronecker, rev(matrices))
#'   })
#'   p <- prod(dims)
#'   J <- matrix(1, p, p) - diag(p)
#'   alpha[1] * J + alpha[2] * Reduce("+", decay) + alpha[3] * diag(p)
#' }
#'
#' ## Example 1
#' set.seed(1)
#' n <- 150
#' coords <- 1:10
#' beta <- c(3, 7)
#' alpha <- list(c(1, -0.30, 1.2), c(1, -0.45, 1.4))
#'
#' Xi <- list()
#' for (g in 1:2) Xi[[g]] <- spatial_cov(coords, beta[g], alpha[[g]])
#'
#' X <- rbind(
#'   MASS::mvrnorm(n, rep(0, 10), Xi[[1]]),
#'   MASS::mvrnorm(n, rep(2.5, 10), Xi[[2]])
#' )
#' truth <- rep(1:2, each = n)
#'
#' fit <- spatmix(
#'   X, G = 2, coords = coords, spatial_decay = "sigmoid",
#'   sigmoid_ctrl = list(init = 4, lower = 0.1, upper = 12),
#'   max_iter = 50, tol = 0.01, verbose = FALSE
#' )
#'
#' ord <- order(colMeans(fit$means[, 1, ]))
#' table(truth, fitted = match(fit$cluster, ord))
#' round(rbind(truth = beta, fitted = fit$sigmoid[1, ord]), 2)
#' round(cbind(
#'   truth.1 = alpha[[1]], fitted.1 = fit$alpha[, ord[1]],
#'   truth.2 = alpha[[2]], fitted.2 = fit$alpha[, ord[2]]
#' ), 2)
#'
#' ## Example 2
#' set.seed(16)
#' n <- 60
#' coords <- list(x = c(0, 0.5, 2), y = c(0, 0.5, 2))
#' p <- prod(lengths(coords))
#' q <- 2
#' beta <- c(2, 4)
#' alpha <- list(c(1, -0.30, 1.2), c(1, -0.45, 1.4))
#'
#' Xi <- list()
#' for (g in 1:2) Xi[[g]] <- spatial_cov(coords, beta[g], alpha[[g]])
#'
#' Lambda <- list(
#'   matrix(c(1, 0.6), ncol = 1),
#'   matrix(c(-0.5, 0.9), ncol = 1)
#' )
#' Psi <- list(c(0.5, 0.4), c(0.4, 0.7))
#' Omega <- list()
#' for (g in 1:2) Omega[[g]] <- tcrossprod(Lambda[[g]]) + diag(Psi[[g]])
#'
#' M <- list(matrix(0, p, q), matrix(3, p, q))
#' z <- rbind(
#'   MASS::mvrnorm(n, as.vector(M[[1]]), kronecker(Omega[[1]], Xi[[1]])),
#'   MASS::mvrnorm(n, as.vector(M[[2]]), kronecker(Omega[[2]], Xi[[2]]))
#' )
#' X <- array(t(z), dim = c(p, q, 2 * n))
#'
#' fit <- spatmix(
#'   X, G = 2, r = 1, coords = coords, nknots = 1, degree = 2,
#'   spatial_decay = "ispline", decay_group = 1,
#'   mean_structure = "constrained", max_iter = 15, tol = 0.01,
#'   spatial_max_iter = 20, verbose = FALSE
#' )
#'
#' ord <- order(sapply(1:2, function(g) mean(fit$means[, , g])))
#' table(truth = rep(1:2, each = n), fitted = match(fit$cluster, ord))
#'
#' d <- seq(0, 2, length.out = 200)
#' basis <- splines2::iSpline(
#'   d, knots = fit$knots[[1]], degree = 2, intercept = TRUE,
#'   Boundary.knots = range(d)
#' )
#' fitted_decay <- basis %*% fit$coordinate_beta[[1]][, ord]
#' plot(d, sigmoid(d, beta[1]), type = "l", lwd = 2,
#'      col = "firebrick", ylim = c(0, 1),
#'      xlab = "Normalized distance", ylab = "Decay")
#' lines(d, fitted_decay[, 1], col = "firebrick", lwd = 2, lty = 2)
#' lines(d, sigmoid(d, beta[2]), col = "navy", lwd = 2)
#' lines(d, fitted_decay[, 2], col = "navy", lwd = 2, lty = 2)
#' legend("bottomright",
#'        c("Component 1: sigmoid", "Component 1: I-spline",
#'          "Component 2: sigmoid", "Component 2: I-spline"),
#'        col = c("firebrick", "firebrick", "navy", "navy"),
#'        lty = c(1, 2, 1, 2), lwd = 2, bty = "n")
#' @export
spatmix <- function(X, G, r = NULL, coords, nknots = 6L, degree = 3L,
                 common_noise = TRUE,
                 mean_structure = c("unconstrained", "constrained"),
                 max_iter = 1000L, tol = 0.1,
                 seed = 1L, init = c("kmeans", "covariance", "random"),
                 mpem_sweeps = 1L,
                 cov_warm_start = 1,
                 spatial_max_iter = 5000L, spatial_tol = 1e-6,
                 verbose = interactive(),
                 spatial_decay = c("ispline", "sigmoid"),
                 decay_group = 1,
                 sigmoid_ctrl = list()) {
  call <- match.call()
  init <- match.arg(init)
  mean_structure <- match.arg(mean_structure)
  spatial_decay <- match.arg(spatial_decay)
  if (!is.numeric(X) || !(length(dim(X)) %in% c(2L, 3L))) {
    stop("X must be a numeric matrix or three-dimensional array.", call. = FALSE)
  }
  if (any(is.infinite(X))) stop("X may contain NA but not infinite values.", call. = FALSE)
  if (all(is.na(X))) stop("X must contain observed values.", call. = FALSE)
  input_is_matrix <- length(dim(X)) == 2L
  if (input_is_matrix) {
    X <- array(t(X), dim = c(ncol(X), 1L, nrow(X)))
  }
  dimensions <- dim(X)
  p <- dimensions[1L]
  q <- dimensions[2L]
  n <- dimensions[3L]
  spatial_only <- q == 1L
  G <- .positive_integer(G, "G")
  if (spatial_only) {
    if (!is.null(r)) {
      r <- .nonnegative_integer(r, "r")
      if (r != 0L) stop("r must be NULL or 0 for spatial-only data.", call. = FALSE)
    }
    r <- 0L
  } else {
    if (is.null(r)) stop("r is required when X has a non-spatial dimension.",
                            call. = FALSE)
    r <- .positive_integer(r, "r")
  }
  max_iter <- .positive_integer(max_iter, "max_iter")
  mpem_sweeps <- .positive_integer(mpem_sweeps, "mpem_sweeps")
  spatial_max_iter <- .positive_integer(spatial_max_iter, "spatial_max_iter")
  if (G >= n) stop("G must be smaller than the number of observations.", call. = FALSE)
  if (!spatial_only && r >= q) {
    stop("r must be smaller than dim(X)[2].", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("tol must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(spatial_tol) || length(spatial_tol) != 1L ||
      !is.finite(spatial_tol) || spatial_tol <= 0) {
    stop("spatial_tol must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(cov_warm_start) || length(cov_warm_start) != 1L ||
      !is.finite(cov_warm_start) || cov_warm_start < 0 || cov_warm_start > 1) {
    stop("cov_warm_start must be in [0, 1].", call. = FALSE)
  }
  if (!is.logical(common_noise) || length(common_noise) != 1L || is.na(common_noise)) {
    stop("common_noise must be TRUE or FALSE.", call. = FALSE)
  }

  coordinate <- .coordinate_basis(
    coords, p, nknots, degree, spatial_decay, decay_group, sigmoid_ctrl
  )
  number_axes <- coordinate$axes
  number_decay_groups <- coordinate$decay_groups
  J <- matrix(1, p, p)
  diag(J) <- 0
  I <- diag(p)
  has_missing <- anyNA(X)
  missing_info <- .missing_information(X)

  initial <- X
  cell_mean <- apply(X, c(1L, 2L), mean, na.rm = TRUE)
  overall_mean <- mean(X, na.rm = TRUE)
  cell_mean[!is.finite(cell_mean)] <- overall_mean
  initial[is.na(initial)] <- rep(cell_mean, n)[is.na(initial)]
  completed <- array(rep(initial, G), dim = c(p, q, n, G))

  if (init %in% c("kmeans", "covariance")) {
    flattened <- t(matrix(initial, nrow = p * q, ncol = n))
    if (init == "covariance") {
      flattened <- sweep(flattened, 2L, colMeans(flattened), "-")^2
    }
    cluster <- withr::with_seed(seed, stats::kmeans(
      flattened, centers = G, nstart = 5L
    )$cluster)
    responsibility <- matrix(0, n, G)
    responsibility[cbind(seq_len(n), cluster)] <- 1
  } else {
    responsibility <- withr::with_seed(seed, matrix(stats::runif(n * G), n, G))
    responsibility <- responsibility / rowSums(responsibility)
  }
  if (any(colSums(responsibility) <= 0)) {
    stop("Initialization produced an empty component; try another init method or seed.",
         call. = FALSE)
  }

  means <- array(NA_real_, c(p, q, G))
  row_covariance <- row_precision <- array(NA_real_, c(p, p, G))
  col_covariance <- col_precision <- array(NA_real_, c(q, q, G))
  loading <- array(NA_real_, c(q, r, G))
  uniqueness <- matrix(NA_real_, q, G)
  alpha <- matrix(NA_real_,
                  if (common_noise) number_decay_groups + 2L else
                    number_decay_groups + 1L + p,
                  G)
  beta <- lapply(coordinate$beta_initial, function(initial_beta) {
    matrix(rep(initial_beta, G), nrow = length(initial_beta), ncol = G)
  })
  conditional_cov <- if (has_missing) {
    lapply(seq_len(G), function(g) vector("list", n))
  } else NULL

  group_size <- colSums(responsibility)
  for (g in seq_len(G)) {
    weights <- responsibility[, g]
    completed_g <- array(completed[, , , g], dim = c(p, q, n))
    means[, , g] <- .component_mean(completed_g, weights, mean_structure)
    mean_g <- matrix(means[, , g], nrow = p, ncol = q)
    row_scatter <- .weighted_scatter(
      completed_g, mean_g, weights, diag(q), "row"
    ) / (q * group_size[g])
    if (!spatial_only) {
      normalizer <- row_scatter[1L, 1L]
      if (!is.finite(normalizer) || abs(normalizer) < 1e-10) normalizer <- 1
      row_scatter <- row_scatter / normalizer
    }

    beta_g <- lapply(beta, function(x) x[, g])
    model_matrices <- .spatial_model_matrices(beta_g, coordinate, J, I, common_noise)
    alpha[, g] <- .alpha_update(model_matrices, row_scatter, I, common_noise)
    if (spatial_decay == "sigmoid") {
      beta_g <- .beta_update(
        beta_g, alpha[, g], coordinate, I, row_scatter, common_noise,
        spatial_max_iter, spatial_tol
      )
    }
    for (group in seq_len(number_decay_groups)) beta[[group]][, g] <- beta_g[[group]]
    model_matrices <- .spatial_model_matrices(beta_g, coordinate, J, I, common_noise)
    row_covariance[, , g] <- .stabilize_covariance(
      .build_spatial_covariance(alpha[, g], model_matrices, p, common_noise)
    )
    row_precision[, , g] <- .safe_solve(row_covariance[, , g])

    if (spatial_only) {
      uniqueness[, g] <- 1
      col_covariance[, , g] <- col_precision[, , g] <- diag(q)
    } else {
      col_scatter <- .weighted_scatter(
        completed_g, mean_g, weights,
        matrix(row_precision[, , g], nrow = p, ncol = p), "col"
      ) / (p * group_size[g])
      factor <- .factor_initialize(col_scatter, r)
      covariance <- tcrossprod(factor$loading) + diag(factor$uniqueness)
      loading[, , g] <- factor$loading
      uniqueness[, g] <- factor$uniqueness
      col_covariance[, , g] <- covariance
      loading_g <- matrix(loading[, , g], nrow = q, ncol = r)
      col_precision[, , g] <- .woodbury_precision(loading_g, uniqueness[, g])
    }

    if (has_missing) {
      update <- .mpem_component_moments(
        X, mean_g,
        matrix(row_precision[, , g], nrow = p, ncol = p),
        matrix(col_precision[, , g], nrow = q, ncol = q),
        completed_g, conditional_cov[[g]], weights, mpem_sweeps,
        cov_warm_start, update_moments = TRUE
      )
      completed[, , , g] <- update$completed
      conditional_cov[[g]] <- update$conditional_cov
    }
  }

  log_likelihood <- numeric(max_iter)
  converged <- FALSE
  iterations <- max_iter
  for (iteration in seq_len(max_iter)) {
    group_size <- colSums(responsibility)
    if (any(group_size <= 1e-8)) {
      stop(
        "A mixture component became empty (effective sizes: ",
        paste(signif(group_size, 4L), collapse = ", "),
        "); reduce G or try another initialization method or seed.",
        call. = FALSE
      )
    }
    proportions <- group_size / n

    for (g in seq_len(G)) {
      weights <- responsibility[, g]
      completed_g <- array(completed[, , , g], dim = c(p, q, n))
      means[, , g] <- .component_mean(completed_g, weights, mean_structure)
      mean_g <- matrix(means[, , g], nrow = p, ncol = q)
      correction <- if (has_missing) {
        moments <- .mpem_component_moments(
          X, mean_g,
          matrix(row_precision[, , g], nrow = p, ncol = p),
          matrix(col_precision[, , g], nrow = q, ncol = q),
          completed_g, conditional_cov[[g]], weights, mpem_sweeps,
          cov_warm_start, update_moments = FALSE
        )
        list(row = moments$row_correction, col = moments$col_correction)
      } else list(row = matrix(0, p, p), col = matrix(0, q, q))

      row_scatter <- (
        .weighted_scatter(completed_g, mean_g, weights,
                          matrix(col_precision[, , g], q, q), "row") + correction$row
      ) / (q * group_size[g])
      if (!spatial_only) {
        normalizer <- row_scatter[1L, 1L]
        if (!is.finite(normalizer) || abs(normalizer) < 1e-10) normalizer <- 1
        row_scatter <- row_scatter / normalizer
      }

      beta_g <- lapply(beta, function(x) x[, g])
      precision_g <- matrix(row_precision[, , g], p, p)
      model_matrices <- .spatial_model_matrices(beta_g, coordinate, J, I, common_noise)
      alpha[, g] <- .alpha_update(
        model_matrices, row_scatter, precision_g, common_noise
      )
      # Let the Kronecker scale and column factors stabilize before releasing
      # the higher-dimensional I-spline shape parameters.
      if (spatial_decay == "sigmoid" || iteration > 5L) {
        beta_g <- .beta_update(
          beta_g, alpha[, g], coordinate, precision_g, row_scatter,
          common_noise, spatial_max_iter, spatial_tol
        )
      }
      for (group in seq_len(number_decay_groups)) beta[[group]][, g] <- beta_g[[group]]
      model_matrices <- .spatial_model_matrices(beta_g, coordinate, J, I, common_noise)
      row_covariance[, , g] <- .stabilize_covariance(
        .build_spatial_covariance(alpha[, g], model_matrices, p, common_noise)
      )
      row_precision[, , g] <- .safe_solve(row_covariance[, , g])

      if (!spatial_only) {
        col_scatter <- (
          .weighted_scatter(completed_g, mean_g, weights,
                            matrix(row_precision[, , g], p, p), "col") + correction$col
        ) / (p * group_size[g])
        loading_g <- matrix(loading[, , g], nrow = q, ncol = r)
        factor <- .factor_update(col_scatter, loading_g, uniqueness[, g])
        loading[, , g] <- factor$loading
        uniqueness[, g] <- factor$uniqueness
        col_covariance[, , g] <- factor$covariance
        col_precision[, , g] <- .woodbury_precision(factor$loading, factor$uniqueness)
      }

      if (has_missing) {
        update <- .mpem_component_moments(
          X, mean_g,
          matrix(row_precision[, , g], nrow = p, ncol = p),
          matrix(col_precision[, , g], nrow = q, ncol = q),
          completed_g, conditional_cov[[g]], weights, mpem_sweeps,
          cov_warm_start, update_moments = TRUE
        )
        completed[, , , g] <- update$completed
        conditional_cov[[g]] <- update$conditional_cov
      }
    }

    expectation <- .matrix_normal_estep(
      completed, proportions, means, row_precision, col_precision
    )
    responsibility <- expectation$responsibility
    log_likelihood[iteration] <- expectation$log_likelihood
    if (verbose) {
      cat("SpatMix iteration", iteration, ": log-likelihood =",
          format(log_likelihood[iteration], digits = 9L), "\n")
    }
    if (iteration >= 6L &&
        abs(log_likelihood[iteration] - log_likelihood[iteration - 1L]) < tol) {
      converged <- TRUE
      iterations <- iteration
      break
    }
  }
  log_likelihood <- log_likelihood[seq_len(iterations)]
  group_size <- colSums(responsibility)
  proportions <- group_size / n
  classification <- max.col(responsibility, ties.method = "first")
  imputed <- initial
  if (has_missing) {
    for (i in seq_len(n)) {
      if (!length(missing_info[[i]]$index)) next
      mixture_completed <- matrix(completed[, , i, ], nrow = p * q, ncol = G)
      values <- as.vector(mixture_completed %*% responsibility[i, ])
      imputed[, , i][missing_info[[i]]$index] <- values[missing_info[[i]]$index]
    }
  }
  if (input_is_matrix) imputed <- t(matrix(imputed, nrow = p, ncol = n))

  mean_parameters <- if (mean_structure == "unconstrained") p * q else q
  spatial_parameters <- if (common_noise) number_decay_groups + 2L else
    number_decay_groups + 1L + p
  decay_parameters <- if (spatial_decay == "ispline") {
    sum(coordinate$group_basis_size - 1L)
  } else {
    number_decay_groups
  }
  factor_parameters <- if (spatial_only) 0 else
    q * r + q - r * (r - 1) / 2
  number_parameters <- G - 1L + G *
    (mean_parameters + spatial_parameters + decay_parameters + factor_parameters)
  bic <- 2 * utils::tail(log_likelihood, 1L) - log(n) * number_parameters

  sigmoid <- if (spatial_decay == "sigmoid") {
    if (number_decay_groups == 1L) beta[[1L]] else beta
  } else NULL

  out <- list(
    call = call, converged = converged, iterations = iterations,
    log_likelihood = log_likelihood, BIC = bic, n_parameters = number_parameters,
    proportions = proportions, means = means,
    row_covariance = row_covariance, row_precision = row_precision,
    col_covariance = col_covariance, col_precision = col_precision,
    alpha = alpha, coordinate_beta = beta, sigmoid = sigmoid, loading = loading,
    uniqueness = uniqueness, responsibility = responsibility,
    cluster = classification, imputation = imputed, has_missing = has_missing,
    knots = coordinate$knots, coordinate_dimensions = coordinate$dimensions,
    settings = list(G = G, r = r, init = init, seed = seed,
                    common_noise = common_noise,
                    mean_structure = mean_structure, spatial_only = spatial_only,
                    input_is_matrix = input_is_matrix,
                    spatial_decay = spatial_decay,
                    decay_group = coordinate$decay_group,
                    decay_group_labels = coordinate$decay_group_labels,
                    sigmoid_ctrl = coordinate$sigmoid_ctrl,
                    mpem_sweeps = mpem_sweeps,
                    cov_warm_start = cov_warm_start)
  )
  class(out) <- "spatmixfit"
  out
}

#' @export
print.spatmixfit <- function(x, ...) {
  cat(if (x$settings$spatial_only) "Spatial Gaussian mixture model\n" else
    "Spatial mixture model with a factor-analyzer covariance\n")
  cat("  components:", x$settings$G, "\n")
  if (!x$settings$spatial_only) cat("  factors:", x$settings$r, "\n")
  cat("  spatial decay:", x$settings$spatial_decay, "\n")
  cat("  axis decay groups:", paste(x$settings$decay_group, collapse = ", "), "\n")
  cat("  mean structure:", x$settings$mean_structure, "\n")
  cat("  observations:", nrow(x$responsibility), "\n")
  cat("  missing-data MPEM:", if (x$has_missing) "yes" else "no", "\n")
  cat("  iterations:", x$iterations,
      if (x$converged) "(converged)" else "(maximum reached)", "\n")
  cat("  BIC (higher is better):", format(x$BIC, digits = 7L), "\n")
  invisible(x)
}

#' @export
summary.spatmixfit <- function(object, ...) {
  list(
    converged = object$converged,
    iterations = object$iterations,
    G = object$settings$G,
    r = object$settings$r,
    spatial_only = object$settings$spatial_only,
    spatial_decay = object$settings$spatial_decay,
    decay_group = object$settings$decay_group,
    sigmoid = object$sigmoid,
    mean_structure = object$settings$mean_structure,
    proportions = object$proportions,
    log_likelihood = utils::tail(object$log_likelihood, 1L),
    BIC = object$BIC,
    n_parameters = object$n_parameters,
    cluster_sizes = tabulate(object$cluster, nbins = object$settings$G),
    used_missing_data_MPEM = object$has_missing
  )
}
