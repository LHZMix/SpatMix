.simplex_project <- function(y) {
  u <- sort(y, decreasing = TRUE)
  index <- seq_along(u)
  cumulative <- cumsum(u)
  candidates <- index[(u + (1 - cumulative) / index) > 0]
  rho <- if (length(candidates)) max(candidates) else 1L
  pmax(y + (1 - cumulative[rho]) / rho, 0)
}

.quantile_knots <- function(distance, number) {
  number <- as.integer(number)
  if (number <= 0L) return(numeric())
  values <- distance[is.finite(distance)]
  lower <- min(values)
  upper <- max(values)
  if (!is.finite(lower) || !is.finite(upper) || upper <= lower) {
    stop("Each coordinate axis must contain at least two distinct locations.", call. = FALSE)
  }
  positive <- values[values > lower]
  source <- if (length(unique(positive)) >= 2L) positive else values
  probability <- seq(0, 1, length.out = number + 2L)[-c(1L, number + 2L)]
  knots <- sort(unique(as.numeric(stats::quantile(
    source, probability, na.rm = TRUE, names = FALSE, type = 7
  ))))
  knots <- knots[knots > lower & knots < upper]
  if (!length(knots)) {
    knots <- seq(lower, upper, length.out = number + 2L)[-c(1L, number + 2L)]
  }
  knots
}

.sigmoid_decay <- function(distance, beta, shift = 3) {
  baseline <- stats::plogis(-shift)
  denominator <- stats::plogis(2 * beta - shift) - baseline
  if (!is.finite(denominator) || abs(denominator) < 1e-12) {
    return(distance / 2)
  }
  (stats::plogis(distance * beta - shift) - baseline) / denominator
}

.automatic_sigmoid_upper <- function(distance, lower, shift) {
  positive <- distance[upper.tri(distance) & distance > 0]
  if (!length(positive)) return(max(10, lower * 10))
  nearest <- min(positive)
  midpoint <- function(beta) .sigmoid_decay(nearest, beta, shift) - 0.5
  if (midpoint(lower) >= 0) return(max(10, lower * 10))
  upper <- max(1, lower * 2)
  while (upper < 1e4 && midpoint(upper) < 0) upper <- min(upper * 2, 1e4)
  if (midpoint(upper) < 0) return(upper)
  stats::uniroot(midpoint, c(lower, upper))$root
}

.decay_groups <- function(decay_group, axes) {
  if (is.null(decay_group)) {
    return(list(axis = seq_len(axes), labels = as.character(seq_len(axes)),
                count = axes))
  }
  if (is.list(decay_group) || is.matrix(decay_group)) {
    stop("decay_group must be NULL or an atomic vector.", call. = FALSE)
  }
  if (length(decay_group) == 1L) decay_group <- rep(decay_group, axes)
  if (length(decay_group) != axes) {
    stop("decay_group must have length 1 or one value per coordinate axis.",
         call. = FALSE)
  }
  if (anyNA(decay_group)) stop("decay_group must not contain NA.", call. = FALSE)
  labels <- as.character(decay_group)
  if (any(!nzchar(labels))) {
    stop("decay_group labels must not be empty.", call. = FALSE)
  }
  group_labels <- unique(labels)
  list(axis = match(labels, group_labels), labels = group_labels,
       count = length(group_labels))
}

.sigmoid_ctrls <- function(control, distances, axis_group) {
  if (!is.list(control)) stop("sigmoid_ctrl must be a list.", call. = FALSE)
  allowed <- c("init", "lower", "upper", "shift")
  unknown <- setdiff(names(control), allowed)
  if (length(unknown)) {
    stop("Unknown sigmoid_ctrl field: ", paste(unknown, collapse = ", "), ".",
         call. = FALSE)
  }
  axes <- length(distances)
  groups <- max(axis_group)
  expand_axis <- function(value, name) {
    if (!is.numeric(value) || any(!is.finite(value)) ||
        !(length(value) %in% c(1L, axes, groups))) {
      stop("sigmoid_ctrl$", name,
           " must contain finite values and have length 1, one value per ",
           "coordinate axis, or one value per decay group.",
           call. = FALSE)
    }
    if (length(value) == groups && groups != axes) return(value[axis_group])
    rep(value, length.out = axes)
  }

  shift <- if (is.null(control$shift)) 3 else control$shift
  if (!is.numeric(shift) || length(shift) != 1L || !is.finite(shift)) {
    stop("sigmoid_ctrl$shift must be a finite scalar.", call. = FALSE)
  }
  lower_axis <- expand_axis(if (is.null(control$lower)) 1e-3 else control$lower,
                            "lower")
  if (any(lower_axis <= 0)) {
    stop("sigmoid_ctrl$lower must be positive.", call. = FALSE)
  }
  upper_axis <- if (is.null(control$upper)) {
    vapply(seq_len(axes), function(j) {
      .automatic_sigmoid_upper(distances[[j]], lower_axis[j], shift)
    }, numeric(1L))
  } else {
    expand_axis(control$upper, "upper")
  }
  lower <- vapply(seq_len(groups), function(j) {
    max(lower_axis[axis_group == j])
  }, numeric(1L))
  upper <- vapply(seq_len(groups), function(j) {
    min(upper_axis[axis_group == j])
  }, numeric(1L))
  if (any(upper <= lower)) {
    stop("sigmoid_ctrl$upper must exceed sigmoid_ctrl$lower.", call. = FALSE)
  }
  initial_axis <- expand_axis(if (is.null(control$init)) 1 else control$init,
                              "init")
  initial <- vapply(seq_len(groups), function(j) {
    mean(initial_axis[axis_group == j])
  }, numeric(1L))
  initial <- pmin(pmax(initial, lower), upper)
  list(init = initial, lower = lower, upper = upper, shift = shift)
}

.axis_term <- function(axis_matrix, axis, dimensions) {
  matrices <- lapply(seq_along(dimensions), function(j) {
    if (j == axis) axis_matrix else matrix(1, dimensions[j], dimensions[j])
  })
  Reduce(kronecker, rev(matrices))
}

.coordinate_basis <- function(coords, p, nknots, degree,
                              spatial_decay = c("ispline", "sigmoid"),
                              decay_group = 1, sigmoid_ctrl = list()) {
  spatial_decay <- match.arg(spatial_decay)
  if (is.list(coords)) {
    dimensions <- vapply(coords, length, integer(1L))
    if (prod(dimensions) != p) {
      stop("For grid coordinates, prod(lengths(coords)) must equal dim(X)[1].",
           call. = FALSE)
    }
    axes <- coords
    lifted <- TRUE
  } else {
    coordinate_matrix <- as.matrix(coords)
    if (nrow(coordinate_matrix) != p) {
      if (is.vector(coords) && length(coords) == p) {
        coordinate_matrix <- matrix(coords, ncol = 1L)
      } else {
        stop("Point coordinates must have one row per spatial location.", call. = FALSE)
      }
    }
    axes <- list(coordinate_matrix)
    dimensions <- p
    lifted <- FALSE
  }

  number_axes <- length(axes)
  decay <- .decay_groups(decay_group, number_axes)
  axis_group <- decay$axis
  number_decay_groups <- decay$count
  if (spatial_decay == "ispline") {
    if (length(nknots) == 1L) nknots <- rep(nknots, number_axes)
    if (length(degree) == 1L) degree <- rep(degree, number_axes)
    if (length(nknots) != number_axes || length(degree) != number_axes) {
      stop("nknots and degree must have length 1 or one value per coordinate axis.",
           call. = FALSE)
    }
    if (any(nknots < 0 | nknots != as.integer(nknots))) {
      stop("nknots must contain non-negative integers.", call. = FALSE)
    }
    if (any(degree < 1 | degree != as.integer(degree))) {
      stop("degree must contain positive integers.", call. = FALSE)
    }
  } else {
    nknots <- degree <- integer(number_axes)
  }

  distance_matrix <- lapply(axes, function(axis) as.matrix(stats::dist(axis)))
  if (spatial_decay == "ispline") {
    distance <- lapply(distance_matrix, as.vector)
    knots <- lapply(seq_len(number_axes), function(j) {
      .quantile_knots(distance[[j]], nknots[j])
    })
    splines <- lapply(seq_len(number_axes), function(j) {
      range_j <- range(distance[[j]])
      splines2::iSpline(
        distance[[j]], knots = knots[[j]], degree = degree[j], intercept = TRUE,
        Boundary.knots = range_j
      )
    })
    basis_size <- vapply(splines, ncol, integer(1L))
    group_basis_size <- vapply(seq_len(number_decay_groups), function(group) {
      sizes <- unique(basis_size[axis_group == group])
      if (length(sizes) != 1L) {
        stop("Axes that share an I-spline decay must have the same number ",
             "of basis functions. Adjust nknots, degree, or decay_group.",
             call. = FALSE)
      }
      sizes
    }, integer(1L))

    if (lifted) {
      basis <- lapply(seq_len(number_axes), function(axis) {
        size <- dimensions[axis]
        lapply(seq_len(ncol(splines[[axis]])), function(k) {
          .axis_term(matrix(splines[[axis]][, k], size, size), axis, dimensions)
        })
      })
    } else {
      basis <- list(lapply(seq_len(ncol(splines[[1L]])), function(k) {
        matrix(splines[[1L]][, k], p, p)
      }))
    }
    beta_initial <- lapply(group_basis_size, function(size) rep(1 / size, size))
    sigmoid <- NULL
  } else {
    distance <- lapply(distance_matrix, function(x) {
      maximum <- max(x)
      if (!is.finite(maximum) || maximum <= 0) {
        stop("Each coordinate axis must contain at least two distinct locations.",
             call. = FALSE)
      }
      2 * x / maximum
    })
    sigmoid <- .sigmoid_ctrls(sigmoid_ctrl, distance, axis_group)
    knots <- replicate(number_axes, numeric(), simplify = FALSE)
    splines <- basis <- NULL
    group_basis_size <- rep(1L, number_decay_groups)
    beta_initial <- lapply(sigmoid$init, function(x) x)
  }

  list(
    splines = splines, basis = basis, dimensions = dimensions,
    knots = knots, actual_nknots = vapply(knots, length, integer(1L)),
    degree = as.integer(degree), axes = number_axes, lifted = lifted,
    decay_group = axis_group, decay_group_labels = decay$labels,
    decay_groups = number_decay_groups, group_basis_size = group_basis_size,
    spatial_decay = spatial_decay, distance = distance,
    sigmoid_ctrl = sigmoid, beta_initial = beta_initial
  )
}

.weighted_spatial_term <- function(beta, basis) {
  out <- matrix(0, nrow(basis[[1L]]), ncol(basis[[1L]]))
  for (k in seq_along(beta)) out <- out + beta[k] * basis[[k]]
  out
}

.decay_terms <- function(beta, coordinate) {
  if (coordinate$spatial_decay == "ispline") {
    return(lapply(seq_len(coordinate$axes), function(axis) {
      group <- coordinate$decay_group[axis]
      .weighted_spatial_term(beta[[group]], coordinate$basis[[axis]])
    }))
  }
  lapply(seq_len(coordinate$axes), function(axis) {
    group <- coordinate$decay_group[axis]
    term <- .sigmoid_decay(
      coordinate$distance[[axis]], beta[[group]][1L],
      coordinate$sigmoid_ctrl$shift
    )
    if (coordinate$lifted) {
      .axis_term(term, axis, coordinate$dimensions)
    } else {
      term
    }
  })
}

.group_decay_terms <- function(beta, coordinate) {
  axis_terms <- .decay_terms(beta, coordinate)
  lapply(seq_len(coordinate$decay_groups), function(group) {
    Reduce(`+`, axis_terms[coordinate$decay_group == group])
  })
}

.spatial_model_matrices <- function(beta, coordinate, J, I, common_noise) {
  matrices <- c(list(J), .group_decay_terms(beta, coordinate))
  if (common_noise) matrices <- c(matrices, list(I))
  matrices
}

.build_spatial_covariance <- function(alpha, matrices, p, common_noise) {
  covariance <- matrix(0, p, p)
  for (j in seq_along(matrices)) covariance <- covariance + alpha[j] * matrices[[j]]
  if (!common_noise) {
    start <- length(matrices) + 1L
    diag(covariance) <- diag(covariance) + alpha[start:(start + p - 1L)]
  }
  (covariance + t(covariance)) / 2
}

.alpha_update <- function(matrices, scatter, precision, common_noise) {
  p <- nrow(precision)
  base <- length(matrices)
  number <- if (common_noise) base else base + p
  transformed_scatter <- precision %*% scatter %*% precision
  rhs <- numeric(number)
  theta <- matrix(0, number, number)
  transformed <- vector("list", base)

  for (j in seq_len(base)) {
    transformed[[j]] <- precision %*% matrices[[j]] %*% precision
    rhs[j] <- sum(matrices[[j]] * transformed_scatter)
    for (i in seq_len(j)) {
      theta[i, j] <- theta[j, i] <- sum(matrices[[i]] * transformed[[j]])
    }
  }
  if (!common_noise) {
    diagonal_index <- (base + 1L):(base + p)
    for (j in seq_len(base)) {
      values <- diag(transformed[[j]])
      theta[j, diagonal_index] <- values
      theta[diagonal_index, j] <- values
    }
    theta[diagonal_index, diagonal_index] <- precision^2
    rhs[diagonal_index] <- diag(transformed_scatter)
  }
  as.vector(MASS::ginv(theta) %*% rhs)
}

.group_beta_gradient <- function(beta, basis, axes, alpha_group, precision, target) {
  group_term <- Reduce(`+`, lapply(axes, function(axis) {
    .weighted_spatial_term(beta, basis[[axis]])
  }))
  residual <- target - alpha_group * group_term
  middle <- precision %*% residual %*% precision
  -vapply(seq_along(beta), function(k) {
    alpha_group * sum(vapply(
      axes, function(axis) sum(basis[[axis]][[k]] * middle), numeric(1L)
    ))
  }, numeric(1L))
}

.group_beta_objective <- function(beta, basis, axes, alpha_group, precision,
                                  target) {
  group_term <- Reduce(`+`, lapply(axes, function(axis) {
    .weighted_spatial_term(beta, basis[[axis]])
  }))
  residual <- target - alpha_group * group_term
  0.5 * sum(residual * (precision %*% residual %*% precision))
}

.spatial_covariance_objective <- function(covariance, scatter) {
  covariance <- (covariance + t(covariance)) / 2
  log_determinant <- determinant(covariance, logarithm = TRUE)
  if (log_determinant$sign <= 0 || !is.finite(as.numeric(log_determinant$modulus))) {
    return(.Machine$double.xmax / 100)
  }
  precision <- tryCatch(solve(covariance), error = function(e) NULL)
  if (is.null(precision) || any(!is.finite(precision))) {
    return(.Machine$double.xmax / 100)
  }
  as.numeric(log_determinant$modulus) + sum(scatter * precision)
}

.sigmoid_beta_update <- function(beta, alpha, coordinate, scatter,
                                 common_noise, tolerance) {
  p <- nrow(scatter)
  J <- matrix(1, p, p)
  diag(J) <- 0
  I <- diag(p)
  for (group in seq_along(beta)) {
    objective <- function(value) {
      candidate <- beta
      candidate[[group]] <- value
      matrices <- .spatial_model_matrices(candidate, coordinate, J, I, common_noise)
      covariance <- .build_spatial_covariance(alpha, matrices, p, common_noise)
      .spatial_covariance_objective(covariance, scatter)
    }
    result <- tryCatch(stats::optimize(
      objective,
      interval = c(coordinate$sigmoid_ctrl$lower[group],
                   coordinate$sigmoid_ctrl$upper[group]),
      tol = max(tolerance, 1e-8)
    ), error = function(e) NULL)
    if (!is.null(result) && is.finite(result$objective)) beta[[group]] <- result$minimum
  }
  beta
}

.beta_update <- function(beta, alpha, coordinate, precision, scatter,
                         common_noise, max_iter, tolerance) {
  if (coordinate$spatial_decay == "sigmoid") {
    return(.sigmoid_beta_update(
      beta, alpha, coordinate, scatter, common_noise, tolerance
    ))
  }
  basis <- coordinate$basis
  p <- nrow(precision)
  number_groups <- length(beta)
  J <- matrix(1, p, p)
  diag(J) <- 0
  I <- diag(p)

  for (group in seq_len(number_groups)) {
    terms <- .group_decay_terms(beta, coordinate)
    group_axes <- which(coordinate$decay_group == group)
    target <- scatter - alpha[1L] * J
    for (other in setdiff(seq_len(number_groups), group)) {
      target <- target - alpha[other + 1L] * terms[[other]]
    }
    if (common_noise) {
      target <- target - alpha[number_groups + 2L] * I
    } else {
      diag(target) <- diag(target) -
        alpha[(number_groups + 2L):(number_groups + 1L + p)]
    }
    current <- beta[[group]]
    spatial_scale <- alpha[group + 1L]
    if (!is.finite(spatial_scale) || abs(spatial_scale) < 1e-10) next
    gradient <- .group_beta_gradient(
      current, basis, group_axes, spatial_scale, precision, target
    )
    objective <- .group_beta_objective(
      current, basis, group_axes, spatial_scale, precision, target
    )
    basis_norm <- sum(vapply(group_axes, function(axis) {
      sum(vapply(basis[[axis]], function(x) sum(x^2), numeric(1L)))
    }, numeric(1L)))
    lipschitz <- max(sum(diag(precision))^2 * spatial_scale^2 * basis_norm, 1e-8)
    step <- 1 / lipschitz
    for (iteration in seq_len(max_iter)) {
      candidate <- .simplex_project(current - step * gradient)
      candidate_objective <- .group_beta_objective(
        candidate, basis, group_axes, spatial_scale, precision, target
      )
      backtracks <- 0L
      while ((!is.finite(candidate_objective) || candidate_objective > objective) &&
             backtracks < 30L) {
        step <- step / 2
        candidate <- .simplex_project(current - step * gradient)
        candidate_objective <- .group_beta_objective(
          candidate, basis, group_axes, spatial_scale, precision, target
        )
        backtracks <- backtracks + 1L
      }
      difference <- candidate - current
      new_gradient <- .group_beta_gradient(
        candidate, basis, group_axes, spatial_scale, precision, target
      )
      gradient_difference <- new_gradient - gradient
      denominator <- sum(difference * gradient_difference)
      step <- if (is.finite(denominator) && denominator > .Machine$double.eps) {
        sum(difference^2) / denominator
      } else {
        1 / lipschitz
      }
      current <- candidate
      gradient <- new_gradient
      objective <- candidate_objective
      if (iteration > 1L && sqrt(sum(difference^2)) < tolerance) break
    }
    beta[[group]] <- current
  }
  beta
}

.matrix_normal_estep <- function(data, proportions, means, row_precision,
                                 col_precision) {
  p <- dim(data)[1L]
  q <- dim(data)[2L]
  n <- dim(data)[3L]
  G <- length(proportions)
  log_weight <- matrix(NA_real_, n, G)
  for (g in seq_len(G)) {
    row_precision_g <- matrix(row_precision[, , g], nrow = p, ncol = p)
    col_precision_g <- matrix(col_precision[, , g], nrow = q, ncol = q)
    det_row <- determinant(row_precision_g, logarithm = TRUE)
    det_col <- determinant(col_precision_g, logarithm = TRUE)
    if (det_row$sign <= 0 || det_col$sign <= 0) {
      stop("A precision matrix is not positive definite.", call. = FALSE)
    }
    mean_g <- matrix(means[, , g], nrow = p, ncol = q)
    centered <- sweep(data[, , , g, drop = FALSE], c(1L, 2L), mean_g, "-")
    centered <- array(centered, c(p, q, n))
    quadratic <- vapply(seq_len(n), function(i) {
      centered_i <- matrix(centered[, , i], nrow = p, ncol = q)
      sum((row_precision_g %*% centered_i %*% col_precision_g) * centered_i)
    }, numeric(1L))
    log_density <- -0.5 * (p * q * log(2 * pi) -
      q * as.numeric(det_row$modulus) - p * as.numeric(det_col$modulus) + quadratic)
    log_weight[, g] <- log(proportions[g]) + log_density
  }
  likelihood <- sum(matrixStats::rowLogSumExps(log_weight))
  maximum <- apply(log_weight, 1L, max)
  responsibility <- exp(log_weight - maximum)
  responsibility <- responsibility / rowSums(responsibility)
  list(responsibility = responsibility, log_likelihood = likelihood)
}

.factor_initialize <- function(scatter, r, epsilon = 1e-8) {
  decomposition <- eigen((scatter + t(scatter)) / 2, symmetric = TRUE)
  loading <- decomposition$vectors[, seq_len(r), drop = FALSE] %*%
    diag(sqrt(pmax(decomposition$values[seq_len(r)], epsilon)), r)
  uniqueness <- pmax(diag(scatter - tcrossprod(loading)), epsilon)
  list(loading = loading, uniqueness = uniqueness)
}

.factor_update <- function(scatter, loading, uniqueness, epsilon = 1e-8) {
  r <- ncol(loading)
  scaled_loading <- loading * (1 / uniqueness)
  inverse_w <- .safe_solve(diag(r) + crossprod(loading, scaled_loading))
  regression <- scaled_loading %*% inverse_w
  scatter_regression <- scatter %*% regression
  new_loading <- scatter_regression %*%
    .safe_solve(inverse_w + crossprod(regression, scatter_regression))
  new_uniqueness <- pmax(diag(scatter) - rowSums(new_loading * scatter_regression), epsilon)
  covariance <- tcrossprod(new_loading) + diag(new_uniqueness)
  list(
    loading = new_loading,
    uniqueness = new_uniqueness,
    covariance = covariance
  )
}

.woodbury_precision <- function(loading, uniqueness) {
  r <- ncol(loading)
  inverse_uniqueness <- 1 / uniqueness
  scaled <- loading * inverse_uniqueness
  diag(inverse_uniqueness) - scaled %*%
    .safe_solve(diag(r) + crossprod(loading, scaled)) %*% t(scaled)
}

.missing_information <- function(X) {
  p <- dim(X)[1L]
  lapply(seq_len(dim(X)[3L]), function(i) {
    missing <- which(is.na(X[, , i]))
    list(
      index = missing,
      row = ((missing - 1L) %% p) + 1L,
      col = ((missing - 1L) %/% p) + 1L
    )
  })
}

.mpem_component_moments <- function(X, mean, row_precision, col_precision,
                                    completed, conditional_cov, weights,
                                    sweeps, warm_start,
                                    update_moments = TRUE) {
  p <- dim(X)[1L]
  q <- dim(X)[2L]
  n <- dim(X)[3L]
  result <- mpem_moments(
    x = t(matrix(X, nrow = p * q, ncol = n)),
    rows = p,
    cols = q,
    mean = mean,
    row_precision = row_precision,
    col_precision = col_precision,
    precision_multiplier = 1,
    completed = t(matrix(completed, nrow = p * q, ncol = n)),
    conditional_cov = conditional_cov,
    weights = weights,
    sweeps = sweeps,
    warm_start = warm_start,
    update_moments = update_moments
  )
  list(
    completed = array(t(result$completed), dim = c(p, q, n)),
    conditional_cov = result$conditional_cov,
    row_correction = result$row_correction,
    col_correction = result$col_correction
  )
}

.weighted_mean <- function(data, weights) {
  total <- matrix(0, dim(data)[1L], dim(data)[2L])
  for (i in seq_len(dim(data)[3L])) total <- total + weights[i] * data[, , i]
  total / sum(weights)
}

.component_mean <- function(data, weights, structure) {
  estimate <- .weighted_mean(data, weights)
  if (structure == "constrained") {
    estimate <- matrix(rep(colMeans(estimate), each = nrow(estimate)),
                       nrow(estimate), ncol(estimate))
  }
  estimate
}

.weighted_scatter <- function(data, mean, weights, precision, side = c("row", "col")) {
  side <- match.arg(side)
  size <- if (side == "row") nrow(mean) else ncol(mean)
  total <- matrix(0, size, size)
  for (i in seq_len(dim(data)[3L])) {
    residual <- data[, , i] - mean
    term <- if (side == "row") residual %*% precision %*% t(residual) else
      t(residual) %*% precision %*% residual
    total <- total + weights[i] * term
  }
  total
}
