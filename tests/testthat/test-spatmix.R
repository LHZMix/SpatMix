test_that("spatmix fits complete data and grid coordinates", {
  set.seed(5)
  x <- array(rnorm(4 * 3 * 30), c(4, 3, 30))
  x[, , 16:30] <- x[, , 16:30] + 1.5
  fit <- spatmix(
    x, G = 2, r = 1, coords = list(x = 1:2, y = 1:2),
    nknots = 0, degree = 1, max_iter = 6, spatial_max_iter = 20,
    verbose = FALSE
  )

  expect_s3_class(fit, "spatmixfit")
  expect_identical(fit$imputation, x)
  expect_equal(rowSums(fit$responsibility), rep(1, 30), tolerance = 1e-8)
  expect_length(fit$cluster, 30)
  expect_false(fit$has_missing)
})

test_that("spatmix MPEM completes missing values", {
  set.seed(6)
  x <- array(rnorm(4 * 3 * 30), c(4, 3, 30))
  x[, , 16:30] <- x[, , 16:30] + 1.5
  incomplete <- x
  incomplete[1, 1, seq(1, 30, by = 5)] <- NA_real_
  incomplete[4, 3, seq(3, 30, by = 6)] <- NA_real_
  fit <- spatmix(
    incomplete, G = 2, r = 1, coords = 1:4,
    nknots = 1, degree = 1, max_iter = 6, spatial_max_iter = 20,
    verbose = FALSE
  )

  expect_true(fit$has_missing)
  expect_false(anyNA(fit$imputation))
  expect_equal(fit$imputation[!is.na(incomplete)], incomplete[!is.na(incomplete)])
})

test_that("constrained mean structure is enforced", {
  set.seed(8)
  x <- array(rnorm(4 * 3 * 24), c(4, 3, 24))
  fit <- spatmix(
    x, G = 2, r = 1, coords = 1:4, mean_structure = "constrained",
    nknots = 0, degree = 1, max_iter = 6, spatial_max_iter = 10,
    verbose = FALSE
  )
  for (g in 1:2) {
    expect_equal(fit$means[, , g],
                 matrix(rep(fit$means[1, , g], each = 4), 4, 3))
  }
})

test_that("covariance initialization handles equal-mean variance clusters", {
  set.seed(81)
  x <- array(rnorm(6 * 3 * 80), c(6, 3, 80))
  x[, , 41:80] <- 3 * x[, , 41:80]
  fit <- spatmix(
    x, G = 2, r = 1, coords = 1:6, init = "covariance",
    nknots = 0, degree = 1, mean_structure = "constrained",
    max_iter = 8, spatial_max_iter = 20, verbose = FALSE
  )

  expect_identical(fit$settings$init, "covariance")
  expect_equal(rowSums(fit$responsibility), rep(1, 80), tolerance = 1e-8)
  expect_true(all(diff(fit$log_likelihood) >= -1e-7))
  tab <- table(rep(1:2, each = 40), fit$cluster)
  expect_equal(max(sum(diag(tab)), sum(diag(tab[, 2:1]))), 80)
})

test_that("spatial-only matrix data skips the column factor analyzer", {
  set.seed(9)
  x <- matrix(rnorm(40 * 4), nrow = 40, ncol = 4)
  x[21:40, ] <- x[21:40, ] + 1.5
  fit <- spatmix(
    x, G = 2, coords = 1:4, spatial_decay = "sigmoid",
    sigmoid_ctrl = list(init = 1, lower = 0.001, upper = 10, shift = 3),
    max_iter = 6, verbose = FALSE
  )

  expect_true(fit$settings$spatial_only)
  expect_identical(fit$settings$r, 0L)
  expect_identical(dim(fit$loading), c(1L, 0L, 2L))
  expect_identical(dim(fit$imputation), dim(x))
  expect_equal(fit$imputation, x)
  expect_identical(dim(fit$sigmoid), c(1L, 2L))
  expect_identical(fit$settings$sigmoid_ctrl$upper, 10)
  expect_true(all(fit$sigmoid > 0))
  expect_true(all(vapply(seq_len(2), function(g) {
    min(eigen(fit$row_covariance[, , g], symmetric = TRUE,
              only.values = TRUE)$values) > 0
  }, logical(1L))))
  expect_identical(fit$n_parameters, 17)
  expect_error(
    spatmix(x, G = 2, r = 1, coords = 1:4, max_iter = 1, verbose = FALSE),
    "NULL or 0"
  )
})

test_that("spatial-only sigmoid data uses MPEM for missing entries", {
  set.seed(10)
  x <- matrix(rnorm(36 * 4), nrow = 36, ncol = 4)
  x[19:36, ] <- x[19:36, ] + 1
  incomplete <- x
  incomplete[seq(1, 36, by = 4), 1] <- NA_real_
  incomplete[seq(2, 36, by = 5), 3] <- NA_real_
  fit <- spatmix(
    incomplete, G = 2, coords = 1:4, spatial_decay = "sigmoid",
    max_iter = 6, verbose = FALSE
  )

  expect_true(fit$has_missing)
  expect_identical(dim(fit$imputation), dim(x))
  expect_false(anyNA(fit$imputation))
  expect_equal(fit$imputation[!is.na(incomplete)], incomplete[!is.na(incomplete)])
})

test_that("sigmoid decay supports matrix-variate data and grid axes", {
  set.seed(12)
  x <- array(rnorm(4 * 3 * 30), c(4, 3, 30))
  x[, , 16:30] <- x[, , 16:30] + 1.2
  fit <- spatmix(
    x, G = 2, r = 1, coords = list(x = 1:2, y = 1:2),
    spatial_decay = "sigmoid", max_iter = 6, verbose = FALSE
  )

  expect_false(fit$settings$spatial_only)
  expect_identical(fit$settings$spatial_decay, "sigmoid")
  expect_identical(fit$settings$decay_group, c(1L, 1L))
  expect_length(fit$coordinate_beta, 1L)
  expect_identical(dim(fit$alpha), c(3L, 2L))
  expect_identical(dim(fit$sigmoid), c(1L, 2L))
  expect_true(all(fit$sigmoid > 0))
})

test_that("sigmoid decay can be made axis-specific", {
  set.seed(13)
  x <- array(rnorm(4 * 3 * 30), c(4, 3, 30))
  x[, , 16:30] <- x[, , 16:30] + 1.2
  fit <- spatmix(
    x, G = 2, r = 1, coords = list(x = 1:2, y = 1:2),
    spatial_decay = "sigmoid", decay_group = NULL,
    sigmoid_ctrl = list(init = 1, lower = 0.001, upper = 10),
    max_iter = 6, verbose = FALSE
  )

  expect_identical(fit$settings$decay_group, c(1L, 2L))
  expect_identical(fit$settings$decay_group_labels, c("1", "2"))
  expect_length(fit$coordinate_beta, 2L)
  expect_length(fit$sigmoid, 2L)
  expect_true(all(vapply(fit$sigmoid, function(x) all(x > 0), logical(1L))))
  expect_identical(fit$n_parameters, 49)
})

test_that("sigmoid decay defaults to sharing across grid axes", {
  set.seed(13)
  x <- array(rnorm(4 * 3 * 30), c(4, 3, 30))
  x[, , 16:30] <- x[, , 16:30] + 1.2
  fit <- spatmix(
    x, G = 2, r = 1, coords = list(x = 1:2, y = 1:2),
    spatial_decay = "sigmoid",
    sigmoid_ctrl = list(init = 1, lower = 0.001, upper = 10),
    max_iter = 6, verbose = FALSE
  )

  expect_identical(fit$settings$decay_group, c(1L, 1L))
  expect_identical(fit$settings$decay_group_labels, "1")
  expect_length(fit$coordinate_beta, 1L)
  expect_identical(dim(fit$alpha), c(3L, 2L))
  expect_identical(dim(fit$sigmoid), c(1L, 2L))
  expect_identical(fit$n_parameters, 45)
})

test_that("I-spline decay supports partial sharing across grid axes", {
  set.seed(14)
  x <- array(rnorm(8 * 2 * 32), c(8, 2, 32))
  x[, , 17:32] <- x[, , 17:32] + 1
  fit <- spatmix(
    x, G = 2, r = 1, coords = list(x = 1:2, y = 1:2, z = 1:2),
    nknots = 0, degree = 1, decay_group = c("xy", "xy", "z"),
    max_iter = 6, spatial_max_iter = 20, verbose = FALSE
  )

  expect_identical(fit$settings$decay_group, c(1L, 1L, 2L))
  expect_identical(fit$settings$decay_group_labels, c("xy", "z"))
  expect_length(fit$coordinate_beta, 2L)
  expect_identical(dim(fit$alpha), c(4L, 2L))
  expect_identical(fit$n_parameters, 53)
  expect_true(all(vapply(seq_len(2), function(g) {
    min(eigen(fit$row_covariance[, , g], symmetric = TRUE,
              only.values = TRUE)$values) > 0
  }, logical(1L))))
})

test_that("shared I-spline axes require compatible bases", {
  set.seed(15)
  x <- array(rnorm(6 * 2 * 20), c(6, 2, 20))
  expect_error(
    spatmix(
      x, G = 2, r = 1, coords = list(x = 1:2, y = 1:3),
      nknots = c(0, 1), degree = 1, decay_group = 1,
      max_iter = 1, verbose = FALSE
    ),
    "same number of basis functions"
  )
})
