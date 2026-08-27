#' Conditional moments for MPEM updates
#'
#' Internal bridge to the compiled conditional-moment core used by `spatmix()`.
#'
#' @noRd
mpem_moments <- function(x, rows, cols, mean, row_precision, col_precision,
                         precision_multiplier, completed, conditional_cov,
                         weights, sweeps = 1L, warm_start = 1,
                         update_moments = TRUE) {
  mpem_moments_cpp(
    x = x,
    rows = rows,
    cols = cols,
    mean = mean,
    row_precision = row_precision,
    col_precision = col_precision,
    precision_multiplier = precision_multiplier,
    completed = completed,
    conditional_cov = conditional_cov,
    weights = weights,
    sweeps = sweeps,
    warm_start = warm_start,
    update_moments = update_moments
  )
}
