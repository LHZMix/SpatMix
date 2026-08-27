// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

struct MpemMoments {
  arma::mat completed;
  std::vector<arma::mat> conditional_cov;
  arma::mat row_correction;
  arma::mat col_correction;
  double scale_correction;
};

static arma::uvec finite_indices(const arma::rowvec& x, bool finite) {
  std::vector<arma::uword> index;
  index.reserve(x.n_elem);
  for (arma::uword j = 0; j < x.n_elem; ++j) {
    if (std::isfinite(x[j]) == finite) index.push_back(j);
  }
  return arma::uvec(index);
}

static bool valid_covariance_block(const arma::mat& x,
                                   arma::uword dimension,
                                   double eps) {
  if (x.n_rows != dimension || x.n_cols != dimension || !x.is_finite())
    return false;
  for (arma::uword j = 0; j < dimension; ++j) {
    if (x(j, j) <= eps) return false;
  }
  return true;
}

static std::vector<arma::uvec> missing_indices(const arma::mat& x) {
  std::vector<arma::uvec> out(x.n_rows);
  for (arma::uword i = 0; i < x.n_rows; ++i)
    out[i] = finite_indices(x.row(i), false);
  return out;
}

static MpemMoments mpem_moments_core(
    const arma::mat& x, int rows, int cols, const arma::mat& mean,
    const arma::mat& row_precision, const arma::mat& col_precision,
    double precision_multiplier, arma::mat completed,
    std::vector<arma::mat> conditional_cov,
    const std::vector<arma::uvec>& missing_index,
    const arma::vec& weights, int sweeps, double warm_start,
    bool update_moments) {
  const double eps = 1e-10;
  const arma::vec row_diagonal = row_precision.diag();
  const arma::vec col_diagonal = col_precision.diag();
  const arma::uword observations = x.n_rows;
  warm_start = std::max(0.0, std::min(1.0, warm_start));

  MpemMoments out;
  out.completed = completed;
  out.conditional_cov = conditional_cov;
  out.row_correction.zeros(rows, rows);
  out.col_correction.zeros(cols, cols);
  out.scale_correction = 0.0;

  for (arma::uword i = 0; i < observations; ++i) {
    const arma::uvec& missing = missing_index[i];
    const arma::uword number_missing = missing.n_elem;
    if (number_missing == 0) continue;

    arma::uvec missing_rows(number_missing);
    arma::uvec missing_cols(number_missing);
    for (arma::uword j = 0; j < number_missing; ++j) {
      missing_rows[j] = missing[j] % rows;
      missing_cols[j] = missing[j] / rows;
    }

    if (update_moments) {
      arma::vec y = out.completed.row(i).t();
      arma::mat observation(y.memptr(), rows, cols, false, true);
      arma::mat residual_precision =
        row_precision * (observation - mean) * col_precision;

      for (int sweep = 0; sweep < sweeps; ++sweep) {
        for (arma::uword j = 0; j < number_missing; ++j) {
          arma::uword row = missing_rows[j];
          arma::uword col = missing_cols[j];
          double denominator = row_diagonal[row] * col_diagonal[col];
          if (!std::isfinite(denominator) || std::abs(denominator) < eps)
            denominator = eps;
          double delta = -residual_precision(row, col) / denominator;
          y[missing[j]] += delta;
          residual_precision += delta *
            row_precision.col(row) * col_precision.row(col);
        }
      }
      out.completed.row(i) = y.t();

      arma::vec diagonal(number_missing);
      for (arma::uword j = 0; j < number_missing; ++j) {
        double denominator = precision_multiplier *
          row_diagonal[missing_rows[j]] * col_diagonal[missing_cols[j]];
        if (!std::isfinite(denominator) || std::abs(denominator) < eps)
          denominator = eps;
        diagonal[j] = 1.0 / denominator;
      }
      arma::mat diagonal_start = arma::diagmat(diagonal);
      arma::mat covariance = out.conditional_cov[i];
      if (valid_covariance_block(covariance, number_missing, eps)) {
        covariance = warm_start * 0.5 * (covariance + covariance.t()) +
          (1.0 - warm_start) * diagonal_start;
      } else {
        covariance = diagonal_start;
      }

      arma::rowvec coefficient(number_missing, arma::fill::zeros);
      for (int sweep = 0; sweep < sweeps; ++sweep) {
        for (arma::uword aa = 0; aa < number_missing; ++aa) {
          arma::uword row_a = missing_rows[aa];
          arma::uword col_a = missing_cols[aa];
          double denominator = row_diagonal[row_a] * col_diagonal[col_a];
          if (!std::isfinite(denominator) || std::abs(denominator) < eps)
            denominator = eps;
          if (number_missing == 1) {
            covariance(0, 0) = 1.0 / (precision_multiplier * denominator);
            continue;
          }

          coefficient.zeros();
          for (arma::uword jj = 0; jj < number_missing; ++jj) {
            if (jj != aa) {
              coefficient[jj] =
                row_precision(row_a, missing_rows[jj]) *
                col_precision(col_a, missing_cols[jj]) / denominator;
            }
          }
          arma::rowvec off_diagonal = -coefficient * covariance;
          double updated_diagonal = 1.0 / (precision_multiplier * denominator) -
            arma::dot(coefficient, off_diagonal);
          for (arma::uword jj = 0; jj < number_missing; ++jj) {
            if (jj != aa) {
              covariance(aa, jj) = off_diagonal[jj];
              covariance(jj, aa) = off_diagonal[jj];
            }
          }
          covariance(aa, aa) =
            (std::isfinite(updated_diagonal) && updated_diagonal > eps)
              ? updated_diagonal
              : 1.0 / (precision_multiplier * denominator);
        }
      }
      out.conditional_cov[i] = 0.5 * (covariance + covariance.t());
    }

    arma::mat covariance = out.conditional_cov[i];
    if (!valid_covariance_block(covariance, number_missing, eps)) {
      arma::vec diagonal(number_missing);
      for (arma::uword j = 0; j < number_missing; ++j) {
        double denominator = precision_multiplier *
          row_diagonal[missing_rows[j]] * col_diagonal[missing_cols[j]];
        if (!std::isfinite(denominator) || std::abs(denominator) < eps)
          denominator = eps;
        diagonal[j] = 1.0 / denominator;
      }
      covariance = arma::diagmat(diagonal);
      out.conditional_cov[i] = covariance;
    }

    for (arma::uword jj = 0; jj < number_missing; ++jj) {
      arma::uword col_j = missing_cols[jj];
      arma::uword row_j = missing_rows[jj];
      for (arma::uword kk = 0; kk < number_missing; ++kk) {
        arma::uword col_k = missing_cols[kk];
        arma::uword row_k = missing_rows[kk];
        double value = weights[i] * covariance(kk, jj);
        out.col_correction(col_k, col_j) +=
          value * row_precision(row_k, row_j);
        out.row_correction(row_k, row_j) +=
          value * col_precision(col_k, col_j);
        out.scale_correction += value *
          row_precision(row_k, row_j) * col_precision(col_k, col_j);
      }
    }
  }
  return out;
}

static Rcpp::List covariance_list(const std::vector<arma::mat>& covariance) {
  Rcpp::List out(covariance.size());
  for (std::size_t i = 0; i < covariance.size(); ++i) out[i] = covariance[i];
  return out;
}

// [[Rcpp::export]]
Rcpp::List mpem_moments_cpp(
    const arma::mat& x, int rows, int cols, const arma::mat& mean,
    const arma::mat& row_precision, const arma::mat& col_precision,
    double precision_multiplier, const arma::mat& completed,
    const Rcpp::List& conditional_cov, const arma::vec& weights,
    int sweeps = 1, double warm_start = 1.0, bool update_moments = true) {
  const arma::uword observations = x.n_rows;
  const arma::uword dimension = rows * cols;
  if (x.n_cols != dimension || completed.n_rows != observations ||
      completed.n_cols != dimension)
    Rcpp::stop("x and completed must have compatible dimensions.");
  if (mean.n_rows != static_cast<arma::uword>(rows) ||
      mean.n_cols != static_cast<arma::uword>(cols))
    Rcpp::stop("mean has incompatible dimensions.");
  if (row_precision.n_rows != static_cast<arma::uword>(rows) ||
      row_precision.n_cols != static_cast<arma::uword>(rows) ||
      col_precision.n_rows != static_cast<arma::uword>(cols) ||
      col_precision.n_cols != static_cast<arma::uword>(cols))
    Rcpp::stop("precision matrices have incompatible dimensions.");
  if (weights.n_elem != observations)
    Rcpp::stop("weights must have one value per observation.");
  if (conditional_cov.size() != static_cast<R_xlen_t>(observations))
    Rcpp::stop("conditional_cov must have one element per observation.");
  if (update_moments && sweeps < 1)
    Rcpp::stop("sweeps must be positive when moments are updated.");
  if (!std::isfinite(precision_multiplier) || precision_multiplier <= 0)
    Rcpp::stop("precision_multiplier must be positive.");

  std::vector<arma::mat> covariance(observations);
  for (arma::uword i = 0; i < observations; ++i) {
    SEXP element = conditional_cov[i];
    if (!Rf_isNull(element)) covariance[i] = Rcpp::as<arma::mat>(element);
  }
  std::vector<arma::uvec> missing = missing_indices(x);
  MpemMoments moments = mpem_moments_core(
    x, rows, cols, mean, row_precision, col_precision,
    precision_multiplier, completed, covariance, missing, weights,
    sweeps, warm_start, update_moments
  );
  return Rcpp::List::create(
    Rcpp::Named("completed") = moments.completed,
    Rcpp::Named("conditional_cov") = covariance_list(moments.conditional_cov),
    Rcpp::Named("row_correction") = moments.row_correction,
    Rcpp::Named("col_correction") = moments.col_correction,
    Rcpp::Named("scale_correction") = moments.scale_correction
  );
}
