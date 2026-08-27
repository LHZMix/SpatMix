## Test environments

- Local: macOS Sequoia 15.7.7, R 4.5.2
- GitHub Actions: Windows, macOS, and Ubuntu with R-release
- GitHub Actions: Ubuntu with R-devel and R-oldrel-1

## R CMD check results

- `R CMD check --no-manual`: 0 errors | 0 warnings | 0 notes
- `R CMD check --as-cran`: 0 errors | 0 warnings | 3 notes
- GitHub Actions: 0 errors | 0 warnings on all five configurations

The `--as-cran` notes identify this as a new submission, skip the README check
because Pandoc is not installed locally, and skip HTML validation because the
local HTML Tidy is not recent enough. The PDF manual, examples, and tests pass.

## Submission notes

This is a new submission.

`SpatMix` fits spatial Gaussian mixtures and mixtures with non-spatial
factor-analyzer covariance structures. Spatial decay may be represented by
monotone I-splines or a normalized sigmoid function. Missing values are
represented by `NA` and handled by internal matrix partial EM updates.

The examples demonstrate a spatial-only sigmoid fit and a matrix-variate grid
fit in which data generated with a shared sigmoid decay across axes are fitted
using a shared I-spline decay. Missing-data fitting is exercised in the test
suite.
