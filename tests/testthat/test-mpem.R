test_that("MPEM remains an internal SpatMix engine", {
  expect_setequal(
    getNamespaceExports("SpatMix"),
    "spatmix"
  )
  expect_false(exists("introduce_missing", envir = asNamespace("SpatMix"),
                      inherits = FALSE))
  expect_false(exists("impute", envir = asNamespace("SpatMix"), inherits = FALSE))
  expect_false(exists("mpem", envir = asNamespace("SpatMix"), inherits = FALSE))
  expect_false(exists("mpem_cpp", envir = asNamespace("SpatMix"), inherits = FALSE))
  expect_true(is.function(getFromNamespace("mpem_moments", "SpatMix")))
  expect_true(is.function(getFromNamespace("mpem_moments_cpp", "SpatMix")))
})
