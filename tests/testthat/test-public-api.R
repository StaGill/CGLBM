test_that("the exported surface is the Stage-4 API", {
  expected <- sort(c(
    "apply_detection_limits", "center_cglbm_data", "censored_lbm",
    "cglbm_icl_bic", "comparison_cari", "complete_boundary",
    "complete_halfdouble", "complete_mi", "complete_truncnorm",
    "exact_cglbm_data", "new_cglbm_data", "pool_mi_fits",
    "select_cglbm_model"
  ))
  expect_identical(sort(getNamespaceExports("CGLBM")), expected)
  expect_false(exists("fit_cglbm_hh", asNamespace("CGLBM"), inherits = FALSE))
  expect_false(exists("fit_cglbm_hs", asNamespace("CGLBM"), inherits = FALSE))
  expect_false(exists("fit_cglbm_sh", asNamespace("CGLBM"), inherits = FALSE))
  expect_false(exists("fit_cglbm_ss", asNamespace("CGLBM"), inherits = FALSE))
  expect_false(exists("%||%", asNamespace("CGLBM"), inherits = FALSE))
})

test_that("S3 methods are registered", {
  expect_true(is.function(getS3method("print", "cglbm_data")))
  expect_true(is.function(getS3method("print", "cglbm_fit")))
  expect_true(is.function(getS3method("predict", "cglbm_fit")))
  expect_true(is.function(getS3method("summary", "cglbm_fit")))
  expect_true(is.function(getS3method("print", "cglbm_model_selection")))
})

test_that("compiled routines are registered", {
  dll <- getLoadedDLLs()[["CGLBM"]]
  expect_false(is.null(dll))
  expect_true(is.loaded("_CGLBM_cglbm_cpp_all_cell_loglik_raw",
                        PACKAGE = "CGLBM"))
  expect_true(is.loaded("_CGLBM_cglbm_cpp_weighted_gaussian_em_step_raw",
                        PACKAGE = "CGLBM"))
})
