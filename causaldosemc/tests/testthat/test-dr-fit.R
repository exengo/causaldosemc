test_that("cross-fitted DR estimator returns linear lag coefficients", {
  panel <- simulate_cdmc_data(
    n_units = 20,
    n_times = 12,
    rank = 2,
    beta = 1,
    lag_beta = 0.3,
    n_covariates = 1,
    noise_sd = 0.03,
    switch_on_prob = 0.15,
    switch_off_prob = 0.4,
    seed = 1111
  )
  panel$w <- rep(seq(0.9, 1.1, length.out = 12), times = 20)

  fit <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weights = "w",
    n_folds = 2,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 1,
    seed = 1111
  )

  expect_s3_class(fit, "cdmc_dr_fit")
  expect_equal(fit$n_folds, 2)
  expect_true(all(c("dose_lag0", "dose_lag1") %in% names(fit$effect$coefficients)))
  expect_equal(unname(fit$effect$coefficients["dose_lag0"]), 1, tolerance = 0.8)
  expect_equal(unname(fit$effect$coefficients["dose_lag1"]), 0.3, tolerance = 0.8)
  expect_true(sum(fit$effect$sample_mask, na.rm = TRUE) > 0)
})

test_that("cross-fitted DR estimator reuses supplied fold assignments deterministically", {
  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 8,
    rank = 2,
    beta = 0.9,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.18,
    switch_off_prob = 0.42,
    seed = 1122
  )

  fit_seeded <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weight_method = "gaussian_gps",
    n_folds = 2,
    lambda = 0.2,
    rank_max = 2,
    washout = 0,
    lag_order = 1,
    seed = 1122
  )

  fit_reused <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weight_method = "gaussian_gps",
    fold_assignments = fit_seeded$fold_assignments,
    lambda = 0.2,
    rank_max = 2,
    washout = 0,
    lag_order = 1,
    seed = 9999
  )

  expect_s3_class(fit_reused, "cdmc_dr_fit")
  expect_equal(fit_reused$fold_assignments, fit_seeded$fold_assignments)
  expect_equal(fit_reused$data$.cdmc_fold_id, fit_seeded$data$.cdmc_fold_id)
  expect_equal(fit_reused$effect$coefficients, fit_seeded$effect$coefficients, tolerance = 1e-10)
  expect_equal(fit_reused$effect$tau_dr, fit_seeded$effect$tau_dr, tolerance = 1e-10)
})

test_that("cross-fitted DR estimator accepts CBPS weight objects", {
  testthat::skip_if_not_installed("CBPS")

  panel <- simulate_cdmc_data(
    n_units = 20,
    n_times = 10,
    rank = 2,
    beta = 0.8,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.18,
    switch_off_prob = 0.42,
    seed = 1212
  )
  cbps_weights <- cdmc_cbps_weights(
    data = panel,
    dose = "dose",
    covariates = "x1",
    iterations = 250L
  )

  fit <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weights = cbps_weights,
    n_folds = 2,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 1,
    seed = 1212
  )

  expect_s3_class(fit, "cdmc_dr_fit")
  expect_true(all(is.finite(fit$data$.cdmc_weight)))
  expect_true(any(fit$data$.cdmc_dr_sample))
})

test_that("cross-fitted DR estimator can estimate internal CBPS weights", {
  testthat::skip_if_not_installed("CBPS")

  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 8,
    rank = 2,
    beta = 0.85,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.18,
    switch_off_prob = 0.42,
    seed = 1262
  )

  fit <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weight_method = "cbps",
    cbps_iterations = 250L,
    n_folds = 2,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 1,
    seed = 1262
  )

  expect_s3_class(fit, "cdmc_dr_fit")
  expect_identical(fit$weight_method, "cbps")
  expect_identical(fit$cbps_method, "over")
  expect_true(all(is.finite(fit$data$.cdmc_weight[fit$data$.cdmc_observed])))
  expect_true(any(abs(fit$data$.cdmc_weight[fit$data$.cdmc_observed] - 1) > 1e-6))
  expect_true(any(fit$data$.cdmc_dr_sample))
  expect_true(all(vapply(
    fit$baseline$fold_summaries,
    function(summary) !is.null(summary$cbps_train_method) && !is.null(summary$cbps_holdout_method),
    logical(1)
  )))
})

test_that("cross-fitted DR estimator can estimate internal gaussian GPS weights", {
  panel <- simulate_cdmc_data(
    n_units = 20,
    n_times = 10,
    rank = 2,
    beta = 0.9,
    lag_beta = 0.25,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.18,
    switch_off_prob = 0.42,
    seed = 1313
  )

  fit <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weight_method = "gaussian_gps",
    n_folds = 2,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 1,
    seed = 1313
  )

  expect_s3_class(fit, "cdmc_dr_fit")
  expect_equal(fit$weight_method, "gaussian_gps")
  expect_true(all(is.finite(fit$data$.cdmc_weight)))
  expect_true(any(abs(fit$data$.cdmc_weight - 1) > 1e-6))
  expect_true(all(c("dose_lag0", "dose_lag1") %in% names(fit$effect$coefficients)))
})

test_that("internal gaussian GPS weighting works without covariates", {
  panel <- simulate_cdmc_data(
    n_units = 16,
    n_times = 10,
    rank = 2,
    beta = 0.7,
    lag_beta = 0.15,
    n_covariates = 0,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 1414
  )

  fit <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    weight_method = "gaussian_gps",
    n_folds = 2,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 1,
    seed = 1414
  )

  expect_s3_class(fit, "cdmc_dr_fit")
  expect_true(all(is.finite(fit$data$.cdmc_weight)))
  expect_true(any(fit$data$.cdmc_dr_sample))
})

test_that("internal spline Gaussian GPS weights support nonlinear assignment patterns", {
  panel <- simulate_cdmc_data(
    n_units = 20,
    n_times = 10,
    rank = 2,
    beta = 0.6,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.18,
    switch_off_prob = 0.42,
    seed = 1515
  )
  original_dose <- panel$dose
  active <- abs(original_dose) > 0
  panel$dose[active] <- sign(original_dose[active]) * (0.25 + panel$x1[active]^2 + panel$time[active] / max(panel$time))
  panel$y <- panel$y + 0.6 * (panel$dose - original_dose)

  fit_linear <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weight_method = "gaussian_gps",
    gps_model = "linear",
    n_folds = 2,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 0,
    seed = 1515
  )

  fit_spline <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weight_method = "gaussian_gps",
    gps_model = "spline",
    gps_df = 3,
    n_folds = 2,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 0,
    seed = 1515
  )

  spline_formulas <- vapply(
    fit_spline$baseline$fold_summaries,
    function(summary) paste(deparse(summary$gps_formula), collapse = " "),
    character(1)
  )

  expect_s3_class(fit_spline, "cdmc_dr_fit")
  expect_identical(fit_spline$gps_model, "spline")
  expect_identical(fit_spline$gps_df, 3L)
  expect_true(all(is.finite(fit_spline$data$.cdmc_weight)))
  expect_true(any(grepl("splines::ns\\(", spline_formulas)))
  expect_true(any(abs(fit_spline$data$.cdmc_weight - fit_linear$data$.cdmc_weight) > 1e-6, na.rm = TRUE))
  expect_true(any(fit_spline$data$.cdmc_dr_sample))
})

test_that("internal GAM Gaussian GPS weights support smooth nonlinear assignment patterns", {
  testthat::skip_if_not_installed("mgcv")

  panel <- simulate_cdmc_data(
    n_units = 20,
    n_times = 10,
    rank = 2,
    beta = 0.65,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.18,
    switch_off_prob = 0.42,
    seed = 1545
  )
  original_dose <- panel$dose
  active <- abs(original_dose) > 0
  panel$dose[active] <- sign(original_dose[active]) * (
    0.2 + sin(2 * panel$x1[active]) + (panel$time[active] / max(panel$time))^2
  )
  panel$y <- panel$y + 0.5 * (panel$dose - original_dose)

  fit_linear <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weight_method = "gaussian_gps",
    gps_model = "linear",
    n_folds = 2,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 0,
    seed = 1545
  )

  fit_gam <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weight_method = "gaussian_gps",
    gps_model = "gam",
    gps_df = 4,
    n_folds = 2,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 0,
    seed = 1545
  )

  gam_formulas <- vapply(
    fit_gam$baseline$fold_summaries,
    function(summary) paste(deparse(summary$gps_formula), collapse = " "),
    character(1)
  )

  expect_s3_class(fit_gam, "cdmc_dr_fit")
  expect_identical(fit_gam$gps_model, "gam")
  expect_true(all(is.finite(fit_gam$data$.cdmc_weight)))
  expect_true(any(grepl("s\\(", gam_formulas)))
  expect_true(any(abs(fit_gam$data$.cdmc_weight - fit_linear$data$.cdmc_weight) > 1e-6, na.rm = TRUE))
  expect_true(any(fit_gam$data$.cdmc_dr_sample))
})

test_that("internal tree Gaussian GPS weights support threshold assignment patterns", {
  testthat::skip_if_not_installed("rpart")

  panel <- simulate_cdmc_data(
    n_units = 20,
    n_times = 10,
    rank = 2,
    beta = 0.6,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.18,
    switch_off_prob = 0.42,
    seed = 1555
  )
  original_dose <- panel$dose
  active <- abs(original_dose) > 0
  x_threshold <- stats::median(panel$x1)
  t_threshold <- stats::median(panel$time)
  panel$dose[active] <- sign(original_dose[active]) * (
    0.2 +
      ifelse(panel$x1[active] > x_threshold, 1.0, 0.2) +
      ifelse(panel$time[active] > t_threshold, 0.7, 0.1)
  )
  panel$y <- panel$y + 0.5 * (panel$dose - original_dose)

  fit_linear <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weight_method = "gaussian_gps",
    gps_model = "linear",
    n_folds = 2,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 0,
    seed = 1555
  )

  fit_tree <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weight_method = "gaussian_gps",
    gps_model = "tree",
    n_folds = 2,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 0,
    seed = 1555
  )

  expect_s3_class(fit_tree, "cdmc_dr_fit")
  expect_identical(fit_tree$gps_model, "tree")
  expect_true(all(is.finite(fit_tree$data$.cdmc_weight)))
  expect_true(any(abs(fit_tree$data$.cdmc_weight - fit_linear$data$.cdmc_weight) > 1e-6, na.rm = TRUE))
  expect_true(any(fit_tree$data$.cdmc_dr_sample))
})

test_that("internal forest Gaussian GPS weights support interaction-heavy assignment patterns", {
  testthat::skip_if_not_installed("ranger")

  panel <- simulate_cdmc_data(
    n_units = 20,
    n_times = 10,
    rank = 2,
    beta = 0.6,
    lag_beta = NULL,
    n_covariates = 2,
    noise_sd = 0.04,
    switch_on_prob = 0.18,
    switch_off_prob = 0.42,
    seed = 1561
  )
  original_dose <- panel$dose
  active <- abs(original_dose) > 0
  x1_threshold <- stats::median(panel$x1)
  interaction_score <- panel$x2 + panel$time / max(panel$time)
  panel$dose[active] <- sign(original_dose[active]) * (
    0.15 +
      ifelse(panel$x1[active] > x1_threshold, 0.9, 0.2) +
      ifelse(interaction_score[active] > stats::median(interaction_score), 0.7, -0.1)
  )
  panel$y <- panel$y + 0.5 * (panel$dose - original_dose)

  fit_linear <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = c("x1", "x2"),
    weight_method = "gaussian_gps",
    gps_model = "linear",
    n_folds = 2,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 0,
    seed = 1561
  )

  fit_forest <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = c("x1", "x2"),
    weight_method = "gaussian_gps",
    gps_model = "forest",
    gps_forest_trees = 50L,
    gps_forest_mtry = 2L,
    gps_forest_min_node_size = 3L,
    n_folds = 2,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 0,
    seed = 1561
  )

  expect_s3_class(fit_forest, "cdmc_dr_fit")
  expect_identical(fit_forest$gps_model, "forest")
  expect_identical(fit_forest$gps_forest_trees, 50L)
  expect_identical(fit_forest$gps_forest_mtry, 2L)
  expect_identical(fit_forest$gps_forest_min_node_size, 3L)
  expect_true(all(is.finite(fit_forest$data$.cdmc_weight)))
  expect_true(any(abs(fit_forest$data$.cdmc_weight - fit_linear$data$.cdmc_weight) > 1e-6, na.rm = TRUE))
  expect_true(all(vapply(
    fit_forest$baseline$fold_summaries,
    function(summary) identical(summary$gps_forest_trees, 50L),
    logical(1)
  )))
  expect_true(any(fit_forest$data$.cdmc_dr_sample))
})

test_that("internal kernel GPS weights support non-Gaussian assignment patterns", {
  panel <- simulate_cdmc_data(
    n_units = 20,
    n_times = 10,
    rank = 2,
    beta = 0.7,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.18,
    switch_off_prob = 0.42,
    seed = 1565
  )
  original_dose <- panel$dose
  active <- abs(original_dose) > 0
  panel$dose[active] <- sign(original_dose[active]) * (
    abs(original_dose[active]) + 0.25 + panel$x1[active]^2 + rexp(sum(active), rate = 2)
  )
  panel$y <- panel$y + 0.7 * (panel$dose - original_dose)

  fit_gaussian <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weight_method = "gaussian_gps",
    gps_model = "spline",
    gps_df = 3,
    n_folds = 2,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 0,
    seed = 1565
  )

  fit_kernel <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weight_method = "kernel_gps",
    gps_model = "spline",
    gps_df = 3,
    gps_bandwidth = 0.3,
    n_folds = 2,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 0,
    seed = 1565
  )

  expect_s3_class(fit_kernel, "cdmc_dr_fit")
  expect_identical(fit_kernel$weight_method, "kernel_gps")
  expect_equal(fit_kernel$gps_bandwidth, 0.3)
  expect_true(all(is.finite(fit_kernel$data$.cdmc_weight)))
  expect_true(any(abs(fit_kernel$data$.cdmc_weight - fit_gaussian$data$.cdmc_weight) > 1e-6, na.rm = TRUE))
  expect_true(any(fit_kernel$data$.cdmc_dr_sample))
})

test_that("cross-fitted DR estimator supports unbalanced panels", {
  panel <- simulate_cdmc_data(
    n_units = 18,
    n_times = 10,
    rank = 2,
    beta = 0.8,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.18,
    switch_off_prob = 0.42,
    seed = 1616
  )
  drop_index <- head(which(panel$dose != 0), 3L)
  panel <- panel[-drop_index, ]

  fit <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weight_method = "gaussian_gps",
    n_folds = 2,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 1,
    seed = 1616
  )

  expect_s3_class(fit, "cdmc_dr_fit")
  expect_true(any(fit$data$.cdmc_dr_sample))
  expect_true(all(is.na(fit$data$.cdmc_tau_dr_linear[!fit$data$.cdmc_observed])))
  expect_true(all(fit$data$.cdmc_weight[!fit$data$.cdmc_observed] == 0))
})