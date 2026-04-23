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

test_that("cross-fitted DR estimator validates cv_workers", {
  panel <- simulate_cdmc_data(
    n_units = 8,
    n_times = 8,
    rank = 2,
    beta = 0.8,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.18,
    switch_off_prob = 0.42,
    seed = 1121
  )

  expect_error(
    cdmc_dr_fit(
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
      cv_workers = 0,
      seed = 1121
    ),
    "cv_workers must be a positive integer"
  )
})

test_that("cross-fitted DR estimator validates dr_workers", {
  panel <- simulate_cdmc_data(
    n_units = 8,
    n_times = 8,
    rank = 2,
    beta = 0.8,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.18,
    switch_off_prob = 0.42,
    seed = 1123
  )

  expect_error(
    cdmc_dr_fit(
      data = panel,
      outcome = "y",
      dose = "dose",
      unit = "unit",
      time = "time",
      covariates = "x1",
      weight_method = "gaussian_gps",
      n_folds = 2,
      dr_workers = 0,
      lambda = 0.2,
      rank_max = 2,
      washout = 0,
      lag_order = 1,
      seed = 1123
    ),
    "dr_workers must be a positive integer"
  )
})

test_that("cross-fitted DR estimator parallel folds match sequential folds", {
  skip_on_os("windows")

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
    seed = 1125
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
    seed = 1125
  )

  fit_seq <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weight_method = "gaussian_gps",
    fold_assignments = fit_seeded$fold_assignments,
    dr_workers = 1,
    lambda = 0.2,
    rank_max = 2,
    washout = 0,
    lag_order = 1,
    seed = 999
  )

  fit_par <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weight_method = "gaussian_gps",
    fold_assignments = fit_seeded$fold_assignments,
    dr_workers = 2,
    lambda = 0.2,
    rank_max = 2,
    washout = 0,
    lag_order = 1,
    seed = 999
  )

  expect_false(isTRUE(fit_seq$dr_parallel))
  expect_true(isTRUE(fit_par$dr_parallel))
  expect_equal(fit_seq$dr_workers, 1L)
  expect_equal(fit_par$dr_workers, 2L)
  expect_equal(fit_seq$effect$coefficients, fit_par$effect$coefficients, tolerance = 1e-10)
  expect_equal(fit_seq$effect$tau_dr, fit_par$effect$tau_dr, tolerance = 1e-10)
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

test_that("adaptive internal weight clipping caps extreme weights deterministically", {
  raw_weights <- c(rep(1, 99), 100)

  cap <- cdmc_resolve_internal_max_weight(raw_weights, max_weight = "adaptive")
  clipped <- cdmc_cap_internal_weights(raw_weights, max_weight = "adaptive")

  expect_equal(cap, 10, tolerance = 1e-8)
  expect_equal(max(clipped), cap, tolerance = 1e-8)
  expect_equal(cdmc_cap_internal_weights(raw_weights, max_weight = NULL), raw_weights)
  expect_equal(cdmc_cap_internal_weights(raw_weights, max_weight = 5), pmin(raw_weights, 5))
})

test_that("cross-fitted DR estimator uses adaptive internal clipping by default", {
  panel <- simulate_cdmc_data(
    n_units = 16,
    n_times = 10,
    rank = 2,
    beta = 0.8,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.18,
    switch_off_prob = 0.42,
    seed = 1148
  )
  original_dose <- panel$dose
  active <- abs(original_dose) > 0
  panel$dose[active] <- sign(original_dose[active]) * (
    0.15 + 2.5 * panel$x1[active]^3 + ifelse(panel$time[active] > stats::median(panel$time), 1.2, -0.2)
  )
  panel$y <- panel$y + 0.4 * (panel$dose - original_dose)

  fit_adaptive <- cdmc_dr_fit(
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
    seed = 1148
  )

  fit_none <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weight_method = "gaussian_gps",
    gps_model = "linear",
    max_weight = NULL,
    fold_assignments = fit_adaptive$fold_assignments,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 0,
    seed = 1148
  )

  applied_caps <- vapply(
    fit_adaptive$baseline$fold_summaries,
    function(summary) summary$applied_max_weight %||% NA_real_,
    numeric(1)
  )
  applied_caps <- applied_caps[is.finite(applied_caps)]

  expect_s3_class(fit_adaptive, "cdmc_dr_fit")
  expect_identical(fit_adaptive$max_weight, "adaptive")
  expect_null(fit_none$max_weight)
  expect_true(length(applied_caps) > 0L)
  expect_true(all(is.finite(applied_caps)))
  expect_equal(nrow(fit_adaptive$weight_diagnostics$by_fold), fit_adaptive$n_folds)
  expect_true(all(c(
    "observed_ess",
    "observed_clipped_share",
    "dr_sample_ess",
    "dr_sample_ess_fraction"
  ) %in% names(fit_adaptive$weight_diagnostics$by_fold)))
  expect_true(is.finite(fit_adaptive$weight_diagnostics$overall$observed$ess))
  expect_true(is.finite(fit_adaptive$weight_diagnostics$overall$dr_sample$ess))
  expect_true(fit_adaptive$weight_diagnostics$overall$observed$ess <= fit_adaptive$weight_diagnostics$overall$observed$n_weights + 1e-8)
  expect_true(fit_adaptive$weight_diagnostics$overall$dr_sample$ess <= fit_adaptive$weight_diagnostics$overall$dr_sample$n_weights + 1e-8)
  expect_true(max(fit_adaptive$data$.cdmc_weight, na.rm = TRUE) <= max(applied_caps) + 1e-8)
  expect_equal(fit_none$weight_diagnostics$overall$observed$clipped_count, 0L)
  expect_true(max(fit_none$data$.cdmc_weight, na.rm = TRUE) >= max(fit_adaptive$data$.cdmc_weight, na.rm = TRUE))
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
    function(summary) !is.null(summary$cbps_train_method) && !is.null(summary$weight_transport_model),
    logical(1)
  )))
})

test_that("cross-fitted DR estimator can estimate internal entropy-balance weights", {
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
    seed = 1268
  )
  original_dose <- panel$dose
  active <- abs(original_dose) > 0
  panel$dose[active] <- sign(original_dose[active]) * (
    0.2 + panel$x1[active]^2 + ifelse(panel$time[active] > stats::median(panel$time), 0.7, -0.1)
  )

  fit <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weight_method = "entropy_balance",
    gps_model = "linear",
    entropy_balance_iterations = 500L,
    n_folds = 2,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 1,
    seed = 1268
  )

  expect_s3_class(fit, "cdmc_dr_fit")
  expect_identical(fit$weight_method, "entropy_balance")
  expect_identical(fit$entropy_balance_degree, 1L)
  expect_true(all(is.finite(fit$data$.cdmc_weight[fit$data$.cdmc_observed])))
  expect_true(any(abs(fit$data$.cdmc_weight[fit$data$.cdmc_observed] - 1) > 1e-6))
  expect_true(any(fit$data$.cdmc_dr_sample))
  expect_true(all(vapply(
    fit$baseline$fold_summaries,
    function(summary) !is.null(summary$entropy_balance_max_abs_balance),
    logical(1)
  )))
})

test_that("cross-fitted DR estimator can estimate internal kernel-balance weights", {
  panel <- simulate_cdmc_data(
    n_units = 30,
    n_times = 10,
    rank = 2,
    beta = 0.85,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.1,
    switch_off_prob = 0.5,
    seed = 1274
  )
  original_dose <- panel$dose
  active <- abs(original_dose) > 0
  panel$dose[active] <- sign(original_dose[active]) * (
    0.25 + sin(panel$x1[active]) + ifelse(panel$time[active] > stats::median(panel$time), 0.6, -0.1)
  )

  fit <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weight_method = "kernel_balance",
    gps_model = "linear",
    kernel_balance_centers = 6L,
    kernel_balance_iterations = 400L,
    n_folds = 2,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 1,
    seed = 1274
  )

  expect_s3_class(fit, "cdmc_dr_fit")
  expect_identical(fit$weight_method, "kernel_balance")
  expect_identical(fit$kernel_balance_degree, 1L)
  expect_identical(fit$kernel_balance_centers, 6L)
  expect_true(all(is.finite(fit$data$.cdmc_weight[fit$data$.cdmc_observed])))
  expect_true(any(abs(fit$data$.cdmc_weight[fit$data$.cdmc_observed] - 1) > 1e-6))
  expect_true(any(fit$data$.cdmc_dr_sample))
  expect_true(all(vapply(
    fit$baseline$fold_summaries,
    function(summary) !is.null(summary$kernel_balance_max_abs_balance),
    logical(1)
  )))
})

test_that("cross-fitted DR estimator can estimate internal adaptive-balance weights", {
  panel <- simulate_cdmc_data(
    n_units = 30,
    n_times = 10,
    rank = 2,
    beta = 0.85,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.1,
    switch_off_prob = 0.5,
    seed = 1284
  )
  original_dose <- panel$dose
  active <- abs(original_dose) > 0
  panel$dose[active] <- sign(original_dose[active]) * (
    0.25 + sin(panel$x1[active]) + ifelse(panel$time[active] > stats::median(panel$time), 0.6, -0.1)
  )

  fit <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weight_method = "adaptive_balance",
    adaptive_balance_methods = c("entropy_balance", "kernel_balance"),
    gps_model = "linear",
    kernel_balance_centers = 6L,
    entropy_balance_iterations = 400L,
    kernel_balance_iterations = 400L,
    n_folds = 2,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 1,
    seed = 1284
  )

  expect_s3_class(fit, "cdmc_dr_fit")
  expect_identical(fit$weight_method, "adaptive_balance")
  expect_identical(fit$adaptive_balance_methods, c("entropy_balance", "kernel_balance"))
  expect_true(all(is.finite(fit$data$.cdmc_weight[fit$data$.cdmc_observed])))
  expect_true(any(abs(fit$data$.cdmc_weight[fit$data$.cdmc_observed] - 1) > 1e-6))
  expect_true(any(fit$data$.cdmc_dr_sample))
  expect_true(all(vapply(
    fit$baseline$fold_summaries,
    function(summary) !is.null(summary$adaptive_balance_selected_method),
    logical(1)
  )))
  expect_true(all(vapply(
    fit$baseline$fold_summaries,
    function(summary) is.data.frame(summary$adaptive_balance_candidate_scores),
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

test_that("cross-fitted DR print reports fold lambda tuning guidance", {
  panel <- simulate_cdmc_data(
    n_units = 10,
    n_times = 8,
    rank = 2,
    beta = 0.85,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.18,
    switch_off_prob = 0.42,
    seed = 12680
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
    lambda = NULL,
    lambda_selection = "heuristic",
    rank_max = 2,
    washout = 0,
    lag_order = 0,
    seed = 12680
  )

  output <- paste(capture.output(print(fit)), collapse = "\n")

  expect_true(all(vapply(
    fit$baseline$fold_summaries,
    function(summary) identical(summary$lambda_method, "heuristic"),
    logical(1)
  )))
  expect_match(output, 'lambda selection: heuristic', fixed = TRUE)
  expect_match(output, 'empirical workflows should prefer lambda_selection = "cv"', fixed = TRUE)
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

test_that("internal stacked Gaussian GPS weights combine available nuisance learners", {
  testthat::skip_if_not_installed("rpart")

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
    seed = 1567
  )
  original_dose <- panel$dose
  active <- abs(original_dose) > 0
  interaction_score <- panel$x2 + panel$time / max(panel$time)
  panel$dose[active] <- sign(original_dose[active]) * (
    0.1 +
      panel$x1[active]^2 +
      ifelse(interaction_score[active] > stats::median(interaction_score), 0.6, -0.2)
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
    seed = 1567
  )

  fit_stack <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = c("x1", "x2"),
    weight_method = "gaussian_gps",
    gps_model = "stack",
    gps_stack_models = c("linear", "spline", "tree"),
    n_folds = 2,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 0,
    seed = 1567
  )

  expect_s3_class(fit_stack, "cdmc_dr_fit")
  expect_identical(fit_stack$gps_model, "stack")
  expect_identical(fit_stack$gps_stack_models, c("linear", "spline", "tree"))
  expect_true(all(is.finite(fit_stack$data$.cdmc_weight)))
  expect_true(any(vapply(
    fit_stack$baseline$fold_summaries,
    function(summary) identical(summary$gps_stack_models, c("linear", "spline", "tree")),
    logical(1)
  )))
  expect_true(any(abs(fit_stack$data$.cdmc_weight - fit_linear$data$.cdmc_weight) > 1e-6, na.rm = TRUE))
  expect_true(any(fit_stack$data$.cdmc_dr_sample))
})

test_that("internal boost Gaussian GPS weights support additive nonlinear assignment patterns", {
  testthat::skip_if_not_installed("gbm")

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
    seed = 1573
  )
  original_dose <- panel$dose
  active <- abs(original_dose) > 0
  interaction_score <- panel$x2 + panel$time / max(panel$time)
  panel$dose[active] <- sign(original_dose[active]) * (
    0.15 +
      panel$x1[active]^2 +
      0.5 * interaction_score[active] +
      0.25 * pmax(panel$x1[active], 0)
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
    seed = 1573
  )

  fit_boost <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = c("x1", "x2"),
    weight_method = "gaussian_gps",
    gps_model = "boost",
    gps_boost_trees = 40L,
    gps_boost_depth = 2L,
    gps_boost_shrinkage = 0.05,
    gps_boost_min_obs_node = 3L,
    n_folds = 2,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 0,
    seed = 1573
  )

  expect_s3_class(fit_boost, "cdmc_dr_fit")
  expect_identical(fit_boost$gps_model, "boost")
  expect_identical(fit_boost$gps_boost_trees, 40L)
  expect_identical(fit_boost$gps_boost_depth, 2L)
  expect_equal(fit_boost$gps_boost_shrinkage, 0.05)
  expect_identical(fit_boost$gps_boost_min_obs_node, 3L)
  expect_true(all(is.finite(fit_boost$data$.cdmc_weight)))
  expect_true(any(abs(fit_boost$data$.cdmc_weight - fit_linear$data$.cdmc_weight) > 1e-6, na.rm = TRUE))
  expect_true(all(vapply(
    fit_boost$baseline$fold_summaries,
    function(summary) identical(summary$gps_boost_trees, 40L),
    logical(1)
  )))
  expect_true(any(fit_boost$data$.cdmc_dr_sample))
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
test_that("PLR (Robinson) score returns linear lag coefficients", {
  panel <- simulate_cdmc_data(
    n_units = 25,
    n_times = 12,
    rank = 2,
    beta = 1,
    lag_beta = 0.3,
    n_covariates = 1,
    noise_sd = 0.03,
    switch_on_prob = 0.15,
    switch_off_prob = 0.4,
    seed = 2222
  )

  fit <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    n_folds = 2,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 1,
    dr_score = "plr",
    seed = 2222
  )

  expect_s3_class(fit, "cdmc_dr_fit")
  expect_identical(fit$dr_score, "plr")
  expect_true(all(c("dose_lag0", "dose_lag1") %in% names(fit$effect$coefficients)))
  expect_equal(unname(fit$effect$coefficients["dose_lag0"]), 1, tolerance = 0.8)
})

test_that("baseline_weighting='none' fits without GPS-weighted baseline", {
  panel <- simulate_cdmc_data(
    n_units = 20,
    n_times = 10,
    rank = 2,
    beta = 1,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 3333
  )

  fit <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    n_folds = 2,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 1,
    baseline_weighting = "none",
    seed = 3333
  )

  expect_s3_class(fit, "cdmc_dr_fit")
  expect_identical(fit$baseline_weighting, "none")
  expect_true("dose_lag0" %in% names(fit$effect$coefficients))
})

test_that("summary.cdmc_dr_fit returns a structured table", {
  panel <- simulate_cdmc_data(
    n_units = 20,
    n_times = 10,
    rank = 2,
    beta = 1,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 4444
  )

  fit <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    n_folds = 2,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 1,
    seed = 4444
  )

  s <- summary(fit)
  expect_s3_class(s, "summary.cdmc_dr_fit")
  expect_true(all(c("coefficients", "folds", "weight", "score", "dimensions") %in% names(s)))
  expect_true(nrow(s$coefficients) >= 1L)
  expect_identical(s$score$dr_score, "aipw")
  expect_silent(invisible(capture.output(print(s))))
})

test_that("cdmc_dr_score_contrast aligns AIPW and PLR coefficients", {
  panel <- simulate_cdmc_data(
    n_units = 25,
    n_times = 12,
    rank = 2,
    beta = 1,
    lag_beta = 0.3,
    n_covariates = 1,
    noise_sd = 0.03,
    switch_on_prob = 0.15,
    switch_off_prob = 0.4,
    seed = 5555
  )
  common <- list(
    data = panel, outcome = "y", dose = "dose", unit = "unit", time = "time",
    covariates = "x1", n_folds = 2, lambda = 0.2, rank_max = 3,
    washout = 0, lag_order = 1, seed = 5555
  )
  fit_aipw <- do.call(cdmc_dr_fit, c(common, list(dr_score = "aipw")))
  fit_plr  <- do.call(cdmc_dr_fit, c(common, list(dr_score = "plr")))

  contrast <- cdmc_dr_score_contrast(fit_aipw, fit_plr)
  expect_s3_class(contrast, "data.frame")
  expect_true(all(c("term", "aipw_estimate", "plr_estimate", "abs_diff",
                    "rel_diff", "sign_agree") %in% names(contrast)))
  expect_true("dose_lag0" %in% contrast$term)
  expect_equal(contrast$abs_diff, contrast$plr_estimate - contrast$aipw_estimate)

  expect_error(cdmc_dr_score_contrast(fit_aipw, "not a fit"))
})
