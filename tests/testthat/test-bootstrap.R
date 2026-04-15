test_that("bootstrap inference returns coefficient summaries", {
  panel <- simulate_cdmc_data(
    n_units = 8,
    n_times = 8,
    rank = 2,
    beta = 1,
    lag_beta = 0.3,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 101
  )

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = 0.2,
    rank_max = 3,
    washout = 1,
    lag_order = 1,
    seed = 101
  )

  boot <- cdmc_bootstrap(
    fit,
    n_boot = 4,
    statistics = c("coefficients", "average_tau"),
    seed = 101
  )

  expect_s3_class(boot, "cdmc_bootstrap")
  expect_equal(boot$n_boot, 4)
  expect_true(all(c("coef_dose_lag0", "coef_dose_lag1", "mean_tau_active_dose") %in% colnames(boot$draws)))
  expect_true(all(c("coef_dose_lag0", "coef_dose_lag1", "mean_tau_active_dose") %in% boot$summary$statistic))
})

test_that("bootstrap preserves joint-objective fit specifications", {
  panel <- simulate_cdmc_data(
    n_units = 8,
    n_times = 8,
    rank = 2,
    beta = 0.8,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 111
  )

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = 0.2,
    rank_max = 3,
    lag_order = 1,
    objective = "joint",
    seed = 111
  )

  boot <- cdmc_bootstrap(
    fit,
    n_boot = 2,
    statistics = c("coefficients", "average_tau"),
    seed = 111
  )

  expect_s3_class(boot, "cdmc_bootstrap")
  expect_equal(boot$lambda_method, "fixed")
  expect_true(all(c("coef_dose_lag0", "coef_dose_lag1", "mean_tau_active_dose") %in% colnames(boot$draws)))
})

test_that("bootstrap can rerun cross-validated lambda selection", {
  panel <- simulate_cdmc_data(
    n_units = 8,
    n_times = 7,
    rank = 2,
    beta = 0.8,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.45,
    seed = 202
  )

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda_selection = "cv",
    lambda_grid = c(0.4, 0.2),
    cv_rounds = 2,
    cv_block_size = 1,
    rank_max = 3,
    washout = 1,
    lag_order = 1,
    seed = 202
  )

  boot <- cdmc_bootstrap(
    fit,
    n_boot = 2,
    statistics = "coefficients",
    rerun_tuning = TRUE,
    seed = 202
  )

  expect_true(boot$rerun_tuning)
  expect_equal(boot$lambda_method, "cv")
  expect_true(boot$n_failures <= boot$n_boot)
})

test_that("bootstrap preserves weighted fit specifications", {
  panel <- simulate_cdmc_data(
    n_units = 8,
    n_times = 8,
    rank = 2,
    beta = 1,
    lag_beta = 0.3,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 3030
  )
  panel$w <- 0.5 + abs(panel$x1)

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weights = "w",
    lambda = 0.2,
    rank_max = 3,
    washout = 1,
    lag_order = 1,
    seed = 3030
  )

  boot <- cdmc_bootstrap(
    fit,
    n_boot = 2,
    statistics = "coefficients",
    seed = 3030
  )

  expect_true(fit$weight_supplied)
  expect_equal(fit$fit_control$weights, ".cdmc_weight")
  expect_true(boot$n_failures <= boot$n_boot)
})

test_that("bootstrap supports unbalanced panel fits", {
  panel <- simulate_cdmc_data(
    n_units = 8,
    n_times = 8,
    rank = 2,
    beta = 0.9,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 4040
  )
  drop_index <- head(which(panel$dose != 0), 2L)
  panel <- panel[-drop_index, ]

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = 0.2,
    rank_max = 3,
    washout = 1,
    lag_order = 1,
    seed = 4040
  )

  boot <- cdmc_bootstrap(
    fit,
    n_boot = 2,
    statistics = "coefficients",
    seed = 4040
  )

  expect_s3_class(boot, "cdmc_bootstrap")
  expect_true(boot$n_failures <= boot$n_boot)
})

test_that("bootstrap preserves spline second-stage fit specifications", {
  panel <- simulate_cdmc_data(
    n_units = 8,
    n_times = 8,
    rank = 2,
    beta = 0.3,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 3131
  )
  panel$y <- panel$y + 0.8 * panel$dose^2

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = 0.2,
    rank_max = 3,
    washout = 1,
    lag_order = 0,
    effect_model = "spline",
    effect_df = 3,
    seed = 3131
  )

  boot <- cdmc_bootstrap(
    fit,
    n_boot = 2,
    statistics = "coefficients",
    seed = 3131
  )

  expect_s3_class(boot, "cdmc_bootstrap")
  expect_true(any(startsWith(colnames(boot$draws), "coef_dose_lag0_ns")))
  expect_true(any(startsWith(boot$summary$statistic, "coef_dose_lag0_ns")))
  expect_true(boot$n_failures <= boot$n_boot)
})

test_that("bootstrap supports dose-response fits and prediction functionals", {
  panel <- simulate_cdmc_data(
    n_units = 10,
    n_times = 10,
    rank = 2,
    beta = 0.2,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    signed_dose = TRUE,
    negative_dose_prob = 0.5,
    seed = 7171
  )
  panel$y <- panel$y + 0.8 * panel$dose^2

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 0,
    seed = 7171
  )
  dr <- cdmc_dose_response(
    fit,
    model = "spline",
    lag_order = 0,
    df = 3,
    include_zero_dose = TRUE
  )

  boot <- cdmc_bootstrap(
    dr,
    n_boot = 2,
    statistics = c("coefficients", "prediction"),
    prediction_dose = c(-0.5, 0.5),
    prediction_type = "response",
    seed = 7171
  )

  expect_s3_class(boot, "cdmc_bootstrap")
  expect_true(any(startsWith(colnames(boot$draws), "coef_dose_lag0_ns")))
  expect_true(all(c("dose_response_response_1", "dose_response_response_2") %in% colnames(boot$draws)))
  expect_true(all(c("dose_response_response_1", "dose_response_response_2") %in% boot$summary$statistic))
  expect_true(boot$n_failures <= boot$n_boot)
})

test_that("bootstrap supports DR fits with internal gaussian GPS weights", {
  panel <- simulate_cdmc_data(
    n_units = 8,
    n_times = 8,
    rank = 2,
    beta = 0.9,
    lag_beta = 0.25,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 4040
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
    seed = 4040
  )

  boot <- cdmc_bootstrap(
    fit,
    n_boot = 3,
    statistics = c("coefficients", "average_tau_dr", "average_tau_linear"),
    seed = 4040
  )

  expect_s3_class(boot, "cdmc_bootstrap")
  expect_equal(boot$lambda_method, "fixed")
  expect_true(all(c(
    "coef_dose_lag0",
    "coef_dose_lag1",
    "mean_tau_dr",
    "mean_tau_dr_linear"
  ) %in% colnames(boot$draws)))
  expect_true(all(c(
    "coef_dose_lag0",
    "coef_dose_lag1",
    "mean_tau_dr",
    "mean_tau_dr_linear"
  ) %in% boot$summary$statistic))
})

test_that("bootstrap supports DR fits with internal spline Gaussian GPS weights", {
  panel <- simulate_cdmc_data(
    n_units = 8,
    n_times = 8,
    rank = 2,
    beta = 0.6,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 4141
  )
  original_dose <- panel$dose
  active <- abs(original_dose) > 0
  panel$dose[active] <- sign(original_dose[active]) * (0.25 + panel$x1[active]^2 + panel$time[active] / max(panel$time))
  panel$y <- panel$y + 0.6 * (panel$dose - original_dose)

  fit <- cdmc_dr_fit(
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
    seed = 4141
  )

  boot <- cdmc_bootstrap(
    fit,
    n_boot = 2,
    statistics = c("coefficients", "average_tau_dr"),
    seed = 4141
  )

  expect_s3_class(boot, "cdmc_bootstrap")
  expect_identical(fit$gps_model, "spline")
  expect_true(all(c("coef_dose_lag0", "mean_tau_dr") %in% colnames(boot$draws)))
  expect_true(all(c("coef_dose_lag0", "mean_tau_dr") %in% boot$summary$statistic))
  expect_true(boot$n_failures <= boot$n_boot)
})

test_that("bootstrap supports DR fits with internal GAM Gaussian GPS weights", {
  testthat::skip_if_not_installed("mgcv")

  panel <- simulate_cdmc_data(
    n_units = 8,
    n_times = 8,
    rank = 2,
    beta = 0.55,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 4191
  )
  original_dose <- panel$dose
  active <- abs(original_dose) > 0
  panel$dose[active] <- sign(original_dose[active]) * (
    0.25 + cos(panel$x1[active]) + (panel$time[active] / max(panel$time))^2
  )
  panel$y <- panel$y + 0.55 * (panel$dose - original_dose)

  fit <- cdmc_dr_fit(
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
    seed = 4191
  )

  boot <- cdmc_bootstrap(
    fit,
    n_boot = 2,
    statistics = c("coefficients", "average_tau_dr"),
    seed = 4191
  )

  expect_s3_class(boot, "cdmc_bootstrap")
  expect_identical(fit$gps_model, "gam")
  expect_true(all(c("coef_dose_lag0", "mean_tau_dr") %in% colnames(boot$draws)))
  expect_true(boot$n_failures <= boot$n_boot)
})

test_that("bootstrap supports DR fits with internal tree Gaussian GPS weights", {
  testthat::skip_if_not_installed("rpart")

  panel <- simulate_cdmc_data(
    n_units = 8,
    n_times = 8,
    rank = 2,
    beta = 0.55,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 4292
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
  panel$y <- panel$y + 0.55 * (panel$dose - original_dose)

  fit <- cdmc_dr_fit(
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
    seed = 4292
  )

  boot <- cdmc_bootstrap(
    fit,
    n_boot = 2,
    statistics = c("coefficients", "average_tau_dr"),
    seed = 4292
  )

  expect_s3_class(boot, "cdmc_bootstrap")
  expect_identical(fit$gps_model, "tree")
  expect_true(all(c("coef_dose_lag0", "mean_tau_dr") %in% colnames(boot$draws)))
  expect_true(boot$n_failures <= boot$n_boot)
})

test_that("bootstrap supports DR fits with internal forest Gaussian GPS weights", {
  testthat::skip_if_not_installed("ranger")

  panel <- simulate_cdmc_data(
    n_units = 8,
    n_times = 8,
    rank = 2,
    beta = 0.55,
    lag_beta = NULL,
    n_covariates = 2,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 4328
  )
  original_dose <- panel$dose
  active <- abs(original_dose) > 0
  interaction_score <- panel$x2 + panel$time / max(panel$time)
  panel$dose[active] <- sign(original_dose[active]) * (
    0.15 +
      ifelse(panel$x1[active] > stats::median(panel$x1), 0.9, 0.2) +
      ifelse(interaction_score[active] > stats::median(interaction_score), 0.7, -0.1)
  )
  panel$y <- panel$y + 0.55 * (panel$dose - original_dose)

  fit <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = c("x1", "x2"),
    weight_method = "gaussian_gps",
    gps_model = "forest",
    gps_forest_trees = 40L,
    gps_forest_mtry = 2L,
    gps_forest_min_node_size = 3L,
    n_folds = 2,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 0,
    seed = 4328
  )

  boot <- cdmc_bootstrap(
    fit,
    n_boot = 2,
    statistics = c("coefficients", "average_tau_dr"),
    seed = 4328
  )

  expect_s3_class(boot, "cdmc_bootstrap")
  expect_identical(fit$gps_model, "forest")
  expect_true(all(c("coef_dose_lag0", "mean_tau_dr") %in% colnames(boot$draws)))
  expect_true(boot$n_failures <= boot$n_boot)
})

test_that("bootstrap supports DR fits with internal kernel GPS weights", {
  panel <- simulate_cdmc_data(
    n_units = 8,
    n_times = 8,
    rank = 2,
    beta = 0.7,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 4343
  )
  original_dose <- panel$dose
  active <- abs(original_dose) > 0
  panel$dose[active] <- sign(original_dose[active]) * (
    abs(original_dose[active]) + 0.25 + panel$x1[active]^2 + rexp(sum(active), rate = 2)
  )
  panel$y <- panel$y + 0.6 * (panel$dose - original_dose)

  fit <- cdmc_dr_fit(
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
    seed = 4343
  )

  boot <- cdmc_bootstrap(
    fit,
    n_boot = 2,
    statistics = c("coefficients", "average_tau_dr"),
    seed = 4343
  )

  expect_s3_class(boot, "cdmc_bootstrap")
  expect_true(all(c("coef_dose_lag0", "mean_tau_dr") %in% colnames(boot$draws)))
  expect_true(boot$n_failures <= boot$n_boot)
})

test_that("bootstrap supports DR fits with internal CBPS weights", {
  testthat::skip_if_not_installed("CBPS")

  panel <- simulate_cdmc_data(
    n_units = 10,
    n_times = 8,
    rank = 2,
    beta = 0.75,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 4494
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
    lag_order = 0,
    seed = 4494
  )

  boot <- cdmc_bootstrap(
    fit,
    n_boot = 2,
    statistics = c("coefficients", "average_tau_dr"),
    seed = 4494
  )

  expect_s3_class(boot, "cdmc_bootstrap")
  expect_true(all(c("coef_dose_lag0", "mean_tau_dr") %in% colnames(boot$draws)))
  expect_true(boot$n_failures <= boot$n_boot)
})

test_that("bootstrap supports DR fits with external weights", {
  panel <- simulate_cdmc_data(
    n_units = 8,
    n_times = 8,
    rank = 2,
    beta = 0.85,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 5050
  )
  panel$w <- 0.5 + abs(panel$x1)

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
    seed = 5050
  )

  boot <- cdmc_bootstrap(
    fit,
    n_boot = 2,
    statistics = c("coefficients", "average_tau_dr"),
    seed = 5050
  )

  expect_s3_class(boot, "cdmc_bootstrap")
  expect_equal(boot$lambda_method, "fixed")
  expect_true(all(c("coef_dose_lag0", "coef_dose_lag1", "mean_tau_dr") %in% colnames(boot$draws)))
  expect_true(boot$n_failures <= boot$n_boot)
})

test_that("bootstrap supports richer DR lag-path and contrast summaries", {
  panel <- simulate_cdmc_data(
    n_units = 8,
    n_times = 8,
    rank = 2,
    beta = 0.9,
    lag_beta = 0.25,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 6060
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
    seed = 6060
  )

  boot <- cdmc_bootstrap(
    fit,
    n_boot = 2,
    statistics = c("lag_average_tau_dr", "dose_contrast_dr"),
    contrast_history = c(dose_lag0 = 1, dose_lag1 = 0.5),
    reference_history = c(dose_lag0 = 0, dose_lag1 = 0),
    seed = 6060
  )

  expect_s3_class(boot, "cdmc_bootstrap")
  expect_true(all(c(
    "mean_tau_dr_dose_lag0",
    "mean_tau_dr_dose_lag1",
    "dr_path_contrast"
  ) %in% colnames(boot$draws)))
  expect_true(all(c(
    "mean_tau_dr_dose_lag0",
    "mean_tau_dr_dose_lag1",
    "dr_path_contrast"
  ) %in% boot$summary$statistic))
})

test_that("bootstrap supports placebo diagnostics", {
  panel <- simulate_cdmc_data(
    n_units = 10,
    n_times = 10,
    rank = 2,
    beta = 1,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.03,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 7070
  )

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = 0.2,
    rank_max = 3,
    washout = 1,
    lag_order = 1,
    seed = 7070
  )
  placebo <- cdmc_placebo_test(fit, periods = -2:0)

  boot <- cdmc_bootstrap(placebo, n_boot = 2, seed = 7070)

  expect_s3_class(boot, "cdmc_bootstrap")
  expect_equal(boot$lambda_method, "fixed")
  expect_true("mean_placebo_tau" %in% colnames(boot$draws))
  expect_true("mean_placebo_tau" %in% boot$summary$statistic)
})

test_that("bootstrap supports refit carryover diagnostics", {
  panel <- simulate_cdmc_data(
    n_units = 10,
    n_times = 10,
    rank = 2,
    beta = 1,
    lag_beta = 0.4,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.2,
    switch_off_prob = 0.45,
    seed = 8080
  )

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 1,
    seed = 8080
  )
  diagnostic <- cdmc_carryover_refit_test(fit, periods = 1)

  boot <- cdmc_bootstrap(diagnostic, n_boot = 2, seed = 8080)

  expect_s3_class(boot, "cdmc_bootstrap")
  expect_true("mean_refit_carryover_tau" %in% colnames(boot$draws))
  expect_true("mean_refit_carryover_tau" %in% boot$summary$statistic)
})

test_that("bootstrap supports residual carryover diagnostics", {
  panel <- simulate_cdmc_data(
    n_units = 10,
    n_times = 10,
    rank = 2,
    beta = 1,
    lag_beta = 0.35,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.2,
    switch_off_prob = 0.45,
    seed = 8181
  )

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 1,
    seed = 8181
  )
  diagnostic <- cdmc_carryover_test(fit, periods = 1)

  boot <- cdmc_bootstrap(diagnostic, n_boot = 2, seed = 8181)

  expect_s3_class(boot, "cdmc_bootstrap")
  expect_true("mean_carryover_tau" %in% colnames(boot$draws))
  expect_true("mean_carryover_tau" %in% boot$summary$statistic)
})

test_that("bootstrap supports joint placebo diagnostics", {
  panel <- simulate_cdmc_data(
    n_units = 10,
    n_times = 12,
    rank = 2,
    beta = 1,
    lag_beta = 0.15,
    n_covariates = 1,
    noise_sd = 0.02,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 9090
  )

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = 0.2,
    rank_max = 3,
    washout = 1,
    lag_order = 1,
    seed = 9090
  )
  joint <- cdmc_joint_placebo_test(fit, periods = -2:-1)

  boot <- cdmc_bootstrap(joint, n_boot = 2, seed = 9090)

  expect_s3_class(boot, "cdmc_bootstrap")
  expect_true(all(c(
    "joint_placebo_mean_tau_period_m2",
    "joint_placebo_mean_tau_period_m1",
    "joint_placebo_p_value"
  ) %in% colnames(boot$draws)))
  expect_true(all(c(
    "joint_placebo_mean_tau_period_m2",
    "joint_placebo_mean_tau_period_m1",
    "joint_placebo_p_value"
  ) %in% boot$summary$statistic))
})

test_that("bootstrap supports equivalence diagnostics", {
  panel <- simulate_cdmc_data(
    n_units = 10,
    n_times = 12,
    rank = 2,
    beta = 1,
    lag_beta = 0.1,
    n_covariates = 1,
    noise_sd = 0.02,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 9141
  )

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = 0.2,
    rank_max = 3,
    washout = 1,
    lag_order = 1,
    seed = 9141
  )
  placebo <- cdmc_placebo_test(fit, periods = -2:-1)
  equivalence <- cdmc_equivalence_test(placebo, margin = 0.2)

  boot <- cdmc_bootstrap(equivalence, n_boot = 2, seed = 9141)

  expect_s3_class(boot, "cdmc_bootstrap")
  expect_true(all(c("equivalence_mean_tau", "equivalence_p_value") %in% colnames(boot$draws)))
  expect_true(all(c("equivalence_mean_tau", "equivalence_p_value") %in% boot$summary$statistic))
})

test_that("bootstrap supports SCIA diagnostics", {
  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 10,
    rank = 2,
    beta = 0.8,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 9191
  )

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 1,
    seed = 9191
  )
  diagnostic <- cdmc_scia_test(fit, lags = 1)

  boot <- cdmc_bootstrap(diagnostic, n_boot = 2, seed = 9191)

  expect_s3_class(boot, "cdmc_bootstrap")
  expect_true(all(c("scia_f_statistic", "scia_p_value") %in% colnames(boot$draws)))
  expect_true(all(c("scia_f_statistic", "scia_p_value") %in% boot$summary$statistic))
})
