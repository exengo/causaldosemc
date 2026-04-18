test_that("cdmc_fit can tune lambda with blocked cross-validation", {
  panel <- simulate_cdmc_data(
    n_units = 10,
    n_times = 9,
    rank = 2,
    beta = 0.9,
    lag_beta = 0.25,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 99
  )

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = NULL,
    lambda_selection = "cv",
    lambda_grid = c(0.4, 0.2, 0.1),
    cv_rounds = 2,
    cv_block_size = 1,
    rank_max = 3,
    washout = 1,
    lag_order = 1,
    seed = 99
  )

  expect_s3_class(fit, "cdmc_fit")
  expect_equal(fit$lambda_tuning$method, "cv")
  expect_true(fit$lambda %in% c(0.4, 0.2, 0.1))
  expect_equal(fit$lambda, fit$lambda_tuning$selected_lambda)
  expect_equal(length(fit$lambda_tuning$lambda_grid), 3)
  expect_true(all(is.finite(fit$lambda_tuning$mean_scores)))
  expect_true(all(fit$lambda_tuning$holdout_counts > 0))
})

test_that("cdmc_fit validates cv_workers", {
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
    seed = 109
  )

  expect_error(
    cdmc_fit(
      data = panel,
      outcome = "y",
      dose = "dose",
      unit = "unit",
      time = "time",
      covariates = "x1",
      lambda = NULL,
      lambda_selection = "cv",
      lambda_grid = c(0.4, 0.2),
      cv_rounds = 2,
      cv_block_size = 1,
      cv_workers = 0,
      rank_max = 3,
      washout = 1,
      lag_order = 1,
      seed = 109
    ),
    "cv_workers must be a positive integer"
  )

  expect_error(
    cdmc_fit(
      data = panel,
      outcome = "y",
      dose = "dose",
      unit = "unit",
      time = "time",
      covariates = "x1",
      lambda = NULL,
      lambda_selection = "cv",
      lambda_grid = c(0.4, 0.2),
      cv_rounds = 2,
      cv_block_size = 1,
      cv_workers = 1.5,
      rank_max = 3,
      washout = 1,
      lag_order = 1,
      seed = 109
    ),
    "cv_workers must be a positive integer"
  )

  expect_error(
    cdmc_fit(
      data = panel,
      outcome = "y",
      dose = "dose",
      unit = "unit",
      time = "time",
      covariates = "x1",
      lambda = NULL,
      lambda_selection = "cv",
      lambda_grid = c(0.4, 0.2),
      cv_rounds = 2,
      cv_block_size = 1,
      cv_top_k = 0,
      rank_max = 3,
      washout = 1,
      lag_order = 1,
      seed = 109
    ),
    "cv_top_k must be NULL or a positive integer"
  )

  expect_error(
    cdmc_fit(
      data = panel,
      outcome = "y",
      dose = "dose",
      unit = "unit",
      time = "time",
      covariates = "x1",
      lambda = NULL,
      lambda_selection = "cv",
      lambda_grid = c(0.4, 0.3, 0.2, 0.1),
      cv_rounds = 2,
      cv_block_size = 1,
      cv_coarse_to_fine = TRUE,
      cv_coarse_nlambda = 2,
      rank_max = 3,
      washout = 1,
      lag_order = 1,
      seed = 109
    ),
    "cv_coarse_nlambda must be NULL or an integer >= 3"
  )
})

test_that("cdmc_fit parallel cv tuning matches sequential tuning", {
  skip_on_os("windows")

  panel <- simulate_cdmc_data(
    n_units = 10,
    n_times = 9,
    rank = 2,
    beta = 0.9,
    lag_beta = 0.25,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 209
  )

  fit_seq <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = NULL,
    lambda_selection = "cv",
    lambda_grid = c(0.4, 0.2, 0.1),
    cv_rounds = 2,
    cv_block_size = 1,
    cv_workers = 1,
    rank_max = 3,
    washout = 1,
    lag_order = 1,
    seed = 209
  )

  fit_par <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = NULL,
    lambda_selection = "cv",
    lambda_grid = c(0.4, 0.2, 0.1),
    cv_rounds = 2,
    cv_block_size = 1,
    cv_workers = 2,
    rank_max = 3,
    washout = 1,
    lag_order = 1,
    seed = 209
  )

  expect_equal(fit_seq$lambda, fit_par$lambda)
  expect_equal(fit_seq$lambda_tuning$scores, fit_par$lambda_tuning$scores)
  expect_equal(fit_seq$lambda_tuning$mean_scores, fit_par$lambda_tuning$mean_scores)
  expect_false(isTRUE(fit_seq$lambda_tuning$cv_parallel))
  expect_true(isTRUE(fit_par$lambda_tuning$cv_parallel))
  expect_false(isTRUE(fit_seq$lambda_tuning$cv_warm_starts))
  expect_false(isTRUE(fit_par$lambda_tuning$cv_warm_starts))
  expect_equal(fit_seq$lambda_tuning$cv_workers, 1L)
  expect_equal(fit_par$lambda_tuning$cv_workers, 2L)
})

test_that("cdmc_fit can opt into sequential cv warm starts", {
  panel <- simulate_cdmc_data(
    n_units = 10,
    n_times = 9,
    rank = 2,
    beta = 0.9,
    lag_beta = 0.25,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 239
  )

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = NULL,
    lambda_selection = "cv",
    lambda_grid = c(0.4, 0.2, 0.1),
    cv_rounds = 2,
    cv_block_size = 1,
    cv_workers = 1,
    cv_warm_starts = TRUE,
    rank_max = 3,
    washout = 1,
    lag_order = 1,
    seed = 239
  )

  expect_equal(fit$lambda_tuning$method, "cv")
  expect_true(isTRUE(fit$lambda_tuning$cv_warm_starts))
  expect_true(is.finite(fit$lambda_tuning$selected_lambda))
})

test_that("cdmc_fit can screen lambda candidates during cv", {
  panel <- simulate_cdmc_data(
    n_units = 10,
    n_times = 9,
    rank = 2,
    beta = 0.9,
    lag_beta = 0.25,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 249
  )

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = NULL,
    lambda_selection = "cv",
    lambda_grid = c(0.5, 0.35, 0.25, 0.18, 0.12),
    cv_rounds = 3,
    cv_block_size = 1,
    cv_top_k = 2,
    cv_workers = 1,
    rank_max = 3,
    washout = 1,
    lag_order = 1,
    seed = 249
  )

  expect_equal(fit$lambda_tuning$method, "cv")
  expect_true(isTRUE(fit$lambda_tuning$cv_screening))
  expect_equal(fit$lambda_tuning$cv_top_k, 2L)
  expect_true(any(fit$lambda_tuning$cv_rounds_evaluated < fit$lambda_tuning$cv_rounds))
  expect_equal(sum(is.finite(fit$lambda_tuning$mean_scores)), 2L)
})

test_that("cdmc_fit can run coarse-to-fine cv tuning", {
  panel <- simulate_cdmc_data(
    n_units = 10,
    n_times = 9,
    rank = 2,
    beta = 0.9,
    lag_beta = 0.25,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 259
  )

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = NULL,
    lambda_selection = "cv",
    lambda_grid = c(0.6, 0.45, 0.35, 0.27, 0.2, 0.15, 0.11),
    cv_rounds = 2,
    cv_block_size = 1,
    cv_coarse_to_fine = TRUE,
    cv_coarse_nlambda = 4,
    cv_workers = 1,
    rank_max = 3,
    washout = 1,
    lag_order = 1,
    seed = 259
  )

  expect_true(isTRUE(fit$lambda_tuning$cv_coarse_to_fine))
  expect_equal(fit$lambda_tuning$cv_coarse_nlambda, 4L)
  expect_true(length(fit$lambda_tuning$cv_coarse_grid) <= 4L)
  expect_true(length(fit$lambda_tuning$cv_fine_grid) <= 3L)
  expect_true(fit$lambda %in% fit$lambda_tuning$cv_fine_grid)
})

test_that("cdmc_fit can tune lambda with observation weights", {
  panel <- simulate_cdmc_data(
    n_units = 10,
    n_times = 9,
    rank = 2,
    beta = 0.9,
    lag_beta = 0.25,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 199
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
    lambda = NULL,
    lambda_selection = "cv",
    lambda_grid = c(0.4, 0.2, 0.1),
    cv_rounds = 2,
    cv_block_size = 1,
    rank_max = 3,
    washout = 1,
    lag_order = 1,
    seed = 199
  )

  expect_s3_class(fit, "cdmc_fit")
  expect_true(fit$weight_supplied)
  expect_equal(fit$lambda_tuning$method, "cv")
  expect_true(all(is.finite(fit$lambda_tuning$mean_scores)))
})
