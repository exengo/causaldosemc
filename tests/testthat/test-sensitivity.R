manual_sensitivity_refit_source_fit <- function(source_fit, data) {
  if (inherits(source_fit, "cdmc_dr_fit")) {
    return(do.call(
      cdmc_dr_fit,
      list(
        data = data,
        outcome = source_fit$outcome,
        dose = source_fit$dose,
        unit = source_fit$unit,
        time = source_fit$time,
        covariates = source_fit$covariates,
        weights = source_fit$fit_control$weights,
        weight_method = source_fit$weight_method,
        weight_covariates = source_fit$weight_covariates,
        gps_time_effects = source_fit$gps_time_effects,
        gps_model = source_fit$gps_model,
        gps_df = if (is.null(source_fit$gps_df)) 4L else source_fit$gps_df,
        gps_spline_covariates = source_fit$gps_spline_covariates,
        gps_stack_models = source_fit$gps_stack_models,
        gps_bandwidth = source_fit$gps_bandwidth,
        gps_forest_trees = source_fit$gps_forest_trees,
        gps_forest_mtry = source_fit$gps_forest_mtry,
        gps_forest_min_node_size = source_fit$gps_forest_min_node_size,
        gps_boost_trees = source_fit$gps_boost_trees,
        gps_boost_depth = source_fit$gps_boost_depth,
        gps_boost_shrinkage = source_fit$gps_boost_shrinkage,
        gps_boost_min_obs_node = source_fit$gps_boost_min_obs_node,
        adaptive_balance_methods = source_fit$adaptive_balance_methods,
        entropy_balance_degree = source_fit$entropy_balance_degree,
        entropy_balance_standardize = source_fit$entropy_balance_standardize,
        entropy_balance_iterations = source_fit$entropy_balance_iterations,
        entropy_balance_reltol = source_fit$entropy_balance_reltol,
        kernel_balance_degree = source_fit$kernel_balance_degree,
        kernel_balance_centers = source_fit$kernel_balance_centers,
        kernel_balance_bandwidth = source_fit$kernel_balance_bandwidth,
        kernel_balance_standardize = source_fit$kernel_balance_standardize,
        kernel_balance_iterations = source_fit$kernel_balance_iterations,
        kernel_balance_reltol = source_fit$kernel_balance_reltol,
        stabilize_weights = source_fit$stabilize_weights,
        max_weight = source_fit$max_weight,
        n_folds = source_fit$n_folds,
        fold_assignments = source_fit$fold_assignments,
        lambda = source_fit$lambda,
        rank_max = source_fit$rank_max,
        lambda_fraction = source_fit$fit_control$lambda_fraction,
        lambda_selection = source_fit$lambda_selection,
        lambda_grid = source_fit$fit_control$lambda_grid,
        nlambda = source_fit$fit_control$nlambda,
        lambda_min_ratio = source_fit$fit_control$lambda_min_ratio,
        cv_rounds = source_fit$fit_control$cv_rounds,
        cv_block_size = source_fit$fit_control$cv_block_size,
        washout = source_fit$washout,
        lag_order = source_fit$lag_order,
        outer_maxit = source_fit$fit_control$outer_maxit,
        fe_maxit = source_fit$fit_control$fe_maxit,
        soft_maxit = source_fit$fit_control$soft_maxit,
        tol = source_fit$fit_control$tol,
        fe_tol = source_fit$fit_control$fe_tol,
        zero_tolerance = source_fit$zero_tolerance,
        verbose = FALSE
      )
    ))
  }

  cdmc_fit(
    data = data,
    outcome = source_fit$outcome,
    dose = source_fit$dose,
    unit = source_fit$unit,
    time = source_fit$time,
    covariates = source_fit$covariates,
    weights = source_fit$fit_control$weights,
    lambda = source_fit$lambda,
    rank_max = source_fit$rank_max,
    washout = source_fit$washout,
    lag_order = source_fit$lag_order,
    effect_model = source_fit$effect_model,
    effect_df = if (is.null(source_fit$effect_df)) 4L else source_fit$effect_df,
    outer_maxit = source_fit$fit_control$outer_maxit,
    fe_maxit = source_fit$fit_control$fe_maxit,
    soft_maxit = source_fit$fit_control$soft_maxit,
    tol = source_fit$fit_control$tol,
    fe_tol = source_fit$fit_control$fe_tol,
    zero_tolerance = source_fit$zero_tolerance,
    verbose = FALSE,
    objective = if (is.null(source_fit$objective)) "staged" else source_fit$objective
  )
}

manual_sensitivity_target <- function(object, refit_fit) {
  if (inherits(object, "cdmc_fit")) {
    return(refit_fit)
  }

  if (inherits(object, "cdmc_dr_fit")) {
    return(refit_fit)
  }

  if (inherits(object, "cdmc_dose_response")) {
    return(cdmc_dose_response(
      refit_fit,
      model = object$model,
      lag_order = object$lag_order,
      df = if (is.null(object$design_info$df)) 4L else object$design_info$df,
      forest_trees = object$forest_trees %||% 200L,
      forest_mtry = object$forest_mtry,
      forest_min_node_size = object$forest_min_node_size,
      include_zero_dose = object$include_zero_dose,
      weights = object$weights
    ))
  }

  if (inherits(object, "cdmc_dynamic_estimand")) {
    source_target <- if (inherits(object$source_object, "cdmc_dose_response")) {
      cdmc_dose_response(
        refit_fit,
        model = object$source_object$model,
        lag_order = object$source_object$lag_order,
        df = if (is.null(object$source_object$design_info$df)) 4L else object$source_object$design_info$df,
        forest_trees = object$source_object$forest_trees %||% 200L,
        forest_mtry = object$source_object$forest_mtry,
        forest_min_node_size = object$source_object$forest_min_node_size,
        include_zero_dose = object$source_object$include_zero_dose,
        weights = object$source_object$weights
      )
    } else {
      refit_fit
    }

    return(cdmc_dynamic_estimand(
      source_target,
      history = object$history,
      reference_history = object$reference_history,
      type = object$type,
      aggregate = object$aggregate,
      path_weights = object$path_weights,
      eps = object$eps,
      labels = object$labels
    ))
  }

  stop(
    "manual sensitivity helper only supports cdmc_fit, cdmc_dr_fit, cdmc_dose_response, and cdmc_dynamic_estimand objects.",
    call. = FALSE
  )
}

manual_sensitivity_statistics <- function(object, prediction_history = NULL, prediction_type = "response") {
  if (inherits(object, "cdmc_fit")) {
    out <- numeric(0)
    coefficients <- object$effect$coefficients
    if (length(coefficients) > 0L) {
      names(coefficients) <- paste0("coef_", names(coefficients))
      out <- c(out, coefficients)
    }

    active_mask <- abs(object$dose_matrix) > object$zero_tolerance
    out <- c(
      out,
      mean_tau_active_dose = if (any(active_mask)) {
        if (isTRUE(object$weight_supplied)) {
          stats::weighted.mean(object$effect$tau[active_mask], w = object$weight_matrix[active_mask])
        } else {
          mean(object$effect$tau[active_mask])
        }
      } else {
        NA_real_
      }
    )

    return(out)
  }

  if (inherits(object, "cdmc_dr_fit")) {
    out <- numeric(0)
    coefficients <- object$effect$coefficients
    if (length(coefficients) > 0L) {
      names(coefficients) <- paste0("coef_", names(coefficients))
      out <- c(out, coefficients)
    }

    dr_mask <- object$effect$sample_mask
    out <- c(
      out,
      mean_tau_dr = if (any(dr_mask, na.rm = TRUE)) mean(object$effect$tau_dr[dr_mask], na.rm = TRUE) else NA_real_,
      mean_tau_dr_linear = if (any(dr_mask, na.rm = TRUE)) mean(object$effect$fitted[dr_mask], na.rm = TRUE) else NA_real_
    )

    return(out)
  }

  if (inherits(object, "cdmc_dose_response")) {
    prediction <- predict(object, history = prediction_history, type = prediction_type)
    estimates <- prediction$estimate
    names(estimates) <- paste0("dose_response_", prediction_type, "_", seq_along(estimates))
    return(estimates)
  }

  if (inherits(object, "cdmc_dynamic_estimand")) {
    return(object$estimates)
  }

  stop(
    "manual sensitivity helper only supports cdmc_fit, cdmc_dr_fit, cdmc_dose_response, and cdmc_dynamic_estimand objects.",
    call. = FALSE
  )
}

manual_sensitivity_optimization_radius <- function(
  object,
  statistic_name,
  step,
  prediction_history = NULL,
  prediction_type = "response"
) {
  source_fit <- cdmc_sensitivity_source_object(object)
  full_data <- source_fit$data
  observed_rows <- rep(TRUE, nrow(full_data))
  if (".cdmc_observed" %in% names(full_data)) {
    observed_rows <- !is.na(full_data$.cdmc_observed) & full_data$.cdmc_observed
  }

  data_columns <- cdmc_sensitivity_original_columns(object)
  fit_data <- full_data[observed_rows, data_columns, drop = FALSE]
  row_lookup <- integer(length(observed_rows))
  row_lookup[observed_rows] <- seq_len(sum(observed_rows))

  optimization_mask <- cdmc_sensitivity_bounds_optimization_mask(source_fit)
  perturbation_rows <- row_lookup[observed_rows & as.vector(t(optimization_mask))]
  perturbation_rows <- perturbation_rows[perturbation_rows > 0L]

  base_statistics <- manual_sensitivity_statistics(
    object,
    prediction_history = prediction_history,
    prediction_type = prediction_type
  )
  base_outcome <- fit_data[[source_fit$outcome]]
  derivatives <- numeric(length(perturbation_rows))

  for (index in seq_along(perturbation_rows)) {
    perturbed_data <- fit_data
    perturbed_data[[source_fit$outcome]][perturbation_rows[[index]]] <-
      base_outcome[[perturbation_rows[[index]]]] + step
    refit_fit <- manual_sensitivity_refit_source_fit(source_fit, perturbed_data)
    refit_target <- manual_sensitivity_target(object, refit_fit)
    refit_statistics <- manual_sensitivity_statistics(
      refit_target,
      prediction_history = prediction_history,
      prediction_type = prediction_type
    )
    derivatives[[index]] <- (refit_statistics[[statistic_name]] - base_statistics[[statistic_name]]) / step
  }

  list(
    multiplier = sum(abs(derivatives)),
    n_cells = length(perturbation_rows)
  )
}

test_that("sensitivity scan summarizes washout and weight scenarios", {
  panel <- simulate_cdmc_data(
    n_units = 10,
    n_times = 10,
    rank = 2,
    beta = 1,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 1111
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
    seed = 1111
  )

  scan <- cdmc_sensitivity_scan(
    fit,
    washout_grid = 0:1,
    zero_tolerance_grid = fit$zero_tolerance,
    include_unweighted = TRUE
  )

  expect_s3_class(scan, "cdmc_sensitivity_scan")
  expect_equal(nrow(scan$results), 4)
  expect_true(all(c("current", "unweighted") %in% scan$results$weight_scenario))
  expect_true(all(scan$results$support_ok))
  expect_true(all(scan$results$fit_success))
  expect_true(all(c("coef_dose_lag0", "coef_dose_lag1", "mean_tau_active_dose") %in% names(scan$results)))
  expect_true(all(c(
    "support_fraction",
    "support_loss",
    "optimization_cells",
    "optimization_fraction",
    "delta_eligible_controls",
    "delta_support_fraction",
    "delta_optimization_cells",
    "delta_optimization_fraction",
    "delta_mean_tau_active_dose"
  ) %in% names(scan$results)))
  expect_true(all(scan$results$objective == "staged"))
  expect_true(all(scan$results$optimization_sample == "eligible_zero_dose"))
  expect_equal(scan$reference$washout, fit$washout)
  expect_equal(scan$reference$weight_scenario, "current")
  expect_true(all(c("support_summary", "washout_summary", "weight_summary") %in% names(scan)))
  expect_equal(nrow(scan$support_summary), 4)
  expect_equal(nrow(scan$washout_summary), 2)
  expect_equal(nrow(scan$weight_summary), 2)
  expect_true(all(scan$washout_summary$n_washout == 2))
  expect_true(all(scan$weight_summary$n_weight_scenarios == 2))
  expect_s3_class(summary(scan), "summary.cdmc_sensitivity_scan")
})

test_that("sensitivity scan validates workers argument", {
  panel <- simulate_cdmc_data(
    n_units = 8,
    n_times = 8,
    rank = 2,
    beta = 1,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 1131
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
    seed = 1131
  )

  expect_error(
    cdmc_sensitivity_scan(fit, washout_grid = 0:1, workers = 0),
    "workers must be a positive integer"
  )
  expect_error(
    cdmc_sensitivity_scan(fit, washout_grid = 0:1, workers = 1.5),
    "workers must be a positive integer"
  )
})

test_that("sensitivity scan parallel replay matches sequential replay", {
  skip_on_os("windows")

  panel <- simulate_cdmc_data(
    n_units = 8,
    n_times = 8,
    rank = 2,
    beta = 1,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 1141
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
    seed = 1141
  )

  scan_seq <- cdmc_sensitivity_scan(
    fit,
    washout_grid = 0:1,
    zero_tolerance_grid = fit$zero_tolerance,
    include_unweighted = TRUE,
    workers = 1
  )
  scan_par <- cdmc_sensitivity_scan(
    fit,
    washout_grid = 0:1,
    zero_tolerance_grid = fit$zero_tolerance,
    include_unweighted = TRUE,
    workers = 2
  )

  expect_false(isTRUE(scan_seq$parallel))
  expect_equal(scan_seq$workers, 1L)
  expect_true(isTRUE(scan_par$parallel))
  expect_equal(scan_par$workers, 2L)
  expect_equal(scan_seq$results, scan_par$results)
  expect_equal(scan_seq$support_summary, scan_par$support_summary)
  expect_equal(scan_seq$washout_summary, scan_par$washout_summary)
  expect_equal(scan_seq$weight_summary, scan_par$weight_summary)
})

test_that("sensitivity scan records support failures without stopping", {
  panel <- expand.grid(
    unit = paste0("u", 1:4),
    time = 1:4,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  panel$dose <- c(
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1
  )
  panel$x1 <- rep(c(-1, -0.5, 0.5, 1), each = 4)
  unit_effect <- c(u1 = -0.2, u2 = 0, u3 = 0.2, u4 = 0.4)
  time_effect <- c(`1` = -0.1, `2` = 0, `3` = 0.1, `4` = 0.2)
  panel$y <- 1 + unit_effect[panel$unit] + time_effect[as.character(panel$time)] + 0.5 * panel$dose

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = 0.1,
    rank_max = 1,
    washout = 0,
    lag_order = 0,
    seed = 123
  )

  scan <- cdmc_sensitivity_scan(
    fit,
    washout_grid = c(0, 3),
    zero_tolerance_grid = fit$zero_tolerance,
    include_unweighted = FALSE
  )

  failing_row <- scan$results[scan$results$washout == 3, , drop = FALSE]

  expect_s3_class(scan, "cdmc_sensitivity_scan")
  expect_equal(nrow(scan$results), 2)
  expect_false(failing_row$support_ok)
  expect_false(failing_row$fit_success)
  expect_match(failing_row$fit_error, "eligible zero-dose observation")
  expect_true("support_loss" %in% names(scan$support_summary))
  expect_equal(scan$washout_summary$successful_washout, 1)
  expect_true(is.na(failing_row$delta_mean_tau_active_dose))
})

test_that("sensitivity scan distinguishes joint-objective optimization support", {
  panel <- expand.grid(
    unit = paste0("u", 1:4),
    time = 1:4,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  panel$dose <- c(
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1
  )
  panel$x1 <- rep(c(-1, -0.5, 0.5, 1), each = 4)
  unit_effect <- c(u1 = -0.2, u2 = 0, u3 = 0.2, u4 = 0.4)
  time_effect <- c(`1` = -0.1, `2` = 0, `3` = 0.1, `4` = 0.2)
  panel$y <- 1 + unit_effect[panel$unit] + time_effect[as.character(panel$time)] + 0.5 * panel$dose

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = 0.1,
    rank_max = 1,
    washout = 0,
    lag_order = 0,
    objective = "joint",
    seed = 123
  )

  scan <- cdmc_sensitivity_scan(
    fit,
    washout_grid = c(0, 3),
    zero_tolerance_grid = fit$zero_tolerance,
    include_unweighted = FALSE
  )

  washout_three <- scan$results[scan$results$washout == 3, , drop = FALSE]

  expect_s3_class(scan, "cdmc_sensitivity_scan")
  expect_true(all(scan$results$objective == "joint"))
  expect_true(all(scan$results$optimization_sample == "observed_panel"))
  expect_true(all(scan$results$support_ok))
  expect_true(all(scan$results$fit_success))
  expect_false(washout_three$control_support_ok)
  expect_true(washout_three$support_ok)
  expect_true(washout_three$fit_success)
  expect_equal(washout_three$optimization_cells, scan$reference$optimization_cells)
  expect_equal(washout_three$delta_optimization_cells, 0)
  expect_equal(washout_three$delta_optimization_fraction, 0)
})

test_that("sensitivity scan supports DR fits", {
  panel <- simulate_cdmc_data(
    n_units = 16,
    n_times = 10,
    rank = 2,
    beta = 0.9,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.18,
    switch_off_prob = 0.42,
    seed = 1717
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
    seed = 1717
  )

  scan <- cdmc_sensitivity_scan(
    fit,
    washout_grid = 0:1,
    zero_tolerance_grid = fit$zero_tolerance
  )

  expect_s3_class(scan, "cdmc_sensitivity_scan")
  expect_equal(scan$reference$source_class, "cdmc_dr_fit")
  expect_equal(scan$reference$target_class, "cdmc_dr_fit")
  expect_true(all(c(
    "coef_dose_lag0",
    "coef_dose_lag1",
    "mean_tau_dr",
    "mean_tau_dr_linear",
    "dr_sample_cells",
    "delta_mean_tau_dr"
  ) %in% names(scan$results)))
  expect_true(all(scan$results$weighted_fit))
  expect_true(all(scan$results$weight_mode == "gaussian_gps"))
  expect_equal(nrow(scan$washout_summary), 1)
})

test_that("sensitivity scan can compare kernel and Gaussian GPS nuisance models", {
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
    seed = 1767
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
    seed = 1767
  )

  scan <- cdmc_sensitivity_scan(
    fit,
    washout_grid = fit$washout,
    zero_tolerance_grid = fit$zero_tolerance,
    weight_specs = list(gaussian = list(weight_method = "gaussian_gps", gps_model = "spline", gps_df = 3))
  )

  expect_s3_class(scan, "cdmc_sensitivity_scan")
  expect_true(all(c("current", "gaussian") %in% scan$results$weight_scenario))
  expect_true(all(c("kernel_gps", "gaussian_gps") %in% scan$results$weight_mode))
  expect_true("delta_mean_tau_dr" %in% names(scan$results))
})

test_that("sensitivity scan can compare linear and GAM Gaussian GPS nuisance models", {
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
    seed = 1768
  )
  original_dose <- panel$dose
  active <- abs(original_dose) > 0
  panel$dose[active] <- sign(original_dose[active]) * (
    0.2 + sin(2 * panel$x1[active]) + (panel$time[active] / max(panel$time))^2
  )
  panel$y <- panel$y + 0.6 * (panel$dose - original_dose)

  fit <- cdmc_dr_fit(
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
    seed = 1768
  )

  scan <- cdmc_sensitivity_scan(
    fit,
    washout_grid = fit$washout,
    zero_tolerance_grid = fit$zero_tolerance,
    weight_specs = list(gam = list(weight_method = "gaussian_gps", gps_model = "gam", gps_df = 4))
  )

  expect_s3_class(scan, "cdmc_sensitivity_scan")
  expect_true(all(c("current", "gam") %in% scan$results$weight_scenario))
  expect_true(all(scan$results$weight_mode == "gaussian_gps"))
  expect_true(any(scan$results$weight_scenario == "gam" & scan$results$fit_success))
  expect_true("delta_mean_tau_dr" %in% names(scan$results))
})

test_that("sensitivity scan can compare linear and tree Gaussian GPS nuisance models", {
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
    seed = 1779
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

  fit <- cdmc_dr_fit(
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
    seed = 1779
  )

  scan <- cdmc_sensitivity_scan(
    fit,
    washout_grid = fit$washout,
    zero_tolerance_grid = fit$zero_tolerance,
    weight_specs = list(tree = list(weight_method = "gaussian_gps", gps_model = "tree"))
  )

  expect_s3_class(scan, "cdmc_sensitivity_scan")
  expect_true(all(c("current", "tree") %in% scan$results$weight_scenario))
  expect_true(all(scan$results$weight_mode == "gaussian_gps"))
  expect_true(any(scan$results$weight_scenario == "tree" & scan$results$fit_success))
  expect_true("delta_mean_tau_dr" %in% names(scan$results))
})

test_that("sensitivity scan can compare linear and forest Gaussian GPS nuisance models", {
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
    seed = 1786
  )
  original_dose <- panel$dose
  active <- abs(original_dose) > 0
  interaction_score <- panel$x2 + panel$time / max(panel$time)
  panel$dose[active] <- sign(original_dose[active]) * (
    0.15 +
      ifelse(panel$x1[active] > stats::median(panel$x1), 0.9, 0.2) +
      ifelse(interaction_score[active] > stats::median(interaction_score), 0.7, -0.1)
  )
  panel$y <- panel$y + 0.5 * (panel$dose - original_dose)

  fit <- cdmc_dr_fit(
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
    seed = 1786
  )

  scan <- cdmc_sensitivity_scan(
    fit,
    washout_grid = fit$washout,
    zero_tolerance_grid = fit$zero_tolerance,
    weight_specs = list(
      forest = list(
        weight_method = "gaussian_gps",
        gps_model = "forest",
        gps_forest_trees = 50L,
        gps_forest_mtry = 2L,
        gps_forest_min_node_size = 3L
      )
    )
  )

  expect_s3_class(scan, "cdmc_sensitivity_scan")
  expect_true(all(c("current", "forest") %in% scan$results$weight_scenario))
  expect_true(all(scan$results$weight_mode == "gaussian_gps"))
  expect_true(any(scan$results$weight_scenario == "forest" & scan$results$fit_success))
  expect_true("delta_mean_tau_dr" %in% names(scan$results))
})

test_that("sensitivity scan can replay stacked Gaussian GPS nuisance models", {
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
    seed = 1788
  )
  original_dose <- panel$dose
  active <- abs(original_dose) > 0
  interaction_score <- panel$x2 + panel$time / max(panel$time)
  panel$dose[active] <- sign(original_dose[active]) * (
    0.15 + panel$x1[active]^2 + ifelse(interaction_score[active] > stats::median(interaction_score), 0.5, -0.2)
  )
  panel$y <- panel$y + 0.5 * (panel$dose - original_dose)

  fit <- cdmc_dr_fit(
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
    seed = 1788
  )

  scan <- cdmc_sensitivity_scan(
    fit,
    washout_grid = fit$washout,
    zero_tolerance_grid = fit$zero_tolerance,
    weight_specs = list(
      stack = list(
        weight_method = "gaussian_gps",
        gps_model = "stack",
        gps_stack_models = c("linear", "spline", "tree")
      )
    )
  )

  expect_s3_class(scan, "cdmc_sensitivity_scan")
  expect_true(all(c("current", "stack") %in% scan$results$weight_scenario))
  expect_true(all(scan$results$weight_mode == "gaussian_gps"))
  expect_true(any(scan$results$weight_scenario == "stack" & scan$results$fit_success))
  expect_true("delta_mean_tau_dr" %in% names(scan$results))
})

test_that("sensitivity scan can replay boosted Gaussian GPS nuisance models", {
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
    seed = 1791
  )
  original_dose <- panel$dose
  active <- abs(original_dose) > 0
  interaction_score <- panel$x2 + panel$time / max(panel$time)
  panel$dose[active] <- sign(original_dose[active]) * (
    0.15 + panel$x1[active]^2 + 0.4 * interaction_score[active]
  )
  panel$y <- panel$y + 0.5 * (panel$dose - original_dose)

  fit <- cdmc_dr_fit(
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
    seed = 1791
  )

  scan <- cdmc_sensitivity_scan(
    fit,
    washout_grid = fit$washout,
    zero_tolerance_grid = fit$zero_tolerance,
    weight_specs = list(
      boost = list(
        weight_method = "gaussian_gps",
        gps_model = "boost",
        gps_boost_trees = 40L,
        gps_boost_depth = 2L,
        gps_boost_shrinkage = 0.05,
        gps_boost_min_obs_node = 3L
      )
    )
  )

  expect_s3_class(scan, "cdmc_sensitivity_scan")
  expect_true(all(c("current", "boost") %in% scan$results$weight_scenario))
  expect_true(all(scan$results$weight_mode == "gaussian_gps"))
  expect_true(any(scan$results$weight_scenario == "boost" & scan$results$fit_success))
  expect_true("delta_mean_tau_dr" %in% names(scan$results))
})

test_that("sensitivity scan can replay internal CBPS weighting", {
  testthat::skip_if_not_installed("CBPS")

  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 8,
    rank = 2,
    beta = 0.8,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.18,
    switch_off_prob = 0.42,
    seed = 1778
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
    lag_order = 0,
    seed = 1778
  )

  scan <- cdmc_sensitivity_scan(
    fit,
    washout_grid = fit$washout,
    zero_tolerance_grid = fit$zero_tolerance,
    weight_specs = list(cbps = list(weight_method = "cbps", cbps_iterations = 250L))
  )

  expect_s3_class(scan, "cdmc_sensitivity_scan")
  expect_true(all(c("current", "cbps") %in% scan$results$weight_scenario))
  expect_true(all(c("gaussian_gps", "cbps") %in% scan$results$weight_mode))
  expect_true("delta_mean_tau_dr" %in% names(scan$results))
})

test_that("sensitivity scan can replay internal entropy-balance weighting", {
  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 8,
    rank = 2,
    beta = 0.8,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.18,
    switch_off_prob = 0.42,
    seed = 1782
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
    lag_order = 0,
    seed = 1782
  )

  scan <- cdmc_sensitivity_scan(
    fit,
    washout_grid = fit$washout,
    zero_tolerance_grid = fit$zero_tolerance,
    weight_specs = list(entropy_balance = list(
      weight_method = "entropy_balance",
      entropy_balance_iterations = 500L
    ))
  )

  expect_s3_class(scan, "cdmc_sensitivity_scan")
  expect_true(all(c("current", "entropy_balance") %in% scan$results$weight_scenario))
  expect_true(all(c("gaussian_gps", "entropy_balance") %in% scan$results$weight_mode))
  expect_true("delta_mean_tau_dr" %in% names(scan$results))
})

test_that("sensitivity scan can replay internal kernel-balance weighting", {
  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 8,
    rank = 2,
    beta = 0.8,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.18,
    switch_off_prob = 0.42,
    seed = 1788
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
    lag_order = 0,
    seed = 1788
  )

  scan <- cdmc_sensitivity_scan(
    fit,
    washout_grid = fit$washout,
    zero_tolerance_grid = fit$zero_tolerance,
    weight_specs = list(kernel_balance = list(
      weight_method = "kernel_balance",
      kernel_balance_centers = 6L,
      kernel_balance_iterations = 400L
    ))
  )

  expect_s3_class(scan, "cdmc_sensitivity_scan")
  expect_true(all(c("current", "kernel_balance") %in% scan$results$weight_scenario))
  expect_true(all(c("gaussian_gps", "kernel_balance") %in% scan$results$weight_mode))
  expect_true("delta_mean_tau_dr" %in% names(scan$results))
})

test_that("sensitivity scan can replay internal adaptive-balance weighting", {
  panel <- simulate_cdmc_data(
    n_units = 20,
    n_times = 10,
    rank = 2,
    beta = 0.8,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.12,
    switch_off_prob = 0.5,
    seed = 1798
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
    lag_order = 0,
    seed = 1798
  )

  scan <- cdmc_sensitivity_scan(
    fit,
    washout_grid = fit$washout,
    zero_tolerance_grid = fit$zero_tolerance,
    weight_specs = list(adaptive_balance = list(
      weight_method = "adaptive_balance",
      adaptive_balance_methods = c("entropy_balance", "kernel_balance"),
      kernel_balance_centers = 6L,
      entropy_balance_iterations = 400L,
      kernel_balance_iterations = 400L
    ))
  )

  expect_s3_class(scan, "cdmc_sensitivity_scan")
  expect_true(all(c("current", "adaptive_balance") %in% scan$results$weight_scenario))
  expect_true(all(c("gaussian_gps", "adaptive_balance") %in% scan$results$weight_mode))
  expect_true("delta_mean_tau_dr" %in% names(scan$results))
})

test_that("sensitivity scan supports dynamic estimands built from dose-response objects", {
  panel <- simulate_cdmc_data(
    n_units = 12,
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
    seed = 1818
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
    seed = 1818
  )
  dose_response <- cdmc_dose_response(
    fit,
    model = "spline",
    lag_order = 0,
    df = 3,
    include_zero_dose = TRUE
  )
  dynamic <- cdmc_dynamic_estimand(
    dose_response,
    history = data.frame(dose_lag0 = c(-0.5, 0.5)),
    type = "slope",
    aggregate = "both"
  )

  scan <- cdmc_sensitivity_scan(
    dynamic,
    washout_grid = 0:1,
    zero_tolerance_grid = fit$zero_tolerance
  )

  expect_s3_class(scan, "cdmc_sensitivity_scan")
  expect_equal(scan$reference$target_class, "cdmc_dynamic_estimand")
  expect_true(all(c(
    "dynamic_slope_path1",
    "dynamic_slope_path2",
    "mean_dynamic_slope",
    "dynamic_path_count",
    "delta_mean_dynamic_slope"
  ) %in% names(scan$results)))
  expect_equal(nrow(scan$support_summary), 2)
  expect_s3_class(summary(scan), "summary.cdmc_sensitivity_scan")
})

test_that("formal sensitivity bounds match staged coefficient and average-tau maps", {
  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 10,
    rank = 2,
    beta = 1,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 1919
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
    seed = 1919
  )

  bounds <- cdmc_sensitivity_bounds(
    fit,
    statistics = c("coefficients", "average_tau"),
    gamma_grid = c(0, 0.5),
    scale = "manual",
    scale_value = 1
  )

  operator <- qr.coef(
    qr(fit$effect$design_info$design),
    diag(1, nrow = nrow(fit$effect$design_info$design), ncol = nrow(fit$effect$design_info$design))
  )
  operator <- as.matrix(operator)
  operator[is.na(operator)] <- 0

  coef_row <- bounds$bounds[
    bounds$bounds$statistic == "coef_dose_lag0" & bounds$bounds$gamma == 0.5,
    ,
    drop = FALSE
  ]
  coef_multiplier <- sum(abs(operator["dose_lag0", ]))

  mean_row <- bounds$bounds[
    bounds$bounds$statistic == "mean_tau_active_dose" & bounds$bounds$gamma == 0.5,
    ,
    drop = FALSE
  ]
  mean_reference <- bounds$reference[
    bounds$reference$statistic == "mean_tau_active_dose",
    ,
    drop = FALSE
  ]

  expect_s3_class(bounds, "cdmc_sensitivity_bounds")
  expect_equal(coef_row$radius, 0.5 * coef_multiplier, tolerance = 1e-8)
  expect_equal(
    coef_row$lower,
    unname(fit$effect$coefficients[["dose_lag0"]]) - 0.5 * coef_multiplier,
    tolerance = 1e-8
  )
  expect_equal(mean_row$radius, 0.5, tolerance = 1e-8)
  expect_equal(mean_reference$sensitivity_multiplier, 1, tolerance = 1e-8)
  expect_equal(mean_reference$breakdown_gamma, abs(mean_reference$estimate), tolerance = 1e-8)
  expect_s3_class(summary(bounds), "summary.cdmc_sensitivity_bounds")
})

test_that("formal sensitivity bounds support sparse contamination budgets", {
  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 10,
    rank = 2,
    beta = 1,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 1927
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
    seed = 1927
  )

  bounds <- cdmc_sensitivity_bounds(
    fit,
    statistics = "coefficients",
    gamma_grid = 0.5,
    contamination_fraction_grid = c(0.1, 1),
    scale = "manual",
    scale_value = 1
  )

  operator <- qr.coef(
    qr(fit$effect$design_info$design),
    diag(1, nrow = nrow(fit$effect$design_info$design), ncol = nrow(fit$effect$design_info$design))
  )
  operator <- as.matrix(operator)
  operator[is.na(operator)] <- 0
  loading <- abs(operator["dose_lag0", ])
  sparse_cells <- max(1L, ceiling(0.1 * length(loading)))
  sparse_multiplier <- sum(sort(loading, decreasing = TRUE)[seq_len(sparse_cells)])
  dense_multiplier <- sum(loading)

  sparse_row <- bounds$bounds[
    bounds$bounds$statistic == "coef_dose_lag0" &
      abs(bounds$bounds$contamination_fraction - 0.1) < 1e-12,
    ,
    drop = FALSE
  ]
  dense_row <- bounds$bounds[
    bounds$bounds$statistic == "coef_dose_lag0" &
      abs(bounds$bounds$contamination_fraction - 1) < 1e-12,
    ,
    drop = FALSE
  ]

  expect_s3_class(bounds, "cdmc_sensitivity_bounds")
  expect_equal(sparse_row$contamination_cells, sparse_cells)
  expect_equal(sparse_row$effective_multiplier, sparse_multiplier, tolerance = 1e-8)
  expect_equal(sparse_row$radius, 0.5 * sparse_multiplier, tolerance = 1e-8)
  expect_equal(dense_row$effective_multiplier, dense_multiplier, tolerance = 1e-8)
  expect_lt(sparse_row$radius, dense_row$radius)
})

test_that("formal sensitivity bounds support energy-budget perturbations", {
  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 10,
    rank = 2,
    beta = 1,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 1933
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
    seed = 1933
  )

  bounds <- cdmc_sensitivity_bounds(
    fit,
    statistics = "coefficients",
    gamma_grid = 0.5,
    perturbation_constraint = "energy",
    scale = "manual",
    scale_value = 1
  )

  operator <- qr.coef(
    qr(fit$effect$design_info$design),
    diag(1, nrow = nrow(fit$effect$design_info$design), ncol = nrow(fit$effect$design_info$design))
  )
  operator <- as.matrix(operator)
  operator[is.na(operator)] <- 0
  energy_multiplier <- sqrt(sum(operator["dose_lag0", ] ^ 2))

  energy_row <- bounds$bounds[
    bounds$bounds$statistic == "coef_dose_lag0" & bounds$bounds$gamma == 0.5,
    ,
    drop = FALSE
  ]

  expect_s3_class(bounds, "cdmc_sensitivity_bounds")
  expect_identical(bounds$perturbation_constraint, "energy")
  expect_equal(energy_row$effective_multiplier, energy_multiplier, tolerance = 1e-8)
  expect_equal(energy_row$radius, 0.5 * energy_multiplier, tolerance = 1e-8)
})

test_that("formal sensitivity bounds support unit- and time-scope perturbations", {
  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 10,
    rank = 2,
    beta = 1,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 1937
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
    seed = 1937
  )

  bounds_unit <- cdmc_sensitivity_bounds(
    fit,
    statistics = "coefficients",
    gamma_grid = 0.5,
    perturbation_scope = "unit",
    scale = "manual",
    scale_value = 1
  )
  bounds_time <- cdmc_sensitivity_bounds(
    fit,
    statistics = "coefficients",
    gamma_grid = 0.5,
    perturbation_scope = "time",
    scale = "manual",
    scale_value = 1
  )

  operator <- qr.coef(
    qr(fit$effect$design_info$design),
    diag(1, nrow = nrow(fit$effect$design_info$design), ncol = nrow(fit$effect$design_info$design))
  )
  operator <- as.matrix(operator)
  operator[is.na(operator)] <- 0
  sample_indices <- which(fit$effect$sample_mask, arr.ind = TRUE)
  unit_ids <- as.character(fit$unit_levels[sample_indices[, 1L]])
  time_ids <- as.character(fit$time_levels[sample_indices[, 2L]])
  unit_multiplier <- sum(abs(tapply(operator["dose_lag0", ], unit_ids, sum)))
  time_multiplier <- sum(abs(tapply(operator["dose_lag0", ], time_ids, sum)))

  unit_row <- bounds_unit$bounds[
    bounds_unit$bounds$statistic == "coef_dose_lag0" & bounds_unit$bounds$gamma == 0.5,
    ,
    drop = FALSE
  ]
  time_row <- bounds_time$bounds[
    bounds_time$bounds$statistic == "coef_dose_lag0" & bounds_time$bounds$gamma == 0.5,
    ,
    drop = FALSE
  ]

  expect_s3_class(bounds_unit, "cdmc_sensitivity_bounds")
  expect_identical(bounds_unit$perturbation_scope, "unit")
  expect_equal(unit_row$effective_multiplier, unit_multiplier, tolerance = 1e-8)
  expect_equal(unit_row$radius, 0.5 * unit_multiplier, tolerance = 1e-8)

  expect_s3_class(bounds_time, "cdmc_sensitivity_bounds")
  expect_identical(bounds_time$perturbation_scope, "time")
  expect_equal(time_row$effective_multiplier, time_multiplier, tolerance = 1e-8)
  expect_equal(time_row$radius, 0.5 * time_multiplier, tolerance = 1e-8)
})

test_that("formal sensitivity bounds support dose-response predictions", {
  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 10,
    rank = 2,
    beta = 1,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 1919
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
    seed = 1919
  )
  dose_response <- cdmc_dose_response(
    fit,
    model = "linear",
    lag_order = 1,
    include_zero_dose = TRUE
  )
  history <- data.frame(
    dose_lag0 = c(1, 0.5),
    dose_lag1 = c(0.5, 1)
  )

  bounds <- cdmc_sensitivity_bounds(
    dose_response,
    statistics = "prediction",
    prediction_history = history,
    gamma_grid = c(0, 0.25),
    scale = "manual",
    scale_value = 2
  )

  operator <- qr.coef(
    qr(dose_response$design_info$design),
    diag(1, nrow = nrow(dose_response$design_info$design), ncol = nrow(dose_response$design_info$design))
  )
  operator <- as.matrix(operator)
  operator[is.na(operator)] <- 0

  prediction_rows <- as.matrix(history[, c("dose_lag0", "dose_lag1"), drop = FALSE])
  loading_one <- prediction_rows[1, , drop = FALSE] %*% operator
  response_row <- bounds$bounds[
    bounds$bounds$statistic == "dose_response_response_1" & bounds$bounds$gamma == 0.25,
    ,
    drop = FALSE
  ]

  expect_s3_class(bounds, "cdmc_sensitivity_bounds")
  expect_equal(response_row$radius, 0.5 * sum(abs(loading_one)), tolerance = 1e-8)
})

test_that("formal sensitivity bounds reject nonlinear gam dose-response fits", {
  testthat::skip_if_not_installed("mgcv")

  panel <- simulate_cdmc_data(
    n_units = 12,
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
    seed = 1941
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
    seed = 1941
  )
  dose_response <- cdmc_dose_response(
    fit,
    model = "gam",
    lag_order = 0,
    df = 4,
    include_zero_dose = TRUE
  )

  expect_error(
    cdmc_sensitivity_bounds(
      dose_response,
      statistics = "prediction",
      prediction_dose = c(-0.5, 0.5),
      gamma_grid = c(0, 0.25),
      scale = "manual",
      scale_value = 1
    ),
    "Use perturbation_layer = 'optimization'"
  )
})

test_that("formal sensitivity bounds reject nonlinear tree dose-response fits", {
  testthat::skip_if_not_installed("rpart")

  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 10,
    rank = 2,
    beta = 0,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    signed_dose = TRUE,
    negative_dose_prob = 0.5,
    seed = 1947
  )
  panel$y <- panel$y + ifelse(panel$dose > 0.25, 0.8, 0)

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
    seed = 1947
  )
  dose_response <- cdmc_dose_response(
    fit,
    model = "tree",
    lag_order = 0,
    include_zero_dose = TRUE
  )

  expect_error(
    cdmc_sensitivity_bounds(
      dose_response,
      statistics = "prediction",
      prediction_dose = c(0, 0.8),
      gamma_grid = c(0, 0.25),
      scale = "manual",
      scale_value = 1
    ),
    "Use perturbation_layer = 'optimization'"
  )
})

test_that("formal sensitivity bounds reject nonlinear forest dose-response fits", {
  testthat::skip_if_not_installed("ranger")

  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 10,
    rank = 2,
    beta = 0,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    signed_dose = TRUE,
    negative_dose_prob = 0.5,
    seed = 1949
  )
  panel$y <- panel$y + ifelse(panel$dose > 0.25, 0.8, 0)

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
    seed = 1949
  )
  dose_response <- cdmc_dose_response(
    fit,
    model = "forest",
    lag_order = 0,
    forest_trees = 50,
    include_zero_dose = TRUE
  )

  expect_error(
    cdmc_sensitivity_bounds(
      dose_response,
      statistics = "prediction",
      prediction_dose = c(0, 0.8),
      gamma_grid = c(0, 0.25),
      scale = "manual",
      scale_value = 1
    ),
    "Use perturbation_layer = 'optimization'"
  )
})

test_that("formal sensitivity bounds support optimization-layer nonlinear dose-response predictions", {
  model_specs <- list()
  if (requireNamespace("mgcv", quietly = TRUE)) {
    model_specs$gam <- list(model = "gam", df = 4)
  }
  if (requireNamespace("rpart", quietly = TRUE)) {
    model_specs$tree <- list(model = "tree")
  }
  if (requireNamespace("ranger", quietly = TRUE)) {
    model_specs$forest <- list(model = "forest", forest_trees = 50)
  }
  if (length(model_specs) < 1L) {
    testthat::skip("requires mgcv, rpart, or ranger")
  }

  panel <- simulate_cdmc_data(
    n_units = 5,
    n_times = 6,
    rank = 2,
    beta = 0.2,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.03,
    switch_on_prob = 0.18,
    switch_off_prob = 0.45,
    signed_dose = TRUE,
    negative_dose_prob = 0.5,
    seed = 2651
  )
  panel$y <- panel$y + 0.8 * panel$dose^2

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = 0.15,
    rank_max = 2,
    washout = 0,
    lag_order = 0,
    seed = 2651
  )
  history <- data.frame(dose_lag0 = c(-0.5, 0.5))
  step <- 1e-3

  for (model_name in names(model_specs)) {
    args <- c(
      list(
        object = fit,
        lag_order = 0,
        include_zero_dose = TRUE
      ),
      model_specs[[model_name]]
    )
    dose_response <- do.call(cdmc_dose_response, args)

    expected <- manual_sensitivity_optimization_radius(
      dose_response,
      statistic_name = "dose_response_response_1",
      step = step,
      prediction_history = history
    )

    bounds <- cdmc_sensitivity_bounds(
      dose_response,
      statistics = "prediction",
      prediction_history = history,
      gamma_grid = c(0, 0.5),
      perturbation_layer = "optimization",
      scale = "manual",
      scale_value = 1,
      refit_step = step
    )

    response_row <- bounds$bounds[
      bounds$bounds$statistic == "dose_response_response_1" & bounds$bounds$gamma == 0.5,
      ,
      drop = FALSE
    ]

    expect_s3_class(bounds, "cdmc_sensitivity_bounds")
    expect_identical(bounds$perturbation_layer, "optimization")
    expect_identical(bounds$bound_type, "local_refit")
    expect_equal(bounds$n_perturbation_cells, expected$n_cells)
    expect_equal(response_row$radius, 0.5 * expected$multiplier, tolerance = 1e-6)
  }
})

test_that("formal sensitivity bounds support optimization-layer dynamic estimands backed by nonlinear dose-response models", {
  testthat::skip_if_not_installed("mgcv")

  panel <- simulate_cdmc_data(
    n_units = 5,
    n_times = 6,
    rank = 2,
    beta = 0.2,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.03,
    switch_on_prob = 0.18,
    switch_off_prob = 0.45,
    signed_dose = TRUE,
    negative_dose_prob = 0.5,
    seed = 2653
  )
  panel$y <- panel$y + 0.8 * panel$dose^2

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = 0.15,
    rank_max = 2,
    washout = 0,
    lag_order = 0,
    seed = 2653
  )
  dose_response <- cdmc_dose_response(
    fit,
    model = "gam",
    lag_order = 0,
    df = 4,
    include_zero_dose = TRUE
  )
  dynamic <- cdmc_dynamic_estimand(
    dose_response,
    history = data.frame(dose_lag0 = c(-0.5, 0.5)),
    type = "response",
    aggregate = "mean"
  )

  step <- 1e-3
  expected <- manual_sensitivity_optimization_radius(
    dynamic,
    statistic_name = "mean_dynamic_response",
    step = step
  )

  bounds <- cdmc_sensitivity_bounds(
    dynamic,
    gamma_grid = c(0, 0.5),
    perturbation_layer = "optimization",
    scale = "manual",
    scale_value = 1,
    refit_step = step
  )

  mean_row <- bounds$bounds[
    bounds$bounds$statistic == "mean_dynamic_response" & bounds$bounds$gamma == 0.5,
    ,
    drop = FALSE
  ]

  expect_s3_class(bounds, "cdmc_sensitivity_bounds")
  expect_identical(bounds$perturbation_layer, "optimization")
  expect_identical(bounds$bound_type, "local_refit")
  expect_equal(bounds$n_perturbation_cells, expected$n_cells)
  expect_equal(mean_row$radius, 0.5 * expected$multiplier, tolerance = 1e-6)
})

test_that("formal sensitivity bounds support DR dynamic estimands", {
  panel <- simulate_cdmc_data(
    n_units = 16,
    n_times = 10,
    rank = 2,
    beta = 0.9,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.18,
    switch_off_prob = 0.42,
    seed = 2121
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
    seed = 2121
  )
  history <- data.frame(
    dose_lag0 = c(1, 0.25),
    dose_lag1 = c(0.5, 1)
  )
  dynamic <- cdmc_dynamic_estimand(
    fit,
    history = history,
    type = "contrast",
    aggregate = "both"
  )

  bounds <- cdmc_sensitivity_bounds(
    dynamic,
    gamma_grid = c(0, 0.25),
    scale = "manual",
    scale_value = 2
  )

  sample_indices <- which(fit$effect$sample_mask, arr.ind = TRUE)
  design <- cbind(
    dose_lag0 = fit$dose_matrix[sample_indices],
    dose_lag1 = fit$dose_matrix[cbind(sample_indices[, 1], sample_indices[, 2] - 1L)]
  )
  operator <- qr.coef(
    qr(design),
    diag(1, nrow = nrow(design), ncol = nrow(design))
  )
  operator <- as.matrix(operator)
  operator[is.na(operator)] <- 0

  mean_loading <- colMeans(as.matrix(history)) %*% operator
  mean_row <- bounds$bounds[
    bounds$bounds$statistic == "mean_dynamic_contrast" & bounds$bounds$gamma == 0.25,
    ,
    drop = FALSE
  ]

  expect_s3_class(bounds, "cdmc_sensitivity_bounds")
  expect_true(all(c(
    "dynamic_contrast_path1",
    "dynamic_contrast_path2",
    "mean_dynamic_contrast"
  ) %in% bounds$reference$statistic))
  expect_equal(mean_row$radius, 0.5 * sum(abs(mean_loading)), tolerance = 1e-8)
})

test_that("formal sensitivity bounds support direct joint-objective main fits", {
  panel <- simulate_cdmc_data(
    n_units = 10,
    n_times = 8,
    rank = 2,
    beta = 1,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 2222
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
    lag_order = 0,
    objective = "joint",
    seed = 2222
  )

  bounds <- cdmc_sensitivity_bounds(
    fit,
    statistics = c("coefficients", "average_tau"),
    gamma_grid = c(0, 0.25),
    scale = "manual",
    scale_value = 1
  )

  operator <- qr.coef(
    qr(fit$effect$design_info$design),
    diag(1, nrow = nrow(fit$effect$design_info$design), ncol = nrow(fit$effect$design_info$design))
  )
  operator <- as.matrix(operator)
  operator[is.na(operator)] <- 0

  coef_row <- bounds$bounds[
    bounds$bounds$statistic == "coef_dose_lag0" & bounds$bounds$gamma == 0.25,
    ,
    drop = FALSE
  ]
  mean_row <- bounds$bounds[
    bounds$bounds$statistic == "mean_tau_active_dose" & bounds$bounds$gamma == 0.25,
    ,
    drop = FALSE
  ]

  expect_s3_class(bounds, "cdmc_sensitivity_bounds")
  expect_equal(bounds$perturbation_layer, "stage2")
  expect_equal(coef_row$radius, 0.25 * sum(abs(operator["dose_lag0", ])), tolerance = 1e-8)
  expect_equal(mean_row$radius, 0.25, tolerance = 1e-8)
})

test_that("formal sensitivity bounds support dynamic estimands from joint fits", {
  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 10,
    rank = 2,
    beta = 1,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 2323
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
    objective = "joint",
    seed = 2323
  )
  history <- data.frame(
    dose_lag0 = c(1, 0.25),
    dose_lag1 = c(0.5, 1)
  )
  dynamic <- cdmc_dynamic_estimand(
    fit,
    history = history,
    type = "contrast",
    aggregate = "both"
  )

  bounds <- cdmc_sensitivity_bounds(
    dynamic,
    gamma_grid = c(0, 0.25),
    scale = "manual",
    scale_value = 2
  )

  operator <- qr.coef(
    qr(fit$effect$design_info$design),
    diag(1, nrow = nrow(fit$effect$design_info$design), ncol = nrow(fit$effect$design_info$design))
  )
  operator <- as.matrix(operator)
  operator[is.na(operator)] <- 0

  mean_loading <- colMeans(as.matrix(history)) %*% operator
  mean_row <- bounds$bounds[
    bounds$bounds$statistic == "mean_dynamic_contrast" & bounds$bounds$gamma == 0.25,
    ,
    drop = FALSE
  ]

  expect_s3_class(bounds, "cdmc_sensitivity_bounds")
  expect_true(all(c(
    "dynamic_contrast_path1",
    "dynamic_contrast_path2",
    "mean_dynamic_contrast"
  ) %in% bounds$reference$statistic))
  expect_equal(mean_row$radius, 0.5 * sum(abs(mean_loading)), tolerance = 1e-8)
})

test_that("formal sensitivity bounds support baseline-layer DR violations", {
  panel <- simulate_cdmc_data(
    n_units = 16,
    n_times = 10,
    rank = 2,
    beta = 0.9,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.18,
    switch_off_prob = 0.42,
    seed = 2424
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
    seed = 2424
  )

  stage2_bounds <- cdmc_sensitivity_bounds(
    fit,
    statistics = "coefficients",
    gamma_grid = c(0, 0.25),
    scale = "manual",
    scale_value = 1
  )
  baseline_bounds <- cdmc_sensitivity_bounds(
    fit,
    statistics = "coefficients",
    gamma_grid = c(0, 0.25),
    perturbation_layer = "baseline",
    scale = "manual",
    scale_value = 1
  )

  sample_indices <- which(fit$effect$sample_mask, arr.ind = TRUE)
  design <- cbind(
    dose_lag0 = fit$dose_matrix[sample_indices],
    dose_lag1 = fit$dose_matrix[cbind(sample_indices[, 1], sample_indices[, 2] - 1L)]
  )
  operator <- qr.coef(
    qr(design),
    diag(1, nrow = nrow(design), ncol = nrow(design))
  )
  operator <- as.matrix(operator)
  operator[is.na(operator)] <- 0
  sample_weights <- as.numeric(fit$weight_matrix[fit$effect$sample_mask])
  weighted_operator <- sweep(operator, 2, sample_weights, FUN = "*")

  stage2_row <- stage2_bounds$bounds[
    stage2_bounds$bounds$statistic == "coef_dose_lag0" & stage2_bounds$bounds$gamma == 0.25,
    ,
    drop = FALSE
  ]
  baseline_row <- baseline_bounds$bounds[
    baseline_bounds$bounds$statistic == "coef_dose_lag0" & baseline_bounds$bounds$gamma == 0.25,
    ,
    drop = FALSE
  ]

  expect_s3_class(baseline_bounds, "cdmc_sensitivity_bounds")
  expect_equal(baseline_bounds$perturbation_layer, "baseline")
  expect_gt(max(abs(sample_weights - 1)), 1e-6)
  expect_equal(stage2_row$radius, 0.25 * sum(abs(operator["dose_lag0", ])), tolerance = 1e-8)
  expect_equal(baseline_row$radius, 0.25 * sum(abs(weighted_operator["dose_lag0", ])), tolerance = 1e-8)
  expect_false(isTRUE(all.equal(stage2_row$radius, baseline_row$radius, tolerance = 1e-8)))
})

test_that("formal sensitivity bounds support optimization-layer local refits for joint main fits", {
  panel <- simulate_cdmc_data(
    n_units = 6,
    n_times = 6,
    rank = 2,
    beta = 1,
    lag_beta = 0.15,
    n_covariates = 1,
    noise_sd = 0.03,
    switch_on_prob = 0.18,
    switch_off_prob = 0.45,
    seed = 2525
  )

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = 0.15,
    rank_max = 2,
    washout = 0,
    lag_order = 0,
    objective = "joint",
    seed = 2525
  )

  step <- 1e-3
  expected <- manual_sensitivity_optimization_radius(
    fit,
    statistic_name = "coef_dose_lag0",
    step = step
  )

  bounds <- cdmc_sensitivity_bounds(
    fit,
    statistics = "coefficients",
    gamma_grid = c(0, 0.25),
    perturbation_layer = "optimization",
    scale = "manual",
    scale_value = 1,
    refit_step = step
  )

  coef_row <- bounds$bounds[
    bounds$bounds$statistic == "coef_dose_lag0" & bounds$bounds$gamma == 0.25,
    ,
    drop = FALSE
  ]

  expect_s3_class(bounds, "cdmc_sensitivity_bounds")
  expect_equal(bounds$perturbation_layer, "optimization")
  expect_equal(bounds$bound_type, "local_refit")
  expect_equal(bounds$refit_step, step)
  expect_equal(bounds$n_perturbation_cells, expected$n_cells)
  expect_equal(coef_row$radius, 0.25 * expected$multiplier, tolerance = 1e-6)
})

test_that("formal sensitivity bounds support optimization-layer dose-response predictions", {
  panel <- simulate_cdmc_data(
    n_units = 6,
    n_times = 7,
    rank = 2,
    beta = 1,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.03,
    switch_on_prob = 0.18,
    switch_off_prob = 0.45,
    seed = 2626
  )

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = 0.15,
    rank_max = 2,
    washout = 0,
    lag_order = 1,
    seed = 2626
  )
  dose_response <- cdmc_dose_response(
    fit,
    model = "linear",
    lag_order = 1,
    include_zero_dose = TRUE
  )
  history <- data.frame(
    dose_lag0 = c(1, 0.5),
    dose_lag1 = c(0.5, 1)
  )

  step <- 1e-3
  expected <- manual_sensitivity_optimization_radius(
    dose_response,
    statistic_name = "dose_response_response_1",
    step = step,
    prediction_history = history
  )

  bounds <- cdmc_sensitivity_bounds(
    dose_response,
    statistics = "prediction",
    prediction_history = history,
    gamma_grid = c(0, 0.5),
    perturbation_layer = "optimization",
    scale = "manual",
    scale_value = 1,
    refit_step = step
  )

  response_row <- bounds$bounds[
    bounds$bounds$statistic == "dose_response_response_1" & bounds$bounds$gamma == 0.5,
    ,
    drop = FALSE
  ]

  expect_s3_class(bounds, "cdmc_sensitivity_bounds")
  expect_equal(bounds$bound_type, "local_refit")
  expect_equal(bounds$n_perturbation_cells, expected$n_cells)
  expect_equal(response_row$radius, 0.5 * expected$multiplier, tolerance = 1e-6)
})

test_that("formal sensitivity bounds support optimization-layer DR coefficients", {
  panel <- simulate_cdmc_data(
    n_units = 5,
    n_times = 6,
    rank = 2,
    beta = 0.9,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.03,
    switch_on_prob = 0.18,
    switch_off_prob = 0.45,
    seed = 2727
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
    lambda = 0.15,
    rank_max = 2,
    washout = 0,
    lag_order = 1,
    seed = 2727
  )

  step <- 1e-3
  expected <- manual_sensitivity_optimization_radius(
    fit,
    statistic_name = "coef_dose_lag0",
    step = step
  )

  bounds <- cdmc_sensitivity_bounds(
    fit,
    statistics = "coefficients",
    gamma_grid = c(0, 0.25),
    perturbation_layer = "optimization",
    scale = "manual",
    scale_value = 1,
    refit_step = step
  )

  coef_row <- bounds$bounds[
    bounds$bounds$statistic == "coef_dose_lag0" & bounds$bounds$gamma == 0.25,
    ,
    drop = FALSE
  ]

  expect_s3_class(bounds, "cdmc_sensitivity_bounds")
  expect_equal(bounds$source_class, "cdmc_dr_fit")
  expect_equal(bounds$perturbation_layer, "optimization")
  expect_equal(bounds$bound_type, "local_refit")
  expect_equal(bounds$refit_step, step)
  expect_equal(bounds$n_perturbation_cells, expected$n_cells)
  expect_equal(coef_row$radius, 0.25 * expected$multiplier, tolerance = 1e-6)
})

test_that("formal sensitivity bounds support optimization-layer DR dynamic estimands", {
  panel <- simulate_cdmc_data(
    n_units = 5,
    n_times = 6,
    rank = 2,
    beta = 0.9,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.03,
    switch_on_prob = 0.18,
    switch_off_prob = 0.45,
    seed = 2828
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
    lambda = 0.15,
    rank_max = 2,
    washout = 0,
    lag_order = 1,
    seed = 2828
  )
  history <- data.frame(
    dose_lag0 = c(1, 0.25),
    dose_lag1 = c(0.5, 1)
  )
  dynamic <- cdmc_dynamic_estimand(
    fit,
    history = history,
    type = "contrast",
    aggregate = "both"
  )

  step <- 1e-3
  expected <- manual_sensitivity_optimization_radius(
    dynamic,
    statistic_name = "mean_dynamic_contrast",
    step = step
  )

  bounds <- cdmc_sensitivity_bounds(
    dynamic,
    gamma_grid = c(0, 0.25),
    perturbation_layer = "optimization",
    scale = "manual",
    scale_value = 1,
    refit_step = step
  )

  mean_row <- bounds$bounds[
    bounds$bounds$statistic == "mean_dynamic_contrast" & bounds$bounds$gamma == 0.25,
    ,
    drop = FALSE
  ]

  expect_s3_class(bounds, "cdmc_sensitivity_bounds")
  expect_equal(bounds$source_class, "cdmc_dr_fit")
  expect_equal(bounds$bound_type, "local_refit")
  expect_equal(bounds$n_perturbation_cells, expected$n_cells)
  expect_equal(mean_row$radius, 0.25 * expected$multiplier, tolerance = 1e-6)
})

test_that("formal sensitivity bounds support sparse optimization-layer contamination budgets", {
  panel <- simulate_cdmc_data(
    n_units = 8,
    n_times = 8,
    rank = 2,
    beta = 0.9,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.03,
    switch_on_prob = 0.16,
    switch_off_prob = 0.5,
    seed = 2836
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
    lambda = 0.15,
    rank_max = 3,
    washout = 0,
    lag_order = 1,
    seed = 2836
  )

  step <- 1e-3
  source_fit <- cdmc_sensitivity_source_object(fit)
  full_data <- source_fit$data
  observed_rows <- rep(TRUE, nrow(full_data))
  if (".cdmc_observed" %in% names(full_data)) {
    observed_rows <- !is.na(full_data$.cdmc_observed) & full_data$.cdmc_observed
  }
  optimization_mask <- cdmc_sensitivity_bounds_optimization_mask(source_fit)
  row_lookup <- integer(length(observed_rows))
  row_lookup[observed_rows] <- seq_len(sum(observed_rows))
  perturbation_rows <- row_lookup[observed_rows & as.vector(t(optimization_mask))]
  perturbation_rows <- perturbation_rows[perturbation_rows > 0L]
  fit_data <- full_data[observed_rows, cdmc_sensitivity_original_columns(fit), drop = FALSE]
  base_statistics <- manual_sensitivity_statistics(fit)
  base_outcome <- fit_data[[source_fit$outcome]]
  derivatives <- numeric(length(perturbation_rows))

  for (index in seq_along(perturbation_rows)) {
    perturbed_data <- fit_data
    perturbed_data[[source_fit$outcome]][perturbation_rows[[index]]] <-
      base_outcome[[perturbation_rows[[index]]]] + step
    refit_fit <- manual_sensitivity_refit_source_fit(source_fit, perturbed_data)
    refit_statistics <- manual_sensitivity_statistics(refit_fit)
    derivatives[[index]] <- (refit_statistics[["coef_dose_lag0"]] - base_statistics[["coef_dose_lag0"]]) / step
  }

  sparse_fraction <- 0.1
  sparse_cells <- max(1L, ceiling(sparse_fraction * length(derivatives)))
  sparse_multiplier <- sum(sort(abs(derivatives), decreasing = TRUE)[seq_len(sparse_cells)])

  bounds <- cdmc_sensitivity_bounds(
    fit,
    statistics = "coefficients",
    gamma_grid = 0.25,
    contamination_fraction_grid = c(sparse_fraction, 1),
    perturbation_layer = "optimization",
    scale = "manual",
    scale_value = 1,
    refit_step = step
  )

  sparse_row <- bounds$bounds[
    bounds$bounds$statistic == "coef_dose_lag0" &
      abs(bounds$bounds$contamination_fraction - sparse_fraction) < 1e-12,
    ,
    drop = FALSE
  ]
  dense_row <- bounds$bounds[
    bounds$bounds$statistic == "coef_dose_lag0" &
      abs(bounds$bounds$contamination_fraction - 1) < 1e-12,
    ,
    drop = FALSE
  ]

  expect_s3_class(bounds, "cdmc_sensitivity_bounds")
  expect_equal(sparse_row$contamination_cells, sparse_cells)
  expect_equal(sparse_row$effective_multiplier, sparse_multiplier, tolerance = 1e-6)
  expect_equal(sparse_row$radius, 0.25 * sparse_multiplier, tolerance = 1e-6)
  expect_lt(sparse_row$radius, dense_row$radius)
})

test_that("formal sensitivity bounds support optimization-layer energy budgets", {
  panel <- simulate_cdmc_data(
    n_units = 8,
    n_times = 8,
    rank = 2,
    beta = 0.9,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.03,
    switch_on_prob = 0.16,
    switch_off_prob = 0.5,
    seed = 2844
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
    lambda = 0.15,
    rank_max = 3,
    washout = 0,
    lag_order = 1,
    seed = 2844
  )

  step <- 1e-3
  source_fit <- cdmc_sensitivity_source_object(fit)
  full_data <- source_fit$data
  observed_rows <- rep(TRUE, nrow(full_data))
  if (".cdmc_observed" %in% names(full_data)) {
    observed_rows <- !is.na(full_data$.cdmc_observed) & full_data$.cdmc_observed
  }
  optimization_mask <- cdmc_sensitivity_bounds_optimization_mask(source_fit)
  row_lookup <- integer(length(observed_rows))
  row_lookup[observed_rows] <- seq_len(sum(observed_rows))
  perturbation_rows <- row_lookup[observed_rows & as.vector(t(optimization_mask))]
  perturbation_rows <- perturbation_rows[perturbation_rows > 0L]
  fit_data <- full_data[observed_rows, cdmc_sensitivity_original_columns(fit), drop = FALSE]
  base_statistics <- manual_sensitivity_statistics(fit)
  base_outcome <- fit_data[[source_fit$outcome]]
  derivatives <- numeric(length(perturbation_rows))

  for (index in seq_along(perturbation_rows)) {
    perturbed_data <- fit_data
    perturbed_data[[source_fit$outcome]][perturbation_rows[[index]]] <-
      base_outcome[[perturbation_rows[[index]]]] + step
    refit_fit <- manual_sensitivity_refit_source_fit(source_fit, perturbed_data)
    refit_statistics <- manual_sensitivity_statistics(refit_fit)
    derivatives[[index]] <- (refit_statistics[["coef_dose_lag0"]] - base_statistics[["coef_dose_lag0"]]) / step
  }

  energy_multiplier <- sqrt(sum(derivatives ^ 2))

  bounds <- cdmc_sensitivity_bounds(
    fit,
    statistics = "coefficients",
    gamma_grid = 0.25,
    perturbation_constraint = "energy",
    perturbation_layer = "optimization",
    scale = "manual",
    scale_value = 1,
    refit_step = step
  )

  energy_row <- bounds$bounds[
    bounds$bounds$statistic == "coef_dose_lag0" & bounds$bounds$gamma == 0.25,
    ,
    drop = FALSE
  ]

  expect_s3_class(bounds, "cdmc_sensitivity_bounds")
  expect_identical(bounds$perturbation_constraint, "energy")
  expect_equal(energy_row$effective_multiplier, energy_multiplier, tolerance = 1e-6)
  expect_equal(energy_row$radius, 0.25 * energy_multiplier, tolerance = 1e-6)
})

test_that("formal sensitivity bounds support optimization-layer unit-scope perturbations", {
  panel <- simulate_cdmc_data(
    n_units = 8,
    n_times = 8,
    rank = 2,
    beta = 0.9,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.03,
    switch_on_prob = 0.16,
    switch_off_prob = 0.5,
    seed = 2844
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
    lambda = 0.15,
    rank_max = 3,
    washout = 0,
    lag_order = 1,
    seed = 2844
  )

  step <- 1e-3
  source_fit <- cdmc_sensitivity_source_object(fit)
  full_data <- source_fit$data
  observed_rows <- rep(TRUE, nrow(full_data))
  if (".cdmc_observed" %in% names(full_data)) {
    observed_rows <- !is.na(full_data$.cdmc_observed) & full_data$.cdmc_observed
  }
  optimization_mask <- cdmc_sensitivity_bounds_optimization_mask(source_fit)
  row_lookup <- integer(length(observed_rows))
  row_lookup[observed_rows] <- seq_len(sum(observed_rows))
  perturbation_rows <- row_lookup[observed_rows & as.vector(t(optimization_mask))]
  perturbation_rows <- perturbation_rows[perturbation_rows > 0L]
  fit_data <- full_data[observed_rows, cdmc_sensitivity_original_columns(fit), drop = FALSE]
  base_statistics <- manual_sensitivity_statistics(fit)
  base_outcome <- fit_data[[source_fit$outcome]]
  derivatives <- numeric(length(perturbation_rows))

  for (index in seq_along(perturbation_rows)) {
    perturbed_data <- fit_data
    perturbed_data[[source_fit$outcome]][perturbation_rows[[index]]] <-
      base_outcome[[perturbation_rows[[index]]]] + step
    refit_fit <- manual_sensitivity_refit_source_fit(source_fit, perturbed_data)
    refit_statistics <- manual_sensitivity_statistics(refit_fit)
    derivatives[[index]] <- (refit_statistics[["coef_dose_lag0"]] - base_statistics[["coef_dose_lag0"]]) / step
  }

  unit_ids <- as.character(fit_data[[source_fit$unit]][perturbation_rows])
  unit_multiplier <- sum(abs(tapply(derivatives, unit_ids, sum)))

  bounds <- cdmc_sensitivity_bounds(
    fit,
    statistics = "coefficients",
    gamma_grid = 0.25,
    perturbation_scope = "unit",
    perturbation_layer = "optimization",
    scale = "manual",
    scale_value = 1,
    refit_step = step
  )

  unit_row <- bounds$bounds[
    bounds$bounds$statistic == "coef_dose_lag0" & bounds$bounds$gamma == 0.25,
    ,
    drop = FALSE
  ]

  expect_s3_class(bounds, "cdmc_sensitivity_bounds")
  expect_identical(bounds$perturbation_scope, "unit")
  expect_equal(unit_row$effective_multiplier, unit_multiplier, tolerance = 1e-6)
  expect_equal(unit_row$radius, 0.25 * unit_multiplier, tolerance = 1e-6)
})