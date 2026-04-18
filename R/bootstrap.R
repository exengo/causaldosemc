cdmc_bootstrap_sample_data <- function(data, unit, original_columns) {
  if (is.data.frame(data) && ".cdmc_observed" %in% names(data)) {
    data <- data[!is.na(data$.cdmc_observed) & data$.cdmc_observed, , drop = FALSE]
  }
  unit_ids <- unique(data[[unit]])
  if (length(unit_ids) == 0L) {
    return(data[, original_columns, drop = FALSE])
  }

  sampled_units <- as.character(sample(unit_ids, length(unit_ids), replace = TRUE))
  unit_chunks <- split(data[, original_columns, drop = FALSE], data[[unit]], drop = TRUE)
  bootstrap_chunks <- unname(unit_chunks[sampled_units])

  bootstrap_chunks <- lapply(seq_along(bootstrap_chunks), function(sample_index) {
    bootstrap_chunk <- bootstrap_chunks[[sample_index]]
    bootstrap_chunk[[unit]] <- paste0("boot_unit_", sample_index)
    bootstrap_chunk
  })

  do.call(rbind, bootstrap_chunks)
}

cdmc_supported_bootstrap_classes <- function() {
  c(
    "cdmc_fit",
    "cdmc_dynamic_estimand",
    "cdmc_dose_response",
    "cdmc_dr_fit",
    "cdmc_placebo_test",
    "cdmc_carryover_test",
    "cdmc_carryover_refit_test",
    "cdmc_joint_placebo_test",
    "cdmc_equivalence_test",
    "cdmc_scia_test"
  )
}

cdmc_bootstrap_object_class <- function(object) {
  class_name <- cdmc_first_inherited_class(object, cdmc_supported_bootstrap_classes())
  if (is.null(class_name)) {
    stop(
      "object must inherit from 'cdmc_fit', 'cdmc_dose_response', 'cdmc_dr_fit', or a supported diagnostic class.",
      call. = FALSE
    )
  }

  class_name
}

cdmc_bootstrap_statistics_map <- function(default = TRUE) {
  statistics <- list(
    cdmc_fit = c("coefficients", "average_tau"),
    cdmc_dynamic_estimand = "estimate",
    cdmc_dose_response = c("coefficients", "prediction"),
    cdmc_dr_fit = c(
      "coefficients",
      "average_tau_dr",
      "average_tau_linear",
      "lag_average_tau_dr",
      "dose_contrast_dr"
    ),
    cdmc_placebo_test = "mean_tau",
    cdmc_carryover_test = "mean_tau",
    cdmc_carryover_refit_test = "mean_tau",
    cdmc_joint_placebo_test = c("period_mean_tau", "window_mean_tau", "joint_p_value"),
    cdmc_equivalence_test = c("mean_tau", "equivalence_p_value"),
    cdmc_scia_test = c("f_statistic", "p_value", "restriction_p_value", "restriction_adj_p_value")
  )

  if (isTRUE(default)) {
    statistics$cdmc_dose_response <- "coefficients"
    statistics$cdmc_dr_fit <- c("coefficients", "average_tau_dr", "average_tau_linear")
    statistics$cdmc_scia_test <- c("f_statistic", "p_value")
  }

  statistics
}

cdmc_default_bootstrap_statistics <- function(object) {
  cdmc_bootstrap_statistics_map(default = TRUE)[[cdmc_bootstrap_object_class(object)]]
}

cdmc_supported_bootstrap_statistics <- function(object) {
  cdmc_bootstrap_statistics_map(default = FALSE)[[cdmc_bootstrap_object_class(object)]]
}

cdmc_bootstrap_source_object <- function(object) {
  if (cdmc_inherits_any(object, c("cdmc_fit", "cdmc_dr_fit"))) {
    return(object)
  }

  if (inherits(object, "cdmc_dynamic_estimand")) {
    source_object <- object$source_object
    if (inherits(source_object, "cdmc_dose_response")) {
      if (!inherits(source_object$fit_object, "cdmc_fit")) {
        stop(
          "Dynamic estimand bootstrap requires a stored parent cdmc_fit object for cdmc_dose_response sources.",
          call. = FALSE
        )
      }
      return(source_object$fit_object)
    }

    if (cdmc_inherits_any(source_object, c("cdmc_fit", "cdmc_dr_fit"))) {
      return(source_object)
    }

    stop(
      "Dynamic estimand bootstrap requires a source object inheriting from 'cdmc_fit', 'cdmc_dr_fit', or 'cdmc_dose_response'.",
      call. = FALSE
    )
  }

  if (inherits(object, "cdmc_dose_response")) {
    if (!inherits(object$fit_object, "cdmc_fit")) {
      stop(
        "Dose-response bootstrap requires a stored parent cdmc_fit object. Recreate the dose-response fit with the current package version.",
        call. = FALSE
      )
    }
    return(object$fit_object)
  }

  if (cdmc_inherits_any(object, c("cdmc_placebo_test", "cdmc_carryover_test", "cdmc_carryover_refit_test", "cdmc_joint_placebo_test", "cdmc_equivalence_test", "cdmc_scia_test"))) {
    if (!inherits(object$fit_object, "cdmc_fit")) {
      stop(
        "Diagnostic bootstrap requires a stored parent cdmc_fit object. Recreate the diagnostic with the current package version.",
        call. = FALSE
      )
    }
    return(object$fit_object)
  }

  stop(
    "object must inherit from 'cdmc_fit', 'cdmc_dose_response', 'cdmc_dr_fit', or a supported diagnostic class.",
    call. = FALSE
  )
}

cdmc_bootstrap_period_label <- function(period) {
  if (period < 0L) {
    return(sprintf("m%d", abs(as.integer(period))))
  }
  if (period > 0L) {
    return(sprintf("p%d", as.integer(period)))
  }

  "0"
}

cdmc_bootstrap_safe_label <- function(label) {
  label <- gsub("[^A-Za-z0-9]+", "_", as.character(label))
  label <- gsub("_+", "_", label)
  label <- gsub("^_+|_+$", "", label)
  if (!nzchar(label)) {
    return("window")
  }

  label
}

cdmc_joint_placebo_mean_names <- function(object) {
  if (!is.data.frame(object$tests) || nrow(object$tests) < 1L) {
    return(character(0))
  }

  if (all(c("n_periods", "period", "window_name") %in% names(object$tests)) &&
      all(object$tests$n_periods == 1L) &&
      all(!is.na(object$tests$period)) &&
      all(object$tests$window_name == paste0("period_", vapply(object$tests$period, cdmc_bootstrap_period_label, character(1))))) {
    return(paste0(
      "joint_placebo_mean_tau_period_",
      vapply(object$tests$period, cdmc_bootstrap_period_label, character(1))
    ))
  }

  paste0(
    "joint_placebo_mean_tau_window_",
    vapply(
      if ("window_name" %in% names(object$tests)) object$tests$window_name else seq_len(nrow(object$tests)),
      cdmc_bootstrap_safe_label,
      character(1)
    )
  )
}

cdmc_dr_design_columns <- function(object) {
  object$effect$design_columns %||% paste0("dose_lag", seq.int(0L, object$lag_order))
}

cdmc_dr_coefficient_vector <- function(object) {
  design_columns <- cdmc_dr_design_columns(object)
  coefficient_vector <- numeric(length(design_columns))
  names(coefficient_vector) <- design_columns
  if (length(object$effect$coefficients) > 0L) {
    common_columns <- intersect(names(object$effect$coefficients), design_columns)
    coefficient_vector[common_columns] <- object$effect$coefficients[common_columns]
  }

  coefficient_vector
}

cdmc_dr_history_table <- function(object) {
  lag_array <- cdmc_build_lagged_doses(object$dose_matrix, lag_order = object$lag_order)
  sample_indices <- which(object$effect$sample_mask, arr.ind = TRUE)
  data.frame(cdmc_sample_array_matrix(lag_array, sample_indices), check.names = FALSE)
}

cdmc_resolve_dr_contrast_history <- function(object, history = NULL, default = c("ones", "zero")) {
  design_columns <- cdmc_dr_design_columns(object)
  default <- match.arg(default)

  if (is.null(history)) {
    values <- if (identical(default, "ones")) rep(1, length(design_columns)) else rep(0, length(design_columns))
    names(values) <- design_columns
    return(values)
  }

  resolved <- if (is.data.frame(history) || is.matrix(history)) {
    history <- as.data.frame(history)
    if (nrow(history) != 1L) {
      stop("contrast_history and reference_history must contain exactly one row when supplied as data frames or matrices.", call. = FALSE)
    }
    missing_columns <- setdiff(design_columns, names(history))
    if (length(missing_columns) > 0L) {
      stop(
        sprintf("History input is missing required lag columns: %s.", paste(missing_columns, collapse = ", ")),
        call. = FALSE
      )
    }
    as.numeric(history[1L, design_columns, drop = TRUE])
  } else if (is.numeric(history)) {
    if (length(history) == 1L) {
      values <- rep(0, length(design_columns))
      values[[1L]] <- history
      values
    } else if (length(history) == length(design_columns)) {
      as.numeric(history)
    } else if (!is.null(names(history))) {
      values <- rep(0, length(design_columns))
      names(values) <- design_columns
      unknown_names <- setdiff(names(history), design_columns)
      if (length(unknown_names) > 0L) {
        stop(
          sprintf("Unknown lag names in history input: %s.", paste(unknown_names, collapse = ", ")),
          call. = FALSE
        )
      }
      values[names(history)] <- history
      unname(values)
    } else {
      stop(
        sprintf("Numeric history inputs must have length 1 or %d.", length(design_columns)),
        call. = FALSE
      )
    }
  } else {
    stop("History inputs must be NULL, numeric, a matrix, or a data frame.", call. = FALSE)
  }

  names(resolved) <- design_columns
  resolved
}

cdmc_compute_dr_contrast <- function(object, contrast_history = NULL, reference_history = NULL) {
  coefficient_vector <- cdmc_dr_coefficient_vector(object)
  contrast_values <- cdmc_resolve_dr_contrast_history(
    object = object,
    history = contrast_history,
    default = "ones"
  )
  reference_values <- cdmc_resolve_dr_contrast_history(
    object = object,
    history = reference_history,
    default = "zero"
  )

  sum((contrast_values - reference_values) * coefficient_vector)
}

cdmc_bootstrap_evaluate_target <- function(object, bootstrap_fit) {
  if (inherits(object, "cdmc_fit") || inherits(object, "cdmc_dr_fit")) {
    return(bootstrap_fit)
  }

  if (inherits(object, "cdmc_dynamic_estimand")) {
    dynamic_estimand <- get("cdmc_dynamic_estimand", mode = "function")
    source_target <- object$source_object
    if (inherits(source_target, "cdmc_dose_response")) {
      source_target <- cdmc_dose_response(
        bootstrap_fit,
        model = source_target$model,
        lag_order = source_target$lag_order,
        df = source_target$design_info$df %||% 4L,
        forest_trees = source_target$forest_trees %||% 200L,
        forest_mtry = source_target$forest_mtry,
        forest_min_node_size = source_target$forest_min_node_size,
        include_zero_dose = source_target$include_zero_dose,
        weights = source_target$weights
      )
    } else if (!inherits(source_target, "cdmc_fit") && !inherits(source_target, "cdmc_dr_fit")) {
      stop(
        "Dynamic estimand bootstrap requires a source object inheriting from 'cdmc_fit', 'cdmc_dr_fit', or 'cdmc_dose_response'.",
        call. = FALSE
      )
    } else {
      source_target <- bootstrap_fit
    }

    return(dynamic_estimand(
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

  if (inherits(object, "cdmc_dose_response")) {
    return(cdmc_dose_response(
      bootstrap_fit,
      model = object$model,
      lag_order = object$lag_order,
      df = object$design_info$df %||% 4L,
      forest_trees = object$forest_trees %||% 200L,
      forest_mtry = object$forest_mtry,
      forest_min_node_size = object$forest_min_node_size,
      include_zero_dose = object$include_zero_dose,
      weights = object$weights
    ))
  }

  if (inherits(object, "cdmc_placebo_test")) {
    return(cdmc_placebo_test(
      bootstrap_fit,
      periods = object$periods,
      rerun_tuning = object$rerun_tuning %||% FALSE,
      verbose = FALSE
    ))
  }

  if (inherits(object, "cdmc_carryover_test")) {
    return(cdmc_carryover_test(bootstrap_fit, periods = object$periods))
  }

  if (inherits(object, "cdmc_carryover_refit_test")) {
    return(cdmc_carryover_refit_test(
      bootstrap_fit,
      periods = object$periods,
      rerun_tuning = object$rerun_tuning %||% FALSE,
      verbose = FALSE
    ))
  }

  if (inherits(object, "cdmc_joint_placebo_test")) {
    return(cdmc_joint_placebo_test(
      bootstrap_fit,
      periods = object$periods,
      placebo_windows = object$placebo_windows %||% NULL,
      alpha = object$alpha,
      equivalence_margin = object$equivalence_margin,
      rerun_tuning = object$rerun_tuning %||% FALSE,
      verbose = FALSE
    ))
  }

  if (inherits(object, "cdmc_equivalence_test")) {
    diagnostic_target <- cdmc_bootstrap_evaluate_target(object$diagnostic_object, bootstrap_fit)
    return(cdmc_equivalence_test(
      diagnostic_target,
      margin = object$margin,
      alpha = object$alpha
    ))
  }

  if (inherits(object, "cdmc_scia_test")) {
    return(cdmc_scia_test(
      bootstrap_fit,
      lags = object$lags,
      outcome_proxy = object$outcome_proxy,
      covariates = object$covariates,
      include_current_covariates = object$include_current_covariates,
      include_covariate_lags = object$include_covariate_lags,
      include_unit_effects = object$include_unit_effects,
      include_time_effects = object$include_time_effects,
      restriction_blocks = object$restriction_blocks,
      p_adjust_method = object$p_adjust_method %||% "holm",
      alpha = object$alpha
    ))
  }

  stop(
    "object must inherit from 'cdmc_fit', 'cdmc_dose_response', 'cdmc_dr_fit', or a supported diagnostic class.",
    call. = FALSE
  )
}

cdmc_bootstrap_prediction_statistics <- function(
  object,
  prediction_dose = NULL,
  prediction_history = NULL,
  prediction_type = c("response", "slope")
) {
  prediction_type <- match.arg(prediction_type)
  prediction <- stats::predict(
    object,
    dose = prediction_dose,
    history = prediction_history,
    type = prediction_type
  )
  estimates <- prediction$estimate
  names(estimates) <- paste0("dose_response_", prediction_type, "_", seq_along(estimates))
  estimates
}

cdmc_bootstrap_original_columns <- function(object) {
  source_object <- cdmc_bootstrap_source_object(object)

  if (inherits(source_object, "cdmc_fit")) {
    fit_control <- source_object$fit_control %||% list(weights = source_object$weights %||% NULL)
    return(unique(c(
      source_object$outcome,
      source_object$dose,
      source_object$unit,
      source_object$time,
      source_object$covariates,
      fit_control$weights
    )))
  }

  if (inherits(source_object, "cdmc_dr_fit")) {
    fit_control <- source_object$fit_control %||% list(
      weights = if (identical(source_object$weight_method, "external")) ".cdmc_weight" else NULL,
      weight_covariates = source_object$weight_covariates %||% NULL
    )
    return(unique(c(
      source_object$outcome,
      source_object$dose,
      source_object$unit,
      source_object$time,
      source_object$covariates,
      fit_control$weights,
      fit_control$weight_covariates
    )))
  }

  stop("Bootstrap source object must inherit from 'cdmc_fit' or 'cdmc_dr_fit'.", call. = FALSE)
}

cdmc_bootstrap_lambda_method <- function(object) {
  source_object <- cdmc_bootstrap_source_object(object)

  if (inherits(source_object, "cdmc_fit")) {
    return(source_object$lambda_tuning$method)
  }
  if (inherits(source_object, "cdmc_dr_fit")) {
    return(if (is.null(source_object$lambda)) source_object$lambda_selection else "fixed")
  }

  stop("Bootstrap source object must inherit from 'cdmc_fit' or 'cdmc_dr_fit'.", call. = FALSE)
}

cdmc_bootstrap_fit_function <- function(object) {
  source_object <- cdmc_bootstrap_source_object(object)

  if (inherits(source_object, "cdmc_fit")) {
    return(cdmc_fit)
  }
  if (inherits(source_object, "cdmc_dr_fit")) {
    return(cdmc_dr_fit)
  }

  stop("Bootstrap source object must inherit from 'cdmc_fit' or 'cdmc_dr_fit'.", call. = FALSE)
}

cdmc_collect_bootstrap_statistics <- function(
  object,
  statistics,
  contrast_history = NULL,
  reference_history = NULL,
  prediction_dose = NULL,
  prediction_history = NULL,
  prediction_type = c("response", "slope")
) {
  statistics <- unique(statistics)
  out <- numeric(0)
  prediction_type <- match.arg(prediction_type)

  if (inherits(object, "cdmc_dynamic_estimand")) {
    if ("estimate" %in% statistics) {
      out <- c(out, object$estimates)
    }

    return(out)
  }

  if (inherits(object, "cdmc_dose_response")) {
    if ("coefficients" %in% statistics) {
      coefficients <- object$coefficients
      if (length(coefficients) > 0L) {
        names(coefficients) <- paste0("coef_", names(coefficients))
        out <- c(out, coefficients)
      }
    }

    if ("prediction" %in% statistics) {
      out <- c(
        out,
        cdmc_bootstrap_prediction_statistics(
          object,
          prediction_dose = prediction_dose,
          prediction_history = prediction_history,
          prediction_type = prediction_type
        )
      )
    }

    return(out)
  }

  if ("coefficients" %in% statistics) {
    coefficients <- object$effect$coefficients
    if (length(coefficients) > 0L) {
      names(coefficients) <- paste0("coef_", names(coefficients))
      out <- c(out, coefficients)
    }
  }

  if (inherits(object, "cdmc_fit") && "average_tau" %in% statistics) {
    active_mask <- cdmc_active_dose_mask(object$dose_matrix, zero_tolerance = object$zero_tolerance)
    if (any(active_mask)) {
      out <- c(
        out,
        mean_tau_active_dose = if (isTRUE(object$weight_supplied)) {
          stats::weighted.mean(object$effect$tau[active_mask], w = object$weight_matrix[active_mask])
        } else {
          mean(object$effect$tau[active_mask])
        }
      )
    } else {
      out <- c(out, mean_tau_active_dose = NA_real_)
    }
  }

  if (inherits(object, "cdmc_dr_fit")) {
    dr_mask <- object$effect$sample_mask

    if ("average_tau_dr" %in% statistics) {
      out <- c(
        out,
        mean_tau_dr = if (any(dr_mask, na.rm = TRUE)) {
          mean(object$effect$tau_dr[dr_mask], na.rm = TRUE)
        } else {
          NA_real_
        }
      )
    }

    if ("average_tau_linear" %in% statistics) {
      out <- c(
        out,
        mean_tau_dr_linear = if (any(dr_mask, na.rm = TRUE)) {
          mean(object$effect$fitted[dr_mask], na.rm = TRUE)
        } else {
          NA_real_
        }
      )
    }

    if ("lag_average_tau_dr" %in% statistics) {
      history <- cdmc_dr_history_table(object)
      coefficient_vector <- cdmc_dr_coefficient_vector(object)
      lag_averages <- if (nrow(history) > 0L) {
        colMeans(sweep(as.matrix(history), 2, coefficient_vector, FUN = "*"), na.rm = TRUE)
      } else {
        rep(NA_real_, length(coefficient_vector))
      }
      names(lag_averages) <- paste0("mean_tau_dr_", names(coefficient_vector))
      out <- c(out, lag_averages)
    }

    if ("dose_contrast_dr" %in% statistics) {
      out <- c(
        out,
        dr_path_contrast = cdmc_compute_dr_contrast(
          object = object,
          contrast_history = contrast_history,
          reference_history = reference_history
        )
      )
    }
  }

  if (inherits(object, "cdmc_placebo_test") && "mean_tau" %in% statistics) {
    out <- c(out, mean_placebo_tau = object$mean_tau)
  }

  if (inherits(object, "cdmc_carryover_test") && "mean_tau" %in% statistics) {
    out <- c(out, mean_carryover_tau = object$mean_tau)
  }

  if (inherits(object, "cdmc_carryover_refit_test") && "mean_tau" %in% statistics) {
    out <- c(out, mean_refit_carryover_tau = object$mean_tau)
  }

  if (inherits(object, "cdmc_joint_placebo_test")) {
    if ("period_mean_tau" %in% statistics || "window_mean_tau" %in% statistics) {
      period_means <- object$tests$mean_tau
      names(period_means) <- cdmc_joint_placebo_mean_names(object)
      out <- c(out, period_means)
    }

    if ("joint_p_value" %in% statistics) {
      out <- c(out, joint_placebo_p_value = object$joint_p_value)
    }
  }

  if (inherits(object, "cdmc_equivalence_test")) {
    if ("mean_tau" %in% statistics) {
      out <- c(out, equivalence_mean_tau = object$mean_tau)
    }
    if ("equivalence_p_value" %in% statistics) {
      out <- c(out, equivalence_p_value = object$p_value)
    }
  }

  if (inherits(object, "cdmc_scia_test")) {
    if ("f_statistic" %in% statistics) {
      out <- c(out, scia_f_statistic = object$f_statistic)
    }
    if ("p_value" %in% statistics) {
      out <- c(out, scia_p_value = object$p_value)
    }
    if ("restriction_p_value" %in% statistics && nrow(object$restriction_table) > 0L) {
      block_p_values <- object$restriction_table$p_value
      names(block_p_values) <- paste0("scia_restriction_p_value_", object$restriction_table$restriction_name)
      out <- c(out, block_p_values)
    }
    if ("restriction_adj_p_value" %in% statistics && nrow(object$restriction_table) > 0L) {
      block_adj_p_values <- object$restriction_table$adjusted_p_value
      names(block_adj_p_values) <- paste0("scia_restriction_adj_p_value_", object$restriction_table$restriction_name)
      out <- c(out, block_adj_p_values)
    }
  }

  out
}

cdmc_build_bootstrap_fit_spec <- function(object, rerun_tuning = NULL) {
  source_object <- cdmc_bootstrap_source_object(object)
  if (!identical(source_object, object)) {
    return(cdmc_build_bootstrap_fit_spec(source_object, rerun_tuning = rerun_tuning))
  }

  if (inherits(object, "cdmc_fit")) {
    fit_control <- object$fit_control
    if (is.null(fit_control)) {
      fit_control <- list(
        lambda = if (identical(object$lambda_tuning$method, "fixed")) object$lambda else NULL,
        weights = object$weights %||% NULL,
        rank_max = object$rank_max,
        lambda_fraction = object$lambda_tuning$lambda_fraction %||% 0.25,
        lambda_selection = if (identical(object$lambda_tuning$method, "cv")) "cv" else "heuristic",
        lambda_grid = object$lambda_tuning$lambda_grid,
        nlambda = length(object$lambda_tuning$lambda_grid %||% numeric(0)),
        lambda_min_ratio = 0.05,
        cv_rounds = object$lambda_tuning$cv_rounds %||% 5L,
        cv_block_size = object$lambda_tuning$cv_block_size_requested %||% 2L,
        cv_workers = object$fit_control$cv_workers %||% 1L,
        cv_top_k = object$fit_control$cv_top_k %||% NULL,
        cv_coarse_to_fine = object$fit_control$cv_coarse_to_fine %||% FALSE,
        cv_coarse_nlambda = object$fit_control$cv_coarse_nlambda %||% NULL,
        cv_warm_starts = object$fit_control$cv_warm_starts %||% FALSE,
        washout = object$washout,
        lag_order = object$lag_order,
        effect_model = object$effect_model,
        effect_df = object$effect_df %||% 4L,
        objective = object$objective %||% "staged",
        outer_maxit = 20L,
        fe_maxit = 200L,
        soft_maxit = 100L,
        tol = 1e-5,
        fe_tol = 1e-8,
        zero_tolerance = object$zero_tolerance
      )
    }

    if (is.null(rerun_tuning)) {
      rerun_tuning <- !identical(object$lambda_tuning$method, "fixed")
    }
    rerun_tuning <- isTRUE(rerun_tuning)

    fit_spec <- fit_control
    if (rerun_tuning && !identical(object$lambda_tuning$method, "fixed")) {
      fit_spec$lambda <- NULL
    } else {
      fit_spec$lambda <- object$lambda
    }

    fit_spec$verbose <- FALSE
    return(list(fit_spec = fit_spec, rerun_tuning = rerun_tuning))
  }

  if (inherits(object, "cdmc_dr_fit")) {
    fit_control <- object$fit_control %||% list(
      weights = if (identical(object$weight_method, "external")) ".cdmc_weight" else NULL,
      weight_method = object$weight_method,
      weight_covariates = object$weight_covariates %||% NULL,
      gps_time_effects = object$gps_time_effects,
      gps_model = object$gps_model %||% "linear",
      gps_df = object$gps_df %||% 4L,
      gps_spline_covariates = object$gps_spline_covariates %||% object$weight_covariates %||% NULL,
      gps_stack_models = object$gps_stack_models %||% NULL,
      gps_bandwidth = object$gps_bandwidth %||% NULL,
      gps_forest_trees = object$gps_forest_trees %||% 200L,
      gps_forest_mtry = object$gps_forest_mtry %||% NULL,
      gps_forest_min_node_size = object$gps_forest_min_node_size %||% NULL,
      gps_boost_trees = object$gps_boost_trees %||% 200L,
      gps_boost_depth = object$gps_boost_depth %||% 2L,
      gps_boost_shrinkage = object$gps_boost_shrinkage %||% 0.05,
      gps_boost_min_obs_node = object$gps_boost_min_obs_node %||% 10L,
      cbps_standardize = object$cbps_standardize %||% TRUE,
      cbps_method = object$cbps_method %||% "over",
      cbps_iterations = object$cbps_iterations %||% 1000L,
      cbps_twostep = object$cbps_twostep %||% TRUE,
      adaptive_balance_methods = object$adaptive_balance_methods %||% NULL,
      entropy_balance_degree = object$entropy_balance_degree %||% 1L,
      entropy_balance_standardize = object$entropy_balance_standardize %||% TRUE,
      entropy_balance_iterations = object$entropy_balance_iterations %||% 1000L,
      entropy_balance_reltol = object$entropy_balance_reltol %||% 1e-8,
      kernel_balance_degree = object$kernel_balance_degree %||% 1L,
      kernel_balance_centers = object$kernel_balance_centers %||% 25L,
      kernel_balance_bandwidth = object$kernel_balance_bandwidth %||% NULL,
      kernel_balance_standardize = object$kernel_balance_standardize %||% TRUE,
      kernel_balance_iterations = object$kernel_balance_iterations %||% 1000L,
      kernel_balance_reltol = object$kernel_balance_reltol %||% 1e-8,
      stabilize_weights = object$stabilize_weights,
      max_weight = object$max_weight,
      n_folds = object$n_folds,
      dr_workers = object$dr_workers %||% object$fit_control$dr_workers %||% 1L,
      lambda = object$lambda,
      rank_max = object$rank_max %||% NULL,
      lambda_fraction = 0.25,
      lambda_selection = object$lambda_selection,
      lambda_grid = NULL,
      nlambda = 5L,
      lambda_min_ratio = 0.05,
      cv_rounds = 5L,
      cv_block_size = 2L,
      cv_workers = object$fit_control$cv_workers %||% 1L,
      cv_top_k = object$fit_control$cv_top_k %||% NULL,
      cv_coarse_to_fine = object$fit_control$cv_coarse_to_fine %||% FALSE,
      cv_coarse_nlambda = object$fit_control$cv_coarse_nlambda %||% NULL,
      cv_warm_starts = object$fit_control$cv_warm_starts %||% FALSE,
      washout = object$washout,
      lag_order = object$lag_order,
      outer_maxit = 20L,
      fe_maxit = 200L,
      soft_maxit = 100L,
      tol = 1e-5,
      fe_tol = 1e-8,
      zero_tolerance = object$zero_tolerance %||% 1e-8
    )

    if (is.null(rerun_tuning)) {
      rerun_tuning <- is.null(object$lambda)
    }
    rerun_tuning <- isTRUE(rerun_tuning)

    if (!rerun_tuning && is.null(object$lambda)) {
      stop(
        "Bootstrap for cdmc_dr_fit objects with automatic foldwise lambda selection requires rerun_tuning = TRUE.",
        call. = FALSE
      )
    }

    fit_spec <- fit_control
    fit_spec$lambda <- if (rerun_tuning) NULL else object$lambda
    fit_spec$verbose <- FALSE
    return(list(fit_spec = fit_spec, rerun_tuning = rerun_tuning))
  }

  stop("Bootstrap source object must inherit from 'cdmc_fit' or 'cdmc_dr_fit'.", call. = FALSE)
}

cdmc_resolve_bootstrap_simultaneous <- function(object, simultaneous, n_statistics) {
  if (is.null(simultaneous)) {
    return(inherits(object, "cdmc_dynamic_estimand") && n_statistics > 1L)
  }

  if (!is.logical(simultaneous) || length(simultaneous) != 1L || is.na(simultaneous)) {
    stop("simultaneous must be TRUE, FALSE, or NULL.", call. = FALSE)
  }

  if (isTRUE(simultaneous) && !inherits(object, "cdmc_dynamic_estimand")) {
    stop(
      "simultaneous bands are currently supported only for cdmc_dynamic_estimand bootstrap objects.",
      call. = FALSE
    )
  }

  isTRUE(simultaneous)
}

cdmc_bootstrap_simultaneous_bands <- function(estimates, draws, conf_level) {
  joint_draws <- draws[, names(estimates), drop = FALSE]
  complete_rows <- stats::complete.cases(joint_draws)
  joint_draws <- joint_draws[complete_rows, , drop = FALSE]
  estimates <- unname(estimates)
  names(estimates) <- colnames(joint_draws)

  if (nrow(joint_draws) == 0L) {
    return(list(
      enabled = TRUE,
      method = "max-t",
      critical_value = NA_real_,
      n_success = 0L,
      conf_low = rep(NA_real_, length(estimates)),
      conf_high = rep(NA_real_, length(estimates))
    ))
  }

  bootstrap_sd <- vapply(seq_along(estimates), function(index) {
    stats::sd(joint_draws[, index], na.rm = TRUE)
  }, numeric(1))
  stable <- is.finite(bootstrap_sd) & bootstrap_sd > sqrt(.Machine$double.eps)

  conf_low <- estimates
  conf_high <- estimates
  critical_value <- 0

  if (any(stable)) {
    standardized <- sweep(
      abs(sweep(joint_draws[, stable, drop = FALSE], 2, estimates[stable], FUN = "-")),
      2,
      bootstrap_sd[stable],
      FUN = "/"
    )
    max_t <- vapply(seq_len(nrow(standardized)), function(index) {
      max(standardized[index, ], na.rm = TRUE)
    }, numeric(1))
    critical_value <- stats::quantile(max_t, probs = conf_level, names = FALSE, na.rm = TRUE)
    conf_low[stable] <- estimates[stable] - critical_value * bootstrap_sd[stable]
    conf_high[stable] <- estimates[stable] + critical_value * bootstrap_sd[stable]
  }

  list(
    enabled = TRUE,
    method = "max-t",
    critical_value = unname(critical_value),
    n_success = nrow(joint_draws),
    conf_low = conf_low,
    conf_high = conf_high
  )
}

cdmc_bootstrap_summary_table <- function(estimates, draws, conf_level, simultaneous = NULL) {
  alpha <- (1 - conf_level) / 2
  summary_table <- data.frame(
    statistic = names(estimates),
    estimate = unname(estimates),
    bootstrap_mean = vapply(seq_along(estimates), function(index) {
      mean(draws[, index], na.rm = TRUE)
    }, numeric(1)),
    bootstrap_sd = vapply(seq_along(estimates), function(index) {
      stats::sd(draws[, index], na.rm = TRUE)
    }, numeric(1)),
    conf_low = vapply(seq_along(estimates), function(index) {
      stats::quantile(draws[, index], probs = alpha, na.rm = TRUE, names = FALSE)
    }, numeric(1)),
    conf_high = vapply(seq_along(estimates), function(index) {
      stats::quantile(draws[, index], probs = 1 - alpha, na.rm = TRUE, names = FALSE)
    }, numeric(1)),
    n_success = vapply(seq_along(estimates), function(index) {
      sum(is.finite(draws[, index]))
    }, numeric(1))
  )

  if (!is.null(simultaneous) && isTRUE(simultaneous$enabled)) {
    summary_table$simultaneous_conf_low <- simultaneous$conf_low
    summary_table$simultaneous_conf_high <- simultaneous$conf_high
  }

  summary_table
}

cdmc_bootstrap <- function(
  object,
  n_boot = 100L,
  statistics = c("coefficients", "average_tau"),
  rerun_tuning = NULL,
  conf_level = 0.95,
  simultaneous = NULL,
  contrast_history = NULL,
  reference_history = NULL,
  prediction_dose = NULL,
  prediction_history = NULL,
  prediction_type = c("response", "slope"),
  workers = 1L,
  seed = NULL,
  verbose = FALSE
) {
  missing_statistics <- missing(statistics)

  cdmc_bootstrap_object_class(object)

  prediction_type <- match.arg(prediction_type)

  if (!inherits(object, "cdmc_dr_fit") && (!is.null(contrast_history) || !is.null(reference_history))) {
    stop("contrast_history and reference_history are only supported for cdmc_dr_fit bootstrap summaries.", call. = FALSE)
  }

  if (!inherits(object, "cdmc_dose_response") && (!is.null(prediction_dose) || !is.null(prediction_history))) {
    stop("prediction_dose and prediction_history are only supported for cdmc_dose_response bootstrap summaries.", call. = FALSE)
  }

  n_boot <- as.integer(n_boot)
  if (n_boot < 1L) {
    stop("n_boot must be a positive integer.", call. = FALSE)
  }

  if (!is.numeric(conf_level) || length(conf_level) != 1L || conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must be a scalar in (0, 1).", call. = FALSE)
  }

  if (!is.numeric(workers) || length(workers) != 1L || !is.finite(workers) || workers < 1 || workers != floor(workers)) {
    stop("workers must be a positive integer.", call. = FALSE)
  }
  workers <- as.integer(workers)

  statistics <- if (missing_statistics) {
    cdmc_default_bootstrap_statistics(object)
  } else {
    unique(match.arg(
      statistics,
      cdmc_supported_bootstrap_statistics(object),
      several.ok = TRUE
    ))
  }
  if (length(statistics) == 0L) {
    stop("statistics must request at least one bootstrap target.", call. = FALSE)
  }

  if (inherits(object, "cdmc_dose_response") &&
      "prediction" %in% statistics &&
      is.null(prediction_dose) &&
      is.null(prediction_history)) {
    stop(
      "Provide prediction_dose or prediction_history when requesting the 'prediction' bootstrap statistic for cdmc_dose_response objects.",
      call. = FALSE
    )
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  use_parallel <- workers > 1L
  if (use_parallel && identical(.Platform$OS.type, "windows")) {
    warning(
      "Parallel bootstrap currently uses multicore execution and is not available on Windows. Falling back to sequential execution.",
      call. = FALSE
    )
    use_parallel <- FALSE
    workers <- 1L
  }
  worker_count <- if (use_parallel) min(workers, n_boot) else 1L
  replicate_seeds <- if (use_parallel) sample.int(.Machine$integer.max, n_boot, replace = TRUE) else integer(0)

  source_object <- cdmc_bootstrap_source_object(object)
  original_columns <- cdmc_bootstrap_original_columns(object)
  bootstrap_spec <- cdmc_build_bootstrap_fit_spec(object, rerun_tuning = rerun_tuning)
  reference_statistics <- cdmc_collect_bootstrap_statistics(
    object,
    statistics = statistics,
    contrast_history = contrast_history,
    reference_history = reference_history,
    prediction_dose = prediction_dose,
    prediction_history = prediction_history,
    prediction_type = prediction_type
  )
  fit_function <- cdmc_bootstrap_fit_function(object)

  if (length(reference_statistics) == 0L) {
    stop(
      "No bootstrap statistics are available for the current fit and requested statistic set.",
      call. = FALSE
    )
  }

  simultaneous <- cdmc_resolve_bootstrap_simultaneous(
    object = object,
    simultaneous = simultaneous,
    n_statistics = length(reference_statistics)
  )

  bootstrap_draws <- matrix(
    NA_real_,
    nrow = n_boot,
    ncol = length(reference_statistics),
    dimnames = list(paste0("boot", seq_len(n_boot)), names(reference_statistics))
  )
  errors <- character(0)

  run_bootstrap_replicate <- function(bootstrap_index, replicate_seed = NULL) {
    if (!is.null(replicate_seed)) {
      set.seed(replicate_seed)
    }

    if (verbose && !use_parallel) {
      message(sprintf("bootstrap replicate %d/%d", bootstrap_index, n_boot))
    }

    bootstrap_data <- cdmc_bootstrap_sample_data(
      data = source_object$data,
      unit = source_object$unit,
      original_columns = original_columns
    )

    bootstrap_fit <- tryCatch(
      do.call(
        fit_function,
        c(
          list(
            data = bootstrap_data,
            outcome = source_object$outcome,
            dose = source_object$dose,
            unit = source_object$unit,
            time = source_object$time,
            covariates = source_object$covariates
          ),
          bootstrap_spec$fit_spec
        )
      ),
      error = function(error) error
    )

    if (inherits(bootstrap_fit, "error")) {
      return(list(index = bootstrap_index, statistics = NULL, error = conditionMessage(bootstrap_fit)))
    }

    bootstrap_target <- tryCatch(
      cdmc_bootstrap_evaluate_target(object, bootstrap_fit),
      error = function(error) error
    )

    if (inherits(bootstrap_target, "error")) {
      return(list(index = bootstrap_index, statistics = NULL, error = conditionMessage(bootstrap_target)))
    }

    bootstrap_statistics <- cdmc_collect_bootstrap_statistics(
      bootstrap_target,
      statistics = statistics,
      contrast_history = contrast_history,
      reference_history = reference_history,
      prediction_dose = prediction_dose,
      prediction_history = prediction_history,
      prediction_type = prediction_type
    )

    list(index = bootstrap_index, statistics = bootstrap_statistics, error = NULL)
  }

  if (use_parallel && verbose) {
    message(sprintf("running bootstrap in parallel with %d workers", worker_count))
  }

  replicate_results <- if (use_parallel) {
    parallel::mclapply(
      seq_len(n_boot),
      function(bootstrap_index) {
        run_bootstrap_replicate(
          bootstrap_index = bootstrap_index,
          replicate_seed = replicate_seeds[[bootstrap_index]]
        )
      },
      mc.cores = worker_count,
      mc.set.seed = FALSE
    )
  } else {
    lapply(seq_len(n_boot), function(bootstrap_index) {
      run_bootstrap_replicate(bootstrap_index = bootstrap_index)
    })
  }

  for (replicate_result in replicate_results) {
    if (!is.null(replicate_result$error)) {
      errors <- c(errors, replicate_result$error)
      next
    }

    bootstrap_draws[replicate_result$index, names(replicate_result$statistics)] <- replicate_result$statistics
  }

  simultaneous_result <- if (isTRUE(simultaneous)) {
    cdmc_bootstrap_simultaneous_bands(
      estimates = reference_statistics,
      draws = bootstrap_draws,
      conf_level = conf_level
    )
  } else {
    NULL
  }

  summary_table <- cdmc_bootstrap_summary_table(
    estimates = reference_statistics,
    draws = bootstrap_draws,
    conf_level = conf_level,
    simultaneous = simultaneous_result
  )

  result <- list(
    call = match.call(),
    statistics = statistics,
    estimates = reference_statistics,
    draws = bootstrap_draws,
    summary = summary_table,
    n_boot = n_boot,
    conf_level = conf_level,
    simultaneous = simultaneous_result,
    workers = worker_count,
    parallel = use_parallel,
    rerun_tuning = bootstrap_spec$rerun_tuning,
    lambda_method = cdmc_bootstrap_lambda_method(object),
    n_failures = length(errors),
    errors = errors
  )

  class(result) <- "cdmc_bootstrap"
  result
}

print.cdmc_bootstrap <- function(x, ...) {
  cat("causaldosemc bootstrap inference\n")
  cat(sprintf("  bootstrap replications: %d\n", x$n_boot))
  cat(sprintf("  successful replications: %d\n", x$n_boot - x$n_failures))
  cat(sprintf("  confidence level: %.3f\n", x$conf_level))
  cat(sprintf("  execution mode: %s (%d worker%s)\n", if (isTRUE(x$parallel)) "parallel" else "sequential", x$workers %||% 1L, if ((x$workers %||% 1L) == 1L) "" else "s"))
  if (!is.null(x$simultaneous) && isTRUE(x$simultaneous$enabled)) {
    cat(sprintf("  simultaneous bands: %s\n", x$simultaneous$method))
    cat(sprintf("  joint successful replications: %d\n", x$simultaneous$n_success))
  }
  cat(sprintf("  lambda method in original fit: %s\n", x$lambda_method))
  cat(sprintf("  rerun tuning in bootstrap: %s\n", if (x$rerun_tuning) "yes" else "no"))
  if (x$n_failures > 0L) {
    cat(sprintf("  failures: %d\n", x$n_failures))
  }
  print(x$summary, row.names = FALSE)
  invisible(x)
}
