cdmc_sensitivity_supported_classes <- function() {
  c("cdmc_fit", "cdmc_dr_fit", "cdmc_dose_response", "cdmc_dynamic_estimand")
}

cdmc_sensitivity_first_class <- function(object, classes) {
  matches <- vapply(classes, function(class_name) inherits(object, class_name), logical(1))
  if (!any(matches)) {
    return(NULL)
  }

  classes[[which(matches)[1L]]]
}

cdmc_sensitivity_target_class <- function(object) {
  class_name <- cdmc_sensitivity_first_class(object, cdmc_sensitivity_supported_classes())
  if (is.null(class_name)) {
    stop(
      "object must inherit from 'cdmc_fit', 'cdmc_dr_fit', 'cdmc_dose_response', or 'cdmc_dynamic_estimand'.",
      call. = FALSE
    )
  }

  class_name
}

cdmc_sensitivity_source_class <- function(source_object) {
  class_name <- cdmc_sensitivity_first_class(source_object, c("cdmc_fit", "cdmc_dr_fit"))
  if (is.null(class_name)) {
    stop(
      "source_object must inherit from 'cdmc_fit' or 'cdmc_dr_fit'.",
      call. = FALSE
    )
  }

  class_name
}

cdmc_sensitivity_source_object <- function(object) {
  cdmc_sensitivity_target_class(object)

  cdmc_bootstrap_source_object(object)
}

cdmc_default_washout_grid <- function(object) {
  source_object <- cdmc_sensitivity_source_object(object)
  unique(pmax(0L, c(source_object$washout - 1L, source_object$washout, source_object$washout + 1L)))
}

cdmc_default_zero_tolerance_grid <- function(object) {
  source_object <- cdmc_sensitivity_source_object(object)
  zero_tolerance <- source_object$zero_tolerance %||% 1e-8
  sort(unique(c(zero_tolerance / 10, zero_tolerance, zero_tolerance * 10)))
}

cdmc_sensitivity_support_fraction <- function(observed_controls, eligible_controls) {
  if (!is.finite(observed_controls) || observed_controls <= 0 || !is.finite(eligible_controls)) {
    return(NA_real_)
  }

  as.numeric(eligible_controls) / as.numeric(observed_controls)
}

cdmc_sensitivity_fit_objective <- function(source_object) {
  source_class <- cdmc_sensitivity_source_class(source_object)

  if (identical(source_class, "cdmc_fit")) {
    return(source_object$objective %||% "staged")
  }

  "dr"
}

cdmc_sensitivity_optimization_sample <- function(source_object) {
  source_class <- cdmc_sensitivity_source_class(source_object)

  if (identical(source_class, "cdmc_fit") && identical(cdmc_sensitivity_fit_objective(source_object), "joint")) {
    return("observed_panel")
  }

  "eligible_zero_dose"
}

cdmc_sensitivity_named_weight_specs <- function(weight_specs) {
  if (is.null(weight_specs)) {
    return(NULL)
  }

  if (!is.list(weight_specs)) {
    weight_specs <- list(custom = weight_specs)
  }

  weight_names <- names(weight_specs)
  if (is.null(weight_names)) {
    weight_names <- rep("", length(weight_specs))
  }
  empty_names <- weight_names == ""
  weight_names[empty_names] <- paste0("custom", seq_len(sum(empty_names)))
  names(weight_specs) <- weight_names
  weight_specs
}

cdmc_sensitivity_fit_weight_scenario <- function(weights) {
  list(
    weights = weights,
    weight_mode = if (is.null(weights)) "unweighted" else "weighted"
  )
}

cdmc_sensitivity_fit_weight_value <- function(spec) {
  if (is.list(spec) && !is.null(names(spec)) && "weights" %in% names(spec)) {
    return(spec$weights)
  }

  spec
}

cdmc_sensitivity_default_dr_weight_scenario <- function(source_object) {
  list(
    weights = if (identical(source_object$weight_method, "external")) {
      source_object$fit_control$weights %||% source_object$weights
    } else {
      NULL
    },
    weight_method = source_object$weight_method,
    weight_covariates = source_object$weight_covariates,
    gps_time_effects = source_object$gps_time_effects,
    gps_model = source_object$gps_model %||% "linear",
    gps_df = source_object$gps_df %||% 4L,
    gps_spline_covariates = source_object$gps_spline_covariates,
    gps_stack_models = source_object$gps_stack_models %||% NULL,
    gps_bandwidth = source_object$gps_bandwidth %||% NULL,
    gps_forest_trees = source_object$gps_forest_trees %||% 200L,
    gps_forest_mtry = source_object$gps_forest_mtry %||% NULL,
    gps_forest_min_node_size = source_object$gps_forest_min_node_size %||% NULL,
    gps_boost_trees = source_object$gps_boost_trees %||% 200L,
    gps_boost_depth = source_object$gps_boost_depth %||% 2L,
    gps_boost_shrinkage = source_object$gps_boost_shrinkage %||% 0.05,
    gps_boost_min_obs_node = source_object$gps_boost_min_obs_node %||% 10L,
    cbps_standardize = source_object$cbps_standardize %||% TRUE,
    cbps_method = source_object$cbps_method %||% "over",
    cbps_iterations = source_object$cbps_iterations %||% 1000L,
    cbps_twostep = source_object$cbps_twostep %||% TRUE,
    adaptive_balance_methods = source_object$adaptive_balance_methods %||% NULL,
    entropy_balance_degree = source_object$entropy_balance_degree %||% 1L,
    entropy_balance_standardize = source_object$entropy_balance_standardize %||% TRUE,
    entropy_balance_iterations = source_object$entropy_balance_iterations %||% 1000L,
    entropy_balance_reltol = source_object$entropy_balance_reltol %||% 1e-8,
    kernel_balance_degree = source_object$kernel_balance_degree %||% 1L,
    kernel_balance_centers = source_object$kernel_balance_centers %||% 25L,
    kernel_balance_bandwidth = source_object$kernel_balance_bandwidth %||% NULL,
    kernel_balance_standardize = source_object$kernel_balance_standardize %||% TRUE,
    kernel_balance_iterations = source_object$kernel_balance_iterations %||% 1000L,
    kernel_balance_reltol = source_object$kernel_balance_reltol %||% 1e-8,
    stabilize_weights = source_object$stabilize_weights,
    max_weight = source_object$max_weight,
    weight_mode = source_object$weight_method
  )
}

cdmc_sensitivity_parse_dr_weight_spec <- function(source_object, spec) {
  scenario <- cdmc_sensitivity_default_dr_weight_scenario(source_object)
  known_names <- c(
    "weights",
    "weight_method",
    "weight_covariates",
    "gps_time_effects",
    "gps_model",
    "gps_df",
    "gps_spline_covariates",
    "gps_stack_models",
    "gps_bandwidth",
    "gps_forest_trees",
    "gps_forest_mtry",
    "gps_forest_min_node_size",
    "gps_boost_trees",
    "gps_boost_depth",
    "gps_boost_shrinkage",
    "gps_boost_min_obs_node",
    "cbps_standardize",
    "cbps_method",
    "cbps_iterations",
    "cbps_twostep",
    "adaptive_balance_methods",
    "entropy_balance_degree",
    "entropy_balance_standardize",
    "entropy_balance_iterations",
    "entropy_balance_reltol",
    "kernel_balance_degree",
    "kernel_balance_centers",
    "kernel_balance_bandwidth",
    "kernel_balance_standardize",
    "kernel_balance_iterations",
    "kernel_balance_reltol",
    "stabilize_weights",
    "max_weight"
  )

  if (is.list(spec) && !is.null(names(spec)) && any(names(spec) %in% known_names)) {
    for (argument_name in intersect(names(spec), known_names)) {
      scenario[[argument_name]] <- spec[[argument_name]]
    }
  } else {
    scenario$weights <- spec
    scenario$weight_method <- "external"
  }

  scenario$weight_method <- match.arg(scenario$weight_method, c("external", "gaussian_gps", "kernel_gps", "cbps", "entropy_balance", "kernel_balance", "adaptive_balance"))
  scenario$gps_model <- match.arg(scenario$gps_model, c("linear", "spline", "gam", "tree", "forest", "stack", "boost"))
  scenario$gps_df <- as.integer(scenario$gps_df %||% 4L)
  resolve_gps_stack_models <- get("cdmc_resolve_gps_stack_models", mode = "function")
  scenario$gps_stack_models <- if (identical(scenario$gps_model, "stack")) {
    resolve_gps_stack_models(scenario$gps_stack_models %||% NULL)
  } else {
    NULL
  }
  resolve_max_weight_spec <- get("cdmc_resolve_max_weight_spec", mode = "function")
  resolve_adaptive_balance_methods <- get("cdmc_resolve_adaptive_balance_methods", mode = "function")
  resolve_entropy_balance_degree <- get("cdmc_resolve_entropy_balance_degree", mode = "function")
  resolve_entropy_balance_iterations <- get("cdmc_resolve_entropy_balance_iterations", mode = "function")
  resolve_entropy_balance_reltol <- get("cdmc_resolve_entropy_balance_reltol", mode = "function")
  resolve_kernel_balance_centers <- get("cdmc_resolve_kernel_balance_centers", mode = "function")
  resolve_kernel_balance_bandwidth <- get("cdmc_resolve_kernel_balance_bandwidth", mode = "function")
  scenario$adaptive_balance_methods <- if (identical(scenario$weight_method, "adaptive_balance")) {
    resolve_adaptive_balance_methods(scenario$adaptive_balance_methods %||% NULL)
  } else {
    NULL
  }
  scenario$max_weight <- resolve_max_weight_spec(scenario$max_weight)
  scenario$entropy_balance_degree <- resolve_entropy_balance_degree(scenario$entropy_balance_degree %||% 1L)
  scenario$entropy_balance_iterations <- resolve_entropy_balance_iterations(scenario$entropy_balance_iterations %||% 1000L)
  scenario$entropy_balance_reltol <- resolve_entropy_balance_reltol(scenario$entropy_balance_reltol %||% 1e-8)
  scenario$kernel_balance_degree <- resolve_entropy_balance_degree(scenario$kernel_balance_degree %||% 1L)
  scenario$kernel_balance_centers <- resolve_kernel_balance_centers(scenario$kernel_balance_centers %||% 25L)
  scenario$kernel_balance_bandwidth <- resolve_kernel_balance_bandwidth(scenario$kernel_balance_bandwidth %||% NULL)
  scenario$kernel_balance_iterations <- resolve_entropy_balance_iterations(scenario$kernel_balance_iterations %||% 1000L)
  scenario$kernel_balance_reltol <- resolve_entropy_balance_reltol(scenario$kernel_balance_reltol %||% 1e-8)
  if (!is.logical(scenario$entropy_balance_standardize) || length(scenario$entropy_balance_standardize) != 1L || is.na(scenario$entropy_balance_standardize)) {
    stop("entropy_balance_standardize must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(scenario$kernel_balance_standardize) || length(scenario$kernel_balance_standardize) != 1L || is.na(scenario$kernel_balance_standardize)) {
    stop("kernel_balance_standardize must be TRUE or FALSE.", call. = FALSE)
  }
  if (identical(scenario$gps_model, "boost")) {
    get("cdmc_assert_installed", mode = "function")("gbm")
  }
  if ((identical(scenario$gps_model, "gam") ||
       (identical(scenario$gps_model, "stack") && "gam" %in% scenario$gps_stack_models)) &&
      scenario$gps_df < 3L) {
    stop("gps_df must be at least 3 when gps_model uses a GAM nuisance learner.", call. = FALSE)
  }
  resolve_gps_bandwidth <- get("cdmc_resolve_gps_bandwidth", mode = "function")
  scenario$gps_bandwidth <- resolve_gps_bandwidth(scenario$gps_bandwidth %||% NULL)
  resolve_gps_forest_trees <- get("cdmc_resolve_gps_forest_trees", mode = "function")
  resolve_gps_forest_mtry <- get("cdmc_resolve_gps_forest_mtry", mode = "function")
  resolve_gps_forest_min_node_size <- get("cdmc_resolve_gps_forest_min_node_size", mode = "function")
  resolve_gps_boost_trees <- get("cdmc_resolve_gps_boost_trees", mode = "function")
  resolve_gps_boost_depth <- get("cdmc_resolve_gps_boost_depth", mode = "function")
  resolve_gps_boost_shrinkage <- get("cdmc_resolve_gps_boost_shrinkage", mode = "function")
  resolve_gps_boost_min_obs_node <- get("cdmc_resolve_gps_boost_min_obs_node", mode = "function")
  scenario$gps_forest_trees <- resolve_gps_forest_trees(scenario$gps_forest_trees %||% 200L)
  scenario$gps_forest_mtry <- resolve_gps_forest_mtry(scenario$gps_forest_mtry %||% NULL)
  scenario$gps_forest_min_node_size <- resolve_gps_forest_min_node_size(scenario$gps_forest_min_node_size %||% NULL)
  scenario$gps_boost_trees <- resolve_gps_boost_trees(scenario$gps_boost_trees %||% 200L)
  scenario$gps_boost_depth <- resolve_gps_boost_depth(scenario$gps_boost_depth %||% 2L)
  scenario$gps_boost_shrinkage <- resolve_gps_boost_shrinkage(scenario$gps_boost_shrinkage %||% 0.05)
  scenario$gps_boost_min_obs_node <- resolve_gps_boost_min_obs_node(scenario$gps_boost_min_obs_node %||% 10L)

  if (scenario$weight_method %in% c("cbps", "entropy_balance", "kernel_balance", "adaptive_balance") && !scenario$gps_model %in% c("linear", "spline")) {
    stop(
      "gps_model is not available for balancing-only replay scenarios.",
      call. = FALSE
    )
  }

  if (identical(scenario$weight_method, "external") && is.null(scenario$weights)) {
    stop(
      "External DR replay scenarios must provide weights through weight_specs.",
      call. = FALSE
    )
  }

  if (!identical(scenario$weight_method, "external")) {
    scenario$weights <- NULL
  }

  scenario$weight_mode <- scenario$weight_method
  scenario
}

cdmc_resolve_sensitivity_weight_specs <- function(object, source_object, weight_specs = NULL, include_unweighted = NULL) {
  if (!is.null(include_unweighted) && (!is.logical(include_unweighted) || length(include_unweighted) != 1L)) {
    stop("include_unweighted must be TRUE, FALSE, or NULL.", call. = FALSE)
  }

  named_specs <- cdmc_sensitivity_named_weight_specs(weight_specs)
  source_class <- cdmc_sensitivity_source_class(source_object)

  if (identical(source_class, "cdmc_fit")) {
    scenarios <- list(current = cdmc_sensitivity_fit_weight_scenario(source_object$weights %||% NULL))

    if (is.null(include_unweighted)) {
      include_unweighted <- !is.null(source_object$weights)
    }
    if (isTRUE(include_unweighted) && !is.null(source_object$weights)) {
      scenarios$unweighted <- cdmc_sensitivity_fit_weight_scenario(NULL)
    }

    if (is.null(named_specs)) {
      return(scenarios)
    }

    duplicated_names <- intersect(names(scenarios), names(named_specs))
    if (length(duplicated_names) > 0L) {
      stop(
        sprintf("weight_specs contains names that are already reserved: %s.", paste(duplicated_names, collapse = ", ")),
        call. = FALSE
      )
    }

    for (scenario_name in names(named_specs)) {
      scenarios[[scenario_name]] <- cdmc_sensitivity_fit_weight_scenario(
        cdmc_sensitivity_fit_weight_value(named_specs[[scenario_name]])
      )
    }

    return(scenarios)
  }

  scenarios <- list(current = cdmc_sensitivity_default_dr_weight_scenario(source_object))

  if (is.null(include_unweighted)) {
    include_unweighted <- FALSE
  }
  if (isTRUE(include_unweighted)) {
    if (!identical(source_object$weight_method, "external")) {
      stop(
        "include_unweighted is only available for cdmc_dr_fit objects that currently use external weights.",
        call. = FALSE
      )
    }
    scenarios$gaussian_gps <- cdmc_sensitivity_parse_dr_weight_spec(
      source_object,
      list(weight_method = "gaussian_gps")
    )
  }

  if (is.null(named_specs)) {
    return(scenarios)
  }

  duplicated_names <- intersect(names(scenarios), names(named_specs))
  if (length(duplicated_names) > 0L) {
    stop(
      sprintf("weight_specs contains names that are already reserved: %s.", paste(duplicated_names, collapse = ", ")),
      call. = FALSE
    )
  }

  for (scenario_name in names(named_specs)) {
    scenarios[[scenario_name]] <- cdmc_sensitivity_parse_dr_weight_spec(source_object, named_specs[[scenario_name]])
  }

  scenarios
}

cdmc_sensitivity_original_columns <- function(object) {
  cdmc_bootstrap_original_columns(object)
}

cdmc_sensitivity_weight_mode <- function(source_object, weight_scenario) {
  if (!is.null(weight_scenario$weight_mode)) {
    return(weight_scenario$weight_mode)
  }

  if (inherits(source_object, "cdmc_dr_fit")) {
    return(weight_scenario$weight_method %||% source_object$weight_method)
  }

  if (is.null(weight_scenario$weights)) "unweighted" else "weighted"
}

cdmc_sensitivity_support_summary <- function(object, data, weight_scenario, washout, zero_tolerance) {
  prepared <- cdmc_prepare_panel(
    data = data,
    outcome = object$outcome,
    dose = object$dose,
    unit = object$unit,
    time = object$time,
    covariates = object$covariates,
    zero_tolerance = zero_tolerance
  )

  eligible_mask <- cdmc_build_eligible_mask(
    dose_matrix = prepared$dose_matrix,
    zero_tolerance = zero_tolerance,
    washout = washout
  )
  objective <- cdmc_sensitivity_fit_objective(object)
  optimization_sample <- cdmc_sensitivity_optimization_sample(object)
  observed_controls <- sum(cdmc_zero_dose_mask(prepared$dose_matrix, zero_tolerance = zero_tolerance))
  observed_optimization_cells <- if (identical(optimization_sample, "observed_panel")) {
    sum(prepared$observed_mask)
  } else {
    observed_controls
  }
  optimization_mask <- if (identical(optimization_sample, "observed_panel")) {
    prepared$observed_mask
  } else {
    eligible_mask
  }
  weight_mode <- cdmc_sensitivity_weight_mode(object, weight_scenario)

  if (!(inherits(object, "cdmc_dr_fit") && identical(weight_mode, "gaussian_gps"))) {
    weight_info <- cdmc_prepare_panel_weights(
      weights = weight_scenario$weights,
      data = prepared$data,
      n_units = prepared$n_units,
      n_times = prepared$n_times,
      eligible_mask = optimization_mask
    )

    if (weight_info$supplied) {
      eligible_mask <- eligible_mask & weight_info$matrix > 0
      optimization_mask <- optimization_mask & weight_info$matrix > 0
    }
  }

  control_support_error <- tryCatch(
    {
      cdmc_validate_support(eligible_mask)
      NULL
    },
    error = function(error) conditionMessage(error)
  )

  optimization_support_error <- tryCatch(
    {
      if (inherits(object, "cdmc_fit") && identical(objective, "joint")) {
        validate_joint_support <- get("cdmc_validate_joint_support", mode = "function")
        validate_joint_support(optimization_mask)
      } else {
        cdmc_validate_support(optimization_mask)
      }
      NULL
    },
    error = function(error) conditionMessage(error)
  )

  list(
    objective = objective,
    optimization_sample = optimization_sample,
    observed_controls = observed_controls,
    eligible_controls = sum(eligible_mask),
    min_unit_support = min(rowSums(eligible_mask)),
    min_time_support = min(colSums(eligible_mask)),
    control_support_ok = is.null(control_support_error),
    control_support_error = control_support_error,
    observed_optimization_cells = observed_optimization_cells,
    optimization_cells = sum(optimization_mask),
    min_unit_optimization_support = min(rowSums(optimization_mask)),
    min_time_optimization_support = min(colSums(optimization_mask)),
    support_ok = is.null(optimization_support_error),
    support_error = optimization_support_error,
    weight_mode = weight_mode
  )
}

cdmc_sensitivity_resolve_statistics <- function(object, statistics = NULL) {
  if (is.null(statistics)) {
    return(cdmc_default_bootstrap_statistics(object))
  }

  unique(match.arg(
    statistics,
    cdmc_supported_bootstrap_statistics(object),
    several.ok = TRUE
  ))
}

cdmc_sensitivity_target_statistics <- function(
  object,
  statistics,
  contrast_history = NULL,
  reference_history = NULL,
  prediction_dose = NULL,
  prediction_history = NULL,
  prediction_type = c("response", "slope")
) {
  prediction_type <- match.arg(prediction_type)
  object_class <- cdmc_sensitivity_target_class(object)

  if (identical(object_class, "cdmc_dose_response") &&
      "prediction" %in% statistics &&
      is.null(prediction_dose) &&
      is.null(prediction_history)) {
    stop(
      "Provide prediction_dose or prediction_history when requesting the 'prediction' sensitivity statistic for cdmc_dose_response objects.",
      call. = FALSE
    )
  }

  out <- list(
    fit_success = TRUE,
    fit_error = NA_character_
  )

  if (identical(object_class, "cdmc_fit")) {
    out$lambda <- object$lambda
    out$effective_rank <- object$baseline$effective_rank
    out$baseline_converged <- object$baseline$converged
  } else if (identical(object_class, "cdmc_dr_fit")) {
    out$lambda <- object$lambda
    out$dr_sample_cells <- sum(object$effect$sample_mask, na.rm = TRUE)
    out$n_folds <- object$n_folds
  } else if (identical(object_class, "cdmc_dose_response")) {
    out$dose_response_sample_size <- nrow(object$history)
  } else if (identical(object_class, "cdmc_dynamic_estimand")) {
    out$dynamic_path_count <- nrow(object$estimate_table)
  }

  collected <- cdmc_collect_bootstrap_statistics(
    object,
    statistics = statistics,
    contrast_history = contrast_history,
    reference_history = reference_history,
    prediction_dose = prediction_dose,
    prediction_history = prediction_history,
    prediction_type = prediction_type
  )

  for (statistic_name in names(collected)) {
    out[[statistic_name]] <- unname(collected[[statistic_name]])
  }

  out
}

cdmc_sensitivity_reference_summary <- function(
  object,
  source_object,
  data,
  weight_scenario,
  statistics,
  contrast_history = NULL,
  reference_history = NULL,
  prediction_dose = NULL,
  prediction_history = NULL,
  prediction_type = c("response", "slope")
) {
  prediction_type <- match.arg(prediction_type)
  source_class <- cdmc_sensitivity_source_class(source_object)
  target_class <- cdmc_sensitivity_target_class(object)
  support <- cdmc_sensitivity_support_summary(
    object = source_object,
    data = data,
    weight_scenario = weight_scenario,
    washout = source_object$washout,
    zero_tolerance = source_object$zero_tolerance
  )

  reference <- c(
    list(
      source_class = source_class,
      target_class = target_class,
      objective = support$objective,
      optimization_sample = support$optimization_sample,
      weight_scenario = "current",
      weight_mode = support$weight_mode,
      weighted_fit = !identical(support$weight_mode, "unweighted"),
      washout = source_object$washout,
      zero_tolerance = source_object$zero_tolerance,
      observed_controls = support$observed_controls,
      eligible_controls = support$eligible_controls,
      control_support_ok = support$control_support_ok,
      observed_optimization_cells = support$observed_optimization_cells,
      optimization_cells = support$optimization_cells,
      min_unit_support = support$min_unit_support,
      min_time_support = support$min_time_support,
      min_unit_optimization_support = support$min_unit_optimization_support,
      min_time_optimization_support = support$min_time_optimization_support,
      support_ok = support$support_ok
    ),
    cdmc_sensitivity_target_statistics(
      object,
      statistics = statistics,
      contrast_history = contrast_history,
      reference_history = reference_history,
      prediction_dose = prediction_dose,
      prediction_history = prediction_history,
      prediction_type = prediction_type
    )
  )

  reference$support_fraction <- cdmc_sensitivity_support_fraction(
    observed_controls = reference$observed_controls,
    eligible_controls = reference$eligible_controls
  )
  reference$optimization_fraction <- cdmc_sensitivity_support_fraction(
    observed_controls = reference$observed_optimization_cells,
    eligible_controls = reference$optimization_cells
  )
  reference$support_loss <- if (is.finite(reference$observed_controls) && is.finite(reference$eligible_controls)) {
    reference$observed_controls - reference$eligible_controls
  } else {
    NA_real_
  }
  reference$optimization_loss <- if (is.finite(reference$observed_optimization_cells) && is.finite(reference$optimization_cells)) {
    reference$observed_optimization_cells - reference$optimization_cells
  } else {
    NA_real_
  }

  reference
}

cdmc_sensitivity_non_metric_columns <- function() {
  c(
    "source_class",
    "target_class",
    "objective",
    "optimization_sample",
    "weight_scenario",
    "weight_mode",
    "weighted_fit",
    "washout",
    "zero_tolerance",
    "observed_controls",
    "observed_optimization_cells",
    "eligible_controls",
    "min_unit_support",
    "min_time_support",
    "min_unit_optimization_support",
    "min_time_optimization_support",
    "control_support_ok",
    "control_support_error",
    "support_ok",
    "fit_success",
    "fit_error",
    "support_fraction",
    "support_loss",
    "optimization_loss"
  )
}

cdmc_sensitivity_target_stat_columns <- function(data) {
  candidate_columns <- setdiff(
    names(data),
    c(cdmc_sensitivity_non_metric_columns(), grep("^delta_", names(data), value = TRUE))
  )

  candidate_columns[vapply(candidate_columns, function(column_name) {
    is.numeric(data[[column_name]]) || is.integer(data[[column_name]])
  }, logical(1))]
}

cdmc_sensitivity_augment_results <- function(results, reference) {
  results$support_fraction <- vapply(seq_len(nrow(results)), function(index) {
    cdmc_sensitivity_support_fraction(
      observed_controls = results$observed_controls[[index]],
      eligible_controls = results$eligible_controls[[index]]
    )
  }, numeric(1))
  results$optimization_fraction <- vapply(seq_len(nrow(results)), function(index) {
    cdmc_sensitivity_support_fraction(
      observed_controls = results$observed_optimization_cells[[index]],
      eligible_controls = results$optimization_cells[[index]]
    )
  }, numeric(1))
  results$support_loss <- vapply(seq_len(nrow(results)), function(index) {
    observed_controls <- results$observed_controls[[index]]
    eligible_controls <- results$eligible_controls[[index]]
    if (!is.finite(observed_controls) || !is.finite(eligible_controls)) {
      return(NA_real_)
    }

    observed_controls - eligible_controls
  }, numeric(1))
  results$optimization_loss <- vapply(seq_len(nrow(results)), function(index) {
    observed_cells <- results$observed_optimization_cells[[index]]
    optimization_cells <- results$optimization_cells[[index]]
    if (!is.finite(observed_cells) || !is.finite(optimization_cells)) {
      return(NA_real_)
    }

    observed_cells - optimization_cells
  }, numeric(1))

  delta_columns <- intersect(
    c(
      "eligible_controls",
      "support_fraction",
      "optimization_cells",
      "optimization_fraction",
      cdmc_sensitivity_target_stat_columns(results)
    ),
    names(reference)
  )
  for (column_name in delta_columns) {
    results[[paste0("delta_", column_name)]] <- vapply(seq_len(nrow(results)), function(index) {
      value <- results[[column_name]][[index]]
      reference_value <- reference[[column_name]]
      if (!is.finite(value) || !is.finite(reference_value)) {
        return(NA_real_)
      }

      value - reference_value
    }, numeric(1))
  }

  results
}

cdmc_sensitivity_support_table <- function(results) {
  summary_columns <- intersect(
    c(
      "source_class",
      "target_class",
      "objective",
      "optimization_sample",
      "weight_scenario",
      "weight_mode",
      "weighted_fit",
      "washout",
      "zero_tolerance",
      "observed_controls",
      "eligible_controls",
      "control_support_ok",
      "observed_optimization_cells",
      "optimization_cells",
      "optimization_fraction",
      "optimization_loss",
      "support_fraction",
      "support_loss",
      "min_unit_support",
      "min_time_support",
      "min_unit_optimization_support",
      "min_time_optimization_support",
      "delta_eligible_controls",
      "delta_support_fraction",
      "delta_optimization_cells",
      "delta_optimization_fraction",
      "support_ok",
      "fit_success",
      "fit_error"
    ),
    names(results)
  )

  results[, summary_columns, drop = FALSE]
}

cdmc_sensitivity_metric_range <- function(values) {
  finite_values <- values[is.finite(values)]
  if (length(finite_values) == 0L) {
    return(c(min = NA_real_, max = NA_real_, range = NA_real_))
  }

  c(
    min = min(finite_values),
    max = max(finite_values),
    range = max(finite_values) - min(finite_values)
  )
}

cdmc_sensitivity_group_mask <- function(data, row, columns) {
  mask <- rep(TRUE, nrow(data))
  for (column_name in columns) {
    mask <- mask & data[[column_name]] == row[[column_name]]
  }

  mask
}

cdmc_sensitivity_dependence_summary <- function(results, group_columns, varying_column) {
  grouped_values <- unique(results[, group_columns, drop = FALSE])
  metric_columns <- intersect(
    c(
      "eligible_controls",
      "support_fraction",
      "optimization_cells",
      "optimization_fraction",
      cdmc_sensitivity_target_stat_columns(results)
    ),
    names(results)
  )
  delta_columns <- intersect(
    c(
      "delta_eligible_controls",
      "delta_support_fraction",
      "delta_optimization_cells",
      "delta_optimization_fraction",
      paste0("delta_", cdmc_sensitivity_target_stat_columns(results))
    ),
    names(results)
  )

  count_name <- if (identical(varying_column, "washout")) "n_washout" else "n_weight_scenarios"
  success_name <- if (identical(varying_column, "washout")) "successful_washout" else "successful_weight_scenarios"
  values_name <- if (identical(varying_column, "washout")) "washout_values" else "weight_scenarios"

  rows <- vector("list", nrow(grouped_values))
  for (row_index in seq_len(nrow(grouped_values))) {
    grouping_row <- grouped_values[row_index, , drop = FALSE]
    group_rows <- results[cdmc_sensitivity_group_mask(results, grouping_row, group_columns), , drop = FALSE]
    successful_rows <- group_rows[!is.na(group_rows$fit_success) & group_rows$fit_success, , drop = FALSE]

    varying_values <- unique(group_rows[[varying_column]])
    if (is.numeric(varying_values) || is.integer(varying_values)) {
      varying_values <- sort(varying_values)
    } else {
      varying_values <- sort(as.character(varying_values))
    }

    row <- as.list(grouping_row)
    row[[values_name]] <- paste(varying_values, collapse = ", ")
    row[[count_name]] <- length(varying_values)
    row[[success_name]] <- nrow(successful_rows)

    for (column_name in metric_columns) {
      range_summary <- cdmc_sensitivity_metric_range(successful_rows[[column_name]])
      for (suffix in names(range_summary)) {
        row[[paste0(column_name, "_", suffix)]] <- range_summary[[suffix]]
      }
    }

    for (column_name in delta_columns) {
      successful_deltas <- abs(successful_rows[[column_name]])
      successful_deltas <- successful_deltas[is.finite(successful_deltas)]
      row[[paste0("max_abs_", column_name)]] <- if (length(successful_deltas) > 0L) {
        max(successful_deltas)
      } else {
        NA_real_
      }
    }

    rows[[row_index]] <- row
  }

  cdmc_sensitivity_bind_rows(rows)
}

cdmc_sensitivity_bind_rows <- function(rows) {
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  row_data <- lapply(rows, function(row) {
    missing_names <- setdiff(all_names, names(row))
    for (missing_name in missing_names) {
      row[[missing_name]] <- NA
    }
    row[all_names]
  })

  do.call(
    rbind,
    lapply(row_data, function(row) {
      as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
    })
  )
}

cdmc_sensitivity_fit_exclusions <- function(source_object) {
  if (inherits(source_object, "cdmc_fit")) {
    return(c("weights", "washout", "zero_tolerance", "verbose"))
  }

  c(
    "weights",
    "weight_method",
    "weight_covariates",
    "gps_time_effects",
    "gps_model",
    "gps_df",
    "gps_spline_covariates",
    "gps_stack_models",
    "gps_bandwidth",
    "gps_forest_trees",
    "gps_forest_mtry",
    "gps_forest_min_node_size",
    "gps_boost_trees",
    "gps_boost_depth",
    "gps_boost_shrinkage",
    "gps_boost_min_obs_node",
    "cbps_standardize",
    "cbps_method",
    "cbps_iterations",
    "cbps_twostep",
    "adaptive_balance_methods",
    "entropy_balance_degree",
    "entropy_balance_standardize",
    "entropy_balance_iterations",
    "entropy_balance_reltol",
    "kernel_balance_degree",
    "kernel_balance_centers",
    "kernel_balance_bandwidth",
    "kernel_balance_standardize",
    "kernel_balance_iterations",
    "kernel_balance_reltol",
    "stabilize_weights",
    "max_weight",
    "washout",
    "zero_tolerance",
    "verbose"
  )
}

cdmc_sensitivity_refit_source <- function(source_object, data, fit_spec, weight_scenario, washout, zero_tolerance) {
  common_arguments <- list(
    data = data,
    outcome = source_object$outcome,
    dose = source_object$dose,
    unit = source_object$unit,
    time = source_object$time,
    covariates = source_object$covariates,
    washout = washout,
    zero_tolerance = zero_tolerance
  )

  if (inherits(source_object, "cdmc_fit")) {
    return(do.call(
      cdmc_fit,
      c(
        common_arguments,
        list(weights = weight_scenario$weights),
        fit_spec$fit_spec[!names(fit_spec$fit_spec) %in% cdmc_sensitivity_fit_exclusions(source_object)],
        list(verbose = FALSE)
      )
    ))
  }

  do.call(
    cdmc_dr_fit,
    c(
      common_arguments,
      list(
        weights = weight_scenario$weights,
        weight_method = weight_scenario$weight_method,
        weight_covariates = weight_scenario$weight_covariates,
        gps_time_effects = weight_scenario$gps_time_effects,
        gps_model = weight_scenario$gps_model,
        gps_df = weight_scenario$gps_df,
        gps_spline_covariates = weight_scenario$gps_spline_covariates,
        gps_stack_models = weight_scenario$gps_stack_models,
        gps_bandwidth = weight_scenario$gps_bandwidth,
        gps_forest_trees = weight_scenario$gps_forest_trees,
        gps_forest_mtry = weight_scenario$gps_forest_mtry,
        gps_forest_min_node_size = weight_scenario$gps_forest_min_node_size,
        gps_boost_trees = weight_scenario$gps_boost_trees,
        gps_boost_depth = weight_scenario$gps_boost_depth,
        gps_boost_shrinkage = weight_scenario$gps_boost_shrinkage,
        gps_boost_min_obs_node = weight_scenario$gps_boost_min_obs_node,
        cbps_standardize = weight_scenario$cbps_standardize,
        cbps_method = weight_scenario$cbps_method,
        cbps_iterations = weight_scenario$cbps_iterations,
        cbps_twostep = weight_scenario$cbps_twostep,
        adaptive_balance_methods = weight_scenario$adaptive_balance_methods,
        entropy_balance_degree = weight_scenario$entropy_balance_degree,
        entropy_balance_standardize = weight_scenario$entropy_balance_standardize,
        entropy_balance_iterations = weight_scenario$entropy_balance_iterations,
        entropy_balance_reltol = weight_scenario$entropy_balance_reltol,
        kernel_balance_degree = weight_scenario$kernel_balance_degree,
        kernel_balance_centers = weight_scenario$kernel_balance_centers,
        kernel_balance_bandwidth = weight_scenario$kernel_balance_bandwidth,
        kernel_balance_standardize = weight_scenario$kernel_balance_standardize,
        kernel_balance_iterations = weight_scenario$kernel_balance_iterations,
        kernel_balance_reltol = weight_scenario$kernel_balance_reltol,
        stabilize_weights = weight_scenario$stabilize_weights,
        max_weight = weight_scenario$max_weight,
        fold_assignments = source_object$fold_assignments %||% NULL
      ),
      fit_spec$fit_spec[!names(fit_spec$fit_spec) %in% cdmc_sensitivity_fit_exclusions(source_object)],
      list(verbose = FALSE)
    )
  )
}

cdmc_sensitivity_pick_pattern_columns <- function(data, patterns, limit) {
  out <- character(0)
  for (pattern in patterns) {
    matches <- grep(pattern, names(data), value = TRUE)
    out <- c(out, setdiff(matches, out))
    if (length(out) >= limit) {
      break
    }
  }

  out[seq_len(min(length(out), limit))]
}

cdmc_sensitivity_display_range_columns <- function(data) {
  base_columns <- intersect(
    c(
      "eligible_controls_range",
      "support_fraction_range",
      "optimization_cells_range",
      "optimization_fraction_range"
    ),
    names(data)
  )
  target_columns <- cdmc_sensitivity_pick_pattern_columns(
    data,
    patterns = c(
      "^mean_tau_active_dose_range$",
      "^mean_tau_dr_range$",
      "^mean_tau_dr_linear_range$",
      "^dose_response_.*_range$",
      "^mean_dynamic_.*_range$",
      "^dynamic_.*_range$",
      "^coef_.*_range$"
    ),
    limit = 3L
  )

  c(base_columns, target_columns)
}

cdmc_sensitivity_display_delta_columns <- function(data) {
  cdmc_sensitivity_pick_pattern_columns(
    data,
    patterns = c(
      "^max_abs_delta_optimization_cells$",
      "^max_abs_delta_optimization_fraction$",
      "^max_abs_delta_mean_tau_active_dose$",
      "^max_abs_delta_mean_tau_dr$",
      "^max_abs_delta_mean_tau_dr_linear$",
      "^max_abs_delta_dose_response_.*$",
      "^max_abs_delta_mean_dynamic_.*$",
      "^max_abs_delta_dynamic_.*$",
      "^max_abs_delta_coef_.*$"
    ),
    limit = 2L
  )
}

cdmc_sensitivity_scan <- function(
  object,
  washout_grid = cdmc_default_washout_grid(object),
  zero_tolerance_grid = cdmc_default_zero_tolerance_grid(object),
  weight_specs = NULL,
  include_unweighted = NULL,
  statistics = NULL,
  rerun_tuning = NULL,
  contrast_history = NULL,
  reference_history = NULL,
  prediction_dose = NULL,
  prediction_history = NULL,
  prediction_type = c("response", "slope"),
  workers = 1L,
  verbose = FALSE
) {
  object_class <- cdmc_sensitivity_target_class(object)

  prediction_type <- match.arg(prediction_type)

  if (!identical(object_class, "cdmc_dr_fit") && (!is.null(contrast_history) || !is.null(reference_history))) {
    stop("contrast_history and reference_history are only supported for cdmc_dr_fit sensitivity summaries.", call. = FALSE)
  }

  if (!identical(object_class, "cdmc_dose_response") && (!is.null(prediction_dose) || !is.null(prediction_history))) {
    stop("prediction_dose and prediction_history are only supported for cdmc_dose_response sensitivity summaries.", call. = FALSE)
  }

  washout_grid <- unique(as.integer(washout_grid))
  if (length(washout_grid) < 1L || any(!is.finite(washout_grid)) || any(washout_grid < 0L)) {
    stop("washout_grid must contain one or more nonnegative integers.", call. = FALSE)
  }

  zero_tolerance_grid <- unique(as.numeric(zero_tolerance_grid))
  if (length(zero_tolerance_grid) < 1L || any(!is.finite(zero_tolerance_grid)) || any(zero_tolerance_grid < 0)) {
    stop("zero_tolerance_grid must contain one or more nonnegative numeric values.", call. = FALSE)
  }

  if (!is.numeric(workers) || length(workers) != 1L || !is.finite(workers) || workers < 1 || workers != floor(workers)) {
    stop("workers must be a positive integer.", call. = FALSE)
  }
  workers <- as.integer(workers)

  source_object <- cdmc_sensitivity_source_object(object)
  source_class <- cdmc_sensitivity_source_class(source_object)
  statistics <- cdmc_sensitivity_resolve_statistics(object, statistics = statistics)
  weight_scenarios <- cdmc_resolve_sensitivity_weight_specs(
    object = object,
    source_object = source_object,
    weight_specs = weight_specs,
    include_unweighted = include_unweighted
  )
  fit_spec <- cdmc_build_bootstrap_fit_spec(object, rerun_tuning = rerun_tuning)

  fit_data <- source_object$data
  if (".cdmc_observed" %in% names(fit_data)) {
    fit_data <- fit_data[!is.na(fit_data$.cdmc_observed) & fit_data$.cdmc_observed, , drop = FALSE]
  }
  fit_data <- fit_data[, cdmc_sensitivity_original_columns(object), drop = FALSE]

  reference <- cdmc_sensitivity_reference_summary(
    object = object,
    source_object = source_object,
    data = fit_data,
    weight_scenario = weight_scenarios$current,
    statistics = statistics,
    contrast_history = contrast_history,
    reference_history = reference_history,
    prediction_dose = prediction_dose,
    prediction_history = prediction_history,
    prediction_type = prediction_type
  )

  scenario_grid <- expand.grid(
    weight_scenario = names(weight_scenarios),
    washout = washout_grid,
    zero_tolerance = zero_tolerance_grid,
    stringsAsFactors = FALSE
  )
  use_parallel <- workers > 1L
  if (use_parallel && identical(.Platform$OS.type, "windows")) {
    warning(
      "Parallel sensitivity scan currently uses multicore execution and is not available on Windows. Falling back to sequential execution.",
      call. = FALSE
    )
    use_parallel <- FALSE
    workers <- 1L
  }
  worker_count <- if (use_parallel) min(workers, nrow(scenario_grid)) else 1L

  run_sensitivity_scenario <- function(scenario_index) {
    weight_name <- scenario_grid$weight_scenario[[scenario_index]]
    washout <- scenario_grid$washout[[scenario_index]]
    zero_tolerance <- scenario_grid$zero_tolerance[[scenario_index]]
    weight_scenario <- weight_scenarios[[weight_name]]

    if (verbose) {
      message(sprintf(
        "sensitivity scenario %d/%d: weights=%s, washout=%d, zero_tolerance=%.6g",
        scenario_index,
        nrow(scenario_grid),
        weight_name,
        washout,
        zero_tolerance
      ))
    }

    support <- tryCatch(
      cdmc_sensitivity_support_summary(
        object = source_object,
        data = fit_data,
        weight_scenario = weight_scenario,
        washout = washout,
        zero_tolerance = zero_tolerance
      ),
      error = function(error) list(
        objective = cdmc_sensitivity_fit_objective(source_object),
        optimization_sample = cdmc_sensitivity_optimization_sample(source_object),
        observed_controls = NA_integer_,
        eligible_controls = NA_integer_,
        control_support_ok = NA,
        observed_optimization_cells = NA_integer_,
        optimization_cells = NA_integer_,
        min_unit_support = NA_integer_,
        min_time_support = NA_integer_,
        min_unit_optimization_support = NA_integer_,
        min_time_optimization_support = NA_integer_,
        support_ok = FALSE,
        support_error = conditionMessage(error),
        weight_mode = cdmc_sensitivity_weight_mode(source_object, weight_scenario)
      )
    )

    row <- list(
      source_class = source_class,
      target_class = object_class,
      objective = support$objective,
      optimization_sample = support$optimization_sample,
      weight_scenario = weight_name,
      weight_mode = support$weight_mode,
      weighted_fit = !identical(support$weight_mode, "unweighted"),
      washout = washout,
      zero_tolerance = zero_tolerance,
      observed_controls = support$observed_controls,
      eligible_controls = support$eligible_controls,
      control_support_ok = support$control_support_ok,
      observed_optimization_cells = support$observed_optimization_cells,
      optimization_cells = support$optimization_cells,
      min_unit_support = support$min_unit_support,
      min_time_support = support$min_time_support,
      min_unit_optimization_support = support$min_unit_optimization_support,
      min_time_optimization_support = support$min_time_optimization_support,
      support_ok = support$support_ok
    )

    if (!isTRUE(support$support_ok)) {
      row$fit_success <- FALSE
      row$fit_error <- support$support_error
      return(list(index = scenario_index, row = row))
    }

    current_fit <- tryCatch(
      cdmc_sensitivity_refit_source(
        source_object = source_object,
        data = fit_data,
        fit_spec = fit_spec,
        weight_scenario = weight_scenario,
        washout = washout,
        zero_tolerance = zero_tolerance
      ),
      error = function(error) error
    )

    if (inherits(current_fit, "error")) {
      row$fit_success <- FALSE
      row$fit_error <- conditionMessage(current_fit)
      return(list(index = scenario_index, row = row))
    }

    target_object <- tryCatch(
      cdmc_bootstrap_evaluate_target(object, current_fit),
      error = function(error) error
    )

    if (inherits(target_object, "error")) {
      row$fit_success <- FALSE
      row$fit_error <- conditionMessage(target_object)
      return(list(index = scenario_index, row = row))
    }

    target_statistics <- tryCatch(
      cdmc_sensitivity_target_statistics(
        object = target_object,
        statistics = statistics,
        contrast_history = contrast_history,
        reference_history = reference_history,
        prediction_dose = prediction_dose,
        prediction_history = prediction_history,
        prediction_type = prediction_type
      ),
      error = function(error) error
    )

    if (inherits(target_statistics, "error")) {
      row$fit_success <- FALSE
      row$fit_error <- conditionMessage(target_statistics)
      return(list(index = scenario_index, row = row))
    }

    list(index = scenario_index, row = c(row, target_statistics))
  }

  if (use_parallel && verbose) {
    message(sprintf("running sensitivity scan in parallel with %d workers", worker_count))
  }

  scenario_results <- if (use_parallel) {
    parallel::mclapply(
      seq_len(nrow(scenario_grid)),
      run_sensitivity_scenario,
      mc.cores = worker_count,
      mc.set.seed = TRUE
    )
  } else {
    lapply(seq_len(nrow(scenario_grid)), run_sensitivity_scenario)
  }

  rows <- vector("list", nrow(scenario_grid))
  for (scenario_result in scenario_results) {
    rows[[scenario_result$index]] <- scenario_result$row
  }

  results <- cdmc_sensitivity_augment_results(cdmc_sensitivity_bind_rows(rows), reference = reference)

  result <- list(
    call = match.call(),
    statistics = statistics,
    rerun_tuning = fit_spec$rerun_tuning,
    workers = worker_count,
    parallel = use_parallel,
    results = results,
    reference = reference,
    support_summary = cdmc_sensitivity_support_table(results),
    washout_summary = cdmc_sensitivity_dependence_summary(
      results = results,
      group_columns = c("objective", "optimization_sample", "weight_scenario", "weight_mode", "zero_tolerance"),
      varying_column = "washout"
    ),
    weight_summary = cdmc_sensitivity_dependence_summary(
      results = results,
      group_columns = c("objective", "optimization_sample", "washout", "zero_tolerance"),
      varying_column = "weight_scenario"
    ),
    base_object = object,
    source_object = source_object
  )

  class(result) <- "cdmc_sensitivity_scan"
  result
}

summary.cdmc_sensitivity_scan <- function(object, ...) {
  out <- list(
    reference = object$reference,
    support = object$support_summary,
    washout = object$washout_summary,
    weight = object$weight_summary
  )

  class(out) <- "summary.cdmc_sensitivity_scan"
  out
}

print.summary.cdmc_sensitivity_scan <- function(x, ...) {
  cat("causaldosemc sensitivity summary\n")
  cat(sprintf("  source class: %s\n", x$reference$source_class))
  cat(sprintf("  target class: %s\n", x$reference$target_class))
  cat(sprintf(
    "  reference scenario: objective=%s, optimization=%s, weights=%s (%s), washout=%d, zero_tolerance=%.6g\n",
    x$reference$objective,
    x$reference$optimization_sample,
    x$reference$weight_scenario,
    x$reference$weight_mode,
    x$reference$washout,
    x$reference$zero_tolerance
  ))

  if (nrow(x$support) > 0L) {
    cat("\nSupport summary:\n")
    print(
      x$support[, intersect(
        c(
          "objective",
          "optimization_sample",
          "weight_scenario",
          "weight_mode",
          "washout",
          "zero_tolerance",
          "eligible_controls",
          "optimization_cells",
          "support_fraction",
          "optimization_fraction",
          "control_support_ok",
          "delta_eligible_controls",
          "delta_support_fraction",
          "support_ok",
          "fit_success"
        ),
        names(x$support)
      ), drop = FALSE],
      row.names = FALSE
    )
  }

  if (nrow(x$washout) > 0L) {
    cat("\nWashout dependence:\n")
    washout_columns <- intersect(
      c(
        "objective",
        "optimization_sample",
        "weight_scenario",
        "weight_mode",
        "zero_tolerance",
        "washout_values",
        "n_washout",
        "successful_washout"
      ),
      names(x$washout)
    )
    washout_columns <- c(
      washout_columns,
      setdiff(cdmc_sensitivity_display_range_columns(x$washout), washout_columns),
      setdiff(cdmc_sensitivity_display_delta_columns(x$washout), washout_columns)
    )
    print(x$washout[, washout_columns, drop = FALSE], row.names = FALSE)
  }

  if (nrow(x$weight) > 0L) {
    cat("\nWeight dependence:\n")
    weight_columns <- intersect(
      c(
        "objective",
        "optimization_sample",
        "washout",
        "zero_tolerance",
        "weight_scenarios",
        "n_weight_scenarios",
        "successful_weight_scenarios"
      ),
      names(x$weight)
    )
    weight_columns <- c(
      weight_columns,
      setdiff(cdmc_sensitivity_display_range_columns(x$weight), weight_columns),
      setdiff(cdmc_sensitivity_display_delta_columns(x$weight), weight_columns)
    )
    print(x$weight[, weight_columns, drop = FALSE], row.names = FALSE)
  }

  invisible(x)
}

print.cdmc_sensitivity_scan <- function(x, ...) {
  cat("causaldosemc sensitivity scan\n")
  cat(sprintf("  scenarios evaluated: %d\n", nrow(x$results)))
  cat(sprintf("  successful fits: %d\n", sum(as.logical(x$results$fit_success), na.rm = TRUE)))
  cat(sprintf("  execution mode: %s (%d worker%s)\n", if (isTRUE(x$parallel)) "parallel" else "sequential", x$workers %||% 1L, if ((x$workers %||% 1L) == 1L) "" else "s"))
  cat(sprintf("  rerun tuning: %s\n", if (x$rerun_tuning) "yes" else "no"))

  print(summary(x), ...)
  invisible(x)
}