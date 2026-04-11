cdmc_sensitivity_bounds_operator <- function(design, weights = NULL) {
  design <- as.matrix(design)
  if (nrow(design) < 1L || ncol(design) < 1L) {
    stop("design must have at least one row and one column.", call. = FALSE)
  }

  if (is.null(weights)) {
    qr_fit <- qr(design)
    rhs <- diag(1, nrow = nrow(design), ncol = nrow(design))
  } else {
    sqrt_weights <- sqrt(as.numeric(weights))
    qr_fit <- qr(design * sqrt_weights)
    rhs <- diag(sqrt_weights, nrow = length(sqrt_weights), ncol = length(sqrt_weights))
  }

  operator <- qr.coef(qr_fit, rhs)
  operator <- as.matrix(operator)
  operator[is.na(operator)] <- 0
  rownames(operator) <- colnames(design)
  operator
}

cdmc_sensitivity_bounds_normalized_loading <- function(weights) {
  weights <- as.numeric(weights)
  if (length(weights) < 1L) {
    return(numeric(0))
  }

  total_weight <- sum(weights)
  if (!is.finite(total_weight) || total_weight <= 0) {
    stop("weights must sum to a positive finite value.", call. = FALSE)
  }

  weights / total_weight
}

cdmc_sensitivity_bounds_align_coefficient_rows <- function(rows, coefficient_names) {
  rows <- as.matrix(rows)
  aligned <- matrix(0, nrow = nrow(rows), ncol = length(coefficient_names))
  rownames(aligned) <- rownames(rows)
  colnames(aligned) <- coefficient_names

  common_columns <- intersect(colnames(rows), coefficient_names)
  if (length(common_columns) > 0L) {
    aligned[, common_columns] <- rows[, common_columns, drop = FALSE]
  }

  aligned
}

cdmc_sensitivity_bounds_resolve_perturbation_layer <- function(perturbation_layer) {
  match.arg(perturbation_layer, c("stage2", "baseline", "optimization"))
}

cdmc_sensitivity_bounds_layer_label <- function(perturbation_layer) {
  perturbation_layer <- cdmc_sensitivity_bounds_resolve_perturbation_layer(perturbation_layer)

  if (identical(perturbation_layer, "optimization")) {
    return("optimization_outcome")
  }

  if (identical(perturbation_layer, "baseline")) {
    return("baseline_counterfactual")
  }

  "stage2_response"
}

cdmc_sensitivity_bounds_optimization_source <- function(object) {
  source_object <- cdmc_sensitivity_source_object(object)

  if (!inherits(source_object, "cdmc_fit") && !inherits(source_object, "cdmc_dr_fit")) {
    stop(
      paste(
        "Optimization-layer sensitivity bounds are currently available only for targets backed by 'cdmc_fit' or 'cdmc_dr_fit'.",
        "Objects backed by other source classes are not supported."
      ),
      call. = FALSE
    )
  }

  source_object
}

cdmc_sensitivity_bounds_optimization_mask <- function(object) {
  if (inherits(object, "cdmc_fit")) {
    return(object$optimization_mask %||% object$eligible_mask)
  }

  if (inherits(object, "cdmc_dr_fit")) {
    return(object$observed_mask)
  }

  stop("object must inherit from 'cdmc_fit' or 'cdmc_dr_fit'.", call. = FALSE)
}

cdmc_sensitivity_bounds_optimization_basis <- function(object) {
  if (inherits(object, "cdmc_dr_fit")) {
    return("dr_optimization_local_refit")
  }

  if (!inherits(object, "cdmc_fit")) {
    stop("object must inherit from 'cdmc_fit' or 'cdmc_dr_fit'.", call. = FALSE)
  }

  if (identical(object$objective %||% "staged", "joint")) {
    return("joint_optimization_local_refit")
  }

  "fit_optimization_local_refit"
}

cdmc_sensitivity_bounds_optimization_scale_source <- function(object) {
  source_object <- cdmc_sensitivity_bounds_optimization_source(object)
  optimization_mask <- cdmc_sensitivity_bounds_optimization_mask(source_object)
  if (inherits(source_object, "cdmc_dr_fit")) {
    return(list(
      response = as.numeric(source_object$y_matrix[optimization_mask]),
      residuals = as.numeric(source_object$effect$tau[optimization_mask]),
      label = cdmc_sensitivity_bounds_optimization_basis(source_object)
    ))
  }

  optimization_fitted <- source_object$baseline$joint_full_fitted %||% source_object$baseline$baseline_hat

  list(
    response = as.numeric(source_object$y_matrix[optimization_mask]),
    residuals = as.numeric(source_object$y_matrix[optimization_mask] - optimization_fitted[optimization_mask]),
    label = cdmc_sensitivity_bounds_optimization_basis(source_object)
  )
}

cdmc_sensitivity_bounds_resolve_refit_step <- function(refit_step, scale_info) {
  if (is.null(refit_step)) {
    return(max(1e-4, 1e-3 * scale_info$value))
  }

  if (!is.numeric(refit_step) || length(refit_step) != 1L || !is.finite(refit_step) || refit_step <= 0) {
    stop("refit_step must be NULL or a single positive numeric value.", call. = FALSE)
  }

  as.numeric(refit_step)
}

cdmc_sensitivity_bounds_collect_aligned_statistics <- function(
  object,
  statistics,
  statistic_names,
  contrast_history = NULL,
  reference_history = NULL,
  prediction_dose = NULL,
  prediction_history = NULL,
  prediction_type = c("response", "slope")
) {
  prediction_type <- match.arg(prediction_type)
  collected <- cdmc_collect_bootstrap_statistics(
    object,
    statistics = statistics,
    contrast_history = contrast_history,
    reference_history = reference_history,
    prediction_dose = prediction_dose,
    prediction_history = prediction_history,
    prediction_type = prediction_type
  )

  aligned <- rep(NA_real_, length(statistic_names))
  names(aligned) <- statistic_names
  common_names <- intersect(names(collected), statistic_names)
  if (length(common_names) > 0L) {
    aligned[common_names] <- collected[common_names]
  }

  aligned
}

cdmc_sensitivity_bounds_refit_specs <- function(
  object,
  collected,
  statistics,
  refit_step,
  contrast_history = NULL,
  reference_history = NULL,
  prediction_dose = NULL,
  prediction_history = NULL,
  prediction_type = c("response", "slope")
) {
  source_object <- cdmc_sensitivity_bounds_optimization_source(object)
  prediction_type <- match.arg(prediction_type)

  if (any(!is.finite(collected))) {
    stop(
      "Optimization-layer sensitivity bounds require finite reference statistics for all requested summaries.",
      call. = FALSE
    )
  }

  full_data <- source_object$data
  observed_rows <- rep(TRUE, nrow(full_data))
  if (".cdmc_observed" %in% names(full_data)) {
    observed_rows <- !is.na(full_data$.cdmc_observed) & full_data$.cdmc_observed
  }

  original_columns <- cdmc_sensitivity_original_columns(object)
  fit_data <- full_data[observed_rows, original_columns, drop = FALSE]
  row_lookup <- integer(length(observed_rows))
  row_lookup[observed_rows] <- seq_len(sum(observed_rows))

  perturbation_mask <- cdmc_sensitivity_bounds_optimization_mask(source_object)
  perturbation_rows <- row_lookup[observed_rows & cdmc_flatten_matrix(perturbation_mask)]
  perturbation_rows <- perturbation_rows[perturbation_rows > 0L]
  if (length(perturbation_rows) < 1L) {
    stop("No optimization cells were available for optimization-layer sensitivity bounds.", call. = FALSE)
  }

  perturbation_labels <- paste0("cell_", seq_along(perturbation_rows))
  if (all(c(source_object$unit, source_object$time) %in% names(fit_data))) {
    perturbation_labels <- make.names(
      paste0(
        as.character(fit_data[[source_object$unit]][perturbation_rows]),
        "@",
        as.character(fit_data[[source_object$time]][perturbation_rows])
      ),
      unique = TRUE
    )
  }

  fit_spec <- tryCatch(
    cdmc_build_bootstrap_fit_spec(object, rerun_tuning = FALSE),
    error = function(error) error
  )
  if (inherits(fit_spec, "error")) {
    stop(
      paste(
        "Optimization-layer sensitivity bounds require a source fit that can be rebuilt with fixed tuning.",
        conditionMessage(fit_spec)
      ),
      call. = FALSE
    )
  }

  weight_scenario <- cdmc_resolve_sensitivity_weight_specs(
    object = source_object,
    source_object = source_object,
    weight_specs = NULL,
    include_unweighted = FALSE
  )$current
  outcome_name <- source_object$outcome
  base_outcome <- fit_data[[outcome_name]]
  jacobian <- matrix(NA_real_, nrow = length(collected), ncol = length(perturbation_rows))
  rownames(jacobian) <- names(collected)
  colnames(jacobian) <- perturbation_labels

  for (index in seq_along(perturbation_rows)) {
    perturbed_data <- fit_data
    perturbed_data[[outcome_name]][perturbation_rows[[index]]] <- base_outcome[[perturbation_rows[[index]]]] + refit_step

    refit_source <- tryCatch(
      cdmc_sensitivity_refit_source(
        source_object = source_object,
        data = perturbed_data,
        fit_spec = fit_spec,
        weight_scenario = weight_scenario,
        washout = source_object$washout,
        zero_tolerance = source_object$zero_tolerance
      ),
      error = function(error) error
    )
    if (inherits(refit_source, "error")) {
      stop(
        sprintf(
          "Optimization-layer sensitivity refit failed after perturbing %s: %s.",
          perturbation_labels[[index]],
          conditionMessage(refit_source)
        ),
        call. = FALSE
      )
    }

    target_object <- tryCatch(
      cdmc_bootstrap_evaluate_target(object, refit_source),
      error = function(error) error
    )
    if (inherits(target_object, "error")) {
      stop(
        sprintf(
          "Optimization-layer sensitivity target evaluation failed after perturbing %s: %s.",
          perturbation_labels[[index]],
          conditionMessage(target_object)
        ),
        call. = FALSE
      )
    }

    perturbed_statistics <- cdmc_sensitivity_bounds_collect_aligned_statistics(
      object = target_object,
      statistics = statistics,
      statistic_names = names(collected),
      contrast_history = contrast_history,
      reference_history = reference_history,
      prediction_dose = prediction_dose,
      prediction_history = prediction_history,
      prediction_type = prediction_type
    )
    if (any(!is.finite(perturbed_statistics))) {
      stop(
        sprintf(
          "Optimization-layer sensitivity produced nonfinite statistics after perturbing %s.",
          perturbation_labels[[index]]
        ),
        call. = FALSE
      )
    }

    jacobian[, index] <- (perturbed_statistics - collected) / refit_step
  }

  specs <- vector("list", length(collected))
  names(specs) <- names(collected)
  basis <- cdmc_sensitivity_bounds_optimization_basis(source_object)
  unit_ids <- if (source_object$unit %in% names(fit_data)) {
    as.character(fit_data[[source_object$unit]][perturbation_rows])
  } else {
    NULL
  }
  time_ids <- if (source_object$time %in% names(fit_data)) {
    as.character(fit_data[[source_object$time]][perturbation_rows])
  } else {
    NULL
  }
  for (statistic_name in names(collected)) {
    specs[[statistic_name]] <- list(
      loading = as.numeric(jacobian[statistic_name, , drop = TRUE]),
      basis = basis,
      unit_ids = unit_ids,
      time_ids = time_ids
    )
  }

  structure(
    specs,
    bound_type = "local_refit",
    refit_step = refit_step,
    n_perturbation_cells = length(perturbation_rows)
  )
}

cdmc_sensitivity_bounds_fit_control_residuals <- function(object) {
  if (!inherits(object, "cdmc_fit")) {
    stop("object must inherit from 'cdmc_fit'.", call. = FALSE)
  }

  candidate_masks <- list(
    object$eligible_mask,
    cdmc_zero_dose_mask(object$dose_matrix, zero_tolerance = object$zero_tolerance) & object$observed_mask
  )

  for (mask in candidate_masks) {
    values <- as.numeric(object$effect$tau[mask])
    values <- values[is.finite(values)]
    if (length(values) > 0L) {
      return(values)
    }
  }

  numeric(0)
}

cdmc_sensitivity_bounds_dr_control_residuals <- function(object) {
  if (!inherits(object, "cdmc_dr_fit")) {
    stop("object must inherit from 'cdmc_dr_fit'.", call. = FALSE)
  }

  eligible_mask <- cdmc_build_eligible_mask(
    dose_matrix = object$dose_matrix,
    zero_tolerance = object$zero_tolerance,
    washout = object$washout
  )
  candidate_masks <- list(
    eligible_mask & object$observed_mask,
    cdmc_zero_dose_mask(object$dose_matrix, zero_tolerance = object$zero_tolerance) & object$observed_mask
  )

  for (mask in candidate_masks) {
    values <- as.numeric(object$effect$tau[mask])
    values <- values[is.finite(values)]
    if (length(values) > 0L) {
      return(values)
    }
  }

  numeric(0)
}

cdmc_sensitivity_bounds_history_design <- function(design_info, history) {
  built <- cdmc_build_response_design(
    history = cdmc_prepare_dynamic_design_history(history),
    model = design_info$model,
    df = design_info$df %||% 4L,
    basis_spec = design_info$basis_spec
  )

  rows <- matrix(0, nrow = nrow(built$design), ncol = length(design_info$column_names))
  colnames(rows) <- design_info$column_names
  common_columns <- intersect(colnames(built$design), design_info$column_names)
  if (length(common_columns) > 0L) {
    rows[, common_columns] <- as.matrix(built$design[, common_columns, drop = FALSE])
  }

  rows
}

cdmc_sensitivity_bounds_fit_map <- function(object, perturbation_layer = c("stage2", "baseline")) {
  if (!inherits(object, "cdmc_fit")) {
    stop("object must inherit from 'cdmc_fit'.", call. = FALSE)
  }

  perturbation_layer <- cdmc_sensitivity_bounds_resolve_perturbation_layer(perturbation_layer)
  objective <- object$objective %||% "staged"

  design_info <- object$effect$design_info
  if (is.null(design_info) || is.null(design_info$design) || nrow(design_info$design) < 1L) {
    return(NULL)
  }
  sample_indices <- which(object$effect$sample_mask, arr.ind = TRUE)
  panel_ids <- cdmc_sensitivity_bounds_panel_ids(
    unit_levels = object$unit_levels,
    time_levels = object$time_levels,
    sample_indices = sample_indices
  )

  list(
    response = as.numeric(object$effect$tau[object$effect$sample_mask]),
    residuals = as.numeric(object$effect$residuals %||% numeric(0)),
    design = as.matrix(design_info$design),
    weights = object$effect$weights,
    operator = cdmc_sensitivity_bounds_operator(design_info$design, weights = object$effect$weights),
    coefficient_names = design_info$column_names,
    design_info = design_info,
    unit_ids = panel_ids$unit_ids,
    time_ids = panel_ids$time_ids,
    basis = if (identical(perturbation_layer, "baseline")) {
      if (identical(objective, "joint")) "joint_baseline_counterfactual" else "fit_baseline_counterfactual"
    } else {
      if (identical(objective, "joint")) "joint_effect_regression" else "fit_effect_regression"
    }
  )
}

cdmc_sensitivity_bounds_dr_map <- function(object, perturbation_layer = c("stage2", "baseline")) {
  if (!inherits(object, "cdmc_dr_fit")) {
    stop("object must inherit from 'cdmc_dr_fit'.", call. = FALSE)
  }

  perturbation_layer <- cdmc_sensitivity_bounds_resolve_perturbation_layer(perturbation_layer)

  history <- cdmc_dr_history_table(object)
  design_columns <- cdmc_dr_design_columns(object)
  if (nrow(history) < 1L || length(design_columns) < 1L) {
    return(NULL)
  }
  sample_indices <- which(object$effect$sample_mask, arr.ind = TRUE)
  panel_ids <- cdmc_sensitivity_bounds_panel_ids(
    unit_levels = object$unit_levels,
    time_levels = object$time_levels,
    sample_indices = sample_indices
  )

  design <- as.matrix(history[, design_columns, drop = FALSE])
  sample_weights <- if (is.null(object$weight_matrix)) {
    rep(1, nrow(design))
  } else {
    as.numeric(object$weight_matrix[object$effect$sample_mask])
  }

  list(
    response = as.numeric(object$effect$tau_dr[object$effect$sample_mask]),
    residuals = as.numeric(object$effect$residuals %||% numeric(0)),
    design = design,
    weights = NULL,
    operator = cdmc_sensitivity_bounds_operator(design),
    coefficient_names = design_columns,
    history = history,
    sample_weights = sample_weights,
    unit_ids = panel_ids$unit_ids,
    time_ids = panel_ids$time_ids,
    basis = if (identical(perturbation_layer, "baseline")) {
      "dr_baseline_counterfactual"
    } else {
      "dr_pseudo_outcome_regression"
    }
  )
}

cdmc_sensitivity_bounds_dr_operator <- function(map, perturbation_layer = c("stage2", "baseline")) {
  perturbation_layer <- cdmc_sensitivity_bounds_resolve_perturbation_layer(perturbation_layer)

  if (!identical(perturbation_layer, "baseline")) {
    return(map$operator)
  }

  sweep(map$operator, 2, map$sample_weights, FUN = "*")
}

cdmc_sensitivity_bounds_dose_response_weights <- function(object) {
  if (!inherits(object, "cdmc_dose_response")) {
    stop("object must inherit from 'cdmc_dose_response'.", call. = FALSE)
  }

  prepared <- cdmc_prepare_effect_sample(
    object = object$fit_object,
    lag_order = object$lag_order,
    include_zero_dose = object$include_zero_dose
  )

  cdmc_resolve_effect_weights(
    weights = object$weights,
    object = object$fit_object,
    sample_mask = prepared$sample_mask
  )
}

cdmc_sensitivity_bounds_dose_response_map <- function(object, perturbation_layer = c("stage2", "baseline")) {
  if (!inherits(object, "cdmc_dose_response")) {
    stop("object must inherit from 'cdmc_dose_response'.", call. = FALSE)
  }

  perturbation_layer <- cdmc_sensitivity_bounds_resolve_perturbation_layer(perturbation_layer)

  design_info <- object$design_info
  if (is.null(design_info) || is.null(design_info$design) || nrow(design_info$design) < 1L) {
    return(NULL)
  }

  resolved_weights <- cdmc_sensitivity_bounds_dose_response_weights(object)
  sample_indices <- cbind(object$history$row, object$history$col)
  panel_ids <- cdmc_sensitivity_bounds_panel_ids(
    unit_levels = object$fit_object$unit_levels,
    time_levels = object$fit_object$time_levels,
    sample_indices = sample_indices
  )

  list(
    response = as.numeric(object$history$tau),
    residuals = as.numeric(object$residuals %||% numeric(0)),
    design = as.matrix(design_info$design),
    weights = resolved_weights,
    operator = cdmc_sensitivity_bounds_operator(design_info$design, weights = resolved_weights),
    coefficient_names = design_info$column_names,
    design_info = design_info,
    unit_ids = panel_ids$unit_ids,
    time_ids = panel_ids$time_ids,
    basis = if (identical(perturbation_layer, "baseline")) {
      "dose_response_baseline_counterfactual"
    } else {
      "dose_response_regression"
    }
  )
}

cdmc_sensitivity_bounds_primary_scale_source <- function(object, perturbation_layer = c("stage2", "baseline")) {
  perturbation_layer <- cdmc_sensitivity_bounds_resolve_perturbation_layer(perturbation_layer)

  if (identical(perturbation_layer, "optimization")) {
    return(cdmc_sensitivity_bounds_optimization_scale_source(object))
  }

  if (inherits(object, "cdmc_dynamic_estimand")) {
    return(cdmc_sensitivity_bounds_primary_scale_source(object$source_object, perturbation_layer = perturbation_layer))
  }

  if (inherits(object, "cdmc_dose_response")) {
    if (identical(perturbation_layer, "baseline") && inherits(object$fit_object, "cdmc_fit")) {
      return(cdmc_sensitivity_bounds_primary_scale_source(object$fit_object, perturbation_layer = perturbation_layer))
    }

    map <- cdmc_sensitivity_bounds_dose_response_map(object, perturbation_layer = perturbation_layer)
    if (is.null(map)) {
      stop("No stage-2 dose-response sample is available for sensitivity scaling.", call. = FALSE)
    }
    return(list(
      response = map$response,
      residuals = map$residuals,
      label = map$basis
    ))
  }

  if (inherits(object, "cdmc_dr_fit")) {
    if (identical(perturbation_layer, "baseline")) {
      control_residuals <- cdmc_sensitivity_bounds_dr_control_residuals(object)
      if (length(control_residuals) > 0L) {
        return(list(
          response = control_residuals,
          residuals = control_residuals,
          label = "dr_baseline_control_residuals"
        ))
      }
    }

    map <- cdmc_sensitivity_bounds_dr_map(object, perturbation_layer = perturbation_layer)
    if (is.null(map)) {
      stop("No DR pseudo-outcome sample is available for sensitivity scaling.", call. = FALSE)
    }
    return(list(
      response = map$response,
      residuals = map$residuals,
      label = map$basis
    ))
  }

  if (inherits(object, "cdmc_fit")) {
    if (identical(perturbation_layer, "baseline")) {
      control_residuals <- cdmc_sensitivity_bounds_fit_control_residuals(object)
      if (length(control_residuals) > 0L) {
        return(list(
          response = control_residuals,
          residuals = control_residuals,
          label = if (identical(object$objective %||% "staged", "joint")) {
            "joint_baseline_control_residuals"
          } else {
            "fit_baseline_control_residuals"
          }
        ))
      }
    }

    map <- cdmc_sensitivity_bounds_fit_map(object, perturbation_layer = perturbation_layer)
    if (!is.null(map)) {
      return(list(
        response = map$response,
        residuals = map$residuals,
        label = map$basis
      ))
    }

    active_mask <- cdmc_active_dose_mask(object$dose_matrix, zero_tolerance = object$zero_tolerance)
    response <- as.numeric(object$effect$tau[active_mask])
    if (length(response) < 1L) {
      stop("No active-dose treatment-effect sample is available for sensitivity scaling.", call. = FALSE)
    }

    return(list(
      response = response,
      residuals = numeric(0),
      label = "active_treatment_mean"
    ))
  }

  stop(
    "object must inherit from 'cdmc_fit', 'cdmc_dr_fit', 'cdmc_dose_response', or 'cdmc_dynamic_estimand'.",
    call. = FALSE
  )
}

cdmc_sensitivity_bounds_resolve_scale <- function(
  object,
  perturbation_layer = c("stage2", "baseline"),
  scale = c("residual_sd", "response_sd", "manual"),
  scale_value = NULL
) {
  perturbation_layer <- cdmc_sensitivity_bounds_resolve_perturbation_layer(perturbation_layer)
  scale <- match.arg(scale)

  if (identical(scale, "manual")) {
    if (!is.numeric(scale_value) || length(scale_value) != 1L || !is.finite(scale_value) || scale_value <= 0) {
      stop("scale_value must be a single positive numeric value when scale = 'manual'.", call. = FALSE)
    }

    return(list(
      method = scale,
      value = as.numeric(scale_value),
      source = "manual"
    ))
  }

  source <- cdmc_sensitivity_bounds_primary_scale_source(object, perturbation_layer = perturbation_layer)
  values <- if (identical(scale, "residual_sd")) source$residuals else source$response
  values <- values[is.finite(values)]

  if (length(values) < 2L) {
    stop(
      sprintf("Scale '%s' requires at least two finite values in the underlying sensitivity sample.", scale),
      call. = FALSE
    )
  }

  resolved <- stats::sd(values)
  if (!is.finite(resolved) || resolved <= sqrt(.Machine$double.eps)) {
    stop(
      sprintf("Scale '%s' resolved to a nonpositive or degenerate value.", scale),
      call. = FALSE
    )
  }

  list(
    method = scale,
    value = resolved,
    source = source$label
  )
}

cdmc_sensitivity_bounds_resolve_null <- function(null, statistic_names) {
  if (length(statistic_names) < 1L) {
    return(numeric(0))
  }

  if (!is.numeric(null) || any(!is.finite(null))) {
    stop("null must be numeric and finite.", call. = FALSE)
  }

  if (length(null) == 1L && is.null(names(null))) {
    out <- rep(as.numeric(null), length(statistic_names))
    names(out) <- statistic_names
    return(out)
  }

  if (is.null(names(null))) {
    if (length(null) != length(statistic_names)) {
      stop(
        "Unnamed null inputs must have length 1 or match the number of reported statistics.",
        call. = FALSE
      )
    }

    out <- as.numeric(null)
    names(out) <- statistic_names
    return(out)
  }

  unknown_names <- setdiff(names(null), statistic_names)
  if (length(unknown_names) > 0L) {
    stop(
      sprintf("null contains unknown statistic names: %s.", paste(unknown_names, collapse = ", ")),
      call. = FALSE
    )
  }

  out <- rep(0, length(statistic_names))
  names(out) <- statistic_names
  out[names(null)] <- as.numeric(null)
  out
}

cdmc_sensitivity_bounds_breakdown <- function(estimate, null_value, multiplier) {
  if (!is.finite(estimate) || !is.finite(null_value) || !is.finite(multiplier)) {
    return(NA_real_)
  }

  distance <- abs(estimate - null_value)
  if (multiplier <= sqrt(.Machine$double.eps)) {
    return(if (distance <= sqrt(.Machine$double.eps)) 0 else Inf)
  }

  distance / multiplier
}

cdmc_sensitivity_bounds_resolve_contamination_fraction_grid <- function(contamination_fraction_grid) {
  resolved <- sort(unique(as.numeric(contamination_fraction_grid)))
  if (length(resolved) < 1L || any(!is.finite(resolved)) || any(resolved < 0) || any(resolved > 1)) {
    stop("contamination_fraction_grid must contain one or more numeric values in [0, 1].", call. = FALSE)
  }

  resolved
}

cdmc_sensitivity_bounds_resolve_perturbation_constraint <- function(perturbation_constraint) {
  match.arg(perturbation_constraint, c("cellwise", "energy"))
}

cdmc_sensitivity_bounds_resolve_perturbation_scope <- function(perturbation_scope) {
  match.arg(perturbation_scope, c("cell", "unit", "time"))
}

cdmc_sensitivity_bounds_panel_ids <- function(unit_levels, time_levels, sample_indices) {
  if (is.null(sample_indices) || nrow(sample_indices) < 1L) {
    return(list(unit_ids = character(0), time_ids = character(0)))
  }

  list(
    unit_ids = as.character(unit_levels[sample_indices[, 1L]]),
    time_ids = as.character(time_levels[sample_indices[, 2L]])
  )
}

cdmc_sensitivity_bounds_loading_groups <- function(
  loading,
  perturbation_scope = c("cell", "unit", "time"),
  unit_ids = NULL,
  time_ids = NULL
) {
  perturbation_scope <- cdmc_sensitivity_bounds_resolve_perturbation_scope(perturbation_scope)
  loading <- as.numeric(loading)

  if (identical(perturbation_scope, "cell")) {
    return(loading)
  }

  group_ids <- if (identical(perturbation_scope, "unit")) unit_ids else time_ids
  if (is.null(group_ids)) {
    stop(
      sprintf("%s-scope sensitivity requires retained %s identifiers for the perturbed sample.", perturbation_scope, perturbation_scope),
      call. = FALSE
    )
  }
  if (length(group_ids) != length(loading)) {
    stop(
      sprintf("%s identifiers must have the same length as the loading vector.", perturbation_scope),
      call. = FALSE
    )
  }

  aggregated <- stats::aggregate(
    x = loading,
    by = list(group = as.character(group_ids)),
    FUN = sum,
    na.rm = FALSE
  )
  stats::setNames(as.numeric(aggregated$x), aggregated$group)
}

cdmc_sensitivity_bounds_contamination_cells <- function(contamination_fraction, n_cells) {
  n_cells <- as.integer(n_cells)
  if (!is.finite(contamination_fraction) || !is.finite(n_cells) || n_cells < 0L) {
    return(NA_integer_)
  }
  if (n_cells == 0L || contamination_fraction <= 0) {
    return(0L)
  }

  min(n_cells, max(1L, ceiling(contamination_fraction * n_cells)))
}

cdmc_sensitivity_bounds_effective_multiplier <- function(
  loading,
  contamination_fraction = 1,
  perturbation_constraint = c("cellwise", "energy"),
  perturbation_scope = c("cell", "unit", "time"),
  unit_ids = NULL,
  time_ids = NULL
) {
  loading <- as.numeric(loading)
  perturbation_constraint <- cdmc_sensitivity_bounds_resolve_perturbation_constraint(perturbation_constraint)
  perturbation_scope <- cdmc_sensitivity_bounds_resolve_perturbation_scope(perturbation_scope)
  if (length(loading) < 1L || any(!is.finite(loading))) {
    return(NA_real_)
  }

  grouped_loadings <- cdmc_sensitivity_bounds_loading_groups(
    loading = loading,
    perturbation_scope = perturbation_scope,
    unit_ids = unit_ids,
    time_ids = time_ids
  )
  active_cells <- cdmc_sensitivity_bounds_contamination_cells(contamination_fraction, length(grouped_loadings))
  if (!is.finite(active_cells) || active_cells <= 0L) {
    return(0)
  }

  top_loadings <- sort(abs(grouped_loadings), decreasing = TRUE)[seq_len(active_cells)]
  if (identical(perturbation_constraint, "energy")) {
    return(sqrt(sum(top_loadings ^ 2)))
  }

  sum(top_loadings)
}

cdmc_sensitivity_bounds_fit_specs <- function(object, collected, perturbation_layer = c("stage2", "baseline")) {
  if (!inherits(object, "cdmc_fit")) {
    stop("object must inherit from 'cdmc_fit'.", call. = FALSE)
  }

  perturbation_layer <- cdmc_sensitivity_bounds_resolve_perturbation_layer(perturbation_layer)

  specs <- list()
  map <- cdmc_sensitivity_bounds_fit_map(object, perturbation_layer = perturbation_layer)

  if (!is.null(map) && length(map$coefficient_names) > 0L) {
    for (coefficient_name in map$coefficient_names) {
      statistic_name <- paste0("coef_", coefficient_name)
      if (statistic_name %in% names(collected)) {
        specs[[statistic_name]] <- list(
          loading = as.numeric(map$operator[coefficient_name, , drop = TRUE]),
          basis = map$basis,
          unit_ids = map$unit_ids,
          time_ids = map$time_ids
        )
      }
    }
  }

  if ("mean_tau_active_dose" %in% names(collected)) {
    active_mask <- cdmc_active_dose_mask(object$dose_matrix, zero_tolerance = object$zero_tolerance)
    active_tau <- as.numeric(object$effect$tau[active_mask])
    if (length(active_tau) > 0L) {
      loading <- if (isTRUE(object$weight_supplied)) {
        cdmc_sensitivity_bounds_normalized_loading(object$weight_matrix[active_mask])
      } else {
        rep(1 / length(active_tau), length(active_tau))
      }
      specs[["mean_tau_active_dose"]] <- list(
        loading = loading,
        basis = if (identical(perturbation_layer, "baseline")) {
          "active_treatment_baseline_counterfactual"
        } else {
          "active_treatment_mean"
        },
        unit_ids = as.character(object$unit_levels[which(active_mask, arr.ind = TRUE)[, 1L]]),
        time_ids = as.character(object$time_levels[which(active_mask, arr.ind = TRUE)[, 2L]])
      )
    }
  }

  specs
}

cdmc_sensitivity_bounds_dr_specs <- function(
  object,
  collected,
  contrast_history = NULL,
  reference_history = NULL,
  perturbation_layer = c("stage2", "baseline")
) {
  if (!inherits(object, "cdmc_dr_fit")) {
    stop("object must inherit from 'cdmc_dr_fit'.", call. = FALSE)
  }

  perturbation_layer <- cdmc_sensitivity_bounds_resolve_perturbation_layer(perturbation_layer)

  specs <- list()
  map <- cdmc_sensitivity_bounds_dr_map(object, perturbation_layer = perturbation_layer)
  if (is.null(map)) {
    return(specs)
  }
  operator <- cdmc_sensitivity_bounds_dr_operator(map, perturbation_layer = perturbation_layer)

  for (coefficient_name in map$coefficient_names) {
    statistic_name <- paste0("coef_", coefficient_name)
    if (statistic_name %in% names(collected)) {
      specs[[statistic_name]] <- list(
        loading = as.numeric(operator[coefficient_name, , drop = TRUE]),
        basis = map$basis,
        unit_ids = map$unit_ids,
        time_ids = map$time_ids
      )
    }
  }

  if ("mean_tau_dr" %in% names(collected) && length(map$response) > 0L) {
    specs[["mean_tau_dr"]] <- list(
      loading = if (identical(perturbation_layer, "baseline")) {
        map$sample_weights / length(map$sample_weights)
      } else {
        rep(1 / length(map$response), length(map$response))
      },
      basis = if (identical(perturbation_layer, "baseline")) {
        "dr_baseline_counterfactual_mean"
      } else {
        "dr_pseudo_outcome_mean"
      },
      unit_ids = map$unit_ids,
      time_ids = map$time_ids
    )
  }

  if ("mean_tau_dr_linear" %in% names(collected)) {
    coefficient_row <- matrix(colMeans(map$design), nrow = 1L)
    colnames(coefficient_row) <- map$coefficient_names
    loading <- cdmc_sensitivity_bounds_align_coefficient_rows(coefficient_row, map$coefficient_names) %*% operator
    specs[["mean_tau_dr_linear"]] <- list(
      loading = as.numeric(loading[1L, , drop = TRUE]),
      basis = map$basis,
      unit_ids = map$unit_ids,
      time_ids = map$time_ids
    )
  }

  lag_means <- colMeans(map$design)
  for (coefficient_name in map$coefficient_names) {
    statistic_name <- paste0("mean_tau_dr_", coefficient_name)
    if (statistic_name %in% names(collected)) {
      coefficient_row <- matrix(0, nrow = 1L, ncol = length(map$coefficient_names))
      colnames(coefficient_row) <- map$coefficient_names
      coefficient_row[1L, coefficient_name] <- lag_means[[coefficient_name]]
      loading <- coefficient_row %*% operator
      specs[[statistic_name]] <- list(
        loading = as.numeric(loading[1L, , drop = TRUE]),
        basis = map$basis,
        unit_ids = map$unit_ids,
        time_ids = map$time_ids
      )
    }
  }

  if ("dr_path_contrast" %in% names(collected)) {
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
    coefficient_row <- matrix(contrast_values - reference_values, nrow = 1L)
    colnames(coefficient_row) <- map$coefficient_names
    loading <- coefficient_row %*% operator
    specs[["dr_path_contrast"]] <- list(
      loading = as.numeric(loading[1L, , drop = TRUE]),
      basis = map$basis,
      unit_ids = map$unit_ids,
      time_ids = map$time_ids
    )
  }

  specs
}

cdmc_sensitivity_bounds_dose_response_prediction_rows <- function(
  object,
  prediction_dose = NULL,
  prediction_history = NULL,
  prediction_type = c("response", "slope"),
  eps = 1e-4
) {
  prediction_type <- match.arg(prediction_type)
  history <- cdmc_prepare_prediction_history(object, dose = prediction_dose, history = prediction_history)

  if (identical(prediction_type, "response")) {
    rows <- cdmc_sensitivity_bounds_history_design(object$design_info, history)
  } else if (identical(object$model, "linear")) {
    rows <- matrix(0, nrow = nrow(history), ncol = length(object$design_info$column_names))
    colnames(rows) <- object$design_info$column_names
    if ("dose_lag0" %in% colnames(rows)) {
      rows[, "dose_lag0"] <- 1
    }
  } else {
    history_up <- history
    history_down <- history
    history_up$dose_lag0 <- history_up$dose_lag0 + eps
    history_down$dose_lag0 <- history_down$dose_lag0 - eps
    rows_up <- cdmc_sensitivity_bounds_history_design(object$design_info, history_up)
    rows_down <- cdmc_sensitivity_bounds_history_design(object$design_info, history_down)
    rows <- (rows_up - rows_down) / (2 * eps)
  }

  rownames(rows) <- paste0("dose_response_", prediction_type, "_", seq_len(nrow(rows)))
  rows
}

cdmc_sensitivity_bounds_dose_response_specs <- function(
  object,
  collected,
  prediction_dose = NULL,
  prediction_history = NULL,
  prediction_type = c("response", "slope"),
  perturbation_layer = c("stage2", "baseline")
) {
  if (!inherits(object, "cdmc_dose_response")) {
    stop("object must inherit from 'cdmc_dose_response'.", call. = FALSE)
  }

  perturbation_layer <- cdmc_sensitivity_bounds_resolve_perturbation_layer(perturbation_layer)
  prediction_type <- match.arg(prediction_type)
  specs <- list()
  map <- cdmc_sensitivity_bounds_dose_response_map(object, perturbation_layer = perturbation_layer)
  if (is.null(map)) {
    return(specs)
  }

  for (coefficient_name in map$coefficient_names) {
    statistic_name <- paste0("coef_", coefficient_name)
    if (statistic_name %in% names(collected)) {
      specs[[statistic_name]] <- list(
        loading = as.numeric(map$operator[coefficient_name, , drop = TRUE]),
        basis = map$basis,
        unit_ids = map$unit_ids,
        time_ids = map$time_ids
      )
    }
  }

  prediction_rows <- names(collected)[grepl("^dose_response_(response|slope)_", names(collected))]
  if (length(prediction_rows) > 0L) {
    coefficient_rows <- cdmc_sensitivity_bounds_dose_response_prediction_rows(
      object = object,
      prediction_dose = prediction_dose,
      prediction_history = prediction_history,
      prediction_type = prediction_type
    )
    coefficient_rows <- coefficient_rows[prediction_rows, , drop = FALSE]
    loading_matrix <- coefficient_rows %*% map$operator
    for (statistic_name in rownames(loading_matrix)) {
      specs[[statistic_name]] <- list(
        loading = as.numeric(loading_matrix[statistic_name, , drop = TRUE]),
        basis = map$basis,
        unit_ids = map$unit_ids,
        time_ids = map$time_ids
      )
    }
  }

  specs
}

cdmc_sensitivity_bounds_fit_dynamic_rows <- function(object, history, reference_history = NULL, type, eps) {
  if (!inherits(object, "cdmc_fit")) {
    stop("object must inherit from 'cdmc_fit'.", call. = FALSE)
  }

  design_info <- object$effect$design_info
  if (is.null(design_info)) {
    stop("Dynamic sensitivity bounds require a retained fit effect design.", call. = FALSE)
  }

  if (identical(type, "response")) {
    return(cdmc_sensitivity_bounds_history_design(design_info, history))
  }

  if (identical(type, "contrast")) {
    response_rows <- cdmc_sensitivity_bounds_history_design(design_info, history)
    reference_rows <- cdmc_sensitivity_bounds_history_design(design_info, reference_history)
    return(response_rows - reference_rows)
  }

  if (identical(object$effect_model, "linear")) {
    rows <- matrix(0, nrow = nrow(history), ncol = length(design_info$column_names))
    colnames(rows) <- design_info$column_names
    if ("dose_lag0" %in% colnames(rows)) {
      rows[, "dose_lag0"] <- 1
    }
    return(rows)
  }

  history_up <- history
  history_down <- history
  history_up$dose_lag0 <- history_up$dose_lag0 + eps
  history_down$dose_lag0 <- history_down$dose_lag0 - eps
  rows_up <- cdmc_sensitivity_bounds_history_design(design_info, history_up)
  rows_down <- cdmc_sensitivity_bounds_history_design(design_info, history_down)
  (rows_up - rows_down) / (2 * eps)
}

cdmc_sensitivity_bounds_dr_dynamic_rows <- function(object, history, reference_history = NULL, type) {
  if (!inherits(object, "cdmc_dr_fit")) {
    stop("object must inherit from 'cdmc_dr_fit'.", call. = FALSE)
  }

  coefficient_names <- cdmc_dr_design_columns(object)
  response_rows <- as.matrix(history[, coefficient_names, drop = FALSE])

  if (identical(type, "response")) {
    return(response_rows)
  }

  if (identical(type, "contrast")) {
    reference_rows <- as.matrix(reference_history[, coefficient_names, drop = FALSE])
    return(response_rows - reference_rows)
  }

  rows <- matrix(0, nrow = nrow(history), ncol = length(coefficient_names))
  colnames(rows) <- coefficient_names
  if ("dose_lag0" %in% coefficient_names) {
    rows[, "dose_lag0"] <- 1
  }
  rows
}

cdmc_sensitivity_bounds_dynamic_estimate_rows <- function(object, coefficient_rows) {
  prefix <- cdmc_dynamic_stat_prefix(object$type)
  path_names <- paste0(prefix, "_", make.names(object$estimate_table$label, unique = TRUE))
  rownames(coefficient_rows) <- path_names

  if (identical(object$aggregate, "individual")) {
    return(coefficient_rows)
  }

  mean_row <- matrix(
    as.numeric((object$path_weights / sum(object$path_weights)) %*% coefficient_rows),
    nrow = 1L
  )
  colnames(mean_row) <- colnames(coefficient_rows)
  rownames(mean_row) <- paste0("mean_", prefix)

  if (identical(object$aggregate, "mean")) {
    return(mean_row)
  }

  rbind(coefficient_rows, mean_row)
}

cdmc_sensitivity_bounds_dynamic_specs <- function(object, collected, perturbation_layer = c("stage2", "baseline")) {
  if (!inherits(object, "cdmc_dynamic_estimand")) {
    stop("object must inherit from 'cdmc_dynamic_estimand'.", call. = FALSE)
  }

  perturbation_layer <- cdmc_sensitivity_bounds_resolve_perturbation_layer(perturbation_layer)
  source_object <- object$source_object
  operator <- NULL
  if (inherits(source_object, "cdmc_dose_response")) {
    map <- cdmc_sensitivity_bounds_dose_response_map(source_object, perturbation_layer = perturbation_layer)
    coefficient_rows <- cdmc_sensitivity_bounds_dose_response_prediction_rows(
      object = source_object,
      prediction_history = object$history,
      prediction_type = if (identical(object$type, "slope")) "slope" else "response",
      eps = object$eps
    )
    if (identical(object$type, "contrast")) {
      reference_rows <- cdmc_sensitivity_bounds_dose_response_prediction_rows(
        object = source_object,
        prediction_history = object$reference_history,
        prediction_type = "response",
        eps = object$eps
      )
      coefficient_rows <- coefficient_rows - reference_rows
    }
    operator <- map$operator
  } else if (inherits(source_object, "cdmc_dr_fit")) {
    map <- cdmc_sensitivity_bounds_dr_map(source_object, perturbation_layer = perturbation_layer)
    coefficient_rows <- cdmc_sensitivity_bounds_dr_dynamic_rows(
      object = source_object,
      history = object$history,
      reference_history = object$reference_history,
      type = object$type
    )
    operator <- cdmc_sensitivity_bounds_dr_operator(map, perturbation_layer = perturbation_layer)
  } else if (inherits(source_object, "cdmc_fit")) {
    map <- cdmc_sensitivity_bounds_fit_map(source_object, perturbation_layer = perturbation_layer)
    coefficient_rows <- cdmc_sensitivity_bounds_fit_dynamic_rows(
      object = source_object,
      history = object$history,
      reference_history = object$reference_history,
      type = object$type,
      eps = object$eps
    )
    operator <- map$operator
  } else {
    stop(
      "Dynamic sensitivity bounds require a source object inheriting from 'cdmc_fit', 'cdmc_dr_fit', or 'cdmc_dose_response'.",
      call. = FALSE
    )
  }

  if (is.null(map)) {
    return(list())
  }

  coefficient_rows <- cdmc_sensitivity_bounds_dynamic_estimate_rows(object, coefficient_rows)
  coefficient_rows <- coefficient_rows[names(collected), , drop = FALSE]
  loading_matrix <- cdmc_sensitivity_bounds_align_coefficient_rows(
    coefficient_rows,
    map$coefficient_names
  ) %*% operator

  specs <- vector("list", nrow(loading_matrix))
  names(specs) <- rownames(loading_matrix)
  for (statistic_name in rownames(loading_matrix)) {
    specs[[statistic_name]] <- list(
      loading = as.numeric(loading_matrix[statistic_name, , drop = TRUE]),
      basis = map$basis,
      unit_ids = map$unit_ids,
      time_ids = map$time_ids
    )
  }
  specs
}

cdmc_sensitivity_bounds_build_specs <- function(
  object,
  collected,
  statistics,
  contrast_history = NULL,
  reference_history = NULL,
  prediction_dose = NULL,
  prediction_history = NULL,
  prediction_type = c("response", "slope"),
  perturbation_layer = c("stage2", "baseline", "optimization"),
  refit_step = NULL
) {
  perturbation_layer <- cdmc_sensitivity_bounds_resolve_perturbation_layer(perturbation_layer)
  prediction_type <- match.arg(prediction_type)

  if (identical(perturbation_layer, "optimization")) {
    return(cdmc_sensitivity_bounds_refit_specs(
      object = object,
      collected = collected,
      statistics = statistics,
      refit_step = refit_step,
      contrast_history = contrast_history,
      reference_history = reference_history,
      prediction_dose = prediction_dose,
      prediction_history = prediction_history,
      prediction_type = prediction_type
    ))
  }

  if (inherits(object, "cdmc_dynamic_estimand")) {
    return(cdmc_sensitivity_bounds_dynamic_specs(
      object,
      collected = collected,
      perturbation_layer = perturbation_layer
    ))
  }

  if (inherits(object, "cdmc_dose_response")) {
    return(cdmc_sensitivity_bounds_dose_response_specs(
      object = object,
      collected = collected,
      prediction_dose = prediction_dose,
      prediction_history = prediction_history,
      prediction_type = prediction_type,
      perturbation_layer = perturbation_layer
    ))
  }

  if (inherits(object, "cdmc_dr_fit")) {
    return(cdmc_sensitivity_bounds_dr_specs(
      object = object,
      collected = collected,
      contrast_history = contrast_history,
      reference_history = reference_history,
      perturbation_layer = perturbation_layer
    ))
  }

  if (inherits(object, "cdmc_fit")) {
    return(cdmc_sensitivity_bounds_fit_specs(
      object,
      collected = collected,
      perturbation_layer = perturbation_layer
    ))
  }

  stop(
    "object must inherit from 'cdmc_fit', 'cdmc_dr_fit', 'cdmc_dose_response', or 'cdmc_dynamic_estimand'.",
    call. = FALSE
  )
}

cdmc_sensitivity_bounds_reference_table <- function(
  estimates,
  null_values,
  specs,
  scale_info,
  perturbation_layer,
  perturbation_constraint,
  perturbation_scope
) {
  rows <- vector("list", length(estimates))
  statistic_names <- names(estimates)

  for (index in seq_along(estimates)) {
    statistic_name <- statistic_names[[index]]
    spec <- specs[[statistic_name]]
    loading <- if (!is.null(spec)) spec$loading %||% NULL else NULL
    unit_ids <- if (!is.null(spec)) spec$unit_ids %||% NULL else NULL
    time_ids <- if (!is.null(spec)) spec$time_ids %||% NULL else NULL
    multiplier <- if (!is.null(loading)) {
      cdmc_sensitivity_bounds_effective_multiplier(
        loading,
        contamination_fraction = 1,
        perturbation_constraint = perturbation_constraint,
        perturbation_scope = perturbation_scope,
        unit_ids = unit_ids,
        time_ids = time_ids
      )
    } else {
      NA_real_
    }
    grouped_loadings <- if (!is.null(loading)) {
      cdmc_sensitivity_bounds_loading_groups(
        loading = loading,
        perturbation_scope = perturbation_scope,
        unit_ids = unit_ids,
        time_ids = time_ids
      )
    } else {
      numeric(0)
    }
    breakdown_bias <- cdmc_sensitivity_bounds_breakdown(
      estimate = estimates[[index]],
      null_value = null_values[[statistic_name]],
      multiplier = multiplier
    )

    rows[[index]] <- list(
      statistic = statistic_name,
      estimate = unname(estimates[[index]]),
      null_value = unname(null_values[[statistic_name]]),
      perturbation_layer = perturbation_layer,
      perturbation_constraint = perturbation_constraint,
      perturbation_scope = perturbation_scope,
      basis = spec$basis %||% NA_character_,
      n_loading_cells = if (length(grouped_loadings) > 0L) length(grouped_loadings) else NA_integer_,
      sensitivity_multiplier = multiplier,
      scale_method = scale_info$method,
      scale_value = scale_info$value,
      scale_source = scale_info$source,
      breakdown_bias_bound = breakdown_bias,
      breakdown_gamma = if (is.finite(breakdown_bias)) breakdown_bias / scale_info$value else breakdown_bias
    )
  }

  cdmc_sensitivity_bind_rows(rows)
}

cdmc_sensitivity_bounds_table <- function(
  reference_table,
  specs,
  gamma_grid,
  contamination_fraction_grid,
  perturbation_constraint,
  perturbation_scope
) {
  rows <- vector("list", nrow(reference_table) * length(gamma_grid) * length(contamination_fraction_grid))
  row_index <- 0L

  for (reference_index in seq_len(nrow(reference_table))) {
    reference_row <- reference_table[reference_index, , drop = FALSE]
    spec <- specs[[reference_row$statistic[[1L]]]]
    loading <- if (!is.null(spec)) spec$loading %||% NULL else NULL
    unit_ids <- if (!is.null(spec)) spec$unit_ids %||% NULL else NULL
    time_ids <- if (!is.null(spec)) spec$time_ids %||% NULL else NULL

    grouped_loadings <- if (!is.null(loading)) {
      cdmc_sensitivity_bounds_loading_groups(
        loading = loading,
        perturbation_scope = perturbation_scope,
        unit_ids = unit_ids,
        time_ids = time_ids
      )
    } else {
      numeric(0)
    }
    loading_count <- if (length(grouped_loadings) > 0L) length(grouped_loadings) else NA_integer_

    for (contamination_fraction in contamination_fraction_grid) {
      active_cells <- cdmc_sensitivity_bounds_contamination_cells(contamination_fraction, loading_count)
      multiplier <- if (!is.null(loading)) {
        cdmc_sensitivity_bounds_effective_multiplier(
          loading,
          contamination_fraction = contamination_fraction,
          perturbation_constraint = perturbation_constraint,
          perturbation_scope = perturbation_scope,
          unit_ids = unit_ids,
          time_ids = time_ids
        )
      } else {
        NA_real_
      }

      for (gamma in gamma_grid) {
        row_index <- row_index + 1L
        bias_bound <- gamma * reference_row$scale_value[[1L]]
        radius <- if (is.finite(multiplier)) {
          bias_bound * multiplier
        } else {
          NA_real_
        }

        lower <- if (is.finite(radius) && is.finite(reference_row$estimate[[1L]])) {
          reference_row$estimate[[1L]] - radius
        } else {
          NA_real_
        }
        upper <- if (is.finite(radius) && is.finite(reference_row$estimate[[1L]])) {
          reference_row$estimate[[1L]] + radius
        } else {
          NA_real_
        }

        rows[[row_index]] <- c(as.list(reference_row), list(
          contamination_fraction = contamination_fraction,
          contamination_cells = active_cells,
          effective_multiplier = multiplier,
          gamma = gamma,
          bias_bound = bias_bound,
          radius = radius,
          lower = lower,
          upper = upper,
          excludes_null = if (is.finite(lower) && is.finite(upper) && is.finite(reference_row$null_value[[1L]])) {
            reference_row$null_value[[1L]] < lower || reference_row$null_value[[1L]] > upper
          } else {
            NA
          }
        ))
      }
    }
  }

  cdmc_sensitivity_bind_rows(rows)
}

cdmc_sensitivity_bounds <- function(
  object,
  statistics = NULL,
  gamma_grid = c(0, 0.25, 0.5, 1),
  contamination_fraction_grid = 1,
  perturbation_constraint = c("cellwise", "energy"),
  perturbation_scope = c("cell", "unit", "time"),
  perturbation_layer = c("stage2", "baseline", "optimization"),
  scale = c("residual_sd", "response_sd", "manual"),
  scale_value = NULL,
  null = 0,
  contrast_history = NULL,
  reference_history = NULL,
  prediction_dose = NULL,
  prediction_history = NULL,
  prediction_type = c("response", "slope"),
  refit_step = NULL
) {
  if (!inherits(object, "cdmc_fit") &&
      !inherits(object, "cdmc_dr_fit") &&
      !inherits(object, "cdmc_dose_response") &&
      !inherits(object, "cdmc_dynamic_estimand")) {
    stop(
      "object must inherit from 'cdmc_fit', 'cdmc_dr_fit', 'cdmc_dose_response', or 'cdmc_dynamic_estimand'.",
      call. = FALSE
    )
  }

  prediction_type <- match.arg(prediction_type)
  perturbation_layer <- cdmc_sensitivity_bounds_resolve_perturbation_layer(perturbation_layer)
  perturbation_constraint <- cdmc_sensitivity_bounds_resolve_perturbation_constraint(perturbation_constraint)
  perturbation_scope <- cdmc_sensitivity_bounds_resolve_perturbation_scope(perturbation_scope)

  if (!inherits(object, "cdmc_dr_fit") && (!is.null(contrast_history) || !is.null(reference_history))) {
    stop("contrast_history and reference_history are only supported for cdmc_dr_fit sensitivity bounds.", call. = FALSE)
  }

  if (!inherits(object, "cdmc_dose_response") && (!is.null(prediction_dose) || !is.null(prediction_history))) {
    stop("prediction_dose and prediction_history are only supported for cdmc_dose_response sensitivity bounds.", call. = FALSE)
  }

  gamma_grid <- unique(as.numeric(gamma_grid))
  if (length(gamma_grid) < 1L || any(!is.finite(gamma_grid)) || any(gamma_grid < 0)) {
    stop("gamma_grid must contain one or more nonnegative numeric values.", call. = FALSE)
  }
  contamination_fraction_grid <- cdmc_sensitivity_bounds_resolve_contamination_fraction_grid(contamination_fraction_grid)

  statistics <- cdmc_sensitivity_resolve_statistics(object, statistics = statistics)
  collected <- cdmc_collect_bootstrap_statistics(
    object,
    statistics = statistics,
    contrast_history = contrast_history,
    reference_history = reference_history,
    prediction_dose = prediction_dose,
    prediction_history = prediction_history,
    prediction_type = prediction_type
  )

  if (length(collected) < 1L) {
    stop("No scalar statistics were available for the requested object and statistics.", call. = FALSE)
  }

  scale_info <- cdmc_sensitivity_bounds_resolve_scale(
    object = object,
    perturbation_layer = perturbation_layer,
    scale = scale,
    scale_value = scale_value
  )
  resolved_refit_step <- if (identical(perturbation_layer, "optimization")) {
    cdmc_sensitivity_bounds_resolve_refit_step(refit_step, scale_info)
  } else {
    NA_real_
  }

  specs <- cdmc_sensitivity_bounds_build_specs(
    object = object,
    collected = collected,
    statistics = statistics,
    contrast_history = contrast_history,
    reference_history = reference_history,
    prediction_dose = prediction_dose,
    prediction_history = prediction_history,
    prediction_type = prediction_type,
    perturbation_layer = perturbation_layer,
    refit_step = resolved_refit_step
  )
  bound_type <- attr(specs, "bound_type") %||% "exact_stored_map"
  n_perturbation_cells <- attr(specs, "n_perturbation_cells") %||% NA_integer_
  resolved_refit_step <- attr(specs, "refit_step") %||% resolved_refit_step
  missing_specs <- setdiff(names(collected), names(specs))
  if (length(missing_specs) > 0L) {
    stop(
      sprintf(
        "Formal %s sensitivity bounds are not available for the requested statistics: %s.",
        cdmc_sensitivity_bounds_layer_label(perturbation_layer),
        paste(missing_specs, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  null_values <- cdmc_sensitivity_bounds_resolve_null(null, names(collected))
  reference <- cdmc_sensitivity_bounds_reference_table(
    estimates = collected,
    null_values = null_values,
    specs = specs,
    scale_info = scale_info,
    perturbation_layer = perturbation_layer,
    perturbation_constraint = perturbation_constraint,
    perturbation_scope = perturbation_scope
  )
  bounds <- cdmc_sensitivity_bounds_table(
    reference_table = reference,
    specs = specs,
    gamma_grid = gamma_grid,
    contamination_fraction_grid = contamination_fraction_grid,
    perturbation_constraint = perturbation_constraint,
    perturbation_scope = perturbation_scope
  )

  result <- list(
    call = match.call(),
    target_class = class(object)[1L],
    source_class = class(cdmc_sensitivity_source_object(object))[1L],
    statistics = statistics,
    contamination_fraction_grid = contamination_fraction_grid,
    perturbation_constraint = perturbation_constraint,
    perturbation_scope = perturbation_scope,
    perturbation_layer = perturbation_layer,
    bound_type = bound_type,
    refit_step = resolved_refit_step,
    n_perturbation_cells = n_perturbation_cells,
    scale_method = scale_info$method,
    scale_value = scale_info$value,
    scale_source = scale_info$source,
    reference = reference,
    bounds = bounds,
    base_object = object
  )

  class(result) <- "cdmc_sensitivity_bounds"
  result
}

summary.cdmc_sensitivity_bounds <- function(object, ...) {
  max_rows <- do.call(
    rbind,
    lapply(split(object$bounds, object$bounds$statistic), function(data) {
      ordered <- data[order(data$gamma, data$contamination_fraction), , drop = FALSE]
      ordered[nrow(ordered), c(
        "statistic",
        "contamination_fraction",
        "contamination_cells",
        "gamma",
        "bias_bound",
        "radius",
        "lower",
        "upper",
        "excludes_null"
      ), drop = FALSE]
    })
  )

  out <- list(
    target_class = object$target_class,
    source_class = object$source_class,
    contamination_fraction_grid = object$contamination_fraction_grid,
    perturbation_constraint = object$perturbation_constraint,
    perturbation_scope = object$perturbation_scope,
    perturbation_layer = object$perturbation_layer,
    bound_type = object$bound_type,
    refit_step = object$refit_step,
    n_perturbation_cells = object$n_perturbation_cells,
    scale_method = object$scale_method,
    scale_value = object$scale_value,
    scale_source = object$scale_source,
    reference = object$reference,
    max_bounds = max_rows
  )

  class(out) <- "summary.cdmc_sensitivity_bounds"
  out
}

print.summary.cdmc_sensitivity_bounds <- function(x, ...) {
  cat("causaldosemc formal sensitivity summary\n")
  cat(sprintf("  target class: %s\n", x$target_class))
  cat(sprintf("  source class: %s\n", x$source_class))
  cat(sprintf("  contamination fractions: %d\n", length(x$contamination_fraction_grid)))
  cat(sprintf("  perturbation constraint: %s\n", x$perturbation_constraint))
  cat(sprintf("  perturbation scope: %s\n", x$perturbation_scope))
  cat(sprintf("  perturbation layer: %s\n", x$perturbation_layer))
  cat(sprintf("  propagation: %s\n", x$bound_type))
  if (identical(x$bound_type, "local_refit")) {
    cat(sprintf("  local refit step: %.6g\n", x$refit_step))
    cat(sprintf("  perturbation cells: %d\n", x$n_perturbation_cells))
  }
  cat(sprintf(
    "  perturbation scale: %s = %.6g (%s)\n",
    x$scale_method,
    x$scale_value,
    x$scale_source
  ))

  if (nrow(x$reference) > 0L) {
    cat("\nReference statistics:\n")
    print(
      x$reference[, c(
        "statistic",
        "estimate",
        "null_value",
        "basis",
        "sensitivity_multiplier",
        "breakdown_bias_bound",
        "breakdown_gamma"
      ), drop = FALSE],
      row.names = FALSE
    )
  }

  if (nrow(x$max_bounds) > 0L) {
    cat("\nLargest requested perturbation bounds:\n")
    print(x$max_bounds, row.names = FALSE)
  }

  invisible(x)
}

print.cdmc_sensitivity_bounds <- function(x, ...) {
  cat("causaldosemc formal sensitivity bounds\n")
  cat(sprintf("  statistics: %d\n", nrow(x$reference)))
  cat(sprintf("  gamma values: %d\n", length(unique(x$bounds$gamma))))
  cat(sprintf("  contamination fractions: %d\n", length(x$contamination_fraction_grid)))
  cat(sprintf("  perturbation constraint: %s\n", x$perturbation_constraint))
  cat(sprintf("  perturbation scope: %s\n", x$perturbation_scope))
  cat(sprintf("  perturbation layer: %s\n", x$perturbation_layer))
  cat(sprintf("  propagation: %s\n", x$bound_type))
  if (identical(x$bound_type, "local_refit")) {
    cat(sprintf("  local refit step: %.6g\n", x$refit_step))
    cat(sprintf("  perturbation cells: %d\n", x$n_perturbation_cells))
  }
  cat(sprintf(
    "  perturbation scale: %s = %.6g (%s)\n",
    x$scale_method,
    x$scale_value,
    x$scale_source
  ))

  print(summary(x), ...)
  invisible(x)
}