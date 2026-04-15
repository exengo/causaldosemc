cdmc_prepare_effect_sample <- function(
  object,
  lag_order = object$lag_order,
  include_zero_dose = TRUE
) {
  lag_array <- cdmc_build_lagged_doses(object$dose_matrix, lag_order = lag_order)
  valid_history <- apply(!is.na(lag_array), c(1, 2), all) & !is.na(object$effect$tau)
  exposure_magnitude <- cdmc_lag_exposure_magnitude(lag_array)

  sample_mask <- if (include_zero_dose) {
    valid_history
  } else {
    valid_history & exposure_magnitude > 0
  }

  sample_indices <- which(sample_mask, arr.ind = TRUE)
  history <- cdmc_build_sample_history(
    lag_array = lag_array,
    sample_indices = sample_indices,
    tau = object$effect$tau[sample_indices]
  )

  list(
    lag_order = lag_order,
    lag_names = dimnames(lag_array)[[3L]],
    sample_mask = sample_mask,
    history = history
  )
}

cdmc_resolve_effect_weights <- function(weights, object, sample_mask) {
  if (is.null(weights) && !is.null(object$fit_control$weights)) {
    weights <- object$fit_control$weights
  }

  if (is.null(weights)) {
    return(NULL)
  }

  resolved <- cdmc_resolve_weight_vector(
    weights = weights,
    data = object$data,
    n_units = nrow(object$dose_matrix),
    n_times = ncol(object$dose_matrix),
    argument_name = "weights"
  )

  resolved[cdmc_flatten_matrix(sample_mask)]
}

cdmc_build_spline_basis_spec <- function(x, df) {
  boundary_knots <- range(x)
  interior_candidates <- sort(unique(x[x > boundary_knots[1L] & x < boundary_knots[2L]]))
  n_knots <- min(max(df - 1L, 0L), length(interior_candidates))

  if (n_knots < 1L) {
    return(list(
      type = "linear",
      zero_basis = 0,
      ncol = 1L
    ))
  }

  knot_probs <- seq(0, 1, length.out = n_knots + 2L)[seq_len(n_knots) + 1L]
  knots <- as.numeric(stats::quantile(interior_candidates, probs = knot_probs, names = FALSE, type = 1L))
  knots <- unique(knots[knots > boundary_knots[1L] & knots < boundary_knots[2L]])

  if (length(knots) < 1L) {
    return(list(
      type = "linear",
      zero_basis = 0,
      ncol = 1L
    ))
  }

  zero_basis <- as.numeric(splines::ns(0, knots = knots, Boundary.knots = boundary_knots))

  list(
    type = "ns",
    knots = knots,
    boundary_knots = boundary_knots,
    zero_basis = zero_basis,
    ncol = length(zero_basis)
  )
}

cdmc_apply_spline_basis <- function(x, basis_spec) {
  if (identical(basis_spec$type, "linear")) {
    transformed <- matrix(x, ncol = 1L)
  } else {
    transformed <- splines::ns(
      x,
      knots = basis_spec$knots,
      Boundary.knots = basis_spec$boundary_knots
    )
  }

  sweep(transformed, 2, basis_spec$zero_basis, FUN = "-")
}

cdmc_build_response_design <- function(history, model = c("linear", "spline"), df = 4L, basis_spec = NULL) {
  model <- match.arg(model)
  lag_names <- setdiff(names(history), c("row", "col", "tau"))

  if (model == "linear") {
    design <- as.matrix(history[, lag_names, drop = FALSE])
    return(list(
      design = design,
      basis_spec = NULL,
      column_names = colnames(design),
      model = model,
      lag_names = lag_names
    ))
  }

  df <- as.integer(df)
  if (df < 1L) {
    stop("df must be a positive integer for spline dose-response models.", call. = FALSE)
  }

  basis_spec <- basis_spec %||% vector("list", length(lag_names))
  names(basis_spec) <- lag_names
  design_parts <- vector("list", length(lag_names))
  column_names <- character(0)

  for (index in seq_along(lag_names)) {
    lag_name <- lag_names[[index]]
    x <- history[[lag_name]]

    if (is.null(basis_spec[[lag_name]])) {
      basis_spec[[lag_name]] <- cdmc_build_spline_basis_spec(x, df = df)
    }

    centered <- cdmc_apply_spline_basis(x, basis_spec[[lag_name]])
    current_names <- paste0(lag_name, "_ns", seq_len(ncol(centered)))
    colnames(centered) <- current_names

    design_parts[[index]] <- centered
    column_names <- c(column_names, current_names)
  }

  design <- do.call(cbind, design_parts)

  list(
    design = design,
    basis_spec = basis_spec,
    column_names = column_names,
    model = model,
    lag_names = lag_names,
    df = df
  )
}

cdmc_resolve_gam_response_df <- function(df) {
  if (!is.numeric(df) || length(df) != 1L || !is.finite(df) || df < 3L) {
    stop("df must be a single integer greater than or equal to 3 for gam dose-response models.", call. = FALSE)
  }

  as.integer(df)
}

cdmc_resolve_gam_response_k <- function(x, df) {
  n_unique <- length(unique(stats::na.omit(x)))
  if (n_unique < 4L) {
    return(NULL)
  }

  k <- min(as.integer(df), n_unique - 1L)
  if (k < 3L) {
    return(NULL)
  }

  k
}

cdmc_build_gam_response_formula <- function(history, lag_names, df = 4L) {
  cdmc_assert_installed("mgcv")

  termlabels <- vapply(lag_names, function(lag_name) {
    k <- cdmc_resolve_gam_response_k(history[[lag_name]], df = df)
    if (is.null(k)) {
      lag_name
    } else {
      sprintf("s(%s, k = %d)", lag_name, k)
    }
  }, character(1))

  formula_environment <- new.env(parent = baseenv())
  formula_environment$s <- mgcv::s
  stats::as.formula(
    sprintf("tau ~ 0 + %s", paste(termlabels, collapse = " + ")),
    env = formula_environment
  )
}

cdmc_fit_gam_dose_response <- function(history, lag_names, df = 4L, weights = NULL) {
  cdmc_assert_installed("mgcv")
  df <- cdmc_resolve_gam_response_df(df)

  formula <- cdmc_build_gam_response_formula(history, lag_names = lag_names, df = df)
  gam_data <- history
  fit <- if (is.null(weights)) {
    mgcv::gam(
      formula = formula,
      data = gam_data,
      method = "REML"
    )
  } else {
    gam_data$.cdmc_dose_response_weight <- as.numeric(weights)
    .cdmc_dose_response_weight <- gam_data$.cdmc_dose_response_weight
    mgcv::gam(
      formula = formula,
      data = gam_data,
      weights = .cdmc_dose_response_weight,
      method = "REML"
    )
  }

  coefficients <- stats::coef(fit)
  coefficients[is.na(coefficients)] <- 0
  zero_history <- cdmc_named_numeric_data_frame(lag_names)
  zero_reference <- as.numeric(stats::predict(fit, newdata = zero_history, type = "response"))
  fitted_values <- as.numeric(stats::predict(fit, newdata = history[, lag_names, drop = FALSE], type = "response")) - zero_reference

  list(
    coefficients = coefficients,
    fitted_values = fitted_values,
    residuals = history$tau - fitted_values,
    weights = weights,
    model_fit = fit,
    zero_reference = zero_reference,
    design_info = list(
      design = NULL,
      basis_spec = NULL,
      column_names = names(coefficients),
      model = "gam",
      lag_names = lag_names,
      df = df,
      formula = formula
    )
  )
}

cdmc_resolve_response_forest_trees <- function(forest_trees) {
  if (!is.numeric(forest_trees) || length(forest_trees) != 1L || !is.finite(forest_trees) || forest_trees < 1) {
    stop("forest_trees must be a positive integer.", call. = FALSE)
  }

  as.integer(forest_trees)
}

cdmc_resolve_response_forest_mtry <- function(forest_mtry) {
  if (is.null(forest_mtry)) {
    return(NULL)
  }

  if (!is.numeric(forest_mtry) || length(forest_mtry) != 1L || !is.finite(forest_mtry) || forest_mtry < 1) {
    stop("forest_mtry must be NULL or a positive integer.", call. = FALSE)
  }

  as.integer(forest_mtry)
}

cdmc_resolve_response_forest_min_node_size <- function(forest_min_node_size) {
  if (is.null(forest_min_node_size)) {
    return(NULL)
  }

  if (!is.numeric(forest_min_node_size) || length(forest_min_node_size) != 1L || !is.finite(forest_min_node_size) || forest_min_node_size < 1) {
    stop("forest_min_node_size must be NULL or a positive integer.", call. = FALSE)
  }

  as.integer(forest_min_node_size)
}

cdmc_forest_response_importance <- function(fit, lag_names) {
  importance <- numeric(length(lag_names))
  names(importance) <- lag_names

  fitted_importance <- fit$variable.importance
  if (is.null(fitted_importance) || length(fitted_importance) == 0L) {
    return(importance)
  }

  common_names <- intersect(names(fitted_importance), lag_names)
  importance[common_names] <- as.numeric(fitted_importance[common_names])
  importance
}

cdmc_fit_forest_dose_response <- function(
  history,
  lag_names,
  weights = NULL,
  forest_trees = 200L,
  forest_mtry = NULL,
  forest_min_node_size = NULL
) {
  cdmc_assert_installed("ranger")

  forest_trees <- cdmc_resolve_response_forest_trees(forest_trees)
  forest_mtry <- cdmc_resolve_response_forest_mtry(forest_mtry)
  forest_min_node_size <- cdmc_resolve_response_forest_min_node_size(forest_min_node_size)

  forest_data <- history[, c("tau", lag_names), drop = FALSE]
  fit <- ranger::ranger(
    dependent.variable.name = "tau",
    data = forest_data,
    case.weights = weights,
    num.trees = forest_trees,
    mtry = forest_mtry,
    min.node.size = forest_min_node_size,
    importance = "impurity",
    respect.unordered.factors = "order",
    seed = 1L,
    num.threads = 1L
  )

  zero_history <- cdmc_named_numeric_data_frame(lag_names)
  zero_reference <- as.numeric(stats::predict(fit, data = zero_history)$predictions)
  fitted_values <- as.numeric(stats::predict(fit, data = history[, lag_names, drop = FALSE])$predictions) - zero_reference

  list(
    coefficients = cdmc_forest_response_importance(fit, lag_names = lag_names),
    fitted_values = fitted_values,
    residuals = history$tau - fitted_values,
    weights = weights,
    model_fit = fit,
    zero_reference = zero_reference,
    forest_trees = forest_trees,
    forest_mtry = forest_mtry,
    forest_min_node_size = forest_min_node_size,
    design_info = list(
      design = NULL,
      basis_spec = NULL,
      column_names = lag_names,
      model = "forest",
      lag_names = lag_names,
      df = NULL,
      forest_trees = forest_trees,
      forest_mtry = forest_mtry,
      forest_min_node_size = forest_min_node_size
    )
  )
}

cdmc_build_tree_response_formula <- function(lag_names) {
  stats::reformulate(termlabels = lag_names, response = "tau")
}

cdmc_tree_response_importance <- function(fit, lag_names) {
  importance <- numeric(length(lag_names))
  names(importance) <- lag_names

  fitted_importance <- fit$variable.importance
  if (is.null(fitted_importance) || length(fitted_importance) == 0L) {
    return(importance)
  }

  common_names <- intersect(names(fitted_importance), lag_names)
  importance[common_names] <- as.numeric(fitted_importance[common_names])
  importance
}

cdmc_fit_tree_dose_response <- function(history, lag_names, weights = NULL) {
  cdmc_assert_installed("rpart")

  formula <- cdmc_build_tree_response_formula(lag_names)
  tree_data <- history[, c("tau", lag_names), drop = FALSE]
  fit <- if (is.null(weights)) {
    rpart::rpart(
      formula = formula,
      data = tree_data,
      method = "anova"
    )
  } else {
    rpart::rpart(
      formula = formula,
      data = tree_data,
      weights = weights,
      method = "anova"
    )
  }

  zero_history <- cdmc_named_numeric_data_frame(lag_names)
  zero_reference <- as.numeric(stats::predict(fit, newdata = zero_history))
  fitted_values <- as.numeric(stats::predict(fit, newdata = history[, lag_names, drop = FALSE])) - zero_reference

  list(
    coefficients = cdmc_tree_response_importance(fit, lag_names = lag_names),
    fitted_values = fitted_values,
    residuals = history$tau - fitted_values,
    weights = weights,
    model_fit = fit,
    zero_reference = zero_reference,
    design_info = list(
      design = NULL,
      basis_spec = NULL,
      column_names = lag_names,
      model = "tree",
      lag_names = lag_names,
      df = NULL,
      formula = formula
    )
  )
}

cdmc_fit_weighted_regression <- function(design, response, weights = NULL) {
  if (is.null(weights)) {
    fit <- stats::lm.fit(x = design, y = response)
  } else {
    fit <- stats::lm.wfit(x = design, y = response, w = weights)
  }

  coefficients <- fit$coefficients
  coefficients[is.na(coefficients)] <- 0
  names(coefficients) <- colnames(design)
  fitted_values <- as.vector(design %*% coefficients)

  list(
    coefficients = coefficients,
    fitted_values = fitted_values,
    residuals = response - fitted_values,
    weights = weights
  )
}

cdmc_dose_response <- function(
  object,
  model = c("linear", "spline", "gam", "tree", "forest"),
  lag_order = object$lag_order,
  df = 4L,
  forest_trees = 200L,
  forest_mtry = NULL,
  forest_min_node_size = NULL,
  include_zero_dose = TRUE,
  weights = NULL
) {
  if (!inherits(object, "cdmc_fit")) {
    stop("object must inherit from 'cdmc_fit'.", call. = FALSE)
  }

  model <- match.arg(model)
  prepared <- cdmc_prepare_effect_sample(
    object = object,
    lag_order = lag_order,
    include_zero_dose = include_zero_dose
  )
  observation_weights <- cdmc_resolve_effect_weights(
    weights = weights,
    object = object,
    sample_mask = prepared$sample_mask
  )
  resolved_forest_trees <- NULL
  resolved_forest_mtry <- NULL
  resolved_forest_min_node_size <- NULL
  if (identical(model, "gam")) {
    fit <- cdmc_fit_gam_dose_response(
      history = prepared$history,
      lag_names = prepared$lag_names,
      df = df,
      weights = observation_weights
    )
    design_info <- fit$design_info
  } else if (identical(model, "tree")) {
    fit <- cdmc_fit_tree_dose_response(
      history = prepared$history,
      lag_names = prepared$lag_names,
      weights = observation_weights
    )
    design_info <- fit$design_info
  } else if (identical(model, "forest")) {
    fit <- cdmc_fit_forest_dose_response(
      history = prepared$history,
      lag_names = prepared$lag_names,
      weights = observation_weights,
      forest_trees = forest_trees,
      forest_mtry = forest_mtry,
      forest_min_node_size = forest_min_node_size
    )
    design_info <- fit$design_info
    resolved_forest_trees <- fit$forest_trees
    resolved_forest_mtry <- fit$forest_mtry
    resolved_forest_min_node_size <- fit$forest_min_node_size
  } else {
    design_info <- cdmc_build_response_design(
      history = prepared$history,
      model = model,
      df = df
    )
    fit <- cdmc_fit_weighted_regression(
      design = design_info$design,
      response = prepared$history$tau,
      weights = observation_weights
    )
  }

  result <- list(
    call = match.call(),
    model = model,
    lag_order = lag_order,
    include_zero_dose = include_zero_dose,
    weights = weights,
    fit_object = object,
    history = prepared$history,
    design_info = design_info,
    coefficients = fit$coefficients,
    fitted_values = fit$fitted_values,
    residuals = fit$residuals,
    model_fit = fit$model_fit %||% NULL,
    zero_reference = fit$zero_reference %||% 0,
    forest_trees = resolved_forest_trees,
    forest_mtry = resolved_forest_mtry,
    forest_min_node_size = resolved_forest_min_node_size
  )

  class(result) <- "cdmc_dose_response"
  result
}

cdmc_prepare_prediction_history <- function(object, dose = NULL, history = NULL) {
  lag_names <- object$design_info$lag_names

  if (!is.null(history)) {
    history <- as.data.frame(history)
    missing_columns <- setdiff(lag_names, names(history))
    if (length(missing_columns) > 0L) {
      stop(
        sprintf("history is missing required lag columns: %s.", paste(missing_columns, collapse = ", ")),
        call. = FALSE
      )
    }
    return(history[, lag_names, drop = FALSE])
  }

  if (is.null(dose)) {
    stop("Provide either dose or history when predicting from a dose-response fit.", call. = FALSE)
  }

  prediction_history <- cdmc_named_numeric_data_frame(lag_names, n_rows = length(dose))
  prediction_history[["dose_lag0"]] <- dose
  prediction_history
}

cdmc_predict_dose_response_response <- function(object, prediction_history) {
  if (identical(object$model, "gam")) {
    return(
      as.numeric(stats::predict(object$model_fit, newdata = prediction_history, type = "response")) -
        (object$zero_reference %||% 0)
    )
  }

  if (identical(object$model, "tree")) {
    return(
      as.numeric(stats::predict(object$model_fit, newdata = prediction_history)) -
        (object$zero_reference %||% 0)
    )
  }

  if (identical(object$model, "forest")) {
    return(
      as.numeric(stats::predict(object$model_fit, data = prediction_history)$predictions) -
        (object$zero_reference %||% 0)
    )
  }

  design_info <- cdmc_build_response_design(
    history = transform(prediction_history, row = seq_len(nrow(prediction_history)), col = seq_len(nrow(prediction_history)), tau = 0),
    model = object$model,
    df = object$design_info$df %||% 1L,
    basis_spec = object$design_info$basis_spec
  )

  as.vector(design_info$design %*% object$coefficients)
}

predict.cdmc_dose_response <- function(object, dose = NULL, history = NULL, type = c("response", "slope"), eps = 1e-4, ...) {
  type <- match.arg(type)
  prediction_history <- cdmc_prepare_prediction_history(object, dose = dose, history = history)
  response <- cdmc_predict_dose_response_response(object, prediction_history)
  if (type == "response") {
    output <- prediction_history
    output$estimate <- response
    return(output)
  }

  history_up <- prediction_history
  history_down <- prediction_history
  history_up$dose_lag0 <- history_up$dose_lag0 + eps
  history_down$dose_lag0 <- history_down$dose_lag0 - eps

  estimate_up <- cdmc_predict_dose_response_response(object, history_up)
  estimate_down <- cdmc_predict_dose_response_response(object, history_down)
  slope <- (estimate_up - estimate_down) / (history_up$dose_lag0 - history_down$dose_lag0)

  output <- prediction_history
  output$estimate <- slope
  output
}

print.cdmc_dose_response <- function(x, ...) {
  cat("causaldosemc dose-response fit\n")
  cat(sprintf("  model: %s\n", x$model))
  cat(sprintf("  lag order: %d\n", x$lag_order))
  cat(sprintf("  sample size: %d\n", nrow(x$history)))
  cat(sprintf("  include zero-dose sample: %s\n", if (x$include_zero_dose) "yes" else "no"))
  if (length(x$coefficients) > 0L) {
    cat(sprintf("  coefficient count: %d\n", length(x$coefficients)))
  }

  invisible(x)
}
