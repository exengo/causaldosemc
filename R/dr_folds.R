cdmc_assign_unit_folds <- function(unit_levels, n_folds) {
  n_folds <- as.integer(n_folds)
  if (!is.numeric(n_folds) || length(n_folds) != 1L || n_folds < 2L) {
    stop("n_folds must be an integer greater than or equal to 2.", call. = FALSE)
  }
  if (n_folds > length(unit_levels)) {
    stop("n_folds cannot exceed the number of units.", call. = FALSE)
  }

  shuffled_units <- sample(unit_levels)
  fold_sizes <- rep(length(unit_levels) %/% n_folds, n_folds)
  fold_sizes[seq_len(length(unit_levels) %% n_folds)] <- fold_sizes[seq_len(length(unit_levels) %% n_folds)] + 1L
  fold_id <- integer(length(unit_levels))
  names(fold_id) <- unit_levels
  start_index <- 1L

  for (fold_index in seq_len(n_folds)) {
    end_index <- start_index + fold_sizes[fold_index] - 1L
    current_units <- shuffled_units[seq.int(start_index, end_index)]
    fold_id[current_units] <- fold_index
    start_index <- end_index + 1L
  }

  data.frame(
    unit = unit_levels,
    fold = unname(fold_id[unit_levels]),
    stringsAsFactors = FALSE
  )
}

cdmc_resolve_unit_folds <- function(unit_levels, n_folds, fold_assignments = NULL) {
  if (is.null(fold_assignments)) {
    return(cdmc_assign_unit_folds(unit_levels, n_folds = n_folds))
  }

  if (!is.data.frame(fold_assignments)) {
    stop("fold_assignments must be NULL or a data frame with columns 'unit' and 'fold'.", call. = FALSE)
  }

  required_columns <- c("unit", "fold")
  missing_columns <- setdiff(required_columns, names(fold_assignments))
  if (length(missing_columns) > 0L) {
    stop(
      sprintf(
        "fold_assignments is missing required columns: %s.",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  provided_units <- as.character(fold_assignments$unit)
  if (anyNA(provided_units) || anyDuplicated(provided_units)) {
    stop("fold_assignments must contain each unit exactly once with no missing values.", call. = FALSE)
  }

  unit_keys <- as.character(unit_levels)
  missing_units <- setdiff(unit_keys, provided_units)
  extra_units <- setdiff(provided_units, unit_keys)
  if (length(missing_units) > 0L || length(extra_units) > 0L) {
    message_parts <- character(0)
    if (length(missing_units) > 0L) {
      message_parts <- c(message_parts, sprintf("missing units: %s", paste(missing_units, collapse = ", ")))
    }
    if (length(extra_units) > 0L) {
      message_parts <- c(message_parts, sprintf("unknown units: %s", paste(extra_units, collapse = ", ")))
    }

    stop(
      sprintf("fold_assignments must match the current unit set exactly (%s).", paste(message_parts, collapse = "; ")),
      call. = FALSE
    )
  }

  unit_index <- match(unit_keys, provided_units)
  fold_values_raw <- fold_assignments$fold[unit_index]
  if (anyNA(fold_values_raw) || any(!is.finite(fold_values_raw))) {
    stop("fold_assignments$fold must be finite for every unit.", call. = FALSE)
  }
  if (any(abs(fold_values_raw - round(fold_values_raw)) > sqrt(.Machine$double.eps))) {
    stop("fold_assignments$fold must contain integer-valued fold labels.", call. = FALSE)
  }

  fold_values <- as.integer(round(fold_values_raw))
  if (any(fold_values < 1L)) {
    stop("fold_assignments$fold must contain positive fold labels.", call. = FALSE)
  }

  unique_folds <- sort(unique(fold_values))
  if (!identical(unique_folds, seq_along(unique_folds))) {
    stop("fold_assignments$fold must use consecutive fold labels 1, 2, ..., K.", call. = FALSE)
  }
  if (length(unique_folds) < 2L) {
    stop("fold_assignments must define at least two folds.", call. = FALSE)
  }

  if (!is.null(n_folds)) {
    n_folds <- as.integer(n_folds)
    if (!is.numeric(n_folds) || length(n_folds) != 1L || n_folds < 2L) {
      stop("n_folds must be an integer greater than or equal to 2.", call. = FALSE)
    }
    if (length(unique_folds) != n_folds) {
      stop(
        sprintf(
          "Provided fold_assignments define %d folds, but n_folds = %d.",
          length(unique_folds),
          n_folds
        ),
        call. = FALSE
      )
    }
  }

  data.frame(
    unit = unit_levels,
    fold = fold_values,
    stringsAsFactors = FALSE
  )
}
