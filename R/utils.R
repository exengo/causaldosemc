cdmc_assert_installed <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      sprintf("Package '%s' must be installed to use causaldosemc.", pkg),
      call. = FALSE
    )
  }
}

cdmc_relative_change <- function(new_matrix, old_matrix) {
  numerator <- sqrt(sum((new_matrix - old_matrix) ^ 2))
  denominator <- sqrt(sum(old_matrix ^ 2))

  if (denominator < .Machine$double.eps) {
    return(numerator)
  }

  numerator / denominator
}

cdmc_reconstruct_low_rank <- function(fit, n_units, n_times) {
  if (length(fit$d) == 0L || is.null(fit$u) || is.null(fit$v)) {
    return(matrix(0, nrow = n_units, ncol = n_times))
  }

  diagonal <- diag(fit$d, nrow = length(fit$d), ncol = length(fit$d))
  fit$u %*% diagonal %*% t(fit$v)
}

cdmc_flatten_matrix <- function(x) {
  as.vector(t(x))
}

cdmc_inherits_any <- function(object, classes) {
  any(vapply(classes, function(class_name) inherits(object, class_name), logical(1)))
}

cdmc_first_inherited_class <- function(object, classes) {
  matches <- vapply(classes, function(class_name) inherits(object, class_name), logical(1))
  if (!any(matches)) {
    return(NULL)
  }

  classes[[which(matches)[1L]]]
}

cdmc_named_numeric_data_frame <- function(column_names, n_rows = 1L, fill = 0) {
  out <- data.frame(
    matrix(fill, nrow = as.integer(n_rows), ncol = length(column_names)),
    check.names = FALSE
  )
  names(out) <- column_names
  out
}

cdmc_sample_array_matrix <- function(array_3d, sample_indices) {
  slice_names <- dimnames(array_3d)[[3L]]
  n_units <- dim(array_3d)[1L]
  n_times <- dim(array_3d)[2L]
  n_slices <- dim(array_3d)[3L]

  if (is.null(sample_indices) || nrow(sample_indices) == 0L) {
    out <- matrix(numeric(0), nrow = 0L, ncol = n_slices)
    colnames(out) <- slice_names
    return(out)
  }

  flat_indices <- sample_indices[, 1L] + (sample_indices[, 2L] - 1L) * n_units
  out <- matrix(array_3d, nrow = n_units * n_times, ncol = n_slices)[flat_indices, , drop = FALSE]
  colnames(out) <- slice_names
  out
}

cdmc_build_sample_history <- function(lag_array, sample_indices, tau = NULL, include_indices = TRUE) {
  if (is.null(sample_indices)) {
    sample_indices <- matrix(integer(0), ncol = 2L)
  }

  lag_history <- data.frame(
    cdmc_sample_array_matrix(lag_array, sample_indices),
    check.names = FALSE
  )

  components <- list()
  if (include_indices) {
    components$row <- sample_indices[, 1L]
    components$col <- sample_indices[, 2L]
  }
  if (!is.null(tau)) {
    components$tau <- tau
  }

  if (length(components) == 0L) {
    return(lag_history)
  }

  history <- data.frame(components, stringsAsFactors = FALSE, check.names = FALSE)
  cbind(history, lag_history)
}

`%||%` <- function(x, y) {
  if (is.null(x)) {
    y
  } else {
    x
  }
}

cdmc_progress_bar <- function(label, total) {
  if (!requireNamespace("cli", quietly = TRUE)) return(NULL)
  cli::cli_progress_bar(label, total = total, clear = FALSE, .auto_close = FALSE)
}

cdmc_progress_update <- function(pb) {
  if (!is.null(pb)) cli::cli_progress_update(id = pb)
  invisible(NULL)
}

cdmc_progress_done <- function(pb) {
  if (!is.null(pb)) cli::cli_progress_done(id = pb)
  invisible(NULL)
}
