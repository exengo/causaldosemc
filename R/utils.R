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

`%||%` <- function(x, y) {
  if (is.null(x)) {
    y
  } else {
    x
  }
}
