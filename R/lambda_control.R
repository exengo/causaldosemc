#' Bundle lambda and cross-validation tuning options
#'
#' Constructs a list of low-rank-penalty (`lambda`) and blocked-CV options
#' suitable for passing to [cdmc_fit()] or [cdmc_dr_fit()] via the
#' `lambda_control` argument. This is purely a packaging convenience: every
#' field maps one-to-one to an existing scalar argument on those functions.
#'
#' When supplied to `cdmc_fit()` / `cdmc_dr_fit()`, the resolved control
#' overrides any individual `lambda*`/`cv_*`/`nlambda` scalars passed
#' alongside it. The legacy scalars remain on the public API for
#' back-compatibility and continue to work when `lambda_control` is `NULL`
#' (the default).
#'
#' @param lambda Optional explicit lambda. If `NULL`, lambda is selected via
#'   the rule chosen by `selection`.
#' @param fraction Heuristic-mode shrinkage fraction in `(0, 1]` applied to
#'   the data-driven lambda upper bound. Used only when `selection = "heuristic"`.
#' @param selection One of `"cv"` (blocked cross-validation, default) or
#'   `"heuristic"`.
#' @param grid Optional explicit lambda grid for the CV search. When `NULL`,
#'   a log-spaced grid of length `nlambda` is constructed between
#'   `min_ratio * lambda_max` and `lambda_max`.
#' @param nlambda Length of the auto-generated lambda grid.
#' @param min_ratio Smallest lambda in the auto grid as a fraction of
#'   `lambda_max`.
#' @param cv_rounds Number of CV folds (rounds) drawn over the masked panel.
#' @param cv_block_size Block size used when constructing held-out CV cells.
#' @param cv_top_k If non-`NULL`, after the first CV pass refit only the top-k
#'   lambda candidates on the full mask.
#' @param cv_coarse_to_fine Whether to run a coarse-then-fine CV grid.
#' @param cv_coarse_nlambda Length of the coarse grid in coarse-to-fine mode.
#' @param cv_warm_starts If `TRUE`, reuse baseline fit state across lambda
#'   candidates within each fold (single worker only).
#'
#' @return An object of class `cdmc_lambda_control` (a validated list).
#'
#' @seealso [cdmc_fit()], [cdmc_dr_fit()].
#' @export
cdmc_lambda_control <- function(
  lambda = NULL,
  fraction = 0.25,
  selection = c("cv", "heuristic"),
  grid = NULL,
  nlambda = 5L,
  min_ratio = 0.05,
  cv_rounds = 5L,
  cv_block_size = 2L,
  cv_top_k = NULL,
  cv_coarse_to_fine = FALSE,
  cv_coarse_nlambda = NULL,
  cv_warm_starts = FALSE
) {
  selection <- match.arg(selection)

  if (!is.null(lambda)) {
    if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
      stop("lambda must be NULL or a single non-negative numeric value.", call. = FALSE)
    }
    lambda <- as.numeric(lambda)
  }
  if (!is.numeric(fraction) || length(fraction) != 1L || !is.finite(fraction) ||
      fraction <= 0 || fraction > 1) {
    stop("fraction must be a scalar in (0, 1].", call. = FALSE)
  }
  if (!is.null(grid)) {
    if (!is.numeric(grid) || length(grid) < 1L || any(!is.finite(grid)) || any(grid < 0)) {
      stop("grid must be NULL or a vector of non-negative finite numerics.", call. = FALSE)
    }
    grid <- as.numeric(grid)
  }
  nlambda <- cdmc_resolve_scalar(nlambda, "nlambda", type = "integer", min = 1L)
  if (!is.numeric(min_ratio) || length(min_ratio) != 1L || !is.finite(min_ratio) ||
      min_ratio <= 0 || min_ratio >= 1) {
    stop("min_ratio must be a scalar in (0, 1).", call. = FALSE)
  }
  cv_rounds <- cdmc_resolve_scalar(cv_rounds, "cv_rounds", type = "integer", min = 1L)
  cv_block_size <- cdmc_resolve_scalar(cv_block_size, "cv_block_size", type = "integer", min = 1L)
  if (!is.null(cv_top_k)) {
    cv_top_k <- cdmc_resolve_scalar(cv_top_k, "cv_top_k", type = "integer", min = 1L,
                                    allow_null = TRUE)
  }
  if (!is.logical(cv_coarse_to_fine) || length(cv_coarse_to_fine) != 1L ||
      is.na(cv_coarse_to_fine)) {
    stop("cv_coarse_to_fine must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(cv_coarse_nlambda)) {
    cv_coarse_nlambda <- cdmc_resolve_scalar(cv_coarse_nlambda, "cv_coarse_nlambda",
                                             type = "integer", min = 3L,
                                             allow_null = TRUE)
  }
  if (!is.logical(cv_warm_starts) || length(cv_warm_starts) != 1L ||
      is.na(cv_warm_starts)) {
    stop("cv_warm_starts must be TRUE or FALSE.", call. = FALSE)
  }

  structure(
    list(
      lambda = lambda,
      fraction = as.numeric(fraction),
      selection = selection,
      grid = grid,
      nlambda = nlambda,
      min_ratio = as.numeric(min_ratio),
      cv_rounds = cv_rounds,
      cv_block_size = cv_block_size,
      cv_top_k = cv_top_k,
      cv_coarse_to_fine = isTRUE(cv_coarse_to_fine),
      cv_coarse_nlambda = cv_coarse_nlambda,
      cv_warm_starts = isTRUE(cv_warm_starts)
    ),
    class = "cdmc_lambda_control"
  )
}

# Internal: returns a fully-validated cdmc_lambda_control. If `lambda_control`
# is non-NULL, it is validated/coerced and returned as-is. Otherwise a control
# is built from the legacy scalars passed by the caller. Never errors on the
# legacy path beyond what cdmc_lambda_control() itself already validates.
cdmc_resolve_lambda_control <- function(
  lambda_control,
  lambda,
  fraction,
  selection,
  grid,
  nlambda,
  min_ratio,
  cv_rounds,
  cv_block_size,
  cv_top_k,
  cv_coarse_to_fine,
  cv_coarse_nlambda,
  cv_warm_starts
) {
  if (is.null(lambda_control)) {
    return(cdmc_lambda_control(
      lambda = lambda,
      fraction = fraction,
      selection = selection,
      grid = grid,
      nlambda = nlambda,
      min_ratio = min_ratio,
      cv_rounds = cv_rounds,
      cv_block_size = cv_block_size,
      cv_top_k = cv_top_k,
      cv_coarse_to_fine = cv_coarse_to_fine,
      cv_coarse_nlambda = cv_coarse_nlambda,
      cv_warm_starts = cv_warm_starts
    ))
  }

  if (!inherits(lambda_control, "cdmc_lambda_control")) {
    if (!is.list(lambda_control)) {
      stop("lambda_control must be NULL or an object built by cdmc_lambda_control().",
           call. = FALSE)
    }
    # Permit a bare list with the same field names by re-running the constructor.
    return(do.call(cdmc_lambda_control, lambda_control))
  }

  lambda_control
}

#' @export
print.cdmc_lambda_control <- function(x, ...) {
  cat("<cdmc_lambda_control>\n")
  cat(sprintf("  selection       : %s\n", x$selection))
  cat(sprintf("  lambda          : %s\n",
              if (is.null(x$lambda)) "NULL (data-driven)" else format(x$lambda)))
  cat(sprintf("  fraction        : %s\n", format(x$fraction)))
  cat(sprintf("  grid            : %s\n",
              if (is.null(x$grid)) "NULL (auto)"
              else sprintf("user-supplied length %d", length(x$grid))))
  cat(sprintf("  nlambda         : %d\n", x$nlambda))
  cat(sprintf("  min_ratio       : %s\n", format(x$min_ratio)))
  cat(sprintf("  cv_rounds       : %d\n", x$cv_rounds))
  cat(sprintf("  cv_block_size   : %d\n", x$cv_block_size))
  cat(sprintf("  cv_top_k        : %s\n",
              if (is.null(x$cv_top_k)) "NULL" else as.character(x$cv_top_k)))
  cat(sprintf("  cv_coarse_to_fine : %s\n", as.character(x$cv_coarse_to_fine)))
  cat(sprintf("  cv_coarse_nlambda : %s\n",
              if (is.null(x$cv_coarse_nlambda)) "NULL"
              else as.character(x$cv_coarse_nlambda)))
  cat(sprintf("  cv_warm_starts  : %s\n", as.character(x$cv_warm_starts)))
  invisible(x)
}
