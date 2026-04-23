\name{cdmc_lambda_control}
\alias{cdmc_lambda_control}
\alias{print.cdmc_lambda_control}
\title{Bundle lambda and cross-validation tuning options}
\usage{
cdmc_lambda_control(
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
)

\method{print}{cdmc_lambda_control}(x, ...)
}
\arguments{
  \item{lambda}{Optional explicit lambda. If \code{NULL}, lambda is selected
  via the rule chosen by \code{selection}.}

  \item{fraction}{Heuristic-mode shrinkage fraction in \eqn{(0, 1]} applied
  to the data-driven lambda upper bound. Used only when
  \code{selection = "heuristic"}.}

  \item{selection}{One of \code{"cv"} (blocked cross-validation, default) or
  \code{"heuristic"}.}

  \item{grid}{Optional explicit lambda grid for the CV search. When
  \code{NULL}, a log-spaced grid of length \code{nlambda} is constructed
  between \code{min_ratio * lambda_max} and \code{lambda_max}.}

  \item{nlambda}{Length of the auto-generated lambda grid.}

  \item{min_ratio}{Smallest lambda in the auto grid as a fraction of
  \code{lambda_max}.}

  \item{cv_rounds}{Number of CV folds (rounds) drawn over the masked panel.}

  \item{cv_block_size}{Block size used when constructing held-out CV cells.}

  \item{cv_top_k}{If non-\code{NULL}, after the first CV pass refit only the
  top-k lambda candidates on the full mask.}

  \item{cv_coarse_to_fine}{Whether to run a coarse-then-fine CV grid.}

  \item{cv_coarse_nlambda}{Length of the coarse grid in coarse-to-fine mode.}

  \item{cv_warm_starts}{If \code{TRUE}, reuse baseline fit state across
  lambda candidates within each fold (single worker only).}

  \item{x}{A \code{cdmc_lambda_control} object.}

  \item{...}{Unused additional arguments for the print method.}
}
\description{
Constructs a list of low-rank-penalty (\code{lambda}) and blocked-CV options
suitable for passing to \code{\link{cdmc_fit}} or \code{\link{cdmc_dr_fit}}
via the \code{lambda_control} argument. This is purely a packaging
convenience: every field maps one-to-one to an existing scalar argument on
those functions.

When supplied to \code{cdmc_fit()} / \code{cdmc_dr_fit()}, the resolved
control overrides any individual \code{lambda*} / \code{cv_*} / \code{nlambda}
scalars passed alongside it. The legacy scalars remain on the public API for
back-compatibility and continue to work when \code{lambda_control = NULL}
(the default).
}
\value{
An object of class \code{cdmc_lambda_control} (a validated list).
}
\seealso{
\code{\link{cdmc_fit}}, \code{\link{cdmc_dr_fit}}.
}
