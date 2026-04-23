\name{summary.cdmc_dr_fit}
\alias{summary.cdmc_dr_fit}
\alias{print.summary.cdmc_dr_fit}
\title{Summarise a cross-fitted DR fit}
\description{
Returns a structured summary of \code{\link{cdmc_dr_fit}} output: the
dose-lag coefficients, fold-level lambda choices, weight diagnostics, and
the score (\code{aipw} / \code{plr}) plus train baseline weighting setting.
For inferential intervals, pair this with \code{\link{cdmc_bootstrap}}.
}
\usage{
\method{summary}{cdmc_dr_fit}(object, ...)
}
\arguments{
  \item{object}{A \code{cdmc_dr_fit} object.}
  \item{...}{Unused.}
}
\value{
A \code{summary.cdmc_dr_fit} object (a list with \code{coefficients},
\code{folds}, \code{weight}, \code{score}, and \code{dimensions}
components) printed via the registered \code{print} method.
}
\seealso{\code{\link{cdmc_dr_fit}}, \code{\link{cdmc_dr_score_contrast}}}
