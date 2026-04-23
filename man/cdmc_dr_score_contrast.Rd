\name{cdmc_dr_score_contrast}
\alias{cdmc_dr_score_contrast}
\title{Side-by-side comparison of two DR scores}
\description{
Compares dose-lag coefficients across two \code{\link{cdmc_dr_fit}}
objects fitted on the same panel with different \code{dr_score} values
(typically \code{"aipw"} and \code{"plr"}). Both fits should otherwise
share data, lag order, weight method, and folds; differing scores let the
user assess robustness of the headline estimates to the orthogonal score
choice.
}
\usage{
cdmc_dr_score_contrast(aipw, plr)
}
\arguments{
  \item{aipw, plr}{\code{cdmc_dr_fit} objects to contrast. The names refer
  to the typical use; the function does not enforce a particular
  \code{dr_score} value, so it can also be used to compare any two
  compatible fits.}
}
\value{
A data frame with one row per dose-lag term containing both estimates,
their absolute and relative differences, and a flag for sign agreement.
Sign agreement uses a \code{0} tolerance equal to
\code{sqrt(.Machine$double.eps) * max(|estimates|)}.
}
\seealso{\code{\link{cdmc_dr_fit}}, \code{\link{summary.cdmc_dr_fit}}}
