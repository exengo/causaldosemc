\name{simulate_cdmc_data}
\alias{simulate_cdmc_data}
\title{Simulate balanced panel data for causaldosemc}
\description{
Generate a synthetic balanced panel with a low-rank untreated outcome surface,
reversible doses, optional covariates, and a linear
contemporaneous or distributed-lag treatment effect.
}
\usage{
simulate_cdmc_data(
  n_units = 20L,
  n_times = 12L,
  rank = 2L,
  beta = 1,
  lag_beta = NULL,
  n_covariates = 1L,
  noise_sd = 0.1,
  switch_on_prob = 0.2,
  switch_off_prob = 0.35,
  signed_dose = FALSE,
  negative_dose_prob = 0.5,
  seed = NULL
)
}
\arguments{
  \item{n_units}{Number of units in the simulated balanced panel.}

  \item{n_times}{Number of time periods in the simulated balanced panel.}

  \item{rank}{Rank of the latent untreated factor structure.}

  \item{beta}{Coefficient on the contemporaneous dose, \code{dose_lag0}.}

  \item{lag_beta}{Optional numeric vector of lag coefficients. For example,
  \code{lag_beta = c(0.4, 0.2)} adds coefficients for \code{dose_lag1} and
  \code{dose_lag2}.}

  \item{n_covariates}{Number of simulated covariates to append as columns named
  \code{x1}, \code{x2}, and so on.}

  \item{noise_sd}{Standard deviation of the additive outcome noise.}

  \item{switch_on_prob}{Probability of switching from zero dose to active
  exposure at each time point when the unit is currently untreated.}

  \item{switch_off_prob}{Probability of switching from active exposure back to 0
  at each time point when the unit is currently treated.}

  \item{signed_dose}{Logical indicating whether active doses may be negative as
  well as positive.}

  \item{negative_dose_prob}{When \code{signed_dose = TRUE}, the probability
  that an active dose draw is assigned a negative sign.}

  \item{seed}{Optional random seed used before simulation.}
}
\details{
The function is intended for package examples, smoke tests, and development
checks. It produces exact 0 doses when a unit is untreated and nonzero doses
when a unit is exposed. The first period is always simulated at dose 0 for
every unit so washout-sensitive workflows retain at least one eligible baseline
control cell per unit. When \code{signed_dose = TRUE}, active doses may be
positive or negative. The outcome is the sum of a low-rank untreated baseline,
additive unit and time effects, optional observed covariates, a linear dose and
lag effect, and Gaussian noise.
}
\value{
A data frame with columns \code{unit}, \code{time}, \code{y}, \code{dose}, and
optional covariate columns \code{x1}, \code{x2}, and so on.
}
\examples{
panel <- simulate_cdmc_data(
  n_units = 6,
  n_times = 5,
  beta = 1.1,
  lag_beta = 0.3,
  seed = 1
)

head(panel)
}
