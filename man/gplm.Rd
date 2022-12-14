% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gplm.R
\name{gplm}
\alias{gplm}
\title{Fit a gplm model (spBayes spLM)}
\usage{
gplm(
  formula,
  data,
  coords,
  chains = 4,
  n_batch = 200,
  batch_len = 100,
  cov_model = c("gaussian", "exponential", "matern", "spherical"),
  starting = list(phi = sqrt(3)/1, sigma.sq = 1, tau.sq = 1),
  priors = list(phi.unif = c(sqrt(3)/52, sqrt(3)/1), sigma.sq.ig = c(2, 1), tau.sq.ig =
    c(2, 1)),
  tuning = list(phi = 1, sigma.sq = 1, tau.sq = 1),
  burnin_frac = 0.5,
  accept_rate = 0.43,
  thin = 1,
  n_report = 50,
  verbose = FALSE
)
}
\arguments{
\item{formula}{A formula object. A symbolic description of the model to be fitted.}

\item{data}{A data.frame object containing data of all variables used in the model.}

\item{coords}{Either a numeric matrix of coordinates or quoted column names from \code{data} to use for distance calculations.}

\item{chains}{Numeric. Number of MCMC chains to fit.}

\item{n_batch}{Numeric. Number of adaptive batches to use when fitting each chain.}

\item{batch_len}{Numeric. Number of iterations per batch.}

\item{cov_model}{Character. Name of the covariance model to use, supported values are:
"exponential", "matern", "spherical", and "gaussian"}

\item{starting}{Named list of parameter starting values, allowed names: \code{beta}, \code{sigma.sq}, \code{tau.sq}, \code{phi}, and \code{nu}.}

\item{tuning}{Names list of variance values for the MH sampler.}

\item{burnin_frac}{Numeric. Proportion of iterations to discard as burnin.}

\item{accept_rate}{Numeric. Desired acceptance rate used by the adaptive MCMC algorithm}

\item{thin}{Numeric. Amount of thinning to apply to posterior samples.}

\item{n_report}{Numeric. The interval to report Metropolis sampler acceptance and MCMC progress.}

\item{verbose}{Logical. Should verbose output (sampling progress) be printed.}

\item{prior}{Named list of priors, each value is a vector of prior hyperparameter values. See \code{\link[spBayes:spLM]{spBayes::spLM()}} for details.}
}
\description{
A more user friendly (and limited) wrapper around the spBayes' spLM
model.
}
\examples{
m = gplm(
  avg_temp~1, data = avg_temp, coords = "week",
  starting=list(
    "phi"=sqrt(3)/4, "sigma.sq"=1, "tau.sq"=1
  ),
  tuning=list(
    "phi"=1, "sigma.sq"=1, "tau.sq"=1
  ),
  priors=list(
    "phi.unif"=c(sqrt(3)/52, sqrt(3)/1),
    "sigma.sq.ig"=c(2, 1),
    "tau.sq.ig"=c(2, 1)
  ),
  thin=10,
  n_batch = 100,
  batch_len = 50
)

newdata = data.frame(
  week = seq(0,3.5*52) |> jitter()
)

(pred = predict(m, newdata=newdata, coords = "week"))

}
