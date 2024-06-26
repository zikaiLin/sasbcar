% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sas_bcar_fit.R
\name{fit_model}
\alias{fit_model}
\title{Fit Model using Rcpp}
\usage{
fit_model(
  X,
  y,
  partition = NULL,
  radius_candidats = c(1, 2),
  radius_update_thinning = 50,
  burnin = 2000,
  thinning = 1,
  mcmc_sample = 2000,
  prior_max = 1,
  prior_min = 0,
  threshold = 0.6,
  verbose = 200,
  kernel_poly_degree = 20,
  kernel_concentration = 0.01,
  kernel_smoothness = 50,
  seed = 1234
)
}
\arguments{
\item{X}{A matrix of predictors.}

\item{y}{A vector of response.}

\item{partition}{A vector of partition, defaults to NULL.}

\item{radius_candidats}{A vector of candidate radius parameters, defaults to c(1,2).}

\item{radius_update_thinning}{An integer specifying thinning interval for radius updates, defaults to 50.}

\item{burnin}{An integer specifying the number of burn-in iterations, defaults to 2000.}

\item{thinning}{An integer specifying the thinning interval for MCMC sampling, defaults to 1.}

\item{mcmc_sample}{An integer specifying the number of MCMC samples to draw, defaults to 2000.}

\item{prior_max}{A numeric value specifying the maximum prior value of probability map, defaults to 1.}

\item{prior_min}{A numeric value specifying the minimum prior value of probability map, defaults to 0.}

\item{threshold}{A numeric value specifying the threshold for inclusion probabilities, defaults to 0.6.}

\item{verbose}{An integer specifying the verbosity at MCMC sampling steps, defaults to 200 steps.}

\item{kernel_poly_degree}{An integer specifying the polynomial degree for the squared exponential kernel, defaults to 20.}

\item{kernel_concentration}{A numeric value specifying the concentration parameter for the squared exponential kernel, defaults to 0.01.}

\item{kernel_smoothness}{A numeric value specifying the smoothness parameter for the squared exponential kernel, defaults to 50.}

\item{seed}{An integer specifying the random seed for reproducibility, defaults to 1234.}
}
\description{
Fit Model using Rcpp
}
