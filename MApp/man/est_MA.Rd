% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/individ_ic_fun.R
\name{est_MA}
\alias{est_MA}
\title{Compute model averaged estimates for partial regression
  coefficients.}
\usage{
est_MA(weights, coef_mat, inmat, w_plus = FALSE)
}
\arguments{
\item{weights}{A vector of approximate posterior model probabilities
(commonly referred to as model weights in the IC literature).}

\item{coef_mat}{A \eqn{J} by \eqn{p} matrix of maximum likelihood
estimates for the partial regression coefficients from all individual
models considered.}

\item{inmat}{A matrix defining which variables are in each model.}

\item{w_plus}{A logical statement defining which type of MA to use. 
Default is \code{FALSE}, which allows all models defined in \code{inmat} to 
be used in the averaging for each variable. If set to \code{TRUE}, MA
estimates will be made using only the models for which the coefficient
associated with each variable is not set to exactly zero. 
Weights are normalized conditional on the models meeting this criterion.}
}
\value{
The point estimate of the model averaged partial regression
  coefficients.
}
\description{
This function computes point estimates for model averaged
  partial regression coefficients using a weighted averaged of maximum
  likelihood estimates from individual models according to approximate
  posterior model probabilities.
}

