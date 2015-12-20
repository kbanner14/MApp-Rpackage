% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/plotting_fun_bms.R
\name{MApp_bms}
\alias{MApp_bms}
\title{Model Averaged Posteriors Plot}
\usage{
MApp_bms(x, plot_wind, num_sims = 1000, max_display = NULL,
  mod_names = NULL, include_coef = NULL, ...)
}
\arguments{
\item{x}{an object of class \code{bma} from the \code{BMS::bms()} function}

\item{plot_wind}{A vector of length 2 specifying the number of rows and
columns to partition the plotting window into.}

\item{max_display}{An integer specifying the number of models to display in
the plotting window. Can also take on the argument "common3" to show the
top PMP model, the full model, and the model averaged distribuiton.}

\item{mod_names}{A vector specifying the names of the models. If left blank,
models will be named \eqn{M_1, M_2,..., M_K}, where \eqn{M_1} will have the
largest posterior model probability and \eqn{M_K} has the smallests
posterior model probability.}

\item{include_coef}{A vector specifying which partial regression coefficients
to display.}
}
\value{
The \code{MApp_bms} plot and a table of posterior standard deviations
  for the top \code{max_display} individual models specified, along with the
  posterior standard deviation of the model averaged parameter.
}
\description{
For \code{bma} objects obtained using the \code{bms} function in
  the \code{BMS} package. Creates a graphical display showing all relevant
  output from the model averaging procedure. Used to graphically compare
  posterior distributions of partial regression coefficients from individual
  models in the model set to the posterior distributions resulting from the
  model averaged partial regression coefficients. Also used to compare
  posterior standard deviations of the continuous piece of the model averaged
  distribution to corresponding standard deviations from individual posterior
  distributions in tabular form.
}
