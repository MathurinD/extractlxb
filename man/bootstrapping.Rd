% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_lxb_data.R
\name{bootstrapping}
\alias{bootstrapping}
\title{Bootstrap function}
\usage{
bootstrapping(values, func = median, error = 0.05, min_run = 100,
  max_run = 1000 * 1000, EPS = 0.001)
}
\arguments{
\item{func}{Function to use to aggregate the data}

\item{error}{Error rate of the confidence interval to report}

\item{min_run}{Minimum number of bootstrap steps to take}

\item{max_run}{Maximum number of bootstrap steps to take}

\item{ESP}{Convergence limit}
}
\value{
a vector or length 3 with the mean, and the 'error' confidence interval of this mean. 5\% by default (0.025 percentile on each side)
}
\description{
Converging boostrap, with a minimum number of run to ensure some initial variation
}
