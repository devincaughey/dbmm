% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/post_analysis.R
\name{identify_draws}
\alias{identify_draws}
\title{Identify the sign and rotation of the parameter draws}
\usage{
identify_draws(
  raw_draws,
  rotate = NULL,
  varimax = TRUE,
  targets = NULL,
  normalize = TRUE,
  item_type = NULL,
  sign = NULL,
  check = TRUE
)
}
\arguments{
\item{raw_draws}{(\code{dbmm_draws}) Posterior draws}

\item{rotate}{(logical) Should the factor draws be rotated? If \code{NULL} (the
default), \code{rotate} will be set to \code{TRUE} if and only if the number of
factors is greater than 1.}

\item{varimax}{(logical) Should a varimax rotation be applied within each
draw? Defaults to \code{TRUE}.}

\item{targets}{(matrix) Matrix of factor loadings to use a targets for target
rotation (using \code{GPArotation::targetT()}). Can include \code{NA} elements, in
which case partial rotation will be applied. Rows correspond to items and
columns to factors/dimensions. Must have the same number of rows as the
number of \code{item_type} items.}

\item{normalize}{(logical) Should Kaiser normalization be performed before
varimax rotation? Defaults to \code{TRUE}.}

\item{item_type}{(string) Should "binary", "ordinal", or "metric" loadings be
used to identify the model. If \code{NULL} (the default), the largest set of
items will be chosen.}

\item{sign}{(integer) Should the sign of the average identified loading be
negative (\code{-1}) or positive (\code{+1}, the default).}

\item{check}{(logical) Should the class of \code{dbmm_draws} be checked?  Defaults
to \code{TRUE}.}
}
\value{
A \code{dbmm_identified} object. Identified draws from posterior draws.
}
\description{
Identify the sign and rotation of the parameter draws
}
