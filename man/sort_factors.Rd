% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/post_analysis.R
\name{sort_factors}
\alias{sort_factors}
\title{Order factors in a model based on sums of squared loadings}
\usage{
sort_factors(modgirt_rvar)
}
\arguments{
\item{modgirt_rvar}{A \code{draws_rvar} object from a MODGIRT model}
}
\value{
A \code{draws_rvar} object with factors ordered by explanatory power
}
\description{
This function takes a model based on posterior draws and orders the factors
based on their estimated sums of squares. Factors with larger sums of squares
will be placed first in the sort model.
}
