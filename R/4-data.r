#' Policies of the U.S. states, 2010-2012
#'
#' Long-format panel of 111 state policies observed annually from 2010
#' through 2012. Most policies are binary indicators of whether a state had
#' the policy in force; some are ordinal, and some are continuous (tax
#' rates, per-pupil expenditures, benefit levels).
#'
#' @format A data frame with 16,500 rows and 8 columns:
#' \describe{
#'   \item{year}{Year, 2010-2012.}
#'   \item{state_abb}{Two-letter postal abbreviation.}
#'   \item{state_name}{State name.}
#'   \item{policy_variable}{Policy identifier.}
#'   \item{policy_short_description}{Brief description; \code{NA} for one item.}
#'   \item{policy_longer_description}{Fuller description.}
#'   \item{value}{Response, in nominal terms where monetary.}
#'   \item{value_real}{Response, in constant dollars where monetary.}
#' }
#'
#' @source Devin Caughey and Christopher Warshaw, <dynamicdemocracy.us>. Accessed 24 August 2026.
"state_policies_2010_2012"
