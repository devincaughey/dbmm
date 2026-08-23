#' State-level social outcomes, 2020-2021
#'
#' A panel of 30 social, economic, health, and environmental indicators for the
#' 50 U.S. states in 2020 and 2021, used to illustrate the dynamic mixed factor
#' model in `vignette("dbmm")`. The data are in the long format expected by
#' [shape_mixfac()], with one row per state-year-outcome.
#'
#' All indicators are continuous and are treated as metric items. Their native
#' scales differ substantially, from rates per 100,000 to proportions, so
#' standardization before fitting is advisable.
#'
#' The panel is complete in structure (50 states x 2 years x 30 outcomes =
#' 3,000 rows), but 611 values are missing, in three distinct patterns. Four
#' indicators (`x_firearm_ownership_rate`, `x_marijuana_rate`,
#' `x_smoking_packs_pc`, and `x_union_membership`) have no observed values at
#' all, and are dropped by [shape_mixfac()], leaving 26 items. Four more
#' (`x_abortion_rate_guttmacher_occurence`, `x_index_crime_rate`,
#' `x_state_co2`, and `x_teen_birth_rate`) are missing throughout 2021, having
#' not yet been published when the data were assembled. The remainder are
#' scattered single states. The second pattern illustrates one motivation for
#' the dynamic model: the evolution prior on the factor scores lets sparsely
#' observed periods borrow strength from adjacent ones.
#'
#' @format A [tibble][tibble::tibble] with 3,000 rows and 4 variables:
#' \describe{
#'   \item{st}{(character) Two-letter state postal abbreviation.}
#'   \item{year}{(integer) Calendar year, `2020` or `2021`.}
#'   \item{outcome}{(character) Indicator name. Names encode the source where
#'     applicable, e.g. `cps` for the Current Population Survey and
#'     `guttmacher` for the Guttmacher Institute.}
#'   \item{value}{(numeric) Observed value on the indicator's native scale, or
#'     `NA` if unavailable.}
#' }
#'
#' @source TODO: add citations, or point to the sources listed in
#'     `vignette("dbmm")`.
#'
#' @examples
#' head(social_outcomes_2020_2021)
#'
#' # Indicators with no observed values
#' with(social_outcomes_2020_2021, names(which(tapply(value, outcome,
#'      function(v) all(is.na(v))))))
#'
#' # Items with no observed values are dropped, leaving 26 metric items
#' shaped <- shape_mixfac(
#'     social_outcomes_2020_2021,
#'     unit_var = "st",
#'     time_var = "year",
#'     item_var = "outcome",
#'     value_var = "value"
#' )
#' shaped$I_metric
"social_outcomes_2020_2021"
