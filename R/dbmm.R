
#' dbmm main function
#'
#' Fit a dbmm model using stan
#'
#' @param data (list) Data formatted for Stan, typically created by `shape()`.
#' @param chains (positive integer) The number of Markov chains to run. The
#'     default is 4.
#' @param parallelize_within_chains (logical) Should computations in a given
#'     Markov chain be parallelized using Stan's `reduce_sum()` function?
#'     Defaults to `FALSE`. If `TRUE`, then `threads_per_chain` should be set as
#'     well. Since speed gains/losses from parallelization are highly
#'     problem-specific, it is advisable to first experiment with test runs of a
#'     few iterations.
#' @param threads_per_chain (positive integer) Number of parallel processes to
#'     use for within-chain parallelization. Ignored if
#'     `parallelize_within_chains = FALSE`.
#' @param force_recompile (logical) Should `cmdstanr::compile()` be required to
#'     recompile the Stan model? Defaults to `FALSE`.
#' @param init_kappa (logical) Should initial values for the item thresholds
#'     `kappa` be pre-generated? This may reduce warnings of starting values
#'     being rejected due to incompatibility with the order constraints on
#'     `kappa`. Ignored unless `init = NULL`.
#' @param init (multiple options) The initialization method to use for the
#'     variables declared in the parameters block of the Stan program. For for
#'     details, see the documentation for `cmdstanr::sample()`.
#' @param return_data (logical) Should the data list fed to
#'     `cmdstanr::sample()`, which includes additional elements not included in
#'     the output of `shape()`, be returned with the fitted model? Defaults to
#'     `TRUE`.
#' @param n_dim (positive integer) Number of latent factors (i.e., dimensions)
#'     in the factor model. Defaults to `1`.
#' @param constant_alpha (logical) Should the item thresholds (`kappa`) be held
#'     constant across time periods? If `FALSE` (the default), the thresholds
#'     for a given item will be allowed to shift by a constant amount (governed
#'     by `alpha`) in each time period.
#' @param separate_eta (logical) Should units' factor scores (`eta`) be assigned
#'     the same priors in every time period? Defaults to `TRUE`. If `FALSE`, the
#'     scores in each period will be given priors centered on their value in the
#'     previous period, thus smoothing the estimates across periods.
#' @param lambda_zeros (multiple options) Should some item loadings (`lambda`)
#'     be fixed at 0, and if so which ones? Rotational invariance (label
#'     switching) across latent factors can be avoided by setting, for each
#'     factor $d$, $d-1$ loadings to 0. If `lambda_zeros = NULL` (the default),
#'     factor-specific parameters (e.g., `eta`) will not be identified, and the
#'     draws will have to be rotated after sampling using `identify_draws(.,
#'     rotate = TRUE)`. If `lambda_zeros = TRUE`, rotational identification will
#'     be imposed by automatically choosing $d-1$ loadings to set to 0. Users
#'     can also choose which loadings to restrict by inputting a two-column
#'     character matrix, each row of which corresponds to a restriction. The
#'     first column should be the name of an item, and the second column should
#'     be a dimension number (as a character, e.g., `"2"`).
#' @param df_sigma_metric (positive real) Degrees of freedom of the Student's t
#'     prior for `sigma_metric`, the residual standard deviations of the metric
#'     items. Defaults to 4.
#' @param df_sigma_alpha_evol (positive real) Degrees of freedom of the
#'     Student's t prior for `sigma_alpha_evol`, the standard deviation of the
#'     dynamic prior for `alpha`. Defaults to 4.
#' @param df_sigma_eta_evol (positive real) Degrees of freedom of the (half)
#'     Student's t prior for `sigma_eta_evol`, the standard deviation of the
#'     dynamic prior for `eta`. Defaults to `4`.
#' @param mu_sigma_metric (real) Mean of the (half) Student's t prior for
#'     `sigma_alpha_evol`, the residual standard deviations of the metric
#'     items. Defaults to `0.5`.
#' @param mu_sigma_alpha_evol (real) Mean of the (half) Student's t prior for
#'     `sigma_alpha_evol`, the standard deviation of the dynamic prior for
#'     `alpha`. Defaults to `0.1`.
#' @param mu_sigma_eta_evol (real) Mean of the (half) Student's t prior for
#'     `sigma_eta_evol`, the standard deviation of the dynamic prior for
#'     `eta`. Defaults to `0.1`.
#' @param sd_sigma_metric (positive real) Standard deviation of the (half)
#'     Student's t prior for `sigma_alpha_evol`, the residual standard
#'     deviations of the metric items. Defaults to `0.5`.
#' @param sd_sigma_alpha_evol (positive real) Standard deviation of the
#'     Student's t prior for `sigma_alpha_evol`, the standard deviation of the
#'     dynamic prior for `alpha`. Defaults to `0.1`.
#' @param sd_sigma_eta_evol (positive real) Standard deviation of the (half)
#'     Student's t prior for `sigma_alpha_evol`, the residual standard
#'     deviations of the metric items. Defaults to `0.1`.
#' @param seed (positive integer) An integer vector of length one indicating the
#'     state of Stan’s pseudo-random number generator. Defaults to `123`.
#' @param link (string) Which link function should be used for binary and
#'     ordinal outcomes. One of `"probit"` (the default) and `"logit"`.
#' @param ... Additional arguments to `cmdstanr::sample()`.
#'
#' @return A dbmm object containing
#'  \describe{
#'    \item{unit_labels}{}
#'    \item{time_labels}{}
#'    \item{binary_item_labels}{}
#'    \item{ordinal_item_labels}{}
#'    \item{metric_item_labels}{}
#' }
#'
#' @import cmdstanr
#'
#' @export
fit <- function (data,
              # removed `file` argument for now
                chains = 4,
                parallelize_within_chains = FALSE,
                threads_per_chain = NULL,
                force_recompile = FALSE,
                init_kappa = FALSE,
                init = NULL,
                return_data = TRUE,
                n_dim = 1,
                constant_alpha = FALSE,
                separate_eta = TRUE,
                lambda_zeros = NULL,
                df_sigma_metric = 4,
                df_sigma_alpha_evol = 4,
                df_sigma_eta_evol = 4,
                mu_sigma_metric = 0.5,
                mu_sigma_alpha_evol = 0.5,
                mu_sigma_eta_evol = 0.5,
                sd_sigma_metric = 0.5,
                sd_sigma_alpha_evol = 0.5,
                sd_sigma_eta_evol = 0.5,
                seed = 123,
                link = "probit",
                ...)
{

  check_arg_type(arg = data, typename = "dbmm_data")

  if (parallelize_within_chains & !is.null(threads_per_chain)) {
    check_arg_type(arg = threads_per_chain, typename = "numeric")
    check_arg_type(arg = chains, typename = "numeric")

    specified_cores <- threads_per_chain * chains
    if (parallel::detectCores() <= specified_cores) {
      cli::cli_alert_warning("The number of specified cores exceeds the number of cores in your computer.")
    }
  }

  stopifnot(!parallelize_within_chains || threads_per_chain > 0)

  ## Add model options to input data
  data$parallelize <- as.integer(parallelize_within_chains)
  data$constant_alpha <- as.integer(constant_alpha)
  data$separate_eta <- as.integer(separate_eta)
  data$D <- n_dim
  data$df_sigma_metric <- df_sigma_metric
  data$df_sigma_alpha_evol <- df_sigma_alpha_evol
  data$df_sigma_eta_evol <- df_sigma_eta_evol
  data$mu_sigma_metric <- mu_sigma_metric
  data$mu_sigma_alpha_evol <- mu_sigma_alpha_evol
  data$mu_sigma_eta_evol <- mu_sigma_eta_evol
  data$sd_sigma_metric <- sd_sigma_metric
  data$sd_sigma_alpha_evol <- sd_sigma_alpha_evol
  data$sd_sigma_eta_evol <- sd_sigma_eta_evol

  nonzero_binary <- matrix(
    data = 1,
    nrow = data$I_binary,
    ncol = data$D,
    dimnames = list(attr(data, "binary_item_labels"), as.character(1:data$D))
  )

  nonzero_ordinal <- matrix(
    data = 1,
    nrow = data$I_ordinal,
    ncol = data$D,
    dimnames = list(attr(data, "ordinal_item_labels"), as.character(1:data$D))
  )

  nonzero_metric <- matrix(
    data = 1,
    nrow = data$I_metric,
    ncol = data$D,
    dimnames = list(attr(data, "metric_item_labels"), as.character(1:data$D))
  )

  if (data$D > 1 && isTRUE(lambda_zeros)) {
    most_items <- which.max(c(data$I_binary, data$I_ordinal, data$I_metric))
    id_with <- c("binary", "ordinal", "metric")[most_items]
    for (d in 2:data$D) {
      if (id_with == "binary") nonzero_binary[1:(d-1), d] <- 0
      if (id_with == "ordinal") nonzero_ordinal[1:(d-1), d] <- 0
      if (id_with == "metric") nonzero_metric[1:(d-1), d] <- 0
    }
  }

  if (data$D > 1 && !isTRUE(lambda_zeros) && !is.null(lambda_zeros)) {
    for (i in 1:nrow(lambda_zeros)) {
      if (lambda_zeros[i, 1] %in% attr(data, "binary_item_labels")) {
        nonzero_binary[lambda_zeros[i, 1], lambda_zeros[i, 2]] <- 0
      }
      if (lambda_zeros[i, 1] %in% attr(data, "ordinal_item_labels")) {
        nonzero_ordinal[lambda_zeros[i, 1], lambda_zeros[i, 2]] <- 0
      }
      if (lambda_zeros[i, 1] %in% attr(data, "metric_item_labels")) {
        nonzero_metric[lambda_zeros[i, 1], lambda_zeros[i, 2]] <- 0
      }
    }
  }
  data$nonzero_binary <- nonzero_binary
  data$nonzero_ordinal <- nonzero_ordinal
  data$nonzero_metric <- nonzero_metric

  ## Compile model
  file <- system.file(paste0("stan/", link, ".stan"), package = "dbmm")

  m0 <- cmdstan_model(stan_file = file)

  cpp_opts <- list(stan_threads = as.logical(data$parallelize))

  m1 <- m0$compile(cpp_options = cpp_opts, force_recompile = force_recompile)

  if (init_kappa & !is.null(init)) {
    cat("\nIgnoring `init_kappa = TRUE` because `init` is not `NULL`.\n")
  }

  ## Initial Values
  if (init_kappa & is.null(init)) {
    kappa_inits <- replicate(chains, {
      kappa <- stats::rnorm(data$I_ordinal * (data$K - 1))
      kappa_mat <- apply(matrix(kappa, ncol = data$K - 1), 1, sort)
      list(kappa = kappa_mat)
    },
    simplify = FALSE
    )
    init <- kappa_inits
  }

  ## Fit model
  dynfac_fit <- m1$sample(data, chains = chains, init = init,
                          threads_per_chain = threads_per_chain, seed = seed,
                          ...)

  ## Prepare output
  attr(dynfac_fit, "unit_labels") <- attr(data, "unit_labels")
  attr(dynfac_fit, "time_labels") <- attr(data, "time_labels")
  attr(dynfac_fit, "binary_item_labels") <- attr(data, "binary_item_labels")
  attr(dynfac_fit, "ordinal_item_labels") <- attr(data, "ordinal_item_labels")
  attr(dynfac_fit, "metric_item_labels") <- attr(data, "metric_item_labels")
  out <- list(fit = dynfac_fit)

  if (return_data) {
    out$data <- data
  }

  class(out) <- c("dbmm_fitted", class(out))

  return(out)
}
