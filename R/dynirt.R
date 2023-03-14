#' dynIRT main function
#'
#' Fit a dynIRT model using stan
#'
#' @param data
#' @param chains
#' @param parallelize_within_chains
#' @param threads_per_chain
#' @param force_recompile
#' @param init_kappa
#' @param init
#' @param return_data
#' @param n_dim
#' @param constant_alpha
#' @param separate_eta
#' @param lambda_zeros
#' @param df_sigma_metric
#' @param df_sigma_alpha_evol
#' @param df_sigma_eta_evol
#' @param mu_sigma_metric
#' @param mu_sigma_alpha_evol
#' @param mu_sigma_eta_evol
#' @param sd_sigma_metric
#' @param sd_sigma_alpha_evol
#' @param sd_sigma_eta_evol
#' @param seed An integer vector of length one indicating the state of Stan’s pseudo-random number generator
#' @param ...
#'
#' @return A dynIRT object containing
#'  \describe{
#'    \item{unit_labels}
#'    \item{time_labels}
#'    \item{binary_item_labels}
#'    \item{ordinal_item_labels}
#'    \item{metric_item_labels}
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
                mu_sigma_alpha_evol = 0.1,
                mu_sigma_eta_evol = 0.1,
                sd_sigma_metric = 0.5,
                sd_sigma_alpha_evol = 0.1,
                sd_sigma_eta_evol = 0.1,
                seed = 123,
                ...)
{

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
    dimnames = list(attr(data, "binary_item_labels"), NULL)
  )

  nonzero_ordinal <- matrix(
    data = 1,
    nrow = data$I_ordinal,
    ncol = data$D,
    dimnames = list(attr(data, "ordinal_item_labels"), NULL)
  )

  nonzero_metric <- matrix(
    data = 1,
    nrow = data$I_metric,
    ncol = data$D,
    dimnames = list(attr(data, "metric_item_labels"), NULL)
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
    for (i in seq_along(lambda_zeros)) {
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
  file <- system.file(paste0("stan/model.stan"), package = "dynIRTtest")

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

  class(out) <- c("dynIRT_fitted", class(out))

  return(out)
}
