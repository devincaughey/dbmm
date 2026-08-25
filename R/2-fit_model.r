
#' Fit a multidimensional ordinal group-level IRT model using Stan
#'
#' This function fits a MODGIRT model using Stan.
#' It takes in the following arguments:
#'
#' @param stan_data The data for the model.
#' @param chains The number of chains to use (default is 4).
#' @param return_data A logical value indicating whether to return the input
#' data along with the model fit (default is TRUE).
#' @param n_factor The number of factors in the model.
#' @param force_recompile A logical value indicating whether to force
#' recompilation of the Stan model (default is FALSE).
#' @param signed_loadings The signed loadings matrix.
#' @param nonzero_loadings The nonzero loadings matrix.
#' @param link The link function to use (default is "probit").
#' @param seed The random seed for reproducibility.
#' @param shape_lkj_theta (positive real) Shape of the LKJ prior on the
#'     within-group correlation of `theta`. `1` is uniform over correlation
#'     matrices; larger values shrink toward zero correlation. Defaults to
#'     `2`. See Details.
#' @param shape_lkj_bar_theta_evol (positive real) Shape of the LKJ prior on
#'     the correlation of the random-walk innovations in `bar_theta`.
#'     Defaults to `2`. See Details.
#' @param df_sd_theta (positive real) Degrees of freedom of the half
#'     Student's t prior on `sd_theta`, the within-group standard deviations
#'     of `theta` in each latent dimension. Defaults to `4`. Set to `1` to
#'     recover the half Cauchy prior used in earlier versions; larger values
#'     shrink the tail and, with it, the sampler's excursions into extreme
#'     within-group dispersions that the data cannot rule out.
#' @param scale_sd_theta (positive real) Scale of the half Student's t prior
#'     on `sd_theta`. Defaults to `1`. Because period-1 group means are
#'     whitened, the between-group standard deviation is fixed at 1 in every
#'     dimension, so this default places the within-group dispersion of
#'     `theta` on the same scale as the dispersion between groups.
#' @param df_sd_bar_theta_evol (positive real) Degrees of freedom of the half
#'     Student's t prior on `sd_bar_theta_evol`, the standard deviations of
#'     the random-walk innovations in `bar_theta`. Defaults to `4`. Ignored
#'     when `T == 1`, since there is then no transition to inform it.
#' @param scale_sd_bar_theta_evol (positive real) Scale of the half Student's
#'     t prior on `sd_bar_theta_evol`. Defaults to `0.1`, which shrinks the
#'     group means toward being constant over time: on the whitened scale of
#'     period 1, an innovation standard deviation of `0.1` moves a group by
#'     roughly a tenth of the between-group spread from one period to the
#'     next. Ignored when `T == 1`.
#' @param ... Additional arguments to be passed to the Stan sampling function.
#'
#' @details
#' # Regularizing the latent covariances
#'
#' `Sigma_theta`, the within-group covariance of `theta`, enters the
#' likelihood only through the scalars \eqn{\beta_q' \Sigma_\theta \beta_q},
#' one per item. Its off-diagonal elements are therefore informed only by
#' items that load on more than one dimension, and are weakly identified in
#' most applications. Some shrinkage toward zero correlation is usually
#' advisable, which is what `shape_lkj_theta` controls.
#'
#' Two subtleties are worth knowing before raising it well above its default.
#'
#' First, zero correlation is a property of the basis, not of the covariance
#' operator. A matrix that is diagonal in one orthogonal basis is not
#' diagonal in another unless its variances happen to be equal. Shrinking
#' correlations toward zero therefore expresses a mild preference for the
#' basis in which `Sigma_theta` is diagonal — that is, for its eigenbasis,
#' or equivalently the principal components of the within-group dispersion
#' of `theta`. At the default value this preference is negligible relative
#' to the likelihood; it becomes appreciable only for shape parameters an
#' order of magnitude larger.
#'
#' Second, that preferred basis is not the one you end up reporting.
#' [identify_modgirt()] rotates the draws to a varimax solution defined by
#' the sparsity of the loadings, which is a different criterion. The prior
#' consequently affects where the sampler explores, and hence how well it
#' mixes, more than it affects the identified estimates. Because
#' identification maps \eqn{\Sigma_\theta \to R' \Sigma_\theta R}, the
#' reported `Sigma_theta` should not be expected to look near-diagonal even
#' under strong shrinkage.
#'
#' If you want regularization that is genuinely indifferent to rotation, the
#' target is sphericity (\eqn{\sigma^2 I}) rather than diagonality: shrink
#' the elements of `sd_theta` toward a common value, via `df_sd_theta` and
#' `scale_sd_theta`, as well as shrinking the correlations.
#' 
#' Both scale defaults are interpretable because the period-1 group means are
#' de-meaned and whitened, which fixes the between-group covariance of
#' `bar_theta` to the identity. Latent distances are therefore measured in
#' units of the period-1 between-group standard deviation, and priors on
#' `sd_theta` and `sd_bar_theta_evol` can be set on that common scale.
#' 
#' @return A list containing the model fit and optionally the input data.
#'
#' @export
fit_modgirt <- function(
    stan_data,
    chains = 4,
    return_data = TRUE,
    n_factor = 1,
    force_recompile = FALSE,
    signed_loadings,
    nonzero_loadings,
    link = "probit",
    seed = NULL,
    shape_lkj_theta = 2,
    shape_lkj_bar_theta_evol = 2,
    df_sd_theta = 4,
    scale_sd_theta = 1,
    df_sd_bar_theta_evol = 4,
    scale_sd_bar_theta_evol = 0.1,
    ...
) {
    n_item <- stan_data$Q
    if (missing(signed_loadings)) {
        signed_loadings <- matrix(
            data = 0,
            nrow = n_item,
            ncol = n_factor,
            dimnames = list(dimnames(stan_data$SSSS)$ITEM, seq_len(n_factor))
        )
    }
    if (missing(nonzero_loadings)) {
        nonzero_loadings <- matrix(
            data = 1,
            nrow = n_item,
            ncol = n_factor,
            dimnames = list(dimnames(stan_data$SSSS)$ITEM, seq_len(n_factor))
        )
    }
    stopifnot(isTRUE(n_factor == ncol(signed_loadings)))
    stopifnot(isTRUE(n_factor == ncol(nonzero_loadings)))
    stan_data$D <- n_factor
    stan_data$beta_nonzero <- nonzero_loadings
    stan_data$beta_sign <- signed_loadings
    stan_data$shape_lkj_theta <- shape_lkj_theta
    stan_data$shape_lkj_bar_theta_evol <- shape_lkj_bar_theta_evol
    stan_data$df_sd_theta <- df_sd_theta
    stan_data$scale_sd_theta <- scale_sd_theta
    stan_data$df_sd_bar_theta_evol <- df_sd_bar_theta_evol
    stan_data$scale_sd_bar_theta_evol <- scale_sd_bar_theta_evol
    file <- system.file(
        paste0("stan/modgirt_", link, ".stan"),
        package = "dbmm"
    )
    if (!nzchar(file)) {
        stop("No Stan file found for `link = \"", link, "\"`.")
    }    
    m0 <- cmdstan_model(stan_file = file)
    m1 <- m0$compile(force_recompile = force_recompile)
    modgirt_fit <- m1$sample(stan_data, chains = chains, seed = seed, ...)
    ## Prepare output
    attr(modgirt_fit, "unit_labels") <- attr(stan_data, "unit_labels")
    attr(modgirt_fit, "time_labels") <- attr(stan_data, "time_labels")
    attr(modgirt_fit, "item_labels") <- attr(stan_data, "item_labels")
    out <- list(fit = modgirt_fit)
    if (return_data) {
        out$stan_data <- stan_data
    }
    return(out)
}

make_mixfac_out <- function(fit, shaped_data, return_data = TRUE) {
    if (missing(fit)) stop("`fit` is required.")
    if (missing(shaped_data)) stop("`shaped_data` is required.")
    fit <- copy_dbmm_attrs(fit, shaped_data)
    out <- list(fit = fit)
    if (isTRUE(return_data)) {
        out$shaped_data <- shaped_data
    }
    class(out) <- c("mixfac_fit", class(out))
    out
}


#' Fit a dynamic mixed factor model using Stan.
#'
#' @param shaped_data (list) Data formatted for Stan, typically created by `shape()`.
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
#' @param smooth_eta (logical) Should units' factor scores (`eta`) be smoothed
#'     across time periods? Defaults to `TRUE`, in which case the scores in
#'     each period are given priors centred on their values in the previous
#'     period, with evolution standard deviation `sigma_eta_evol`. If `FALSE`,
#'     each period's scores are estimated independently, no random walk is
#'     fitted, and `sigma_eta_evol`, `Lcorr_eta`, and a free `Omega` are not
#'     estimated.
#' @param separate_eta (logical) Deprecated. Use `smooth_eta` instead, noting
#'     that the sense is reversed: `separate_eta = TRUE` corresponds to
#'     `smooth_eta = FALSE`.
#' @param whiten_eta (logical) Should units' factor scores (`eta`) be whitened
#'     for identification?  Defaults to `TRUE`.
#' @param gen_log_lik (logical) Should the per-observation log likelihood be
#'     computed and returned, for use with the \pkg{loo} package? Defaults to
#'     `FALSE`, since `log_lik` has one element per observation per draw and
#'     can dominate the size of the output. See `log_lik_index()` for the
#'     mapping from positions of `log_lik` to observations.
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
#'     `sigma_metric`, the residual standard deviations of the metric
#'     items. Defaults to `0.5`.
#' @param sd_sigma_metric (positive real) Standard deviation of the (half)
#'     Student's t prior for `sigma_metric`, the residual standard
#'     deviations of the metric items. Defaults to `0.5`.
#' @param mu_sigma_alpha_evol (non-negative real) Location of the half
#'     Student's t prior for `sigma_alpha_evol`, the standard deviation of the
#'     dynamic prior for `alpha`. Defaults to `0`, which shrinks the intercept
#'     evolution toward zero (i.e., toward time-constant intercepts).
#' @param mu_sigma_eta_evol (non-negative real) Location of the half Student's
#'     t prior for `sigma_eta_evol`, the standard deviation of the dynamic
#'     prior for `eta`. Defaults to `0`, which shrinks the factor evolution
#'     toward zero (i.e., toward time-constant factor scores).
#' @param sd_sigma_alpha_evol (positive real) Scale of the half Student's t
#'     prior for `sigma_alpha_evol`. Defaults to `0.1`.
#' @param sd_sigma_eta_evol (positive real) Scale of the half Student's t prior
#'     for `sigma_eta_evol`. Defaults to `0.1`.
#' @param seed (positive integer) An integer vector of length one indicating the
#'     state of Stan’s pseudo-random number generator. Defaults to `123`.
#' @param link (string) Which link function should be used for binary and
#'     ordinal outcomes. Currently only `"probit"` is supported.
#' @param ... Additional arguments to `cmdstanr::sample()`.
#'
#' @return A dbmm object containing
#'  \describe{
#'    \item{unit_labels}{}
#'    \item{time_labels}{}
#'    \item{binary_item_labels}{}
#'    \item{trichotomous_item_labels}{}
#'    \item{ordinal_item_labels}{}
#'    \item{metric_item_labels}{}
#' }
#'
#' @import cmdstanr
#'
#' @export
fit_mixfac <- function(shaped_data,
                       chains = 4,
                       parallelize_within_chains = FALSE,
                       threads_per_chain = NULL,
                       force_recompile = FALSE,
                       init = NULL,
                       return_data = TRUE,
                       n_dim = 1,
                       constant_alpha = FALSE,
                       smooth_eta = TRUE,
                       whiten_eta = TRUE,
                       gen_log_lik = FALSE,
                       lambda_zeros = NULL,
                       df_sigma_metric = 4,
                       df_sigma_alpha_evol = 4,
                       df_sigma_eta_evol = 4,
                       mu_sigma_metric = 0.5,
                       mu_sigma_alpha_evol = 0,  
                       mu_sigma_eta_evol = 0,    
                       sd_sigma_metric = 0.5,
                       sd_sigma_alpha_evol = 0.1,
                       sd_sigma_eta_evol = 0.1,
                       seed = NULL,
                       link = "probit",
                       separate_eta = NULL,
                       ...) {
    check_arg_type(arg = shaped_data, typename = "mixfac_data")
    ## Deprecated `separate_eta`, whose sense was the reverse of `smooth_eta`
    if (!is.null(separate_eta)) {
        if (!missing(smooth_eta)) {
            cli::cli_abort(paste(
                     "Supply either {.arg smooth_eta} or the deprecated",
                     "{.arg separate_eta}, not both."
                 ))
        }
        cli::cli_warn(paste(
                 "{.arg separate_eta} is deprecated; use {.arg smooth_eta} instead.",
                 "The sense is reversed: {.code separate_eta = TRUE} corresponds to",
                 "{.code smooth_eta = FALSE}."
             ))
        smooth_eta <- !check_flag(separate_eta)
    }

    smooth_eta               <- check_flag(smooth_eta)
    constant_alpha           <- check_flag(constant_alpha)
    whiten_eta               <- check_flag(whiten_eta)
    gen_log_lik              <- check_flag(gen_log_lik)
    parallelize_within_chains <- check_flag(parallelize_within_chains)

    check_arg_type(arg = chains, typename = "numeric")

    if (parallelize_within_chains) {
        if (is.null(threads_per_chain) || !is.numeric(threads_per_chain) ||
            length(threads_per_chain) != 1L || is.na(threads_per_chain) ||
            threads_per_chain < 1) {
            cli::cli_abort(
                "{.arg threads_per_chain} must be a positive integer when
                 {.code parallelize_within_chains = TRUE}."
            )
        }
        avail <- parallel::detectCores()
        if (!is.na(avail) && threads_per_chain * chains > avail) {
            cli::cli_warn(c(
                "Requested {threads_per_chain * chains} core{?s} but only
                 {avail} detected.",
                "i" = "Oversubscription usually slows sampling."
            ))
        }
    } else if (!is.null(threads_per_chain)) {
        cli::cli_warn(
            "{.arg threads_per_chain} is ignored when
             {.code parallelize_within_chains = FALSE}."
        )
        threads_per_chain <- NULL
    }

    ## Add model options to input shaped_data
    shaped_data$parallelize <- as.integer(parallelize_within_chains)

    ## Add model options to input shaped_data
    shaped_data$parallelize <- as.integer(parallelize_within_chains)
    shaped_data$constant_alpha <- as.integer(constant_alpha)
    shaped_data$smooth_eta <- as.integer(smooth_eta)
    shaped_data$whiten_eta <- as.integer(whiten_eta)
    shaped_data$gen_log_lik <- as.integer(gen_log_lik)
    shaped_data$D <- n_dim
    shaped_data$df_sigma_metric <- df_sigma_metric
    shaped_data$df_sigma_alpha_evol <- df_sigma_alpha_evol
    shaped_data$df_sigma_eta_evol <- df_sigma_eta_evol
    shaped_data$mu_sigma_metric <- mu_sigma_metric
    shaped_data$mu_sigma_alpha_evol <- mu_sigma_alpha_evol
    shaped_data$mu_sigma_eta_evol <- mu_sigma_eta_evol
    shaped_data$sd_sigma_metric <- sd_sigma_metric
    shaped_data$sd_sigma_alpha_evol <- sd_sigma_alpha_evol
    shaped_data$sd_sigma_eta_evol <- sd_sigma_eta_evol

    nonzero_binary <- matrix(
        data = 1,
        nrow = shaped_data$I_binary,
        ncol = shaped_data$D,
        dimnames = list(
            attr(shaped_data, "binary_item_labels"),
            as.character(1:shaped_data$D)
        )
    )

    nonzero_trichot <- matrix(
        data = 1,
        nrow = shaped_data$I_trichot,
        ncol = shaped_data$D,
        dimnames = list(
            attr(shaped_data, "trichotomous_item_labels"),
            as.character(1:shaped_data$D)
        )
    )

    nonzero_ordinal <- matrix(
        data = 1,
        nrow = shaped_data$I_ordinal,
        ncol = shaped_data$D,
        dimnames = list(
            attr(shaped_data, "ordinal_item_labels"),
            as.character(1:shaped_data$D)
        )
    )

    nonzero_metric <- matrix(
        data = 1,
        nrow = shaped_data$I_metric,
        ncol = shaped_data$D,
        dimnames = list(
            attr(shaped_data, "metric_item_labels"),
            as.character(1:shaped_data$D)
        )
    )

    if (shaped_data$D > 1 && isTRUE(lambda_zeros)) {
        item_ns <- c(
            shaped_data$I_binary,
            shaped_data$I_trichot,
            shaped_data$I_ordinal,
            shaped_data$I_metric
        )
        most_items <- which.max(item_ns)
        id_with <- c("binary", "trichotomous", "ordinal", "metric")[most_items]
        for (d in 2:shaped_data$D) {
            if (id_with == "binary") nonzero_binary[1:(d-1), d] <- 0
            if (id_with == "trichotomous") nonzero_trichot[1:(d-1), d] <- 0
            if (id_with == "ordinal") nonzero_ordinal[1:(d-1), d] <- 0
            if (id_with == "metric") nonzero_metric[1:(d-1), d] <- 0
        }
    }

    if (shaped_data$D > 1 && !isTRUE(lambda_zeros) && !is.null(lambda_zeros)) {
        for (i in seq_len(nrow(lambda_zeros))) {
            if (lambda_zeros[i, 1] %in% attr(shaped_data, "binary_item_labels")) {
                nonzero_binary[lambda_zeros[i, 1], lambda_zeros[i, 2]] <- 0
            }
            if (lambda_zeros[i, 1] %in% attr(shaped_data, "trichotomous_item_labels")) {
                nonzero_trichot[lambda_zeros[i, 1], lambda_zeros[i, 2]] <- 0
            }
            if (lambda_zeros[i, 1] %in% attr(shaped_data, "ordinal_item_labels")) {
                nonzero_ordinal[lambda_zeros[i, 1], lambda_zeros[i, 2]] <- 0
            }
            if (lambda_zeros[i, 1] %in% attr(shaped_data, "metric_item_labels")) {
                nonzero_metric[lambda_zeros[i, 1], lambda_zeros[i, 2]] <- 0
            }
        }
    }
    shaped_data$nonzero_binary <- nonzero_binary
    shaped_data$nonzero_trichot <- nonzero_trichot
    shaped_data$nonzero_ordinal <- nonzero_ordinal
    shaped_data$nonzero_metric <- nonzero_metric

    ## Compile model
    if (!identical(link, "probit")) {
        stop("Only `link = \"probit\"` is currently supported for mixfac models.")
    }
    file <- system.file(paste0("stan/mixfac_", link, ".stan"), package = "dbmm")
    if (!nzchar(file)) {
        stop("No Stan file found for `link = \"", link, "\"`.")
    }
    m0 <- cmdstan_model(stan_file = file)

    cpp_opts <- list(stan_threads = as.logical(shaped_data$parallelize))

    m1 <- m0$compile(cpp_options = cpp_opts, force_recompile = force_recompile)

    ## Fit model
    fit <- m1$sample(shaped_data, chains = chains, init = init,
                            threads_per_chain = threads_per_chain, seed = seed,
                            ...)

    ## Prepare output
    out <- make_mixfac_out(
        fit = fit,
        shaped_data = shaped_data,
        return_data = return_data
    )
    return(out)
}
