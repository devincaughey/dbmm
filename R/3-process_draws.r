#' Extract draws from fitted model
#'
#' @param x (`dbmm_fitted` object) A fitted model produced by `fit()`.
#' @param drop (character vector) A regular expression (or logical scalar).
#'     Parameters that match any of the regular expressions will be dropped.
#' @param format (string) The format of the returned draws or point
#'     estimates. Must be a valid format from the ‘posterior’ package. The
#'     default is `"df"`, which is what other `dbmm` functions will
#'     expect, but other options include `"array"`, `"matrix"`, and `"list"`.
#' @param check (logical)
#'
#' @return Draws from the posterior distribution of the selected parameters.
#'
#' @import magrittr
#' @import cmdstanr
#'
#' @export
extract_mixfac_draws <- function(x, drop = "^z_|^chol", check = TRUE)
{

    if (check) {
        check_arg_type(arg = x, typename = "mixfac_fit")
    }

    draws <- as_draws_rvars(x$fit$draws())

    if (isTRUE(drop)) {
        draws <- posterior::subset_draws(
                                x = draws,
                                variable = "^z_|^chol",
                                regex = TRUE,
                                exclude = TRUE
                            )
    } else if (is.character(drop)) {
        draws <- posterior::subset_draws(
                                x = draws,
                                variable = drop,
                                regex = TRUE,
                                exclude = TRUE
                            )
    }

    attr(draws, "unit_labels") <- attr(x$fit, "unit_labels")
    attr(draws, "time_labels") <- attr(x$fit, "time_labels")
    attr(draws, "binary_item_labels") <- attr(x$fit, "binary_item_labels")
    attr(draws, "trichotomous_item_labels") <-
        attr(x$fit, "trichotomous_item_labels")
    attr(draws, "ordinal_item_labels") <- attr(x$fit, "ordinal_item_labels")
    attr(draws, "metric_item_labels") <- attr(x$fit, "metric_item_labels")

    class(draws) <- c("mixfac_draws", class(draws))

    return(draws)
}


#' Exported function
#' @export
identify_mixfac <- function(x, method = "varimax", identify_with_type) {
    ## Accept either a draws_rvars object or a fitted object / draws_df
    if (posterior::is_draws_rvars(x)) {
        draws_rvar <- x
    } else {
        draws_rvar <- posterior::as_draws_rvars(x$fit$draws())
    }
    ## Subset to lambda draws (rvar)
    lambda_rvar <-
        posterior::subset_draws(draws_rvar, variable = "^lambda", regex = TRUE)
    draws_of_lambda <- lambda_rvar |>
        lapply(draws_of, with_chains = TRUE) |>
        abind::abind(along = 3)
    ## Dimensions (chains / iterations / factors)
    n_chain <- posterior::nchains(draws_rvar)
    n_iter <- posterior::niterations(draws_rvar)
    n_factor <- dim(draws_of_lambda)[4]
    ## Choose which loading variable to use based on item_type
    if (missing(identify_with_type)) {
        identify_with_type <-
            names(lambda_rvar)[which.max(sapply(lambda_rvar, length))][1]
        identify_with_type <- sub("^lambda_", "", identify_with_type)
    }
    varname <- switch(
        identify_with_type,
        binary = "lambda_binary",
        trichot = "lambda_trichot",
        ordinal = "lambda_ordinal",
        metric = "lambda_metric",
        stop("Invalid `identify_with_type` argument; must be 'binary', 'trichot', 'ordinal' or 'metric'")
    )
    ## Make varimax matrices
    if (n_factor > 1) {
        vm_rvar <- make_vm_rvar(
            draws_of_lambda,
            n_iter = n_iter,
            n_chain = n_chain,
            n_factor = n_factor,
            method = method
        )
    } else {
        vm_rvar <- matrix(1)
    }
    ## Apply varimax rotations to lambdas
    for (t in seq_along(lambda_rvar)) {
        lambda_rvar[[t]] <- posterior::`%**%`(lambda_rvar[[t]], vm_rvar)        
    }
    ## Compute signed-permutation matrices (factor-switching alignment)
    ## Convert rotated lambdas to a draws-matrix in the format expected by :
    lambda_matrix <- posterior::as_draws_matrix(t(lambda_rvar[[varname]]))
    lambda_matrix <- rename_loading_matrix(lambda_matrix)
    ## factor.switching::rsp_exact expects rows ordered by (chain then iteration),
    ## and returns a list with sign_vectors and permute_vectors
    rsp_out <- factor.switching::rsp_exact(lambda_matrix, rotate = FALSE)
    ## Create sp_rvar (signed-permutation rvar) and apply to lambdas
    sp_rvar <- make_sp_rvar(rsp_out, n_iter, n_chain, n_factor)
    ## Apply signed permutations to lambdas
    for (t in seq_along(lambda_rvar)) {
        lambda_rvar[[t]] <- posterior::`%**%`(lambda_rvar[[t]], sp_rvar)        
    }
    ## RSP matrices
    vm_sp_rvar <- posterior::`%**%`(vm_rvar, sp_rvar)
    ## Apply rotations to `eta`
    eta_rvar <- posterior::subset_draws(draws_rvar, variable = "eta")
    for (t in seq_len(dim(eta_rvar$eta)[1])) {
        eta_rvar$eta[t, , ] <- posterior::`%**%`(
            as.matrix(eta_rvar$eta[t, , , drop = TRUE]),
            vm_sp_rvar
        )
    }
    omega_rvar <- posterior::subset_draws(draws_rvar, variable = "Omega")
    omega_rvar$Omega <- t(vm_sp_rvar) %**% omega_rvar$Omega %**% vm_sp_rvar
    out <- posterior::draws_rvars(
                          eta = eta_rvar$eta,
                          lambda_binary = lambda_rvar$lambda_binary,
                          lambda_trichot = lambda_rvar$lambda_trichot,
                          lambda_ordinal = lambda_rvar$lambda_ordinal,
                          lambda_metric = lambda_rvar$lambda_metric,
                          alpha_binary = draws_rvar$alpha_binary,
                          alpha_trichot = draws_rvar$alpha_trichot,
                          alpha_ordinal = draws_rvar$alpha_ordinal,
                          alpha_metric = draws_rvar$alpha_metric,
                          kappa_trichot = draws_rvar$kappa_trichot,
                          kappa_ordinal = draws_rvar$kappa_ordinal,
                          sigma_alpha_evol = draws_rvar$sigma_alpha_evol,
                          sigma_metric = draws_rvar$sigma_metric,
                          Omega = omega_rvar$Omega,
                          lp__ = draws_rvar$lp__
                      )
    attr(out, "unit_labels") <- attr(x, "unit_labels")
    attr(out, "time_labels") <- attr(x, "time_labels")
    attr(out, "binary_item_labels") <- attr(x, "binary_item_labels")
    attr(out, "trichotomous_item_labels") <- attr(x, "trichotomous_item_labels")
    attr(out, "ordinal_item_labels") <- attr(x, "ordinal_item_labels")
    attr(out, "metric_item_labels") <- attr(x, "metric_item_labels")
    attr(out, "rotation matrix") <- vm_rvar
    attr(out, "signed-permutation matrix") <- sp_rvar

    class(out) <- c("mixfac_id", class(out))

    return(out)
}


#' Exported function
#' @export
label_mixfac <- function (x, make_long = TRUE, check = TRUE) {
    stopifnot(is_draws_rvars(x))
    n_factor <- dim(x$eta)[[3]]
    if (length(attr(x, "binary_item_labels")) > 0) {
        dimnames(x$alpha_binary) <- list(
            period = attr(x, "time_labels"),
            item = attr(x, "binary_item_labels")
        )
        dimnames(x$lambda_binary) <- list(
            item = attr(x, "binary_item_labels"),
            factor = 1:n_factor
        )
    }
    if (length(attr(x, "trichotomous_item_labels")) > 0) {
        dimnames(x$alpha_trichot) <- list(
            period = attr(x, "time_labels"),
            item = attr(x, "trichotomous_item_labels")
        )
        dimnames(x$kappa_trichot) <- list(
            item = attr(x, "trichotomous_item_labels"),
            threshold = 1:dim(x$kappa_trichot)[2]
        )
        dimnames(x$lambda_trichot) <- list(
            item = attr(x, "trichotomous_item_labels"),
            factor = 1:n_factor
        )
    }
    if (length(attr(x, "ordinal_item_labels")) > 0) {
        dimnames(x$alpha_ordinal) <- list(
            period = attr(x, "time_labels"),
            item = attr(x, "ordinal_item_labels")
        )
        dimnames(x$kappa_ordinal) <- list(
            item = attr(x, "ordinal_item_labels"),
            threshold = 1:dim(x$kappa_ordinal)[2]
        )
        dimnames(x$lambda_ordinal) <- list(
            item = attr(x, "ordinal_item_labels"),
            factor = 1:n_factor
        )
    }
    if (length(attr(x, "metric_item_labels")) > 0) {
        dimnames(x$alpha_metric) <- list(
            period = attr(x, "time_labels"),
            item = attr(x, "metric_item_labels")
        )
        dimnames(x$lambda_metric) <- list(
            item = attr(x, "metric_item_labels"),
            factor = 1:n_factor
        )
        dimnames(x$sigma_metric) <- list(
            item = attr(x, "metric_item_labels")
        )
    }
    dimnames(x$Omega) <- list(
        factor = 1:n_factor,
        factor = 1:n_factor
    )
    dimnames(x$eta) <- list(
        period = attr(x, "time_labels"),
        unit = attr(x, "unit_labels"),
        factor = 1:n_factor
    )
    if (make_long) {
        for (i in seq_along(x)) {
            if (!str_detect(names(x)[i], "^lp|^Omega|^sigma_alpha_evol")) {
                x[[i]] <- as.data.frame.table(x[[i]], responseName = "value")
            }
        }
    }
    if (make_long) {
        class(x) <- c("mixfac_lab", "list")
    } else {
        class(x) <- c("mixfac_lab", class(x))
    }
    return(x)
}


#' Exported function
#' @export
combine_types_mixfac <- function (x) {
    ## lambda
    lambda_idx <- grep("^lambda_", names(x))
    lambda <- dplyr::bind_rows(x[lambda_idx], .id = "item_type")
    lambda$item_type <- sub("^lambda_", "", lambda$item_type)
    alpha_idx <- grep("^alpha_", names(x))
    alpha <- dplyr::bind_rows(x[alpha_idx], .id = "item_type")
    alpha$item_type <- sub("^alpha_", "", alpha$item_type)
    kappa_idx <- grep("^kappa_", names(x))
    kappa <- dplyr::bind_rows(x[kappa_idx], .id = "item_type")
    kappa$item_type <- sub("^kappa_", "", kappa$item_type)
    out <- c(
        list(lambda = lambda, alpha = alpha, kappa = kappa),
        x[-c(lambda_idx, alpha_idx, kappa_idx)]
    )
    class(out) <- c("mixfac_comb", class(out))
    return(out)
}


#' Exported function
#' @export
summarize_mixfac <- function (x, summary_functions) {
    if (missing(summary_functions)) {
        summary_functions <- list(
            mean = ~posterior::E(.),
            median = ~posterior:::median.rvar(.),
            sd = ~posterior::sd(.),
            mad = ~posterior::mad(.),
            q5 = ~as.numeric(posterior:::quantile.rvar(., probs = .05)),
            q95 = ~as.numeric(posterior:::quantile.rvar(., probs = .95)),
            rhat = ~posterior::rhat(.),
            ess_bulk = ~posterior::ess_bulk(.),
            ess_tail = ~posterior::ess_tail(.)
        )
    }
    sfun <- function (y) {
        y |>
            mutate(
                across(value, summary_functions),
            ) |>
            select(-value) |>
            rename_with(~str_remove(., "value_")) |>
            as_tibble()
    }
    is_rvar <- sapply(x, inherits, "rvar")
    out_rvar <- out_df <- NULL
    if (any(is_rvar)) {
        out_rvar <- list()
        for (i in 1:sum(is_rvar)) {
            if (length(x[is_rvar][[i]]) == 0) {
                out_rvar[[i]] <- NA
            } else {
                out_rvar[[i]] <-
                    x[is_rvar][[i]] |>
                    summarise_draws() |>
                    mutate(
                        variable = str_replace(
                            variable,
                            "(^x\\[is_rvar\\]\\[\\[i\\]\\])|(^\\.)",
                            names(x[is_rvar])[i]
                        ))
            }
            names(out_rvar)[i] <- names(x[is_rvar])[i]
        }
    }
    if (any(!is_rvar)) {
        out_df <- lapply(x[!is_rvar], sfun)
    }
    out <- c(out_rvar, out_df)
    return(out)
}

#' order factors in a model based on sums of squared loadings
#'
#' This function takes a model based on posterior draws and orders the factors
#' based on their estimated sums of squares. Factors with larger sums of squares
#' will be placed first in the sort model.
#'
#' @param mixfac_rvar A `draws_rvar` object from a mixed-factor model
#'
#' @return A `draws_rvar` object with factors ordered by explanatory power
#'
#' @export
sort_mixfac <- function(x) {
    check_arg_type(arg = x, typename = "mixfac_comb")
    ss <- posterior::rvar_apply(mixfac_rvar$lambda, 2, function(x) {
        posterior::rvar_sum(x^2)
    })
    fo <- order(-posterior::E(ss))
    out <-
        posterior::draws_rvars(
                       eta = x$eta[, , fo],
                       lambda = x$lambda[, fo],
                       alpha = x$alpha,
                       kappa = x$kappa,
                       sigma_alpha_evol = x$sigma_alpha_evol,
                       sigma_metric = x$sigma_metric,
                       Omega = x$Omega[fo, fo],
                       lp__ = x$lp__
                   )
    return(out)
}

#' Set Signs
#'
#' This function sets the signs of the parameters of a MIXFAC model based on 
#' user-defined signs.
#'
#' @param x The model object containing the parameters.
#' @param signs A vector of signs to be applied to the parameters. 
#' Scalar values are allowed and will be recycled. Default is 1.
#'
#' @return A modified model object with the signs of the parameters updated.
#'
#' @details This function sets the signs of the parameters in the model object
#' \code{x} based on the user-defined signs provided in the 
#' \code{signs} argument. The function applies the sign flips to the parameters 
#' and returns a modified model object with the updated signs.
#'
#' @import posterior
#'
#' @export
sign_mixfac <- function(x, signs = 1) {
    check_arg_type(arg = x, typename = "mixfac_comb")
    n_time <- dim(x$eta)[1]
    n_factor <- dim(x$eta)[3]
    stopifnot(length(signs == 1) || length(signs) == n_factor)
    init_signs <- sign(colMeans(E(x$lambda)))
    sign_flips <- ifelse(init_signs == signs, 1, -1)
    sm <- diag(sign_flips, nrow = n_factor, ncol = n_factor)
    for (t in seq_len(n_time)) {
        x$eta[t, , drop = TRUE] <-
            x$eta[t, , drop = TRUE] %**% sm
    }
    posterior::draws_rvars(
        eta = x$eta,
        lambda = x$lambda %**% sm,
        alpha = x$alpha,
        kappa = x$kappa,
        sigma_alpha_evol = x$sigma_alpha_evol,
        sigma_metric = x$sigma_metric,
        Omega = t(sm) %**% x$Omega %**% sm,
        lp__ = x$lp__
    )
}


#' Identify MODGIRT draws
#'
#' This function identifies the MODGIRT model by postprocessing the draws from
#' the posterior distribution.
#'
#' @param x A fitted MODGIRT model object or `draws_rvars` object
#'
#' @return A list containing the identified MODGIRT model parameters.
#'
#' @import posterior
#'
#' @export
identify_modgirt <- function(x, method = "varimax") {
    ## Store draws in `rvars` object
    if (posterior::is_draws_rvars(x)) {
        modgirt_rvar <- x
    } else {
        modgirt_rvar <- posterior::as_draws_rvars(x$fit$draws())
    }
    beta_rvar <-
        posterior::subset_draws(modgirt_rvar, variable = "beta")
    draws_of_beta <- posterior::draws_of(beta_rvar$beta, with_chains = TRUE)
    bar_theta_rvar <-
        posterior::subset_draws(modgirt_rvar, variable = "bar_theta")
    n_chain <- posterior::nchains(modgirt_rvar)
    n_iter <- posterior::niterations(modgirt_rvar)
    n_factor <- ncol(beta_rvar$beta)
    ## Create draw-specific varimax rotations
    if (n_factor > 1) {
        vm_rvar <- make_vm_rvar(draws_of_beta, n_iter, n_chain, n_factor,
                                method = method)
    } else {
        vm_rvar <- matrix(1)
    }
    ## Apply varimax rotations to `beta`
    beta_rvar$beta <- posterior::`%**%`(beta_rvar$beta, vm_rvar)
    ## Create draw-specific signed permutations
    beta_matrix <- posterior::as_draws_matrix(t(beta_rvar$beta))
    lambda_matrix <- rename_loading_matrix(beta_matrix)
    rsp_out <- factor.switching::rsp_exact(lambda_matrix, rotate = FALSE)
    sp_rvar <- make_sp_rvar(rsp_out, n_iter, n_chain, n_factor)
    ## Apply signed permutations to `beta`
    beta_rvar$beta <- posterior::`%**%`(beta_rvar$beta, sp_rvar)
    ## Make single RSP matrix
    vm_sp_rvar <- posterior::`%**%`(vm_rvar, sp_rvar)
    ## Apply rotations to `bar_theta`
    for (t in seq_len(dim(bar_theta_rvar$bar_theta)[1])) {
        bar_theta_rvar$bar_theta[t, , ] <- posterior::`%**%`(
            as.matrix(bar_theta_rvar$bar_theta[t, , , drop = TRUE]),
            vm_sp_rvar
        )
    }
    sigma_theta_rvar <-
        posterior::subset_draws(modgirt_rvar, variable = "Sigma_theta")
    omega_rvar <-
        posterior::subset_draws(modgirt_rvar, variable = "Omega")
    sigma_theta_rvar$Sigma_theta <-
        t(vm_sp_rvar) %**% sigma_theta_rvar$Sigma_theta %**% vm_sp_rvar
    omega_rvar$Omega <-
        t(vm_sp_rvar) %**% omega_rvar$Omega %**% vm_sp_rvar
    modgirt_rvar_id <- posterior::draws_rvars(
        lp__ = modgirt_rvar$lp__,
        alpha = modgirt_rvar$alpha,
        beta = beta_rvar$beta,
        bar_theta = bar_theta_rvar$bar_theta,
        Sigma_theta = sigma_theta_rvar$Sigma_theta,
        Omega = omega_rvar$Omega
    )
    out_ls <- list(
        modgirt_rvar = modgirt_rvar_id,
        vm_rvar = vm_rvar,
        sp_rvar = sp_rvar
    )
    return(out_ls)
}


#' Apply rotation to MODGIRT draws
#'
#' This function applies the given rotation to each draw from the posterior
#' distribution of the MODGIRT parameters
#'
#' @param modgirt_rvar A `draws_rvar` object from a MODGIRT model
#' @param rotat An I-by-D rotation matrix
#'
#' @return A `draws_rvar` object of rotated draws
#'
#' @examples
#' rotmat <- varimax(E(modgirt_signed$beta))$rotmat
#' modgirt_rotated <- rotate_modgirt(modgirt_signed, rotmat)
#' 
#' @import posterior
#'
#' @export
rotate_modgirt <- function(modgirt_rvar, rotmat) {
    ## inverse of transpose (needed for oblique rotation)
    if (is_rvar(rotmat)) {
        rvar_solve <- rfun(solve)
        G <- rvar_solve(t(rotmat))
    } else {
        G <- solve(t(rotmat))
    }
    ## Create parameter-specific `draws_rvar` objects
    beta_rvar <- subset_draws(modgirt_rvar, variable = "beta")
    bar_theta_rvar <- subset_draws(modgirt_rvar, variable = "bar_theta")
    sigma_theta_rvar <- subset_draws(modgirt_rvar, variable = "Sigma_theta")
    omega_rvar <- subset_draws(modgirt_rvar, variable = "Omega")
    n_time <- dim(bar_theta_rvar$bar_theta)[1]
    ## Apply rotations
    beta_rvar$beta <- beta_rvar$beta %**% G
    for (t in seq_len(n_time)) {
        bar_theta_rvar$bar_theta[t, , ] <- 
            bar_theta_rvar$bar_theta[t, , , drop = TRUE] %**% G
    }
    sigma_theta_rvar$Sigma_theta <-
        t(G) %**% sigma_theta_rvar$Sigma_theta %**% G
    omega_rvar$Omega <- t(G) %**% omega_rvar$Omega %**% G
    modgirt_rvar_rot <- draws_rvars(
        lp__ = modgirt_rvar$lp__,
        alpha = modgirt_rvar$alpha,
        beta = beta_rvar$beta,
        bar_theta = bar_theta_rvar$bar_theta,
        Sigma_theta = sigma_theta_rvar$Sigma_theta,
        Omega = omega_rvar$Omega)
    return(modgirt_rvar_rot)
}

#' Order factors in a model based on sums of squared loadings
#'
#' This function takes a model based on posterior draws and orders the factors
#' based on their estimated sums of squares. Factors with larger sums of squares
#' will be placed first in the sort model.
#'
#' @param modgirt_rvar A `draws_rvar` object from a MODGIRT model
#'
#' @return A `draws_rvar` object with factors ordered by explanatory power
#'
#' @export
sort_modgirt <- function(modgirt_rvar) {
    ss <- posterior::rvar_apply(modgirt_rvar$beta, 2, function(x) {
        posterior::rvar_sum(x^2)
    })
    fo <- order(-posterior::E(ss))
    sorted_rvar <- posterior::draws_rvars(
        lp__ = modgirt_rvar$lp__,
        alpha = modgirt_rvar$alpha,
        beta = modgirt_rvar$beta[, fo],
        bar_theta = modgirt_rvar$bar_theta[, , fo],
        Sigma_theta = modgirt_rvar$Sigma_theta[fo, fo],
        Omega = modgirt_rvar$Omega[fo, fo]
    )
    return(sorted_rvar)
}

#' Set Signs
#'
#' This function sets the signs of the parameters of a MODGIRT model based on 
#' user-defined signs.
#'
#' @param modgirt_rvar The model object containing the parameters.
#' @param signs A vector of signs to be applied to the parameters. 
#' Scalar values are allowed and will be recycled. Default is 1.
#'
#' @return A modified model object with the signs of the parameters updated.
#'
#' @details This function sets the signs of the parameters in the model object
#' \code{modgirt_rvar} based on the user-defined signs provided in the 
#' \code{signs} argument. The function applies the sign flips to the parameters 
#' and returns a modified model object with the updated signs.
#'
#' @import posterior
#'
#' @export
sign_modgirt <- function(modgirt_rvar, signs = 1) {
    n_time <- dim(modgirt_rvar$bar_theta)[1]
    n_factor <- dim(modgirt_rvar$bar_theta)[3]
    stopifnot(length(signs == 1) || length(signs) == n_factor)
    init_signs <- sign(colMeans(E(modgirt_rvar$beta)))
    sign_flips <- ifelse(init_signs == signs, 1, -1)
    sm <- diag(sign_flips, nrow = n_factor, ncol = n_factor)
    for (t in seq_len(n_time)) {
        modgirt_rvar$bar_theta[t, , drop = TRUE] <-
            modgirt_rvar$bar_theta[t, , drop = TRUE] %**% sm
    }
    posterior::draws_rvars(
        lp__ = modgirt_rvar$lp__,
        alpha = modgirt_rvar$alpha,
        beta = modgirt_rvar$beta %**% sm,
        bar_theta = modgirt_rvar$bar_theta,
        Sigma_theta = t(sm) %**% modgirt_rvar$Sigma_theta %**% sm,
        Omega = t(sm) %**% modgirt_rvar$Omega %**% sm
    )
}

