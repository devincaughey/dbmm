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
    draws_of_lambda <- lambda_rvar %>%
        lapply(draws_of, with_chains = TRUE) %>%
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
    dimnames(x$kappa_trichot) <- list(
        item = attr(x, "trichotomous_item_labels"),
        threshold = 1:dim(x$kappa_trichot)[2]
    )
    dimnames(x$kappa_ordinal) <- list(
        item = attr(x, "ordinal_item_labels"),
        threshold = 1:dim(x$kappa_ordinal)[2]
    )
    dimnames(x$sigma_metric) <- list(
        item = attr(x, "metric_item_labels")
    )
    dimnames(x$Omega) <- list(
        factor = 1:n_factor,
        factor = 1:n_factor
    )
    dimnames(x$eta) <- list(
        period = attr(x, "time_labels"),
        unit = attr(x, "unit_labels"),
        factor = 1:n_factor
    )
    dimnames(x$alpha_binary) <- list(
        period = attr(x, "time_labels"),
        item = attr(x, "binary_item_labels")
    )
    dimnames(x$alpha_trichot) <- list(
        period = attr(x, "time_labels"),
        item = attr(x, "trichotomous_item_labels")
    )
    dimnames(x$alpha_ordinal) <- list(
        period = attr(x, "time_labels"),
        item = attr(x, "ordinal_item_labels")
    )
    dimnames(x$alpha_metric) <- list(
        period = attr(x, "time_labels"),
        item = attr(x, "metric_item_labels")
    )
    dimnames(x$lambda_binary) <- list(
        item = attr(x, "binary_item_labels"),
        factor = 1:n_factor
    )
    dimnames(x$lambda_trichot) <- list(
        item = attr(x, "trichotomous_item_labels"),
        factor = 1:n_factor
    )
    dimnames(x$lambda_ordinal) <- list(
        item = attr(x, "ordinal_item_labels"),
        factor = 1:n_factor
    )
    dimnames(x$lambda_metric) <- list(
        item = attr(x, "metric_item_labels"),
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
            mean = ~mean(.),
            median = ~median(.),
            sd = ~sd(.),
            mad = ~mad(.),
            q5 = ~as.numeric(quantile(., probs = .05)),
            q95 = ~as.numeric(quantile(., probs = .95)),
            rhat = ~posterior::rhat(.),
            ess_bulk = ~posterior::ess_bulk(.),
            ess_tail = ~posterior::ess_tail(.)
        )
    }
    sfun <- function (y) {
        y %>%
            mutate(
                across(value, summary_functions),
            ) %>%
            select(-value) %>%
            rename_with(~str_remove(., "value_")) %>%
            as_tibble()
    }
    is_rvar <- sapply(x, inherits, "rvar")
    out_rvar <- out_df <- NULL
    if (any(is_rvar)) {
        out_rvar <- list()
        for (i in 1:sum(is_rvar)) {
            out_rvar[[i]] <- x[is_rvar][[i]] %>%
                summarise_draws() %>%
                mutate(
                    variable = str_replace(
                        variable,
                        "(^x\\[is_rvar\\]\\[\\[i\\]\\])|(^\\.)",
                        names(x[is_rvar])[i]
                    ))
            names(out_rvar)[i] <- names(x[is_rvar])[i]
        }
    }
    if (any(!is_rvar)) {
        out_df <- lapply(x[!is_rvar], sfun)
    }
    out <- c(out_rvar, out_df)
    return(out)
}
