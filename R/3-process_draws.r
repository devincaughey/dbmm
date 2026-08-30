.mixfac_drop_default <- "^z_|^chol|^r_|^L_|^Lcorr|^WW$|^sigma_eta_evol"

#' Extract draws from fitted model
#'
#' @param x (`dbmm_fitted` object) A fitted model produced by `fit()`.
#' @param drop (character vector) A regular expression (or logical scalar).
#'     Parameters that match any of the regular expressions will be dropped.
#' @param check (logical)
#'
#' @return Draws from the posterior distribution of the selected parameters.
#'
#' @import cmdstanr
#' @importFrom posterior %**%
#'
#' @export
extract_mixfac_draws <- function(x, drop = .mixfac_drop_default, check = TRUE) {

    if (check) {
        check_arg_type(arg = x, typename = "mixfac_fit")
    }

    draws <- as_draws_rvars(x$fit$draws())

    if (isTRUE(drop)) {
        draws <- posterior::subset_draws(
                                x = draws,
                                variable = .mixfac_drop_default,
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

    draws <- copy_dbmm_attrs(draws, x$fit)
    class(draws) <- c("mixfac_draws", class(draws))

    return(draws)
}


#' Identify mixfac model
#'
#' @param x Either a draws_rvars object or a fitted object
#' @param method Rotation criterion
#' @param whiten Should `bar_theta` be demeaned, standardized, and made orthogonal?
#' @param ref_t In what time period should `bar_theta` be whitened? May be
#'     an integer in 1:T, 'last', or 'mean'. The latter means that the
#'     within-unit averages `bar_theta` are whitened.
#' @param identify_with_type (string) Which loading matrix to use when aligning
#'     factors across draws: `"lambda"`, `"binary"`, `"trichot"`, `"ordinal"`,
#'     or `"metric"`. Defaults to the item type with the most items. Forced to
#'     `"lambda"` when `x` inherits from `"mixfac_comb"`.
#' @param random_starts (non-negative integer) Number of random starting
#'     rotations tried per draw. Defaults to `0`, which starts from the
#'     identity matrix and is deterministic. Values greater than `0` guard
#'     against local minima of the rotation criterion, at proportionate
#'     computational cost.
#' @param seed (positive integer or `NULL`) Seed for the random starts.
#'     Defaults to `123`, making results reproducible. Ignored when
#'     `random_starts = 0`, since the rotation is then deterministic. Set to
#'     `NULL` to use the ambient random number stream.
#' @param quiet (logical) Suppress the iteration log printed by
#'     factor.switching::rsp_exact() Defaults to `TRUE`.
#' @export
identify_mixfac <- function(x, method = "varimax", whiten = FALSE,
                            ref_t = "last", identify_with_type,
                            random_starts = 0, seed = 123,
                            quiet = TRUE) {
    ## Accept either a draws_rvars object or a fitted object. Label attributes
    ## live on the cmdstanr fit rather than on the wrapper list, so the source
    ## of the draws and the source of the labels differ in the latter case.
    if (posterior::is_draws_rvars(x)) {
        draws_rvar <- x
        label_src <- x
    } else if (!is.null(x$fit)) {
        draws_rvar <- posterior::as_draws_rvars(x$fit$draws())
        label_src <- x$fit
    } else {
        stop("`x` must be a draws_rvars object or a fitted model object ",
             "containing a `fit` element.")
    }

    ## Subset to lambda draws (rvar)
    lambda_rvar <-
        posterior::subset_draws(draws_rvar, variable = "^lambda", regex = TRUE)

    ## Dimensions
    n_chain <- posterior::nchains(draws_rvar)
    n_iter <- posterior::niterations(draws_rvar)
    n_time <- dim(draws_rvar$eta)[1]
    n_unit <- dim(draws_rvar$eta)[2]
    n_factor <- dim(draws_rvar$eta)[3]

    if (identical(ref_t, "last")) ref_t <- n_time

    ## Choose which loading variable to use for rsp alignment
    if (inherits(x, "mixfac_comb")) {
        identify_with_type <- "lambda"
    } else if (missing(identify_with_type)) {
        identify_with_type <-
            names(lambda_rvar)[which.max(sapply(lambda_rvar, length))][1]
        identify_with_type <- sub("^lambda_", "", identify_with_type)
    }
    varname <- switch(
        identify_with_type,
        lambda = "lambda",
        binary = "lambda_binary",
        trichot = "lambda_trichot",
        ordinal = "lambda_ordinal",
        metric = "lambda_metric",
        stop("Invalid `identify_with_type`; must be 'lambda', 'binary', 'trichot', 'ordinal', or 'metric'.")
    )

    has_omega <- "Omega" %in% posterior::variables(draws_rvar)
    omega_rvar <- if (has_omega) {
                      posterior::subset_draws(draws_rvar, variable = "Omega")
                  } else NULL
    
    ## ---------- Optional whitening ----------
    ## eta_norm   = (eta - m) %*% WW
    ## lambda_norm = lambda %*% inv(WW)'
    ## alpha_norm  = alpha + lambda %*% m        (original basis)
    ## Omega_norm  = WW' %*% Omega %*% WW
    ## kappa is unchanged. Together these leave every linear predictor
    ## nu = alpha + lambda %*% eta exactly invariant. The varimax and
    ## signed-permutation rotations applied afterwards are orthogonal, so
    ## they preserve lambda %*% eta and require no further adjustment.
    if (whiten) {
        eta_raw <- draws_rvar$eta

        ## Anchor matrix XX (J x D)
        if (identical(ref_t, "mean")) {
            XX <- posterior::rvar_apply(eta_raw, c(2, 3), posterior::rvar_mean)
        } else {
            if (!is.numeric(ref_t) || length(ref_t) != 1L ||
                ref_t < 1 || ref_t > n_time) {
                stop("`ref_t` must be 'last', 'mean', or an integer in 1:T.")
            }
            XX <- eta_raw[ref_t, , , drop = TRUE]
            if (n_factor == 1) XX <- t(t(XX))
        }

        ## Anchor column means, and demeaned anchor
        m_rvar <- posterior::rvar_apply(XX, 2, posterior::rvar_mean)
        DM <- XX
        for (d in seq_len(n_factor)) {
            DM[, d] <- XX[, d] - m_rvar[d]
        }

        ## Whitening matrix, and its inverse-transpose for the loadings
        ## (WW is triangular, not orthogonal, so inv(WW)' != WW)
        WW_rvar <- rvar_whiten_matrix(DM)
        Ginv_rvar <- t(rvar_solve(WW_rvar))

        ## --- Transform eta ---
        eta_w <- eta_raw
        for (tt in seq_len(n_time)) {
            E_t <- eta_raw[tt, , , drop = TRUE]
            if (n_factor == 1) E_t <- t(t(E_t))
            for (d in seq_len(n_factor)) {
                E_t[, d] <- E_t[, d] - m_rvar[d]
            }
            eta_w[tt, , ] <- posterior::`%**%`(E_t, WW_rvar)
        }
        draws_rvar$eta <- eta_w

        ## --- Shift intercepts, using the ORIGINAL loadings and anchor mean ---
        ## Done on the draws arrays rather than with rvar slice assignment,
        ## which recycles unpredictably across the draws dimension.
        shift_alpha_exact <- function(alpha_mat, lambda_mat, m_vec) {
            if (is.null(alpha_mat) || is.null(lambda_mat)) return(alpha_mat)
            A <- posterior::draws_of(alpha_mat)    # [draws, T, I]
            L <- posterior::draws_of(lambda_mat)   # [draws, I, D]
            M <- posterior::draws_of(m_vec)        # [draws, D]
            if (is.null(dim(M))) M <- matrix(M, ncol = 1L)
            n_draw <- dim(A)[1]
            n_it <- dim(A)[3]
            ## S[s, i] = sum_d L[s, i, d] * M[s, d]
            S <- matrix(0, nrow = n_draw, ncol = n_it)
            for (d in seq_len(dim(L)[3])) {
                Ld <- array(L[, , d], dim = dim(L)[1:2])
                S <- S + Ld * M[, d]
            }
            for (tt in seq_len(dim(A)[2])) {
                A[, tt, ] <- A[, tt, ] + S
            }
            out <- posterior::rvar(A, nchains = posterior::nchains(alpha_mat))
            dimnames(out) <- dimnames(alpha_mat)
            out
        }

        if (inherits(x, "mixfac_comb")) {
            if (!is.null(draws_rvar$alpha)) {
                draws_rvar$alpha <- shift_alpha_exact(
                    draws_rvar$alpha, lambda_rvar$lambda, m_rvar
                )
            }
        } else {
            for (ty in c("binary", "trichot", "ordinal", "metric")) {
                anm <- paste0("alpha_", ty)
                lnm <- paste0("lambda_", ty)
                if (!is.null(draws_rvar[[anm]]) && !is.null(lambda_rvar[[lnm]])) {
                    draws_rvar[[anm]] <- shift_alpha_exact(
                        draws_rvar[[anm]], lambda_rvar[[lnm]], m_rvar
                    )
                }
            }
        }

        ## --- Transform loadings (must come AFTER the intercept shift) ---
        for (k in seq_along(lambda_rvar)) {
            lambda_rvar[[k]] <- posterior::`%**%`(lambda_rvar[[k]], Ginv_rvar)
        }

        ## --- Transform the innovation covariance ---
        if (has_omega) {
            omega_rvar$Omega <-
                t(WW_rvar) %**% omega_rvar$Omega %**% WW_rvar
        }
    }

    ## Recompute draws_of_lambda AFTER optional whitening
    draws_of_lambda <- lambda_rvar |>
        lapply(posterior::draws_of, with_chains = TRUE) |>
        abind::abind(along = 3)

    ## ---------- Varimax ----------
    if (n_factor > 1) {
        vm_rvar <- with_seed(seed, make_vm_rvar(
                                       draws_of_lambda,
                                       n_iter = n_iter, n_chain = n_chain, n_factor = n_factor,
                                       method = method, randomStarts = random_starts
                                   ))
    } else {
        vm_rvar <- matrix(1)
    }

    for (k in seq_along(lambda_rvar)) {
        lambda_rvar[[k]] <- posterior::`%**%`(lambda_rvar[[k]], vm_rvar)
    }

    ## ---------- Signed permutation alignment ----------
    lambda_without_names <- lambda_rvar[[varname]]
    dimnames(lambda_without_names) <- NULL
    lambda_matrix <- posterior::as_draws_matrix(t(lambda_without_names))
    lambda_matrix <- rename_loading_matrix(lambda_matrix)
    if (quiet) {
        invisible(
            utils::capture.output(
                       rsp_out <- factor.switching::rsp_exact(
                                                        lambda_matrix,
                                                        rotate = FALSE
                                                    )
                   ))
    } else {
        rsp_out <- factor.switching::rsp_exact(lambda_matrix, rotate = FALSE)
    }
    sp_rvar <- make_sp_rvar(rsp_out, n_iter, n_chain, n_factor)

    for (k in seq_along(lambda_rvar)) {
        lambda_rvar[[k]] <- posterior::`%**%`(lambda_rvar[[k]], sp_rvar)
        dimnames(lambda_rvar[[k]]) <- list(
            item = dimnames(lambda_rvar[[k]])[[1]],
            factor = seq_len(dim(lambda_rvar[[k]])[2])
        )
    }

    vm_sp_rvar <- posterior::`%**%`(vm_rvar, sp_rvar)

    ## Rotate eta
    eta_rvar <- posterior::subset_draws(draws_rvar, variable = "eta")
    eta_rvar$eta <- rotate_eta_rvar(eta_rvar$eta, vm_sp_rvar)

    ## Rotate Omega
    if (has_omega) {
        omega_rvar$Omega <-
            t(vm_sp_rvar) %**% omega_rvar$Omega %**% vm_sp_rvar
    }

    ## ---------- Output ----------
    if (inherits(x, "mixfac_comb")) {
        out <- as_draws_rvars_safe(
            eta = eta_rvar$eta,
            lambda = lambda_rvar$lambda,
            alpha = draws_rvar$alpha,
            kappa_trichot = draws_rvar$kappa_trichot,
            kappa_ordinal = draws_rvar$kappa_ordinal,
            sigma_alpha_evol = draws_rvar$sigma_alpha_evol,
            sigma_metric = draws_rvar$sigma_metric,
            Omega = if (has_omega) omega_rvar$Omega else NULL,
            lp__ = draws_rvar$lp__
        )
    } else {
        out <- as_draws_rvars_safe(
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
            Omega = if (has_omega) omega_rvar$Omega else NULL,
            lp__ = draws_rvar$lp__
        )
    }

    out <- copy_dbmm_attrs(out, label_src)
    attr(out, "rotation matrix") <- vm_rvar
    attr(out, "signed-permutation matrix") <- sp_rvar

    class(out) <- c("mixfac_id", class(out))
    if (inherits(x, "mixfac_comb")) {
        class(out) <- c("mixfac_comb", class(out))
    }

    return(out)
}


#' Exported function
#' @param x An `draws_rvars` object, such as the output of `identify_mixfac()`.
#' @param make_long Should the output be a long data frame? (Requires dropping variables) 
#' @export
label_mixfac <- function (x, make_long = FALSE) {
    stopifnot(is_draws_rvars(x))
    n_factor <- dim(x$eta)[[3]]
    if (length(attr(x, "binary_item_labels")) > 0) {
        if (any(grepl("alpha_binary", names(x)))) {
            dimnames(x$alpha_binary) <- list(
                period = attr(x, "time_labels"),
                item = attr(x, "binary_item_labels")
            )
        }
        dimnames(x$lambda_binary) <- list(
            item = attr(x, "binary_item_labels"),
            factor = 1:n_factor
        )
    }
    if (length(attr(x, "trichotomous_item_labels")) > 0) {
        if (any(grepl("alpha_trichot", names(x)))) {
            dimnames(x$alpha_trichot) <- list(
                period = attr(x, "time_labels"),
                item = attr(x, "trichotomous_item_labels")
            )
        }
        if (any(grepl("kappa_trichot", names(x)))) {
            dimnames(x$kappa_trichot) <- list(
                item = attr(x, "trichotomous_item_labels"),
                threshold = 1:dim(x$kappa_trichot)[2]
            )
        }
        dimnames(x$lambda_trichot) <- list(
            item = attr(x, "trichotomous_item_labels"),
            factor = 1:n_factor
        )
    }
    if (length(attr(x, "ordinal_item_labels")) > 0) {
        if (any(grepl("alpha_ordinal", names(x)))) {
            dimnames(x$alpha_ordinal) <- list(
                period = attr(x, "time_labels"),
                item = attr(x, "ordinal_item_labels")
            )
        }
        if (any(grepl("kappa_ordinal", names(x)))) {
            dimnames(x$kappa_ordinal) <- list(
                item = attr(x, "ordinal_item_labels"),
                threshold = 1:dim(x$kappa_ordinal)[2]
            )
        }
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
    if (!is.null(x$Omega)) {
        dimnames(x$Omega) <- list(
            factor = 1:n_factor,
            factor = 1:n_factor
        )
    }
    dimnames(x$eta) <- list(
        period = attr(x, "time_labels"),
        unit = attr(x, "unit_labels"),
        factor = 1:n_factor
    )
    if (make_long) {
        for (i in seq_along(x)) {
            if (!stringr::str_detect(names(x)[i], "^lp|^Omega|^sigma_alpha_evol")) {
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


#' Combine item-type-specific parameters into single variables
#'
#' Stacks the type-specific loading and intercept variables (`lambda_binary`,
#' `lambda_trichot`, ...) into single `lambda` and `alpha` variables, prefixing
#' item names with their type. Accepts either a `draws_rvars` object or the
#' long-format list produced by `label_mixfac(make_long = TRUE)`.
#'
#' @param x A `mixfac_id` or `mixfac_lab` object. Item labels are normally
#'     assigned by [label_mixfac()] before calling this function; if they are
#'     absent, item indices are used instead.
#'
#' @return A `mixfac_comb` object of the same format as `x`.
#'
#' @export
combine_types_mixfac <- function (x) {
    ## Decide the format once. label_mixfac(make_long = TRUE) yields data
    ## frames; otherwise elements are rvars. Testing per-variable with
    ## all(sapply(...)) is unreliable: a type absent from the fit gives an
    ## empty selection, for which all() is vacuously TRUE.
    is_long <- any(vapply(x, is.data.frame, logical(1)))

    lambda_idx <- grep("^lambda_", names(x))
    alpha_idx  <- grep("^alpha_",  names(x))
    kappa_idx  <- grep("^kappa_",  names(x))

    if (length(lambda_idx) == 0) {
        stop("No `lambda_*` variables found in `x`. ",
             "Has `combine_types_mixfac()` already been applied?")
    }

    ## Prefix item names with the item type. paste0() with a NULL argument
    ## returns a length-1 string rather than character(0), so an unlabelled
    ## object would otherwise get a single name assigned to a dimension of
    ## extent I. `item_dim` is 1 for lambda_* (item, factor) and 2 for
    ## alpha_* (period, item).
    prefix_items <- function (v, type, item_dim) {
        dn <- dimnames(v)
        if (is.null(dn)) {
            dn <- vector("list", length(dim(v)))
        }
        nms <- dn[[item_dim]]
        if (is.null(nms)) {
            nms <- seq_len(dim(v)[item_dim])
        }
        dn[[item_dim]] <- paste0(type, ": ", nms)
        if (is.null(names(dn))) {
            names(dn) <- rep("", length(dn))
        }
        if (!nzchar(names(dn)[item_dim])) {
            names(dn)[item_dim] <- "item"
        }
        dimnames(v) <- dn
        v
    }

    ## lambda
    if (is_long) {
        lambda <- dplyr::bind_rows(x[lambda_idx], .id = "item_type")
        lambda$item_type <- sub("^lambda_", "", lambda$item_type)
    } else {
        for (k in seq_along(lambda_idx)) {
            type_k <- sub("^lambda_", "", names(x)[lambda_idx[k]])
            x[[lambda_idx[k]]] <-
                prefix_items(x[[lambda_idx[k]]], type_k, item_dim = 1)
        }
        lambda <- do.call(rbind, x[lambda_idx])
    }

    ## alpha
    if (is_long) {
        alpha <- dplyr::bind_rows(x[alpha_idx], .id = "item_type")
        alpha$item_type <- sub("^alpha_", "", alpha$item_type)
    } else {
        for (k in seq_along(alpha_idx)) {
            type_k <- sub("^alpha_", "", names(x)[alpha_idx[k]])
            x[[alpha_idx[k]]] <-
                prefix_items(x[[alpha_idx[k]]], type_k, item_dim = 2)
        }
        alpha <- do.call(cbind, x[alpha_idx])
    }

    ## kappa, and assembly
    if (is_long) {
        keep <- setdiff(seq_along(x), c(lambda_idx, alpha_idx, kappa_idx))
        if (length(kappa_idx) > 0) {
            kappa <- dplyr::bind_rows(x[kappa_idx], .id = "item_type")
            kappa$item_type <- sub("^kappa_", "", kappa$item_type)
            out <- c(list(lambda = lambda, alpha = alpha, kappa = kappa),
                     x[keep])
        } else {
            out <- c(list(lambda = lambda, alpha = alpha), x[keep])
        }
    } else {
        out <- as_draws_rvars_safe(
            eta = x$eta,
            lambda = lambda,
            alpha = alpha,
            kappa_trichot = x$kappa_trichot,
            kappa_ordinal = x$kappa_ordinal,
            sigma_alpha_evol = x$sigma_alpha_evol,
            sigma_metric = x$sigma_metric,
            Omega = x$Omega,
            lp__ = x$lp__
        )
    }

    out <- copy_dbmm_attrs(out, x)

    class(out) <- c("mixfac_comb", class(out))
    return(out)
}

#' Summarize the draws from a mixfac fit
#'
#' @param x mixfac draws
#' @param summary_functions Functions for summarizing. If missing, defaults will be used.
#' @export
summarize_mixfac <- function (x, summary_functions) {
    if (missing(summary_functions)) {
        summary_functions <- list(
            mean = ~posterior::E(.),
            median = ~stats::median(.),
            sd = ~posterior::sd(.),
            mad = ~posterior::mad(.),
            q5  = ~ as.numeric(stats::quantile(., probs = .05)),
            q95 = ~ as.numeric(stats::quantile(., probs = .95)),
            rhat = ~posterior::rhat(.),
            ess_bulk = ~posterior::ess_bulk(.),
            ess_tail = ~posterior::ess_tail(.)
        )
    }
    sfun <- function(y) {
        dplyr::mutate(
            y,
            dplyr::across(
                dplyr::all_of("value"), summary_functions, .names = "{.fn}"
            ),
            .keep = "unused"
        )
    }
    is_rvar <- sapply(x, inherits, "rvar")
    out_rvar <- out_df <- NULL
    if (any(is_rvar)) {
        out_rvar <- list()
        for (i in 1:sum(is_rvar)) {
            if (length(x[is_rvar][[i]]) == 0) {
                out_rvar[[i]] <- NA
            } else {
                nm <- names(x[is_rvar])[i]
                d <- posterior::as_draws_rvars(
                    stats::setNames(list(x[is_rvar][[i]]), nm)
                )
                out_rvar[[i]] <- do.call(
                    posterior::summarise_draws,
                    c(list(d), summary_functions)
                )
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

#' Order factors by sums of squared loadings
#'
#' Reorders the latent factors so that those accounting for more variance in
#' the loadings come first. This is a reparameterization: every linear
#' predictor, and hence the likelihood, is unchanged.
#'
#' @param x A `mixfac_comb` object, as returned by [combine_types_mixfac()].
#'     Must hold `rvar` draws, not the long format produced by
#'     `label_mixfac(make_long = TRUE)`.
#' @param check (logical) Check the class of `x`? Defaults to `TRUE`.
#'
#' @return A `mixfac_sorted` object with factors ordered by explanatory power.
#'     The permutation applied is recorded in the `"factor order"` attribute.
#'
#' @import posterior
#'
#' @export
sort_mixfac <- function(x, check = TRUE) {
    if (check) {
        check_arg_type(arg = x, typename = "mixfac_comb")
    }
    if (!posterior::is_draws_rvars(x)) {
        stop("`sort_mixfac()` requires rvar draws. It cannot be applied to ",
             "the long-format output of `label_mixfac(make_long = TRUE)`.")
    }
    n_factor <- dim(x$eta)[3]

    ## Sum of squared loadings per factor, descending
    ss <- posterior::rvar_apply(x$lambda, 2, function (y) {
        posterior::rvar_sum(y^2)
    })
    fo <- order(-posterior::E(ss))

    ## drop = FALSE throughout: with D == 1 these would otherwise collapse
    eta <- x$eta[, , fo, drop = FALSE]
    lambda <- x$lambda[, fo, drop = FALSE]
    Omega <- x$Omega[fo, fo, drop = FALSE]

    ## Factor labels are positional, so renumber rather than permute them
    if (!is.null(dimnames(eta))) dimnames(eta)[[3]] <- seq_along(fo)
    if (!is.null(dimnames(lambda))) dimnames(lambda)[[2]] <- seq_along(fo)
    if (!is.null(dimnames(Omega))) {
        dimnames(Omega)[[1]] <- dimnames(Omega)[[2]] <- seq_along(fo)
    }

    out <- as_draws_rvars_safe(
        eta = eta,
        lambda = lambda,
        alpha = x$alpha,
        kappa_trichot = x$kappa_trichot,
        kappa_ordinal = x$kappa_ordinal,
        sigma_alpha_evol = x$sigma_alpha_evol,
        sigma_metric = x$sigma_metric,
        Omega = Omega,
        lp__ = x$lp__
    )

    out <- copy_dbmm_attrs(out, x)
    attr(out, "rotation matrix") <- attr(x, "rotation matrix")
    attr(out, "signed-permutation matrix") <-
        attr(x, "signed-permutation matrix")
    attr(out, "factor order") <- fo

    class(out) <- unique(c("mixfac_sorted", class(x)))
    return(out)
}

#' Set the sign of each latent factor
#'
#' Flips the sign of each factor so that the mean loading on it has the
#' requested sign. This is a reparameterization: every linear predictor, and
#' hence the likelihood, is unchanged. Item intercepts and thresholds are
#' unaffected, since the sign matrix is orthogonal.
#'
#' @param x A `mixfac_comb` object, as returned by [combine_types_mixfac()].
#'     Must hold `rvar` draws.
#' @param signs (numeric) Desired sign of the mean loading on each factor,
#'     either a scalar (recycled) or one value per factor. Defaults to `1`.
#' @param check (logical) Check the class of `x`? Defaults to `TRUE`.
#'
#' @return A `mixfac_signed` object. The flips applied are recorded in the
#'     `"sign flips"` attribute.
#'
#' @import posterior
#'
#' @export
sign_mixfac <- function(x, signs = 1, check = TRUE) {
    if (check) {
        check_arg_type(arg = x, typename = "mixfac_comb")
    }
    if (!posterior::is_draws_rvars(x)) {
        stop("`sign_mixfac()` requires rvar draws. It cannot be applied to ",
             "the long-format output of `label_mixfac(make_long = TRUE)`.")
    }
    n_factor <- dim(x$eta)[3]
    stopifnot(length(signs) == 1 || length(signs) == n_factor)
    if (!all(signs %in% c(-1, 1))) {
        stop("`signs` must contain only -1 and 1.")
    }

    init_signs <- sign(colMeans(posterior::E(x$lambda)))
    sign_flips <- ifelse(init_signs == signs, 1, -1)
    sm <- diag(sign_flips, nrow = n_factor, ncol = n_factor)

    eta <- rotate_eta_rvar(x$eta, sm)
    dimnames(eta) <- dimnames(x$eta)
    lambda <- posterior::`%**%`(x$lambda, sm)
    dimnames(lambda) <- dimnames(x$lambda)
    Omega <- t(sm) %**% x$Omega %**% sm
    dimnames(Omega) <- dimnames(x$Omega)

    out <- as_draws_rvars_safe(
        eta = eta,
        lambda = lambda,
        alpha = x$alpha,
        kappa_trichot = x$kappa_trichot,
        kappa_ordinal = x$kappa_ordinal,
        sigma_alpha_evol = x$sigma_alpha_evol,
        sigma_metric = x$sigma_metric,
        Omega = Omega,
        lp__ = x$lp__
    )

    out <- copy_dbmm_attrs(out, x)
    attr(out, "rotation matrix") <- attr(x, "rotation matrix")
    attr(out, "signed-permutation matrix") <-
        attr(x, "signed-permutation matrix")
    attr(out, "sign flips") <- sign_flips
    attr(out, "factor order") <- attr(x, "factor order")

    class(out) <- unique(c("mixfac_signed", class(x)))
    return(out)
}


#' Attach unit, period, and item labels to modgirt draws
#'
#' @param x A `draws_rvars` object, typically from [identify_modgirt()].
#' @param check (logical) Check that `x` is a `draws_rvars` object?
#' @return A `modgirt_lab` object.
#' @export
label_modgirt <- function(x, check = TRUE) {
    if (check) stopifnot(posterior::is_draws_rvars(x))
    n_factor <- dim(x$bar_theta)[3]
    dimnames(x$bar_theta) <- list(
        period = attr(x, "time_labels"),
        unit = attr(x, "unit_labels"),
        factor = seq_len(n_factor)
    )
    dimnames(x$beta) <- list(
        item = attr(x, "item_labels"),
        factor = seq_len(n_factor)
    )
    for (nm in c("Sigma_theta", "Omega")) {
        if (!is.null(x[[nm]])) {
            dimnames(x[[nm]]) <- list(factor = seq_len(n_factor),
                                      factor = seq_len(n_factor))
        }
    }
    class(x) <- unique(c("modgirt_lab", class(x)))
    x
}

#' Identify MODGIRT draws
#'
#' This function identifies the MODGIRT model by postprocessing the draws from
#' the posterior distribution.
#'
#' @param x A fitted MODGIRT model object or `draws_rvars` object
#' @param method (string) Rotation criterion passed to
#'     [GPArotation::GPFRSorth()]. Defaults to `"varimax"`.
#' @param random_starts (non-negative integer) Number of random starting
#'     rotations tried per draw. Defaults to `0`, which starts from the
#'     identity matrix and is deterministic. Values greater than `0` guard
#'     against local minima of the rotation criterion, at proportionate
#'     computational cost.
#' @param seed (positive integer or `NULL`) Seed for the random starts.
#'     Defaults to `123`, making results reproducible. Ignored when
#'     `random_starts = 0`, since the rotation is then deterministic. Set to
#'     `NULL` to use the ambient random number stream.
#' @param quiet (logical) Suppress the iteration log printed by
#'     factor.switching::rsp_exact() Defaults to `TRUE`.
#'
#' @return A list containing the identified MODGIRT model parameters.
#'
#' @import posterior
#'
#' @export
identify_modgirt <- function(x, method = "varimax", random_starts = 0, seed = 123, quiet = TRUE) {
    ## Accept either a draws_rvars object or a fitted object. Label attributes
    ## live on the cmdstanr fit rather than on the wrapper list, so the source
    ## of the draws and the source of the labels differ in the latter case.
    if (posterior::is_draws_rvars(x)) {
        draws_rvar <- x
        label_src <- x
    } else if (!is.null(x$fit)) {
        draws_rvar <- posterior::as_draws_rvars(x$fit$draws())
        label_src <- x$fit
    } else {
        stop("`x` must be a draws_rvars object or a fitted model object ",
             "containing a `fit` element.")
    }
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
        vm_rvar <- with_seed(seed, make_vm_rvar(
            draws_of_beta, n_iter, n_chain, n_factor,
            method = method, randomStarts = random_starts
        ))
    } else {
        vm_rvar <- matrix(1)
    }
    ## Apply varimax rotations to `beta`
    beta_rvar$beta <- posterior::`%**%`(beta_rvar$beta, vm_rvar)
    ## Strip dimnames before flattening: rsp_exact() requires positional
    ## LambdaVd_q names, and rename_loading_matrix()'s pattern matches only
    ## numeric indices. Item labels are present whenever `x` is itself the
    ## output of identify_modgirt() or label_modgirt().
    beta_unnamed <- beta_rvar$beta
    dimnames(beta_unnamed) <- NULL
    beta_matrix <- posterior::as_draws_matrix(t(beta_unnamed))
    lambda_matrix <- rename_loading_matrix(beta_matrix)
    ## Create draw-specific signed permutations
    if (quiet) {
        invisible(
            utils::capture.output(
                       rsp_out <- factor.switching::rsp_exact(
                                                        lambda_matrix,
                                                        rotate = FALSE
                                                    )
                   ))
    } else {
        rsp_out <- factor.switching::rsp_exact(lambda_matrix, rotate = FALSE)
    }
    sp_rvar <- make_sp_rvar(rsp_out, n_iter, n_chain, n_factor)
    ## Apply signed permutations to `beta`
    beta_rvar$beta <- posterior::`%**%`(beta_rvar$beta, sp_rvar)
    dimnames(beta_rvar$beta) <- list(
        item = attr(label_src, "item_labels"),
        factor = seq_len(n_factor)
    )
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
    out <- as_draws_rvars_safe(
        lp__ = modgirt_rvar$lp__,
        alpha = modgirt_rvar$alpha,
        beta = beta_rvar$beta,
        bar_theta = bar_theta_rvar$bar_theta,
        Sigma_theta = sigma_theta_rvar$Sigma_theta,
        Omega = omega_rvar$Omega
    )
    out <- copy_dbmm_attrs(out, label_src)
    attr(out, "rotation matrix") <- vm_rvar
    attr(out, "signed-permutation matrix") <- sp_rvar
    class(out) <- c("modgirt_id", class(out))
    return(out)
}


#' Apply rotation to MODGIRT draws
#'
#' This function applies the given rotation to each draw from the posterior
#' distribution of the MODGIRT parameters
#'
#' @param modgirt_rvar A `draws_rvar` object from a MODGIRT model
#' @param rotmat An I-by-D rotation matrix
#'
#' @return A `draws_rvar` object of rotated draws
#'
#' @examples
#' \dontrun{
#' rotmat <- varimax(posterior::E(modgirt_signed$beta))$rotmat
#' modgirt_rotated <- rotate_modgirt(modgirt_signed, rotmat)
#' }
#' 
#' @import posterior
#'
#' @export
rotate_modgirt <- function(modgirt_rvar, rotmat) {
    ## inverse of transpose (needed for oblique rotation)
    if (is_rvar(rotmat)) {
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
    modgirt_rvar_rot <- posterior::draws_rvars(
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
        beta = modgirt_rvar$beta[, fo, drop = FALSE],
        bar_theta = modgirt_rvar$bar_theta[, , fo, drop = FALSE],
        Sigma_theta = modgirt_rvar$Sigma_theta[fo, fo, drop = FALSE],
        Omega = modgirt_rvar$Omega[fo, fo, drop = FALSE]
    )

    n_f <- length(fo)
    sorted_rvar$beta <- renumber_factor_dimnames(sorted_rvar$beta, 2L, n_f)
    sorted_rvar$bar_theta <- renumber_factor_dimnames(sorted_rvar$bar_theta, 3L, n_f)
    sorted_rvar$Sigma_theta <- renumber_factor_dimnames(sorted_rvar$Sigma_theta, 1:2, n_f)
    sorted_rvar$Omega <- renumber_factor_dimnames(sorted_rvar$Omega, 1:2, n_f)
    sorted_rvar <- copy_dbmm_attrs(sorted_rvar, modgirt_rvar)

    attr(sorted_rvar, "rotation matrix") <- attr(modgirt_rvar, "rotation matrix")
    attr(sorted_rvar, "signed-permutation matrix") <-
        attr(modgirt_rvar, "signed-permutation matrix")
    attr(sorted_rvar, "factor order") <- fo
    class(sorted_rvar) <- unique(c("modgirt_sorted", class(modgirt_rvar)))

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
    stopifnot(length(signs) == 1 || length(signs) == n_factor)
    if (!all(signs %in% c(-1, 1))) {
        stop("`signs` must contain only -1 and 1.")
    }
    init_signs <- sign(colMeans(posterior::E(modgirt_rvar$beta)))
    sign_flips <- ifelse(init_signs == signs, 1, -1)
    sm <- diag(sign_flips, nrow = n_factor, ncol = n_factor)

    out <- posterior::draws_rvars(
        lp__ = modgirt_rvar$lp__,
        alpha = modgirt_rvar$alpha,
        beta = posterior::`%**%`(modgirt_rvar$beta, sm),
        bar_theta = rotate_eta_rvar(modgirt_rvar$bar_theta, sm),
        Sigma_theta = t(sm) %**% modgirt_rvar$Sigma_theta %**% sm,
        Omega = t(sm) %**% modgirt_rvar$Omega %**% sm
    )
    ## Preserve dimnames, which the matrix products drop
    dimnames(out$beta) <- dimnames(modgirt_rvar$beta)
    dimnames(out$bar_theta) <- dimnames(modgirt_rvar$bar_theta)
    dimnames(out$Sigma_theta) <- dimnames(modgirt_rvar$Sigma_theta)
    dimnames(out$Omega) <- dimnames(modgirt_rvar$Omega)

    out <- copy_dbmm_attrs(out, modgirt_rvar)
    attr(out, "rotation matrix") <- attr(modgirt_rvar, "rotation matrix")
    attr(out, "signed-permutation matrix") <-
        attr(modgirt_rvar, "signed-permutation matrix")
    attr(out, "sign flips") <- sign_flips
    class(out) <- unique(c("modgirt_signed", class(modgirt_rvar)))
    out
}

#' Save a fitted dbmm model to disk
#'
#' `cmdstanr` reads posterior draws from CSV files on demand, so a fitted
#' object saved with [saveRDS()] will be unusable in a later session once
#' \R has deleted its temporary directory. This function forces the draws
#' into memory before serializing.
#'
#' @param x A `mixfac_fit` object from [fit_mixfac()], or the list returned
#'     by [fit_modgirt()].
#' @param file (string) Path to write to.
#' @return `file`, invisibly.
#' @export
save_dbmm <- function(x, file) {
    if (is.null(x$fit)) {
        cli::cli_abort("{.arg x} has no {.field fit} element.")
    }
    invisible(x$fit$draws())
    for (m in c("sampler_diagnostics", "init", "profiles")) {
        try(invisible(x$fit[[m]]()), silent = TRUE)
    }
    saveRDS(x, file)
    invisible(file)
}
