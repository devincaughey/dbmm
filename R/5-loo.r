#' Group log-likelihood columns
#'
#' @param idx (data frame) Output of [log_lik_index()].
#' @param group (string) Grouping level.
#' @return A factor, in order of first appearance, or `NULL` for
#'     `"observation"`.
#' @keywords internal
log_lik_group_key <- function(idx, group) {
    key <- switch(
        group,
        observation = return(NULL),
        ## Item labels are not unique across item types, so include the type
        dyad      = paste(idx$item_type, idx$item, idx$unit, sep = " | "),
        item      = paste(idx$item_type, idx$item, sep = " | "),
        unit      = idx$unit,
        period    = idx$period,
        item_type = idx$item_type
    )
    factor(key, levels = unique(key))
}

#' Extract per-observation log-likelihood draws
#'
#' @param x (mixfac_fit) A fitted model from [fit_mixfac()], fit with
#'     `gen_log_lik = TRUE`.
#' @param group (string) Level at which to sum log-likelihoods, giving the
#'     joint predictive density of each group. One of `"dyad"` (the default;
#'     unit-by-item, i.e. all periods of one item for one unit),
#'     `"observation"`, `"unit"`, `"period"`, `"item"`, or `"item_type"`.
#' @param shaped_data (mixfac_data or `NULL`) The data the model was fit to,
#'     needed to map positions of `log_lik` to observations. By default taken
#'     from `x`, which requires that it was fit with `return_data = TRUE`;
#'     supply it explicitly otherwise.
#'
#' @return A `draws_array` with dimensions iteration by chain by group.
#'
#' @export
log_lik_draws <- function(x, group = "dyad", shaped_data = NULL) {
    check_arg_type(arg = x, typename = "mixfac_fit")
    group <- match.arg(group, c("dyad", "observation", "unit", "period",
                                "item", "item_type"))
    if (is.null(shaped_data)) shaped_data <- x$shaped_data
    if (is.null(shaped_data) && group != "observation") {
        cli::cli_abort(c(
            "Cannot group {.field log_lik} without the data {.arg x} was fit to.",
            "i" = "Refit with {.code return_data = TRUE}, or pass
                   {.arg shaped_data} directly.",
            "i" = "{.code group = \"observation\"} needs no data and works
                   without it."
        ))
    }
    if (!"log_lik" %in% x$fit$metadata()$stan_variables) {
        cli::cli_abort(c(
            "{.field log_lik} is not present in the draws.",
            "i" = "Refit with {.code gen_log_lik = TRUE}."
        ))
    }
    ll <- x$fit$draws("log_lik", format = "draws_array")
    if (group == "observation") return(ll)

    check_arg_type(arg = shaped_data, typename = "mixfac_data")
    idx <- log_lik_index(shaped_data)
    if (nrow(idx) != dim(ll)[3L]) {
        cli::cli_abort(
            "{.fn log_lik_index} returned {nrow(idx)} row{?s} but
             {.field log_lik} has {dim(ll)[3L]} element{?s}."
        )
    }
    key <- log_lik_group_key(idx, group)
    if (is.null(key)) return(ll)

    out <- array(
        NA_real_,
        dim = c(dim(ll)[1L], dim(ll)[2L], nlevels(key)),
        dimnames = list(
            iteration = dimnames(ll)[[1L]],
            chain     = dimnames(ll)[[2L]],
            variable  = levels(key)
        )
    )
    for (g in seq_len(nlevels(key))) {
        j <- which(as.integer(key) == g)
        out[, , g] <- apply(ll[, , j, drop = FALSE], c(1L, 2L), sum)
    }
    posterior::as_draws_array(out)
}

#' Approximate leave-one-out cross-validation
#'
#' Computes PSIS-LOO for a fitted mixed factor model, using the \pkg{loo}
#' package. Relative effective sample sizes are computed from the chain
#' structure of the draws, so the Pareto k diagnostics account for
#' autocorrelation.
#'
#' @inheritParams log_lik_draws
#' @param ... Additional arguments passed to [loo::loo()].
#' @param drop_nonfinite (logical) If `log_lik` contains non-finite values,
#'     should the affected iterations be discarded rather than triggering an
#'     error? Defaults to `FALSE`. Discarding iterations conditions on the
#'     log-likelihood being finite and so biases the estimate; it is preferable
#'     to diagnose and fix the cause.
#' @rdname loo_mixfac
#'
#' @return A `loo` object.
#'
#' @details
#' The default, `group = "dyad"`, gives leave-one-unit-by-item-out: all
#' periods of a given item for a given unit are held out together. This is the
#' appropriate comparison for choosing the number of factors. Under
#' `group = "observation"`, a held-out response can be predicted by the same
#' unit and item in an adjacent period, since both the item intercepts and the
#' factor scores evolve smoothly; predictive accuracy then reflects temporal
#' interpolation more than factor structure. Holding out the whole dyad forces
#' prediction from the unit's responses to *other* items, which is what the
#' loadings govern.
#'
#' For comparing `smooth_eta = TRUE` against `FALSE`, `group = "observation"`
#' is the more direct test, since it asks precisely whether neighbouring
#' periods of the same unit and item improve prediction.
#'
#' `group = "unit"` holds out a unit's entire response history. With a typical
#' number of items and periods this creates blocks of many observations, over
#' which importance sampling is unreliable; expect poor Pareto k diagnostics
#' and prefer K-fold cross-validation in that case.
#'
#' Models compared with [loo::loo_compare()] must be fit to identical
#' observations and grouped identically. Comparing models that differ in
#' `n_dim`, `smooth_eta`, or `constant_alpha` is valid; comparing models fit to
#' different items or periods is not.
#'
#' @examples
#' \dontrun{
#' f1 <- fit_mixfac(shaped, n_dim = 1, gen_log_lik = TRUE)
#' f2 <- fit_mixfac(shaped, n_dim = 2, gen_log_lik = TRUE)
#' loo::loo_compare(
#'     d1 = loo_mixfac(f1, group = "unit"),
#'     d2 = loo_mixfac(f2, group = "unit")
#' )
#' }
#'
#' @export
loo_mixfac <- function(x, group = "dyad", shaped_data = NULL,
                       drop_nonfinite = FALSE, ...) {
    if (!requireNamespace("loo", quietly = TRUE)) {
        cli::cli_abort("Package {.pkg loo} is required to use
                        {.fn loo_mixfac}.")
    }
    drop_nonfinite <- check_flag(drop_nonfinite)
    ll <- log_lik_draws(x, group = group, shaped_data = shaped_data)
    a <- check_log_lik_finite(unclass(ll), drop = drop_nonfinite)

    dots <- list(...)
    if (is.null(dots$r_eff)) {
        dots$r_eff <- loo::relative_eff(exp(a))
    }
    out <- do.call(loo::loo, c(list(a), dots))
    attr(out, "dbmm_group") <- group
    attr(out, "dbmm_groups") <- dimnames(a)[[3L]]
    out
}

#' Identify influential groups in a LOO estimate
#'
#' Reports the groups whose Pareto k diagnostic exceeds a threshold, and so
#' whose contribution to the LOO estimate is unreliable. High values often
#' indicate observations the model fits poorly rather than a purely
#' computational problem.
#'
#' @param l (loo) Output of [loo_mixfac()].
#' @param threshold (numeric) Report groups with Pareto k above this value.
#'     Defaults to `0.7`.
#'
#' @return A data frame of group names and Pareto k values, in decreasing
#'     order of k.
#'
#' @export
loo_influential <- function(l, threshold = 0.7) {
    if (!requireNamespace("loo", quietly = TRUE)) {
        cli::cli_abort("Package {.pkg loo} is required.")
    }
    k <- loo::pareto_k_values(l)
    nm <- attr(l, "dbmm_groups")
    if (is.null(nm) || length(nm) != length(k)) {
        nm <- as.character(seq_along(k))
    }
    i <- which(k > threshold)
    out <- data.frame(group = nm[i], pareto_k = k[i],
                      stringsAsFactors = FALSE)
    out[order(-out$pareto_k), , drop = FALSE]
}

#' Check log-likelihood draws for non-finite values
#'
#' @param a (array) Iteration by chain by group array of log-likelihoods.
#' @param drop (logical) Discard offending iterations rather than aborting?
#' @param call (environment) Frame to report errors against.
#'
#' @return `a`, possibly with iterations removed.
#'
#' @keywords internal
check_log_lik_finite <- function(a, drop = FALSE, call = rlang::caller_env()) {
    bad <- !is.finite(a)
    if (!any(bad)) return(a)

    ## Non-finite values are usually all-or-nothing within a draw, since an
    ## exception in Stan's generated quantities block blanks the whole block.
    bad_iter <- which(apply(bad, 1L, any))
    n_iter <- dim(a)[1L]
    n_grp_affected <- sum(apply(bad, 3L, any))

    if (!drop) {
        cli::cli_abort(
            c(
                "{.field log_lik} contains {sum(bad)} non-finite value{?s}, in
                 {length(bad_iter)} of {n_iter} iteration{?s} and
                 {n_grp_affected} of {dim(a)[3L]} group{?s}.",
                "i" = "This normally reflects a problem with the fit rather
                       than with {.fn loo_mixfac}.",
                "i" = "If whole iterations are affected, an exception was
                       probably thrown in the generated quantities block; see
                       {.code x$fit$output()}.",
                "i" = "Set {.code drop_nonfinite = TRUE} to discard the
                       affected iterations, bearing in mind that this
                       conditions on the log-likelihood being finite and so
                       biases the estimate."
            ),
            call = call
        )
    }

    keep <- setdiff(seq_len(n_iter), bad_iter)
    if (length(keep) == 0L) {
        cli::cli_abort(
            "Every iteration contains non-finite log-likelihoods.",
            call = call
        )
    }
    cli::cli_warn(c(
        "Discarding {length(bad_iter)} of {n_iter} iteration{?s} with
         non-finite log-likelihoods.",
        "!" = "The result conditions on the log-likelihood being finite and is
               therefore biased. Prefer fixing the underlying problem."
    ))
    out <- a[keep, , , drop = FALSE]
    ## Guard against the residual case of scattered non-finite values, which
    ## dropping whole iterations would not remove
    if (any(!is.finite(out))) {
        cli::cli_abort(
            "Non-finite log-likelihoods remain after dropping affected
             iterations.",
            call = call
        )
    }
    out
}
