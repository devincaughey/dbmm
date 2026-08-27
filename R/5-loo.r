#' Map positions of `log_lik` to group-item-periods
#'
#' @param stan_data A list from [shape_modgirt()], as returned in the
#'     `stan_data` element of [fit_modgirt()].
#'
#' @return A data frame with one row per element of `log_lik`, giving the item,
#'     unit, period, and total (possibly weighted) number of responses.
#'
#' @details
#' Positions run in the order the Stan program fills them: period, then group,
#' then item, with item varying fastest. Cells with no responses are omitted,
#' since they contribute nothing to the likelihood.
#'
#' @export
log_lik_index_modgirt <- function(stan_data) {
    if (is.null(stan_data$SSSS)) {
        cli::cli_abort("{.arg stan_data} has no {.field SSSS} element.")
    }
    dm <- dim(stan_data$SSSS)
    ## Item varies fastest, matching Stan's for (t) for (g) for (q) ordering
    grid <- expand.grid(item = seq_len(dm[3L]),
                        unit = seq_len(dm[2L]),
                        period = seq_len(dm[1L]),
                        KEEP.OUT.ATTRS = FALSE)
    n_resp <- apply(stan_data$SSSS, c(1L, 2L, 3L), sum)   # [T, G, Q]
    grid$n_responses <-
        n_resp[cbind(grid$period, grid$unit, grid$item)]
    grid <- grid[grid$n_responses > 0, , drop = FALSE]

    lab <- function(nm, i) {
        v <- attr(stan_data, nm)
        if (is.null(v)) as.character(i) else v[i]
    }
    out <- data.frame(
        item        = lab("item_labels", grid$item),
        unit        = lab("unit_labels", grid$unit),
        period      = lab("time_labels", grid$period),
        n_responses = grid$n_responses,
        stringsAsFactors = FALSE
    )
    out$position <- seq_len(nrow(out))
    row.names(out) <- NULL
    out
}

#' Group log-likelihood columns
#'
#' @param idx (data frame) Output of [log_lik_index()].
#' @param group (string) Grouping level.
#' @param call (environment) Frame to report errors against.
#' @return A factor, in order of first appearance, or `NULL` for
#'     `"observation"`.
#' @keywords internal
log_lik_group_key <- function(idx, group, call = rlang::caller_env()) {
    if (identical(group, "observation")) return(NULL)

    ## Item labels are not unique across item types, so include the type
    cols <- switch(
        group,
        dyad      = c("item_type", "item", "unit"),
        item      = c("item_type", "item"),
        unit      = "unit",
        period    = "period",
        item_type = "item_type",
        cli::cli_abort(
            "Unrecognized {.arg group}: {.val {group}}.",
            call = call
        )
    )

    ## Models with a single response type (modgirt) have no item_type column;
    ## item labels are then unique on their own
    if ("item_type" %in% cols && !"item_type" %in% names(idx)) {
        if (identical(group, "item_type")) {
            cli::cli_abort(
                     "{.code group = \"item_type\"} is not available for this model,
                 which has a single response type.",
                 call = call
                 )
        }
        cols <- setdiff(cols, "item_type")
    }

    absent <- setdiff(cols, names(idx))
    if (length(absent)) {
        cli::cli_abort(
                 "The log-likelihood index is missing column{?s} {.field {absent}},
             needed for {.code group = {.val {group}}}.",
             call = call
             )
    }

    parts <- lapply(idx[cols], as.character)
    key <- do.call(paste, c(parts, list(sep = " | ")))

    ## The pasted key must be one-to-one with the underlying columns; a label
    ## containing the separator would otherwise merge two distinct groups
    if (length(unique(key)) != nrow(unique(as.data.frame(parts)))) {
        cli::cli_abort(
            c("Cannot form unique {.arg group} labels.",
              "i" = "One or more of {.field {cols}} contains the separator
                     {.val { | }}."),
            call = call
        )
    }

    factor(key, levels = unique(key))
}

#' Aggregate log-likelihood draws by group
#'
#' @param fit A `CmdStanMCMC` object.
#' @param group (string) Grouping level, already validated.
#' @param idx (data frame or `NULL`) Log-likelihood index, as returned by
#'     [log_lik_index()] or [log_lik_index_modgirt()]. `NULL` is an error
#'     unless `group` is `"observation"`.
#' @param data_hint (string) Advice to give when `idx` is `NULL`.
#' @keywords internal
#' @noRd
log_lik_draws_impl <- function(fit, group, idx, data_hint) {
    if (!"log_lik" %in% fit$metadata()$stan_variables) {
        cli::cli_abort(c(
            "{.field log_lik} is not present in the draws.",
            "i" = "Refit with {.code gen_log_lik = TRUE}."
        ))
    }
    ll <- fit$draws("log_lik", format = "draws_array")
    if (group == "observation") {
        attr(ll, "dbmm_group") <- group
        return(ll)
    }

    if (is.null(idx)) {
        cli::cli_abort(c(
            "Cannot group {.field log_lik} without the data the model was
             fit to.",
            "i" = data_hint,
            "i" = "{.code group = \"observation\"} needs no data and works
                   without it."
        ))
    }
    if (nrow(idx) != dim(ll)[3L]) {
        cli::cli_abort(c(
            "The log-likelihood index has {nrow(idx)} row{?s} but
             {.field log_lik} has {dim(ll)[3L]} element{?s}.",
            "i" = "The supplied data may not be what the model was fit to."
        ))
    }
    key <- log_lik_group_key(idx, group)

    out <- array(
        NA_real_,
        dim = c(dim(ll)[1L], dim(ll)[2L], nlevels(key)),
        dimnames = list(
            iteration = dimnames(ll)[[1L]],
            chain     = dimnames(ll)[[2L]],
            variable  = levels(key)
        )
    )
    j_by_group <- split(seq_len(nrow(idx)), key)
    for (g in seq_along(j_by_group)) {
        out[, , g] <- rowSums(ll[, , j_by_group[[g]], drop = FALSE], dims = 2L)
    }

    out <- posterior::as_draws_array(out)
    attr(out, "dbmm_group") <- group
    out
}

#' Extract log-likelihood draws
#'
#' @param x A fitted model from [fit_mixfac()] or [fit_modgirt()], fit with
#'     `gen_log_lik = TRUE`.
#' @param group (string) Level at which to sum log-likelihoods, giving the
#'     joint predictive density of each group. One of `"dyad"` (the default),
#'     `"observation"`, `"unit"`, `"period"`, or `"item"`; mixfac fits also
#'     accept `"item_type"`.
#' @param ... Passed to methods.
#'
#' @return A `draws_array` with dimensions iteration by chain by group. The
#'     grouping is recorded in the `"dbmm_group"` attribute.
#'
#' @details
#' For a mixfac fit the finest level is one response, so `"dyad"` sums over
#' the periods in which a unit was observed on an item. For a modgirt fit the
#' finest level is a group-item-period cell, and each element is already a sum
#' over response categories weighted by the number of respondents; `"dyad"`
#' then sums over periods as well. The two are therefore on different scales,
#' and their `elpd_loo` values are not comparable.
#'
#' @export
log_lik_draws <- function(x, group = "dyad", ...) {
    UseMethod("log_lik_draws")
}

#' @param shaped_data (mixfac_data or `NULL`) The data the model was fit to.
#'     By default taken from `x`, which requires `return_data = TRUE`.
#' @rdname log_lik_draws
#' @export
log_lik_draws.mixfac_fit <- function(x, group = "dyad", shaped_data = NULL,
                                     ...) {
    group <- match.arg(group, c("dyad", "observation", "unit", "period",
                                "item", "item_type"))
    if (is.null(shaped_data)) shaped_data <- x$shaped_data
    idx <- if (!is.null(shaped_data) && length(shaped_data)) {
               check_arg_type(arg = shaped_data, typename = "mixfac_data")
               log_lik_index(shaped_data)
           } else NULL
    log_lik_draws_impl(
        x$fit, group, idx,
        data_hint = "Refit with {.code return_data = TRUE}, or pass
                     {.arg shaped_data} directly."
    )
}

#' @param stan_data (list or `NULL`) The data the model was fit to. By default
#'     taken from `x`, which requires `return_data = TRUE`.
#' @rdname log_lik_draws
#' @export
log_lik_draws.modgirt_fit <- function(x, group = "dyad", stan_data = NULL,
                                      ...) {
    group <- match.arg(group, c("dyad", "observation", "unit", "period",
                                "item"))
    if (is.null(stan_data)) stan_data <- x$stan_data
    idx <- if (!is.null(stan_data$SSSS)) {
               log_lik_index_modgirt(stan_data)
           } else NULL
    log_lik_draws_impl(
        x$fit, group, idx,
        data_hint = "Refit with {.code return_data = TRUE}, or pass
                     {.arg stan_data} directly."
    )
}

#' Compute PSIS-LOO from already-grouped log-likelihood draws
#'
#' @param ll A `draws_array` from [log_lik_draws()].
#' @param group (string) Grouping, for the returned attributes.
#' @param drop_nonfinite (logical) Discard offending iterations?
#' @param dots (list) Additional arguments for [loo::loo()].
#' @keywords internal
#' @noRd
loo_from_draws <- function(ll, group, drop_nonfinite, dots) {
    n_iter_before <- dim(ll)[1L]
    a <- check_log_lik_finite(unclass(ll), drop = drop_nonfinite)
    dropped <- dim(a)[1L] < n_iter_before

    if (is.null(dots$r_eff)) {
        if (dropped) {
            cli::cli_warn(c(
                "Using {.code r_eff = 1} because iterations were discarded.",
                "i" = "Dropping iterations breaks the chain structure that
                       {.fn loo::relative_eff} assumes, so the Pareto k
                       diagnostics no longer account for autocorrelation."
            ))
            dots$r_eff <- rep(1, dim(a)[3L])
        } else {
            ## Shift each group to a maximum of zero before exponentiating:
            ## summed log-likelihoods can fall below
            ## log(.Machine$double.xmin), where exp() underflows to a column
            ## of zeros and the ESS is undefined. ESS is invariant to
            ## per-group rescaling.
            shift <- apply(a, 3L, max)
            dots$r_eff <- loo::relative_eff(exp(sweep(a, 3L, shift, "-")))
        }
    }

    out <- do.call(loo::loo, c(list(a), dots))
    attr(out, "dbmm_group") <- group
    attr(out, "dbmm_groups") <- dimnames(a)[[3L]]
    out
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
#' Non-finite log-likelihoods are an error by default. With drop_nonfinite = TRUE
#' the affected iterations are removed from every chain, so that the array stays
#' rectangular; because this breaks the chain structure, r_eff is then set to 1
#' rather than estimated.
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
    check_arg_type(arg = x, typename = "mixfac_fit")
    ll <- log_lik_draws(x, group = group, shaped_data = shaped_data)
    loo_from_draws(ll, attr(ll, "dbmm_group"),
                   check_flag(drop_nonfinite), list(...))
}

#' Approximate leave-one-out cross-validation for a MODGIRT model
#'
#' Computes PSIS-LOO for a fitted group-level IRT model, using the \pkg{loo}
#' package.
#'
#' @param x A fitted model from [fit_modgirt()], fit with
#'     `gen_log_lik = TRUE`.
#' @param group (string) Level at which to sum log-likelihoods. One of
#'     `"dyad"` (the default; group-by-item, i.e. all periods of one item for
#'     one group), `"observation"`, `"unit"`, `"period"`, or `"item"`.
#' @param stan_data (list or `NULL`) The data the model was fit to. By default
#'     taken from `x`.
#' @param drop_nonfinite (logical) If `log_lik` contains non-finite values,
#'     should the affected iterations be discarded rather than triggering an
#'     error? Defaults to `FALSE`.
#' @param ... Additional arguments passed to [loo::loo()].
#'
#' @return A `loo` object.
#'
#' @details
#' The finest unit of the MODGIRT likelihood is a group-item-period cell: the
#' responses of one demographic group to one item in one period, aggregated
#' into counts. Holding out one such cell therefore removes every respondent
#' who answered that item in that group and period, so this is
#' leave-a-block-of-respondents-out rather than leave-one-respondent-out. Its
#' `elpd_loo` is not on the same scale as that of [loo_mixfac()] and the two
#' cannot be compared.
#'
#' When `shape_modgirt()` was given `weight_var`, the counts in `SSSS` are
#' weighted and generally non-integer, so the quantity being cross-validated
#' is a pseudo-likelihood. The Pareto k diagnostics remain informative about
#' the reliability of the importance sampling, but the implied effective
#' sample sizes should be read loosely.
#'
#' As with [loo_mixfac()], `group = "unit"` holds out a group's entire
#' response history, creating folds over which importance sampling is
#' unreliable; prefer K-fold cross-validation in that case. Models compared
#' with [loo::loo_compare()] must be fit to identical data and grouped
#' identically.
#'
#' @seealso [loo_mixfac()], [log_lik_index_modgirt()], [loo_influential()]
#'
#' @examples
#' \dontrun{
#' f1 <- fit_modgirt(stan_data, n_factor = 1, gen_log_lik = TRUE)
#' f2 <- fit_modgirt(stan_data, n_factor = 2, gen_log_lik = TRUE)
#' loo::loo_compare(list(d1 = loo_modgirt(f1), d2 = loo_modgirt(f2)))
#' }
#'
#' @export
loo_modgirt <- function(x, group = "dyad", stan_data = NULL,
                        drop_nonfinite = FALSE, ...) {
    if (!requireNamespace("loo", quietly = TRUE)) {
        cli::cli_abort("Package {.pkg loo} is required to use
                        {.fn loo_modgirt}.")
    }
    check_arg_type(arg = x, typename = "modgirt_fit")
    ll <- log_lik_draws(x, group = group, stan_data = stan_data)
    loo_from_draws(ll, attr(ll, "dbmm_group"),
                   check_flag(drop_nonfinite), list(...))
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
#'     Defaults to `0.7`. Groups with a missing k are always reported.
#'
#' @return A data frame with columns `group`, `pareto_k`, and, where
#'     available, `elpd_loo` and `n_eff`, in decreasing order of k.
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

    ## A missing k is at least as serious as a large one, and `k > threshold`
    ## would silently drop it
    i <- which(is.na(k) | k > threshold)
    out <- data.frame(group = nm[i], pareto_k = k[i],
                      stringsAsFactors = FALSE)

    if (!is.null(l$pointwise) && "elpd_loo" %in% colnames(l$pointwise)) {
        out$elpd_loo <- l$pointwise[i, "elpd_loo"]
    }
    if (!is.null(l$diagnostics$n_eff)) {
        out$n_eff <- l$diagnostics$n_eff[i]
    }

    out <- out[order(-out$pareto_k, na.last = FALSE), , drop = FALSE]
    row.names(out) <- NULL
    out
}

#' Check log-likelihood draws for non-finite values
#'
#' @param a (array) Iteration by chain by group array of log-likelihoods.
#' @param drop (logical) Discard offending iterations rather than aborting?
#'     An iteration is removed from *all* chains if any chain has a non-finite
#'     value at that position, so that the returned array stays rectangular.
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
