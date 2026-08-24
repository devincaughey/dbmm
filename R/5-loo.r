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

    absent <- setdiff(cols, names(idx))
    if (length(absent)) {
        cli::cli_abort(
            "{.fn log_lik_index} is missing column{?s} {.field {absent}},
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
#' @return A `draws_array` with dimensions iteration by chain by group. The
#'     grouping is recorded in the `"dbmm_group"` attribute.
#'
#' @export
log_lik_draws <- function(x, group = "dyad", shaped_data = NULL) {
    check_arg_type(arg = x, typename = "mixfac_fit")
    group <- match.arg(group, c("dyad", "observation", "unit", "period",
                                "item", "item_type"))

    ## Check for log_lik first: it needs nothing but `x`, and its absence is
    ## the more fundamental problem to report
    if (!"log_lik" %in% x$fit$metadata()$stan_variables) {
        cli::cli_abort(c(
            "{.field log_lik} is not present in the draws.",
            "i" = "Refit with {.code gen_log_lik = TRUE}."
        ))
    }
    ll <- x$fit$draws("log_lik", format = "draws_array")
    if (group == "observation") {
        attr(ll, "dbmm_group") <- group
        return(ll)
    }

    if (is.null(shaped_data)) shaped_data <- x$shaped_data
    if (is.null(shaped_data) || !length(shaped_data)) {
        cli::cli_abort(c(
            "Cannot group {.field log_lik} without the data {.arg x} was fit to.",
            "i" = "Refit with {.code return_data = TRUE}, or pass
                   {.arg shaped_data} directly.",
            "i" = "{.code group = \"observation\"} needs no data and works
                   without it."
        ))
    }
    check_arg_type(arg = shaped_data, typename = "mixfac_data")

    idx <- log_lik_index(shaped_data)
    if (nrow(idx) != dim(ll)[3L]) {
        cli::cli_abort(c(
            "{.fn log_lik_index} returned {nrow(idx)} row{?s} but
             {.field log_lik} has {dim(ll)[3L]} element{?s}.",
            "i" = "{.arg shaped_data} may not be the data {.arg x} was fit to."
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
        ## rowSums(dims = 2) sums over the third margin, and is much faster
        ## than apply() for the many-column case
        out[, , g] <- rowSums(ll[, , j_by_group[[g]], drop = FALSE], dims = 2L)
    }

    out <- posterior::as_draws_array(out)
    attr(out, "dbmm_group") <- group
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
    drop_nonfinite <- check_flag(drop_nonfinite)
    ll <- log_lik_draws(x, group = group, shaped_data = shaped_data)

    n_iter_before <- dim(ll)[1L]
    a <- check_log_lik_finite(unclass(ll), drop = drop_nonfinite)
    dropped <- dim(a)[1L] < n_iter_before

    dots <- list(...)
    if (is.null(dots$r_eff)) {
        if (dropped) {
            ## relative_eff() estimates autocorrelation from a contiguous
            ## chain; after dropping iterations the series has holes, so any
            ## ESS it returns is meaningless. Fall back to the conservative
            ## assumption of independent draws.
            cli::cli_warn(c(
                "Using {.code r_eff = 1} because iterations were discarded.",
                "i" = "Dropping iterations breaks the chain structure that
                       {.fn loo::relative_eff} assumes, so the Pareto k
                       diagnostics no longer account for autocorrelation."
            ))
            dots$r_eff <- rep(1, dim(a)[3L])
        } else {
            ## Shift each group to a maximum of zero before exponentiating:
            ## summed log-likelihoods can be far below log(.Machine$double.xmin),
            ## and exp() would underflow to a column of zeros, making the ESS
            ## undefined. ESS is invariant to per-group rescaling.
            shift <- apply(a, 3L, max)
            dots$r_eff <- loo::relative_eff(exp(sweep(a, 3L, shift, "-")))
        }
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
