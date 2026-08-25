## Shared fixtures for mixfac tests. Fitting is slow, so each configuration is
## fit at most once per test run and cached in this environment.

## state_policies_2010_2012 contains binary, ordinal, and continuous policies,
## so no discretization is needed to exercise the mixed-type code paths. Only a
## handful of items per type are kept, since fitting all 111 would make the
## test suite far too slow.
test_mixfac_long_data <- function(periods = 2010:2012,
                                  n_per_type = c(binary = 4, trichot = 2,
                                                 ordinal = 1, metric = 3),
                                  n_units = 20) {
    data("state_policies_2010_2012", package = "dbmm", envir = environment())
    d <- state_policies_2010_2012
    d <- d[d$year %in% periods & !is.na(d$value_real), ]

    ## A subset of states keeps the fixture fast; alphabetical for determinism
    d <- d[d$state_abb %in% utils::head(sort(unique(d$state_abb)), n_units), ]

    ## Classify exactly as shape_mixfac() does, on the subset actually used
    n_u <- tapply(d$value_real, d$policy_variable,
                  function(x) length(unique(x)))
    type_of <- cut(n_u, breaks = c(0, 1, 2, 3, 10, Inf),
                   labels = c("dropped", "binary", "trichot",
                              "ordinal", "metric"))
    names(type_of) <- names(n_u)

    pick <- unlist(lapply(names(n_per_type), function(ty) {
        cand <- names(type_of)[!is.na(type_of) & type_of == ty]
        ## Fewest categories first: each extra ordinal level costs a threshold
        cand <- cand[order(n_u[cand], cand)]
        if (length(cand) < n_per_type[[ty]]) {
            stop("Only ", length(cand), " ", ty, " item(s) available in the ",
                 "fixture; ", n_per_type[[ty]], " requested.")
        }
        utils::head(cand, n_per_type[[ty]])
    }), use.names = FALSE)

    d[d$policy_variable %in% pick, ]
}

test_mixfac_shaped_data <- function(periods = 2010:2012) {
    suppressMessages(shape_mixfac(
        long_data = test_mixfac_long_data(periods),
        unit_var = "state_abb",
        time_var = "year",
        item_var = "policy_variable",
        value_var = "value_real",
        standardize = TRUE,
        make_indicator_for_zeros = FALSE,
        periods_to_estimate = periods
    ))
}

## Cache keyed on the flag settings that matter for the patches under test.
.mixfac_cache <- new.env(parent = emptyenv())

test_mixfac_fit <- function(smooth_eta = TRUE,
                     constant_alpha = TRUE,
                     gen_log_lik = FALSE,
                     n_dim = 2,
                     periods = 2010:2012) {
    key <- paste(smooth_eta, constant_alpha, gen_log_lik, n_dim,
                 min(periods), max(periods), sep = "_")
    if (!is.null(.mixfac_cache[[key]])) return(.mixfac_cache[[key]])
    skip_if_no_cmdstan()
    out <- suppressMessages(fit_mixfac(
        shaped_data = test_mixfac_shaped_data(periods),
        n_dim = n_dim,
        chains = 2,
        parallel_chains = 2,
        iter_warmup = 200,
        iter_sampling = 200,
        smooth_eta = smooth_eta,
        constant_alpha = constant_alpha,
        gen_log_lik = gen_log_lik,
        refresh = 0,
        seed = 123
    ))
    .mixfac_cache[[key]] <- out
    out
}

## Item type with the most items, for invariance checks that need several items
richest_item_type <- function(x) {
    types <- c("binary", "trichot", "ordinal", "metric")
    n <- vapply(types, function(ty) {
        l <- x[[paste0("lambda_", ty)]]
        if (is.null(l)) 0L else dim(l)[1]
    }, integer(1))
    if (all(n == 0)) return(NA_character_)
    types[which.max(n)]
}

## nu = alpha + lambda %*% eta for one cell
nu_cell <- function(o, type, t = 1, j = 1, i = 1) {
    a <- o[[paste0("alpha_", type)]][t, i]
    l <- o[[paste0("lambda_", type)]][i, , drop = TRUE]
    e <- o$eta[t, j, , drop = TRUE]
    a + posterior::`%**%`(l, e)
}

## nu = alpha + lambda %*% eta for a combined (mixfac_comb) object
nu_cell_comb <- function(o, t = 1, j = 1, i = 1) {
    a <- o$alpha[t, i]
    l <- o$lambda[i, , drop = TRUE]
    e <- o$eta[t, j, , drop = TRUE]
    a + posterior::`%**%`(l, e)
}

## Fully processed combined object, cached
test_mixfac_comb <- function(n_dim = 2) {
    key <- paste0("comb_", n_dim)
    if (!is.null(.mixfac_cache[[key]])) return(.mixfac_cache[[key]])
    d <- extract_mixfac_draws(test_mixfac_fit(n_dim = n_dim))
    out <- combine_types_mixfac(label_mixfac(identify_mixfac(d)))
    .mixfac_cache[[key]] <- out
    out
}

test_mixfac_loo <- function(group = "dyad") {
    key <- paste0("loo_", group)
    if (!is.null(.mixfac_cache[[key]])) return(.mixfac_cache[[key]])
    ## Elevated Pareto k is expected: each dyad bundles several observations
    ## and the fixture has few draws
    out <- suppressWarnings(
        loo_mixfac(test_mixfac_fit(gen_log_lik = TRUE), group = group)
    )
    .mixfac_cache[[key]] <- out
    out
}
