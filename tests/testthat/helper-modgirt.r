## Shared fixtures for modgirt tests. Unlike mixfac, modgirt consumes
## individual-level ordinal responses, which shape_modgirt() aggregates into a
## weighted TIME x UNIT x ITEM x CATEGORY count array. social_outcomes_2020_2021
## is already aggregated and all-metric, so it cannot be reused here; responses
## are simulated instead.

## Simulate individual ordinal responses from a group-level IRT model with a
## drifting latent mean, so that the dynamic parameters are identified.
test_modgirt_long_data <- function(n_unit = 5, n_time = 3, n_item = 6,
                                   n_cat = 3, n_resp = 50, n_factor = 2,
                                   seed = 20260823) {
    withr::local_seed(seed)
    units <- LETTERS[seq_len(n_unit)]
    items <- paste0("item", seq_len(n_item))

    ## Unit-by-period factor scores: independent random walks
    theta <- array(NA_real_, dim = c(n_time, n_unit, n_factor))
    theta[1, , ] <- stats::rnorm(n_unit * n_factor)
    for (t in seq_len(n_time)[-1]) {
        theta[t, , ] <- theta[t - 1, , ] +
            stats::rnorm(n_unit * n_factor, sd = 0.3)
    }

    ## Item parameters. Loadings are drawn once and held fixed across periods.
    lambda <- matrix(stats::rnorm(n_item * n_factor, sd = 0.8),
                     nrow = n_item, ncol = n_factor)
    cuts <- stats::qnorm(seq_len(n_cat - 1) / n_cat)

    grid <- expand.grid(
        respondent = seq_len(n_resp),
        unit = units,
        time = seq_len(n_time),
        item = items,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    ui <- match(grid$unit, units)
    ii <- match(grid$item, items)
    nu <- numeric(nrow(grid))
    for (d in seq_len(n_factor)) {
        theta_d <- theta[cbind(grid$time, ui, d)]   # length nrow(grid)
        nu <- nu + theta_d * lambda[ii, d]
    }
    for (d in seq_len(n_factor)[-1]) {
        nu <- nu + theta[cbind(grid$time, ui, rep(d, nrow(grid)))] *
            lambda[ii, d]
    }
    grid$value <- as.integer(cut(nu + stats::rnorm(nrow(grid)),
                                 breaks = c(-Inf, cuts, Inf)))
    grid$wt <- 1

    ## Every cell of the count array must be reachable, or Stan will see
    ## structural zeros that the model cannot explain.
    stopifnot(length(unique(grid$value)) == n_cat)
    grid[c("unit", "time", "item", "value", "wt")]
}

test_modgirt_shaped_data <- function(...) {
    suppressMessages(shape_modgirt(
        long_data = test_modgirt_long_data(...),
        unit_var = "unit",
        time_var = "time",
        item_var = "item",
        value_var = "value",
        weight_var = "wt"
    ))
}

.modgirt_cache <- new.env(parent = emptyenv())

test_modgirt_fit <- function(n_factor = 2, n_dim_data = n_factor,
                             gen_log_lik = FALSE) {
    key <- paste0("fit_", n_factor, "_", n_dim_data, "_", gen_log_lik)
    if (!is.null(.modgirt_cache[[key]])) return(.modgirt_cache[[key]])
    skip_if_no_cmdstan()
    out <- suppressMessages(fit_modgirt(
        stan_data = test_modgirt_shaped_data(n_factor = n_dim_data),
        n_factor = n_factor,
        gen_log_lik = gen_log_lik,
        chains = 2,
        parallel_chains = 2,
        iter_warmup = 200,
        iter_sampling = 200,
        refresh = 0,
        seed = 123
    ))
    .modgirt_cache[[key]] <- out
    out
}

test_modgirt_loo <- function(group = "dyad") {
    key <- paste0("loo_", group)
    if (!is.null(.modgirt_cache[[key]])) return(.modgirt_cache[[key]])
    out <- suppressWarnings(
        loo_modgirt(test_modgirt_fit(gen_log_lik = TRUE), group = group)
    )
    .modgirt_cache[[key]] <- out
    out
}

## Identified draws, cached
test_modgirt_id <- function(n_factor = 2) {
    key <- paste0("id_", n_factor)
    if (!is.null(.modgirt_cache[[key]])) return(.modgirt_cache[[key]])
    out <- identify_modgirt(test_modgirt_fit(n_factor))
    .modgirt_cache[[key]] <- out
    out
}
