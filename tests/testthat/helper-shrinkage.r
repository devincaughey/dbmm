## Shared fixtures for shrinkage-prior tests. Simulated data have one true
## factor, so a fit with more factors should shrink the extras toward zero.

.shrink_cache <- new.env(parent = emptyenv())

## Long-format binary responses from a one-factor probit model, one period
sim_one_factor_long <- function(n_units = 150, n_items = 30, seed = 2026) {
    withr::with_seed(seed, {
        theta <- stats::rnorm(n_units)
        beta <- stats::runif(n_items, 0.8, 1.5) *
            sample(c(-1, 1), n_items, replace = TRUE)
        alpha <- stats::rnorm(n_items, 0, 0.5)
        grid <- expand.grid(unit = seq_len(n_units), item = seq_len(n_items))
        p <- stats::pnorm(alpha[grid$item] + beta[grid$item] * theta[grid$unit])
        data.frame(
            unit = sprintf("u%03d", grid$unit),
            period = 1L,
            item = sprintf("q%02d", grid$item),
            value = stats::rbinom(nrow(grid), 1, p)
        )
    })
}

sim_one_factor_shaped <- function() {
    if (!is.null(.shrink_cache$shaped)) return(.shrink_cache$shaped)
    out <- suppressMessages(shape_mixfac(
        long_data = sim_one_factor_long(),
        unit_var = "unit",
        time_var = "period",
        item_var = "item",
        value_var = "value",
        make_indicator_for_zeros = FALSE,
        periods_to_estimate = 1L,
        min_cat_n = 0
    ))
    .shrink_cache$shaped <- out
    out
}

## Fit the simulated data with three factors and a shrinkage prior, cached
test_shrinkage_fit <- function(loading_prior = c("horseshoe", "lasso")) {
    loading_prior <- match.arg(loading_prior)
    key <- paste0("fit_", loading_prior)
    if (!is.null(.shrink_cache[[key]])) return(.shrink_cache[[key]])
    skip_if_no_cmdstan()
    args <- list(
        shaped_data = sim_one_factor_shaped(),
        n_dim = 3,
        loading_prior = loading_prior,
        chains = 2,
        parallel_chains = 2,
        iter_warmup = 300,
        iter_sampling = 300,
        adapt_delta = 0.95,
        refresh = 0,
        show_messages = FALSE,
        seed = 123
    )
    if (loading_prior == "horseshoe") args$p0 <- 1
    out <- suppressMessages(suppressWarnings(do.call(fit_mixfac, args)))
    .shrink_cache[[key]] <- out
    out
}
