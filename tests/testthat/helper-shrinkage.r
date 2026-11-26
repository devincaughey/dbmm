## Shared fixtures for shrinkage-prior tests. Simulated data have one or two
## true factors, so a fit with more factors should shrink the extras toward
## zero.

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

## Two factors, one strong (|loadings| about 1.0) and one weak (about 0.5)
sim_two_factor_long <- function(n_units = 200, n_items = 40, seed = 2026) {
    withr::with_seed(seed, {
        theta <- matrix(stats::rnorm(n_units * 2), n_units, 2)
        sgn <- function() sample(c(-1, 1), n_items, replace = TRUE)
        beta <- cbind(stats::runif(n_items, 0.8, 1.2) * sgn(),
                      stats::runif(n_items, 0.4, 0.6) * sgn())
        alpha <- stats::rnorm(n_items, 0, 0.5)
        grid <- expand.grid(unit = seq_len(n_units), item = seq_len(n_items))
        nu <- alpha[grid$item] + rowSums(beta[grid$item, ] * theta[grid$unit, ])
        data.frame(
            unit = sprintf("u%03d", grid$unit),
            period = 1L,
            item = sprintf("q%02d", grid$item),
            value = stats::rbinom(nrow(grid), 1, stats::pnorm(nu))
        )
    })
}

shape_sim_long <- function(long) {
    suppressMessages(shape_mixfac(
        long_data = long,
        unit_var = "unit",
        time_var = "period",
        item_var = "item",
        value_var = "value",
        make_indicator_for_zeros = FALSE,
        periods_to_estimate = 1L,
        min_cat_n = 0
    ))
}

sim_one_factor_shaped <- function() {
    if (!is.null(.shrink_cache$shaped)) return(.shrink_cache$shaped)
    out <- shape_sim_long(sim_one_factor_long())
    .shrink_cache$shaped <- out
    out
}

## Fit the simulated data with three factors and a shrinkage prior, cached.
## `n_true` picks the data: one true factor, or two (strong and weak).
test_shrinkage_fit <- function(loading_prior = c("horseshoe", "lasso"),
                               n_true = 1) {
    loading_prior <- match.arg(loading_prior)
    key <- paste0("fit_", loading_prior, "_", n_true)
    if (!is.null(.shrink_cache[[key]])) return(.shrink_cache[[key]])
    skip_if_no_cmdstan()
    shaped <- if (n_true == 1) {
        sim_one_factor_shaped()
    } else {
        shape_sim_long(sim_two_factor_long())
    }
    args <- list(
        shaped_data = shaped,
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
    if (loading_prior == "horseshoe") args$p0 <- n_true
    out <- suppressMessages(suppressWarnings(do.call(fit_mixfac, args)))
    .shrink_cache[[key]] <- out
    out
}
