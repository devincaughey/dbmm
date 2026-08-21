## Shared fixtures for mixfac tests. Fitting is slow, so each configuration is
## fit at most once per test run and cached in this environment.

skip_if_no_cmdstan <- function() {
    testthat::skip_if_not_installed("cmdstanr")
    ok <- tryCatch(!is.null(cmdstanr::cmdstan_version()),
                   error = function(e) FALSE)
    testthat::skip_if_not(ok, "CmdStan not available")
}

## social_outcomes_2020_2021 is long (st, year, outcome, value) and every item
## is continuous, so on its own it exercises only the metric code path. Four
## items are replaced by discretized versions -- two binary, one trichotomous,
## one five-category ordinal -- so that the binary, trichotomous, and ordinal
## paths (and hence kappa_trichot and kappa_ordinal) are covered.
test_long_data <- function(periods = 2020:2021) {
    data("social_outcomes_2020_2021", package = "dbmm", envir = environment())
    long <- social_outcomes_2020_2021
    long <- long[long$year %in% periods & !is.na(long$value), ]

    src <- c(binary1 = "x_overall_health",
             binary2 = "x_employment_rate",
             trichot = "x_gini_posttax",
             ordinal = "x_state_co2")
    n_cat <- c(binary1 = 2L, binary2 = 2L, trichot = 3L, ordinal = 5L)
    stopifnot(all(src %in% long$outcome))

    ## Quantile-based cuts, so no bin is empty for skewed items
    qcut <- function(x, k) {
        br <- stats::quantile(x, probs = seq(0, 1, length.out = k + 1),
                              na.rm = TRUE, names = FALSE)
        stopifnot(!anyDuplicated(br))
        br[1] <- -Inf
        br[length(br)] <- Inf
        as.integer(cut(x, breaks = br, include.lowest = TRUE))
    }

    derived <- list()
    for (nm in names(src)) {
        d <- long[long$outcome == src[[nm]], ]
        d$value <- qcut(d$value, n_cat[[nm]])
        d$outcome <- paste0("d_", nm)
        ## Fail loudly here rather than mysteriously in shape_mixfac()
        stopifnot(length(unique(d$value)) == n_cat[[nm]])
        derived[[nm]] <- d
    }

    long <- long[!long$outcome %in% src, ]
    do.call(rbind, c(list(long), derived))
}

## make_indicator_for_zeros = FALSE for determinism: several outcomes are rates
## that may include zeros, which would otherwise spawn unpredictable "_zi"
## binary items and set the corresponding metric observations to missing.
test_shaped_data <- function(periods = 2020:2021) {
    suppressMessages(shape_mixfac(
        long_data = test_long_data(periods),
        unit_var = "st",
        time_var = "year",
        item_var = "outcome",
        value_var = "value",
        standardize = TRUE,
        make_indicator_for_zeros = FALSE,
        periods_to_estimate = periods
    ))
}

## Cache keyed on the flag settings that matter for the patches under test.
.fit_cache <- new.env(parent = emptyenv())

test_fit <- function(separate_eta = FALSE,
                     constant_alpha = TRUE,
                     gen_log_lik = FALSE,
                     n_dim = 2,
                     periods = 2020:2021) {
    key <- paste(separate_eta, constant_alpha, gen_log_lik, n_dim,
                 min(periods), max(periods), sep = "_")
    if (!is.null(.fit_cache[[key]])) return(.fit_cache[[key]])
    skip_if_no_cmdstan()
    out <- suppressMessages(fit_mixfac(
        shaped_data = test_shaped_data(periods),
        n_dim = n_dim,
        chains = 2,
        iter_warmup = 200,
        iter_sampling = 200,
        separate_eta = separate_eta,
        constant_alpha = constant_alpha,
        gen_log_lik = gen_log_lik,
        refresh = 0,
        seed = 123
    ))
    .fit_cache[[key]] <- out
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
test_comb <- function(n_dim = 2) {
    key <- paste0("comb_", n_dim)
    if (!is.null(.fit_cache[[key]])) return(.fit_cache[[key]])
    d <- extract_mixfac_draws(test_fit(n_dim = n_dim))
    out <- combine_types_mixfac(label_mixfac(identify_mixfac(d)))
    .fit_cache[[key]] <- out
    out
}
