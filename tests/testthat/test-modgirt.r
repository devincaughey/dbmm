test_that("identify_modgirt() returns a classed draws_rvars object", {
    id <- test_modgirt_id(2)
    expect_s3_class(id, "modgirt_id")
    expect_true(posterior::is_draws_rvars(id))
    expect_true(all(c("beta", "bar_theta", "Sigma_theta", "Omega") %in%
                    posterior::variables(id)))
})

test_that("identify_modgirt() records the rotation matrices as attributes", {
    id <- test_modgirt_id(2)
    expect_equal(dim(attr(id, "rotation matrix")), c(2L, 2L))
    expect_equal(dim(attr(id, "signed-permutation matrix")), c(2L, 2L))
})

test_that("identify_modgirt() propagates label attributes from the fit", {
    id <- test_modgirt_id(2)
    expect_length(attr(id, "item_labels"), 6L)
    expect_length(attr(id, "unit_labels"), 5L)
    expect_length(attr(id, "time_labels"), 3L)
})

test_that("identify_modgirt() is reproducible with random starts", {
    f <- test_modgirt_fit(2)
    a <- identify_modgirt(f, random_starts = 2, seed = 1)
    b <- identify_modgirt(f, random_starts = 2, seed = 1)
    expect_equal(posterior::E(a$beta), posterior::E(b$beta))
})

test_that("identify_modgirt() does not disturb the caller's RNG state", {
    f <- test_modgirt_fit(2)
    set.seed(99)
    before <- .Random.seed
    invisible(identify_modgirt(f, random_starts = 2, seed = 1))
    expect_identical(.Random.seed, before)
})

test_that("sort_modgirt() preserves dimensions and renumbers factors", {
    id <- test_modgirt_id(2)
    s <- sort_modgirt(id)
    expect_equal(dim(s$beta), dim(id$beta))
    expect_equal(dim(s$bar_theta), dim(id$bar_theta))
    expect_equal(dim(s$Omega), dim(id$Omega))
    expect_equal(dimnames(s$beta)[[2]], as.character(1:2))
})

test_that("sort_modgirt() orders factors by decreasing sum of squares", {
    s <- sort_modgirt(test_modgirt_id(2))
    ss <- posterior::E(posterior::rvar_apply(s$beta, 2, function(x) {
        posterior::rvar_sum(x^2)
    }))
    expect_true(all(diff(ss) <= 0))
})

test_that("sort_modgirt() works with a single factor", {
    id <- test_modgirt_id(1)
    s <- sort_modgirt(id)
    expect_equal(dim(s$bar_theta), dim(id$bar_theta))
    expect_equal(dim(s$beta), dim(id$beta))
})

test_that("sort_modgirt() preserves labels, matrices, and classes", {
    id <- test_modgirt_id(2)
    s <- sort_modgirt(id)
    expect_s3_class(s, "modgirt_sorted")
    expect_s3_class(s, "modgirt_id")
    expect_equal(attr(s, "item_labels"), attr(id, "item_labels"))
    expect_equal(attr(s, "rotation matrix"), attr(id, "rotation matrix"))
    expect_length(attr(s, "factor order"), 2L)
})

test_that("sign_modgirt() rejects invalid signs", {
    id <- test_modgirt_id(1)
    expect_error(sign_modgirt(id, signs = 2))
    expect_error(sign_modgirt(id, signs = 0))
    expect_error(sign_modgirt(id, signs = c(1, -1, 1)))
})

test_that("sign_modgirt() preserves dimensions with two factors", {
    id <- test_modgirt_id(2)
    g <- sign_modgirt(id, signs = c(1, 1))
    expect_equal(dim(g$bar_theta), dim(id$bar_theta))
    expect_equal(dim(g$beta), dim(id$beta))
})

test_that("sign_modgirt() works with a single factor", {
    id <- test_modgirt_id(1)
    g <- sign_modgirt(id, signs = 1)
    expect_equal(dim(g$bar_theta), dim(id$bar_theta))
})

test_that("sign_modgirt() gives mean loadings the requested sign", {
    id <- test_modgirt_id(2)
    for (target in list(c(1, 1), c(-1, -1), c(1, -1))) {
        g <- sign_modgirt(id, signs = target)
        expect_equal(sign(colMeans(posterior::E(g$beta))), target,
                     ignore_attr = TRUE)
    }
})

test_that("sign_modgirt() preserves labels, matrices, and classes", {
    id <- test_modgirt_id(2)
    g <- sign_modgirt(id, signs = c(1, 1))
    expect_s3_class(g, "modgirt_signed")
    expect_equal(attr(g, "item_labels"), attr(id, "item_labels"))
    expect_length(attr(g, "sign flips"), 2L)
})

## The linear predictor bar_theta %*% t(beta) must be invariant: sorting and
## sign flipping are reparameterizations, not changes to the model. These are
## the regression tests for the dimension-collapse and bar_theta indexing bugs.
nu_modgirt <- function(o, t = 1) {
    posterior::E(posterior::`%**%`(
        o$bar_theta[t, , , drop = TRUE], t(o$beta)
    ))
}

test_that("sort_modgirt() leaves the linear predictor invariant", {
    id <- test_modgirt_id(2)
    s <- sort_modgirt(id)
    for (t in seq_len(dim(id$bar_theta)[1])) {
        expect_equal(nu_modgirt(s, t), nu_modgirt(id, t),
                     tolerance = 1e-8, ignore_attr = TRUE)
    }
})

test_that("sign_modgirt() leaves the linear predictor invariant", {
    id <- test_modgirt_id(2)
    g <- sign_modgirt(id, signs = c(1, -1))
    for (t in seq_len(dim(id$bar_theta)[1])) {
        expect_equal(nu_modgirt(g, t), nu_modgirt(id, t),
                     tolerance = 1e-8, ignore_attr = TRUE)
    }
})

test_that("rotations leave the innovation covariance conformable", {
    id <- test_modgirt_id(2)
    for (o in list(sort_modgirt(id), sign_modgirt(id, signs = c(1, 1)))) {
        expect_equal(dim(o$Omega), c(2L, 2L))
        expect_equal(dim(o$Sigma_theta), c(2L, 2L))
    }
})

test_that("label_modgirt() attaches labels to dimnames", {
    lab <- label_modgirt(test_modgirt_id(2))
    expect_s3_class(lab, "modgirt_lab")
    expect_equal(dimnames(lab$beta)$item, attr(lab, "item_labels"))
    expect_equal(dimnames(lab$bar_theta)$unit, attr(lab, "unit_labels"))
    expect_equal(dimnames(lab$bar_theta)$period, attr(lab, "time_labels"))
})

test_that("label_modgirt() output survives sort and sign", {
    lab <- label_modgirt(test_modgirt_id(2))
    s <- sort_modgirt(lab)
    g <- sign_modgirt(lab, signs = c(1, 1))
    expect_equal(dimnames(s$beta)$item, attr(lab, "item_labels"))
    expect_equal(dimnames(g$bar_theta)$unit, attr(lab, "unit_labels"))
})

### LOO

test_that("fit_modgirt() returns a classed object", {
    expect_s3_class(test_modgirt_fit(), "modgirt_fit")
})

test_that("log_lik_index_modgirt() has one row per observed cell", {
    sdat <- test_modgirt_shaped_data()
    idx <- log_lik_index_modgirt(sdat)
    n_resp <- apply(sdat$SSSS, c(1, 2, 3), sum)
    expect_identical(nrow(idx), sum(n_resp > 0))
    expect_identical(idx$position, seq_len(nrow(idx)))
    expect_true(all(idx$n_responses > 0))
    expect_false(anyNA(idx$item))
    expect_false(anyNA(idx$unit))
    expect_false(anyNA(idx$period))
})

test_that("gen_log_lik = TRUE emits one element per observed cell", {
    f <- test_modgirt_fit(gen_log_lik = TRUE)
    ll <- log_lik_draws(f, group = "observation")
    expect_equal(dim(ll)[3L], nrow(log_lik_index_modgirt(f$stan_data)))
    expect_true(all(is.finite(posterior::draws_of(
        posterior::as_draws_rvars(f$fit$draws("log_lik"))$log_lik))))
})

test_that("gen_log_lik = FALSE omits log_lik", {
    f <- test_modgirt_fit(gen_log_lik = FALSE)
    expect_error(log_lik_draws(f), "gen_log_lik")
})

test_that("log_lik sums to the model's ordinal contribution", {
    ## Every cell's contribution is a weighted sum of log probabilities, so
    ## the total must be negative and finite
    f <- test_modgirt_fit(gen_log_lik = TRUE)
    ll <- log_lik_draws(f, group = "observation")
    tot <- apply(ll, c(1, 2), sum)
    expect_true(all(is.finite(tot)))
    expect_true(all(tot < 0))
})

test_that("regrouping preserves the total log-likelihood", {
    f <- test_modgirt_fit(gen_log_lik = TRUE)
    tot <- function(g) sum(apply(log_lik_draws(f, group = g), c(1, 2), sum))
    expect_equal(tot("observation"), tot("dyad"),   tolerance = 1e-10)
    expect_equal(tot("observation"), tot("unit"),   tolerance = 1e-10)
    expect_equal(tot("observation"), tot("period"), tolerance = 1e-10)
    expect_equal(tot("observation"), tot("item"),   tolerance = 1e-10)
})

test_that("dyad grouping has one column per group-item pair", {
    f <- test_modgirt_fit(gen_log_lik = TRUE)
    idx <- log_lik_index_modgirt(f$stan_data)
    ld <- log_lik_draws(f, group = "dyad")
    expect_equal(dim(ld)[3L], nrow(dplyr::distinct(idx, item, unit)))
})

test_that("item_type grouping is rejected for modgirt", {
    f <- test_modgirt_fit(gen_log_lik = TRUE)
    expect_error(log_lik_draws(f, group = "item_type"))
})

test_that("loo_modgirt() returns a loo object", {
    skip_if_not_installed("loo")
    l <- test_modgirt_loo()
    expect_s3_class(l, "loo")
    expect_equal(attr(l, "dbmm_group"), "dyad")
    expect_false(anyNA(l$diagnostics$n_eff))
})

test_that("loo_influential() works on a modgirt loo object", {
    skip_if_not_installed("loo")
    l <- test_modgirt_loo()
    inf <- loo_influential(l, threshold = -Inf)
    expect_equal(nrow(inf), length(loo::pareto_k_values(l)))
    expect_false(is.unsorted(rev(inf$pareto_k)))
})
