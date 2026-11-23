## --------------------------------------------------- Stan data helpers ----

test_that("the horseshoe sets tau0 from p0 and the total item count", {
    d <- mixfac_shrinkage_data("horseshoe", n_dim = 5, n_items = 100, p0 = 2)
    expect_identical(d$loading_prior, 1L)
    expect_equal(d$tau0, 2 / 3 / sqrt(100))
    expect_identical(d$slab_scale, 2)
    expect_identical(d$slab_df, 4)
})

test_that("a supplied tau0 overrides p0", {
    d <- mixfac_shrinkage_data("horseshoe", n_dim = 5, n_items = 100,
                               p0 = 2, tau0 = 0.3)
    expect_equal(d$tau0, 0.3)
})

test_that("the horseshoe validates p0, tau0, and the slab", {
    expect_error(mixfac_shrinkage_data("horseshoe", 3, 10), "p0")
    expect_error(mixfac_shrinkage_data("horseshoe", 3, 10, p0 = 3), "smaller")
    expect_error(mixfac_shrinkage_data("horseshoe", 3, 10, p0 = 0), "positive")
    expect_error(mixfac_shrinkage_data("horseshoe", 3, 10, tau0 = -1), "tau0")
    expect_error(mixfac_shrinkage_data("horseshoe", 3, 10, p0 = 1,
                                       slab_scale = 0), "slab")
})

test_that("normal and lasso priors ignore p0 and tau0", {
    expect_identical(mixfac_shrinkage_data("normal", 3, 10)$loading_prior, 0L)
    expect_identical(mixfac_shrinkage_data("lasso", 3, 10)$loading_prior, 2L)
    expect_warning(mixfac_shrinkage_data("normal", 3, 10, p0 = 1), "ignored")
})

test_that("fit_mixfac() rejects unknown loading priors before fitting", {
    expect_error(
        fit_mixfac(test_mixfac_shaped_data(), loading_prior = "ridge"),
        "should be one of"
    )
})

## ------------------------------------------------- post-processing ----

test_that("permute_load_scale() follows the signed permutation", {
    ## Two draws, one chain: the first keeps the order, the second swaps the
    ## factors and flips a sign
    scale <- posterior::rvar(rbind(c(1, 0.1), c(-2, 0.2)), nchains = 1)
    P <- array(0, dim = c(2, 2, 2))
    P[1, , ] <- diag(2)
    P[2, , ] <- matrix(c(0, 1, -1, 0), 2, 2)
    sp <- posterior::rvar(P, nchains = 1)
    out <- posterior::draws_of(permute_load_scale(scale, sp, n_chain = 1))
    expect_equal(unname(out[1, ]), c(1, 0.1))
    expect_equal(unname(out[2, ]), c(0.2, 2))
})

test_that("summarize_shrinkage() reports norms, shares, and scales", {
    L <- array(0, dim = c(4, 3, 2))   # draws x items x factors
    L[, , 1] <- 1
    L[, , 2] <- 0.1
    x <- posterior::draws_rvars(
        lambda_binary = posterior::rvar(L),
        load_scale = posterior::rvar(matrix(c(1, 0.05), 4, 2, byrow = TRUE))
    )
    s <- summarize_shrinkage(x)
    expect_identical(s$factor, 1:2)
    expect_equal(s$norm_mean, c(sqrt(3), sqrt(0.03)))
    expect_equal(s$share_mean, c(3, 0.03) / 3.03)
    expect_equal(s$scale_mean, c(1, 0.05))
})

test_that("summarize_shrinkage() works without load_scale", {
    x <- posterior::draws_rvars(
        lambda = posterior::rvar(array(1, dim = c(4, 3, 2)))
    )
    s <- summarize_shrinkage(x)
    expect_false("scale_mean" %in% names(s))
    expect_equal(s$share_mean, c(0.5, 0.5))
})

## ------------------------------------------------------------- fits ----

test_that("the default normal prior adds no shrinkage parameters", {
    ## zero-length parameters produce no draws
    draws_vars <- posterior::variables(test_mixfac_fit()$fit$draws())
    expect_false(any(grepl("^load_scale|^hs_|^lasso_", draws_vars)))
})

for (prior in c("horseshoe", "lasso")) {
    test_that(paste("the", prior, "prior shrinks unsupported factors"), {
        fit <- test_shrinkage_fit(prior)
        expect_identical(attr(fit, "loading_prior"), prior)
        d <- extract_mixfac_draws(fit)
        expect_true("load_scale" %in% posterior::variables(d))
        id <- identify_mixfac(d)
        expect_true("load_scale" %in% names(id))
        s <- summarize_shrinkage(id)
        s <- s[order(-s$norm_mean), ]
        ## One true factor carries nearly all of the loadings
        expect_gt(s$share_mean[1], 0.8)
        expect_lt(max(s$norm_mean[-1]), s$norm_mean[1] / 4)
    })
}

test_that("identification without rotation leaves lp__ unchanged", {
    d <- extract_mixfac_draws(test_shrinkage_fit("horseshoe"))
    id <- identify_mixfac(d)
    expect_equal(posterior::draws_of(id$lp__), posterior::draws_of(d$lp__))
    expect_true(isTRUE(all.equal(
        attr(id, "rotation matrix"), diag(3)
    )))
})

test_that("method = 'none' keeps every linear predictor invariant", {
    d <- extract_mixfac_draws(test_shrinkage_fit("horseshoe"))
    id <- identify_mixfac(d)
    before <- nu_cell(d, "binary", t = 1, j = 3, i = 2)
    after <- nu_cell(id, "binary", t = 1, j = 3, i = 2)
    expect_equal(posterior::draws_of(after), posterior::draws_of(before),
                 tolerance = 1e-8)
})

test_that("rotating a shrinkage fit warns and drops load_scale", {
    d <- extract_mixfac_draws(test_shrinkage_fit("horseshoe"))
    ## GPArotation may not converge on near-zero columns; that is beside the
    ## point here
    quiet_gpa <- function(expr) {
        withCallingHandlers(expr, warning = function(w) {
            if (grepl("Convergence not obtained", conditionMessage(w))) {
                invokeRestart("muffleWarning")
            }
        })
    }
    expect_warning(id <- quiet_gpa(identify_mixfac(d, method = "varimax")),
                   "load_scale")
    expect_false("load_scale" %in% names(id))
})

test_that("load_scale survives label, combine, sort, and sign", {
    d <- extract_mixfac_draws(test_shrinkage_fit("horseshoe"))
    comb <- combine_types_mixfac(label_mixfac(identify_mixfac(d)))
    expect_true("load_scale" %in% names(comb))
    srt <- sort_mixfac(comb)
    fo <- attr(srt, "factor order")
    expect_equal(posterior::draws_of(srt$load_scale),
                 posterior::draws_of(comb$load_scale)[, fo, drop = FALSE],
                 ignore_attr = TRUE)
    sgn <- sign_mixfac(srt)
    expect_equal(posterior::draws_of(sgn$load_scale),
                 posterior::draws_of(srt$load_scale))
    ## sorting puts the supported factor first
    s <- summarize_shrinkage(sgn)
    expect_identical(which.max(s$norm_mean), 1L)
})
