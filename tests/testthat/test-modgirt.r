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
