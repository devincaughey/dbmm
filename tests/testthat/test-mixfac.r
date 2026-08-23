## ------------------------------------------------------------- fixture ----

test_that("the test fixture covers all four item types", {
    sdat <- test_shaped_data()
    expect_identical(sdat$I_binary, 2L)
    expect_identical(sdat$I_trichot, 1L)
    expect_identical(sdat$I_ordinal, 1L)
    expect_gt(sdat$I_metric, 20L)
    expect_identical(sdat$K_ordinal, 5L)
    expect_identical(sdat$J, 50L)
    expect_identical(sdat$T, 2L)
})

## --------------------------------------------------------------- shape ----

test_that("shape_mixfac() returns a mixfac_data object with expected structure", {
    sdat <- test_shaped_data()
    expect_s3_class(sdat, "mixfac_data")
    expect_true(all(c("J", "T", "N_binary", "I_binary", "N_trichot",
                      "I_trichot", "N_ordinal", "I_ordinal", "K_ordinal",
                      "N_metric", "I_metric") %in% names(sdat)))
    ## tob_* were removed in P5
    expect_false(any(grepl("^tob_", names(sdat))))
    ## every type's index vectors agree in length with its N
    for (ty in c("binary", "trichot", "ordinal", "metric")) {
        n <- sdat[[paste0("N_", ty)]]
        for (pre in c("yy_", "ii_", "jj_", "tt_")) {
            expect_length(sdat[[paste0(pre, ty)]], n)
        }
    }
    ## indices are within bounds
    expect_true(all(sdat$ii_binary >= 1 & sdat$ii_binary <= sdat$I_binary))
    expect_true(all(sdat$jj_binary >= 1 & sdat$jj_binary <= sdat$J))
    expect_true(all(sdat$tt_binary >= 1 & sdat$tt_binary <= sdat$T))
    ## labels present and consistent with counts
    expect_length(attr(sdat, "unit_labels"), sdat$J)
    expect_length(attr(sdat, "time_labels"), sdat$T)
    expect_length(attr(sdat, "binary_item_labels"), sdat$I_binary)
    expect_length(attr(sdat, "trichotomous_item_labels"), sdat$I_trichot)
    expect_length(attr(sdat, "ordinal_item_labels"), sdat$I_ordinal)
    expect_length(attr(sdat, "metric_item_labels"), sdat$I_metric)
})

test_that("shape_mixfac() codes outcomes as Stan expects", {
    sdat <- test_shaped_data()
    ## binary outcomes are 0/1; ordinal and trichotomous are 1-based
    expect_true(all(sdat$yy_binary %in% c(0L, 1L)))
    expect_identical(min(sdat$yy_trichot), 1L)
    expect_identical(max(sdat$yy_trichot), 3L)
    expect_identical(min(sdat$yy_ordinal), 1L)
    expect_identical(max(sdat$yy_ordinal), sdat$K_ordinal)
    ## metric items are standardized
    expect_lt(abs(mean(sdat$yy_metric)), 0.1)
})

test_that("shape_mixfac() accepts item-type vectors of length > 1", {
    ## Regression: `if (is.na(binary_items))` errored on vectors in R >= 4.2
    long <- test_long_data()
    expect_no_error(suppressMessages(
        shape_mixfac(long, "st", "year", "outcome", "value",
                     binary_items = c("d_binary1", "d_binary2"),
                     make_indicator_for_zeros = FALSE,
                     periods_to_estimate = 2020:2021)
    ))
})

test_that("log_lik_index() rows match the total observation count", {
    sdat <- test_shaped_data()
    idx <- log_lik_index(sdat)
    n_tot <- sum(sdat$N_binary, sdat$N_trichot, sdat$N_ordinal, sdat$N_metric)
    expect_identical(nrow(idx), n_tot)
    expect_identical(idx$position, seq_len(n_tot))
    expect_false(anyNA(idx$item))
    expect_false(anyNA(idx$unit))
    expect_false(anyNA(idx$period))
    ## types appear in the order Stan concatenates them
    expect_identical(unique(idx$item_type),
                     c("binary", "trichot", "ordinal", "metric"))
})

## ----------------------------------------------------------------- fit ----

test_that("fit_mixfac() rejects unsupported links", {
    expect_error(
        fit_mixfac(test_shaped_data(), link = "logit"),
        regexp = "probit"
    )
})

test_that("fit_mixfac() returns a mixfac_fit object carrying labels", {
    out <- test_mixfac_fit()
    expect_s3_class(out, "mixfac_fit")
    expect_true("fit" %in% names(out))
    expect_false(is.null(attr(out$fit, "unit_labels")))
    expect_false(is.null(attr(out$fit, "metric_item_labels")))
})

test_that("evolution parameters are estimated only when identified (P2)", {
    ## smooth_eta = TRUE, T > 1: random walk in eta, so Omega is estimated
    d_rw <- extract_mixfac_draws(test_mixfac_fit(smooth_eta = TRUE), drop = "^z_")
    expect_true("sigma_eta_evol" %in% names(d_rw))
    expect_true("Lcorr_eta" %in% names(d_rw))
    Om <- posterior::E(d_rw$Omega)
    expect_false(isTRUE(all.equal(Om, diag(nrow = nrow(Om)),
                                  check.attributes = FALSE)))

    ## smooth_eta = FALSE: no random walk, so they are absent and Omega = I
    d_sep <- extract_mixfac_draws(test_mixfac_fit(smooth_eta = FALSE), drop = "^z_")
    expect_false("sigma_eta_evol" %in% names(d_sep))
    expect_false("Lcorr_eta" %in% names(d_sep))
    Om_sep <- posterior::E(d_sep$Omega)
    expect_equal(Om_sep, diag(nrow = nrow(Om_sep)), ignore_attr = TRUE)
    expect_lt(max(posterior::sd(d_sep$Omega)), 1e-12)
})

test_that("smooth_eta defaults to TRUE", {
    expect_true(formals(fit_mixfac)$smooth_eta)
})

test_that("deprecated separate_eta maps to the reverse of smooth_eta", {
    skip_if_no_cmdstan()
    sdat <- test_shaped_data()
    expect_warning(
        out <- suppressMessages(fit_mixfac(
            sdat, n_dim = 1, chains = 1, iter_warmup = 50,
            iter_sampling = 50, refresh = 0, seed = 1,
            separate_eta = TRUE
        )),
        regexp = "deprecated"
    )
    ## separate_eta = TRUE means no smoothing, hence no evolution parameters
    d <- extract_mixfac_draws(out, drop = "^z_")
    expect_false("sigma_eta_evol" %in% names(d))
})

test_that("supplying both smooth_eta and separate_eta is an error", {
    expect_error(
        fit_mixfac(test_shaped_data(), smooth_eta = TRUE, separate_eta = FALSE),
        regexp = "not both"
    )
})

test_that("intercept deviates are estimated only when identified (P3)", {
    ## constant_alpha = TRUE: no alpha random walk
    d_const <- extract_mixfac_draws(test_mixfac_fit(constant_alpha = TRUE),
                                    drop = "^nothing$")
    expect_false("sigma_alpha_evol" %in% names(d_const))
    expect_false("z_alpha_trichot" %in% names(d_const))
    expect_false("z_alpha_ordinal" %in% names(d_const))
    ## one period of binary and metric deviates only
    expect_identical(dim(d_const$z_alpha_binary), c(1L, 2L))
    expect_identical(dim(d_const$z_alpha_metric)[1], 1L)
    ## intercepts are nonetheless full [T, I] and constant across periods
    a <- posterior::E(d_const$alpha_binary)
    expect_identical(nrow(a), 2L)
    expect_equal(a[1, ], a[2, ], ignore_attr = TRUE)
    ## trichotomous and ordinal intercepts are pinned to zero
    expect_true(all(abs(posterior::E(d_const$alpha_trichot)) < 1e-12))
    expect_true(all(abs(posterior::E(d_const$alpha_ordinal)) < 1e-12))

    ## constant_alpha = FALSE: drift is estimated
    d_drift <- extract_mixfac_draws(test_mixfac_fit(constant_alpha = FALSE),
                                    drop = "^nothing$")
    expect_true("sigma_alpha_evol" %in% names(d_drift))
    expect_identical(dim(d_drift$z_alpha_binary), c(2L, 2L))
})

test_that("extract_mixfac_draws() drops internal variables by default", {
    d <- extract_mixfac_draws(test_mixfac_fit())
    expect_false(any(grepl("^z_|^r_eta$|^WW$|^L_eta$|^Lcorr|^sigma_eta_evol",
                           names(d))))
    ## r_eta, WW, and r_Omega were removed from the Stan program in P1/P4
    expect_false(any(c("r_eta", "WW", "r_Omega") %in% names(d)))
    ## substantive variables retained
    expect_true(all(c("eta", "Omega", "lp__", "lambda_binary", "lambda_metric",
                      "alpha_binary", "kappa_trichot", "kappa_ordinal",
                      "sigma_metric") %in% names(d)))
    expect_s3_class(d, "mixfac_draws")
})

test_that("draws have the dimensions implied by the data", {
    sdat <- test_shaped_data()
    d <- extract_mixfac_draws(test_mixfac_fit(n_dim = 2))
    expect_identical(dim(d$eta), c(sdat$T, sdat$J, 2L))
    expect_identical(dim(d$lambda_binary), c(sdat$I_binary, 2L))
    expect_identical(dim(d$alpha_metric), c(sdat$T, sdat$I_metric))
    expect_identical(dim(d$kappa_trichot), c(sdat$I_trichot, 2L))
    expect_identical(dim(d$kappa_ordinal),
                     c(sdat$I_ordinal, sdat$K_ordinal - 1L))
    expect_identical(dim(d$Omega), c(2L, 2L))
})

test_that("gen_log_lik = TRUE returns one log_lik element per observation (P6)", {
    sdat <- test_shaped_data()
    d <- extract_mixfac_draws(test_mixfac_fit(gen_log_lik = TRUE))
    expect_true("log_lik" %in% names(d))
    expect_length(d$log_lik, nrow(log_lik_index(sdat)))
    ll <- posterior::E(d$log_lik)
    expect_true(all(is.finite(ll)))
    ## discrete contributions are log probabilities, hence negative
    n_disc <- sum(sdat$N_binary, sdat$N_trichot, sdat$N_ordinal)
    expect_true(all(ll[seq_len(n_disc)] < 0))
})

test_that("gen_log_lik = FALSE omits log_lik", {
    d <- extract_mixfac_draws(test_mixfac_fit(gen_log_lik = FALSE))
    expect_false("log_lik" %in% names(d))
})

## -------------------------------------------------------- process draws ----

test_that("the identify -> label -> combine pipeline runs and preserves labels", {
    ## ESS caps are expected: the test fixture has too few draws for stable
    ## ESS estimates.
        d <- extract_mixfac_draws(test_mixfac_fit())
        id <- identify_mixfac(d, whiten = FALSE)
        expect_s3_class(id, "mixfac_id")
        expect_true(posterior::is_draws_rvars(id))
        expect_false(is.null(attr(id, "time_labels")))
        expect_false(is.null(attr(id, "rotation matrix")))
        expect_false(is.null(attr(id, "signed-permutation matrix")))

        lab <- label_mixfac(id)
        expect_s3_class(lab, "mixfac_lab")
        expect_identical(dimnames(lab$eta)[["unit"]], attr(d, "unit_labels"))
        expect_identical(dimnames(lab$lambda_metric)[["item"]],
                         attr(d, "metric_item_labels"))

        cmb <- combine_types_mixfac(lab)
        expect_s3_class(cmb, "mixfac_comb")
        expect_true(posterior::is_draws_rvars(cmb))
        expect_true(all(c("eta", "lambda", "alpha") %in% names(cmb)))
        expect_false(is.null(attr(cmb, "time_labels")))
        ## combined lambda stacks every item type
        sdat <- test_shaped_data()
        expect_identical(dim(cmb$lambda)[1],
                         sum(sdat$I_binary, sdat$I_trichot,
                             sdat$I_ordinal, sdat$I_metric))

        summ <- summarize_mixfac(
            cmb,
            summary_functions = list(
                mean = ~ posterior::E(.),
                sd   = ~ posterior::sd(.)
            ))
        expect_type(summ, "list")
        expect_true("eta" %in% names(summ))
})

test_that("identification leaves lp__ unchanged", {
    ## Rotation and whitening are reparameterizations, so the posterior of the
    ## log density must be untouched.
    d <- extract_mixfac_draws(test_mixfac_fit())
    id <- identify_mixfac(d, whiten = FALSE)
    expect_equal(posterior::E(id$lp__), posterior::E(d$lp__))
})

test_that("combine_types_mixfac() works on unlabelled input", {
    ## Regression: paste0() on NULL dimnames produced a length-1 name vector
    d <- extract_mixfac_draws(test_mixfac_fit())
    id <- identify_mixfac(d, whiten = FALSE)
    expect_no_error(cmb <- combine_types_mixfac(id))
    expect_true(posterior::is_draws_rvars(cmb))
    expect_true(all(grepl(": ", dimnames(cmb$lambda)[["item"]])))
})

test_that("combine_types_mixfac() errors on already-combined input", {
    d <- extract_mixfac_draws(test_mixfac_fit())
    cmb <- combine_types_mixfac(identify_mixfac(d, whiten = FALSE))
    expect_error(combine_types_mixfac(cmb), regexp = "lambda")
})

test_that("combine_types_mixfac() handles the long format", {
    d <- extract_mixfac_draws(test_mixfac_fit())
    lng <- label_mixfac(identify_mixfac(d, whiten = FALSE), make_long = TRUE)
    cmb <- combine_types_mixfac(lng)
    expect_s3_class(cmb, "mixfac_comb")
    expect_s3_class(cmb$lambda, "data.frame")
    expect_true("item_type" %in% names(cmb$lambda))
    expect_false(is.null(attr(cmb, "time_labels")))
})

test_that("summarize_mixfac() accepts a variable subset", {
    d <- extract_mixfac_draws(test_mixfac_fit())
    cmb <- combine_types_mixfac(label_mixfac(identify_mixfac(d)))
    summ <- summarize_mixfac(cmb[c("eta", "lambda")])
    expect_named(summ, c("eta", "lambda"))
})

test_that("identify_mixfac() agrees whether given a fit or its draws", {
    ## Varimax is identified only up to signed permutation of factors, and
    ## GPFRSorth() uses a random start, so canonicalize before comparing.
    fit <- test_mixfac_fit()
    canon <- function(o) sign_mixfac(sort_mixfac(
        combine_types_mixfac(label_mixfac(o))))
    id_a <- canon(identify_mixfac(fit, whiten = FALSE))
    id_b <- canon(identify_mixfac(extract_mixfac_draws(fit), whiten = FALSE))
    expect_equal(posterior::E(id_a$eta), posterior::E(id_b$eta),
                 tolerance = 1e-6)
    expect_equal(posterior::E(id_a$lambda), posterior::E(id_b$lambda),
                 tolerance = 1e-6)
})

test_that("identify_mixfac() works when Omega is absent", {
    skip_on_cran()
    d <- posterior::as_draws_rvars(test_mixfac_fit()$fit$draws())
    d$Omega <- NULL
    id <- expect_no_error(identify_mixfac(d))
    expect_false("Omega" %in% posterior::variables(id))
})

## ------------------------------------------------------------ whitening ----

test_that("whitening leaves every linear predictor invariant (P7)", {
    d <- extract_mixfac_draws(test_mixfac_fit())
    id_raw <- identify_mixfac(d, whiten = FALSE)
    id_wht <- identify_mixfac(d, whiten = TRUE, ref_t = "mean")
    ty <- richest_item_type(id_raw)
    skip_if(is.na(ty))
    n_item <- dim(id_raw[[paste0("lambda_", ty)]])[1]

    worst <- 0
    for (t in 1:2) for (j in 1:3) for (i in seq_len(min(3L, n_item))) {
        delta <- nu_cell(id_raw, ty, t, j, i) - nu_cell(id_wht, ty, t, j, i)
        worst <- max(worst, max(abs(posterior::draws_of(delta))))
    }
    expect_lt(worst, 1e-8)
})

test_that("whitening actually changes the intercepts and factor scores", {
    ## Guards against the invariance test passing because nothing happened
    d <- extract_mixfac_draws(test_mixfac_fit())
    id_raw <- identify_mixfac(d, whiten = FALSE)
    id_wht <- identify_mixfac(d, whiten = TRUE, ref_t = "mean")
    ty <- richest_item_type(id_raw)
    skip_if(is.na(ty))
    anm <- paste0("alpha_", ty)
    expect_gt(max(abs(posterior::E(id_wht[[anm]]) -
                      posterior::E(id_raw[[anm]]))), 1e-6)
    expect_gt(max(abs(posterior::E(id_wht$eta) -
                      posterior::E(id_raw$eta))), 1e-6)
})

test_that("whitening centres the anchor and retains the thresholds", {
    d <- extract_mixfac_draws(test_mixfac_fit())
    id <- identify_mixfac(d, whiten = TRUE, ref_t = "mean")
    XX <- posterior::rvar_apply(id$eta, c(2, 3), posterior::rvar_mean)
    m <- posterior::E(posterior::rvar_apply(XX, 2, posterior::rvar_mean))
    expect_true(all(abs(m) < 1e-8))
    ## kappa is unchanged by whitening and no longer dropped from output
    expect_true("kappa_ordinal" %in% names(id))
    expect_equal(posterior::E(id$kappa_ordinal),
                 posterior::E(d$kappa_ordinal))
    expect_true("kappa_trichot" %in% names(id))
    expect_equal(posterior::E(id$kappa_trichot),
                 posterior::E(d$kappa_trichot))
})

test_that("whitening accepts other ref_t values", {
    d <- extract_mixfac_draws(test_mixfac_fit())
    expect_no_error(identify_mixfac(d, whiten = TRUE, ref_t = "last"))
    expect_no_error(identify_mixfac(d, whiten = TRUE, ref_t = 1))
    expect_error(identify_mixfac(d, whiten = TRUE, ref_t = 99),
                 regexp = "ref_t")
})

test_that("whitening works with a single factor", {
    d <- extract_mixfac_draws(test_mixfac_fit(n_dim = 1))
    expect_no_error(id <- identify_mixfac(d, whiten = TRUE, ref_t = "mean"))
    expect_identical(dim(id$eta)[3], 1L)
})

## ------------------------------------------------------- sort and sign ----

test_that("sort_mixfac() orders factors by sum of squared loadings", {
    srt <- sort_mixfac(test_comb(n_dim = 2))
    ss <- posterior::E(posterior::rvar_apply(
        srt$lambda, 2, function(y) posterior::rvar_sum(y^2)
    ))
    expect_true(all(diff(ss) <= 0))
    expect_s3_class(srt, "mixfac_sorted")
    expect_s3_class(srt, "mixfac_comb")
    expect_true(posterior::is_draws_rvars(srt))
    expect_length(attr(srt, "factor order"), 2L)
})

test_that("sort_mixfac() and sign_mixfac() retain the thresholds", {
    cmb <- test_comb()
    for (fn in list(sort_mixfac, sign_mixfac)) {
        out <- fn(cmb)
        expect_true(all(c("kappa_trichot", "kappa_ordinal") %in% names(out)))
        expect_equal(posterior::E(out$kappa_trichot),
                     posterior::E(cmb$kappa_trichot))
        expect_equal(posterior::E(out$kappa_ordinal),
                     posterior::E(cmb$kappa_ordinal))
        expect_equal(posterior::E(out$sigma_metric),
                     posterior::E(cmb$sigma_metric))
    }
})

test_that("sort_mixfac() and sign_mixfac() retain the label attributes", {
    cmb <- test_comb()
    for (fn in list(sort_mixfac, sign_mixfac)) {
        out <- fn(cmb)
        expect_identical(attr(out, "unit_labels"), attr(cmb, "unit_labels"))
        expect_identical(attr(out, "time_labels"), attr(cmb, "time_labels"))
        expect_identical(attr(out, "metric_item_labels"),
                         attr(cmb, "metric_item_labels"))
    }
})

test_that("sorting and signing leave every linear predictor invariant", {
    cmb <- test_comb()
    n_item <- dim(cmb$lambda)[1]
    srt <- sort_mixfac(cmb)
    sgn <- sign_mixfac(cmb)
    worst_srt <- worst_sgn <- 0
    for (t in 1:2) for (j in 1:3) for (i in seq_len(min(3L, n_item))) {
        base <- nu_cell_comb(cmb, t, j, i)
        worst_srt <- max(worst_srt, max(abs(posterior::draws_of(
            base - nu_cell_comb(srt, t, j, i)))))
        worst_sgn <- max(worst_sgn, max(abs(posterior::draws_of(
            base - nu_cell_comb(sgn, t, j, i)))))
    }
    expect_lt(worst_srt, 1e-10)
    expect_lt(worst_sgn, 1e-10)
})

test_that("sorting and signing leave lp__ unchanged", {
    cmb <- test_comb()
    expect_equal(posterior::E(sort_mixfac(cmb)$lp__), posterior::E(cmb$lp__))
    expect_equal(posterior::E(sign_mixfac(cmb)$lp__), posterior::E(cmb$lp__))
})

test_that("sign_mixfac() produces the requested loading signs", {
    cmb <- test_comb()
    pos <- sign_mixfac(cmb, signs = 1)
    expect_true(all(colMeans(posterior::E(pos$lambda)) > 0))
    neg <- sign_mixfac(cmb, signs = -1)
    expect_true(all(colMeans(posterior::E(neg$lambda)) < 0))
    ## Omega stays a valid covariance matrix under sign flips
    Om <- posterior::E(neg$Omega)
    expect_equal(Om, t(Om), ignore_attr = TRUE)
    expect_true(all(diag(Om) > 0))
})

test_that("sign_mixfac() validates its signs argument", {
    cmb <- test_comb()
    expect_error(sign_mixfac(cmb, signs = 0), regexp = "signs")
    expect_error(sign_mixfac(cmb, signs = c(1, -1, 1)))
})

test_that("sort and sign work with a single factor", {
    cmb1 <- test_comb(n_dim = 1)
    expect_no_error(srt <- sort_mixfac(cmb1))
    expect_identical(dim(srt$eta)[3], 1L)
    expect_identical(dim(srt$lambda)[2], 1L)
    expect_identical(dim(srt$Omega), c(1L, 1L))
    expect_no_error(sgn <- sign_mixfac(cmb1))
    expect_identical(dim(sgn$eta)[3], 1L)
})

test_that("sort and sign compose", {
    cmb <- test_comb()
    out <- sign_mixfac(sort_mixfac(cmb))
    expect_s3_class(out, "mixfac_signed")
    expect_s3_class(out, "mixfac_sorted")
    expect_true(all(colMeans(posterior::E(out$lambda)) > 0))
    ss <- posterior::E(posterior::rvar_apply(
        out$lambda, 2, function(y) posterior::rvar_sum(y^2)
    ))
    expect_true(all(diff(ss) <= 0))
})

test_that("sort and sign reject long-format input", {
    d <- extract_mixfac_draws(test_mixfac_fit())
    lng <- combine_types_mixfac(
        label_mixfac(identify_mixfac(d), make_long = TRUE)
    )
    expect_error(sort_mixfac(lng), regexp = "rvar")
    expect_error(sign_mixfac(lng), regexp = "rvar")
})

