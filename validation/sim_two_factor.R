## Two-factor simulation for mixfac shrinkage priors (design doc, validation
## step 2). Data come from a probit model with two true factors, one strong
## and one weak; fits use more factors than that. A prior passes if both true
## factors keep non-trivial loadings, the extra factors shrink toward zero,
## and the surviving factors recover the true loadings up to rotation about
## as well as an oracle fit (normal prior, true D = 2) does.
##
## Usage, from the package root:
##   Rscript validation/sim_two_factor.R fit <design> <prior> [p0=2] [chains=4]
##       [warmup=1000] [iter=1000] [dim=<fitted D>]
##   Rscript validation/sim_two_factor.R fit <design> normal dim=2  # oracle
##   Rscript validation/sim_two_factor.R report
##
## Designs:
##   doc    200 units, 60 items, |loadings| about 1.0 and 0.5, fitted D = 4
##   mirt   multiIRT's M = 2 setting from test_hs_run.R (40 units, 400 items,
##          loading sd 2 and 1 on the logit scale), rescaled to probit by
##          dividing by 1.7, fitted D = 5
##
## Outputs go to validation/output/sim2_<design>_<prior>[_p0<p0>].rds.

suppressMessages(devtools::load_all(quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE)
out_dir <- "validation/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

designs <- list(
    doc = list(n_units = 200, n_items = 60, n_dim = 4, seed = 2026),
    mirt = list(n_units = 40, n_items = 400, n_dim = 5, seed = 123)
)

## True parameters and long-format binary responses, one period
sim_two_factor <- function(design) {
    d <- designs[[design]]
    withr::with_seed(d$seed, {
        theta <- matrix(stats::rnorm(d$n_units * 2), d$n_units, 2)
        if (design == "doc") {
            sgn <- function() sample(c(-1, 1), d$n_items, replace = TRUE)
            lambda <- cbind(stats::runif(d$n_items, 0.8, 1.2) * sgn(),
                            stats::runif(d$n_items, 0.4, 0.6) * sgn())
            alpha <- stats::rnorm(d$n_items, 0, 0.5)
        } else {
            lambda <- cbind(stats::rnorm(d$n_items, 0, 2),
                            stats::rnorm(d$n_items, 0, 1)) / 1.7
            alpha <- -stats::rnorm(d$n_items, -0.5, 1) / 1.7
        }
        grid <- expand.grid(unit = seq_len(d$n_units), item = seq_len(d$n_items))
        nu <- alpha[grid$item] +
            rowSums(lambda[grid$item, ] * theta[grid$unit, ])
        long <- data.frame(
            unit = sprintf("u%03d", grid$unit),
            period = 1L,
            item = sprintf("q%03d", grid$item),
            value = stats::rbinom(nrow(grid), 1, stats::pnorm(nu))
        )
    })
    rownames(lambda) <- sprintf("q%03d", seq_len(d$n_items))
    rownames(theta) <- sprintf("u%03d", seq_len(d$n_units))
    list(long = long, lambda = lambda, theta = theta, alpha = alpha)
}

## Orthogonal Procrustes rotation of X onto target Y
procrustes <- function(X, Y) {
    s <- svd(crossprod(X, Y))
    s$u %*% t(s$v)
}

## Per-factor norms of loadings, per draw, as [iteration, chain, factor]
factor_norms <- function(lambda_rvar) {
    L <- posterior::draws_of(lambda_rvar, with_chains = TRUE)
    sqrt(apply(L^2, c(1, 2, 4), sum))
}

if (identical(args[1], "fit")) {
    design <- match.arg(args[2], names(designs))
    prior <- match.arg(args[3], c("normal", "horseshoe", "lasso"))
    opt <- list(p0 = "2", chains = "4", warmup = "1000", iter = "1000",
                dim = as.character(designs[[design]]$n_dim))
    for (kv in args[-(1:3)]) {
        kv <- strsplit(kv, "=", fixed = TRUE)[[1]]
        opt[[kv[1]]] <- kv[2]
    }
    sim <- sim_two_factor(design)
    shaped <- suppressMessages(shape_mixfac(
        long_data = sim$long,
        unit_var = "unit",
        time_var = "period",
        item_var = "item",
        value_var = "value",
        make_indicator_for_zeros = FALSE,
        periods_to_estimate = 1L,
        min_cat_n = 0
    ))
    n_dim <- as.integer(opt$dim)
    fit_args <- list(
        shaped_data = shaped,
        n_dim = n_dim,
        loading_prior = prior,
        ## Numeric, not integer: integer chain counts trip check_arg_type()
        chains = as.numeric(opt$chains),
        parallel_chains = as.numeric(opt$chains),
        iter_warmup = as.integer(opt$warmup),
        iter_sampling = as.integer(opt$iter),
        refresh = 250,
        seed = 1
    )
    if (prior == "horseshoe") fit_args$p0 <- as.numeric(opt$p0)
    t0 <- Sys.time()
    fit <- do.call(fit_mixfac, fit_args)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    ## Shrinkage fits fix only signs and order; the normal prior is rotated
    id <- identify_mixfac(extract_mixfac_draws(fit),
                          method = if (prior == "normal") "varimax" else "none")
    shrink <- summarize_shrinkage(id)
    ord <- order(-shrink$norm_mean)

    ## Convergence after identification: R-hat of each factor's loading norm
    ## and of load_scale, in the fitted order
    N <- factor_norms(id$lambda_binary)
    rhat_norm <- apply(N, 3, posterior::rhat)
    rhat_scale <- if (!is.null(id$load_scale)) {
        S <- posterior::draws_of(id$load_scale, with_chains = TRUE)
        apply(S, 3, posterior::rhat)
    }
    ## Per-chain shares, to see whether chains agree on the active factors
    ss <- N^2
    chain_share <- apply(ss / array(apply(ss, c(1, 2), sum), dim(ss)),
                         c(2, 3), mean)

    ## Recovery: posterior-mean loadings and factors of the two largest
    ## fitted factors, rotated onto the truth
    item_lab <- attr(shaped, "binary_item_labels")
    unit_lab <- attr(shaped, "unit_labels")
    L_hat <- apply(posterior::draws_of(id$lambda_binary), c(2, 3), mean)
    E_hat <- apply(posterior::draws_of(id$eta), c(2, 3, 4), mean)[1, , ]
    top <- ord[1:2]
    L_true <- sim$lambda[item_lab, ]
    E_true <- sim$theta[unit_lab, ]
    R <- procrustes(L_hat[, top], L_true)
    L_rot <- L_hat[, top] %*% R
    E_rot <- E_hat[, top] %*% R
    recovery <- data.frame(
        true_factor = 1:2,
        loading_cor_rotated = diag(stats::cor(L_rot, L_true)),
        factor_cor_rotated = diag(stats::cor(E_rot, E_true)),
        ## Without rotation: best absolute correlation with any top factor
        loading_cor_raw = apply(abs(stats::cor(L_hat[, top], L_true)), 2, max),
        true_norm = sqrt(colSums(L_true^2)),
        fitted_norm_rotated = sqrt(colSums(L_rot^2))
    )

    diag <- fit$fit$diagnostic_summary(quiet = TRUE)
    sdiag <- posterior::as_draws_df(fit$fit$sampler_diagnostics())
    res <- list(
        design = design,
        prior = prior,
        options = opt,
        n_dim = n_dim,
        n_units = shaped$J,
        n_items = shaped$I_binary,
        seconds = secs,
        divergences = sum(diag$num_divergent),
        max_treedepth = sum(diag$num_max_treedepth),
        mean_leapfrog = mean(sdiag$n_leapfrog__),
        shrinkage = cbind(shrink[ord, ], rhat_norm = rhat_norm[ord]),
        chain_share = chain_share[, ord, drop = FALSE],
        recovery = recovery
    )
    if (!is.null(rhat_scale)) res$shrinkage$rhat_scale <- rhat_scale[ord]
    tag <- switch(prior, horseshoe = paste0("_p0", opt$p0),
                  normal = paste0("_d", n_dim), "")
    saveRDS(res, file.path(out_dir,
                           paste0("sim2_", design, "_", prior, tag, ".rds")))
    print(res[setdiff(names(res), "options")], digits = 3)
}

if (identical(args[1], "report")) {
    files <- list.files(out_dir, "^sim2_.*\\.rds$", full.names = TRUE)
    for (f in files) {
        r <- readRDS(f)
        s <- r$shrinkage
        k <- nrow(s)
        ## Pass: top two keep >= 95% of squared loadings together, the third
        ## factor's 95% upper bound is below the second's 5% lower bound, and
        ## rotated loadings correlate with the truth no worse than the oracle
        ## (normal prior at the true D = 2) minus 0.02
        oracle <- file.path(out_dir, paste0("sim2_", r$design, "_normal_d2.rds"))
        cut <- if (file.exists(oracle)) {
            readRDS(oracle)$recovery$loading_cor_rotated - 0.02
        } else {
            0.9
        }
        pass <- if (k < 3) NA else sum(s$share_mean[1:2]) >= 0.95 &&
            s$norm_hi[3] < s$norm_lo[2] &&
            all(r$recovery$loading_cor_rotated >= cut)
        cat(sprintf(paste("\n== %s | %s%s | D = %d | %.1f min | %d divergent",
                          "| %.0f steps/draw | %s\n"),
                    r$design, r$prior,
                    if (r$prior == "horseshoe") paste0(" p0 = ", r$options$p0) else "",
                    r$n_dim, r$seconds / 60, r$divergences, r$mean_leapfrog,
                    if (is.na(pass)) "baseline" else if (pass) "PASS" else "FAIL"))
        print(within(s, factor <- NULL), digits = 3, row.names = FALSE)
        cat("per-chain shares (rows = chains):\n")
        print(round(r$chain_share, 3))
        print(r$recovery, digits = 3, row.names = FALSE)
    }
}
