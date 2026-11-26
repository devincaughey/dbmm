## Validate mixfac shrinkage priors on the 109th Senate (pscl::s109) against
## multiIRT's horseshoe and lasso fits.
##
## Usage, from the package root:
##   Rscript validation/s109_shrinkage.R fit horseshoe   # dbmm fit
##   Rscript validation/s109_shrinkage.R fit lasso       # dbmm fit
##   Rscript validation/s109_shrinkage.R fit horseshoe zeros=TRUE slab=5 \
##       chains=2 warmup=300 iter=300 tag=zeros_slab5    # diagnostic run
##   Rscript validation/s109_shrinkage.R ref <multiIRT senate dir>
##   Rscript validation/s109_shrinkage.R report
##
## Outputs go to validation/output/ (git-ignored). The multiIRT fits use a
## logit link and upper-triangle zeros; dbmm uses probit and no zeros, so
## loading norms differ in scale (logit is about 1.6-1.7 times probit). The
## comparison is the number of factors retained and their shares.

suppressMessages(devtools::load_all(quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE)
out_dir <- "validation/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Votes coded as in multiIRT's test_s109_run.R: 1-3 yea, 4-6 nay, else NA
s109_long <- function() {
    data("s109", package = "pscl", envir = environment())
    v <- s109$votes
    v[v == 0 | v > 6] <- NA
    v[v %in% 1:3] <- 1
    v[v %in% 4:6] <- 0
    long <- data.frame(
        senator = rep(rownames(v), ncol(v)),
        period = 1L,
        vote = rep(sprintf("v%03d", seq_len(ncol(v))), each = nrow(v)),
        value = as.vector(v)
    )
    long[!is.na(long$value), ]
}

## Per-factor summaries of a draws array [draws, items, factors]
factor_summary <- function(L) {
    ss <- apply(L^2, c(1, 3), sum)
    norm <- sqrt(ss)
    data.frame(
        norm_median = apply(norm, 2, stats::median),
        share_mean = colMeans(ss / rowSums(ss))
    )
}

if (identical(args[1], "fit")) {
    prior <- match.arg(args[2], c("horseshoe", "lasso"))
    shaped <- shape_mixfac(
        long_data = s109_long(),
        unit_var = "senator",
        time_var = "period",
        item_var = "vote",
        value_var = "value",
        make_indicator_for_zeros = FALSE,
        periods_to_estimate = 1L,
        min_cat_n = 0
    )
    ## Optional key=value settings, e.g. zeros=TRUE slab=5 chains=2
    ## warmup=300 iter=300 tag=diag
    opt <- list(zeros = "FALSE", slab = "2", chains = "4", warmup = "1000",
                iter = "1000", tag = "")
    for (kv in args[-(1:2)]) {
        kv <- strsplit(kv, "=", fixed = TRUE)[[1]]
        opt[[kv[1]]] <- kv[2]
    }
    fit_args <- list(
        shaped_data = shaped,
        n_dim = 5,
        loading_prior = prior,
        lambda_zeros = if (as.logical(opt$zeros)) TRUE else NULL,
        slab_scale = as.numeric(opt$slab),
        chains = as.numeric(opt$chains),
        parallel_chains = as.numeric(opt$chains),
        iter_warmup = as.integer(opt$warmup),
        iter_sampling = as.integer(opt$iter),
        refresh = 50,
        seed = 1
    )
    if (prior == "horseshoe") fit_args$p0 <- 1
    t0 <- Sys.time()
    fit <- do.call(fit_mixfac, fit_args)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    id <- identify_mixfac(extract_mixfac_draws(fit))
    diag <- fit$fit$diagnostic_summary(quiet = TRUE)
    sdiag <- posterior::as_draws_df(fit$fit$sampler_diagnostics())
    res <- list(
        prior = prior,
        options = opt,
        n_items = shaped$I_binary,
        n_units = shaped$J,
        n_obs = shaped$N_binary,
        seconds = secs,
        divergences = sum(diag$num_divergent),
        max_treedepth = sum(diag$num_max_treedepth),
        mean_leapfrog = mean(sdiag$n_leapfrog__),
        mean_treedepth = mean(sdiag$treedepth__),
        shrinkage = summarize_shrinkage(id),
        load_scale_raw = fit$fit$summary(variables = "load_scale")
    )
    tag <- if (nzchar(opt$tag)) paste0("_", opt$tag) else ""
    saveRDS(res, file.path(out_dir, paste0("dbmm_", prior, tag, ".rds")))
    print(res[setdiff(names(res), "load_scale_raw")])
}

if (identical(args[1], "ref")) {
    dir <- args[2]
    files <- c(
        horseshoe = "fit_hs_3_d5_p1_warm1k_iter1k_seed%d.rds",
        lasso = "fit_bl2_d5_warm1k_iter1k_seed%d.rds"
    )
    ref <- list()
    for (prior in names(files)) {
        rows <- list()
        for (seed in 1:10) {
            f <- file.path(dir, sprintf(files[[prior]], seed))
            if (!file.exists(f)) next
            fit <- readRDS(f)
            ## Raw (unrotated) beta: the upper-triangle zeros fix the rotation
            ## Summarized per chain, since chains can settle in different
            ## modes; pooling them would blur the per-factor norms
            LC <- posterior::draws_of(
                posterior::as_draws_rvars(fit$draws("beta"))$beta,
                with_chains = TRUE
            )
            for (ch in seq_len(dim(LC)[2])) {
                L <- LC[, ch, , , drop = FALSE]
                dim(L) <- dim(L)[-2]
                s <- factor_summary(L)
                s <- s[order(-s$norm_median), ]
                s$rank <- seq_len(nrow(s))
                s$seed <- seed
                s$chain <- ch
                rows[[length(rows) + 1]] <- s
            }
            rm(fit, L, LC); gc()
        }
        ref[[prior]] <- do.call(rbind, rows)
    }
    saveRDS(ref, file.path(out_dir, "multiirt_ref.rds"))
    print(ref)
}

if (identical(args[1], "report")) {
    ref <- readRDS(file.path(out_dir, "multiirt_ref.rds"))
    for (prior in c("horseshoe", "lasso")) {
        f <- file.path(out_dir, paste0("dbmm_", prior, ".rds"))
        if (!file.exists(f)) next
        d <- readRDS(f)
        s <- d$shrinkage[order(-d$shrinkage$norm_mean), ]
        r <- stats::aggregate(cbind(norm_median, share_mean) ~ rank,
                              data = ref[[prior]], FUN = stats::median)
        cat("\n==", prior, "| dbmm:", d$n_items, "items,", d$divergences,
            "divergences,", round(d$seconds / 60, 1), "min\n")
        print(data.frame(
            rank = seq_len(nrow(s)),
            dbmm_norm = round(s$norm_mean, 2),
            dbmm_share = round(s$share_mean, 3),
            multiirt_norm = round(r$norm_median, 2),
            multiirt_share = round(r$share_mean, 3)
        ), row.names = FALSE)
    }
}
