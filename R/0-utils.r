stdize <- function (x) {
    (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
}

create_counts <- function (long_data,
                           unit_var = "UNIT",
                           time_var = "TIME",
                           item_var = "ITEM",
                           value_var = "value",
                           weight_var = NULL) {
    xtab_formula <- stats::reformulate(c(time_var, unit_var, item_var, value_var))
    if (is.null(weight_var)) {
        weight_formula <- NULL
    } else {
        weight_formula <- stats::reformulate(weight_var)
    }
    des <- survey::svydesign(~1, data = long_data, weights = weight_formula)
    xtab <- survey::svytable(formula = xtab_formula, design = des) 
    return(xtab)
}
# a more than usually informative error message for handing in the
# wrong type to a function
# from keyATM package
check_arg_type <- function(arg, typename, message = NULL) {
  argname <- deparse(match.call()[['arg']])
  if (!inherits(arg, typename)) {
    if (is.null(message)) {
      cli::cli_abort(paste0('`', argname, '` is not a ', typename, ' object.'))
    } else {
      cli::cli_abort(message)
    }
  }
}

rename_loading_matrix <- function(loading_matrix) {
    nm <- colnames(loading_matrix)
    new <- gsub(
        pattern = "^x\\[([0-9]+),([0-9]+)\\]$",
        replacement = "LambdaV\\2_\\1",
        x = nm
    )
    if (any(new == nm)) {
        cli::cli_abort(c(
            "Cannot rename loading columns for
             {.fn factor.switching::rsp_exact}.",
            "x" = "Unmatched name{?s}: {.val {utils::head(nm[new == nm], 3)}}",
            "i" = "Column names must look like {.val x[i,d]} with numeric
                   indices. Strip {.fn dimnames} from the loadings before
                   flattening them."
        ))
    }
    colnames(loading_matrix) <- new
    loading_matrix
}

make_vm_rvar <- function(loading_draws, n_iter, n_chain, n_factor,
                         method = "varimax", maxit = 1000, randomStarts = 0) {
    ## `loading_draws` should be a `draws_of` of an `draws_rvar` object
    rotmat_array <- array(dim = c(n_iter, n_chain, n_factor, n_factor))
    for (i in seq_len(n_iter)) {
        for (c in seq_len(n_chain)) {
            vm <- GPArotation::GPFRSorth(loading_draws[i, c, , , drop = TRUE],
                                         method = method,
                                         normalize = TRUE,
                                         maxit = maxit,
                                         randomStarts = randomStarts)
            rotmat_array[i, c, , ] <- vm$Th
        }
    }
    rotmat_rvar <- posterior::rvar(rotmat_array, with_chains = TRUE)
    return(rotmat_rvar)
}
make_sp_rvar <- function (rsp_out, n_iter, n_chain, n_factor) {
    stopifnot(nrow(rsp_out$sign_vectors) == n_iter * n_chain)
    sp_array <- array(dim = c(n_iter, n_chain, n_factor, n_factor))
    for (i in seq_len(n_iter)) {
        for (c in seq_len(n_chain)) {
            ## Assumes vectors are ordered by chain then iteration
            sv <- rsp_out$sign_vectors[(c - 1) * n_iter + i, , drop = TRUE]
            pv <- rsp_out$permute_vectors[(c - 1) * n_iter + i, , drop = FALSE]
            sp_array[i, c, , ] <-
                t(
                    diag(sv, nrow = length(sv)) %*%
                    seriation::permutation_vector2matrix(pv)
                )
        }
    }
    sp_rvar <- posterior::rvar(sp_array, with_chains = TRUE)
    return(sp_rvar)
}

#' Maximum-leading-eigenvalue rotation criterion
#'
#' Rotation criterion for use with \pkg{GPArotation}, via
#' `identify_mixfac(method = "maxvar1")`. Maximizes the leading eigenvalue of
#' `crossprod(L)`, concentrating variance in a single dominant factor.
#'
#' @param L (matrix) A loadings matrix.
#' @param ... Ignored; present for compatibility with \pkg{GPArotation}.
#'
#' @return A list with elements `f` (criterion value), `Gq` (gradient), and
#'     `Method`.
#'
#' @keywords internal
#' @export
vgQ.maxvar1 <- function(L, ...) {
  L <- as.matrix(L)
  ## t(L) %*% L
  M <- crossprod(L)
  eig <- eigen(M, symmetric = TRUE)
  v <- eig$vectors[, 1]      # leading eigenvector
  lambda1 <- eig$values[1]   # leading eigenvalue
  ## objective: negative leading eigenvalue
  f <- -lambda1
  ## gradient: -2 * (L v) v'
  Gq <- -2 * (L %*% v) %*% t(v)
  list(f = f, Gq = Gq, Method = "maxvar1")
}
whiten_matrix <- function(DM, eps = 0) {
  ## DM must be de-meaned (observations × variables)
  if (!is.matrix(DM)) DM <- as.matrix(DM)
  n <- nrow(DM)
  p <- ncol(DM)
  
  if (n <= 1) stop("DM must have at least 2 rows")

  ## sample covariance
  SS <- crossprod(DM) / (n - 1)

  ## optional jitter
  if (eps > 0) SS <- SS + diag(eps, p)

  ## precision = inverse covariance
  PP <- solve(SS)

  ## whitening matrix
  WW <- t(chol(PP))

  WW
}

#' Whitening matrix for an `rvar` data matrix
#'
#' Computes the whitening matrix for each posterior draw of `DM`.
#'
#' @param DM (rvar) A de-meaned observations-by-variables random matrix.
#' @param eps (non-negative real) Optional ridge added to the diagonal of the
#'     covariance matrix before inversion.
#' @return An `rvar` of dimension `c(p, p)` with the same draws and chains as
#'     `DM`.
#' @keywords internal
#' @noRd
rvar_whiten_matrix <- function(DM, eps = 0) {
    stopifnot(posterior::is_rvar(DM), length(dim(DM)) == 2L)
    a <- posterior::draws_of(DM)                 # S x n x p
    S <- dim(a)[1L]
    p <- dim(a)[3L]
    out <- array(NA_real_, dim = c(S, p, p))
    for (s in seq_len(S)) {
        out[s, , ] <- whiten_matrix(a[s, , , drop = TRUE], eps = eps)
    }
    dimnames(out) <- list(NULL, dimnames(a)[[3L]], dimnames(a)[[3L]])
    posterior::rvar(out, nchains = posterior::nchains(DM))
}



#' Build a draws_rvars object, dropping absent variables
#'
#' Wrapper around [posterior::draws_rvars()] that silently drops `NULL` and
#' zero-length arguments. Needed because some Stan variables are absent from
#' the draws under certain flag settings (e.g. `sigma_alpha_evol` when
#' `constant_alpha == 1`, `sigma_eta_evol` when `T == 1`), and because item
#' types may be empty (e.g. `kappa_trichot` when `I_trichot == 0`).
#'
#' @param ... Named arguments passed to [posterior::draws_rvars()].
#'
#' @return A `draws_rvars` object containing only the non-empty arguments.
#'
#' @import posterior
#'
#' @keywords internal
as_draws_rvars_safe <- function(...) {
    args <- list(...)
    keep <- vapply(args, function(v) !is.null(v) && length(v) > 0, logical(1))
    do.call(posterior::draws_rvars, args[keep])
}

#' Map positions of `log_lik` to observations
#'
#' @param shaped_data A `mixfac_data` object from [shape_mixfac()].
#'
#' @return A data frame with one row per element of `log_lik`, giving the
#'     item type, item, unit, period, and response value.
#'
#' @export
log_lik_index <- function(shaped_data) {
    check_arg_type(arg = shaped_data, typename = "mixfac_data")
    types <- c("binary", "trichot", "ordinal", "metric")
    label_attrs <- c("binary_item_labels", "trichotomous_item_labels",
                     "ordinal_item_labels", "metric_item_labels")
    out <- vector("list", length(types))
    for (k in seq_along(types)) {
        n <- shaped_data[[paste0("N_", types[k])]]
        if (is.null(n) || n == 0) next
        out[[k]] <- data.frame(
            item_type = types[k],
            item = attr(shaped_data, label_attrs[k])[
                shaped_data[[paste0("ii_", types[k])]]],
            unit = attr(shaped_data, "unit_labels")[
                shaped_data[[paste0("jj_", types[k])]]],
            period = attr(shaped_data, "time_labels")[
                shaped_data[[paste0("tt_", types[k])]]],
            value = shaped_data[[paste0("yy_", types[k])]],
            stringsAsFactors = FALSE
        )
    }
    out <- do.call(rbind, Filter(Negate(is.null), out))
    out$position <- seq_len(nrow(out))
    out
}
#' Copy dbmm label attributes from one object to another
#'
#' Absent attributes are skipped rather than assigned `NULL`, so that "never
#' set" is distinguishable from "set to nothing".
#'
#' @keywords internal
copy_dbmm_attrs <- function(to, from) {
    nms <- c("unit_labels", "time_labels", "item_labels",
             "binary_item_labels", "trichotomous_item_labels",
             "ordinal_item_labels", "metric_item_labels")
    for (nm in nms) {
        a <- attr(from, nm)
        if (!is.null(a)) attr(to, nm) <- a
    }
    to
}

#' Right-multiply each period's factor-score matrix by a D-by-D matrix
#'
#' `eta` is `[T, J, D]`. Indexing with `drop = TRUE` collapses the factor
#' dimension when `D == 1`, so the slice is forced back to `J`-by-1.
#'
#' @keywords internal
rotate_eta_rvar <- function(eta, mat) {
    n_time <- dim(eta)[1]
    n_factor <- dim(eta)[3]
    out <- eta
    for (t in seq_len(n_time)) {
        E_t <- eta[t, , , drop = TRUE]
        if (n_factor == 1) E_t <- t(t(E_t))
        out[t, , ] <- posterior::`%**%`(E_t, mat)
    }
    out
}

check_flag <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
        cli::cli_abort(
            "{.arg {arg}} must be {.code TRUE} or {.code FALSE}, not {.obj_type_friendly {x}}.",
            call = call
        )
    }
    as.logical(x)
}

rvar_solve <- function(A) {
    a <- posterior::draws_of(A)
    for (s in seq_len(dim(a)[1L])) a[s, , ] <- solve(a[s, , , drop = TRUE])
    posterior::rvar(a, nchains = posterior::nchains(A))
}

#' Evaluate code with a temporary RNG seed
#'
#' @param seed (integer or `NULL`) Seed to set. If `NULL`, `code` is evaluated
#'     with the ambient RNG state untouched.
#' @param code Code to evaluate.
#' @keywords internal
with_seed <- function(seed, code) {
    if (is.null(seed)) return(code)
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
        cli::cli_abort("{.arg seed} must be a single integer or {.code NULL}.")
    }
    has_old <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    if (has_old) {
        old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
        on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    } else {
        on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())),
                add = TRUE)
    }
    set.seed(as.integer(seed))
    code
}

#' Renumber factor dimnames after permuting factors
#'
#' Factor labels are positional, so a permutation should renumber them rather
#' than carry the old labels along.
#'
#' @param x (rvar) A variable with one or more factor dimensions.
#' @param dims (integer vector) Which dimensions index factors.
#' @param n_factor (positive integer) Number of factors.
#' @keywords internal
renumber_factor_dimnames <- function(x, dims, n_factor) {
    if (is.null(x) || is.null(dimnames(x))) return(x)
    for (d in dims) dimnames(x)[[d]] <- seq_len(n_factor)
    x
}

#' Classify items by measurement type
#'
#' @param n_unique (named integer vector) Number of distinct observed values
#'     per item.
#' @param max_cats (positive integer) Largest number of distinct values for
#'     which an item is classified as ordinal.
#' @return A named character vector of `"dropped"`, `"binary"`,
#'     `"trichotomous"`, `"ordinal"`, or `"metric"`.
#' @keywords internal
#' @noRd
classify_mixfac_items <- function(n_unique, max_cats = 5) {
    out <- rep(NA_character_, length(n_unique))
    names(out) <- names(n_unique)
    out[n_unique <  2L]                        <- "dropped"
    out[n_unique == 2L]                        <- "binary"
    out[n_unique == 3L]                        <- "trichotomous"
    out[n_unique >  3L & n_unique <= max_cats] <- "ordinal"
    out[n_unique >  max_cats]                  <- "metric"
    out
}
