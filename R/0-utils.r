stdize <- function (x) {
    (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
}

create_counts <- function (long_data,
                           unit_var = "UNIT",
                           time_var = "TIME",
                           item_var = "ITEM",
                           value_var = "value",
                           weight_var = NULL) {
    xtab_formula <- reformulate(c(time_var, unit_var, item_var, value_var))
    if (is.null(weight_var)) {
        weight_formula <- NULL
    } else {
        weight_formula <- reformulate(weight_var)
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
    colnames(loading_matrix) <- gsub(
        pattern = "x\\[([0-9]+),([0-9]+)\\]",
        replacement = "LambdaV\\2_\\1",
        x = colnames(loading_matrix)
    )
    return(loading_matrix)
}
make_vm_rvar <- function(loading_draws, n_iter, n_chain, n_factor,
                         method = "varimax", maxit = 1000, randomStarts = 1) {
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
make_sp_rvar <- function(rsp_out, n_iter, n_chain, n_factor) {
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
rvar_whiten_matrix <- rfun(whiten_matrix)


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
