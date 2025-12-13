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
