# when using dynIRTtest::identify
# Warning message:
# Returning more (or less) than 1 row per `summarise()` group was deprecated in dplyr 1.1.0.
# ℹ Please use `reframe()` instead.
# ℹ When switching from `summarise()` to `reframe()`, remember that `reframe()` always returns an ungrouped data frame and adjust accordingly.
# ℹ The deprecated feature was likely used in the dynIRTtest package.
#   Please report the issue to the authors.


#' Extract object from a stan output
#'
#' @param fit
#' @param drop_rex
#' @param format
#'
#' @return
#'
#' @import magrittr
#' @import cmdstanr
#'
#' @export
extract_draws <- function (fit, drop_rex = "^z_", format = "df")
{
  draws <- fit$draws(format = format) %>%
    dplyr::select(-dplyr::matches(drop_rex))

  attr(draws, "unit_labels") <- attr(fit, "unit_labels")
  attr(draws, "time_labels") <- attr(fit, "time_labels")
  attr(draws, "binary_item_labels") <- attr(fit, "binary_item_labels")
  attr(draws, "ordinal_item_labels") <- attr(fit, "ordinal_item_labels")
  attr(draws, "metric_item_labels") <- attr(fit, "metric_item_labels")
  return(draws)
}


#' Identify the sign and rotation of the output
#'
#' @param raw_draws
#' @param rotate
#' @param varimax
#' @param normalize
#' @param id_with
#'
#' @return
#'
#' @import magrittr
#'
#' @export
identify_draws <- function(raw_draws, rotate = FALSE, varimax = TRUE,
                          normalize = TRUE, id_with = NULL, sign = NULL)
{
  if (rotate) {
    if (is.null(sign)) {
      sign <- -1
      cat("Using `sign`=", sign, "\n")
    }
    outcomes_id <- identify_rotation(raw_draws, varimax = varimax,
                                     normalize = normalize, id_with = id_with)
    outcomes_id$id_draws <-
      identify_sign(outcomes_id$id_draws, sign = sign)$id_draws
  } else {
    if (is.null(sign)) {
      sign <- c(-1, -1)
      cat("Using `sign`=", sign, "\n")
    }
    outcomes_id <- identify_sign(raw_draws, sign = sign)

  }
  return(outcomes_id)
}


#'
#' @param raw_draws
#' @param varimax
#' @param normalize
#' @param id_with
#'
#' @return
#'
#' @import magrittr
#'
identify_rotation <- function (raw_draws, varimax,
                               normalize, id_with)
{

  # TODO: create a function
  Lb0 <- extract_draws_match(raw_draws = raw_draws,
                            regex_pars = "(^lambda_binary\\[)|(^\\.)")
  Lo0 <- extract_draws_match(raw_draws = raw_draws,
                            regex_pars = "(^lambda_ordinal\\[)|(^\\.)")
  Lm0 <- extract_draws_match(raw_draws = raw_draws,
                            regex_pars = "(^lambda_metric\\[)|(^\\.)")
  E0 <- extract_draws_match(raw_draws = raw_draws,
                            regex_pars = "(^eta\\[)|(^\\.)")
  S0 <- extract_draws_match(raw_draws = raw_draws,
                            regex_pars = "(^sigma_eta_evol\\[)|(^\\.)")

  # Lb0 <- dplyr::select(raw_draws,
  #                     dplyr::matches("(^lambda_binary\\[)|(^\\.)")) %>%
  #   as.data.frame()
  # Lo0 <- dplyr::select(raw_draws,
  #                     dplyr::matches("(^lambda_ordinal\\[)|(^\\.)")) %>%
  #   as.data.frame()
  # Lm0 <- dplyr::select(raw_draws,
  #                     dplyr::matches("(^lambda_metric\\[)|(^\\.)")) %>%
  #   as.data.frame()
  # E0 <- dplyr::select(raw_draws,
  #                     dplyr::matches("(^eta\\[)|(^\\.)")) %>%
  #   as.data.frame()
  # S0 <- dplyr::select(raw_draws,
  #                     dplyr::matches("(^sigma_eta_evol\\[)|(^\\.)")) %>%
  #   as.data.frame()

  S <- max(raw_draws$.iteration)
  C <- max(raw_draws$.chain)
  Ib <- length(attr(raw_draws, "binary_item_labels"))
  Io <- length(attr(raw_draws, "ordinal_item_labels"))
  Im <- length(attr(raw_draws, "metric_item_labels"))
  J <- length(attr(raw_draws, "unit_labels"))
  T <- length(attr(raw_draws, "time_labels"))
  D <- ncol(S0) - 3
  stopifnot(D > 1)
  Q0 <- array(NA, dim = c(S, C, D, D))

  for (s in 1:S) {
    for (c in 1:C) {
      Q0[s, c, 1:D, 1:D] <- diag(rep(1, D))
    }
  }
  Q1 <- Q0
  Lb1 <- Lb0
  Lo1 <- Lo0
  Lm1 <- Lm0
  lcols_b <- stringr::str_which(names(Lb0), "lambda")
  lcols_o <- stringr::str_which(names(Lo0), "lambda")
  lcols_m <- stringr::str_which(names(Lm0), "lambda")
  ecols <- stringr::str_which(names(E0), "eta")
  ecols_t <- stringr::str_split(names(E0)[ecols], "[\\[,\\]]",
                                simplify = TRUE)[, 2]

  if (is.null(id_with)) {
    id_with <- c("binary", "ordinal", "metric")[which.max(c(Ib, Io, Im))]
  }
  stopifnot(id_with %in% c("binary", "ordinal", "metric"))
  cat("\nUsing", id_with, "items to identify the model.\n")

  if (varimax & D > 1) {
    ## 1. Apply varimax rotation
    cat("\n**** Applying varimax rotations...")
    for (s in 1:S) {
      for (c in 1:C) {
        if (id_with == "binary") {
          row <- Lb0$.chain == c & Lb0$.iteration == s
          Lb0_cs <- matrix(
            as.numeric(Lb0[row, lcols_b]), nrow = Ib, ncol = D
          )
          vm <- stats::varimax(Lb0_cs, normalize = normalize)
          Lb1[row, lcols_b] <- t(as.numeric(vm$loadings))
        }
        else if (id_with == "ordinal") {
          row <- Lo0$.chain == c & Lo0$.iteration == s
          Lo0_cs <- matrix(
            as.numeric(Lo0[row, lcols_o]), nrow = Io, ncol = D
          )
          vm <- varimax(Lo0_cs, normalize = normalize)
          Lo1[row, lcols_o] <- t(as.numeric(vm$loadings))
        }
        else if (id_with == "metric") {
          row <- Lm0$.chain == c & Lm0$.iteration == s
          Lm0_cs <- matrix(
            as.numeric(Lm0[row, lcols_m]), nrow = Im, ncol = D
          )
          vm <- varimax(Lm0_cs, normalize = normalize)
          Lm1[row, lcols_m] <- t(as.numeric(vm$loadings))
        }
        else {
          stop("Invalid")
        }
        Q1[s, c, 1:D, 1:D] <- vm$rotmat
      }
    }
    cat("done.\n")
  }

  ## 2. Sign-permute by chain
  cat("**** Applying sign-permute rotations...\n")
  Q2 <- Q0
  sp_out <- vector("list", length = C)
  for (c_cur in 1:C) {
    if (id_with == "binary") {
      sp_out[[c_cur]] <- sign_permute(lambda_item = Lb1,
                                      c = c_cur,
                                      lcols = lcols_b,
                                      id_with = "binary")
    } else if (id_with == "ordinal") {
      sp_out[[c_cur]] <- sign_permute(lambda_item = Lo1,
                                      c = c_cur,
                                      lcols = lcols_o,
                                      id_with = "ordinal")
    } else if (id_with == "metric") {
      sp_out[[c_cur]] <- sign_permute(lambda_item = Lm1,
                                      c = c_cur,
                                      lcols = lcols_m,
                                      id_with = "metric")
    } else {
      stop("Invalid")
    }
    sv <- sp_out[[c_cur]]$sign_vectors
    pv <- sp_out[[c_cur]]$permute_vectors
    sm <- apply(sv, 1, function (row) data.frame(diag(row)))
    pm <- apply(pv, 1, function (row)
      data.frame(seriation::permutation_vector2matrix(row)))

    for (s in 1:S) {
      Q2[s, c_cur, 1:D, 1:D] <- t(as.matrix(sm[[s]]) %*% as.matrix(pm[[s]]))
    }
  }
  Q3 <- Q0

  cat("**** Harmonizing rotations across chains...\n")

  hv_out <- harmonize_varimax(sp_out)

  for (c_cur in 1:C) {
    (sv <- hv_out$sign_vectors[c_cur, , drop = FALSE])
    (pv <- hv_out$permute_vectors[c_cur, , drop = FALSE])
    (sm <- diag(as.numeric(sv)))
    (pm <- seriation::permutation_vector2matrix(pv))
    for (s in 1:S) {
      Q3[s, c_cur, 1:D, 1:D] <- t(as.matrix(sm) %*% as.matrix(pm))
    }
  }
  Lb3 <- Lb0
  Lo3 <- Lo0
  Lm3 <- Lm0
  E3 <- E0
  S3 <- S0
  for (c_cur in 1:C) {
    cat("\n** Rotating chain", c_cur, "...")
    for (s in 1:S) {
      row <- E0$.chain == c_cur & E0$.iteration == s
      if (Ib > 0) {
        Lb0_cs <- matrix(unlist(Lb0[row, lcols_b]), nrow = Ib, ncol = D)
        Lb3[row, lcols_b] <-
          Lb0_cs %*%
          Q1[s, c_cur, 1:D, 1:D] %*%
          Q2[s, c_cur, 1:D, 1:D] %*%
          Q3[s, c_cur, 1:D, 1:D] %>%
          as.numeric() %>%
          t()
      }
      if (Io > 0) {
        Lo0_cs <- matrix(unlist(Lo0[row, lcols_o]), nrow = Io, ncol = D)
        Lo3[row, lcols_o] <-
          Lo0_cs %*%
          Q1[s, c_cur, 1:D, 1:D] %*%
          Q2[s, c_cur, 1:D, 1:D] %*%
          Q3[s, c_cur, 1:D, 1:D] %>%
          as.numeric() %>%
          t()
      }
      if (Im > 0) {
        Lm0_cs <- matrix(unlist(Lm0[row, lcols_m]), nrow = Im, ncol = D)
        Lm3[row, lcols_m] <-
          Lm0_cs %*%
          Q1[s, c_cur, 1:D, 1:D] %*%
          Q2[s, c_cur, 1:D, 1:D] %*%
          Q3[s, c_cur, 1:D, 1:D] %>%
          as.numeric() %>%
          t()
        S0_cs <- matrix(unlist(S0[row, 1:D, drop = FALSE]))
        S3[row, 1:D] <-    # might need to transpose for conformability
          t(S0_cs) %*%
          abs(Q2[s, c_cur, 1:D, 1:D]) %*%
          abs(Q3[s, c_cur, 1:D, 1:D])
      }
      for (t in 1:T) {
        E0_cst <- matrix(
          unlist(E3[row, ecols[as.integer(ecols_t) == t]]),
          nrow = J, ncol = D
        )
        E3[row, ecols[as.integer(ecols_t) == t]] <-
          E0_cst %*%
          Q1[s, c_cur, 1:D, 1:D] %*%
          Q2[s, c_cur, 1:D, 1:D] %*%
          Q3[s, c_cur, 1:D, 1:D] %>%
          as.numeric() %>%
          t()
      }
    }
  }
  id_draws <- raw_draws
  id_draws[, stringr::str_detect(names(id_draws), "^lambda_binary\\[")] <-
    Lb3[, lcols_b]
  id_draws[, stringr::str_detect(names(id_draws), "^lambda_ordinal\\[")] <-
    Lo3[, lcols_o]
  id_draws[, stringr::str_detect(names(id_draws), "^lambda_metric\\[")] <-
    Lm3[, lcols_m]
  id_draws[, stringr::str_detect(names(id_draws), "^eta\\[")] <- E3[, ecols]
  id_draws[, stringr::str_detect(names(id_draws), "^sigma_eta_evol\\[")] <-
    S3[, 1:D]
  result <- list(
    id_draws = id_draws,
    rotmats = list()
  )
  result$rotmats$Q <- Q0
  result$rotmats$Q1 <- Q1
  result$rotmats$Q2 <- Q2
  result$rotmats$Q3 <- Q3
  for (s in 1:S) {
    for (c_cur in 1:C) {
      result$rotmats$Q[s, c_cur, 1:D, 1:D] <-
        Q1[s, c_cur, 1:D, 1:D] %*%
        Q2[s, c_cur, 1:D, 1:D] %*%
        Q3[s, c_cur, 1:D, 1:D]
    }
  }
  return(result)
}


#' @import magrittr
extract_draws_match <- function(raw_draws, regex_pars) {
  dplyr::select(raw_draws, dplyr::matches(regex_pars)) %>%
    as.data.frame()
}


#' @param lambda_item lambda object
#' @param c_cur current row indicator
#' @param lcols column
#' @param id_with
#'
#' @return sign-permuted object
#'
sign_permute <- function(lambda_item, c_cur, lcols,
                        id_with = c("binary", "ordinal", "metric"))
{
  rows <- lambda_item$.chain == c_cur

  L_c <- as.matrix(lambda_item[rows, lcols])

  var <- as.integer(
    stringr::str_replace(colnames(L_c),
                        ".+\\[([0-9]+),([0-9]+)\\]$", "\\1")
    )

  dim <- as.integer(
    stringr::str_replace(colnames(L_c),
                        ".+\\[([0-9]+),([0-9]+)\\]$", "\\2")
    )

  if (id_with == "binary") {
    replace_original <- "lambda_binary\\[([0-9]+),([0-9]+)\\]$"
  } else if (id_with == "ordinal") {
    replace_original <- "lambda_ordinal\\[([0-9]+),([0-9]+)\\]$"
  } else if (id_with == "metric") {
    replace_original <- "lambda_metric\\[([0-9]+),([0-9]+)\\]$"
  } else {
    stop("Inavlud `id_with` argument")
  }

  colord <- order(var, dim)
  colnames(L_c) <- stringr::str_replace(
    colnames(L_c),
    replace_original,
    "LambdaV\\1_\\2"
  )

  out <- factor.switching::rsp_exact(
    lambda_mcmc = L_c[, colord],
    rotate = FALSE,
    maxIter = 100,
    threshold = 1e-6,
    verbose = FALSE
  )

  return(out)

}


#'
#' @param beta_rsp
#'
#' @return list of
#'
harmonize_varimax <- function (beta_rsp) {
  nChains <- length(beta_rsp)
  cnames <- colnames(beta_rsp[[1]]$lambda_reordered_mcmc)
  d <- dim(beta_rsp[[1]]$lambda_reordered_mcmc)[2]
  q <- as.numeric(strsplit(cnames[d], split = "_")[[1]][2])
  p <- d / q
  lambda_hat_values <- matrix(nrow = nChains, ncol = d)
  colnames(lambda_hat_values) <-
    colnames(beta_rsp[[1]]$lambda_reordered_mcmc)
  for (i in 1:nChains) {
    lambda_hat_values[i, ] <- c(t(beta_rsp[[i]]$lambda_hat))
  }
  tankard <- factor.switching::rsp_exact(
    lambda_mcmc = lambda_hat_values,
    maxIter = 100,
    threshold = 1e-06,
    verbose = TRUE,
    rotate = FALSE)

  reorderedChains <- vector("list", length = nChains)
  for (i in 1:nChains) {
    mcmcIterations <- dim(beta_rsp[[i]]$lambda_reordered_mcmc)[1]
    v_vectors <- matrix(
      tankard$permute_vectors[i, ],
      nrow = mcmcIterations,
      ncol = q,
      byrow = TRUE
    )
    c_vectors <- matrix(
      tankard$sign_vectors[i, ],
      nrow = mcmcIterations,
      ncol = q,
      byrow = TRUE
    )
    lrm <- beta_rsp[[i]]$lambda_reordered_mcmc
    reorderedChains[[i]] <-
      coda::as.mcmc(
        factor.switching::switch_and_permute(
          lambda_mcmc = lrm,
          switch_vectors = c_vectors,
          permute_vectors = v_vectors
        )
      )
  }
  reorderedChains <- coda::as.mcmc.list(reorderedChains)
  return(list(mcmc = reorderedChains,
              sign_vectors = tankard$sign_vectors,
              permute_vectors = tankard$permute_vectors))
}


#'
#' @param raw_drawas
#' @param sign
#'
#' @return
#'
#' @import magrittr
#'
identify_sign <- function (raw_draws, sign) {
  raw_draws_df <- posterior::as_draws_df(raw_draws)

  long_draws <- suppressWarnings(
    raw_draws_df %>%
      dplyr::select(dplyr::matches("^\\.|^lambda")) %>%
      tidyr::pivot_longer(cols = -dplyr::matches("^\\."),
                          names_to = "variable") %>%
      dplyr::group_by(.data$.chain, .data$variable) %>%
      dplyr::summarise(est = mean.default(.data$value), .groups = "drop") %>%
      tidyr::separate(
        col = variable,
        into = c("parameter", "type", "item", "dimension", "extra"),
        convert = TRUE
      )
  )

  sign_df <- data.frame(
    dimension = 1:max(long_draws$dimension),
    sign = rep(sign, length.out = max(long_draws$dimension))
  )

  long_draws <- dplyr::left_join(long_draws, sign_df, by = "dimension")

  chain_flips <- long_draws %>%
    dplyr::group_by(.data$.chain, .data$dimension) %>%
    dplyr::summarise(flip = (sign * mean.default(.data$est)) < 0,
                    .groups = "drop") %>%
    dplyr::select(".chain", "dimension", "flip")

  identified_draws <- raw_draws_df

  for (d in 1:max(chain_flips$dimension)) {
    rows <- identified_draws$.chain %in%
      dplyr::filter(chain_flips, .data$flip & .data$dimension == d)$.chain
    cols <- stringr::str_detect(
      names(identified_draws),
      stringr::str_c("^(eta|lambda).+", d, "\\]")
    )
    suppressWarnings(
      identified_draws[rows, cols] <- -1 * identified_draws[rows, cols]
    )
  }

  return(list(id_draws = identified_draws))
}


#' Label the output
#'
#' @param draws
#' @param regex_pars
#'
#' @import magrittr
#'
#' @export
label_draws <- function (draws, regex_pars = NULL)
{
  # TODO: what to do with regex_pars
  if (is.null(regex_pars)) {
    regex_pars <- c(
      eta = "^eta\\[",
      sigma_eta_evol = "^sigma_eta_evol\\[",
      lambda_binary = "^lambda_binary\\[",
      alpha_binary = "^alpha_binary\\[",
      lambda_ordinal = "^lambda_ordinal\\[",
      kappa = "^kappa\\[",
      alpha_ordinal = "^alpha_ordinal\\[",
      lambda_metric = "^lambda_metric\\[",
      sigma_metric = "^sigma_metric\\[",
      alpha_metric = "^alpha_metric\\[",
      sigma_alpha_evol = "^sigma_alpha_evol$"
    )
  }
  draws_df <- posterior::as_draws_df(draws)
  draws_ls <- vector("list", length(regex_pars))
  names(draws_ls) <- names(regex_pars)

  for (p in seq_along(regex_pars)) {
    regex_pars_p <- regex_pars[p]

    if (!any(stringr::str_detect(names(draws_df), regex_pars_p))) next

    draws_ls[[p]] <-
      draws_df %>%
      as.data.frame() %>%
      dplyr::select(dplyr::matches("^\\."), dplyr::matches(regex_pars_p)) %>%
      tidyr::pivot_longer(
        cols = dplyr::matches(regex_pars_p),
        names_to = "name",
        values_to = "value"
      )

    draws_ls_p <- draws_ls[[p]] %>% as.data.frame()

    if (regex_pars_p == "^eta\\[") {
      draws_ls[[p]] <- draws_ls_p %>%
        dplyr::mutate(
          par = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 1],
          time = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 2],
          unit = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 3],
          dim = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 4],
          dplyr::across(.data$time:.data$dim, as.integer),
          TIME = factor(.data$time, labels = attr(draws, "time_labels")),
          UNIT = factor(.data$unit, labels = attr(draws, "unit_labels"))
        ) %>%
        dplyr::select("par", "TIME", "UNIT", "dim", "value",
                      dplyr::everything())
    } else if (regex_pars_p == "sigma_eta_evol\\[") {
      draws_ls[[p]] <- draws_ls_p %>%
        dplyr::mutate(
          par = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 1],
          dim = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 2],
          dplyr::across(.data$dim, as.integer),
        ) %>%
        dplyr::select("par", "dim", "value", dplyr::everything())
    } else if (regex_pars_p == "^lambda_binary\\[") {
      draws_ls[[p]] <- draws_ls_p %>%
        dplyr::mutate(
          par = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 1],
          item = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 2],
          dim = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 3],
          dplyr::across(.data$item:.data$dim, as.integer),
          ITEM = factor(
            x = .data$item,
            labels = attr(draws, "binary_item_labels")
          )
        ) %>%
        dplyr::select("par", "ITEM", "dim", "value", dplyr::everything())
    } else if (regex_pars_p == "^lambda_metric\\[") {
      draws_ls[[p]] <- draws_ls_p %>%
        dplyr::mutate(
          par = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 1],
          item = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 2],
          dim = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 3],
          dplyr::across(.data$item:.data$dim, as.integer),
          ITEM = factor(
            x = .data$item,
            labels = attr(draws, "metric_item_labels")
          )
        ) %>%
        dplyr::select("par", "ITEM", "dim", "value", dplyr::everything())
    } else if (regex_pars_p == "^sigma_metric\\[") {
      draws_ls[[p]] <- draws_ls_p %>%
        dplyr::mutate(
          par = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 1],
          item = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 2],
          dplyr::across(.data$item, as.integer),
          ITEM = factor(
            x = .data$item,
            labels = attr(draws, "metric_item_labels")
          )
        ) %>%
        dplyr::select("par", "ITEM", "value", dplyr::everything())
    } else if (regex_pars_p == "^alpha_metric\\[") {
      draws_ls[[p]] <- draws_ls_p %>%
        dplyr::mutate(
          par = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 1],
          time = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 2],
          item = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 3],
          dplyr::across(.data$time:.data$item, as.integer),
          TIME = factor(.data$time, labels = attr(draws, "time_labels")),
          ITEM = factor(
            x = .data$item,
            labels = attr(draws, "metric_item_labels")
          )
        ) %>%
        dplyr::select("par", "TIME", "ITEM", "value", dplyr::everything())
    } else if (regex_pars_p == "^alpha_binary\\[") {
      draws_ls[[p]] <- draws_ls_p %>%
        dplyr::mutate(
          par = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 1],
          time = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 2],
          item = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 3],
          dplyr::across(.data$time:.data$item, as.integer),
          TIME = factor(.data$time, labels = attr(draws, "time_labels")),
          ITEM = factor(
            x = .data$item,
            labels = attr(draws, "binary_item_labels")
          )
        ) %>%
        dplyr::select("par", "TIME", "ITEM", "value", dplyr::everything())
    } else if (regex_pars_p == "^lambda_ordinal\\[") {
      draws_ls[[p]] <- draws_ls_p %>%
        dplyr::mutate(
          par = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 1],
          item = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 2],
          dim = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 3],
          dplyr::across(.data$item:.data$dim, as.integer),
          ITEM = factor(
            x = .data$item,
            labels = attr(draws, "ordinal_item_labels")
          )
        ) %>%
        dplyr::select("par", "ITEM", "dim", "value", dplyr::everything())
    } else if (regex_pars_p == "^kappa\\[") {
      draws_ls[[p]] <- draws_ls_p %>%
        dplyr::mutate(
          par = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 1],
          item = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 2],
          threshold =
            stringr::str_split(.data$name,
                                    "[\\[,\\]]", simplify = TRUE)[, 3],
          dplyr::across(.data$threshold, as.integer),
          ITEM = factor(
            x = .data$item,
            labels = attr(draws, "ordinal_item_labels")
          )
        ) %>%
        dplyr::select("par", "ITEM", "threshold", "value", dplyr::everything())
    } else if (regex_pars_p == "^alpha_ordinal\\[") {
      draws_ls[[p]] <- draws_ls_p %>%
        dplyr::mutate(
          par = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 1],
          time = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 2],
          item = stringr::str_split(.data$name,
                                  "[\\[,\\]]", simplify = TRUE)[, 3],
          dplyr::across(.data$time:.data$item, as.integer),
          TIME = factor(.data$time, labels = attr(draws, "time_labels")),
          ITEM = factor(
            x = .data$item,
            labels = attr(draws, "ordinal_item_labels")
          )
        ) %>%
        dplyr::select("par", "TIME", "ITEM", "value", dplyr::everything())
    } else if (regex_pars_p == "^sigma_alpha_evol") {
      draws_ls[[p]] <- draws_ls_p %>%
        dplyr::rename(par = "name") %>%
        dplyr::select("par", "value", dplyr::everything())
    } else {
      # pass the element that does not match regex
      next
    }
  }

  attr(draws_ls, "unit_labels") <- attr(draws, "unit_labels")
  attr(draws_ls, "time_labels") <- attr(draws, "time_labels")
  attr(draws_ls, "binary_item_labels") <- attr(draws, "binary_item_labels")
  attr(draws_ls, "ordinal_item_labels") <- attr(draws, "ordinal_item_labels")
  attr(draws_ls, "metric_item_labels") <- attr(draws, "metric_item_labels")

  return(draws_ls)
}
