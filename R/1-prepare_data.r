#' Prepare data for MODGIRT
#'
#' This function prepares the data for the MODGIRT analysis by creating a
#' four-dimensional cross-tabulation and setting up the necessary variables 
#' and parameters.
#'
#' @param long_data The long-format data frame containing the raw data.
#' @param unit_var The name of the variable representing the unit identifier.
#' @param time_var The name of the variable representing the time identifier.
#' @param item_var The name of the variable representing the item identifier.
#' @param value_var The name of the variable representing the repsonse value.
#' @param weight_var The name of the variable representing subject-specific
#'   weights (optional).
#' @param periods_to_estimate The periods to estimate (optional).
#'
#' @return A list containing the dimensions of the data (T, G, Q, K, D), the
#'  four-dimensional cross-tabulation (SSSS), the matrix of nonzero loadings
#'  (beta_nonzero), and the matrix of signed loadings (beta_sign).
#'
#' @export
shape_modgirt <- function(long_data,
                               unit_var,
                               time_var,
                               item_var,
                               value_var,
                               weight_var = NULL,
                               periods_to_estimate) {
    if (missing(periods_to_estimate)) {
        if (is.factor(long_data[[time_var]])) {
            periods_to_estimate <- levels(long_data[[time_var]])
        } else {
            periods_to_estimate <-
                min(long_data[[time_var]]):max(long_data[[time_var]])
        }
    }
    ## Subset data and create identifier variables
    use_data <- long_data[long_data[[time_var]] %in% periods_to_estimate, ]
    use_data$unit <- use_data[[unit_var]]
    use_data$time <- use_data[[time_var]]
    use_data$item <- use_data[[item_var]]
    use_data$value <- as.numeric(use_data[[value_var]])
    stopifnot(!anyNA(use_data$unit))
    stopifnot(!anyNA(use_data$time))
    stopifnot(!anyNA(use_data$item))
    stopifnot(!anyNA(use_data$value))
    ## Create factors (drops absent levels)
    use_data$UNIT <- factor(use_data[[unit_var]])
    use_data$TIME <- factor(use_data[[time_var]], levels = periods_to_estimate)
    use_data$ITEM <- factor(use_data[[item_var]])
    ## Create four-dimensional cross-tabulation
    count_array <- create_counts(long_data = use_data, weight_var = weight_var)
    n_time <- dim(count_array)[1]
    n_unit <- dim(count_array)[2]
    n_item <- dim(count_array)[3]
    n_value <- dim(count_array)[4]
    stan_data <- list(
        T = n_time,
        G = n_unit,
        Q = n_item,
        K = n_value,
        SSSS = count_array
    )
    attr(stan_data, "unit_labels") <- levels(use_data$UNIT)
    attr(stan_data, "time_labels") <- levels(use_data$TIME)
    attr(stan_data, "item_labels") <- levels(use_data$ITEM)
    return(stan_data)
}

#' Check long-format data for problems affecting a mixed factor model
#'
#' Reports data problems that [shape_mixfac()] would otherwise pass over
#' silently, or handle in ways that may surprise. It is called by
#' [shape_mixfac()], and can also be run beforehand to inspect data.
#'
#' @param long_data (data frame) Data in long (unit-period-item) form.
#' @param unit_var,time_var,item_var,value_var (strings) Column names, as for
#'     [shape_mixfac()].
#' @param periods_to_estimate (vector or `NULL`) Periods to retain before
#'     checking. `NULL`, the default, uses every period present.
#' @param binary_items,trichotomous_items,ordinal_items (character vectors)
#'     Items of each type. `NA`, the default, classifies items by their number
#'     of distinct values, as [shape_mixfac()] does.
#' @param max_cats (positive integer) Largest number of distinct values for
#'     which an item is classified as ordinal. Defaults to `5`.
#' @param min_cat_n (non-negative integer) Warn when a trichotomous or ordinal
#'     item has a response category with fewer than this many observations.
#'     Cutpoints bounding sparse categories are weakly identified. `0`
#'     disables the check. Defaults to `3`.
#' @param quiet (logical) Suppress the informational summary? Warnings are
#'     emitted regardless.
#'
#' @return Invisibly, a list of data frames: `duplicates`, `no_variation`,
#'     `single_period`, `sparse_categories`, `missing_levels`, and
#'     `recoded_binary`, each with zero rows if that check found nothing.
#'
#' @details
#' The checks are:
#'
#' * **Duplicated observations.** [shape_mixfac()] takes the data at face
#'   value, so a repeated unit-period-item combination enters the likelihood
#'   more than once.
#' * **Items with no variation**, which [shape_mixfac()] drops.
#' * **Items observed in a single period**, which inform the loadings but
#'   contribute nothing to the estimated dynamics.
#' * **Sparse response categories** in trichotomous and ordinal items.
#' * **Missing intermediate levels.** [shape_mixfac()] recodes ordered
#'   responses with `as.integer(ordered(value))`, which numbers the observed
#'   values consecutively. An item coded 1, 2, 4, 5 therefore becomes a
#'   four-category item, and its second cutpoint separates the original values
#'   2 and 4 rather than 2 and 3.
#' * **Binary items not coded 0/1**, which are recoded to 0 and 1 in the
#'   order of their values.
#'
#' @examples
#' data("state_policies_2010_2012")
#' invisible(check_mixfac_data(
#'     subset(state_policies_2010_2012, !is.na(value_real)),
#'     unit_var = "state_abb", time_var = "year",
#'     item_var = "policy_variable", value_var = "value_real",
#'     min_cat_n = 0, quiet = TRUE
#' ))
#'
#' @export
check_mixfac_data <- function(long_data,
                              unit_var,
                              time_var,
                              item_var,
                              value_var,
                              periods_to_estimate = NULL,
                              binary_items = NA,
                              trichotomous_items = NA,
                              ordinal_items = NA,
                              max_cats = 5,
                              min_cat_n = 3,
                              quiet = FALSE) {

    quiet <- check_flag(quiet)
    absent <- setdiff(c(unit_var, time_var, item_var, value_var),
                      names(long_data))
    if (length(absent)) {
        cli::cli_abort("Column{?s} {.val {absent}} not found in {.arg long_data}.")
    }
    if (!is.numeric(min_cat_n) || length(min_cat_n) != 1L || min_cat_n < 0) {
        cli::cli_abort("{.arg min_cat_n} must be a single non-negative number.")
    }

    d <- data.frame(
        unit  = as.character(long_data[[unit_var]]),
        time  = long_data[[time_var]],
        item  = as.character(long_data[[item_var]]),
        value = suppressWarnings(as.numeric(long_data[[value_var]])),
        stringsAsFactors = FALSE
    )
    if (!is.null(periods_to_estimate)) {
        d <- d[d$time %in% periods_to_estimate, , drop = FALSE]
    }
    if (!nrow(d)) {
        cli::cli_abort("No observations remain after subsetting to
                        {.arg periods_to_estimate}.")
    }

    empty <- function(...) {
        stats::setNames(rep(list(character(0)), length(c(...))), c(...))
    }
    out <- list(
        duplicates        = as.data.frame(empty("unit", "time", "item", "n")),
        no_variation      = as.data.frame(empty("item")),
        single_period     = as.data.frame(empty("item", "period")),
        sparse_categories = as.data.frame(empty("item", "value", "n")),
        missing_levels    = as.data.frame(empty("item", "missing")),
        recoded_binary    = as.data.frame(empty("item", "values"))
    )

    ## ---- Duplicated unit-period-item combinations -----------------------
    key <- paste(d$unit, d$time, d$item, sep = "\r")
    if (anyDuplicated(key)) {
        tab <- table(key)
        tab <- tab[tab > 1L]
        parts <- do.call(rbind, strsplit(names(tab), "\r", fixed = TRUE))
        out$duplicates <- data.frame(
            unit = parts[, 1L], time = parts[, 2L], item = parts[, 3L],
            n = as.integer(tab), stringsAsFactors = FALSE
        )
        n_items_dup <- length(unique(out$duplicates$item))
        cli::cli_warn(c(
            "{nrow(out$duplicates)} unit-period-item combination{?s}
             appear{?s/} more than once.",
            "!" = "{.fn shape_mixfac} does not de-duplicate, so these
                   responses would enter the likelihood more than once.",
            "i" = "{cli::qty(n_items_dup)}Affected item{?s}:
                   {.val {utils::head(unique(out$duplicates$item), 5)}}"
        ))
    }

    ## ---- Classify items, dropping missing responses first ---------------
    obs <- d[!is.na(d$value), , drop = FALSE]
    n_missing <- nrow(d) - nrow(obs)
    n_unique <- vapply(split(obs$value, obs$item),
                       function(x) length(unique(x)), integer(1))
    auto <- classify_mixfac_items(n_unique, max_cats = max_cats)

    use_auto <- function(x) length(x) == 1L && is.na(x)
    b_items <- if (use_auto(binary_items)) {
                   names(auto)[auto == "binary"]
               } else binary_items
    t_items <- if (use_auto(trichotomous_items)) {
                   names(auto)[auto == "trichotomous"]
               } else trichotomous_items
    o_items <- if (use_auto(ordinal_items)) {
                   names(auto)[auto == "ordinal"]
               } else ordinal_items
    ordered_items <- c(t_items, o_items)

    ## ---- Items with no variation ----------------------------------------
    flat <- names(n_unique)[n_unique < 2L]
    if (length(flat)) {
        out$no_variation <- data.frame(item = flat, stringsAsFactors = FALSE)
        cli::cli_warn(c(
            "{length(flat)} item{?s} {?has/have} fewer than two distinct
             values and will be dropped.",
            "i" = "{.val {flat}}"
        ))
    }

    ## ---- Items observed in a single period ------------------------------
    n_period_total <- length(unique(obs$time))
    if (n_period_total > 1L) {
        n_per_item <- vapply(split(obs$time, obs$item),
                             function(x) length(unique(x)), integer(1))
        one <- setdiff(names(n_per_item)[n_per_item == 1L], flat)
        if (length(one)) {
            out$single_period <- data.frame(
                item = one,
                period = vapply(one, function(i) {
                    as.character(unique(obs$time[obs$item == i]))
                }, character(1)),
                stringsAsFactors = FALSE
            )
            cli::cli_warn(c(
                "{length(one)} item{?s} {?is/are} observed in only one of
                 {n_period_total} periods.",
                "i" = "{cli::qty(length(one))}{?It/They} inform{?s/} the
                       loadings but contribute{?s/} nothing to the estimated
                       dynamics.",
                "i" = "{.val {one}}"
            ))
        }
    }

    ## ---- Sparse categories, and missing intermediate levels -------------
    sparse_rows <- gap_rows <- list()
    for (i in intersect(ordered_items, names(n_unique))) {
        v <- obs$value[obs$item == i]
        tab <- table(v)

        if (min_cat_n > 0 && any(tab < min_cat_n)) {
            sparse_rows[[i]] <- data.frame(
                item = i,
                value = names(tab)[tab < min_cat_n],
                n = as.integer(tab[tab < min_cat_n]),
                stringsAsFactors = FALSE
            )
        }

        ## Gaps are only meaningful for integer-coded scales
        u <- sort(unique(v))
        if (length(u) > 1L && all(u == round(u))) {
            gaps <- setdiff(seq(min(u), max(u)), u)
            if (length(gaps)) {
                gap_rows[[i]] <- data.frame(
                    item = i, missing = paste(gaps, collapse = ", "),
                    stringsAsFactors = FALSE
                )
            }
        }
    }

    if (length(sparse_rows)) {
        out$sparse_categories <- do.call(rbind, sparse_rows)
        row.names(out$sparse_categories) <- NULL
        cli::cli_warn(c(
            "{length(sparse_rows)} ordered item{?s} {?has/have} response
             categories with fewer than {min_cat_n} observations.",
            "!" = "Cutpoints bounding sparse categories are weakly identified,
                   and can take implausible values and impede sampling.",
            "i" = "Consider collapsing categories, or naming the genuinely
                   ordinal items in {.arg ordinal_items}.",
            "i" = "{.val {names(sparse_rows)}}"
        ))
    }

    if (length(gap_rows)) {
        out$missing_levels <- do.call(rbind, gap_rows)
        row.names(out$missing_levels) <- NULL
        cli::cli_warn(c(
            "{length(gap_rows)} ordered item{?s} {?has/have} unobserved
             intermediate levels.",
            "!" = "{.fn shape_mixfac} numbers the observed values
                   consecutively, so unobserved levels are silently closed up
                   and the cutpoints do not correspond to the original codes.",
            "i" = paste0("{.val {out$missing_levels$item}}: missing ",
                         "{.val {out$missing_levels$missing}}")
        ))
    }

    ## ---- Binary items not coded 0/1 -------------------------------------
    odd_binary <- list()
    for (i in intersect(b_items, names(n_unique))) {
        u <- sort(unique(obs$value[obs$item == i]))
        if (!identical(as.numeric(u), c(0, 1))) {
            odd_binary[[i]] <- data.frame(
                item = i, values = paste(u, collapse = ", "),
                stringsAsFactors = FALSE
            )
        }
    }
    if (length(odd_binary)) {
        out$recoded_binary <- do.call(rbind, odd_binary)
        row.names(out$recoded_binary) <- NULL
        if (!quiet) {
            cli::cli_alert_info(
                "{length(odd_binary)} binary item{?s} {?is/are} not coded 0/1;
                 the lower value is recoded to 0 and the higher to 1."
            )
        }
    }

    ## ---- Informational summary ------------------------------------------
    if (!quiet) {
        if (n_missing > 0) {
            cli::cli_alert_info(
                "Dropping {n_missing} observation{?s} with a missing response."
            )
        }
        if (sum(vapply(out, nrow, integer(1))) == 0L) {
            cli::cli_alert_success("No data problems detected.")
        }
    }

    invisible(out)
}

#' Prepare data for dynamic mixed factor model
#'
#' @param long_data (data frame) Data in long (unit-period-item) form
#' @param unit_var (string) Name of variable identifying units
#' @param time_var (string) Name of variable identifying time periods
#' @param item_var (string) Name of variable identifying items
#' @param value_var (string) Name of the response variable
#' @param binary_items (character vector) Items that should be treated as
#'     binary. If `NA` (the default), then any item with 2 unique values will be
#'     treated as binary.
#' @param trichotomous_items (character vector) Items that should be treated as
#'     trichotomous. If `NA` (the default), then any item with 3 unique values
#'     will be treated as trichotomous.
#' @param ordinal_items (character vector) Items that should be treated as
#'     ordinal. If `NA` (the default), then any item with more than 3 and no
#'     more than `max_cats` unique values will be treated as ordinal.
#' @param max_cats (positive integer) Maximum number of response categories for
#'     ordinal items. Items with more distinct values than this are treated as
#'     metric. Defaults to `5`.
#' @param standardize (logical) Should metric items be standardized? Defaults to
#'     `TRUE`.
#' @param make_indicator_for_zeros (logical) Should metric items with a lower
#'     bound of zero be modeled as zero-inflated? If `TRUE` (the default), then
#'     for each such variable observations with a response value of zero will be
#'     set to missing and dropped, and a separate dummy variable ending in "_zi"
#'     will be created to indicate observations with responses greater than
#'     zero.
#' @param periods_to_estimate (vector) Values of `time_var` for which to
#'     estimate parameter values. If `NULL` (the default), `periods_to_estimate`
#'     will be set to `min(long_data[[time_var]]):max(long_data[[time_var]])`,
#'     under the assumption that `long_data[[time_var]]` is an integer
#'     vector. Data from periods not in `periods_to_estimate` will be
#'     dropped. If `periods_to_estimate` includes a period for which there is no
#'     data, parameter values for that year will be imputed by the dynamic
#'     model.
#' @param min_cat_n (non-negative integer) Warn when a trichotomous or ordinal
#'     item has a response category with fewer than this many observations,
#'     since the cutpoints bounding sparse categories are weakly identified.
#'     Passed to [check_mixfac_data()]. `0` disables the check. Defaults to
#'     `3`.
#' @param check (logical) Run [check_mixfac_data()] on the input? Defaults to
#'     `TRUE`.
#' @param quiet (logical) Suppress messages reporting how items were
#'     categorized? Defaults to `FALSE`. Warnings about data problems are
#'     emitted regardless.
#'
#' @return A list formatted for Stan
#'
#' @seealso [check_mixfac_data()]
#'
#' @importFrom rlang .data
#'
#' @export
shape_mixfac <- function(long_data,
                         unit_var,
                         time_var,
                         item_var,
                         value_var,
                         binary_items = NA,
                         trichotomous_items = NA,
                         ordinal_items = NA,
                         max_cats = 5,
                         standardize = TRUE,
                         make_indicator_for_zeros = TRUE,
                         periods_to_estimate,
                         min_cat_n = 3,
                         check = TRUE,
                         quiet = FALSE) {

    quiet <- check_flag(quiet)
    check <- check_flag(check)
    standardize <- check_flag(standardize)
    make_indicator_for_zeros <- check_flag(make_indicator_for_zeros)

    absent <- setdiff(c(unit_var, time_var, item_var, value_var),
                      names(long_data))
    if (length(absent)) {
        cli::cli_abort("Column{?s} {.val {absent}} not found in {.arg long_data}.")
    }

    if (missing(periods_to_estimate)) {
        periods_to_estimate <-
            min(long_data[[time_var]]):max(long_data[[time_var]])
    }

    use_data <- long_data[long_data[[time_var]] %in% periods_to_estimate, ]
    if (!nrow(use_data)) {
        cli::cli_abort(c(
            "No observations fall in {.arg periods_to_estimate}.",
            "i" = "Periods present in the data:
                   {.val {sort(unique(long_data[[time_var]]))}}"
        ))
    }
    use_data$unit <- use_data[[unit_var]]
    use_data$UNIT <- factor(use_data[[unit_var]])
    use_data$time <- use_data[[time_var]]
    use_data$TIME <- factor(use_data[[time_var]], periods_to_estimate)
    use_data$item <- use_data[[item_var]]
    use_data$value <- suppressWarnings(as.numeric(use_data[[value_var]]))
    use_data <- dplyr::select(use_data,
                              "unit", "UNIT", "time",
                              "TIME", "item", "value")

    ## Retained before dropping missing responses, so that the checks can
    ## report how many were dropped and can see duplicated observations
    raw_data <- use_data

    use_data <- tidyr::drop_na(use_data)
    if (!nrow(use_data)) {
        cli::cli_abort("Every observation has a missing response.")
    }

    ## Units were factored before missing responses were dropped, so a unit
    ## with no observed response would otherwise survive as a level with no
    ## data, inflating J and leaving its factor scores informed only by their
    ## prior.
    absent_units <- setdiff(levels(use_data$UNIT),
                            unique(as.character(use_data$UNIT)))
    if (length(absent_units)) {
        cli::cli_warn(c(
            "Dropping {length(absent_units)} unit{?s} with no observed
             responses.",
            "i" = "{.val {absent_units}}"
        ))
        use_data$UNIT <- droplevels(use_data$UNIT)
    }

    ## Periods with no data are retained deliberately: the dynamic model
    ## imputes them. Worth saying so, since it is easy to do by accident.
    absent_periods <- setdiff(levels(use_data$TIME),
                              unique(as.character(use_data$TIME)))
    if (!quiet && length(absent_periods)) {
        cli::cli_alert_info(
            "{length(absent_periods)} period{?s} in {.arg periods_to_estimate}
             {?has/have} no data ({.val {absent_periods}}).
             {cli::qty(length(absent_periods))}Parameters for {?it/them} will
             be imputed by the dynamic model."
        )
    }

    items <- sort(unique(use_data$item))
    n_unique <- vapply(split(use_data$value, use_data$item),
                       function(x) length(unique(x)), integer(1))
    auto <- classify_mixfac_items(n_unique, max_cats = max_cats)
    drop_items <- sort(names(auto)[auto == "dropped"])

    if (length(binary_items) == 1L && is.na(binary_items)) {
        binary_items <- sort(names(auto)[auto == "binary"])
    }
    if (length(trichotomous_items) == 1L && is.na(trichotomous_items)) {
        trichotomous_items <- sort(names(auto)[auto == "trichotomous"])
    }
    if (length(ordinal_items) == 1L && is.na(ordinal_items)) {
        ordinal_items <- sort(names(auto)[auto == "ordinal"])
    }
    metric_items <- sort(setdiff(
        items,
        c(binary_items, trichotomous_items, ordinal_items, drop_items)
    ))

    ## Run the checks against the data as supplied, with the item
    ## classification resolved above, so that the two cannot disagree. Items
    ## created below by `make_indicator_for_zeros` are 0/1 by construction and
    ## need no checking.
    if (check) {
        check_mixfac_data(
            long_data          = raw_data,
            unit_var           = "unit",
            time_var           = "time",
            item_var           = "item",
            value_var          = "value",
            binary_items       = binary_items,
            trichotomous_items = trichotomous_items,
            ordinal_items      = ordinal_items,
            max_cats           = max_cats,
            min_cat_n          = min_cat_n,
            quiet              = quiet
        )
    }
    if (make_indicator_for_zeros) {
        for (i in seq_along(metric_items)) {
            obs_i <- which(use_data$item == metric_items[i])
            is_zero <- use_data$value[obs_i] == 0
            if (any(is_zero) && 0 %in% range(use_data$value[obs_i])) {
                ## TODO: Add collinearity check
                use_data$value[obs_i[is_zero]] <- NA_real_
                zi_i <- stringr::str_c(metric_items[i], "_zi")
                newdat <- use_data[obs_i, ]
                newdat$value <- as.integer(!is_zero)
                newdat$item <- zi_i
                use_data <- dplyr::bind_rows(use_data, newdat)
                binary_items <- c(binary_items, zi_i)
            }
        }
    }

    if (!quiet && length(binary_items) > 0) {
        cli::cli_alert_info("Categorizing {length(binary_items)} item{?s} as binary.")
        cli::cli_ul(binary_items)
    }
    if (!quiet && length(trichotomous_items) > 0) {
        cli::cli_alert_info("Categorizing {length(trichotomous_items)} item{?s} as trichotomous.")
        cli::cli_ul(trichotomous_items)
    }
    if (!quiet && length(ordinal_items) > 0) {
        cli::cli_alert_info("Categorizing {length(ordinal_items)} item{?s} as ordinal.")
        cli::cli_ul(ordinal_items)
    }
    if (!quiet && length(metric_items) > 0) {
        cli::cli_alert_info("Categorizing {length(metric_items)} item{?s} as metric.")
        cli::cli_ul(metric_items)
    }

    binary_data <- use_data |>
        dplyr::filter(.data$item %in% binary_items) |>
        dplyr::mutate(ITEM = factor(.data$item, levels = binary_items)) |>
        dplyr::group_by(.data$ITEM) |>
        dplyr::mutate(yy = as.integer(ordered(.data$value)) - 1L) |>
        dplyr::filter(!is.na(.data$yy)) |>
        dplyr::ungroup() |>
        dplyr::arrange(.data$TIME, .data$ITEM, .data$UNIT) # time must vary last

    trichotomous_data <- use_data |>
        dplyr::filter(.data$item %in% trichotomous_items) |>
        dplyr::mutate(ITEM = factor(.data$item, levels = trichotomous_items)) |>
        dplyr::group_by(.data$ITEM) |>
        dplyr::mutate(yy = as.integer(ordered(.data$value))) |>
        dplyr::filter(!is.na(.data$yy)) |>
        dplyr::ungroup() |>
        dplyr::arrange(.data$TIME, .data$ITEM, .data$UNIT) # time must vary last

    ordinal_data <- use_data |>
        dplyr::filter(.data$item %in% ordinal_items) |>
        dplyr::mutate(ITEM = factor(.data$item, levels = ordinal_items)) |>
        dplyr::group_by(.data$ITEM) |>
        dplyr::mutate(yy = as.integer(ordered(.data$value))) |>
        dplyr::filter(!is.na(.data$yy)) |>
        dplyr::ungroup() |>
        dplyr::arrange(.data$TIME, .data$ITEM, .data$UNIT) # time must vary last

    metric_data <- use_data |>
        dplyr::filter(.data$item %in% metric_items) |>
        dplyr::mutate(ITEM = factor(.data$item, levels = metric_items)) |>
        dplyr::mutate(yy = .data$value) |>
        dplyr::filter(!is.na(.data$yy)) |>
        dplyr::arrange(.data$TIME, .data$ITEM, .data$UNIT) # time must vary last

    if (standardize) {
        metric_data <- metric_data |>
            dplyr::group_by(.data$ITEM) |>
            dplyr::mutate(yy = stdize(.data$yy)) |>
            dplyr::ungroup()
    }

    stan_data <- list(
        J = nlevels(use_data$UNIT),
        T = nlevels(use_data$TIME),
        N_binary = nrow(binary_data),
        I_binary = nlevels(binary_data$ITEM),
        yy_binary = as.integer(binary_data$yy),
        ii_binary = as.integer(binary_data$ITEM),
        jj_binary = as.integer(binary_data$UNIT),
        tt_binary = as.integer(binary_data$TIME),
        N_trichot = nrow(trichotomous_data),
        I_trichot = nlevels(trichotomous_data$ITEM),
        yy_trichot = as.integer(trichotomous_data$yy),
        ii_trichot = as.integer(trichotomous_data$ITEM),
        jj_trichot = as.integer(trichotomous_data$UNIT),
        tt_trichot = as.integer(trichotomous_data$TIME),
        N_ordinal = nrow(ordinal_data),
        I_ordinal = nlevels(ordinal_data$ITEM),
        K_ordinal = if (nrow(ordinal_data) > 0) max(ordinal_data$yy) else 1L,
        yy_ordinal = as.integer(ordinal_data$yy),
        ii_ordinal = as.integer(ordinal_data$ITEM),
        jj_ordinal = as.integer(ordinal_data$UNIT),
        tt_ordinal = as.integer(ordinal_data$TIME),
        N_metric = nrow(metric_data),
        I_metric = nlevels(metric_data$ITEM),
        yy_metric = metric_data$yy,
        ii_metric = as.integer(metric_data$ITEM),
        jj_metric = as.integer(metric_data$UNIT),
        tt_metric = as.integer(metric_data$TIME)
    )
    attr(stan_data, "unit_labels") <- levels(use_data$UNIT)
    attr(stan_data, "time_labels") <- levels(use_data$TIME)
    attr(stan_data, "binary_item_labels") <- levels(binary_data$ITEM)
    attr(stan_data, "trichotomous_item_labels") <- levels(trichotomous_data$ITEM)
    attr(stan_data, "ordinal_item_labels") <- levels(ordinal_data$ITEM)
    attr(stan_data, "metric_item_labels") <- levels(metric_data$ITEM)

    class(stan_data) <- c("mixfac_data", class(stan_data))

    return(stan_data)
}


#' Build a complete item-by-period grid for each item type
#'
#' Expands the item and period indices of a `mixfac_data` object into a full
#' grid and joins on the observed responses, so that unobserved item-periods
#' appear explicitly with `value = NA`.
#'
#' @param shaped_data (mixfac_data) Output of [shape_mixfac()].
#' @param types (character vector) Item types to include.
#'
#' @return A data frame with columns `item_type`, `item`, `time`, `ITEM`,
#'     `TIME`, and `value`.
#'
#' @keywords internal
#' @noRd
make_item_time_grid <- function(shaped_data,
                                types = c("binary", "trichotomous",
                                          "ordinal", "metric")) {
    types <- match.arg(types, several.ok = TRUE)
    n_time <- shaped_data$T
    time_labels <- attr(shaped_data, "time_labels")
    grid_ls <- vector("list", length(types))
    names(grid_ls) <- types
    for (n in seq_along(types)) {
        typ <- stringr::str_remove(types[n], "omous")
        n_item <- shaped_data[[paste0("I_", typ)]]
        if (is.null(n_item) || n_item == 0L) next
        item_labels <- attr(shaped_data, paste0(types[n], "_item_labels"))
        observed <- data.frame(
            item  = as.integer(shaped_data[[paste0("ii_", typ)]]),
            time  = as.integer(shaped_data[[paste0("tt_", typ)]]),
            value = as.numeric(shaped_data[[paste0("yy_", typ)]])
        )
        grid <- expand.grid(
            item = seq_len(n_item),
            time = seq_len(n_time)
        )
        grid <- dplyr::mutate(
            grid,
            ITEM = factor(.data$item, levels = seq_len(n_item),
                          labels = item_labels),
            TIME = factor(.data$time, levels = seq_len(n_time),
                          labels = time_labels)
        )
        grid_ls[[n]] <- dplyr::left_join(grid, observed,
                                         by = c("item", "time"))
    }
    dplyr::bind_rows(grid_ls, .id = "item_type")
}
