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
#'     ordinal. If `NA` (the default), then any item with between 3 and
#'     `max_cats` unique values will be treated as ordinal.
#' @param max_cats (positive integer) Maximum number of response categories for
#'     ordinal items. Defaults to `10`.
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
#' @param quiet (logical) Suppress messages reporting how items were
#'     categorized? Defaults to `FALSE`.
#'
#' @return A list formatted for Stan
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
                         max_cats = 10,
                         standardize = TRUE,
                         make_indicator_for_zeros = TRUE,
                         periods_to_estimate,
                         quiet = FALSE) {

    if (missing(periods_to_estimate)) {
        periods_to_estimate <-
            min(long_data[[time_var]]):max(long_data[[time_var]])
    }

    use_data <- long_data[long_data[[time_var]] %in% periods_to_estimate, ]
    use_data$unit <- use_data[[unit_var]]
    use_data$UNIT <- factor(use_data[[unit_var]])
    use_data$time <- use_data[[time_var]]
    use_data$TIME <- factor(use_data[[time_var]], periods_to_estimate)
    use_data$item <- use_data[[item_var]]
    use_data$value <- as.numeric(use_data[[value_var]])
    use_data <- dplyr::select(use_data,
                              "unit", "UNIT", "time",
                              "TIME", "item", "value")
    use_data <- tidyr::drop_na(use_data)
    stopifnot(!anyNA(use_data$unit))
    stopifnot(!anyNA(use_data$time))
    stopifnot(!anyNA(use_data$item))
    stopifnot(!anyNA(use_data$value))

    items <- sort(unique(use_data$item))

    unique_df <- use_data |>
        dplyr::summarise(.by = "item", n = length(unique(.data$value)))

    drop_items <- dplyr::filter(unique_df, .data$n < 2)$item
    if (length(drop_items) > 0) {
        cat("\nDropping the following items due to lack of variation:\n")
        cat(c("  *", paste(drop_items, collapse = "\n  * "), "\n"))
    }
    if (length(binary_items) == 1L && is.na(binary_items)) {
        binary_items <- unique_df$item[unique_df$n == 2]
        binary_items <- sort(setdiff(binary_items, drop_items))
    }
    if (length(trichotomous_items) == 1L && is.na(trichotomous_items)) {
        trichotomous_items <- unique_df$item[unique_df$n == 3]
        trichotomous_items <- sort(setdiff(trichotomous_items, drop_items))
    }
    if (length(ordinal_items) == 1L && is.na(ordinal_items)) {
        ordinal_items <- unique_df$item[unique_df$n <= max_cats]
        ordinal_items <- setdiff(
            ordinal_items,
            c(drop_items, binary_items, trichotomous_items)
        )
        ordinal_items <- sort(ordinal_items)
    }
    metric_items <- setdiff(
        items,
        c(binary_items, trichotomous_items, ordinal_items, drop_items)
    )
    metric_items <- sort(metric_items)
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
