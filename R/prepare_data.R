#' prepare data for dynIRT
#'
#' @param long_data
#' @param unit_var
#' @param time_var
#' @param item_var
#' @param value_var
#' @param ordinal_items
#' @param binary_items
#' @param max_cats
#' @param standardize
#' @param make_indicator_for_zeros
#' @param periods_to_estimate
#'
#' @return A list that is formatted for stan
#'
#' @import magrittr
#'
#' @export
shape_data <- function (long_data,
                        unit_var,
                        time_var,
                        item_var,
                        value_var,
                        ordinal_items = NA,
                        binary_items = NA,
                        max_cats = 10,
                        standardize = TRUE,
                        make_indicator_for_zeros = TRUE,
                        periods_to_estimate = NULL)
{
  # TODO: check the input


  stdize <- function (x) (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)

  if (is.null(periods_to_estimate)) {
    periods_to_estimate <-
      min(long_data[[time_var]]):max(long_data[[time_var]])
  }

  long_data <- long_data[long_data[[time_var]] %in% periods_to_estimate, ]
  long_data$unit <- long_data[[unit_var]]
  long_data$UNIT <- factor(long_data[[unit_var]])
  long_data$time <- long_data[[time_var]]
  long_data$TIME <- factor(long_data[[time_var]], periods_to_estimate)
  long_data$item <- long_data[[item_var]]
  long_data$value <- as.numeric(long_data[[value_var]])
  long_data <- dplyr::select(long_data, unit, UNIT,
                            time, TIME, item, value)
  items <- sort(unique(long_data$item))

  unique_df <- long_data %>%
    dplyr::group_by(.data$item) %>%
    dplyr::summarise(n = length(unique(.data$value)), .groups = "drop")

  drop_items <- dplyr::filter(unique_df, .data$n < 2)$item
  cat("\nDropping the following items due to lack of variation:\n")
  cat(c("  *", paste(drop_items, collapse = "\n  * "), "\n"))
  if (is.na(binary_items)) {
    binary_items <-
      sort(dplyr::filter(unique_df, .data$n == 2 &
                        !.data$item %in% drop_items)$item)
  }
  if (is.na(ordinal_items)) {
    ordinal_items <-
      sort(dplyr::filter(unique_df, .data$n %in% 3:max_cats &
                        !.data$item %in% drop_items)$item)
  }

  metric_items <- setdiff(items, c(binary_items, ordinal_items, drop_items))
  metric_items <- sort(metric_items)
  if (make_indicator_for_zeros) {
    for (i in seq_along(metric_items)) {
      obs_i <- which(long_data$item == metric_items[i])
      is_zero <- long_data$value[obs_i] == 0
      if (any(is_zero) && 0 %in% range(long_data$value[obs_i])) {
        ## TODO: Add collinearity check
        long_data$value[obs_i[is_zero]] <- NA_real_
        zi_i <- stringr::str_c(metric_items[i], "_zi")
        newdat <- long_data[obs_i, ]
        newdat$value <- as.integer(!is_zero)
        newdat$item <- zi_i
        long_data <- dplyr::bind_rows(long_data, newdat)
        binary_items <- c(binary_items, zi_i)
      }
    }
  }

  cat("\nCategorizing the following items as binary:\n")
  cat(c("  *", paste(binary_items, collapse = "\n  * "), "\n"))
  cat("\nCategorizing the following items as ordinal:\n")
  cat(c("  *", paste(ordinal_items, collapse = "\n  * "), "\n"))
  cat("\nCategorizing the following items as metric:\n")
  cat(c("  *", paste(metric_items, collapse = "\n  * "), "\n"))

  # NOTE: use functions?
  binary_data <- long_data %>%
    dplyr::filter(.data$item %in% binary_items) %>%
    dplyr::mutate(ITEM = factor(.data$item, levels = binary_items)) %>%
    dplyr::group_by(.data$ITEM) %>%
    dplyr::mutate(yy = as.integer(ordered(.data$value)) - 1L) %>%
    dplyr::filter(!is.na(.data$yy)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$TIME, .data$ITEM, .data$UNIT)       # time must vary last

  ordinal_data <- long_data %>%
    dplyr::filter(.data$item %in% ordinal_items) %>%
    dplyr::mutate(ITEM = factor(.data$item, levels = ordinal_items)) %>%
    dplyr::group_by(.data$ITEM) %>%
    dplyr::mutate(yy = as.integer(ordered(.data$value))) %>%
    dplyr::filter(!is.na(.data$yy)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$TIME, .data$ITEM, .data$UNIT)       # time must vary last

  metric_data <- long_data %>%
    dplyr::filter(.data$item %in% metric_items) %>%
    dplyr::mutate(ITEM = factor(.data$item, levels = metric_items)) %>%
    dplyr::mutate(yy = .data$value) %>%
    dplyr::filter(!is.na(.data$yy)) %>%
    dplyr::arrange(.data$TIME, .data$ITEM, .data$UNIT)       # time must vary last

  if (standardize) {
    metric_data <- metric_data %>%
      dplyr::group_by(.data$ITEM) %>%
      dplyr::mutate(yy = stdize(.data$yy)) %>%
      dplyr::ungroup()
  }

  tob_b <- sapply(1:nlevels(long_data$TIME), function (t) {
    x <- as.integer(binary_data$TIME) == t
    if (any(x)) c(min(which(x)), max(which(x)))
    else c(0, 0)
  })
  tob_o <- sapply(1:nlevels(long_data$TIME), function (t) {
    x <- as.integer(ordinal_data$TIME) == t
    if (any(x)) c(min(which(x)), max(which(x)))
    else c(0, 0)
  })
  tob_m <- sapply(1:nlevels(long_data$TIME), function (t) {
    x <- as.integer(metric_data$TIME) == t
    if (any(x)) c(min(which(x)), max(which(x)))
    else c(0, 0)
  })

  stan_data <- list(
    J = nlevels(long_data$UNIT),
    T = nlevels(long_data$TIME),
    N_binary = nrow(binary_data),
    I_binary = nlevels(binary_data$ITEM),
    yy_binary = as.integer(binary_data$yy),
    ii_binary = as.integer(binary_data$ITEM),
    jj_binary = as.integer(binary_data$UNIT),
    tt_binary = as.integer(binary_data$TIME),
    tob_b = t(tob_b),
    N_ordinal = nrow(ordinal_data),
    I_ordinal = nlevels(ordinal_data$ITEM),
    K_ordinal = if (nrow(ordinal_data) > 0) max(ordinal_data$yy) else 1L,
    yy_ordinal = as.integer(ordinal_data$yy),
    ii_ordinal = as.integer(ordinal_data$ITEM),
    jj_ordinal = as.integer(ordinal_data$UNIT),
    tt_ordinal = as.integer(ordinal_data$TIME),
    tob_o = t(tob_o),
    N_metric = nrow(metric_data),
    I_metric = nlevels(metric_data$ITEM),
    yy_metric = metric_data$yy,
    ii_metric = as.integer(metric_data$ITEM),
    jj_metric = as.integer(metric_data$UNIT),
    tt_metric = as.integer(metric_data$TIME),
    tob_m = t(tob_m)
  )
  attr(stan_data, "unit_labels") <- levels(long_data$UNIT)
  attr(stan_data, "time_labels") <- levels(long_data$TIME)
  attr(stan_data, "binary_item_labels") <- levels(binary_data$ITEM)
  attr(stan_data, "ordinal_item_labels") <- levels(ordinal_data$ITEM)
  attr(stan_data, "metric_item_labels") <- levels(metric_data$ITEM)

  class(stan_data) <- c("dynIRT_data", class(stan_data))

  return(stan_data)
}
