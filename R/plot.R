#' internal function to check `item_labels` argument
check_item_labels <- function(outcomes_labeled_item, item_labels) {
  original_labels <- sort(as.character(unique(outcomes_labeled_item)))

  if (length(original_labels) != length(item_labels)) {
    stop("The length of original labels and `item_labels` are different.")
  } else if (sum(is.na(item_labels))) {
    stop("New labels should not contain `NA`.")
  } else if (sum(original_labels != sort(names(item_labels))) > 0) {
    stop("Element names of `item_labels` must match original labels")
  }
}


#' test function
#' @param outcomes_labeled  A `dynIRT_labeled` object
#' @param xtitle (string) x-axis label of the plot
#' @param ytitle (string) y-axis label of the plot
#' @param item_type (string) Should "binary", "ordinal", or "metric" loadings be
#' used to visualize? Default is "metric".
#' @param item_labels (string) Named string vector where each element represents the new labels and the names of the elements correspond to the original label in `outcomes_labeled`. Default is NULL.
#' @param maintitle (string) Title of the plot
#'
#' @return A plot showing the estimated intercept
#'
#' @import magrittr ggplot2
#' @importFrom rlang .data
#' @export
plot_intercept_test <- function(outcomes_labeled,
                          xtitle = NULL,
                          ytitle = "Estimated Intercept",
                          maintitle = NULL,
                          item_type = c("metric", "binary", "ordinal"),
                          item_labels = NULL) {


  if (item_type == "metric") {
    dat <- outcomes_labeled$alpha_metric |> dplyr::mutate(type = "metric")
  } else if (item_type == "binary") {
    dat <- outcomes_labeled$alpha_binary |> dplyr::mutate(type = "binary")
  } else if (item_type == "ordinal") {
    dat <- outcomes_labeled$alpha_ordinal |> dplyr::mutate(type = "ordinal")
  } else {
    stop("Invalid `item_type` argument")
  }

  if (!is.null(item_labels)) {
    check_arg_type(item_labels, "character")
    check_item_labels(dat$ITEM, item_labels)
  }

  dat %>%
    dplyr::group_by(.data$TIME, .data$ITEM) %>%
    dplyr::summarise(
      est = mean.default(.data$value),
      err = stats::sd(.data$value),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      year = as.integer(as.character(.data$TIME))
    ) %>%
    {if (!is.null(item_labels)) dplyr::mutate(., ITEM = dplyr::recode(.data$ITEM, !!!item_labels)) else . } %>%
    dplyr::mutate(
      ITEM = stats::reorder(.data$ITEM, .data$est, FUN = stats::sd)
    ) %>%
    ggplot(aes(x = .data$year, y = .data$est)) +
      # facet_wrap(~.data$ITEM, ncol = 5) +
      facet_wrap(~.data$ITEM) +
      geom_line() +
      geom_ribbon(
        aes(ymin = .data$est - 1.96 * .data$err,
            ymax = .data$est + 1.96 * .data$err),
            color = NA, alpha = 1/4
      ) +
      labs(
        title = maintitle,
        y = ytitle,
        x = xtitle
        # color = NULL,
        # fill = NULL
      ) -> p

  class(p) <- c("dynIRTtest_viz", class(p))
  p
}


#' plot item intercepts
#'
#' @param outcomes_labeled  A `dynIRT_labeled` object
#' @param xtitle (string) x-axis label of the plot
#' @param ytitle (string) y-axis label of the plot
#' @param maintitle (string) Title of the plot
#' @param item_labels (string) Named string vector where each element represents the new labels and the names of the elements correspond to the original label in `outcomes_labeled`. Default is NULL.
#'
#' @return A plot showing the estimated intercept
#'
#' @import magrittr ggplot2
#' @importFrom rlang .data
#'
#' @export
plot_intercept <- function(outcomes_labeled,
                          xtitle = NULL,
                          ytitle = "Estimated Intercept",
                          maintitle = NULL,
                          item_labels = NULL) {

  ## function to check if item_labels is properly

  if (!is.null(item_labels)) {
    check_arg_type(item_labels, "character")
    check_item_labels(outcomes_labeled$alpha_metric$ITEM, item_labels)
  }

  outcomes_labeled$alpha_metric %>%
    dplyr::group_by(.data$TIME, .data$ITEM) %>%
    dplyr::summarise(
      est = mean.default(.data$value),
      err = stats::sd(.data$value),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      year = as.integer(as.character(.data$TIME))
    ) %>%
    {if (!is.null(item_labels)) dplyr::mutate(., ITEM = dplyr::recode(.data$ITEM, !!!item_labels)) else . } %>%
    dplyr::mutate(
      ITEM = stats::reorder(.data$ITEM, .data$est, FUN = stats::sd)
    ) %>%
    ggplot(aes(x = .data$year, y = .data$est)) +
      # facet_wrap(~.data$ITEM, ncol = 5) +
      facet_wrap(~.data$ITEM) +
      geom_line() +
      geom_ribbon(
        aes(ymin = .data$est - 1.96 * .data$err,
            ymax = .data$est + 1.96 * .data$err),
            color = NA, alpha = 1/4
      ) +
      labs(
        title = maintitle,
        y = ytitle,
        x = xtitle
        # color = NULL,
        # fill = NULL
      ) -> p

  class(p) <- c("dynIRTtest_viz", class(p))

  p
}


#' plot item loadings
#'
#' @param outcomes_labeled A `dynIRT_labeled` object
#' @param xtitle (string) x-axis label of the plot
#' @param ytitle (string) y-axis label of the plot
#' @param maintitle (string) Title of the plot. Default is ``Item Loadings''
#' @param item_labels (string) Named string vector where each element represents the new labels and the names of the elements correspond to the original label in `outcomes_labeled`. Default is NULL.
#'
#' @return A plot showing item loadings
#'
#' @import magrittr ggplot2
#' @importFrom rlang .data
#'
#' @export
plot_loadings <- function(outcomes_labeled,
                          xtitle = NULL,
                          ytitle = NULL,
                          maintitle = "Item Loadings",
                          item_labels = NULL) {

  if (!is.null(item_labels)) {
    check_arg_type(item_labels, "character")
    check_item_labels(outcomes_labeled$lambda_metric$ITEM, item_labels)
  }

  outcomes_labeled$lambda_metric %>%
    {if (!is.null(item_labels)) dplyr::mutate(., ITEM = dplyr::recode(.data$ITEM, !!!item_labels)) else . } %>%
    dplyr::group_by(.data$ITEM, .data$dim) %>%
    dplyr::summarise(est = mean.default(.data$value),
                     err = stats::sd(.data$value)) %>%
    tidyr::pivot_wider(id_cols = "ITEM", names_from = "dim",
                      values_from = c("est", "err")) %>%
    ggplot() +
      aes(x = .data$est_1, y = .data$est_2, label = .data$ITEM) +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_point() +
      geom_linerange(
        aes(xmin = .data$est_1 - 1.96 * .data$err_1,
            xmax = .data$est_1 + 1.96 * .data$err_1),
            alpha = 1/4, linewidth = 2
      ) +
      geom_linerange(
        aes(ymin = .data$est_2 - 1.96 * .data$err_2,
            ymax = .data$est_2 + 1.96 * .data$err_2),
            alpha = 1/4, linewidth = 2
      ) +
      ggrepel::geom_text_repel() +
      labs(
          title = maintitle,
          x = xtitle,
          y = ytitle
      ) +
      coord_fixed() -> p

  class(p) <- c("dynIRTtest_viz", class(p))
  p
}


#' @import magrittr
#' @importFrom rlang .data
create_factor_scores <- function(outcomes_labeled) {
  eta_ave <- outcomes_labeled$eta %>%
    dplyr::group_by(.data$UNIT, .data$dim, .data$.draw) %>%
    dplyr::summarise(ave = mean.default(.data$value)) %>%
    dplyr::summarise(est = mean.default(.data$ave),
                     err = stats::sd(.data$ave),
                     .groups = "drop") %>%
    dplyr::mutate(DIMENSION = paste("Dimension", .data$dim)) %>%
    tidyr::pivot_wider(id_cols = "UNIT", names_from = "dim",
                      values_from = c("est", "err"))

  eta_ave
}


#' plot average factor scores across different units
#'
#' @param outcomes_labeled  A `dynIRT_labeled` object
#' @param xtitle (string) x-axis label of the plot
#' @param ytitle (string) y-axis label of the plot
#' @param maintitle (string) Title of the plot. Default is "Agerage Factor Scores".
#' @param unit_labels (string) Named string vector where each element represents the new unit labels and the names of the elements correspond to the original unit labels in `outcomes_labeled`. Default is NULL.
#'
#' @return A plot showing average factor scores
#'
#' @import magrittr ggplot2
#' @importFrom rlang .data
#'
#' @export
plot_scores_ave <- function(outcomes_labeled,
                            xtitle = NULL,
                            ytitle = NULL,
                            maintitle = "Average Factor Scores",
                            unit_labels = NULL)
{
  eta_ave <- create_factor_scores(outcomes_labeled = outcomes_labeled)

  if (!is.null(unit_labels)) {
    check_arg_type(unit_labels, "character")
    check_item_labels(eta_ave$UNIT, unit_labels)
  }

  eta_ave %>%
    {if (!is.null(unit_labels)) dplyr::mutate(., UNIT = dplyr::recode(.data$UNIT, !!!unit_labels)) else . } %>%
    ggplot() +
      aes(x = .data$est_1, y = .data$est_2, label = .data$UNIT) +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_linerange(
        aes(xmin = .data$est_1 - 1.96 * .data$err_1,
            xmax = .data$est_1 + 1.96 * .data$err_1),
        alpha = 1/4, linewidth = 2
      ) +
      geom_linerange(
        aes(ymin = .data$est_2 - 1.96 * .data$err_2,
            ymax = .data$est_2 + 1.96 * .data$err_2),
        alpha = 1/4, linewidth = 2
      ) +
      geom_point() +
      ggrepel::geom_text_repel()  +
      labs(
        title = maintitle,
        x = xtitle,
        y = ytitle
      ) +
      coord_fixed() -> p

  class(p) <- c("dynIRTtest_viz", class(p))
  p
}


#' Plot time series factor scores. Currently only support 2-dimensional plot.
#'
#' @param outcomes_labeled  A `dynIRT_labeled` object
#' @param xtitle (string) x-axis label of the plot
#' @param ytitle (string) y-axis label of the plot
#' @param maintitle (string) Title of the plot
#' @param unit_labels (string) Named string vector where each element represents the new unit labels and the names of the elements correspond to the original unit labels in `outcomes_labeled`. Default is NULL.
#'
#' @return A plot plot showing time trend
#'
#' @import magrittr ggplot2
#' @importFrom rlang .data
#'
#' @export
plot_scores_timetrend <- function(outcomes_labeled,
                                  xtitle = "Dimension 1",
                                  ytitle = "Dimension 2",
                                  maintitle = NULL,
                                  unit_labels = NULL)
{

  if (!is.null(unit_labels)) {
    check_arg_type(unit_labels, "character")
    check_item_labels(outcomes_labeled$eta$UNIT, unit_labels)
  }

  outcomes_labeled$eta %>%
    dplyr::group_by(.data$TIME, .data$UNIT, .data$dim) %>%
    dplyr::summarise(est = mean.default(.data$value),
                     err = stats::sd(.data$value),
                     .groups = "drop") %>%
    {if (!is.null(unit_labels)) dplyr::mutate(., UNIT = dplyr::recode(.data$UNIT, !!!unit_labels)) else . } %>%
    dplyr::mutate(
      DIMENSION = dplyr::if_else(.data$dim == 1, xtitle, ytitle),
      UNIT = stats::reorder(.data$UNIT, -.data$est, FUN = mean.default),
      year = as.integer(as.character(.data$TIME))
    ) %>%
    ggplot() +
      # facet_wrap(~ .data$faUNIT, ncol = 5) +
      facet_wrap(~ .data$UNIT) +
      aes(x = .data$year, y = .data$est,
          color = .data$DIMENSION, fill = .data$DIMENSION) +
      geom_ribbon(
        aes(ymin = .data$est - 1.96 * .data$err,
            ymax = .data$est + 1.96 * .data$err),
        color = NA, alpha = 1/4
      ) +
      geom_line() +
      # scale_x_continuous(breaks = xbreaks) +
      scale_color_brewer(type = "qual") +
      scale_fill_brewer(type = "qual") +
      # coord_cartesian(ylim = c(-3.2, 3.2)) +
      labs(
        title = maintitle,
        y = ytitle,
        x = xtitle
        # color = NULL,
        # fill = NULL
      ) +
      theme(legend.position = "bottom") -> p

  class(p) <- c("dynIRTtest_viz", class(p))
  p
}
