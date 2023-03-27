#' @import magrittr
create_metric_label <- function(outcomes_labeled) {
  metric_labels <- attr(outcomes_labeled, "metric_item_labels") %>%
    stringr::str_remove("^[xz]_") %>%
    stringr::str_remove("cps_") %>%
    stringr::str_remove("spm_") %>%
    stringr::str_remove("_guttmacher_occurence") %>%
    stringr::str_replace_all("per_capita", "pc") %>%
    stringr::str_replace_all("_", " ")
  names(metric_labels) <- attr(outcomes_labeled, "metric_item_labels")

  return(metric_labels)
}



#' plot item intercepts
#'
#' @param outcomes_labeled
#'
#' @return plot
#'
#' @import magrittr ggplot2
#' @importFrom rlang .data
#'
#' @export
plot_intercept <- function(outcomes_labeled) {
  metric_labels <- create_metric_label(outcomes_labeled = outcomes_labeled)
  outcomes_labeled$alpha_metric %>%
    dplyr::group_by(.data$TIME, .data$ITEM) %>%
    dplyr::summarise(
      est = mean(.data$value),
      err = stats::sd(.data$value),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      ITEM = dplyr::recode(.data$ITEM, !!!metric_labels),
      # ITEM = stats::reorder(.data$ITEM, .data$est, FUN = stats::sd),
      year = as.integer(as.character(.data$TIME))
    ) %>%
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
      # scale_x_continuous(breaks = seq(1960, 2020, 20)) +
      labs(
        title = "Item Intercepts over Time",
        y = "Estimated Intercept",
        x = NULL,
        color = NULL,
        fill = NULL
      ) -> p

  p
}


#' plot item loadings
#'
#' @param outcomes_labeled
#'
#' @return plot
#'
#' @import magrittr ggplot2
#' @importFrom rlang .data
#'
#' @export
plot_loadings <- function(outcomes_labeled) {
  metric_labels <- create_metric_label(outcomes_labeled = outcomes_labeled)

  outcomes_labeled$lambda_metric %>%
    dplyr::mutate(ITEM = dplyr::recode(.data$ITEM, !!!metric_labels)) %>%
    dplyr::group_by(.data$ITEM, .data$dim) %>%
    dplyr::summarise(est = mean(.data$value), err = stats::sd(.data$value)) %>%
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
      # labs(
      #     title = "Item Loadings",
      #     x = "Dimension 1 (Valence)",
      #     y = "Dimension 2 (Spatial)"
      # ) +
      coord_fixed() -> p

  p
}


#' @import magrittr
#' @importFrom rlang .data
create_factor_scores <- function(outcomes_labeled) {
  eta_ave <- outcomes_labeled$eta %>%
    dplyr::group_by(.data$UNIT, .data$dim, .data$.draw) %>%
    dplyr::summarise(ave = mean(.data$value)) %>%
    dplyr::summarise(est = mean(.data$ave), err = stats::sd(.data$ave), .groups = "drop") %>%
    dplyr::mutate(DIMENSION = paste("Dimension", .data$dim)) %>%
    tidyr::pivot_wider(id_cols = "UNIT", names_from = "dim",
                      values_from = c("est", "err"))

  eta_ave
}


#' plot average factor scores
#'
#' @param outcomes_labeled
#'
#' @return plot
#'
#' @import magrittr ggplot2
#' @importFrom rlang .data
#'
#' @export
plot_scores_ave <- function(outcomes_labeled,
                            xdim = "Dimension 1 (Valence)",
                            ydim = "Dimension 2 (Spatial)")
{
  eta_ave <- create_factor_scores(outcomes_labeled = outcomes_labeled)

  eta_ave %>%
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
        # title = "Average State Outcome Scores",
        x = xdim,
        y = ydim
      ) +
      coord_fixed() -> p

  p
}


#' plot time series factor scores
#'
#' @param outcomes_labeled
#'
#' @return plot
#'
#' @import magrittr ggplot2
#' @importFrom rlang .data
#'
#' @export
plot_scores_timetrend <- function(outcomes_labeled,
                                  xdim = "Dimension 1 (Valence)",
                                  ydim = "Dimension 2 (Spatial)")
{
  outcomes_labeled$eta %>%
    dplyr::group_by(.data$TIME, .data$UNIT, .data$dim) %>%
    dplyr::summarise(est = mean(.data$value),
                    err = stats::sd(.data$value),
                    .groups = "drop") %>%
    dplyr::mutate(
      DIMENSION = dplyr::if_else(.data$dim == 1, xdim,
                                ydim),
      UNIT = stats::reorder(.data$UNIT, -.data$est, FUN = mean),
      year = as.integer(as.character(.data$TIME))
    ) %>%
    ggplot() +
      # facet_wrap(~ .data$faUNIT, ncol = 5) +
      facet_wrap(~ .data$faUNIT) +
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
      # labs(
      #   title = "State Outcome Scores over Time",
      #   y = "Estimated Factor Score",
      #   x = NULL,
      #   color = NULL,
      #   fill = NULL
      # ) +
      theme(legend.position = "bottom") -> p

  p
}