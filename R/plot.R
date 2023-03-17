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
#'
#' @export
plot_intercept <- function(outcomes_labeled) {
  metric_labels <- create_metric_label(outcomes_labeled = outcomes_labeled)
  outcomes_labeled$alpha_metric %>%
    dplyr::group_by(TIME, ITEM) %>%
    dplyr::summarise(est = mean(value), err = sd(value), .groups = "drop") %>%
    dplyr::mutate(
        ITEM = dplyr::recode(ITEM, !!!metric_labels),
        ITEM = reorder(ITEM, est, FUN = sd),
        year = as.integer(as.character(TIME))
    ) %>%
    ggplot(aes(x = year, y = est)) +
      facet_wrap(~ITEM, ncol = 5) +
      geom_line() +
      geom_ribbon(
        aes(ymin = est - 1.96 * err, ymax = est + 1.96 * err),
            color = NA, alpha = 1/4
    ) +
    scale_x_continuous(breaks = seq(1960, 2020, 20)) +
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
#'
#' @export
plot_loadings <- function(outcomes_labeled) {
  metric_labels <- create_metric_label(outcomes_labeled = outcomes_labeled)

  outcomes_labeled$lambda_metric %>%
    dplyr::mutate(ITEM = dplyr::recode(ITEM, !!!metric_labels)) %>%
    dplyr::group_by(ITEM, dim) %>%
    dplyr::summarise(est = mean(value), err = sd(value)) %>%
    tidyr::pivot_wider(id_cols = "ITEM", names_from = "dim",
                      values_from = c("est", "err")) %>%
    ggplot() +
      aes(x = est_1, y = est_2, label = ITEM) +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_point() +
      geom_linerange(
        aes(xmin = est_1 - 1.96*err_1, xmax = est_1 + 1.96*err_1),
            alpha = 1/4, linewidth = 2
      ) +
      geom_linerange(
        aes(ymin = est_2 - 1.96*err_2, ymax = est_2 + 1.96*err_2),
            alpha = 1/4, linewidth = 2
      ) +
      ggrepel::geom_text_repel() +
      labs(
          title = "Item Loadings",
          x = "Dimension 1 (Valence)",
          y = "Dimension 2 (Spatial)"
      ) +
      coord_fixed() -> p

  p
}


#' @import magrittr
create_factor_scores <- function(outcomes_labeled) {
  eta_ave <- outcomes_labeled$eta %>%
    dplyr::group_by(UNIT, dim, .draw) %>%
    dplyr::summarise(ave = mean(value)) %>%
    dplyr::summarise(est = mean(ave), err = sd(ave), .groups = "drop") %>%
    dplyr::mutate(DIMENSION = paste("Dimension", dim)) %>%
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
#'
#' @export
plot_scores_ave <- function(outcomes_labeled) {
  eta_ave <- create_factor_scores(outcomes_labeled = outcomes_labeled)

  eta_ave %>%
    ggplot() +
      aes(x = est_1, y = est_2, label = UNIT) +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_linerange(
        aes(xmin = est_1 - 1.96*err_1, xmax = est_1 + 1.96*err_1),
            alpha = 1/4, linewidth = 2
      ) +
      geom_linerange(
        aes(ymin = est_2 - 1.96*err_2, ymax = est_2 + 1.96*err_2),
            alpha = 1/4, linewidth = 2
      ) +
      geom_point() +
      ggrepel::geom_text_repel()  +
      labs(
        title = "Average State Outcome Scores",
        x = "Dimension 1 (Valence)",
        y = "Dimension 2 (Spatial)"
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
#'
#' @export
plot_scores_timetrend <- function(outcomes_labeled) {
  outcomes_labeled$eta %>%
    dplyr::group_by(TIME, UNIT, dim) %>%
    dplyr::summarise(est = mean(value), err = sd(value), .groups = "drop") %>%
    dplyr::mutate(
      DIMENSION = dplyr::if_else(dim == 1, "Dimension 1 (Valence)",
                                "Dimension 2 (Spatial)"),
      UNIT = reorder(UNIT, -est, FUN = mean),
      year = as.integer(as.character(TIME))
    ) %>%
    ggplot() +
      facet_wrap(~UNIT, ncol = 5) +
      aes(x = year, y = est, color = DIMENSION, fill = DIMENSION) +
      geom_ribbon(
        aes(ymin = est - 1.96*err, ymax = est + 1.96*err),
            color = NA, alpha = 1/4
      ) +
      geom_line() +
      scale_x_continuous(breaks = seq(1960, 2020, 20)) +
      scale_color_brewer(type = "qual") +
      scale_fill_brewer(type = "qual") +
      coord_cartesian(ylim = c(-3.2, 3.2)) +
      labs(
        title = "State Outcome Scores over Time",
        y = "Estimated Factor Score",
        x = NULL,
        color = NULL,
        fill = NULL
      ) +
      theme(legend.position = "bottom") -> p

  p
}