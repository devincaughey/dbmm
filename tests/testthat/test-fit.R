data(social_outcomes_2020_2021)

shaped_data <- shape_data(
    long_data = social_outcomes_2020_2021,
    unit_var = "st",
    time_var = "year",
    item_var = "outcome",
    value_var = "value",
    standardize = TRUE,
    periods_to_estimate = 2020:2021
)

fitted <- fit(
    data = shaped_data,
    lambda_zeros = data.frame(
        item = c("x_cps_spm_poverty"),
        dim = c(2)
    ),
    n_dim = 2,
    chains = 2,
    parallelize_within_chains = TRUE,
    threads_per_chain = 2,
    constant_alpha = FALSE,
    separate_eta = TRUE,
    init_kappa = FALSE,
    force_recompile = FALSE,
    iter_warmup = 50,
    iter_sampling = 50,
    adapt_delta = .9,
    refresh = 10,
    seed = 123
)

fitted_draws <- extract_draws(fitted)

test_that("Extracted draws", {
  expect_equal(names(fitted_draws)[1], "lp__")
  expect_equal(nrow(fitted_draws), 100)
  }
)

test_that("Check convergence", {
  expect_message(check_convergence(fitted_draws),
                "The following variables have not converged")
  }
)


identified <- identify_draws(fitted_draws, rotate = TRUE)

test_that("Identified rotation", {
  expect_equal(length(identified$rotmats), 4)
  expect_s3_class(identified, "dynIRT_identified")
  }
)

labeled <- label_draws(identified)

test_that("Labeled draws", {
  expect_s3_class(labeled, "dynIRT_labeled")
  }
)

test_that("Figures", {
  p <- plot_intercept(labeled)
  expect_s3_class(p, "dynIRTtest_viz")

  create_metric_label <- function(outcomes_labeled) {
    metric_labels <- attr(outcomes_labeled, "metric_item_labels") |>
      stringr::str_remove("^[xz]_") |>
      stringr::str_remove("cps_") |>
      stringr::str_remove("spm_") |>
      stringr::str_remove("_guttmacher_occurence") |>
      stringr::str_replace_all("per_capita", "pc") |>
      stringr::str_replace_all("_", " ")
    names(metric_labels) <- attr(outcomes_labeled, "metric_item_labels")

    return(metric_labels)
  }

  p <- plot_intercept(labeled, item_labels = create_metric_label(labeled))
  expect_s3_class(p, "dynIRTtest_viz")

  p <- plot_loadings(labeled)
  expect_s3_class(p, "dynIRTtest_viz")

  p <- plot_loadings(labeled, item_labels = create_metric_label(labeled))
  expect_s3_class(p, "dynIRTtest_viz")

  # check `check_item_labels` works
  wrong_labels <- create_metric_label(labeled)
  names(wrong_labels)[1] <- "a"

  expect_error(plot_intercept(labeled, item_labels = wrong_labels),
              regexp = "Element names of `item_labels")

  p <- plot_scores_ave(labeled)
  expect_s3_class(p, "dynIRTtest_viz")

  p <- plot_scores_timetrend(labeled)
  expect_s3_class(p, "dynIRTtest_viz")

  new_unit_names <- state.name
  names(new_unit_names) <- state.abb

  p <- plot_scores_ave(labeled, unit_labels = new_unit_names)
  expect_s3_class(p, "dynIRTtest_viz")

  p <- plot_scores_timetrend(labeled, unit_labels = new_unit_names)
  expect_s3_class(p, "dynIRTtest_viz")

  }
)
