data(social_outcomes_2021)

shaped_data <- shape_data(
    long_data = social_outcomes_2021,
    unit_var = "st",
    time_var = "year",
    item_var = "outcome",
    value_var = "value",
    standardize = TRUE,
    periods_to_estimate = 2021
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

fitted_draws <- extract_draws(fitted$fit)

test_that("Main extracted draws", {
  expect_equal(names(fitted_draws)[1], "lp__")
  expect_equal(nrow(fitted_draws), 100)
  }
)

identified <- identify_draws(fitted_draws, rotate = TRUE)

test_that("Identified rotation", {
  expect_equal(length(identified$rotmats), 4)
  }
)