
<!-- README.md is generated from README.Rmd. Please edit that file -->

# (temporary name) dynIRTtest

<!-- badges: start -->
<!-- badges: end -->

We will change the package name later. Please edit the `README.Rmd` file
and do not directly edit the `README.md` file.

## Installation

You can install the development version of dynIRTtest from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("tomoya-sasaki/dynIRTtest")
```

## TODOs

- [ ] Add descriptions
- [ ] Figure specifications

## Explanation

- How to use `dynIRTtest`

### Step 1: load data

- Load state-level public opinion data from 2020 and 2021

``` r
data("social_outcomes_2020_2021")
```

### Step 2: Reshape data

``` r
shaped_data <- shape_data(
    long_data = social_outcomes_2020_2021,
    unit_var = "st",
    time_var = "year",
    item_var = "outcome",
    value_var = "value",
    standardize = TRUE,
    periods_to_estimate = 2020:2021
)
```

### Step 3: Fit the model

``` r
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
```

### Step 4: Post analysis
