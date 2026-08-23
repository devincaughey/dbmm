# dbmm 0.0.1.9000

## Breaking changes

* `fit_mixfac()`'s `separate_eta` argument is renamed `smooth_eta`, **and its
  sense is reversed**: `separate_eta = TRUE` corresponds to
  `smooth_eta = FALSE`. Passing `separate_eta` still works but warns; passing
  both is an error.

* **The default dynamic behaviour has changed.** `separate_eta` defaulted to
  `TRUE`, so a default `fit_mixfac()` call estimated each period's factor
  scores independently — a static factor model fitted separately by period,
  with no dynamic component. `smooth_eta` defaults to `TRUE`, so the default
  is now a random walk in `eta`. Default fits are slower, produce different
  estimates, and include `sigma_eta_evol`, `Lcorr_eta`, and a freely
  estimated `Omega`. Pass `smooth_eta = FALSE` for the previous behaviour.

* **Default priors on the evolution standard deviations have changed.** The
  priors on `sigma_alpha_evol` and `sigma_eta_evol` are now half Student's t
  with location `0` and scale `0.1`, previously location `0.5` and scale
  `0.5`. Because these priors are truncated at zero, the location mattered
  less than its value suggests, but the prior median per-period standard
  deviation was near 0.6. On the whitened scale, where period-1 `eta` has
  cross-sectional standard deviation 1, that implied cumulative drift of
  roughly `0.6 * sqrt(T - 1)` — about 2.6 over 20 periods, i.e. effectively
  no smoothing. The new defaults shrink toward time-constant intercepts and
  factor scores while retaining heavy tails. Pass
  `mu_sigma_alpha_evol = 0.5, sd_sigma_alpha_evol = 0.5` (and likewise for
  `*_eta_evol`) to recover the previous behaviour. Priors on `sigma_metric`
  are unchanged.

* Several variables are now absent from the posterior draws under flag
  settings in which they do not enter the likelihood:

  * `sigma_eta_evol` and `Lcorr_eta` unless `T > 1` and `smooth_eta = TRUE`.
    `Omega` is fixed to the identity in that case.
  * `sigma_alpha_evol` unless `T > 1` and `constant_alpha = FALSE`. It is now
    indexed, i.e. `sigma_alpha_evol[1]`.
  * `z_alpha_trichot` and `z_alpha_ordinal` when `constant_alpha = TRUE`; and
    `z_alpha_binary` and `z_alpha_metric` then hold only one period of
    deviates.

  `sigma_eta_evol` and `Lcorr_eta` are also indexed, i.e.
  `sigma_eta_evol[1, d]` and `Lcorr_eta[1, i, j]`.

* `r_eta`, `r_Omega`, and `WW` are no longer returned. `r_eta` duplicated
  `eta` exactly; `r_Omega` was an intermediate; `WW` was retained only for a
  generated-quantities computation that has been removed. `WW`, `Lcorr_eta`,
  and `sigma_eta_evol` are now excluded by the default `drop` pattern of
  `extract_mixfac_draws()`.

* `dbmm_probit.stan` is renamed `mixfac_probit.stan`, matching the
  mixfac/modgirt naming used throughout the R code.

* `dbmm_logit.stan` is removed, and `fit_mixfac()` rejects `link = "logit"`
  until a logit variant is reintroduced. The file had received none of the
  changes below and would have failed on the data fields that `fit_mixfac()`
  now supplies.

* `dbmm_probit_alt.stan`, an alternative parameterization, is removed. It
  anchored whitening on the across-period mean of `eta`, which created an
  exact scale invariance between `z_eta` and
  `(sigma_eta_evol, Lcorr_eta)`, leaving the evolution parameters identified
  only by their priors. The pooled anchor it was reaching for is now
  available post hoc via `identify_mixfac(whiten = TRUE, ref_t = "mean")`,
  which is deterministic per draw and so carries no identification cost.

* `shape_mixfac()` no longer returns `tob_b`, `tob_t`, `tob_o`, or `tob_m`.
  These per-period index ranges were computed and passed to Stan but never
  used.

* In `modgirt_probit.stan`, the transition covariance `Omega` is fixed to the
  identity and its parameters are not estimated when `T == 1`.
  
## Improvements

* `summarize_mixfac()` honours `summary_functions` for `rvar` elements and no
  longer calls unexported posterior functions.
* Internal whitening and matrix-inversion helpers now operate on
  `posterior::draws_of()` arrays directly, preserving chain information for
  convergence diagnostics.

## Bug fixes

* **`Omega` was reported incorrectly.** It was computed as
  `quad_form(r_Omega, WW)`, rotating the innovation covariance by the
  period-1 whitening matrix — but the random walk is built directly on the
  whitened period-1 configuration and its innovations are never rotated by
  `WW`. `Omega` is now `L_eta L_eta'`. Values reported by earlier versions
  were wrong whenever `whiten_eta = TRUE` and `D > 1`.

* `identify_mixfac()` read `n_time`, the whitening anchor, and the label
  attributes from its `x` argument rather than from the draws and the fitted
  object, so the documented fitted-object input path did not work: with
  `whiten = TRUE` it errored, and otherwise it returned draws stripped of
  their unit, time, and item labels, which `label_mixfac()` then applied as
  `NULL` dimnames without complaint.

* `identify_mixfac()` computed the varimax rotation from the pre-whitening
  loadings and applied it to the whitened ones, so the rotation was not the
  varimax solution for the matrix it was applied to.

* `identify_mixfac()` read `n_time` and the whitening anchor from its `x`
  argument rather than from the extracted draws, so the documented
  fitted-object input path did not work.

* `sort_mixfac()` and `sign_mixfac()` referenced a `kappa` variable that
  `combine_types_mixfac()` does not produce (it emits `kappa_trichot` and
  `kappa_ordinal`), silently discarding the ordered thresholds. Both now
  carry them through.

* `sign_mixfac()` indexed `eta` as `x$eta[t, , drop = TRUE]`, supplying two
  indices for a three-dimensional `[T, J, D]` array.

* `sort_mixfac()` indexed `eta`, `lambda`, and `Omega` without
  `drop = FALSE`, collapsing them in single-factor models.

* `sort_mixfac()` and `sign_mixfac()` did not copy the unit, time, or item
  label attributes, and `sign_mixfac()` returned an unclassed object, so its
  result was not a `mixfac_comb` and could not be passed on. Both now also
  record the transformation applied, in the `"factor order"` and
  `"sign flips"` attributes, and both reject the long-format output of
  `label_mixfac(make_long = TRUE)`, which they cannot process.

* `sign_mixfac()` and `sign_modgirt()` tested `length(signs == 1)` rather
  than `length(signs) == 1`, so the guard never fired. `sign_mixfac()` now
  also validates that `signs` contains only `-1` and `1`.

* `combine_types_mixfac()` chose between its data-frame and rvar code paths
  by testing `all(sapply(x[idx], is.data.frame))` separately for `lambda`,
  `alpha`, and `kappa`. When an item type is absent the selection is empty,
  and `all()` of an empty set is `TRUE`. With no trichotomous or ordinal
  items the data-frame branch was taken on an rvar object, returning a plain
  list with a spurious empty `kappa` and — because attribute copying occurred
  only in the other branch — none of the label attributes.

* `combine_types_mixfac()` prefixed item names via
  `paste0(type, ": ", dimnames(v)[["item"]])`, which returns a length-1
  string when the dimnames are `NULL`, producing a dimnames length error on
  objects that had not been through `label_mixfac()`. Item indices are now
  used as a fallback.

* `shape_mixfac()`'s `binary_items`, `trichotomous_items`, and
  `ordinal_items` arguments were tested with `if (is.na(x))`, which errors on
  vectors of length greater than one in R >= 4.2, so the documented ability
  to specify item types manually was unusable for more than one item per
  type.

* `shape_mixfac()` printed item-categorization headers even when the
  corresponding set was empty, yielding a bare bullet.

* `shape_mixfac()` passed `.data$item` to `dplyr::summarise()`'s `.by`
  argument and bare column names to `dplyr::select()`, both of which are
  tidyselect contexts; these are now strings.

## Sampling efficiency

The changes below remove parameters that were declared and given priors but
did not enter the likelihood, so they were explored by HMC while merely
reproducing their priors. Marginal posteriors of the retained parameters are
unchanged.

* `sigma_eta_evol` and `Lcorr_eta` are no longer estimated unless the
  transition is identified, replacing a `lkj_corr_cholesky(20)` prior that
  approximated the same effect. Removes `D + D(D-1)/2` parameters.

* Intercept deviates are sized by whether they are used. With
  `constant_alpha = TRUE` this removes
  `(T-1)(I_binary + I_metric) + T(I_trichot + I_ordinal) + 1` parameters —
  often several hundred to a few thousand.

* `r_eta`, which duplicated `eta` in every branch, is no longer computed or
  stored, removing `T * J * D` values per draw and one array from the
  autodiff tape. `WW` is now local to `whiten()`.

* The unused functions `p2l_array()`, `num_matches()`, and `which_equal()`,
  and the unused `tob_*` data arrays, are removed from the Stan program.

## New features

* `fit_mixfac()` gains `gen_log_lik`, which adds a per-observation
  `log_lik` to the generated quantities for use with the **loo** package.
  Defaults to `FALSE`, since `log_lik` has one element per observation per
  draw. Note that because `eta[t, j]` is informed by the observations being
  held out, this is conditional leave-one-observation-out cross-validation;
  elevated Pareto k is expected for sparsely observed cells.

* `log_lik_index()` maps each position of `log_lik` to its item type, item,
  unit, period, and response value.

* `identify_mixfac()`'s `ref_t` argument gains a `"mean"` option, anchoring
  the normalization on each unit's average across periods rather than a
  single period. A pooled anchor does not propagate sampling noise from one
  cross-section through the whole series, and it puts periods on a common
  scale when `smooth_eta = FALSE`.

* Both fit functions now error informatively when no Stan file matches the
  requested `link`, rather than passing an empty path to **cmdstanr**.
  
## Reproducibility

* `identify_mixfac()` is now deterministic. The varimax rotation used a
  random start, so factor order and sign varied between calls on identical
  draws. Note that factor orientation is only fully pinned down by following
  `identify_mixfac()` with `sort_mixfac()` and `sign_mixfac()`.

## Documentation

* `@param white_eta` in `fit_mixfac()` corrected to `whiten_eta`.

* Documented defaults for `mu_sigma_*` and `sd_sigma_*` now match the code;
  several had drifted.

* `mu_sigma_metric` and `sd_sigma_metric` were described as priors for
  `sigma_alpha_evol`.

* Note that `link = "probit"` describes the scale on which coefficients are
  interpretable, not the likelihood: the model applies a logit approximation
  to the linear predictor and then uses `bernoulli_logit` and
  `ordered_logistic`.

## Internal

* Adds a **testthat** suite for the mixfac model, covering data preparation,
  the flag configurations affected by the parameter-sizing changes, draw
  dimensions, optional `log_lik`, the
  extract–identify–label–combine–sort–sign pipeline, and post-hoc whitening,
  including that every linear predictor and `lp__` are invariant to
  identification, sorting, and signing. Fixtures discretize four of the
  continuous outcomes in `social_outcomes_2020_2021` so that all four
  item-type code paths are exercised.

* Adds `as_draws_rvars_safe()`, which drops absent variables when
  constructing a `draws_rvars` object, and `copy_mixfac_attrs()`.

## Licensing

* The package is licensed under GPL-2. The repository previously contained an
  MIT `LICENSE` file, which conflicted with the `GPL (>= 3)` declaration in
  `DESCRIPTION`; the intended license is GPL-2.
