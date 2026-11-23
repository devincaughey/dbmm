functions {
  vector p2l_vector (vector x) { // coverts vector from probit to logit scale
    return x .* (1.5976 + 0.07056 * x .* x);
  }
  real p2l (real x) { // converts scalar from probit to logit scale
    return x * (1.5976 + 0.07056 * x * x);
  }
  /* Return whitening matrix for demeaned matrix DM */
  matrix make_whiten_matrix(matrix DM) {
    int K = cols(DM);
    real eps = 1e-8;
    matrix[K, K] SS = crossprod(DM) ./ (rows(DM) - 1.0);
    // small nugget for numerical stability
    SS += eps * diag_matrix(rep_vector(1.0, K));
    matrix[K, K] L = cholesky_decompose(SS);     // SS = L * L'
    // inv(L): solve L * X = I
    matrix[K, K] invL = mdivide_left_tri_low(L, diag_matrix(rep_vector(1.0, K)));
    matrix[K, K] WW = invL';                      // WW = inv(L)'
    return WW;
  }
  // whiten: demean columns of XX then multiply by whitening matrix (DM * WW)
  matrix whiten(matrix XX) {
    int R = rows(XX);
    int C = cols(XX);
    matrix[R, C] DM;
    for (c in 1:C) DM[:, c] = XX[:, c] - mean(XX[:, c]);
    matrix[C, C] WW = make_whiten_matrix(DM);
    return DM * WW;
  }

  /* Linear predictors for a block of observations: alpha_flat holds the
     intercepts flattened period by period, LL the item loadings, and EE the
     factors stacked period by period (row (t - 1) * J + j). ti and tj give
     each observation's position in alpha_flat and its row in EE. */
  vector linear_predictor(vector alpha_flat, matrix LL, matrix EE,
                          array[] int ti, array[] int ii, array[] int tj) {
    return alpha_flat[ti] + rows_dot_product(LL[ii], EE[tj]);
  }

  real bprobit_partial_sum_lpmf(array[] int yy_b_slice,
                               int start,
                               int end,
                               vector alpha_b,
                               matrix lambda_b,
                               matrix eta_mat,
                               array[] int ti_b,
                               array[] int ii_b,
                               array[] int tj_b) {
    vector[end - start + 1] nu_slice =
      linear_predictor(alpha_b, lambda_b, eta_mat, ti_b[start:end],
                       ii_b[start:end], tj_b[start:end]);
    return bernoulli_logit_lupmf(yy_b_slice | p2l_vector(nu_slice));
  }

  real oprobit_partial_sum_lpmf(array[] int yy_o_slice,
                               int start,
                               int end,
                               vector alpha_o,
                               matrix lambda_o,
                               matrix eta_mat,
                               array[] int ti_o,
                               array[] int ii_o,
                               array[] int tj_o,
                               array[] vector kappa) {
    vector[end - start + 1] nu_slice =
      linear_predictor(alpha_o, lambda_o, eta_mat, ti_o[start:end],
                       ii_o[start:end], tj_o[start:end]);
    return ordered_logistic_lupmf(yy_o_slice | p2l_vector(nu_slice),
                                  kappa[ii_o[start:end]]);
  }

  real normal_partial_sum_lpdf(array[] real yy_m_slice,
                               int start,
                               int end,
                               vector alpha_m,
                               matrix lambda_m,
                               matrix eta_mat,
                               array[] int ti_m,
                               array[] int ii_m,
                               array[] int tj_m,
                               vector sigma) {
    vector[end - start + 1] nu_slice =
      linear_predictor(alpha_m, lambda_m, eta_mat, ti_m[start:end],
                       ii_m[start:end], tj_m[start:end]);
    return normal_lupdf(yy_m_slice | nu_slice, sigma[ii_m[start:end]]);
  }
}
data {
  int<lower=0,upper=1> parallelize;    /* parallelize within chains? */
  int<lower=0,upper=1> constant_alpha; /* keep alphas constant? */
  int<lower=0,upper=1> smooth_eta;     /* smooth eta across periods? */
  int<lower=0,upper=1> whiten_eta;     /* whiten eta */
  int<lower=0,upper=1> gen_log_lik; /* compute per-observation log_lik? */
  int<lower=1> D;                      /* number of latent dimensions */
  int<lower=1> J;                      /* number of units */
  int<lower=1> T;                      /* number of time periods */
  // Binary data //
  int<lower=0> N_binary;                          /* number of observations */
  int<lower=0> I_binary;                          /* number of items */
  array[N_binary] int<lower=0,upper=1> yy_binary; /* outcomes */
  array[N_binary] int<lower=1> ii_binary;         /* item indicator */
  array[N_binary] int<lower=1> jj_binary;         /* unit indicator */
  array[N_binary] int<lower=1> tt_binary;         /* time indicator */
  matrix[I_binary, D] nonzero_binary;		  /* nonzero loadings */
  // Trichotomous data //
  int<lower=0> N_trichot;                   /* number of observations */
  int<lower=0> I_trichot;                   /* number of items */
  array[N_trichot] int<lower=1> yy_trichot; /* outcomes */
  array[N_trichot] int<lower=1> ii_trichot; /* item indicator */
  array[N_trichot] int<lower=1> jj_trichot; /* unit indicator */
  array[N_trichot] int<lower=1> tt_trichot; /* period indicator */
  matrix[I_trichot, D] nonzero_trichot;     /* nonzero loadings */
  // Ordinal data //
  int<lower=0> N_ordinal;                   /* number of observations */
  int<lower=0> I_ordinal;                   /* number of items */
  int<lower=1> K_ordinal;                   /* max response categories */
  array[N_ordinal] int<lower=1> yy_ordinal; /* outcomes */
  array[N_ordinal] int<lower=1> ii_ordinal; /* item indicator */
  array[N_ordinal] int<lower=1> jj_ordinal; /* unit indicator */
  array[N_ordinal] int<lower=1> tt_ordinal; /* period indicator */
  matrix[I_ordinal, D] nonzero_ordinal;     /* nonzero loadings */
  // Metric data //
  int<lower=0> N_metric;                    /* number of observations */
  int<lower=0> I_metric;                    /* number of items */
  vector[N_metric] yy_metric;               /* outcomes */
  array[N_metric] int<lower=1> ii_metric;   /* item indicator */
  array[N_metric] int<lower=1> jj_metric;   /* unit indicator */
  array[N_metric] int<lower=1> tt_metric;   /* period indicator */
  matrix[I_metric, D] nonzero_metric;       /* nonzero loadings */
  // Priors //
  real<lower=0> df_sigma_metric;
  real<lower=0> df_sigma_alpha_evol;
  real<lower=0> df_sigma_eta_evol;
  real<lower=0> mu_sigma_metric;
  real<lower=0> mu_sigma_alpha_evol;
  real<lower=0> mu_sigma_eta_evol;
  real<lower=0> sd_sigma_metric;
  real<lower=0> sd_sigma_alpha_evol;
  real<lower=0> sd_sigma_eta_evol;
  // Shrinkage prior on each factor's loadings //
  int<lower=0, upper=2> loading_prior; /* 0 = normal, 1 = reg. horseshoe, 2 = lasso */
  real<lower=0> tau0;                  /* horseshoe global scale */
  real<lower=0> slab_scale;            /* horseshoe slab scale */
  real<lower=0> slab_df;               /* horseshoe slab degrees of freedom */
}
transformed data {
  /* Shrinkage parameters have zero length unless their prior is chosen, so
     fits with the default normal prior are unchanged. */
  int<lower=0, upper=1> est_hs = loading_prior == 1;
  int<lower=0, upper=1> est_lasso = loading_prior == 2;
  int<lower=0, upper=1> est_shrink = loading_prior > 0;
  /* The evolution covariance is identified only if an actual random walk in
     eta is used: T > 1 and periods are not estimated separately. Otherwise
     sigma_eta_evol / Lcorr_eta would be sampled from their priors alone, so
     they are given zero length and Omega is fixed to I_D. */
  int<lower=0, upper=1> est_Omega = (T > 1) && (smooth_eta == 1);
  /* Likewise, the intercept evolution SD enters the model only if the
     intercepts actually follow a random walk. */
  int<lower=0, upper=1> est_sigma_alpha = (T > 1) && (constant_alpha == 0);
  /* Periods of free intercept deviates. With constant_alpha == 1, later
     periods copy period 1, so only one period of deviates is needed. */
  int<lower=1> T_alpha = constant_alpha == 1 ? 1 : T;
  /* Trichotomous and ordinal intercepts are pinned to 0 (their ordered
     thresholds absorb the item level), so they need deviates only when the
     intercepts are allowed to drift. */
  int<lower=0> T_alpha_k = constant_alpha == 1 ? 0 : T;
  /* Each observation's row in the stacked (T * J) x D factor matrix (tj) and
     the position of its intercept in the flattened [T, I] intercepts (ti),
     so linear predictors can be computed in one vectorized step. */
  array[N_binary] int tj_binary;
  array[N_binary] int ti_binary;
  array[N_trichot] int tj_trichot;
  array[N_trichot] int ti_trichot;
  array[N_ordinal] int tj_ordinal;
  array[N_ordinal] int ti_ordinal;
  array[N_metric] int tj_metric;
  array[N_metric] int ti_metric;
  for (n in 1:N_binary) {
    tj_binary[n] = (tt_binary[n] - 1) * J + jj_binary[n];
    ti_binary[n] = (tt_binary[n] - 1) * I_binary + ii_binary[n];
  }
  for (n in 1:N_trichot) {
    tj_trichot[n] = (tt_trichot[n] - 1) * J + jj_trichot[n];
    ti_trichot[n] = (tt_trichot[n] - 1) * I_trichot + ii_trichot[n];
  }
  for (n in 1:N_ordinal) {
    tj_ordinal[n] = (tt_ordinal[n] - 1) * J + jj_ordinal[n];
    ti_ordinal[n] = (tt_ordinal[n] - 1) * I_ordinal + ii_ordinal[n];
  }
  for (n in 1:N_metric) {
    tj_metric[n] = (tt_metric[n] - 1) * J + jj_metric[n];
    ti_metric[n] = (tt_metric[n] - 1) * I_metric + ii_metric[n];
  }
}
parameters {
  array[T, J, D] real z_eta;                     /* latent factors (deviate) */
  array[T_alpha, I_binary] real z_alpha_binary;	/* intercepts (deviate) */
  matrix[I_binary, D] z_lambda_binary;           /* binary loadings */
  array[T_alpha_k, I_trichot] real z_alpha_trichot; /* intercepts (deviate) */
  array[I_trichot] ordered[2] kappa_trichot;	 /* trichot. thresholds */
  matrix[I_trichot, D] z_lambda_trichot;         /* trichot. loadings */
  array[T_alpha_k, I_ordinal] real z_alpha_ordinal; /* intercepts (deviate) */
  array[I_ordinal] ordered[K_ordinal - 1] kappa_ordinal; /* ordinal thresholds */
  matrix[I_ordinal, D] z_lambda_ordinal;         /* ordinal loadings */
  array[est_sigma_alpha] real<lower=0> sigma_alpha_evol; /* evolution SD of alpha */
  array[T_alpha, I_metric] real z_alpha_metric;    /* intercepts (deviate) */
  matrix[I_metric, D] z_lambda_metric;           /* metric loadings */
  vector<lower=0>[I_metric] sigma_metric;        /* metric residual sd */
  array[est_Omega] vector<lower=0>[D] sigma_eta_evol; /* evolution SD of eta */
  array[est_Omega] cholesky_factor_corr[D] Lcorr_eta; /* evolution corr. (chol.) */
  /* Regularized horseshoe: half-Cauchy local and global scales, written as
     tan() of uniforms for better sampling, plus an inverse-gamma slab. */
  vector<lower=0, upper=pi() / 2>[est_hs ? D : 0] hs_local_raw;
  array[est_hs] real<lower=0, upper=pi() / 2> hs_global_raw;
  array[est_hs] real<lower=0> hs_slab_aux;
  /* Bayesian lasso on each factor's loading scale */
  vector[est_lasso ? D : 0] lasso_omega;
  vector<lower=0>[est_lasso ? D : 0] lasso_tau;
  array[est_lasso] real<lower=0> lasso_rate;
}
transformed parameters {
  array[T, J, D] real eta;		   /* latent factors (final) */
  array[T, I_binary] real alpha_binary;	   /* binary intercepts */
  array[T, I_ordinal] real alpha_ordinal;  /* ordinal intercepts */
  array[T, I_trichot] real alpha_trichot;  /* trichot intercepts */
  array[T, I_metric] real alpha_metric;	   /* metric intercepts */
  array[I_binary, D] real lambda_binary;   /* binary loadings */
  array[I_trichot, D] real lambda_trichot; /* trichot loadings */
  array[I_ordinal, D] real lambda_ordinal; /* ordinal loadings */
  array[I_metric, D] real lambda_metric;   /* metric loadings */
  /* Per-factor loading scale. One scale per factor is shared by every item
     type, since a factor is one latent dimension. */
  vector[est_shrink ? D : 0] load_scale;

  if (est_hs) {
    real hs_tau = tau0 * tan(hs_global_raw[1]);
    vector[D] hs_local = tan(hs_local_raw);
    real c2 = square(slab_scale) * hs_slab_aux[1];
    /* Caps each factor's scale near slab_scale (Piironen & Vehtari 2017) */
    load_scale = hs_tau * sqrt(c2 * square(hs_local)
                               ./ (c2 + square(hs_tau) * square(hs_local)));
  }
  if (est_lasso) {
    load_scale = lasso_omega;
  }

  if (est_shrink) {
    lambda_binary = to_array_2d(z_lambda_binary .* nonzero_binary
                                .* rep_matrix(load_scale', I_binary));
    lambda_trichot = to_array_2d(z_lambda_trichot .* nonzero_trichot
                                 .* rep_matrix(load_scale', I_trichot));
    lambda_ordinal = to_array_2d(z_lambda_ordinal .* nonzero_ordinal
                                 .* rep_matrix(load_scale', I_ordinal));
    lambda_metric = to_array_2d(z_lambda_metric .* nonzero_metric
                                .* rep_matrix(load_scale', I_metric));
  } else {
    lambda_binary = to_array_2d(z_lambda_binary .* nonzero_binary);
    lambda_trichot = to_array_2d(z_lambda_trichot .* nonzero_trichot);
    lambda_ordinal = to_array_2d(z_lambda_ordinal .* nonzero_ordinal);
    lambda_metric = to_array_2d(z_lambda_metric .* nonzero_metric);
  }

  // build cholesky factor for evolution: L_eta = diag(sigma) * Lcorr_eta
  matrix[D, D] L_eta;
  if (est_Omega) {
    L_eta = diag_pre_multiply(sigma_eta_evol[1], Lcorr_eta[1]);
  } else {
    /* No random walk in eta; L_eta is reported as I_D but never used below. */
    L_eta = identity_matrix(D);
  }
for (t in 1:T) {
    if (t == 1) {
      if (whiten_eta == 1) {
        /* Period-1 factors are de-meaned and whitened, fixing the location,
           scale, and rotation of the latent space. */
        eta[t, 1:J, 1:D] = to_array_2d(whiten(to_matrix(z_eta[t, 1:J, 1:D])));
      } else {
        eta[t, 1:J, 1:D] = z_eta[t, 1:J, 1:D];
      }
      alpha_metric[t] = z_alpha_metric[1];
      alpha_binary[t] = z_alpha_binary[1];
      alpha_trichot[t, ] = rep_array(0.0, I_trichot);
      alpha_ordinal[t, ] = rep_array(0.0, I_ordinal);
    } else {
      if (smooth_eta == 0) {
	/* Periods are independent: no random walk in eta */
        eta[t, 1:J, 1:D] = z_eta[t, 1:J, 1:D];
      } else {
        /* Random walk, built directly on the period-1 configuration. */
        for (j in 1:J) {
          eta[t, j, 1:D] =
            to_array_1d(to_vector(eta[t - 1, j, 1:D])
                        + L_eta * to_vector(z_eta[t, j, 1:D]));
        }
      }
      if (constant_alpha == 1) {
        alpha_metric[t] = z_alpha_metric[1]; /* copy first period */
        alpha_binary[t] = z_alpha_binary[1];
        for (i in 1:I_trichot) {
          alpha_trichot[t, i] = 0;
        }
        for (i in 1:I_ordinal) {
          alpha_ordinal[t, i] = 0;
        }
      } else {
        for (i in 1:I_binary) {
          alpha_binary[t][i] = alpha_binary[t - 1][i] +
            z_alpha_binary[t][i] * sigma_alpha_evol[1];
        }
        for (i in 1:I_trichot) {
          alpha_trichot[t][i] = alpha_trichot[t - 1][i] +
            z_alpha_trichot[t][i] * sigma_alpha_evol[1];
        }
        for (i in 1:I_ordinal) {
          alpha_ordinal[t][i] = alpha_ordinal[t - 1][i] +
            z_alpha_ordinal[t][i] * sigma_alpha_evol[1];
        }
        for (i in 1:I_metric) {
          alpha_metric[t][i] = alpha_metric[t - 1][i] +
            z_alpha_metric[t][i] * sigma_alpha_evol[1];
        }
      }
    }
  }
}
model {
  // Likelihood //
  {
    /* Matrix copies of the factors, loadings, and intercepts. Arrays are
       kept as the transformed parameters so the output is unchanged; the
       flattening is row-major, matching tj and ti. */
    matrix[T * J, D] eta_mat = to_matrix(to_array_1d(eta), T * J, D, 0);
    matrix[I_binary, D] L_binary =
      to_matrix(to_array_1d(lambda_binary), I_binary, D, 0);
    matrix[I_trichot, D] L_trichot =
      to_matrix(to_array_1d(lambda_trichot), I_trichot, D, 0);
    matrix[I_ordinal, D] L_ordinal =
      to_matrix(to_array_1d(lambda_ordinal), I_ordinal, D, 0);
    matrix[I_metric, D] L_metric =
      to_matrix(to_array_1d(lambda_metric), I_metric, D, 0);
    vector[T * I_binary] a_binary = to_vector(to_array_1d(alpha_binary));
    vector[T * I_trichot] a_trichot = to_vector(to_array_1d(alpha_trichot));
    vector[T * I_ordinal] a_ordinal = to_vector(to_array_1d(alpha_ordinal));
    vector[T * I_metric] a_metric = to_vector(to_array_1d(alpha_metric));

    if (parallelize == 0) {
      // Linear predictors //
      vector[N_binary] nu_binary;
      vector[N_trichot] nu_trichot;
      vector[N_ordinal] nu_ordinal;
      vector[N_metric] nu_metric;
      profile("linear_predictor") {
        nu_binary = linear_predictor(a_binary, L_binary, eta_mat,
                                     ti_binary, ii_binary, tj_binary);
        nu_trichot = linear_predictor(a_trichot, L_trichot, eta_mat,
                                      ti_trichot, ii_trichot, tj_trichot);
        nu_ordinal = linear_predictor(a_ordinal, L_ordinal, eta_mat,
                                      ti_ordinal, ii_ordinal, tj_ordinal);
        nu_metric = linear_predictor(a_metric, L_metric, eta_mat,
                                     ti_metric, ii_metric, tj_metric);
      }
      profile("likelihood") {
        target += bernoulli_logit_lupmf(yy_binary | p2l_vector(nu_binary));
        target += ordered_logistic_lupmf(yy_trichot | p2l_vector(nu_trichot),
                                         kappa_trichot[ii_trichot]);
        target += ordered_logistic_lupmf(yy_ordinal | p2l_vector(nu_ordinal),
                                         kappa_ordinal[ii_ordinal]);
        target += normal_lupdf(yy_metric | nu_metric,
                               sigma_metric[ii_metric]);
      }
    }
    if (parallelize == 1) {
      profile("parallel") {
        int grainsize = 1;
        target += reduce_sum(bprobit_partial_sum_lpmf,
                             to_array_1d(yy_binary),
                             grainsize,
                             a_binary,
                             L_binary,
                             eta_mat,
                             ti_binary,
                             ii_binary,
                             tj_binary);
        target += reduce_sum(oprobit_partial_sum_lpmf,
                             to_array_1d(yy_trichot),
                             grainsize,
                             a_trichot,
                             L_trichot,
                             eta_mat,
                             ti_trichot,
                             ii_trichot,
                             tj_trichot,
                             kappa_trichot);
        target += reduce_sum(oprobit_partial_sum_lpmf,
                             to_array_1d(yy_ordinal),
                             grainsize,
                             a_ordinal,
                             L_ordinal,
                             eta_mat,
                             ti_ordinal,
                             ii_ordinal,
                             tj_ordinal,
                             kappa_ordinal);
        target += reduce_sum(normal_partial_sum_lpdf,
                             to_array_1d(yy_metric),
                             grainsize,
                             a_metric,
                             L_metric,
                             eta_mat,
                             ti_metric,
                             ii_metric,
                             tj_metric,
                             sigma_metric);
      }
    }
  }

  // Priors //
  to_array_1d(z_eta) ~ std_normal();
  to_array_1d(z_alpha_binary) ~ std_normal();
  to_array_1d(z_lambda_binary) ~ std_normal();
  to_array_1d(z_alpha_trichot) ~ std_normal();
  for (i in 1:I_trichot) {
    to_array_1d(kappa_trichot[i]) ~ std_normal();
  }
  to_array_1d(z_lambda_trichot) ~ std_normal();
  to_array_1d(z_alpha_ordinal) ~ std_normal();
  for (i in 1:I_ordinal) {
    to_array_1d(kappa_ordinal[i]) ~ std_normal();
  }
  to_array_1d(z_lambda_ordinal) ~ std_normal();
  to_array_1d(z_alpha_metric) ~ std_normal();
  to_array_1d(z_lambda_metric) ~ std_normal();

  sigma_metric ~ student_t(df_sigma_metric,
                           mu_sigma_metric,
                           sd_sigma_metric);
  if (est_sigma_alpha) {
    sigma_alpha_evol[1] ~ student_t(df_sigma_alpha_evol,
				    mu_sigma_alpha_evol,
				    sd_sigma_alpha_evol);
  }
  if (est_Omega) {
    sigma_eta_evol[1] ~ student_t(df_sigma_eta_evol,
                                  mu_sigma_eta_evol,
                                  sd_sigma_eta_evol);
    Lcorr_eta[1] ~ lkj_corr_cholesky(2);
  }
  if (est_hs) {
    hs_local_raw ~ uniform(0, pi() / 2);   /* tan() of these is half-Cauchy */
    hs_global_raw[1] ~ uniform(0, pi() / 2);
    hs_slab_aux[1] ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  }
  if (est_lasso) {
    /* As in multiIRT. The gamma(0.001, 0.001) hyperprior is very weak. */
    lasso_omega ~ normal(0, sqrt(lasso_tau));
    lasso_tau ~ exponential(lasso_rate[1] / 2);
    lasso_rate[1] ~ gamma(0.001, 0.001);
  }
}
generated quantities {
  /* In this parameterization the random walk is built directly on the
     whitened period-1 configuration, so the innovation covariance in the
     eta space *is* L_eta L_eta'. No WW rotation is applied. */
  matrix[D, D] Omega = multiply_lower_tri_self_transpose(L_eta);

  /* Per-observation log likelihood, for loo/waic. Positions run in the
     order the observations were supplied: all binary, then trichotomous,
     then ordinal, then metric, each block ordered by period, item, unit
     (see shape_mixfac()). Normalized _lpmf/_lpdf forms are used, unlike
     the _lupmf/_lupdf forms in the model block, so that values are
     comparable across models. */
  vector[gen_log_lik ? N_binary + N_trichot + N_ordinal + N_metric : 0] log_lik;

  if (gen_log_lik) {
    int pos = 1;
    for (n in 1 : N_binary) {
      real nu = alpha_binary[tt_binary[n], ii_binary[n]]
        + dot_product(lambda_binary[ii_binary[n], 1:D],
                      eta[tt_binary[n], jj_binary[n], 1:D]);
      log_lik[pos] = bernoulli_logit_lpmf(yy_binary[n] | p2l(nu));
      pos += 1;
    }
    for (n in 1 : N_trichot) {
      real nu = alpha_trichot[tt_trichot[n], ii_trichot[n]]
        + dot_product(lambda_trichot[ii_trichot[n], 1:D],
                      eta[tt_trichot[n], jj_trichot[n], 1:D]);
      log_lik[pos] = ordered_logistic_lpmf(yy_trichot[n] | p2l(nu),
                                           kappa_trichot[ii_trichot[n]]);
      pos += 1;
    }
    for (n in 1 : N_ordinal) {
      real nu = alpha_ordinal[tt_ordinal[n], ii_ordinal[n]]
        + dot_product(lambda_ordinal[ii_ordinal[n], 1:D],
                      eta[tt_ordinal[n], jj_ordinal[n], 1:D]);
      log_lik[pos] = ordered_logistic_lpmf(yy_ordinal[n] | p2l(nu),
                                           kappa_ordinal[ii_ordinal[n]]);
      pos += 1;
    }
    for (n in 1 : N_metric) {
      real nu = alpha_metric[tt_metric[n], ii_metric[n]]
        + dot_product(lambda_metric[ii_metric[n], 1:D],
                      eta[tt_metric[n], jj_metric[n], 1:D]);
      log_lik[pos] = normal_lpdf(yy_metric[n] | nu, sigma_metric[ii_metric[n]]);
      pos += 1;
    }
  }
}

