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

  real bprobit_partial_sum_lpmf(array[] int yy_b_slice,
                               int start,
                               int end,
                               array[,] real alpha_b,
                               array[,] real lambda_b,
                               array[,,] real eta,
                               array[] int tt_b,
                               array[] int ii_b,
                               array[] int jj_b) {
    int D = dims(eta)[3];
    int N_slice = end - start + 1;
    array[N_slice] int tt_slice = tt_b[start:end];
    array[N_slice] int jj_slice = jj_b[start:end];
    array[N_slice] int ii_slice = ii_b[start:end];
    vector[N_slice] nu_slice;
    for (n in 1:N_slice) {
      real a_n = alpha_b[tt_slice[n], ii_slice[n]];
      row_vector[D] l_n = to_row_vector(lambda_b[ii_slice[n], 1:D]);
      vector[D] e_n = to_vector(eta[tt_slice[n], jj_slice[n], 1:D]);
      nu_slice[n] = a_n + l_n * e_n;
    }
    return bernoulli_logit_lupmf(yy_b_slice | p2l_vector(nu_slice));
  }

  real oprobit_partial_sum_lpmf(array[] int yy_o_slice,
                               int start,
                               int end,
                               array[,] real alpha_o,
                               array[,] real lambda_o,
                               array[,,] real eta,
                               array[] int tt_o,
                               array[] int ii_o,
                               array[] int jj_o,
                               array[] vector kappa) {
    int D = dims(eta)[3];
    int K = size(kappa[1]);
    int N_slice = end - start + 1;
    array[N_slice] int tt_slice = tt_o[start:end];
    array[N_slice] int jj_slice = jj_o[start:end];
    array[N_slice] int ii_slice = ii_o[start:end];
    vector[N_slice] nu_slice;
    array[N_slice] vector[K] kappa_slice = kappa[ii_slice];
    for (n in 1:N_slice) {
      real a_n = alpha_o[tt_slice[n], ii_slice[n]];
      row_vector[D] l_n = to_row_vector(lambda_o[ii_slice[n], 1:D]);
      vector[D] e_n = to_vector(eta[tt_slice[n], jj_slice[n], 1:D]);
      nu_slice[n] = a_n + l_n * e_n;
    }
    return ordered_logistic_lupmf(yy_o_slice | p2l_vector(nu_slice), kappa_slice);
  }

  real normal_partial_sum_lpdf(array[] real yy_m_slice,
                               int start,
                               int end,
                               array[,] real alpha_m,
                               array[,] real lambda_m,
                               array[,,] real eta,
                               array[] int tt_m,
                               array[] int ii_m,
                               array[] int jj_m,
                               array[] real sigma) {
    int D = dims(eta)[3];
    int N_slice = end - start + 1;
    array[N_slice] int tt_slice = tt_m[start:end];
    array[N_slice] int jj_slice = jj_m[start:end];
    array[N_slice] int ii_slice = ii_m[start:end];
    array[N_slice] real nu_slice;
    array[N_slice] real sigma_slice = sigma[ii_slice];
    for (n in 1:N_slice) {
      real a_n = alpha_m[tt_slice[n], ii_slice[n]];
      row_vector[D] l_n = to_row_vector(lambda_m[ii_slice[n], 1:D]);
      vector[D] e_n = to_vector(eta[tt_slice[n], jj_slice[n], 1:D]);
      nu_slice[n] = a_n + l_n * e_n;
    }
    return normal_lupdf(yy_m_slice | nu_slice, sigma_slice);
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
}
transformed data {
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

  lambda_binary = to_array_2d(z_lambda_binary .* nonzero_binary);
  lambda_trichot = to_array_2d(z_lambda_trichot .* nonzero_trichot);
  lambda_ordinal = to_array_2d(z_lambda_ordinal .* nonzero_ordinal);
  lambda_metric = to_array_2d(z_lambda_metric .* nonzero_metric);

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
  if (parallelize == 0) {
    // Linear predictors //
    vector[N_binary] nu_binary;
    vector[N_trichot] nu_trichot;
    vector[N_ordinal] nu_ordinal;
    vector[N_metric] nu_metric;
    profile("linear_predictor") {
      // compute each nu once
      for (n in 1:N_binary) {
        int tt = tt_binary[n];
        int ii = ii_binary[n];
        int jj = jj_binary[n];
        nu_binary[n] = alpha_binary[tt, ii] +
          dot_product(lambda_binary[ii, 1:D],
                      eta[tt, jj, 1:D]);
      }
      for (n in 1:N_trichot) {
        int tt = tt_trichot[n];
        int ii = ii_trichot[n];
        int jj = jj_trichot[n];
        nu_trichot[n] = alpha_trichot[tt, ii] +
          dot_product(lambda_trichot[ii, 1:D],
                      eta[tt, jj, 1:D]);
      }
      for (n in 1:N_ordinal) {
        int tt = tt_ordinal[n];
        int ii = ii_ordinal[n];
        int jj = jj_ordinal[n];
        nu_ordinal[n] = alpha_ordinal[tt, ii] +
          dot_product(lambda_ordinal[ii, 1:D],
                      eta[tt, jj, 1:D]);
      }
      for (n in 1:N_metric) {
        int tt = tt_metric[n];
        int ii = ii_metric[n];
        int jj = jj_metric[n];
        nu_metric[n] = alpha_metric[tt, ii] +
          dot_product(lambda_metric[ii, 1:D],
                      eta[tt, jj, 1:D]);
      }
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
                           alpha_binary,
                           lambda_binary,
                           eta,
                           tt_binary,
                           ii_binary,
                           jj_binary);
      target += reduce_sum(oprobit_partial_sum_lpmf,
                           to_array_1d(yy_trichot),
                           grainsize,
                           alpha_trichot,
                           lambda_trichot,
                           eta,
                           tt_trichot,
                           ii_trichot,
                           jj_trichot,
                           kappa_trichot);
      target += reduce_sum(oprobit_partial_sum_lpmf,
                           to_array_1d(yy_ordinal),
                           grainsize,
                           alpha_ordinal,
                           lambda_ordinal,
                           eta,
                           tt_ordinal,
                           ii_ordinal,
                           jj_ordinal,
                           kappa_ordinal);
      target += reduce_sum(normal_partial_sum_lpdf,
                           to_array_1d(yy_metric),
                           grainsize,
                           alpha_metric,
                           lambda_metric,
                           eta,
                           tt_metric,
                           ii_metric,
                           jj_metric,
                           to_array_1d(sigma_metric));
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
}
generated quantities {
  /* In this parameterization the random walk is built directly on the
     whitened period-1 configuration, so the innovation covariance in the
     eta space *is* L_eta L_eta'. No WW rotation is applied. */
  cov_matrix[D] Omega = multiply_lower_tri_self_transpose(L_eta);

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

