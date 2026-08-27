functions {
  /* De-mean and 'whiten' (cov = I) XX */
  matrix whiten(matrix XX) {
    matrix[rows(XX), cols(XX)] DM;
    matrix[cols(XX), cols(XX)] SS;
    matrix[cols(XX), cols(XX)] PP;
    matrix[cols(XX), cols(XX)] WW;
    for (d in 1 : cols(XX)) {
      DM[ : , d] = XX[ : , d] - mean(XX[ : , d]); /* de-mean each column */
    }
    SS = crossprod(DM) ./ (rows(XX) - 1.0); /* covariance of XX */
    PP = inverse_spd(SS); /* precision of XX */
    WW = cholesky_decompose(PP); /* Cholesky decomposition of precision */
    return DM * WW; /* de-meaned and whitened XX */
  }
  /* Number of group-item-periods with at least one response. Declared as a
     function so that transformed data can size its index arrays in a single
     declaration block. */
  int count_observed_cells(array[,,,] real SS) {
    int n = 0;
    for (t in 1 : dims(SS)[1]) {
      for (g in 1 : dims(SS)[2]) {
        for (q in 1 : dims(SS)[3]) {
          real tot = 0;
          for (k in 1 : dims(SS)[4]) {
            tot += SS[t, g, q, k];
          }
          if (tot > 0) {
            n += 1;
          }
        }
      }
    }
    return n;
  }
}
data {
  int<lower=1> T; // number of periods
  int<lower=1> G; // number of groups
  int<lower=1> Q; // number of items
  int<lower=1> K; // max number of answer options
  int<lower=1> D; // number of latent dimensions
  array[T, G, Q, K] real<lower=0> SSSS; // number of responses (possibly non-integer)
  array[Q, D] int beta_nonzero; // loading point restrictions
  array[Q, D] int beta_sign; // loading sign restrictions
  real<lower=0> shape_lkj_theta;            // 1 = uniform, >1 = diagonal
  real<lower=0> shape_lkj_bar_theta_evol;
  real<lower=0> df_sd_theta;
  real<lower=0> scale_sd_theta;
  real<lower=0> df_sd_bar_theta_evol;
  real<lower=0> scale_sd_bar_theta_evol;
  int<lower=0, upper=1> gen_log_lik; // compute per-cell log_lik?
}
transformed data {
  /* The transition covariance Omega is identified only if T > 1. When T == 1
     there is no transition, so Omega is fixed to the identity rather than
     estimated (see transformed parameters). */
  int<lower=0, upper=1> est_Omega = T > 1;
  /* Group-item-periods with at least one response. Cells with no responses
     contribute nothing to the likelihood, so giving them a log_lik entry of
     exactly zero would create degenerate observations for PSIS. */
  int<lower=0> N_obs = count_observed_cells(SSSS);
  array[N_obs] int tt_obs;
  array[N_obs] int gg_obs;
  array[N_obs] int qq_obs;
  {
    int pos = 1;
    for (t in 1 : T) {
      for (g in 1 : G) {
        for (q in 1 : Q) {
          real tot = 0;
          for (k in 1 : K) {
            tot += SSSS[t, g, q, k];
          }
          if (tot > 0) {
            tt_obs[pos] = t;
            gg_obs[pos] = g;
            qq_obs[pos] = q;
            pos += 1;
          }
        }
      }
    }
  }
}
parameters {
  array[Q] ordered[K - 1] z_alpha; // thresholds (difficulties)
  array[Q, D] real beta_free; // unconstrained discriminations
  array[Q, D] real<lower=0> beta_pos; // sign-constrained discriminations
  array[T, G, D] real z_bar_theta;
  vector<lower=0>[D] sd_theta; // within-group SD of theta
  cholesky_factor_corr[D] L_corr_theta;
  /* Zero-length (hence not estimated) when T == 1 */
  array[est_Omega] vector<lower=0>[D] sd_bar_theta_evol;
  array[est_Omega] cholesky_factor_corr[D] L_corr_bar_theta_evol;
}
transformed parameters {
  array[T, Q] vector[K - 1] alpha; // thresholds (difficulty)
  array[T, Q] real alpha_drift;     // question-specific time trends
  matrix[Q, D] beta;
  array[T] matrix[G, D] bar_theta; // group ideal point means
  matrix[D, D] chol_Sigma_theta = diag_pre_multiply(sd_theta, L_corr_theta);
  matrix[D, D] Sigma_theta =
      multiply_lower_tri_self_transpose(chol_Sigma_theta);
  matrix[D, D] chol_Omega;
  matrix[D, D] Omega;
  if (est_Omega) {
    chol_Omega = diag_pre_multiply(sd_bar_theta_evol[1],
                                   L_corr_bar_theta_evol[1]);
    Omega = multiply_lower_tri_self_transpose(chol_Omega);
  } else {
    chol_Omega = identity_matrix(D);
    Omega = identity_matrix(D);
  }
  for (q in 1 : Q) {
    for (d in 1 : D) {
      if (beta_sign[q, d] == 0) {
        beta[q, d] = beta_nonzero[q, d] * beta_free[q, d];
      } else if (beta_sign[q, d] > 0) {
        beta[q, d] = beta_nonzero[q, d] * beta_pos[q, d];
      } else if (beta_sign[q, d] < 0) {
        beta[q, d] = -1.0 * beta_nonzero[q, d] * beta_pos[q, d];
      }
    }
  }
  for (t in 1 : T) {
    if (t == 1) {
      /* Make period 1 ideal points orthogonal and mean zero */
      bar_theta[t][1 : G, 1 : D] = 
        whiten(to_matrix(z_bar_theta[t, 1 : G, 1 : D]));
    }
    if (t > 1) {
      for (g in 1 : G) {
        vector[D] bt_vec_tm1 = to_vector(bar_theta[t-1][g, 1:D]);
        vector[D] zbt_t = to_vector(z_bar_theta[t, g, 1:D]);
        vector[D] bt_vec_t = bt_vec_tm1 + chol_Omega * zbt_t;
        bar_theta[t][g, 1:D] = to_row_vector(bt_vec_t);
      }
    }
    for (q in 1 : Q) {
      alpha_drift[t, q] = 0; // could estimate but for now set to 0.
      alpha[t, q][1 : (K - 1)] = z_alpha[q][1 : (K - 1)] + alpha_drift[t, q];
    }
  }
}
model {
  /* Priors */
  to_array_1d(z_bar_theta[1 : T, 1 : G, 1 : D]) ~ std_normal();
  to_array_1d(beta_free[1 : Q, 1 : D]) ~ std_normal();
  to_array_1d(beta_pos[1 : Q, 1 : D]) ~ std_normal();
  for (q in 1 : Q) {
    z_alpha[q][1 : (K - 1)] ~ std_normal();
  }
  sd_theta ~ student_t(df_sd_theta, 0, scale_sd_theta);
  L_corr_theta ~ lkj_corr_cholesky(shape_lkj_theta);
  if (est_Omega) {
    sd_bar_theta_evol[1] ~ student_t(df_sd_bar_theta_evol, 0,
                                     scale_sd_bar_theta_evol);
    L_corr_bar_theta_evol[1] ~ lkj_corr_cholesky(shape_lkj_bar_theta_evol);
  }
  /* Likelihood */
  if (K > 1) {
    /* ordinal outcomes */
    for (q in 1 : Q) {
      real sbs = dot_self(chol_Sigma_theta' * to_vector(beta[q][1 : D]));
      real denom = sqrt(1 + sbs);
      for (t in 1 : T) {
        vector[K - 1] cuts = alpha[t, q][1 : (K - 1)] / denom;
        for (g in 1 : G) {
          real eta = to_row_vector(beta[q][1 : D])
                     * to_vector(bar_theta[t, g, 1 : D]) / denom;
          for (k in 1 : K) {
            if (SSSS[t, g, q, k] > 0) {
              target += SSSS[t, g, q, k] * ordered_probit_lupmf(k | eta, cuts);
            }
          }
        }
      }
    }
  }
}
generated quantities {
  /* Per-cell log likelihood, for loo/waic. One element per group-item-period
     with at least one response, in the order period, group, item (item varying
     fastest); see log_lik_index_modgirt(). Each element is the weighted sum
     over response categories, so holding one out removes every respondent who
     answered that item in that group and period. Normalized _lpmf is used,
     unlike the _lupmf in the model block, so that values are comparable
     across models. The multinomial coefficient is omitted: it does not depend
     on the parameters, and is not well defined when SSSS is weighted. */
  vector[gen_log_lik ? N_obs : 0] log_lik;

  if (gen_log_lik) {
    vector[Q] denom_q;
    for (q in 1 : Q) {
      denom_q[q] =
        sqrt(1 + dot_self(chol_Sigma_theta' * to_vector(beta[q][1 : D])));
    }
    for (n in 1 : N_obs) {
      int t = tt_obs[n];
      int g = gg_obs[n];
      int q = qq_obs[n];
      vector[K - 1] cuts = alpha[t, q][1 : (K - 1)] / denom_q[q];
      real eta = to_row_vector(beta[q][1 : D])
                 * to_vector(bar_theta[t, g, 1 : D]) / denom_q[q];
      real ll = 0;
      for (k in 1 : K) {
        if (SSSS[t, g, q, k] > 0) {
          ll += SSSS[t, g, q, k] * ordered_probit_lpmf(k | eta, cuts);
        }
      }
      log_lik[n] = ll;
    }
  }
}
