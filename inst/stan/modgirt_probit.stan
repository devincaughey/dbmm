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
}
parameters {
  array[Q] ordered[K - 1] z_alpha; // thresholds (difficulties)
  array[Q, D] real beta_free; // unconstrained discriminations
  array[Q, D] real<lower=0> beta_pos; // sign-constrained discriminations
  array[T, G, D] real z_bar_theta;
  vector<lower=0>[D] sd_theta; // within-group SD of theta
  corr_matrix[D] corr_theta;   // within-group correlation of theta across dimensions
  vector<lower=0>[D] sd_bar_theta_evol; // evolution SD of bar_theta
  corr_matrix[D] corr_bar_theta_evol;   // cross-dimension correlation of transition model */
}
transformed parameters {
  array[T, Q] vector[K - 1] alpha; // thresholds (difficulty)
  array[T, Q] real alpha_drift;     // question-specific time trends
  matrix[Q, D] beta;
  array[T] matrix[G, D] bar_theta; // group ideal point means
  cov_matrix[D] Sigma_theta; // within-group variance-covariance
  cov_matrix[D] Omega; // transition variance-covariance
  Sigma_theta = quad_form_diag(corr_theta, sd_theta);
  Omega = quad_form_diag(corr_bar_theta_evol, sd_bar_theta_evol);
  matrix[D, D] chol_Omega = cholesky_decompose(Omega);
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
  sd_theta ~ cauchy(0, 1);
  corr_theta ~ lkj_corr(2);
  sd_bar_theta_evol ~ cauchy(0, .1);
  corr_bar_theta_evol ~ lkj_corr(2);
  /* Likelihood */
  if (K > 1) {
    /* ordinal outcomes */
    for (t in 1 : T) {
      for (q in 1 : Q) {
        real denom; // denominator of linear predictor
        vector[K - 1] cuts; // ordered probit cutpoints
	      real sbs;
	      sbs = quad_form(Sigma_theta[1 : D, 1 : D], to_vector(beta[q][1 : D]));
        denom = sqrt(1 + sbs);
        cuts = alpha[t, q][1 : (K - 1)] / denom;
        for (g in 1 : G) {
          real eta; // linear predictor
          eta = to_row_vector(beta[q][1 : D])
                * to_vector(bar_theta[t, g, 1 : D]) / denom;
          for (k in 1 : K) {
            if (SSSS[t, g, q, k] > 0) {
              /* Add SSSS[t, g, q, k] log normalized ordinal probit densities */
              target += SSSS[t, g, q, k] * ordered_probit_lupmf(k | eta, cuts);
            }
          }
        }
      }
    }
  }
}
