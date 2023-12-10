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
  array[Q] ordered[K - 1] alpha_raw; // thresholds (difficulties)
  array[Q, D] real beta_free; // unconstrained discriminations
  array[Q, D] real<lower=0> beta_pos; // sign-constrained discriminations
  array[T, G, D] real z_bar_theta;
  vector<lower=0>[D] sd_theta; // within-group SD of theta
  vector<lower=0>[D] sd_bar_theta_evol; // evolution SD of bar_theta
}
transformed parameters {
  array[T, Q] vector[K - 1] alpha; // thresholds (difficulty)
  matrix[Q, D] beta;
  array[T] matrix[G, D] bar_theta; // group ideal point means
  cov_matrix[D] Sigma_theta; // diagonal matrix of within-group variances
  Sigma_theta = diag_matrix(sd_theta .* sd_theta);
  cov_matrix[D] Sigma_bar_theta_evol; // diagonal matrix of within-group variances
  Sigma_bar_theta_evol = diag_matrix(sd_bar_theta_evol .* sd_bar_theta_evol);
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
      bar_theta[t][1 : G, 1 : D] = whiten(to_matrix(z_bar_theta[t, 1 : G, 1 : D]));
    }
    if (t > 1) {
      for (g in 1 : G) {
        for (d in 1 : D) {
          bar_theta[t][g, d] = bar_theta[t - 1][g, d]
                               + z_bar_theta[t, g, d] * sd_bar_theta_evol[d];
        }
      }
    }
    for (q in 1 : Q) {
      /* Use same alpha each period */
      alpha[t, q][1 : (K - 1)] = alpha_raw[q][1 : (K - 1)];
    }
  }
}
model {
  /* Priors */
  to_array_1d(z_bar_theta[1 : T, 1 : G, 1 : D]) ~ std_normal();
  to_array_1d(beta_free[1 : Q, 1 : D]) ~ std_normal();
  to_array_1d(beta_pos[1 : Q, 1 : D]) ~ std_normal();
  for (q in 1 : Q) {
    alpha_raw[q][1 : (K - 1)] ~ std_normal();
  }
  for (d in 1 : D) {
    sd_theta[d] ~ cauchy(0, 1);
    sd_bar_theta_evol[d] ~ cauchy(0, .1);
  }
  /* Likelihood */
  if (K > 1) {
    /* ordinal outcomes */
    for (t in 1 : T) {
      for (q in 1 : Q) {
        real denom; // denominator of linear predictor
        vector[K - 1] cuts; // ordered probit cutpoints
	      real sb;
	      sb = quad_form(Sigma_theta[1 : D, 1 : D], to_vector(beta[q][1 : D]));
        denom = sqrt(1 + sb);
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
