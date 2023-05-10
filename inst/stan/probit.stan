/* This version produced workable results with 2 dimensions. */

functions {
  /* De-mean and 'whiten' (cov = I) XX */
  real p2l_real (real x) {	// coverts scalar from probit to logit scale
    real y;
    y = 0.07056 * pow(x, 3) + 1.5976 * x;
    return y;
  }
  matrix whiten(matrix XX) {
    matrix[rows(XX), cols(XX)] DM;
    matrix[cols(XX), cols(XX)] SS;
    matrix[cols(XX), cols(XX)] PP;
    matrix[cols(XX), cols(XX)] WW;
    for (d in 1:cols(XX)) {
      DM[:, d] = XX[:, d] - mean(XX[:, d]); // de-mean each column
    }
    SS = crossprod(DM) ./ (rows(XX) - 1.0); // covariance of XX
    PP = inverse_spd(SS);                   // precision of XX
    WW = cholesky_decompose(PP);            // Cholesky decomposition
    return DM * WW;                         // de-meaned and whitened XX
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
    int T = dims(eta)[1];
    int J = dims(eta)[2];
    int D = dims(eta)[3];
    int N_slice = end - start + 1;
    array[N_slice] int tt_slice = tt_b[start:end];
    array[N_slice] int jj_slice = jj_b[start:end];
    array[N_slice] int ii_slice = ii_b[start:end];
    array[N_slice] real nu_slice;
    for (n in 1:N_slice) {
      real a_n = alpha_b[tt_slice[n], ii_slice[n]];
      row_vector[D] l_n = to_row_vector(lambda_b[ii_slice[n], 1:D]);
      vector[D] e_n = to_vector(eta[tt_slice[n], jj_slice[n], 1:D]);
      nu_slice[n] = a_n + l_n * e_n;
    }
    return bernoulli_logit_lupmf(yy_b_slice | p2l_real(nu_slice));
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
    int T = dims(eta)[1];
    int J = dims(eta)[2];
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
    return ordered_logistic_lupmf(yy_o_slice | p2l_real(nu_slice), kappa_slice);
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
    int T = dims(eta)[1];
    int J = dims(eta)[2];
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
  int num_matches(array[] int x, int a) {
    int n = 0;
    for (i in 1:size(x))
      if (x[i] == a)
        n += 1;
    return n;
  }
  array[] int which_equal(array[] int x, int a) {
    array[num_matches(x, a)] int match_positions;
    int pos = 1;
    for (i in 1:size(x)) {
      if (x[i] == a) {
        match_positions[pos] = i;
        pos += 1;
      }
    }
    return match_positions;
  }
}
data {
  int<lower=0,upper=1> parallelize;    /* parallelize within chains? */
  int<lower=0,upper=1> constant_alpha; /* keep alphas constant? */
  int<lower=0,upper=1> separate_eta;   /* estimate eta separately by period */
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
  array[T, 2] int<lower=0> tob_b;		  /* time ranges */
  matrix[I_binary, D] nonzero_binary;		  /* nonzero loadings */
  // Ordinal data //
  int<lower=0> N_ordinal;                   /* number of observations */
  int<lower=0> I_ordinal;                   /* number of items */
  int<lower=1> K_ordinal;                   /* max response categories */
  array[N_ordinal] int<lower=1> yy_ordinal; /* outcomes */
  array[N_ordinal] int<lower=1> ii_ordinal; /* item indicator */
  array[N_ordinal] int<lower=1> jj_ordinal; /* unit indicator */
  array[N_ordinal] int<lower=1> tt_ordinal; /* period indicator */
  array[T, 2] int<lower=0> tob_o;	    /* time ranges */
  matrix[I_ordinal, D] nonzero_ordinal;     /* nonzero loadings */
  // Metric data //
  int<lower=0> N_metric;                    /* number of observations */
  int<lower=0> I_metric;                    /* number of items */
  vector[N_metric] yy_metric;               /* outcomes */
  array[N_metric] int<lower=1> ii_metric;   /* item indicator */
  array[N_metric] int<lower=1> jj_metric;   /* unit indicator */
  array[N_metric] int<lower=1> tt_metric;   /* period indicator */
  array[T, 2] int<lower=0> tob_m;	    /* time ranges */
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
}
parameters {
  array[T, J, D] real z_eta;                     /* latent factors (raw) */
  array[T, I_binary] real z_alpha_binary;        /* intercepts (raw) */
  matrix[I_binary, D] z_lambda_binary;           /* binary loadings */
  array[T, I_ordinal] real z_alpha_ordinal;      /* intercepts (raw) */
  array[I_ordinal] ordered[K_ordinal - 1] kappa; /* ordinal thresholds */
  matrix[I_ordinal, D] z_lambda_ordinal;         /* ordinal loadings */
  real<lower=0> sigma_alpha_evol;                /* evolution SD of alpha */
  array[T, I_metric] real z_alpha_metric;        /* intercepts (raw) */
  matrix[I_metric, D] z_lambda_metric;           /* metric loadings */
  vector<lower=0>[I_metric] sigma_metric;        /* metric residual sd */
  vector<lower=0>[D] sigma_eta_evol;             /* evolution SD of eta */
}
transformed parameters {
  array[T, J, D] real eta;                /* latent factors (whitened) */
  array[T, I_binary] real alpha_binary;   /* binary intercepts */
  array[T, I_ordinal] real alpha_ordinal; /* ordinal intercepts */
  array[T, I_metric] real alpha_metric;   /* metric intercepts */
  array[I_binary, D] real lambda_binary;  /* binary loadings */
  array[I_ordinal, D] real lambda_ordinal; /* ordinal loadings */
  array[I_metric, D] real lambda_metric;   /* metric loadings */
  lambda_binary = to_array_2d(z_lambda_binary .* nonzero_binary);
  lambda_ordinal = to_array_2d(z_lambda_ordinal .* nonzero_ordinal);
  lambda_metric = to_array_2d(z_lambda_metric .* nonzero_metric);
  for (t in 1:T) {
    if (t == 1) {
      eta[t, 1:J, 1:D] = to_array_2d(whiten(to_matrix(z_eta[t, 1:J, 1:D])));
      alpha_metric[t] = z_alpha_metric[t];
      alpha_binary[t] = z_alpha_binary[t];
      alpha_ordinal[t, ] = rep_array(0.0, I_ordinal);
    }
    if (t > 1) {
      if (separate_eta == 1) {
	eta[t, 1:J, 1:D] = z_eta[t, 1:J, 1:D];
      } else {
	for (j in 1:J) {
	  for (d in 1:D) {
	    eta[t, j, d] = eta[t - 1, j, d] + z_eta[t, j, d] * sigma_eta_evol[d];
	  }
	}
      } 
      if (constant_alpha == 1) {
        alpha_metric[t] = z_alpha_metric[1]; /* copy first period */
        alpha_binary[t] = z_alpha_binary[1]; /* copy first period */
        for (i in 1:I_ordinal) {
          alpha_ordinal[t, i] = 0;
        }
      } else {
        for (i in 1:I_binary) {
          alpha_binary[t][i] = alpha_binary[t - 1][i] +
            z_alpha_binary[t][i] * sigma_alpha_evol;
        }
        for (i in 1:I_ordinal) {
          alpha_ordinal[t][i] = alpha_ordinal[t - 1][i] +
            z_alpha_ordinal[t][i] * sigma_alpha_evol;
        }
        for (i in 1:I_metric) {
          alpha_metric[t][i] = alpha_metric[t - 1][i] +
            z_alpha_metric[t][i] * sigma_alpha_evol;
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
    vector[N_ordinal] nu_ordinal;
    vector[N_metric] nu_metric;
    profile("linear_predictor") {
      for (t in 1:T) {
	for (n in 1:N_binary) {
	  nu_binary[n] = alpha_binary[tt_binary[n], ii_binary[n]] +
	    to_row_vector(lambda_binary[ii_binary[n], 1:D]) *
	    to_vector(eta[tt_binary[n], jj_binary[n], 1:D]);
	}
	for (n in 1:N_ordinal) {
	  nu_ordinal[n] = alpha_ordinal[tt_ordinal[n], ii_ordinal[n]] +
	    to_row_vector(lambda_ordinal[ii_ordinal[n], 1:D]) *
	    to_vector(eta[tt_ordinal[n], jj_ordinal[n], 1:D]);
	}
	for (n in 1:N_metric) {
	  nu_metric[n] = alpha_metric[tt_metric[n], ii_metric[n]] +
	    to_row_vector(lambda_metric[ii_metric[n], 1:D]) *
	    to_vector(eta[tt_metric[n], jj_metric[n], 1:D]);
	}
      }
    }
    profile("likelihood") {
      target += bernoulli_logit_lupmf(yy_binary | p2l_real(nu_binary));
      target += ordered_logistic_lupmf(yy_ordinal | p2l_real(nu_ordinal),
				       kappa[ii_ordinal]);
      target += normal_lupdf(yy_metric | nu_metric,
			     sigma_metric[ii_metric]);
    }
  }
  if (parallelize == 1) {
    profile("parallel") {
      int grainsize = 1;
      target += reduce_sum(bprobit_partial_sum_lupmf,
                           to_array_1d(yy_binary),
                           grainsize,
                           alpha_binary,
                           lambda_binary,
                           eta,
                           tt_binary,
                           ii_binary,
                           jj_binary);
      target += reduce_sum(oprobit_partial_sum_lupmf,
                           to_array_1d(yy_ordinal),
                           grainsize,
                           alpha_ordinal,
                           lambda_ordinal,
                           eta,
                           tt_ordinal,
                           ii_ordinal,
                           jj_ordinal,                   
                           kappa);
      target += reduce_sum(normal_partial_sum_lupdf,
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
  to_array_1d(z_alpha_ordinal) ~ std_normal();
  for (i in 1:I_ordinal) {
    to_array_1d(kappa[i]) ~ std_normal();
  }
  to_array_1d(z_lambda_ordinal) ~ std_normal();
  to_array_1d(z_alpha_metric) ~ std_normal();
  to_array_1d(z_lambda_metric) ~ std_normal();
  sigma_metric ~ student_t(df_sigma_metric,
			   mu_sigma_metric,
			   sd_sigma_metric);
  sigma_alpha_evol ~ student_t(df_sigma_alpha_evol,
			       mu_sigma_alpha_evol,
			       sd_sigma_alpha_evol);
  sigma_eta_evol ~ student_t(df_sigma_eta_evol,
			     mu_sigma_eta_evol,
			     sd_sigma_eta_evol);
}
generated quantities {
  
}
