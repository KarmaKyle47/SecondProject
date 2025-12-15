functions {
  // Now this function is just a lightweight matrix multiplication wrapper
  vector predict_path(matrix PHI, vector beta, real prior_mean) {
    // PHI is pre-computed [N x M]
    // beta[2:M+1] corresponds to the waves

    return (PHI * beta[2:rows(beta)]) + beta[1] + prior_mean;
  }
}

data {
  int<lower=1> N_data;
  matrix[N_data,3] Data;         // [t, x, y]
  int<lower=1> M_Drifter;
  vector[2] drifter_prior_means;
  vector[2] drifter_boundaries;
}

transformed data {
  real Lt = drifter_boundaries[2] - drifter_boundaries[1];

  // 1. Pre-compute Frequencies
  vector[M_Drifter+1] omega_drifter;
  for(m in 1:(M_Drifter+1)) omega_drifter[m] = (m-1) * pi();

  // 2. PRE-COMPUTE THE BASIS MATRIX (The "Golden" Optimization)
  matrix[N_data, M_Drifter] PHI;
  {
      vector[N_data] t_norm = (Data[,1] - drifter_boundaries[1]) / Lt;
      row_vector[M_Drifter] w_row = to_row_vector(omega_drifter[2:(M_Drifter+1)]);

      // Calculate Cosines ONCE. Never again.
      PHI = sqrt(2) * cos(t_norm * w_row);
  }

  // Hyperparameters
  // real drifter_k_alpha = 4; real drifter_k_beta = 10;
  // real drifter_l_alpha = 20; real drifter_l_beta = 1;
  // real sigma_pos_alpha = 10; real sigma_pos_beta = 1;
}

parameters {
  matrix[M_Drifter+1, 2] drifter_zs;
  vector<lower=0>[2] drifter_ks;
  vector<lower=0>[2] drifter_ls;
  real<lower=0> sigma_pos;
}

transformed parameters {
  matrix[M_Drifter+1, 2] drifter_betas;

  vector[M_Drifter+1] spd_x = sqrt(sqrt(2*pi()) * drifter_ls[1] * exp(-0.5 * square(drifter_ls[1] * omega_drifter)));
  vector[M_Drifter+1] spd_y = sqrt(sqrt(2*pi()) * drifter_ls[2] * exp(-0.5 * square(drifter_ls[2] * omega_drifter)));

  drifter_betas[,1] = drifter_ks[1] * spd_x .* drifter_zs[,1];
  drifter_betas[,2] = drifter_ks[2] * spd_y .* drifter_zs[,2];
}

model {
  // --- Priors ---
  to_vector(drifter_ks) ~ std_normal();

  // CRITICAL: Use LogNormal to avoid hitting the 0.1 lower bound wall
  to_vector(drifter_ls) ~ lognormal(log(0.3), 0.5);

  sigma_pos ~ std_normal();
  to_vector(drifter_zs) ~ std_normal();

  // --- Observation Likelihood ---
  // Now extremely fast matrix-vector products

  vector[N_data] x_pos_pred = predict_path(PHI, drifter_betas[,1], drifter_prior_means[1]);
  vector[N_data] y_pos_pred = predict_path(PHI, drifter_betas[,2], drifter_prior_means[2]);

  Data[, 2] ~ normal(x_pos_pred, sigma_pos);
  Data[, 3] ~ normal(y_pos_pred, sigma_pos);
}
