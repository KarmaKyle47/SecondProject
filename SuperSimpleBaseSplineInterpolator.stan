data {
  int<lower=1> D;
  int<lower=1> K;
  int<lower=1> N_drifter;
  matrix[N_drifter, K] Phi;

  vector[N_drifter] t_grid_data;
  matrix[N_drifter, D] x_pos;
  matrix[N_drifter, D] y_pos;
}

parameters {
  // NEW: Non-centered raw weights
  matrix[K, D] z_x;
  matrix[K, D] z_y;

  // NEW: Changed to row_vector to enable matrix broadcasting
  row_vector<lower=-10, upper=10>[D] start_pos_x;
  row_vector<lower=-10, upper=10>[D] start_pos_y;

  real<lower=0> sigma_pos;
  real<lower=0> sigma_vel;
  real logcoef1;
  real logcoef2;
}

transformed parameters {
  // NEW: Reconstruct the true weights from the raw Z-scores
  matrix[K, D] w_x = z_x * sigma_vel;
  matrix[K, D] w_y = z_y * sigma_vel;
}

model {
  logcoef1 ~ normal(0, 0.35);
  logcoef2 ~ normal(0, 0.35);

  // OPTIMIZED: std_normal() skips unnecessary math operations
  sigma_pos ~ std_normal();
  sigma_vel ~ std_normal();

  // OPTIMIZED: Clean N(0,1) prior for the raw weights
  to_vector(z_x) ~ std_normal();
  to_vector(z_y) ~ std_normal();

  vector[N_drifter] p_x_base = exp(logcoef1) * t_grid_data;
  vector[N_drifter] p_y_base = exp(logcoef2) * t_grid_data;

  // OPTIMIZED: Replaced the for-loop with direct matrix broadcasting
  // rep_matrix(vector, D) copies columns. rep_matrix(row_vector, N) copies rows.
  matrix[N_drifter, D] p_x_base_mat = rep_matrix(p_x_base, D) + rep_matrix(start_pos_x, N_drifter);
  matrix[N_drifter, D] p_y_base_mat = rep_matrix(p_y_base, D) + rep_matrix(start_pos_y, N_drifter);

  matrix[N_drifter, D] p_x_mat = p_x_base_mat + Phi * w_x;
  matrix[N_drifter, D] p_y_mat = p_y_base_mat + Phi * w_y;

  to_vector(x_pos) ~ normal(to_vector(p_x_mat), sigma_pos);
  to_vector(y_pos) ~ normal(to_vector(p_y_mat), sigma_pos);
}
