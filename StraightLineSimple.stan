data {
  int<lower=1> D;
  int<lower=1> K;
  int<lower=1> N_drifter;

  vector[N_drifter] t_grid_data;
  matrix[N_drifter, D] x_pos;
  matrix[N_drifter, D] y_pos;
}

parameters {

  vector[D] a_x;
  vector[D] b_x;

  vector[D] a_y;
  vector[D] b_y;

  real<lower=0> sigma_pos;
  real<lower=0> sigma_vel;
  real logcoef1;
  real logcoef2;
}

model {
  logcoef1 ~ normal(0, 5);
  logcoef2 ~ normal(0, 5);

  // OPTIMIZED: std_normal() skips unnecessary math operations
  sigma_pos ~ std_normal();
  sigma_vel ~ std_normal();

  // // OPTIMIZED: Clean N(0,1) prior for the raw weights
  b_x ~ normal(exp(logcoef1), sigma_vel);
  b_y ~ normal(exp(logcoef2), sigma_vel);

  // b_x ~ normal(0, 100);
  // b_y ~ normal(0, 100);


  matrix[N_drifter, D] p_x_mat;
  matrix[N_drifter, D] p_y_mat;

  for(i in 1:D){

    p_x_mat[,i] = a_x[i] + b_x[i] * t_grid_data;
    p_y_mat[,i] = a_y[i] + b_y[i] * t_grid_data;

  }

  to_vector(x_pos) ~ normal(to_vector(p_x_mat), sigma_pos);
  to_vector(y_pos) ~ normal(to_vector(p_y_mat), sigma_pos);
}
