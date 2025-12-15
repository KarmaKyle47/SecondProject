// File: SurfaceTrajectories.stan

functions {

  real evaluate_HSGP_Surface(vector curPos, int M, array[,] real beta, real prior_mean, vector omega,
                     real Lx, real Ly, vector border) {

    // 2. Build Basis Functions (Neumann / Cosine)
    vector[M+1] phi_x;
    vector[M+1] phi_y;

    // n=0 case
    phi_x[1] = 1.0;
    phi_y[1] = 1.0;

    // n>0 cases
    // Stan supports element-wise cos() on vectors
    phi_x[2:(M+1)] = sqrt(2) * cos(omega[2:(M+1)] * (curPos[1] - border[1]) / Lx);
    phi_y[2:(M+1)] = sqrt(2) * cos(omega[2:(M+1)] * (curPos[2] - border[2]) / Ly);

    // 3. Tensor Contraction
    // (vector transpose) * matrix * vector -> scalar
    return (phi_x' * to_matrix(beta) * phi_y) + prior_mean;
  }

  vector evaluate_surface_vectorized(matrix Phi_x, matrix Phi_y, matrix beta) {
      // PURE BASIS EXPANSION. No shifts.
      matrix[rows(Phi_x), cols(Phi_x)] H = Phi_x * beta;
      return rows_dot_product(H, Phi_y);
  }

  real evaluate_HSGP_Path(real t, int M, vector beta, real prior_mean, vector omega,
                     real Lt, vector boundary) {

    // 2. Build Basis Functions (Neumann / Cosine)
    vector[M+1] phi;

    // n=0 case
    phi[1] = 1.0;

    // n>0 cases
    // Stan supports element-wise cos() on vectors
    phi[2:(M+1)] = sqrt(2) * cos(omega[2:(M+1)] * (t - boundary[1]) / Lt);

    // 3. Tensor Contraction
    return sum(to_vector(beta) * phi) + prior_mean;
  }

  real evaluate_HSGP_Der_Path(real t, int M, array[,] real beta, vector omega,
                     real Lt, vector boundary) {

    // 2. Build Basis Functions (Neumann / Cosine)
    vector[M] phi = cos(omega[2:(M+1)] * (t - boundary[1]) / Lt);;

    // n>0 cases
    // Stan supports element-wise cos() on vectors
    phi[2:(M+1)] = sqrt(2) * sin(omega[2:(M+1)] * (t - boundary[1]) / Lt);

    // 3. Tensor Contraction
    // (vector transpose) * matrix * vector -> scalar
    return sum(beta * phi);
  }

  matrix baseVectorFields(real t, vector curPos){

    matrix[2,2] VF;

    real inv_norm = inv_sqrt(curPos[1] * curPos[1] + curPos[2] * curPos[2]);

    VF[1,1] = curPos[2]*inv_norm;
    VF[2,1] = -1*curPos[1]*inv_norm;
    VF[1,2] = curPos[1]*inv_norm;
    VF[2,2] = curPos[2]*inv_norm;
    // VF[1,3] = sqrt(2);
    // VF[2,3] = sqrt(2);
    // VF[1,4] = -1*curPos[2]*inv_norm;
    // VF[2,4] = curPos[1]*inv_norm;
    // VF[1,5] = -1*curPos[1]*inv_norm;
    // VF[2,5] = -1*curPos[2]*inv_norm;

    return VF;

  }

  vector TrajWeightedBaseVectorFields(real t, vector curPos,
                                      int M, array[,,] real log_betas,
                                      vector omega, int N_models,
                                      real Lx, real Ly, vector border){


    vector[N_models] CoefValues;

    for(i in 1:N_models){

      CoefValues[i] = exp(evaluate_HSGP_Surface(curPos, M, log_betas[,,i], log_prior_mean, omega, Lx, Ly, border));

    }

    vector[N_models] TrajValues = CoefValues;

    matrix[2,N_models] ModelVels = baseVectorFields(t, curPos);

    vector[2] TrajWeightedVel = ModelVels * TrajValues;

    return TrajWeightedVel;

  }

  matrix getMeanCoefsEst(array[,,] real log_betas, array[,] real path_betas, real Lx, real Ly, vector border, real Lt, vector boundary, int t_res,
                         int M_Surface, int M_Drifter, vector drifter_prior_means, vector omega_drifter, vector omega_surface, int N_models){

    vector[t_res+1] t_grid = linspaced_vector(boundary[1], boundary[2], t_res + 1);

    real delta_t = t_res / Lt;

    vector[t_res+1] x_seq;
    vector[t_res+1] y_seq;

    matrix[t_res+1, 2] mu_v;

    for(i in 1:(t_res+1)){

      x_seq[i] = evaluate_HSGP_Path(t_grid[i], M_Drifter, path_betas[,1], drifter_prior_means[1], omega_drifter, Lt, boundary);
      y_seq[i] = evaluate_HSGP_Path(t_grid[i], M_Drifter, path_betas[,2], drifter_prior_means[2], omega_drifter, Lt, boundary);

      mu_v[i,] = TrajWeightedBaseVectorFields(t_grid[i], [x_seq[i], y_seq[i]], M_Surface, log_betas, omega_surface, N_models, Lx, Ly, border)'

    }

    vector[t_res+1] norm_t = (t_grid - boundary[1])/Lt;
    row_vector[M_Drifter] omega_row = linspaced_row_vector(M, 1, M) * pi();

    matrix[t_res+1,M_Drifter] sin_b = sqrt(2) * sin(norm_t * omega_row);

    matrix[M_Drifter,2] c = (mu_v * sin_b)' * delta_t) / Lt;

    return c;

  }


}

data {
  int<lower=1> N_data;      // Number of observations
  int<lower=1> N_drifters; //Number of drifters
  array[N_drifters] int<lower=2> obs_per_drifter; //Number of observations per drifter
  matrix[N_data,3] Data;         // Particle Velocities and time
  int<lower=1> M_Surface; // Number of eigenfunctions to use for the surface
  int<lower=1> M_Drifter; //Number of eigenfunctions to use for the paths
  vector[4] border; // (lower x, lower y, upper x, upper y)
  matrix[N_drifters, 2] drifter_prior_means; // empirical means of each drifter's positions (x,y)
  matrix[N_drifters, 2] drifter_boundaries; // starting and ending time for each drifter
  int<lower=1> t_res;
}

transformed data {

  real Lx = border[3] - border[1];
  real Ly = border[4] - border[2];
  vector Lts = drifter_boundaries[,2] - drifter_boundaries[,1];

  vector[M_Surface+1] omega_surface = 0:M_Surface * pi();
  vector[M_Drifter+1] omega_drifter = 0:M_Drifter * pi();

  int N_models = 2;    // Number of models

  real log_prior_mean = 0.0;

  real log_k_alpha = 9;
  real log_k_beta = 2;

  real log_l_alpha = 4;
  real log_l_beta = 1;

  real path_k_alpha = 4;
  real path_k_beta = 10;

  real path_l_alpha = 20;
  real path_l_beta = 1;

  real sigma_vel_alpha = 4;
  real sigma_vel_beta = 1;

  real sigma_pos_alpha = 10;
  real sigma_pos_beta = 1;

  matrix[N_data, M+1] Phi_x_data;
  matrix[N_data, M+1] Phi_y_data;

  // Fill Column 1 (Mean / n=0)
  Phi_x_data[, 1] = rep_vector(1.0, N_data);
  Phi_y_data[, 1] = rep_vector(1.0, N_data);

  // Fill Columns 2:M+1 (Ripples)
  for (m in 2:(M+1)) {
     // Vectorized cos() on the whole data column
     Phi_x_data[, m] = sqrt(2) * cos(omega[m] * (Data[, 2] - border[1]) / Lx);
     Phi_y_data[, m] = sqrt(2) * cos(omega[m] * (Data[, 3] - border[2]) / Ly);
  }

  matrix[N_data, N_models] Base_Vx;
  matrix[N_data, N_models] Base_Vy;

  for(i in 1:N_data) {
     // Run your physics math ONCE
     matrix[2, 2] vf = baseVectorFields(Data[i,1], Data[i, 2:3]');

     // Store X-velocities
     Base_Vx[i, 1] = vf[1, 1]; // Model 1 X
     Base_Vx[i, 2] = vf[1, 2]; // Model 2 X
     // Base_Vx[i, 3] = vf[1, 3]; // Model 1 X
     // Base_Vx[i, 4] = vf[1, 4]; // Model 2 X
     // Base_Vx[i, 5] = vf[1, 5]; // Model 1 X

     // Store Y-velocities
     Base_Vy[i, 1] = vf[2, 1]; // Model 1 Y
     Base_Vy[i, 2] = vf[2, 2]; // Model 2 Y
     // Base_Vy[i, 3] = vf[2, 3]; // Model 1 Y
     // Base_Vy[i, 4] = vf[2, 4]; // Model 2 Y
     // Base_Vy[i, 5] = vf[2, 5]; // Model 1 Y
  }

}

parameters {
  vector<lower=0>[N_models] log_ks;   // SD of the log surfaces
  vector<lower=0>[N_models] log_ls;   // Length Scales for the log surfaces
  array[N_models] matrix[M_Surface+1, M_Surface+1] log_zs; // zs for the log surfaces
  matrix<lower=0>[N_drifters,2] drifter_ks;
  matrix<lower=0>[N_drifters,2] drifter_ls;
  array[N_drifter] matrix[M_Drifter+1,2] drifter_zs;
  real<lower=0> sigma_vel; // Velocity Error Sigma
  real<lower=0> sigma_pos; // Positional Error Sigma
}

transformed parameters {

  array[N_models] matrix[M_Surface+1, M_Surface+1] log_betas;
  array[N_drifter] matrix[M_Drifter+1,2] drifter_betas;

  profile("Spectral Scaling"){
    for(i in 1:N_models){
      // Standard HSGP Scaling
      vector[M_Surface+1] sqrt_spd_log = sqrt(sqrt(2*pi()) * log_ls[i] * exp(-0.5 * square(log_ls[i] * omega_surface)));

      log_betas[i] = log_ks[i] * diag_post_multiply(diag_pre_multiply(sqrt_spd_log, log_zs[i]), sqrt_spd_log);
    }

    for(i in 1:N_drifters){

      vector[M_Drifters+1] sqrt_spd_x = sqrt(sqrt(2*pi()) * drifter_ls[i,1] * exp(-0.5 * square(drifter_ls[i,1] * omega_drifter)));
      vector[M_Drifters+1] sqrt_spd_y = sqrt(sqrt(2*pi()) * drifter_ls[i,2] * exp(-0.5 * square(drifter_ls[i,2] * omega_drifter)));

      col(drifter_betas[i],1) = drifter_ks[i,1] * sqrt_spd_x .* col(drifter_zs[i],1);
      col(drifter_betas[i],2) = drifter_ks[i,2] * sqrt_spd_y .* col(drifter_zs[i],2);
    }
  }

  real sigma2_vel = sigma_vel*sigma_vel;
  real sigma2_pos = sigma_pos*sigma_pos;


}

model {
  // Priors

  profile("Priors"){
    log_ks ~ inv_gamma(log_k_alpha, log_k_beta); // most of prob between 0.01 and 0.5
    log_ls ~ inv_gamma(log_l_alpha, log_l_beta);  //approx jeffrey's

    drifter_ks ~ inv_gamma(drifter_k_alpha, drifter_k_beta);
    drifter_ls ~ inv_gamma(drifter_l_alpha, drifter_l_beta);

    sigma_vel ~ inv_gamma(sigma_vel_alpha, sigma_vel_beta); //approx jeffrey's
    sigma_pos ~ inv_gamma(sigma_pos_alpha, sigma_pos_beta);

    for(k in 1:N_models){
       // Pure Standard Normal.
       // No shifting. No scaling. No dynamic dependencies.
       to_vector(log_zs[k]) ~ std_normal();
    }

    for(d in 1:N_drifters){
      to_vector(drifter_zs[i]) ~ std_normal();
    }

  }

  // Bridge between Path and Trajectories

  // Getting the means for the derivative coefficients

  for(i in 1:N_drifters){

    matrix[M_Drifter,2] c_hat = getMeanCoefsEst(log_betas, drifter_betas, Lx, Ly, border, Lt, drifter_boundaries[i,]', t_res,
                                                M_Surface, M_Drifter, drifter_prior_means[i,]', omega_drifter, omega_surface, N_models);

    matrix[M_Drifter,2] drifter_der_betas;

    drifter_der_betas[,1] = col(drifter_betas[i],1) .* omega_drifter[1:M_Drifter+1] * pi() * (1/Lt) * -1;
    drifter_der_betas[,2] = col(drifter_betas[i],2) .* omega_drifter[1:M_Drifter+1] * pi() * (1/Lt) * -1;

    to_vector(drifter_der_betas) ~ normal(to_vector(c_hat), sigma2_pos/Lt);

  }


/*  // --- 2. Likelihood ---
  vector[N_data] mu_vx = rep_vector(0, N_data);
  vector[N_data] mu_vy = rep_vector(0, N_data);

  profile("Surface Evaluation"){
    for (k in 1:N_models) {
      // A. Calculate Deviation Surfaces (Centered at 0)
      vector[N_data] log_values = evaluate_surface_vectorized(Phi_x_data, Phi_y_data, log_betas[k]);

      // Log = Base + Deviation
      // If log_dev is 0 (prior mean), Log is 1.0.
      vector[N_data] c = exp(log_prior_mean + log_values);

      // C. Accumulate
      vector[N_data] effective_scale = c;
      mu_vx += effective_scale .* col(Base_Vx, k);
      mu_vy += effective_scale .* col(Base_Vy, k);
    }
  }*/

  profile("Likelihood"){

/*    Data[, 4] ~ normal(mu_vx, sigma_vel);
    Data[, 5] ~ normal(mu_vy, sigma_vel);*/


    int start_index = 1;

    for(i in 1:N_drifters){

      int cur_N_obs = obs_per_drifter[i];

      int end_index = start_index + cur_N_obs - 1;

      matrix[cur_N_obs,3] curData = Data[start_index:end_index,];

      vector[cur_N_obs] x_pos_pred;
      vector[cur_N_obs] y_pos_pred;

      for(i in 1:cur_N_obs){

        x_pos_pred[i] = evaluate_HSGP_Path(curData[i,1], M_drifter, col(drifter_betas[i],1), drifter_prior_means[i,1], omega_drifter, Lt, boundary)
        y_pos_pred[i] = evaluate_HSGP_Path(curData[i,1], M_drifter, col(drifter_betas[i],2), drifter_prior_means[i,2], omega_drifter, Lt, boundary)

      }

      x_resid = curData[,2] - x_pos_pred;
      y_resid = curData[,3] - y_pos_pred;

      x_resid ~ normal(0, sigma2_pos);
      y_resid ~ normal(0, sigma2_pos);


      // vector[cur_N_obs-1] cur_t_present = Data[start_index:(end_index-1),1];
      // vector[cur_N_obs-1] cur_x_pos_present = Data[start_index:(end_index-1),2];
      // vector[cur_N_obs-1] cur_y_pos_present = Data[start_index:(end_index-1),3];
      // vector[cur_N_obs-1] cur_x_vel_present = mu_vx[start_index:(end_index-1)];
      // vector[cur_N_obs-1] cur_y_vel_present = mu_vy[start_index:(end_index-1)];
      //
      // vector[cur_N_obs-1] cur_t_future = Data[(start_index+1):end_index,1];
      // vector[cur_N_obs-1] cur_x_pos_future = Data[(start_index+1):end_index,2];
      // vector[cur_N_obs-1] cur_y_pos_future = Data[(start_index+1):end_index,3];
      //
      // vector[cur_N_obs-1] cur_dt = cur_t_future - cur_t_present;
      //
      // cur_x_pos_future ~ normal(cur_x_pos_present + cur_x_vel_present.*cur_dt, sqrt(2*sigma2_pos + sigma2_vel.*cur_dt));
      // cur_y_pos_future ~ normal(cur_y_pos_present + cur_y_pos_present.*cur_dt, sqrt(2*sigma2_pos + sigma2_vel.*cur_dt));

      start_index += cur_N_obs;

    }



  }



}

generated quantities {
    // Something here
}
