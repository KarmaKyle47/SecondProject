// File: SurfaceTrajectories.stan

functions {

  real evaluate_HSGP(vector curPos, int M, array[,] real beta, real prior_mean, vector omega) {

    // 2. Build Basis Functions (Neumann / Cosine)
    vector[M+1] phi_x;
    vector[M+1] phi_y;

    // n=0 case
    phi_x[1] = 1.0;
    phi_y[1] = 1.0;

    // n>0 cases
    // Stan supports element-wise cos() on vectors
    phi_x[2:(M+1)] = sqrt(2) * cos(omega[2:(M+1)] * curPos[1]);
    phi_y[2:(M+1)] = sqrt(2) * cos(omega[2:(M+1)] * curPos[2]);

    // 3. Tensor Contraction
    // (vector transpose) * matrix * vector -> scalar
    return (phi_x' * to_matrix(beta) * phi_y) + prior_mean;
  }

vector evaluate_surface_vectorized(matrix Phi_x, matrix Phi_y, matrix beta) {
    // PURE BASIS EXPANSION. No shifts.
    matrix[rows(Phi_x), cols(Phi_x)] H = Phi_x * beta;
    return rows_dot_product(H, Phi_y);
}

  vector getCompSpacePos(data matrix GMM_means, array[,,] real GMM_cov, data vector GMM_weights, data vector curPos_Phy,
                         data vector GMM_sd_x, data vector GMM_sd_y_cond, data vector GMM_cond_slope){
    int n_mixtures = rows(GMM_means);
    real x_phy = curPos_Phy[1];
    real y_phy = curPos_Phy[2];

    real x_cdf = 0.0;
    vector[n_mixtures] cond_log_weights_uw;

    for(i in 1:n_mixtures){
      x_cdf += normal_cdf(x_phy | GMM_means[i,1], GMM_sd_x[i])*GMM_weights[i];
      cond_log_weights_uw[i] = normal_lpdf(x_phy | GMM_means[i,1], GMM_sd_x[i]) + log(GMM_weights[i]);
    }

    vector[n_mixtures] cond_log_weights = cond_log_weights_uw - log_sum_exp(cond_log_weights_uw);
    vector[n_mixtures] cond_weights = exp(cond_log_weights);

    real y_given_x_cdf = 0.0;

    for(i in 1:n_mixtures){
      real cur_cond_mean = GMM_means[i,2] + GMM_cond_slope[i]*(x_phy - GMM_means[i,1]);
      real cur_cond_sd = GMM_sd_y_cond[i];

      y_given_x_cdf += normal_cdf(y_phy | cur_cond_mean, cur_cond_sd)*cond_weights[i];
    }

    return [x_cdf, y_given_x_cdf]';

  }

  matrix baseVectorFields(data real t, data vector curPos){

    matrix[2,2] VF;

    real inv_norm = inv_sqrt(curPos[1] * curPos[1] + curPos[2] * curPos[2]);

    VF[1,1] = curPos[2]*inv_norm;
    VF[2,1] = -1*curPos[1]*inv_norm;
    VF[1,2] = curPos[1]*inv_norm;
    VF[2,2] = curPos[2]*inv_norm;
    //VF[1,3] = sqrt(2);
    //VF[2,3] = sqrt(2);

    return VF;

  }

  vector TrajWeightedBaseVectorFields(data real t, data vector phySpacePos, data vector compSpacePos,
                                      int M, array[,,] real log_betas, array[,,] real logit_betas, real log_prior_mean, real logit_prior_mean, vector omega, int N_models){


    vector[N_models] CoefValues;

    for(i in 1:N_models){

      CoefValues[i] = exp(evaluate_HSGP(compSpacePos, M, log_betas[,,i], log_prior_mean, omega));

    }

    vector[N_models] WeightValues;

    for(i in 1:N_models){

      WeightValues[i] = inv_logit(evaluate_HSGP(compSpacePos, M, logit_betas[,,i], logit_prior_mean, omega));

    }

    vector[N_models] TrajValues = CoefValues .* WeightValues;

    matrix[2,N_models] ModelVels = baseVectorFields(t, phySpacePos);

    vector[2] TrajWeightedVel = ModelVels * TrajValues;

    return TrajWeightedVel;

  }


}

data {
  int<lower=1> N_data;      // Number of observations
  int<lower=1> N_drifters; //Number of drifters
  array[N_drifters] int<lower=2> obs_per_drifter; //Number of observations per drifter
  matrix[N_data,3] Data;         // Particle Velocities and time
  int<lower=1> GMM_num; // Number of Gaussian Mixtures in the Transformation
  matrix[GMM_num,2] GMM_means; // Means of all Gaussian Mixtures;
  array[2,2,GMM_num] real GMM_cov; // Covariance matrices of all Mixtures
  simplex[GMM_num] GMM_weights; // Weights for Mixtures
  int<lower=1> M; // Number of eigenfunctions to use
}

transformed data {
  vector[M+1] omega;
  for (m in 1:(M+1)) {
    omega[m] = (m - 1) * pi();
  }

  int N_models = 2;    // Number of models

  real logit_prior_mean;

  if(N_models == 2){
    logit_prior_mean = 10;
  } else{
    logit_prior_mean = log(2.0/(N_models-2));
  }

  real log_prior_mean = 0.0;

  real log_k_alpha = 9;
  real log_k_beta = 2;

  real logit_k_alpha = 4;
  real logit_k_beta = 1;

  real log_l_alpha = 4;
  real log_l_beta = 1;

  real logit_l_alpha = 4;
  real logit_l_beta = 1;

  real sigma_vel_alpha = 4;
  real sigma_vel_beta = 1;

  real sigma_pos_alpha = 10;
  real sigma_pos_beta = 1;

  // --- Pre-calculate GMM values ---
  vector[GMM_num] GMM_sd_x;        // sqrt(cov[1,1])
  vector[GMM_num] GMM_sd_y_cond;   // Sqrt of conditional variance
  vector[GMM_num] GMM_cond_slope;  // cov[2,1] / cov[1,1]

  for(i in 1:GMM_num) {
    real inv_cov_11 = 1.0 / GMM_cov[1,1,i];
    real cond_slope = GMM_cov[2,1,i] * inv_cov_11;
    GMM_sd_x[i] = sqrt(GMM_cov[1,1,i]);
    GMM_cond_slope[i] = cond_slope;
    GMM_sd_y_cond[i] = sqrt(GMM_cov[2,2,i] - cond_slope * GMM_cov[1,2,i]);
  }

  matrix[N_data, 2] compSpacePos_data;
  for (i in 1:N_data) {
    compSpacePos_data[i, 1:2] = getCompSpacePos(GMM_means, GMM_cov, GMM_weights, Data[i, 2:3]',
                                                GMM_sd_x, GMM_sd_y_cond, GMM_cond_slope)';
  }

  matrix[N_data, M+1] Phi_x_data;
  matrix[N_data, M+1] Phi_y_data;

  // Fill Column 1 (Mean / n=0)
  Phi_x_data[, 1] = rep_vector(1.0, N_data);
  Phi_y_data[, 1] = rep_vector(1.0, N_data);

  // Fill Columns 2:M+1 (Ripples)
  for (m in 2:(M+1)) {
     // Vectorized cos() on the whole data column
     Phi_x_data[, m] = sqrt(2) * cos(omega[m] * compSpacePos_data[, 1]);
     Phi_y_data[, m] = sqrt(2) * cos(omega[m] * compSpacePos_data[, 2]);
  }

  matrix[N_data, N_models] Base_Vx;
  matrix[N_data, N_models] Base_Vy;

  for(i in 1:N_data) {
     // Run your physics math ONCE
     matrix[2, 2] vf = baseVectorFields(Data[i,1], Data[i, 2:3]');

     // Store X-velocities
     Base_Vx[i, 1] = vf[1, 1]; // Model 1 X
     Base_Vx[i, 2] = vf[1, 2]; // Model 2 X

     // Store Y-velocities
     Base_Vy[i, 1] = vf[2, 1]; // Model 1 Y
     Base_Vy[i, 2] = vf[2, 2]; // Model 2 Y
  }

  // vector<lower=0>[N_models] log_ks = [0.175, 0.175]';      // SD of the log surfaces
  // vector<lower=0>[N_models] logit_ks = [2.0, 2.0]'; // SD of the logit surfaces
  // vector<lower=0>[N_models] log_ls = [0.2,0.2]';   // Length Scales for the log surfaces
  // vector<lower=0>[N_models] logit_ls = [0.1,0.1]';  // Length Scales for the logit surfaces

}

parameters {
  vector<lower=0>[N_models] log_ks;      // SD of the log surfaces
  vector<lower=0>[N_models] logit_ks; // SD of the logit surfaces
  vector<lower=0>[N_models] log_ls;   // Length Scales for the log surfaces
  vector<lower=0>[N_models] logit_ls;  // Length Scales for the logit surfaces
  array[N_models] matrix[M+1, M+1] log_zs; // zs for the log surfaces
  array[N_models] matrix[M+1, M+1] logit_zs; // zs for the logit surfaces
  real<lower=0> sigma_vel; // Velocity Error Sigma
  real<lower=0> sigma_pos; // Positional Error Sigma
}

transformed parameters {

  array[N_models] matrix[M+1, M+1] log_betas;
  array[N_models] matrix[M+1, M+1] logit_betas;

  profile("Spectral Scaling"){
    for(i in 1:N_models){
      // Standard HSGP Scaling
      vector[M+1] sqrt_spd_log = sqrt(sqrt(2*pi()) * log_ls[i] * exp(-0.5 * square(log_ls[i] * omega)));
      vector[M+1] sqrt_spd_logit = sqrt(sqrt(2*pi()) * logit_ls[i] * exp(-0.5 * square(logit_ls[i] * omega)));

      log_betas[i] = log_ks[i] * diag_post_multiply(diag_pre_multiply(sqrt_spd_log, log_zs[i]), sqrt_spd_log);
      logit_betas[i] = logit_ks[i] * diag_post_multiply(diag_pre_multiply(sqrt_spd_logit, logit_zs[i]), sqrt_spd_logit);
    }
  }

  real sigma2_vel = sigma_vel*sigma_vel;
  real sigma2_pos = sigma_pos*sigma_pos;

}

model {
  // Priors

  profile("Priors"){
    log_ks ~ inv_gamma(log_k_alpha, log_k_beta); // most of prob between 0.01 and 0.5
    logit_ks ~ inv_gamma(logit_k_alpha, logit_k_beta); //approx jeffrey's
    log_ls ~ inv_gamma(log_l_alpha, log_l_beta);  //approx jeffrey's
    logit_ls ~ inv_gamma(logit_l_alpha, logit_l_beta);  //approx jeffrey's

    sigma_vel ~ inv_gamma(sigma_vel_alpha, sigma_vel_beta); //approx jeffrey's
    sigma_pos ~ inv_gamma(sigma_pos_alpha, sigma_pos_beta);

    for(k in 1:N_models){
       // Pure Standard Normal.
       // No shifting. No scaling. No dynamic dependencies.
       to_vector(log_zs[k]) ~ std_normal();
       to_vector(logit_zs[k]) ~ std_normal();
    }

  }


  // --- 2. Likelihood ---
  vector[N_data] mu_vx = rep_vector(0, N_data);
  vector[N_data] mu_vy = rep_vector(0, N_data);

  profile("Surface Evaluation"){
    for (k in 1:N_models) {
      // A. Calculate Deviation Surfaces (Centered at 0)
      vector[N_data] logit_dev = evaluate_surface_vectorized(Phi_x_data, Phi_y_data, logit_betas[k]);
      vector[N_data] log_dev     = evaluate_surface_vectorized(Phi_x_data, Phi_y_data, log_betas[k]);

      // B. Apply Explicit Shifts (The "Residual" Logic)
      // Logit = Base + Deviation
      vector[N_data] w = inv_logit(logit_prior_mean + logit_dev);

      // Log = Base + Deviation
      // If log_dev is 0 (prior mean), Log is 1.0.
      vector[N_data] c = exp(log_prior_mean + log_dev);

      // C. Accumulate
      vector[N_data] effective_scale = w .* c;
      mu_vx += effective_scale .* col(Base_Vx, k);
      mu_vy += effective_scale .* col(Base_Vy, k);
    }
  }

  profile("Likelihood"){

    int start_index = 1;

    for(i in 1:N_drifters){

      int cur_N_obs = obs_per_drifter[i];

      int end_index = start_index + cur_N_obs - 1;

      vector[cur_N_obs-1] cur_t_present = Data[start_index:(end_index-1),1];
      vector[cur_N_obs-1] cur_x_pos_present = Data[start_index:(end_index-1),2];
      vector[cur_N_obs-1] cur_y_pos_present = Data[start_index:(end_index-1),3];
      vector[cur_N_obs-1] cur_x_vel_present = mu_vx[start_index:(end_index-1)];
      vector[cur_N_obs-1] cur_y_vel_present = mu_vy[start_index:(end_index-1)];

      vector[cur_N_obs-1] cur_t_future = Data[(start_index+1):end_index,1];
      vector[cur_N_obs-1] cur_x_pos_future = Data[(start_index+1):end_index,2];
      vector[cur_N_obs-1] cur_y_pos_future = Data[(start_index+1):end_index,3];

      vector[cur_N_obs-1] cur_dt = cur_t_future - cur_t_present;

      cur_x_pos_future ~ normal(cur_x_pos_present + cur_x_vel_present.*cur_dt, sqrt(2*sigma2_pos + sigma2_vel.*cur_dt));
      cur_y_pos_future ~ normal(cur_y_pos_present + cur_y_pos_present.*cur_dt, sqrt(2*sigma2_pos + sigma2_vel.*cur_dt));

      start_index += cur_N_obs;

    }



  }



}

generated quantities {
    // Something here
}
