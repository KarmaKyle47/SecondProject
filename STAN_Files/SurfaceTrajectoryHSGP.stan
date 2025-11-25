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

  vector evaluate_surface_vectorized(matrix Phi_x, matrix Phi_y, matrix beta, real prior_mean) {
    // 1. Project Phi_x through the coefficients
    // (N_data x M) * (M x M) -> (N_data x M)
    matrix[rows(Phi_x), cols(Phi_x)] H = Phi_x * beta;

    // 2. Dot product with Phi_y row-by-row
    // Result is vector of length N_data
    return rows_dot_product(H, Phi_y) + prior_mean;
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
                                      int M, array[,,] real coef_betas, array[,,] real logit_betas, vector coef_prior_means, real logit_prior_mean, vector omega, int N_models){


    vector[N_models] CoefValues;

    for(i in 1:N_models){

      CoefValues[i] = evaluate_HSGP(compSpacePos, M, coef_betas[,,i], coef_prior_means[i], omega);

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
  int<lower=0> N_data;      // Number of observations
  matrix[N_data,5] Data;         // Particle Velocities with Positions for now (t, x, y, v_x, v_y)
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
    logit_prior_mean = 0;
  } else{
    logit_prior_mean = log(2.0/(M-2));
  }

  real coef_k_alpha = 9;
  real coef_k_beta = 2;

  real logit_k_alpha = 4;
  real logit_k_beta = 1;

  real coef_l_alpha = 4;
  real coef_l_beta = 1;

  real logit_l_alpha = 4;
  real logit_l_beta = 1;

  real sigma_vel_alpha = 4;
  real sigma_vel_beta = 1;

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

  array[N_data] matrix[2, 2] Base_Vels;

  for(i in 1:N_data) {
     Base_Vels[i] = baseVectorFields(Data[i,1], Data[i, 2:3]');
  }

}

parameters {
  vector<lower=0>[N_models] coef_ks;      // SD of the coefficient surfaces
  vector<lower=0>[N_models] logit_ks; // SD of the logit surfaces
  vector<lower=0>[N_models] coef_ls;   // Length Scales for the coefficient surfaces
  vector<lower=0>[N_models] logit_ls;  // Length Scales for the logit surfaces
  vector[N_models] coef_prior_means;   // Prior mean of the coefficient surfaces
  array[N_models] matrix[M+1, M+1] coef_zs; // zs for the coefficient surfaces
  array[N_models] matrix[M+1, M+1] logit_zs; // zs for the logit surfaces
  real<lower=0> sigma_vel; // Velocity Sigma
}

transformed parameters {

  array[N_models] matrix[M+1, M+1] coef_betas;
  array[N_models] matrix[M+1, M+1] logit_betas;

  // OPTIMIZATION: No diag_matrix()
  // Use diag_pre/post_multiply for speed
  for(i in 1:N_models){
    // 1. Calculate Spectral Density
    vector[M+1] sqrt_spd_coef = coef_ks[i] * sqrt(sqrt(2*pi()) * coef_ls[i] * exp(-0.5 * square(coef_ls[i] * omega)));
    vector[M+1] sqrt_spd_logit = logit_ks[i] * sqrt(sqrt(2*pi()) * logit_ls[i] * exp(-0.5 * square(logit_ls[i] * omega)));

    // 2. Scale Z -> Beta
    coef_betas[i] = diag_post_multiply(diag_pre_multiply(sqrt_spd_coef, coef_zs[i]), sqrt_spd_coef);
    logit_betas[i] = diag_post_multiply(diag_pre_multiply(sqrt_spd_logit, logit_zs[i]), sqrt_spd_logit);
  }

}

model {
  // Priors
  coef_ks ~ inv_gamma(coef_k_alpha, coef_k_beta); // most of prob between 0.01 and 0.5
  logit_ks ~ inv_gamma(logit_k_alpha, logit_k_beta); //approx jeffrey's
  coef_ls ~ inv_gamma(coef_l_alpha, coef_l_beta);  //approx jeffrey's
  logit_ls ~ inv_gamma(logit_l_alpha, logit_l_beta);  //approx jeffrey's

  coef_prior_means ~ normal(1, 0.5); //centered at 0

  sigma_vel ~ inv_gamma(sigma_vel_alpha, sigma_vel_beta); //approx jeffrey's

  for(k in 1:N_models){
     to_vector(coef_zs[k]) ~ std_normal();
     to_vector(logit_zs[k]) ~ std_normal();
  }

  //Likelihood

  vector[N_data] mu_vx = rep_vector(0, N_data);
  vector[N_data] mu_vy = rep_vector(0, N_data);

  // 1. Calculate contributions from each model
  for (k in 1:N_models) {

      // A. Get Weight Vector (for all data points at once)
      vector[N_data] logits = evaluate_surface_vectorized(Phi_x_data, Phi_y_data, logit_betas[k], logit_prior_mean);
      vector[N_data] w = inv_logit(logits);

      // B. Get Coefficient Vector (for all data points at once)
      vector[N_data] c = evaluate_surface_vectorized(Phi_x_data, Phi_y_data, coef_betas[k], coef_prior_means[k]);

      // C. Accumulate Physics
      // We loop over data here only for the final summation, which is cheap
      // (The expensive surface evaluation is already done)
      for (i in 1:N_data) {
         // w[i] * c[i] * Base_Vel[i]
         // Base_Vels[i] is 2x2. Column k is the vector for model k.
         vector[2] v_k = Base_Vels[i][, k];

         real scalar_mult = w[i] * c[i];

         mu_vx[i] += scalar_mult * v_k[1];
         mu_vy[i] += scalar_mult * v_k[2];
      }
  }

  // 2. Evaluate Normal (Vectorized)
  Data[, 4] ~ normal(mu_vx, sigma_vel);
  Data[, 5] ~ normal(mu_vy, sigma_vel);

}

generated quantities {
    // Something here
}
