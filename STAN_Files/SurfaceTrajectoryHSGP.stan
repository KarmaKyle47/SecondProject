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

    real inv_norm = inv_sqrt(curPos[1] * curPos[1] + curPos[2] * curPos[2] + 1e-6);

    VF[1,1] = curPos[2]*inv_norm;
    VF[2,1] = -1*curPos[1]*inv_norm;
    VF[1,2] = curPos[1]*inv_norm;
    VF[2,2] = curPos[2]*inv_norm;
    //VF[1,3] = sqrt(2);
    //VF[2,3] = sqrt(2);

    return VF;

  }

  vector TrajWeightedBaseVectorFields(data real t, data vector phySpacePos, data vector compSpacePos,
                                      int M, array[,,] real coef_betas, array[,,] real logit_betas, real coef_prior_mean, real logit_prior_mean, vector omega, int N_models){


    vector[N_models] CoefValues;

    for(i in 1:N_models){

      CoefValues[i] = evaluate_HSGP(compSpacePos, M, coef_betas[,,i], coef_prior_mean, omega);

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
    logit_prior_mean = 100;
  } else{
    logit_prior_mean = log(2.0/(M-2));
  }

  real coef_prior_mean = 1.0;

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

}

parameters {
  vector<lower=1e-6>[N_models] coef_ks;      // SD of the coefficient surfaces
  vector<lower=1e-6>[N_models] logit_ks; // SD of the logit surfaces
  vector<lower=0.05>[N_models] coef_ls;   // Length Scales for the coefficient surfaces
  vector<lower=0.05>[N_models] logit_ls;  // Length Scales for the logit surfaces
  array[N_models] matrix[M+1, M+1] coef_zs; // zs for the coefficient surfaces
  array[N_models] matrix[M+1, M+1] logit_zs; // zs for the logit surfaces
  real<lower=1e-6> sigma_vel; // Velocity Sigma
}

transformed parameters {

  array[N_models] matrix[M+1, M+1] coef_betas;
  array[N_models] matrix[M+1, M+1] logit_betas;

  profile("Spectral Scaling"){
    for(i in 1:N_models){
      // Standard HSGP Scaling
      vector[M+1] sqrt_spd_coef = coef_ks[i] * sqrt(sqrt(2*pi()) * coef_ls[i] * exp(-0.5 * square(coef_ls[i] * omega)));
      vector[M+1] sqrt_spd_logit = logit_ks[i] * sqrt(sqrt(2*pi()) * logit_ls[i] * exp(-0.5 * square(logit_ls[i] * omega)));

      coef_betas[i] = diag_post_multiply(diag_pre_multiply(sqrt_spd_coef, coef_zs[i]), sqrt_spd_coef);
      logit_betas[i] = diag_post_multiply(diag_pre_multiply(sqrt_spd_logit, logit_zs[i]), sqrt_spd_logit);
    }
  }





}

model {
  // Priors

  profile("Priors"){
    coef_ks ~ inv_gamma(coef_k_alpha, coef_k_beta); // most of prob between 0.01 and 0.5
    logit_ks ~ inv_gamma(logit_k_alpha, logit_k_beta); //approx jeffrey's
    coef_ls ~ inv_gamma(coef_l_alpha, coef_l_beta);  //approx jeffrey's
    logit_ls ~ inv_gamma(logit_l_alpha, logit_l_beta);  //approx jeffrey's

    sigma_vel ~ inv_gamma(sigma_vel_alpha, sigma_vel_beta); //approx jeffrey's

    for(k in 1:N_models){
       // Pure Standard Normal.
       // No shifting. No scaling. No dynamic dependencies.
       to_vector(coef_zs[k]) ~ std_normal();
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
      vector[N_data] c_dev     = evaluate_surface_vectorized(Phi_x_data, Phi_y_data, coef_betas[k]);

      // B. Apply Explicit Shifts (The "Residual" Logic)
      // Logit = Base + Deviation
      vector[N_data] w = inv_logit(logit_prior_mean + logit_dev);

      // Coeff = Base + Deviation
      // If c_dev is 0 (prior mean), Coeff is 1.0.
      vector[N_data] c = coef_prior_mean + c_dev;

      // C. Accumulate
      vector[N_data] effective_scale = w .* c;
      mu_vx += effective_scale .* col(Base_Vx, k);
      mu_vy += effective_scale .* col(Base_Vy, k);
    }
  }

  profile("Likelihood"){
    Data[, 4] ~ normal(mu_vx, sigma_vel);
    Data[, 5] ~ normal(mu_vy, sigma_vel);
  }



}

generated quantities {
    // Something here
}
