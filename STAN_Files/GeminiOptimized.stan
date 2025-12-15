functions {

  // --- 1. Core HSGP Evaluators ---

  real evaluate_HSGP_Surface(vector curPos, int M, matrix beta, vector omega,
                     real Lx, real Ly, vector border) {
    vector[M+1] phi_x;
    vector[M+1] phi_y;

    phi_x[1] = 1.0;
    phi_y[1] = 1.0;

    // Vectorized Cosine
    phi_x[2:(M+1)] = sqrt(2) * cos(omega[2:(M+1)] * (curPos[1] - border[1]) / Lx);
    phi_y[2:(M+1)] = sqrt(2) * cos(omega[2:(M+1)] * (curPos[2] - border[2]) / Ly);

    return (phi_x' * beta * phi_y);
  }

  // 1D Path Evaluator (Time -> Position)
  real evaluate_HSGP_Path(real t, int M, vector beta, real prior_mean, vector omega,
                     real Lt, data vector boundary) {
    vector[M+1] phi;
    phi[1] = 1.0;
    phi[2:(M+1)] = sqrt(2) * cos(omega[2:(M+1)] * (t - boundary[1]) / Lt);

    return dot_product(beta, phi) + prior_mean;
  }

  vector evaluate_HSGP_Path_Vectorized(vector t, int M, vector beta, real prior_mean,
                                       vector omega, real Lt, vector boundary) {
    int N = num_elements(t);

    // 1. Normalize Time [0, 1]
    vector[N] t_norm = (t - boundary[1]) / Lt;

    // 2. Build Basis Matrix for n=1..M
    // We skip n=0 because it's just the intercept
    row_vector[M] w_row = to_row_vector(omega[2:(M+1)]); // Row vector [1 x M]

    // Outer Product: (N x 1) * (1 x M) -> (N x M)
    matrix[N, M] PHI = sqrt(2) * cos(t_norm * w_row);

    // 3. Project
    // beta[1] is the coefficient for phi_0(t)=1
    // beta[2:M+1] are coefficients for the cosine waves
    vector[N] waves = PHI * beta[2:(M+1)];

    return waves + beta[1] + prior_mean;
  }

  vector predict_path(matrix PHI, vector beta, real prior_mean) {
    // PHI is pre-computed [N x M]
    // beta[2:M+1] corresponds to the waves

    return (PHI * beta[2:rows(beta)]) + beta[1] + prior_mean;
  }

  // --- 2. Physics & Linking ---

  matrix baseVectorFields(real t, vector curPos){
    matrix[2,2] VF;
    real inv_norm = inv_sqrt(curPos[1] * curPos[1] + curPos[2] * curPos[2]);
    // Model 1: Rotation
    VF[1,1] = curPos[2]*inv_norm;
    VF[2,1] = -1*curPos[1]*inv_norm;
    // Model 2: Expansion
    VF[1,2] = curPos[1]*inv_norm;
    VF[2,2] = curPos[2]*inv_norm;
    return VF;
  }

  // Calculates Velocity at a point given the Surface Parameters
  vector TrajWeightedBaseVectorFields(real t, vector curPos,
                                      int M, array[] matrix log_betas,
                                      vector omega, int N_models,
                                      real Lx, real Ly, vector border){
    vector[N_models] CoefValues;
    for(i in 1:N_models){
      CoefValues[i] = exp(evaluate_HSGP_Surface(curPos, M, log_betas[i], omega, Lx, Ly, border));
    }
    // Assuming Weights are 1.0 for now as per your request

    matrix[2,N_models] ModelVels = baseVectorFields(t, curPos);
    return ModelVels * CoefValues;
  }

  // --- 3. The FDA Projector (Optimized) ---

  matrix getMeanCoefsEst(array[] matrix log_betas, matrix path_betas,
                         real Lx, real Ly, vector border, real Lt,
                         int M_Surface, int M_Drifter, vector drifter_prior_means,
                         vector omega_drifter, vector omega_surface, int N_models,
                         vector t_grid, real delta_t, matrix PHI_EstMean, matrix sin_b, int t_res){

    // 2. Evaluate Path Positions (Reconstruct X(t) and Y(t))
    // We need these to query the vector field
    vector[t_res+1] x_seq = predict_path(PHI_EstMean, path_betas[,1], drifter_prior_means[1]);
    vector[t_res+1] y_seq = predict_path(PHI_EstMean, path_betas[,2], drifter_prior_means[2]);

    // 3. Evaluate Vector Field along the Path
    matrix[2, t_res+1] mu_v; // [2 x Time]

    for(i in 1:(t_res+1)){
       // Query the surface physics at the reconstructed position
       mu_v[, i] = TrajWeightedBaseVectorFields(t_grid[i], [x_seq[i], y_seq[i]]', M_Surface, log_betas, omega_surface, N_models, Lx, Ly, border);
    }

    // 4. Project onto Sine Basis (The "Derivative" Basis)
    // We only project onto n=1..M (The mean n=0 derivative is 0)


    // Projection: Integral(Velocity * Basis)
    // (2 x Time) * (Time x M) -> (2 x M)
    matrix[2, M_Drifter] proj = mu_v * sin_b;

    // Scale by dt/Lt (Riemann Sum Normalization)
    return proj' * (delta_t / Lt); // Returns [M x 2]
  }
}

data {
  int<lower=1> N_data;
  int<lower=1> N_drifters;
  int<lower=1> N_models;
  array[N_drifters] int<lower=2> obs_per_drifter;
  matrix[N_data,3] Data;         // [t, x, y]

  int<lower=1> M_Surface;
  int<lower=1> M_Drifter;

  vector[N_models] log_ks;
  vector[N_models] log_ls;
  vector[4] border;

  matrix[N_drifters,2] drifter_ks;
  matrix[N_drifters,2] drifter_ls;

  matrix[N_drifters, 2] drifter_prior_means;
  matrix[N_drifters, 2] drifter_boundaries;
  int<lower=1> t_res;
}

transformed data {
  real Lx = border[3] - border[1];
  real Ly = border[4] - border[2];

  vector[N_drifters] Lts;
  for(i in 1:N_drifters){
    Lts[i] = drifter_boundaries[i,2] - drifter_boundaries[i,1];
  }

  vector[N_drifters] delta_ts = Lts/t_res;


  vector[M_Surface+1] omega_surface;
  for(m in 1:(M_Surface+1)) omega_surface[m] = (m-1) * pi();

  vector[M_Drifter+1] omega_drifter;
  for(m in 1:(M_Drifter+1)) omega_drifter[m] = (m-1) * pi();

  matrix[N_drifters, t_res+1] drifter_t_grids;
  for(i in 1:N_drifters){
    drifter_t_grids[i,] = linspaced_vector(t_res+1, drifter_boundaries[i,1], drifter_boundaries[i,2])';
  }

  matrix[N_data, M_Drifter] PHI_Data;
  array[N_drifters] matrix[t_res+1, M_Drifter] PHI_MeanEst;
  array[N_drifters] matrix[t_res+1, M_Drifter] sin_bs;
  int cur_index = 0;
  for(i in 1:N_drifters){
    int cur_N_obs = obs_per_drifter[i];
    matrix[cur_N_obs,3] curData = Data[(cur_index + 1) : (cur_index + cur_N_obs),];
    vector[cur_N_obs] t_norm_data = (curData[,1] - drifter_boundaries[i,1]) / Lts[i];
    row_vector[M_Drifter] w_row = to_row_vector(omega_drifter[2:(M_Drifter+1)]);
    PHI_Data[(cur_index + 1) : (cur_index + cur_N_obs),] = sqrt(2) * cos(t_norm_data * w_row);

    vector[t_res+1] t_norm_est = (drifter_t_grids[i,]' - drifter_boundaries[i,1]) / Lts[i];
    PHI_MeanEst[i] = sqrt(2) * cos(t_norm_est * w_row);
    sin_bs[i] = sqrt(2) * sin(t_norm_est * w_row);

    cur_index += cur_N_obs;

  }

  real sigma_vel_alpha = 4; real sigma_vel_beta = 1;
  real sigma_pos_alpha = 10; real sigma_pos_beta = 1;
}

parameters {
  // Surface Params
  array[N_models] matrix[M_Surface+1, M_Surface+1] log_zs;

  // Array of Matrices for Drifter Zs
  array[N_drifters] matrix[M_Drifter+1, 2] drifter_zs;

  real<lower=0> sigma_vel;
  real<lower=0> sigma_pos;
}

transformed parameters {
  array[N_models] matrix[M_Surface+1, M_Surface+1] log_betas;
  array[N_drifters] matrix[M_Drifter+1, 2] drifter_betas;

  // 1. Surface Spectral Scaling
  for(i in 1:N_models){
      vector[M_Surface+1] sqrt_spd = sqrt(sqrt(2*pi()) * log_ls[i] * exp(-0.5 * square(log_ls[i] * omega_surface)));
      log_betas[i] = log_ks[i] * diag_post_multiply(diag_pre_multiply(sqrt_spd, log_zs[i]), sqrt_spd);
  }

  // 2. Drifter Path Spectral Scaling
  for(i in 1:N_drifters){
      // Independent X and Y scaling
      vector[M_Drifter+1] spd_x = sqrt(sqrt(2*pi()) * drifter_ls[i,1] * exp(-0.5 * square(drifter_ls[i,1] * omega_drifter)));
      vector[M_Drifter+1] spd_y = sqrt(sqrt(2*pi()) * drifter_ls[i,2] * exp(-0.5 * square(drifter_ls[i,2] * omega_drifter)));

      // Access columns of the matrix inside the array
      drifter_betas[i][:,1] = drifter_ks[i,1] * spd_x .* drifter_zs[i][:,1];
      drifter_betas[i][:,2] = drifter_ks[i,2] * spd_y .* drifter_zs[i][:,2];
  }

  array[N_drifters] matrix[M_Drifter, 2] drifter_betas_der;
  for(i in 1:N_drifters){

    // Access columns of the matrix inside the array
    drifter_betas_der[i][,1] = drifter_betas[i][2:(M_Drifter+1), 1] .* omega_drifter[2:(M_Drifter+1)] * (-1.0/Lts[i]);
    drifter_betas_der[i][,2] = drifter_betas[i][2:(M_Drifter+1), 2] .* omega_drifter[2:(M_Drifter+1)] * (-1.0/Lts[i]);
  }


  real sigma2_pos = sigma_pos*sigma_pos;
}

model {
  // --- Priors ---

  // Use to_vector to flatten matrices for efficient sampling
  for(k in 1:N_models) to_vector(log_zs[k]) ~ std_normal();

  // Corrected Drifter Priors

  sigma_vel ~ inv_gamma(sigma_vel_alpha, sigma_vel_beta);
  sigma_pos ~ inv_gamma(sigma_pos_alpha, sigma_pos_beta);

  for(i in 1:N_drifters) to_vector(drifter_zs[i]) ~ std_normal();

  // --- Collocation Likelihood (Physics Link) ---

  for(i in 1:N_drifters){

    // A. Get Target Coefficients from Vector Field (Physics)
    // Returns [M x 2] matrix for indices n=1..M
    matrix[M_Drifter, 2] c_hat = getMeanCoefsEst(log_betas, drifter_betas[i], Lx, Ly, border, Lts[i],
                                                  M_Surface, M_Drifter, drifter_prior_means[i]', omega_drifter, omega_surface, N_models,
                                                  drifter_t_grids[i,]', delta_ts[i], PHI_MeanEst[i], sin_bs[i], t_res);

    // C. The Physics Mismatch Penalty
    // We treat the Physics Prediction (c_hat) as the "Mean" and the Path Derivative as the "Observation"
    // Variance is scaled by 1/Lt because it's an integral norm
    to_vector(drifter_betas_der[i]) ~ normal(to_vector(c_hat), sigma_vel / sqrt(Lts[i]));
    // Note: sigma_vel is SD. If using variance in normal(), check parameterization. Stan uses SD.
  }

  // --- Observation Likelihood (Data Link) ---

  int start_index = 1;
  for(i in 1:N_drifters){
    int cur_N_obs = obs_per_drifter[i];
    int end_index = start_index + cur_N_obs - 1;

    matrix[cur_N_obs, 3] curData = Data[start_index:end_index,];
  // 1. Extract Time Vector for this drifter

    // 2. Evaluate X and Y paths in one shot
    vector[cur_N_obs] x_pos_pred = predict_path(PHI_Data[start_index:end_index,], drifter_betas[i][,1], drifter_prior_means[i,1]);
    vector[cur_N_obs] y_pos_pred = predict_path(PHI_Data[start_index:end_index,], drifter_betas[i][,2], drifter_prior_means[i,2]);

    // 3. Likelihood
    curData[, 2] ~ normal(x_pos_pred, sigma_pos);
    curData[, 3] ~ normal(y_pos_pred, sigma_pos);
    start_index += cur_N_obs;
  }
}
