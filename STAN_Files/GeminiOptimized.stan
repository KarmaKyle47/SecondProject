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

  // --- PARALLEL LIKELIHOOD FUNCTION ---
  real partial_sum_likelihood(
      array[] int slice_drifters, // The chunk of drifter indices (e.g., 1..5)
      int start, int end,         // Stan internal indices
      // --- Parameters ---
      array[] matrix drifter_betas,
      array[] matrix drifter_betas_der,
      array[] matrix log_betas,
      real sigma_vel,
      real sigma_pos,
      // --- Data ---
      matrix Data,
      matrix PHI_Data,
      array[] matrix PHI_MeanEst,
      array[] matrix sin_bs,
      matrix drifter_t_grids,
      array[] int obs_per_drifter,
      array[] int start_indices,    // Pre-computed index map
      matrix drifter_prior_means,
      matrix drifter_boundaries,
      vector delta_ts,
      vector Lts,
      // --- Geometry/Constants ---
      real Lx, real Ly, vector border,
      int M_Surface, int M_Drifter, int N_models,
      vector omega_drifter, vector omega_surface, int t_res
  ) {
      real lp = 0;

      // Loop over the drifters in this specific chunk
      for (i in 1:size(slice_drifters)) {
          int d = slice_drifters[i]; // Get the actual Drifter ID

          // --- 1. Physics Likelihood ---
          matrix[M_Drifter, 2] c_hat = getMeanCoefsEst(
              log_betas, drifter_betas[d], Lx, Ly, border, Lts[d],
              M_Surface, M_Drifter, drifter_prior_means[d]', omega_drifter, omega_surface, N_models,
              drifter_t_grids[d,]', delta_ts[d], PHI_MeanEst[d], sin_bs[d], t_res
          );

          // Accumulate Physics Log-Prob
          lp += normal_lpdf(to_vector(drifter_betas_der[d]) | to_vector(c_hat), sigma_vel / sqrt(Lts[d]));

          // --- 2. Observation Likelihood ---
          int cur_N = obs_per_drifter[d];
          int s_idx = start_indices[d]; // O(1) Lookup instead of loop

          // Slice Data and Basis for this drifter
          // Note: using block() or slicing syntax is fast in Stan
          vector[cur_N] x_pred = predict_path(PHI_Data[s_idx : s_idx + cur_N - 1, ], col(drifter_betas[d], 1), drifter_prior_means[d,1]);
          vector[cur_N] y_pred = predict_path(PHI_Data[s_idx : s_idx + cur_N - 1, ], col(drifter_betas[d], 2), drifter_prior_means[d,2]);

          lp += normal_lpdf(Data[s_idx : s_idx + cur_N - 1, 2] | x_pred, sigma_pos);
          lp += normal_lpdf(Data[s_idx : s_idx + cur_N - 1, 3] | y_pred, sigma_pos);
      }
      return lp;
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

  array[N_drifters] int start_indices;
  int current_idx = 1;
  for(i in 1:N_drifters) {
      start_indices[i] = current_idx;
      current_idx += obs_per_drifter[i];
  }

  // We create a sequence 1, 2, ... N_drifters
  array[N_drifters] int seq_drifters;
  for(i in 1:N_drifters) seq_drifters[i] = i;

  // Grainsize: 1 implies "split as finely as possible" (good for heavy per-item cost)
  int grainsize = 1;

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
  array[N_drifters] matrix[M_Drifter, 2] drifter_betas_der;

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

    drifter_betas_der[i][,1] = drifter_betas[i][2:(M_Drifter+1), 1] .* omega_drifter[2:(M_Drifter+1)] * (-1.0/Lts[i]);
    drifter_betas_der[i][,2] = drifter_betas[i][2:(M_Drifter+1), 2] .* omega_drifter[2:(M_Drifter+1)] * (-1.0/Lts[i]);
  }

}

model {
  // --- Priors ---

  // Use to_vector to flatten matrices for efficient sampling
  for(k in 1:N_models) to_vector(log_zs[k]) ~ std_normal();
  for(i in 1:N_drifters) to_vector(drifter_zs[i]) ~ std_normal();

  // Corrected Drifter Priors
  sigma_vel ~ std_normal();
  sigma_pos ~ std_normal();

    // --- PARALLEL LIKELIHOOD ---

  target += reduce_sum(
      partial_sum_likelihood, seq_drifters, grainsize,
      // Pass Parameters
      drifter_betas, drifter_betas_der, log_betas, sigma_vel, sigma_pos,
      // Pass Data
      Data, PHI_Data, PHI_MeanEst, sin_bs, drifter_t_grids,
      obs_per_drifter, start_indices, drifter_prior_means, drifter_boundaries,
      delta_ts, Lts,
      // Pass Constants
      Lx, Ly, border, M_Surface, M_Drifter, N_models,
      omega_drifter, omega_surface, t_res
  );

}
