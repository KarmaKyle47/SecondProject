functions {
  // --- 1. Vectorized HSGP Evaluator ---
  vector evaluate_HSGP_Surface_vec(int N_quad, vector p_x, vector p_y, int M,
                                   matrix beta, vector omega,
                                   real Lx, real Ly, vector border) {
    matrix[N_quad, M+1] phi_x;
    matrix[N_quad, M+1] phi_y;

    phi_x[,1] = rep_vector(1.0, N_quad);
    phi_y[,1] = rep_vector(1.0, N_quad);

    for (m in 2:(M+1)) {
      phi_x[,m] = sqrt(2.0) * cos(omega[m] * (p_x - border[1]) / Lx);
      phi_y[,m] = sqrt(2.0) * cos(omega[m] * (p_y - border[2]) / Ly);
    }
    return rows_dot_product(phi_x * beta, phi_y);
  }

  // --- Vectorized Neural Network Emulator Evaluator ---
  array[] matrix evaluate_VF_vec(vector t, vector x, vector y,
                                 int N_points, int N_models){


    vector[N_points] window = cos(pi()*x/20) .* cos(pi()*y/20);

    array[N_models] matrix[N_points, 2] VFs;

    VFs[1][,1] = y.*window;
    VFs[1][,2] = -1*x.*window;

    VFs[2][,1] = x.*window.*sin(t/8);
    VFs[2][,2] = y.*window.*sin(t/8);

    return VFs;
  }

  // --- 2. Parallel Likelihood Function ---
  // Notice K, c_part, Phi, and Z are completely gone from the arguments!
  real partial_sum(array[] int drifter_slice, int start, int end,
                   vector w_x, vector w_y, real sigma_vel,
                   array[] vector p_x_base, array[] vector p_y_base,
                   array[] vector v_x_base, array[] vector v_y_base,
                   matrix Phi_Z_flat, matrix Phi_d_Z_flat,
                   array[] int M, array[] int M_starts,
                   array[] vector t_quad, vector dt_weight,
                   int M_Surface, array[] matrix log_betas, vector omega_surface,
                   real Lx, real Ly, vector border,
                   int N_models) {

    real log_lik = 0.0;
    int N_quad = rows(Phi_Z_flat);

    for (i in 1:size(drifter_slice)) {
      int d = drifter_slice[i];
      real sigma_scaled = sigma_vel / dt_weight[d];

      int idx_start = M_starts[d];
      int idx_end = idx_start + M[d] - 1;

      vector[M[d]] w_x_d = w_x[idx_start:idx_end];
      vector[M[d]] w_y_d = w_y[idx_start:idx_end];

      // Compute pos/vel exclusively from precomputed matrices
      vector[N_quad] p_x = p_x_base[d] + Phi_Z_flat[, idx_start:idx_end] * w_x_d;
      vector[N_quad] p_y = p_y_base[d] + Phi_Z_flat[, idx_start:idx_end] * w_y_d;
      vector[N_quad] v_x = v_x_base[d] + Phi_d_Z_flat[, idx_start:idx_end] * w_x_d;
      vector[N_quad] v_y = v_y_base[d] + Phi_d_Z_flat[, idx_start:idx_end] * w_y_d;
      vector[N_quad] t_val = t_quad[d];

      array[N_models] matrix[N_quad, 2] VFs = evaluate_VF_vec(
        t_val, p_x, p_y, N_quad, N_models
      );

      vector[N_quad] coef1 = exp(evaluate_HSGP_Surface_vec(N_quad, p_x, p_y, M_Surface, log_betas[1], omega_surface, Lx, Ly, border));
      vector[N_quad] target_v_x = VFs[1][, 1] .* coef1;
      vector[N_quad] target_v_y = VFs[1][, 2] .* coef1;

      vector[N_quad] coef2 = exp(evaluate_HSGP_Surface_vec(N_quad, p_x, p_y, M_Surface, log_betas[2], omega_surface, Lx, Ly, border));
      target_v_x += VFs[2][, 1] .* coef2;
      target_v_y += VFs[2][, 2] .* coef2;

      log_lik += normal_lpdf(v_x | target_v_x, sigma_scaled) +
                 normal_lpdf(v_y | target_v_y, sigma_scaled);
    }

    return log_lik;
  }
}

data {
  int<lower=1> D;
  int<lower=1> N_quad;

  // --- Ragged Dimensions ---
  array[D] int<lower=1> M;   // null space columns per drifter
  int<lower=1> total_M;

  array[D] int<lower=1> K;   // NEW: basis functions per drifter
  int<lower=1> total_K;      // NEW: sum of K array
  int<lower=1> total_Z;      // NEW: sum of (K[d] * M[d])

  // --- Flattened Ragged Arrays ---
  vector[total_K] c_part_x_flat;
  vector[total_K] c_part_y_flat;
  // matrix[N_quad, total_K] Phi_flat;
  // matrix[N_quad, total_K] Phi_d_flat;
  // vector[total_Z] Z_elements;         // Z matrices flattened into 1 long vector

  //array[D] vector[N_quad] t_grid;
  //vector[D] Lt;

  matrix[D,2] drifter_boundaries;
  vector[total_K - total_M] drifter_times;


  // --- Surface Physics Data ---
  int<lower=1> N_models;
  int<lower=1> M_Surface;
  vector[4] border;
  vector<lower=0>[N_models] fixed_ks;
  //vector<lower=0>[N_models] fixed_ls;

}

transformed data {
  real Lx = border[3] - border[1];
  real Ly = border[4] - border[2];

  // --- Ragged Index Trackers ---
  array[D] int drifter_idxs;
  array[D] int M_starts;
  array[D] int K_starts;
  array[D] int Z_starts;
  array[D] int N_starts;
  array[D] vector[N_quad] t_grid;
  vector[D] Lt;

  int cur_m = 1;
  int cur_k = 1;
  int cur_z = 1;
  int cur_n = 1;

  for (d in 1:D) {
    drifter_idxs[d] = d;
    M_starts[d] = cur_m;
    K_starts[d] = cur_k;
    Z_starts[d] = cur_z;
    N_starts[d] = cur_n;

    cur_m += M[d];
    cur_k += K[d];
    cur_z += K[d] * M[d]; // Z is size K x M
    cur_n += K[d] - M[d];

    t_grid[d] = linspaced_vector(N_quad, drifter_boundaries[d,1], drifter_boundaries[d,2]);
    Lt[d] = drifter_boundaries[d,2] - drifter_boundaries[d,1];

  }

  vector[D] dt_weight;
  for (d in 1:D) {
    dt_weight[d] = sqrt(Lt[d] / N_quad);
  }

  // Compute Phi_flat, Phi_d_flat, Z_elements

  matrix[N_quad, total_K] Phi_flat;
  matrix[N_quad, total_K] Phi_d_flat;
  vector[total_Z] Z_elements;

  for(d in 1:D){

    vector[N_quad] cur_quad_grid = t_grid[d];
    real cur_Lt = Lt[d];
    int cur_K = K[d];
    int cur_M = M[d];
    int cur_N = cur_K - cur_M;

    int n_start = N_starts[d];
    int n_end = n_start + cur_N - 1;

    int k_start = K_starts[d];
    int k_end = k_start + cur_K - 1;

    vector[cur_N] cur_drifter_times = drifter_times[n_start:n_end];

    matrix[cur_N, cur_K] cur_P;
    matrix[N_quad, cur_K] cur_Phi;
    matrix[N_quad, cur_K] cur_Phi_d;

    for(k in 1:cur_K){
      cur_P[,k] = cos(((cur_drifter_times - drifter_boundaries[d,1])/cur_Lt) * pi() * (k-1));
      cur_Phi[,k] = cos(((cur_quad_grid - drifter_boundaries[d,1])/cur_Lt) * pi() * (k-1));
      cur_Phi_d[,k] = sin(((cur_quad_grid - drifter_boundaries[d,1])/cur_Lt) * pi() * (k-1)) * -1 * pi() * (k-1) / cur_Lt;
    }

    matrix[cur_K, cur_K] cur_V_full = qr_Q(cur_P');
    matrix[cur_K, cur_M] cur_Z = cur_V_full[,(cur_K - cur_M + 1):cur_K];

    int z_start = Z_starts[d];
    int z_end = z_start + (cur_K * cur_M) - 1;

    Z_elements[z_start:z_end] = to_vector(cur_Z);

    Phi_flat[,k_start:k_end] = cur_Phi;
    Phi_d_flat[,k_start:k_end] = cur_Phi_d;

  }

  // --- Precomputations (Solves Tape Bloat AND Ragged Geometry) ---
  array[D] vector[N_quad] p_x_base;
  array[D] vector[N_quad] p_y_base;
  array[D] vector[N_quad] v_x_base;
  array[D] vector[N_quad] v_y_base;

  matrix[N_quad, total_M] Phi_Z_flat;
  matrix[N_quad, total_M] Phi_d_Z_flat;

  for (d in 1:D) {
    int k_start = K_starts[d];
    int k_end = k_start + K[d] - 1;

    int m_start = M_starts[d];
    int m_end = m_start + M[d] - 1;

    int z_start = Z_starts[d];
    int z_end = z_start + (K[d] * M[d]) - 1;

    // Extract ragged data blocks
    vector[K[d]] c_x_d = c_part_x_flat[k_start:k_end];
    vector[K[d]] c_y_d = c_part_y_flat[k_start:k_end];
    matrix[N_quad, K[d]] Phi_d_mat = Phi_flat[, k_start:k_end];
    matrix[N_quad, K[d]] Phi_d_deriv_mat = Phi_d_flat[, k_start:k_end];

    // Reconstruct 2D Z matrix from the flat vector
    matrix[K[d], M[d]] Z_d = to_matrix(Z_elements[z_start:z_end], K[d], M[d]);

    // Compute the static matrix algebra
    p_x_base[d] = Phi_d_mat * c_x_d;
    p_y_base[d] = Phi_d_mat * c_y_d;
    v_x_base[d] = Phi_d_deriv_mat * c_x_d;
    v_y_base[d] = Phi_d_deriv_mat * c_y_d;

    Phi_Z_flat[, m_start:m_end] = Phi_d_mat * Z_d;
    Phi_d_Z_flat[, m_start:m_end] = Phi_d_deriv_mat * Z_d;
  }

  vector[M_Surface+1] omega_surface;
  for(m in 1:(M_Surface+1)) omega_surface[m] = (m-1) * pi();

}

parameters {
  vector[total_M] w_x;
  vector[total_M] w_y;
  real<lower=0> sigma_vel;
  array[N_models] matrix[M_Surface+1, M_Surface+1] z_surface_betas;
  real<lower=0> fixed_ls;
}

transformed parameters {
  array[N_models] matrix[M_Surface+1, M_Surface+1] prior_std_surface;
  for (i in 1:N_models) {
    vector[M_Surface+1] sqrt_spd = sqrt(sqrt(2*pi()) * fixed_ls * exp(-0.5 * square(fixed_ls * omega_surface)));
    prior_std_surface[i] = fixed_ks[i] * (sqrt_spd * sqrt_spd');
  }

  array[N_models] matrix[M_Surface+1, M_Surface+1] surface_betas;
  for (i in 1:N_models) {
    surface_betas[i] = z_surface_betas[i] .* prior_std_surface[i];
  }
}

model {
  w_x ~ normal(0, 1);
  w_y ~ normal(0, 1);
  sigma_vel ~ normal(0, 1);

  for (i in 1:N_models) {
    // FIX: Target z_surface_betas for the clean N(0,1) prior
    to_vector(z_surface_betas[i]) ~ std_normal();
  }

  int grainsize = 1;

  target += reduce_sum(partial_sum, drifter_idxs, grainsize,
                       w_x, w_y, sigma_vel,
                       p_x_base, p_y_base, v_x_base, v_y_base,
                       Phi_Z_flat, Phi_d_Z_flat,
                       M, M_starts, t_grid, dt_weight,
                       M_Surface, surface_betas, omega_surface,
                       Lx, Ly, border,
                       N_models);
}
