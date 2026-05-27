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

  // --- 2. Vectorized Base VF Evaluator ---
  array[] matrix evaluate_VF_vec(vector t, vector x, vector y,
                                 int N_points, int N_models){

    vector[N_points] window = cos(pi()*x/20.0) .* cos(pi()*y/20.0);

    array[N_models] matrix[N_points, 2] VFs;

    VFs[1][,1] = y .* window;
    VFs[1][,2] = -1.0 * x .* window;

    VFs[2][,1] = x .* window .* sin(t/8.0);
    VFs[2][,2] = y .* window .* sin(t/8.0);

    return VFs;
  }

  // --- 3. Parallel Likelihood Function ---
  real partial_sum(array[] int point_slice, int start, int end,
                   vector x, vector y, vector v_x, vector v_y, vector t,
                   real sigma_vel, int M_Surface, array[] matrix surface_betas,
                   vector omega_surface, real Lx, real Ly, vector border, int N_models) {

    int N_slice = end - start + 1;

    // Slice the data vectors for this specific thread
    vector[N_slice] x_slice = x[start:end];
    vector[N_slice] y_slice = y[start:end];
    vector[N_slice] vx_slice = v_x[start:end];
    vector[N_slice] vy_slice = v_y[start:end];
    vector[N_slice] t_slice = t[start:end];

    // Evaluate the Base Vector Fields for the slice
    array[N_models] matrix[N_slice, 2] VFs = evaluate_VF_vec(t_slice, x_slice, y_slice, N_slice, N_models);

    // Evaluate the HSGP Surfaces for the slice
    vector[N_slice] coef1 = exp(evaluate_HSGP_Surface_vec(N_slice, x_slice, y_slice, M_Surface, surface_betas[1], omega_surface, Lx, Ly, border));
    vector[N_slice] target_v_x = VFs[1][, 1] .* coef1;
    vector[N_slice] target_v_y = VFs[1][, 2] .* coef1;

    vector[N_slice] coef2 = exp(evaluate_HSGP_Surface_vec(N_slice, x_slice, y_slice, M_Surface, surface_betas[2], omega_surface, Lx, Ly, border));
    target_v_x += VFs[2][, 1] .* coef2;
    target_v_y += VFs[2][, 2] .* coef2;

    // Return the sum of the log-likelihoods for this slice
    return normal_lpdf(vx_slice | target_v_x, sigma_vel) +
           normal_lpdf(vy_slice | target_v_y, sigma_vel);
  }
}

data {
  int<lower=1> N_points;
  vector[N_points] x;
  vector[N_points] y;
  vector[N_points] v_x;
  vector[N_points] v_y;
  vector[N_points] t;

  int<lower=1> N_models;
  int<lower=1> M_Surface;
  vector[4] border;
  vector<lower=0>[N_models] fixed_ks;
  real<lower=0> fixed_ls;
}

transformed data {
  real Lx = border[3] - border[1];
  real Ly = border[4] - border[2];

  vector[M_Surface+1] omega_surface;
  for(m in 1:(M_Surface+1)) omega_surface[m] = (m-1) * pi();

  // --- Index array for reduce_sum slicing ---
  array[N_points] int point_idxs;
  for (n in 1:N_points) {
    point_idxs[n] = n;
  }

  array[N_models] matrix[M_Surface+1, M_Surface+1] prior_std_surface;
  for (i in 1:N_models) {
    vector[M_Surface+1] sqrt_spd = sqrt(sqrt(2.0*pi()) * fixed_ls * exp(-0.5 * square(fixed_ls * omega_surface)));
    prior_std_surface[i] = fixed_ks[i] * (sqrt_spd * sqrt_spd');
  }
}

parameters {
  real<lower=0> sigma_vel;
  array[N_models] matrix[M_Surface+1, M_Surface+1] z_surface_betas;
}

transformed parameters {

  array[N_models] matrix[M_Surface+1, M_Surface+1] surface_betas;
  for (i in 1:N_models) {
    surface_betas[i] = z_surface_betas[i] .* prior_std_surface[i];
  }
}

model {
  // Priors
  sigma_vel ~ normal(0, 1);
  for (i in 1:N_models) {
    to_vector(z_surface_betas[i]) ~ std_normal();
  }

  // Parallel Likelihood via reduce_sum
  int grainsize = 1; // Stan's TBB will auto-tune batching if set to 1
  target += reduce_sum(partial_sum, point_idxs, grainsize,
                       x, y, v_x, v_y, t,
                       sigma_vel, M_Surface, surface_betas,
                       omega_surface, Lx, Ly, border, N_models);
}
