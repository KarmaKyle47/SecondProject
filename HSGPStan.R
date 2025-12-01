library(cmdstanr)
options(mc.cores = parallel::detectCores())
library(MASS)

sampledGMM = sampleGMM(c(-5,-5,5,5))

baseGridBoundaries = generate_grid_tree_boundaries(10)
tree = generate_grid_tree(0.1, c(0,0,1,1))
ModelLogits = sample_new_models_one_pass(tree, 5, 5)
basePatches = sample_NewPatches(ModelLogits, k = 0.5)

library(abind)
baseGridCornerQuantities_Array = do.call(abind, c(basePatches, along = 4))
dim(baseGridCornerQuantities_Array)

# --- 2. Assemble the Stan Data List ---
# Names must match the 'data' block in the .stan file
stan_data <- list(
  # N_mixtures = length(sampledGMM$Cov),
  # GMM_means_data = sampledGMM$Mean,
  # GMM_cov_data = array(unlist(sampledGMM$Cov), dim = c(2, 2, length(sampledGMM$Cov))),
  # GMM_weights_data = sampledGMM$Weights,
  # curPos_Phy_data = c(-3.2,1.4),
  # baseBoundaries_Data = as.matrix(baseTree$boundaries),
  # baseModels_Data = as.matrix(sampleNewModels),
  # t_data = 0,
  # fullBoundaries_data = as.matrix(updatedTree$boundaries),
  # fullCoefs1 = updatedTree$coefs[[2]],
  # fullCoefs2 = updatedTree$coefs[[3]],
  comp_res = 10,
  trans_prop = 8/9,
  baseGridCornerQuantities = baseGridCornerQuantities_Array,
  baseGridBoundaries = baseGridBoundaries,
  model_num = 4,
  selfPenalty = 10,
  ModelLogits = ModelLogits,

  y_dummy = 0.0 # Just a placeholder
)

baseVectorFields = function(t, curPos){

  f1 = c(curPos[2],-1*curPos[1])/sqrt(sum(c(curPos[1],curPos[2])^2))
  f2 = c(curPos[1],curPos[2])/sqrt(sum(c(curPos[1],curPos[2])^2))
  #f3 = c(sqrt(2),sqrt(2))
  matrix(c(f1,f2), nrow = 2, byrow = F)

}

visualizeNewModelExistence(tree$boundaries, sampled_ModelProbs, model = 2)

sampleGMM(c(-2,-2,2,2), lambda = 1)

GMM_means = matrix(c(0,0,1,1), nrow = 2, byrow = T)
GMM_cov = array(c(0.1,0,0,0.1,0.1,0,0,0.1), dim = c(2,2,2))
GMM_weights = c(0.5,0.5)

sampled_GMM = list(Mean = GMM_means,
                   Cov = asplit(GMM_cov, MARGIN = 3),
                   Weights = GMM_weights)

M=10
z_log1 = matrix(rnorm((M+1)*(M+1)), nrow = M+1, ncol = M+1)
z_log2 = matrix(rnorm((M+1)*(M+1)), nrow = M+1, ncol = M+1)
k_log = 0.175
l_log = 0.2
z_logit1 = matrix(rnorm((M+1)*(M+1)), nrow = M+1, ncol = M+1)
z_logit2 = matrix(rnorm((M+1)*(M+1)), nrow = M+1, ncol = M+1)
k_logit = 2
l_logit = 0.1
logit_prior_mean = 10


Data_t = rep(0,1000)
Data_pos = matrix(nrow = 1000, ncol = 2)
Data_vel = matrix(nrow = 1000, ncol = 2)

for(i in 1:1000){

  which_normal = sample(c(1,2), 1, prob = GMM_weights)
  Data_pos[i,] = mvrnorm(mu = GMM_means[which_normal,], Sigma = GMM_cov[,,which_normal])

  cur_compPos = get_compSpace_pos(sampled_GMM, Data_pos[i,])

  curLogSurfaceVal_1 = evaluateHSGP(z = z_log1, k = k_log, l = l_log, M = 10, curPos = cur_compPos)
  curLogitSurfaceVal_1 = evaluateHSGP(z = z_logit1, k = k_logit, l = l_logit, M = 10, curPos = cur_compPos) + 10
  curLogSurfaceVal_2 = evaluateHSGP(z = z_log2, k = k_log, l = l_log, M = 10, curPos = cur_compPos)
  curLogitSurfaceVal_2 = evaluateHSGP(z = z_logit2, k = k_logit, l = l_logit, M = 10, curPos = cur_compPos) + 10

  curTrajVal = c(exp(curLogSurfaceVal_1)*invlogit(curLogitSurfaceVal_1),
                 exp(curLogSurfaceVal_2)*invlogit(curLogitSurfaceVal_2))

  Data_vel[i,] = t(baseVectorFields(0,Data_pos[i,]) %*% matrix(curTrajVal)) + rnorm(2,0,0.1)

}


plot(Data_pos)
plot(Data_vel)

plot(t(apply(Data_pos, MARGIN = 1, FUN = get_compSpace_pos, GMM = sampled_GMM)))

Stan_Data = matrix(c(Data_t, Data_pos, Data_vel), nrow = 1000, ncol = 5, byrow = F)

stan_data = list(
  N_data = 1000,                 # Number of data points
  Data = Stan_Data,                 # Particle Velocities with Positions for now (t, x, y, v_x, v_y)
  GMM_num = 2,              # Number of Gaussian Mixtures in the Transformation
  GMM_means = GMM_means,     # Means of all Gaussian Mixtures;
  GMM_cov = GMM_cov,        # Covariance matrices of all Mixtures
  GMM_weights = GMM_weights, #Weights for Mixtures
  M = 10 #Number of eigenfunctions
)

init_fun <- function() {
  list(
    # Start Magnitudes small (Assume physics is correct, no scaling needed yet)
    coef_ks = rep(0.1, N_models),
    logit_ks = rep(0.5, N_models),

    # Start Length Scales at a "safe" middle ground (not 0.05, not 100)
    coef_ls = rep(0.2, N_models),
    logit_ls = rep(0.15, N_models),

    # Start Surfaces FLAT (Crucial!)
    coef_zs = array(0, dim = c(N_models, M+1, M+1)),
    logit_zs = array(0, dim = c(N_models, M+1, M+1)),

    # Start noise reasonable
    sigma_vel = 0.5
  )
}

cat("--- STARTING GPU BENCHMARK ---\n")
t_start <- Sys.time()
fit_gpu <- mod$sample(
  data = stan_data,
  chains = 1,
  iter_warmup = 10,
  iter_sampling = 0,
  init = init_fun, # Use your init function
  refresh = 1
)
t_end <- Sys.time()
print(t_end - t_start)

# --- 3. Compile the Stan Model ---
# This will take a minute, but it uses the REAL compiler
mod <- rstan::stan_model("STAN_Files/SurfaceTrajectoryOnlyFunctions.stan")
mod <- rstan::stan_model("STAN_Files/test_functions.stan")
mod <- rstan::stan_model("STAN_Files/SurfaceTrajectory.stan")

cuda_path <- "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v13.0"

opencl_flags <- list(
  stan_opencl = TRUE,
  opencl_platform = 1,
  opencl_device = 0,
  LDFLAGS = paste0("-L\"", cuda_path, "/lib/x64\" -lOpenCL"),
  CXXFLAGS = paste0("-I\"", cuda_path, "/include\"")
)

opencl_flags <- list(
  stan_opencl = TRUE,
  opencl_platform = 1,
  opencl_device = 0
)

mod <- cmdstan_model(
  "STAN_Files/SurfaceTrajectoryHSGP.stan",
  cpp_options = opencl_flags,
  force_recompile = TRUE
)

mod <- cmdstan_model(
  "STAN_Files/SurfaceTrajectoryHSGP.stan",
  force_recompile = TRUE
)


mod <- cmdstan_model("STAN_Files/SurfaceTrajectoryHSGP.stan")

fit <- mod$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 1000
)

fit$profiles()
fit$cmdstan_diagnose()

fit_sum = fit$summary()
max(fit_sum$rhat)

# 1. Extract all draws for the coefficient Zs
# format: [iterations, chains, flat_index]
log_draws <- fit$draws("log_zs", format = "draws_matrix")
logit_draws <- fit$draws("logit_zs", format = "draws_matrix")

log_draws_array_1 = array(dim = c(M+1, M+1, 4000))
log_draws_array_2 = array(dim = c(M+1, M+1, 4000))

logit_draws_array_1 = array(dim = c(M+1, M+1, 4000))
logit_draws_array_2 = array(dim = c(M+1, M+1, 4000))

for(i in 1:4000){

  log_draws_array_1[,,i] = matrix(log_draws[i,1:121*2 - 1], nrow = M+1, ncol = M+1, byrow = F)
  log_draws_array_2[,,i] = matrix(log_draws[i,1:121*2], nrow = M+1, ncol = M+1, byrow = F)

  logit_draws_array_1[,,i] = matrix(logit_draws[i,1:121*2 - 1], nrow = M+1, ncol = M+1, byrow = F)
  logit_draws_array_2[,,i] = matrix(logit_draws[i,1:121*2], nrow = M+1, ncol = M+1, byrow = F)

}

grid_res = 100
n_iters = 4000

getGridStanHSGP = function(fit, grid_res, M, n_iters){

  log_draws <- fit$draws("log_zs", format = "draws_matrix")
  logit_draws <- fit$draws("logit_zs", format = "draws_matrix")

  log_draws_array_1 = array(dim = c(M+1, M+1, n_iters))
  log_draws_array_2 = array(dim = c(M+1, M+1, n_iters))

  logit_draws_array_1 = array(dim = c(M+1, M+1, n_iters))
  logit_draws_array_2 = array(dim = c(M+1, M+1, n_iters))

  for(i in 1:n_iters){

    log_draws_array_1[,,i] = matrix(log_draws[i,1:((M+1)^2)*2 - 1], nrow = M+1, ncol = M+1, byrow = F)
    log_draws_array_2[,,i] = matrix(log_draws[i,1:((M+1)^2)*2], nrow = M+1, ncol = M+1, byrow = F)

    logit_draws_array_1[,,i] = matrix(logit_draws[i,1:((M+1)^2)*2 - 1], nrow = M+1, ncol = M+1, byrow = F)
    logit_draws_array_2[,,i] = matrix(logit_draws[i,1:((M+1)^2)*2], nrow = M+1, ncol = M+1, byrow = F)

  }

  x_seq = y_seq = seq(0,1,length.out = grid_res+1)

  plot_grid = expand.grid(x_seq, y_seq)

  coef_grid_draws_1 = array(dim = c(grid_res+1, grid_res+1, n_iters))
  weight_grid_draws_1 = array(dim = c(grid_res+1, grid_res+1, n_iters))

  coef_grid_draws_2 = array(dim = c(grid_res+1, grid_res+1, n_iters))
  weight_grid_draws_2 = array(dim = c(grid_res+1, grid_res+1, n_iters))


  for(i in 1:n_iters){
    coef_grid_draws_1[,,i] = matrix(exp(apply(plot_grid, MARGIN = 1, FUN = evaluateHSGP, z = log_draws_array_1[,,i], k= k_log, l=l_log, M=M)), byrow = F, nrow = grid_res+1)
    weight_grid_draws_1[,,i] = matrix(invlogit(apply(plot_grid, MARGIN = 1, FUN = evaluateHSGP, z = logit_draws_array_1[,,i], k= k_logit, l=l_logit, M=M)), byrow = F, nrow = grid_res+1)
    coef_grid_draws_2[,,i] = matrix(exp(apply(plot_grid, MARGIN = 1, FUN = evaluateHSGP, z = log_draws_array_2[,,i], k= k_log, l=l_log, M=M)), byrow = F, nrow = grid_res+1)
    weight_grid_draws_2[,,i] = matrix(invlogit(apply(plot_grid, MARGIN = 1, FUN = evaluateHSGP, z = logit_draws_array_2[,,i], k= k_logit, l=l_logit, M=M)), byrow = F, nrow = grid_res+1)
    svMisc::progress(i,4000)
  }



}
n_iters = 4000
getGridStanHSGP = function(fit, grid_res, M, n_iters){

  # --- 1. Setup Constants ---
  # Number of coefficients per model
  N_basis = (M+1)^2
  omega = (0:M) * pi

  # Create Grid Vectors (0 to 1)
  x_seq = seq(0, 1, length.out = grid_res + 1)
  y_seq = seq(0, 1, length.out = grid_res + 1)

  # --- 2. Pre-calculate Basis Matrices (The Speedup) ---
  # We do this ONCE, not inside the loop.
  # Phi_x: [(grid_res+1) x (M+1)]

  # Helper to make cosine matrix
  make_phi <- function(coords) {
    # n=0 is column of 1s
    # n>0 is sqrt(2)*cos(...)
    # outer() creates the matrix of coords vs frequencies
    waves <- sqrt(2) * cos(outer(coords, omega[2:(M+1)]))
    return(cbind(1, waves))
  }

  Phi_x <- make_phi(x_seq)
  Phi_y <- make_phi(y_seq)

  # --- 3. Extract Hyperparameters ---
  # We need the specific k and l for every draw to scale the Zs correctly
  # Assuming N_models = 2

  # Coefficient (Log) Hypers
  log_ks_1 <- as.vector(fit$draws("log_ks[1]", format = "draws_matrix")[1:n_iters])
  log_ls_1 <- as.vector(fit$draws("log_ls[1]", format = "draws_matrix")[1:n_iters])
  log_ks_2 <- as.vector(fit$draws("log_ks[2]", format = "draws_matrix")[1:n_iters])
  log_ls_2 <- as.vector(fit$draws("log_ls[2]", format = "draws_matrix")[1:n_iters])

  # log_ks_1 <- rep(k_log, n_iters)
  # log_ls_1 <- rep(l_log, n_iters)
  # log_ks_2 <- rep(k_log, n_iters)
  # log_ls_2 <- rep(l_log, n_iters)

  # Weight (Logit) Hypers
  logit_ks_1 <- as.vector(fit$draws("logit_ks[1]", format = "draws_matrix")[1:n_iters])
  logit_ls_1 <- as.vector(fit$draws("logit_ls[1]", format = "draws_matrix")[1:n_iters])
  logit_ks_2 <- as.vector(fit$draws("logit_ks[2]", format = "draws_matrix")[1:n_iters])
  logit_ls_2 <- as.vector(fit$draws("logit_ls[2]", format = "draws_matrix")[1:n_iters])

  # logit_ks_1 <- rep(k_logit, n_iters)
  # logit_ls_1 <- rep(l_logit, n_iters)
  # logit_ks_2 <- rep(k_logit, n_iters)
  # logit_ls_2 <- rep(l_logit, n_iters)

  # --- 4. Extract Raw Z Matrices ---
  # Safer extraction using variable names

  log_draws <- fit$draws("log_zs", format = "draws_matrix")
  logit_draws <- fit$draws("logit_zs", format = "draws_matrix")

  log_draws_1 = array(dim = c(M+1, M+1, n_iters))
  log_draws_2 = array(dim = c(M+1, M+1, n_iters))

  logit_draws_1 = array(dim = c(M+1, M+1, n_iters))
  logit_draws_2 = array(dim = c(M+1, M+1, n_iters))

  for(i in 1:n_iters){

    log_draws_1[,,i] = matrix(log_draws[i,1:((M+1)^2)*2 - 1], nrow = M+1, ncol = M+1, byrow = F)
    log_draws_2[,,i] = matrix(log_draws[i,1:((M+1)^2)*2], nrow = M+1, ncol = M+1, byrow = F)

    logit_draws_1[,,i] = matrix(logit_draws[i,1:((M+1)^2)*2 - 1], nrow = M+1, ncol = M+1, byrow = F)
    logit_draws_2[,,i] = matrix(logit_draws[i,1:((M+1)^2)*2], nrow = M+1, ncol = M+1, byrow = F)

  }

  # --- 5. Initialize Output Arrays ---
  coef_grid_draws_1 = array(dim = c(grid_res+1, grid_res+1, n_iters))
  weight_grid_draws_1 = array(dim = c(grid_res+1, grid_res+1, n_iters))
  traj_grid_draws_1 = array(dim = c(grid_res+1, grid_res+1, n_iters))

  coef_grid_draws_2 = array(dim = c(grid_res+1, grid_res+1, n_iters))
  weight_grid_draws_2 = array(dim = c(grid_res+1, grid_res+1, n_iters))
  traj_grid_draws_2 = array(dim = c(grid_res+1, grid_res+1, n_iters))

  # Helper for Spectral Density calculation
  get_beta_matrix <- function(z_flat, k, l, M) {
    # Calculate Spectral Density for this specific length scale
    # Squared Exponential Kernel
    # Note: Using Unit Variance first (k=1.0 inside), multiply by k at end
    spd_unit <- sqrt(2*pi) * l * exp(-0.5 * l^2 * omega^2)
    scale_vec <- sqrt(spd_unit)

    # Reshape Z to matrix
    z_mat <- matrix(z_flat, nrow = M+1, ncol = M+1)

    # Apply Scaling: D * Z * D
    beta_unit <- diag(scale_vec) %*% z_mat %*% diag(scale_vec)

    return(k * beta_unit)
  }

  # --- 6. The Main Loop (Now Fast) ---
  # Uses matrix multiplication instead of 'apply'

  print("Reconstructing surfaces...")
  pb <- txtProgressBar(min = 0, max = n_iters, style = 3)

  for(i in 1:n_iters){

    # --- Model 1 ---
    # A. Coefficients
    beta_log_1 <- get_beta_matrix(log_draws_1[,,i], log_ks_1[i], log_ls_1[i], M)
    # Tensor Projection: Phi_x * Beta * Phi_y'
    surf_log_1 <- Phi_x %*% beta_log_1 %*% t(Phi_y)
    coef_grid_draws_1[,,i] <- exp(surf_log_1) # Exponential Link

    # B. Weights
    beta_logit_1 <- get_beta_matrix(logit_draws_1[,,i], logit_ks_1[i], logit_ls_1[i], M)
    surf_logit_1 <- Phi_x %*% beta_logit_1 %*% t(Phi_y)
    # Add prior mean shift for logits if needed (e.g., -2 or +0.69)
    # assuming prior_mean is handled implicitly or added here
    weight_grid_draws_1[,,i] <- invlogit(surf_logit_1 + 10) # Sigmoid Link

    traj_grid_draws_1[,,i] = coef_grid_draws_1[,,i]*weight_grid_draws_1[,,i]

    # --- Model 2 ---
    # A. Coefficients
    beta_log_2 <- get_beta_matrix(log_draws_2[,,i], log_ks_2[i], log_ls_2[i], M)
    surf_log_2 <- Phi_x %*% beta_log_2 %*% t(Phi_y)
    coef_grid_draws_2[,,i] <- exp(surf_log_2)

    # B. Weights
    beta_logit_2 <- get_beta_matrix(logit_draws_2[,,i], logit_ks_2[i], logit_ls_2[i], M)
    surf_logit_2 <- Phi_x %*% beta_logit_2 %*% t(Phi_y)
    weight_grid_draws_2[,,i] <- invlogit(surf_logit_2 + 10)

    traj_grid_draws_2[,,i] = coef_grid_draws_2[,,i]*weight_grid_draws_2[,,i]

    setTxtProgressBar(pb, i)
  }
  close(pb)

  # Return list of results
  return(list(
    coef_1 = coef_grid_draws_1,
    weight_1 = weight_grid_draws_1,
    traj_1 = traj_grid_draws_1,
    coef_2 = coef_grid_draws_2,
    weight_2 = weight_grid_draws_2,
    traj_2 = traj_grid_draws_2,
    x = x_seq,
    y = y_seq
  ))
}


fast_mean_3d <- function(arr_3d) {
  d <- dim(arr_3d)
  dim(arr_3d) <- c(d[1] * d[2], d[3])
  means <- rowMeans(arr_3d)
  dim(means) <- c(d[1], d[2])
  return(means)
}

testSTANGrid = getGridStanHSGP(fit, 100, M = 10, n_iters = 4000)
truth_1 = plotHSGP(grid_res = 100, z_log = z_log1, k_log = k_log, l_log = l_log,
                   z_logit = z_logit1, k_logit = k_logit, l_logit = l_logit, logit_prior_mean = 10, M = 10, color_limits = c(0,2))
truth_2 = plotHSGP(grid_res = 100, z_log = z_log2, k_log = k_log, l_log = l_log,
                   z_logit = z_logit2, k_logit = k_logit, l_logit = l_logit, logit_prior_mean = 10, M = 10, color_limits = c(0,2))



# Usage
x_seq = y_seq = seq(0,1,length.out = 101)
final_mean_traj1 <- fast_mean_3d(testSTANGrid$traj_1)
final_mean_traj2 <- fast_mean_3d(testSTANGrid$traj_2)

plot_ly(x = x_seq, y = y_seq, z = final_mean_traj1, type = "surface", color_scale = "Viridis", cmin = 0, cmax = 2)
truth_1$Trajectory

plot_ly(x = x_seq, y = y_seq, z = final_mean_traj2, type = "surface", color_scale = "Viridis", cmin = 0, cmax = 2)
truth_2$Trajectory


log_means_flat <- colMeans(log_draws)
log_means_flat_1 = log_means_flat[1:121 * 2 -1]
log_means_flat_2 = log_means_flat[1:121 * 2]

z_log_STAN1 = matrix(log_means_flat_1, nrow = 11, ncol = 11, byrow = F)
z_log_STAN1 = matrix(log_means_flat_2, nrow = 11, ncol = 11, byrow = F)

# Do the same for Logits
logit_draws <- fit$draws("logit_zs", format = "draws_matrix")
logit_means_flat <- colMeans(logit_draws)
logit_means_flat_1 = logit_means_flat[1:121 * 2 -1]
logit_means_flat_2 = logit_means_flat[1:121 * 2]

z_logit_STAN1 = matrix(logit_means_flat_1, nrow = 11, ncol = 11, byrow = F)
z_logit_STAN1 = matrix(logit_means_flat_2, nrow = 11, ncol = 11, byrow = F)

log_ks_draws = fit$draws("log_ks", format = "draws_matrix")
log_ks <- colMeans(log_ks_draws)

log_ks = rep(k_log,2)

logit_ks_draws = fit$draws("logit_ks", format = "draws_matrix")
logit_ks <- colMeans(logit_ks_draws)

logit_ks = rep(k_logit,2)

log_ls_draws = fit$draws("log_ls", format = "draws_matrix")
log_ls <- colMeans(log_ls_draws)

log_ls = rep(l_log, 2)

logit_ls_draws = fit$draws("logit_ls", format = "draws_matrix")
logit_ls <- colMeans(logit_ls_draws)

logit_ls = rep(l_logit, 2)

sigma_vel_draws = fit$draws("sigma_vel", format = "draws_matrix")
sigma_vel <- colMeans(sigma_vel_draws)

bayes_test_1 = plotHSGP(grid_res = 100, z_log = z_log_STAN1, k_log = log_ks[1], l_log = log_ls[1],
                        z_logit = z_logit_array[,,1], k_logit = z_logit_STAN1, l_logit = logit_ls[1], logit_prior_mean = 10, M = 10, color_limits = c(0,2))
bayes_test_2 = plotHSGP(grid_res = 100, z_log = z_log_array[,,2], k_log = log_ks[2], l_log = log_ls[2],
                        z_logit = z_logit_array[,,2], k_logit = logit_ks[2], l_logit = logit_ls[2], logit_prior_mean = 10, M = 10, color_limits = c(0,2))

truth_1 = plotHSGP(grid_res = 100, z_log = z_log1, k_log = k_log, l_log = l_log,
                   z_logit = z_logit1, k_logit = k_logit, l_logit = l_logit, logit_prior_mean = 10, M = 10, color_limits = c(0,2))
truth_2 = plotHSGP(grid_res = 100, z_log = z_log2, k_log = k_log, l_log = l_log,
                   z_logit = z_logit2, k_logit = k_logit, l_logit = l_logit, logit_prior_mean = 10, M = 10, color_limits = c(0,2))

truth_1$Trajectory %>% layout(
  scene = list(
    zaxis = list(range = c(0,2)) # <--- THIS IS THE FIX
  ))
bayes_test_1$Trajectory %>% layout(
  scene = list(
    zaxis = list(range = c(0,2)) # <--- THIS IS THE FIX
  ))

truth_1$Coefficient
truth_2$Trajectory
bayes_test_2$Trajectory
bayes_test_2$Weight






M <- 10
N_models <- 2


fit <- mod$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  refresh = 10, # print update every 500 iters
  init = init_fun
)

fit$profiles()
fit$cmdstan_diagnose()

fit_sum = fit$summary()
max(fit_sum$rhat)
library(bayesplot)
mcmc_trace(fit$draws(c("coef_ls", "coef_ks")))

# --- 4. Run the Model to Test the Function ---
# This part is very fast
fit <- rstan::sampling(
  object = mod,
  data = stan_data,
  algorithm = "NUTS",  # <-- The magic trick
  iter = 1,                   # We only need one "draw"
  warmup = 0,
  chains = 4
)

fit <- rstan::sampling(
  object = mod,
  data = stan_data,
  algorithm = "NUTS",
  chains = 8,
  warmup = 1000,
  iter = 5000
)

?sampling

# --- 5. Extract and Inspect the Result ---
fit_extract <- rstan::extract(fit)
rstan::extract()
print(fit_extract)
Rhat(fit)
plot(fit_extract$modelProbs[,1,], type = 'l')
plot(fit_extract$modelProbs[,4,3], type = 'l')
hist(fit_extract$modelProbs[,4,2])
hist(fit_extract$sigma_vel)
hist(fit_extract$UpdatedCornerQuantities[,4,1,1])
median(fit_extract$UpdatedCornerQuantities[,1,1,1])

sampledTreePlot$ModelRegions[[1]]

save(fit_extract, stan_data, sampledTree, sampledTreePlot, file = "FirstTestSTAN.RData")

plot(fit_extract$modelLogits[,4,2])
plot(fit_extract$baseCornerQuanities[,1,1,1,2])
dim(fit_extract$baseCornerQuanities)

# 'test_output' is what we named our variable in 'generated quantities'
fit_extract$test_baseBoundaries[1,,] == baseTreeSampled[,1:4]# Get the first (and only) row
fit_extract$test_transBoundaries[1,,] == baseTreeSampled[,1:4]# Get the first (and only) row
fit_extract$test_updatedQuantities[1,,]
fit_extract$test_updatedCoefs[1,,]
fit_extract$test_BaseEnergy
print("Function test output:")
print(final_result)

get_compSpace_pos(GMM = sampledGMM, curPos_phy = c(-3.2,1.4))
baseTree$models = sampleNewModels
baseTreeSampled = baseTree$boundaries
baseTreeSampled$model1 = sampleNewModels[,1]
baseTreeSampled$model2 = sampleNewModels[,2]
calculate_tree_energy(baseTreeSampled, self_penalty = 100000, temperature = 1)

updatedTree_V2 = updatedTree
updatedTree_V2$coefs = list(updatedTree$coefs[[2]], updatedTree$coefs[[3]])

TrajWeightedBaseVectorFields(0, c(-3.2,1.4), baseVectorFields, compPatchTree = updatedTree_V2, GMM = sampledGMM)

plotTransitionRegions(baseTreeSampled)
newBoundaries

median(fit$draws("coef_ls")[,,2])
system("clinfo")

# Generate synthetic data for scaling test
run_benchmark <- function(N_test) {

  cat(paste0("\n--- TESTING N = ", N_test, " ---\n"))

  # Fake data structures
  bench_data <- list(
    N_data = N_test,
    Data = matrix(rnorm(N_test * 5), ncol = 5),
    GMM_num = 2,
    GMM_means = matrix(0, 2, 2),
    GMM_cov = array(diag(2), dim=c(2,2,2)),
    GMM_weights = c(0.5, 0.5),
    M = 10 # Keep M constant
  )

  # Run just 10 iterations to check speed
  t_start <- Sys.time()
  fit <- mod$sample(
    data = bench_data,
    chains = 1,
    iter_warmup = 10,
    iter_sampling = 0,
    fixed_param = TRUE, # Just measure likelihood eval speed
    init = init_fun
  )
  print(Sys.time() - t_start)
}

# 1. The "Slow" Zone (Where you are now)
run_benchmark(1000)

# 2. The Break-Even Zone
run_benchmark(10000)

# 3. The "GPU Shine" Zone
run_benchmark(1000000)



fit <- mod$sample(
  data = stan_data,
  init = init_fun,
  chains = 1,           # Just 1 chain to save time/heat
  iter_warmup = 200,    # Short warmup
  iter_sampling = 200,  # Short sampling
  fixed_param = FALSE,  # <--- TURN THE ENGINE ON
  refresh = 1
)

summary_stats <- fit$summary(
  variables = c("coef_ks", "logit_ks", "coef_ls", "logit_ls", "sigma_vel"),
  "mean", "sd", "rhat" # What columns you want
)

print(summary_stats, n = Inf)









