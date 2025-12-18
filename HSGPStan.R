library(cmdstanr)
options(mc.cores = parallel::detectCores())
library(MASS)
library(mclust)
library(stringr)

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

  Data_pos[i,] = runif(2, min = -5, max = 5)
  Data_vel[i,] = TrajWeightedBaseVectorFields_HSGP(0, Data_pos[i,], baseVectorFields, sampledHSGP, logit_prior_mean = 10, M = 10, border = testBorder) + rnorm(2,0,0.1)

}


plot(Data_pos)
plot(Data_vel)

plot(t(apply(Data_pos, MARGIN = 1, FUN = get_compSpace_pos, GMM = sampled_GMM)))

Stan_Data = matrix(c(Data_t, Data_pos, Data_vel), nrow = 1000, ncol = 5, byrow = F)


dataGMM = Mclust(sampledParticles_Subset[,2:3], G = 10)

dataGMM_STAN = list(Mean = t(dataGMM$parameters$mean),
                    Cov = asplit(dataGMM$parameters$variance$sigma, MARGIN = 3),
                    Weights = dataGMM$parameters$pro)

plot(t(apply(sampledParticles[,2:3], MARGIN = 1, FUN = get_compSpace_pos, GMM = dataGMM_STAN)))
plot(sampledParticles[,2:3])

testGMM$parameters

obs_per_drifter = c()

for(i in 1:100){

  cur_drifter = str_c("Particle",i)

  obs_per_drifter = c(obs_per_drifter, sum(sampledParticles_Subset$Particle == cur_drifter))


}

testBorder = border

Lx = testBorder[3] - testBorder[1]
Ly = testBorder[4] - testBorder[2]



stan_data = list(
  N_data = nrow(sampledParticles_Subset[,1:5]),                 # Number of data points
  # N_drifters = 100, # Number of drifters
  # obs_per_drifter = obs_per_drifter,
  Data = sampledParticles_Subset[,1:5],                 # Particle Positions (t, x, y, vx, vy)
  M = 10, #Number of eigenfunctions
  border = c(-5,-5,5,5)
)

init_fun <- function() {
  list(
    # Start Magnitudes small (Assume physics is correct, no scaling needed yet)
    log_ks = rep(0.1, N_models),
    logit_ks = rep(0.5, N_models),

    # Start Length Scales at a "safe" middle ground (not 0.05, not 100)
    log_ls = rep(0.2, N_models),
    logit_ls = rep(0.15, N_models),

    # Start Surfaces FLAT (Crucial!)
    log_zs = array(0, dim = c(N_models, M+1, M+1)),
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

N_models = 2
M = 10

mod <- cmdstan_model("STAN_Files/SurfaceTrajectoryHSGP.stan")
mod <- cmdstan_model("STAN_Files/GeminiOptimized.stan")
mod <- cmdstan_model("STAN_Files/PathHSGPRegression.stan")
mod_Path <- cmdstan_model("STAN_Files/PathHSGPRegressionwLandK.stan")

N_drifters = 100
drifter_prior_means = matrix(nrow = 100, ncol = 2)
drifter_boundaries = matrix(nrow = 100, ncol = 2)
drifter_ks = matrix(nrow = 100, ncol = 2)
drifter_ls = matrix(nrow = 100, ncol = 2)
i=1

for(i in 1:N_drifters){

  curDrifter = sampledParticles[sampledParticles$Particle == str_c("Particle",i),]

  cur_prior_means = c(mean(curDrifter[,2]), mean(curDrifter[,3]))
  cur_boundary = c(min(curDrifter[,1]), max(curDrifter[,1]))

  curDrifter_data = list(
    N_data = nrow(curDrifter),
    Data = curDrifter[,1:3],
    M_Drifter = 10,
    drifter_prior_means = cur_prior_means,
    drifter_boundaries = cur_boundary
  )

  cur_fit <- mod_Path$optimize(data = curDrifter_data, iter = 10000)

  drifter_ks[i,] = as.numeric(cur_fit$draws('drifter_ks'))
  drifter_ls[i,] = as.numeric(cur_fit$draws('drifter_ls'))

  drifter_prior_means[i,] = cur_prior_means
  drifter_boundaries[i,] = cur_boundary

}

4/(pi*0.05)

load(file = "FirstSuccessfulHSGP.RData")

stan_data = list(

  N_data = nrow(sampledParticles),
  N_drifters = 100,
  N_models = 2,
  obs_per_drifter = rep(100,100),
  Data = sampledParticles[,1:3],

  M_Surface = 25,
  M_Drifter = 10,

  log_ks = rep(0.5,2),
  log_ls = rep(0.05,2),
  border = c(-10,-10,10,10),

  drifter_ks = drifter_ks,
  drifter_ls = drifter_ls,

  drifter_prior_means = drifter_prior_means,
  drifter_boundaries = drifter_boundaries,
  t_res = 20

)

mod <- cmdstan_model("STAN_Files/GeminiOptimized.stan", cpp_options = list(stan_threads = TRUE))

fit <- mod$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 7,
  iter_warmup = 500,
  iter_sampling = 1000, refresh = 10
)


stan_data = list(
  N_data = nrow(sampledParticles),                 # Number of data points
  Data = sampledParticles[,1:3],                 # Particle Positions (t, x, y, vx, vy)
  M_Drifter = 35, #Number of eigenfunctions
  drifter_prior_means = c(mean(sampledParticles[,2]), mean(sampledParticles[,3])),
  drifter_boundaries = c(min(sampledParticles[,1]), max(sampledParticles[,1])),
  drifter_ks = c(6.85,8.78),
  drifter_ls = c(0.04,0.04)
)

stan_data = list(
  N_data = nrow(sampledParticles),                 # Number of data points
  Data = sampledParticles[,1:3],                 # Particle Positions (t, x, y, vx, vy)
  M_Drifter = 50, #Number of eigenfunctions
  drifter_prior_means = c(mean(sampledParticles[,2]), mean(sampledParticles[,3])),
  drifter_boundaries = c(min(sampledParticles[,1]), max(sampledParticles[,1]))
)

fit <- mod$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 1000, refresh = 100
)

fit <- mod$optimize(data = stan_data, iter = 10000)
fit <- mod$variational(data = stan_data)

fit$draws()
fit$output()
fit$profiles()
fit$cmdstan_diagnose()

sampledHSGP$log_z[[1]]

fit_sum = fit$summary()
max(fit_sum$rhat, na.rm = T)

fit_sum[3555,]
which(fit_sum$variable == "sigma_pos")


save(stan_data, sampledHSGP, fit, file = "FirstSuccessfulHSGP.RData")


fit_sum$rhat

fit$draws("drifter_zs", format = 'draws_matrix')

M=35
t_grid_res = 1000
N_draws = 4000
t_boundary = c(min(sampledParticles[,1]), max(sampledParticles[,1]))
x_mean = mean(sampledParticles[,2])
y_mean = mean(sampledParticles[,3])
TrueData = sampledParticles[,1:3]

plotPosteriorPath = function(fit, M, t_grid_res, N_draws, t_boundary, x_mean, y_mean, TrueData, drifter_ks, drifter_ls){

  path_z_draws = fit$draws("drifter_zs", format = "draws_matrix")

  x_z_draws = path_z_draws[,1:(M+1)]
  y_z_draws = path_z_draws[,1:(M+1) + M+1]

  # path_l_draws = fit$draws("drifter_ls", format = "draws_matrix")
  #
  # x_l_draws = path_l_draws[,1]
  # y_l_draws = path_l_draws[,2]
  #
  # path_k_draws = fit$draws("drifter_ks", format = "draws_matrix")
  #
  # x_k_draws = path_k_draws[,1]
  # y_k_draws = path_k_draws[,2]

  t_grid = seq(min(sampledParticles[,1]), max(sampledParticles[,1]), length.out = t_grid_res+1)

  x_pos_draws = y_pos_draws = matrix(nrow = N_draws, ncol = t_grid_res+1)

  for(i in 1:4000){

    x_pos_draws[i,] = apply(matrix(t_grid), MARGIN = 1, FUN = evaluateHSGP_1D, z = as.numeric(x_z_draws[i,]), k = drifter_ks[1],
                            l = drifter_ls[1], M=M, boundary = t_boundary, prior_mean = x_mean)

    y_pos_draws[i,] = apply(matrix(t_grid), MARGIN = 1, FUN = evaluateHSGP_1D, z = as.numeric(y_z_draws[i,]), k = drifter_ks[2],
                            l = drifter_ls[2], M=M, boundary = t_boundary, prior_mean = y_mean)

    svMisc::progress(i, 4000)

  }

  x_post_mean = colMeans(x_pos_draws)
  y_post_mean = colMeans(y_pos_draws)

  PostMean = data.frame(x = x_post_mean, y = y_post_mean)

  x_lower_q = apply(x_pos_draws, MARGIN = 2, FUN = quantile, probs = 0.025)
  y_lower_q = apply(y_pos_draws, MARGIN = 2, FUN = quantile, probs = 0.025)

  x_upper_q = apply(x_pos_draws, MARGIN = 2, FUN = quantile, probs = 0.975)
  y_upper_q = apply(y_pos_draws, MARGIN = 2, FUN = quantile, probs = 0.975)

  ConfBand = data.frame(x_lower = x_lower_q, x_upper = x_upper_q, y_lower = y_lower_q, y_upper = y_upper_q)

  ggplot() +
    # geom_rect(data = ConfBand,
    #           aes(xmin = x_lower, xmax = x_upper, ymin = y_lower, ymax = y_upper),
    #           fill = "blue", alpha = 0.2) +
    geom_path(data = PostMean, aes(x = x, y = y), size = 1.5, color = "blue") +
    geom_point(data = TrueData, aes(x = X1, y = X2), color = "red", alpha = 0.25)



}

plotPosteriorPath(fit, M, t_grid_res = 1000, N_draws, t_boundary, x_mean, y_mean, TrueData, drifter_ks = c(6.85,8.78), drifter_ls = c(0.04,0.04))

path_z_draws = fit$draws("drifter_zs", format = "draws_matrix")

x_z_draws = path_z_draws[,1:(M+1)]
y_z_draws = path_z_draws[,1:(M+1) + M+1]

path_l_draws = fit$draws("drifter_ls", format = "draws_matrix")

x_l_draws = path_l_draws[,1]
y_l_draws = path_l_draws[,2]

path_k_draws = fit$draws("drifter_ks", format = "draws_matrix")

x_k_draws = path_k_draws[,1]
y_k_draws = path_k_draws[,2]

t_grid = seq(min(sampledParticles[,1]), max(sampledParticles[,1]), length.out = t_grid_res+1)

x_pos_draws = y_pos_draws = matrix(nrow = N_draws, ncol = t_grid_res+1)

mean(y_k_draws)

sd(TrueData[,2])
sd(TrueData[,3])

median(dist(TrueData[,1])/100)

4/(0.04*pi)

for(i in 1:4000){

  x_pos_draws[i,] = apply(matrix(t_grid), MARGIN = 1, FUN = evaluateHSGP_1D, z = as.numeric(x_z_draws[i,]), k = as.numeric(x_k_draws[i,]),
                          l = as.numeric(x_l_draws[i,]), M=M, boundary = t_boundary, prior_mean = x_mean)

  y_pos_draws[i,] = apply(matrix(t_grid), MARGIN = 1, FUN = evaluateHSGP_1D, z = as.numeric(y_z_draws[i,]), k = as.numeric(y_k_draws[i,]),
                          l = as.numeric(y_l_draws[i,]), M=M, boundary = t_boundary, prior_mean = y_mean)

  svMisc::progress(i, 4000)

}

x_post_mean = colMeans(x_pos_draws)
y_post_mean = colMeans(y_pos_draws)

plot(x_post_mean, y_post_mean)

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
n_iters = 8000
M=10
GMM = dataGMM_STAN
getGridStanHSGP = function(fit, grid_res, phySpaceBorder, M, n_iters){

  log_draws <- fit$draws("log_zs", format = "draws_matrix")
  logit_draws <- fit$draws("logit_zs", format = "draws_matrix")

  log_ks_1 <- as.vector(fit$draws("log_ks[1]", format = "draws_matrix")[1:n_iters])
  log_ls_1 <- as.vector(fit$draws("log_ls[1]", format = "draws_matrix")[1:n_iters])
  log_ks_2 <- as.vector(fit$draws("log_ks[2]", format = "draws_matrix")[1:n_iters])
  log_ls_2 <- as.vector(fit$draws("log_ls[2]", format = "draws_matrix")[1:n_iters])

  logit_ks_1 <- as.vector(fit$draws("logit_ks[1]", format = "draws_matrix")[1:n_iters])
  logit_ls_1 <- as.vector(fit$draws("logit_ls[1]", format = "draws_matrix")[1:n_iters])
  logit_ks_2 <- as.vector(fit$draws("logit_ks[2]", format = "draws_matrix")[1:n_iters])
  logit_ls_2 <- as.vector(fit$draws("logit_ls[2]", format = "draws_matrix")[1:n_iters])

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

  x_seq_phy = seq(phySpaceBorder[1], phySpaceBorder[3], length.out = grid_res+1)
  y_seq_phy = seq(phySpaceBorder[2], phySpaceBorder[4], length.out = grid_res+1)

  phy_grid = expand.grid(x_seq_phy, y_seq_phy)
  phy_grid_comp = compSpaceData(GMM, data.frame(X = phy_grid[,1], Y = phy_grid[,2], Model = NA))[,c(1,2)]


  coef_grid_draws_1 = array(dim = c(grid_res+1, grid_res+1, n_iters))
  weight_grid_draws_1 = array(dim = c(grid_res+1, grid_res+1, n_iters))

  coef_grid_draws_2 = array(dim = c(grid_res+1, grid_res+1, n_iters))
  weight_grid_draws_2 = array(dim = c(grid_res+1, grid_res+1, n_iters))
i=1
curPos = as.numeric(phy_grid_comp[1,])
  for(i in 1:n_iters){
    coef_grid_draws_1[,,i] = matrix(exp(apply(phy_grid_comp, MARGIN = 1, FUN = evaluateHSGP, z = log_draws_array_1[,,i], k= log_ks_1[i], l=log_ls_1[i], M=M)), byrow = F, nrow = grid_res+1)
    weight_grid_draws_1[,,i] = matrix(invlogit(apply(phy_grid_comp, MARGIN = 1, FUN = evaluateHSGP, z = logit_draws_array_1[,,i], k= k_logit, l=l_logit, M=M)), byrow = F, nrow = grid_res+1)
    coef_grid_draws_2[,,i] = matrix(exp(apply(phy_grid_comp, MARGIN = 1, FUN = evaluateHSGP, z = log_draws_array_2[,,i], k= k_log, l=l_log, M=M)), byrow = F, nrow = grid_res+1)
    weight_grid_draws_2[,,i] = matrix(invlogit(apply(phy_grid_comp, MARGIN = 1, FUN = evaluateHSGP, z = logit_draws_array_2[,,i], k= k_logit, l=l_logit, M=M)), byrow = F, nrow = grid_res+1)
    svMisc::progress(i,4000)
  }



}


grid_res

n_iters = 4000
getGridStanHSGP = function(fit, grid_res, M, n_iters, border){

  Lx = border[3] - border[1]
  Ly = border[4] - border[2]

  # --- 1. Setup Constants ---
  # Number of coefficients per model
  N_basis = (M+1)^2
  omega = (0:M) * pi

  # Create Grid Vectors (0 to 1)
  x_seq = seq(border[1], border[3], length.out = grid_res + 1)
  y_seq = seq(border[2], border[4], length.out = grid_res + 1)

  # --- 2. Pre-calculate Basis Matrices (The Speedup) ---
  # We do this ONCE, not inside the loop.
  # Phi_x: [(grid_res+1) x (M+1)]

  # Helper to make cosine matrix
  make_phi <- function(coords, start, L) {
    # n=0 is column of 1s
    # n>0 is sqrt(2)*cos(...)
    # outer() creates the matrix of coords vs frequencies
    waves <- sqrt(2) * cos(outer((coords - start)/L, omega[2:(M+1)]))
    return(cbind(1, waves))
  }

  Phi_x <- make_phi(x_seq, border[1], Lx)
  Phi_y <- make_phi(y_seq, border[2], Ly)

  # --- 3. Extract Hyperparameters ---
  # We need the specific k and l for every draw to scale the Zs correctly
  # Assuming N_models = 2

  # Coefficient (Log) Hypers

  log_ks_1 <- rep(0.35, n_iters)
  log_ls_1 <- rep(0.2, n_iters)
  log_ks_2 <- rep(0.35, n_iters)
  log_ls_2 <- rep(0.2, n_iters)

  # --- 4. Extract Raw Z Matrices ---
  # Safer extraction using variable names

  log_draws <- fit$draws("log_zs", format = "draws_matrix")

  log_draws_1 = array(dim = c(M+1, M+1, n_iters))
  log_draws_2 = array(dim = c(M+1, M+1, n_iters))

  for(i in 1:n_iters){

    log_draws_1[,,i] = matrix(log_draws[i,1:((M+1)^2)*2 - 1], nrow = M+1, ncol = M+1, byrow = F)
    log_draws_2[,,i] = matrix(log_draws[i,1:((M+1)^2)*2], nrow = M+1, ncol = M+1, byrow = F)

  }

  # --- 5. Initialize Output Arrays ---
  traj_grid_draws_1 = array(dim = c(grid_res+1, grid_res+1, n_iters))
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
    traj_grid_draws_1[,,i] <- exp(surf_log_1) # Exponential Lin

    # --- Model 2 ---
    # A. Coefficients
    beta_log_2 <- get_beta_matrix(log_draws_2[,,i], log_ks_2[i], log_ls_2[i], M)
    surf_log_2 <- Phi_x %*% beta_log_2 %*% t(Phi_y)
    traj_grid_draws_2[,,i] <- exp(surf_log_2)

    setTxtProgressBar(pb, i)
  }
  close(pb)

  # Return list of results
  return(list(
    traj_1 = traj_grid_draws_1,
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

fast_quantile_3d = function(arr_3d, q){
  d <- dim(arr_3d)
  dim(arr_3d) <- c(d[1] * d[2], d[3])
  qs <- apply(arr_3d, 1, FUN = quantile, q)
  dim(qs) <- c(d[1], d[2])
  return(qs)
}

testSTANGrid = getGridStanHSGP(fit = fit, grid_res = 100, M = 25, n_iters = 4000, border = c(-10,-10,10,10))

truth = plotFullTrajectoriesHSGP(sampledHSGP, grid_res = 100, color_limits = c(0,4), border = c(-10,-10,10,10))
truth[[1]]$Coefficient
grid_res = 100
# Usage
x_seq = testSTANGrid$x
y_seq = testSTANGrid$y
fast_mean_3d(testSTANGrid$traj_1)
final_mean_traj1 <- fast_mean_3d(testSTANGrid$traj_1)
final_median_traj1 <- fast_quantile_3d(testSTANGrid$traj_1, 0.5)
final_lower_traj1 = fast_quantile_3d(testSTANGrid$traj_1, 0.025)
final_upper_traj1 = fast_quantile_3d(testSTANGrid$traj_1, 0.975)

c_min = 0
c_max = 4

plot_ly() %>%

  # Layer 1: The Lower Bound (Floor)
  add_surface(x = x_seq, y = y_seq, z = final_lower_traj1,
              opacity = 0.3,           # <--- Make it ghostly
              colorscale = "Viridis",
              cmin = c_min, cmax = c_max,
              showscale = FALSE) %>%   # Hide legend for this layer

  # Layer 2: The Upper Bound (Ceiling)
  add_surface(x = x_seq, y = y_seq, z = final_upper_traj1,
              opacity = 0.3,           # <--- Make it ghostly
              colorscale = "Viridis",
              cmin = c_min, cmax = c_max,
              showscale = FALSE) %>%

  # # Layer 3: The Posterior Mean (The Core)
  # add_surface(x = x_seq, y = y_seq, z = final_mean_traj1,
  #             opacity = 0.8,           # <--- Solid
  #             colorscale = "Viridis",
  #             cmin = c_min, cmax = c_max,
  #             colorbar = list(title = "Trajectory")) %>%

  # Layer 4: The Truth
  add_surface(x = x_seq, y = y_seq, z = truth[[1]]$CoefValues,
              opacity = 1.0,           # <--- Solid
              colorscale = "Viridis",
              cmin = c_min, cmax = c_max,
              colorbar = list(title = "Trajectory")) %>%

  layout(title = "Truth with 95% Credible Envelope",
         scene = list(
           zaxis = list(title = "Trajectory Value")
         ))

final_mean_traj2 <- fast_mean_3d(testSTANGrid$traj_2)
final_median_traj2 <- fast_quantile_3d(testSTANGrid$traj_2, 0.5)
final_lower_traj2 = fast_quantile_3d(testSTANGrid$traj_2, 0.025)
final_upper_traj2 = fast_quantile_3d(testSTANGrid$traj_2, 0.975)


plot_ly() %>%

  # Layer 1: The Lower Bound (Floor)
  add_surface(x = x_seq, y = y_seq, z = final_lower_traj2,
              opacity = 0.3,           # <--- Make it ghostly
              colorscale = "Viridis",
              cmin = c_min, cmax = c_max,
              showscale = FALSE) %>%   # Hide legend for this layer

  # Layer 2: The Upper Bound (Ceiling)
  add_surface(x = x_seq, y = y_seq, z = final_upper_traj2,
              opacity = 0.3,           # <--- Make it ghostly
              colorscale = "Viridis",
              cmin = c_min, cmax = c_max,
              showscale = FALSE) %>%

  # # Layer 3: The Posterior Mean (The Core)
  # add_surface(x = x_seq, y = y_seq, z = final_mean_traj2,
  #             opacity = 0.8,           # <--- Solid
  #             colorscale = "Viridis",
  #             cmin = c_min, cmax = c_max,
  #             colorbar = list(title = "Trajectory")) %>%

  # Layer 4: The Truth
  add_surface(x = x_seq, y = y_seq, z = truth[[2]]$CoefValues,
              opacity = 1.0,           # <--- Solid
              colorscale = "Viridis",
              cmin = c_min, cmax = c_max,
              colorbar = list(title = "Trajectory")) %>%

  layout(title = "Posterior Mean with 95% Credible Envelope",
         scene = list(
           zaxis = list(title = "Trajectory Value")
         ))

plot_ly() %>%

  # Layer 3: The Posterior Mean (The Core)
  add_surface(x = x_seq, y = y_seq, z = final_median_traj1,
              opacity = 0.8,           # <--- Solid
              colorscale = "Viridis",
              cmin = c_min, cmax = c_max,
              colorbar = list(title = "Trajectory")) %>%

  layout(title = "Posterior Mean with 95% Credible Envelope",
         scene = list(
           zaxis = list(title = "Trajectory Value")
         ))








plot_ly(x = x_seq, y = y_seq, z = final_mean_traj1, type = "surface", color_scale = "Viridis", cmin = 0, cmax = 2)
truth[[1]]$Physical$Trajectory

plot_ly(x = x_seq, y = y_seq, z = final_mean_traj2, type = "surface", color_scale = "Viridis", cmin = 0, cmax = 2)
truth[[2]]$Physical$Trajectory


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









