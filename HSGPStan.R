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

Data_t = rep(0,1000)
Data_pos = matrix(nrow = 1000, ncol = 2)
Data_vel = matrix(nrow = 1000, ncol = 2)

for(i in 1:1000){

  which_normal = sample(c(1,2), 1, prob = GMM_weights)
  Data_pos[i,] = mvrnorm(mu = GMM_means[which_normal,], Sigma = GMM_cov[,,which_normal])

  Data_vel[i,] = rowSums(baseVectorFields(0,Data_pos[i,])) + rnorm(2,0,0.1)

}


plot(Data_pos)
plot(Data_vel)

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


mod <- cmdstan_model("STAN_Files/SurfaceTrajectoryHSGP.stan",
                     force_recompile = TRUE)

fit <- mod$sample(data = stan_data, chains = 1, iter_warmup = 10, iter_sampling = 10)

fit <- mod$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 5000
)



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

fit$summary()

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


# 1. Extract all draws for the coefficient Zs
# format: [iterations, chains, flat_index]
coef_draws <- fit$draws("coef_zs", format = "draws_matrix")

# 2. Calculate the column means (averaging over all MCMC samples)
coef_means_flat <- colMeans(coef_draws)

# 3. Reshape back to [M+1, M+1, N_models]
# Note: Stan flattens Column-Major (like R), so strict reshaping works.
# Dimensions: (M+1, M+1, N_models)
z_phys_array <- array(coef_means_flat, dim = c(M+1, M+1, N_models))

# Do the same for Logits
logit_draws <- fit$draws("logit_zs", format = "draws_matrix")
logit_means_flat <- colMeans(logit_draws)
z_logit_array <- array(logit_means_flat, dim = c(M+1, M+1, N_models))

plotHSGP(grid_res = 100, coefs = z_phys_array[,,1], prior_mean = 1, k = summary_stats$mean[1], l = summary_stats$mean[5], M = 10)

