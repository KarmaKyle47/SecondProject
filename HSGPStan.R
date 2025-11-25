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

Data_t = rep(0,100)
Data_pos = matrix(nrow = 100, ncol = 2)
Data_vel = matrix(nrow = 100, ncol = 2)

for(i in 1:100){

  which_normal = sample(c(1,2), 1, prob = GMM_weights)
  Data_pos[i,] = mvrnorm(mu = GMM_means[which_normal,], Sigma = GMM_cov[,,which_normal])

  Data_vel[i,] = rowSums(baseVectorFields(0,Data_pos[i,])) + rnorm(2,0,0.01)

}



plot(Data_pos)
plot(Data_vel)

Stan_Data = matrix(c(Data_t, Data_pos, Data_vel), nrow = 100, ncol = 5, byrow = F)

stan_data = list(
  N_data = 100,                 # Number of data points
  Data = Stan_Data,                 # Particle Velocities with Positions for now (t, x, y, v_x, v_y)
  GMM_num = 2,              # Number of Gaussian Mixtures in the Transformation
  GMM_means = GMM_means,     # Means of all Gaussian Mixtures;
  GMM_cov = GMM_cov,        # Covariance matrices of all Mixtures
  GMM_weights = GMM_weights, #Weights for Mixtures
  M = 10 #Number of eigenfunctions
)

# --- 3. Compile the Stan Model ---
# This will take a minute, but it uses the REAL compiler
mod <- rstan::stan_model("STAN_Files/SurfaceTrajectoryOnlyFunctions.stan")
mod <- rstan::stan_model("STAN_Files/test_functions.stan")
mod <- rstan::stan_model("STAN_Files/SurfaceTrajectory.stan")
mod <- cmdstan_model("STAN_Files/SurfaceTrajectoryHSGP.stan")





fit <- mod$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  refresh = 10 # print update every 500 iters
)

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
