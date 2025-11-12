source("AllRequiredFunctions.R")
library(mclust)

testData = sampleFromFullPrior(vel_error_params = c(1,10), pos_error_params = c(1,100000), n_particles = 100)
testData = testPrior
testData$Hyperparameters
ggplot(testData$ObsData, aes(x = X1, y = X2, color = Particle)) + geom_point()

testObsData_GMM = testObsData = testData$ObsData[,c(2,3)]
names(testObsData) = c("X","Y")
testObsData$Model = NA

testGMM = Mclust(testObsData_GMM, G = 100)
mclustModelNames

testGMM

testGMM$modelName
plot(log(testGMM$parameters$pro))

testGMM_Right = list(Mean = t(testGMM$parameters$mean), Cov = apply(testGMM$parameters$variance$sigma, MARGIN = 3, FUN = identity, simplify = FALSE), Weights = testGMM$parameters$pro)

testCompData = compSpaceData(GMM = testGMM_Right, phySpace_Data = testObsData)

testCompData$Particle = testData$ObsData$Particle

ggplot(testCompData, aes(x = X, y = Y, color = Particle)) + geom_point() + theme(legend.position = "none")


library(rstan)

# This compiles the functions and loads them into your R session
expose_stan_functions("STAN_Files/SurfaceTrajectoryOnlyFunctions.stan")

evaluateCubicPatchParY(coef = tree$coefs[[1]][2,], border = c(0,0,1,1), curPos = c(0.5,0.4))
evaluateCubicPatchParXY()

evaluateSampledTreeValue(tree = tree, coefs = tree$coefs[[1]], curPos = c(0.5,0.4))
evaluateSampledSurfaceValue(as.matrix(tree$boundaries), coefs = tree$coefs[[1]], curPos = c(0.5,0.4))

calculatePatch_KnownDerivates(c(0,0,1,1), c(1,1,1,1), c(1,1,1,1), c(2,2,2,2), c(-1,-1,-1,-1))

knownCorners = rnorm(49, mean = 1, sd = 0.1)
knownParX = rnorm(49, mean = 0, sd = 0.1)
knownParY = rnorm(49, mean = 0, sd = 0.1)
knownParXY = rnorm(49, mean = 0, sd = 0.1)

t1 = Sys.time()
baseTree = generate_grid_tree(0.5, c(0,0,1,1))
t2 = Sys.time()
baseModels = sample_models_one_pass(tree = baseTree, num_models = 4, baseWeight = 0.1)[,c(6,7)]

t1 = Sys.time()
stan_surface = calculateSurface_KnownCorners(updatedBoundaries[,1:4], knownCorners, knownParX, knownParY, knownParXY)
t2 = Sys.time()
t2-t1
r_surface = calculateSurface_KnownCorners(data.frame(L1 = updatedBoundaries[,1],
                                                     L2 = updatedBoundaries[,2],
                                                     U1 = updatedBoundaries[,3],
                                                     U2 = updatedBoundaries[,4]), knownCorners, knownParX, knownParY, knownParXY)
stan_updated_surface = updateCornerQuantities()

t3 = Sys.time()

stan_surface == r_surface
updatedBoundaries = updateCompGridTransitions(8/9, baseCompGridBoundaries = as.matrix(baseTree$boundaries), comp_res = 2, models = as.matrix(baseModels), numModels = 4)

ShadeData = data.frame(xmin = updatedBoundaries[,1],
                       xmax = updatedBoundaries[,3],
                       ymin = updatedBoundaries[,2],
                       ymax = updatedBoundaries[,4],
                       Region = str_c("Region", 1:nrow(updatedBoundaries)))

plot = ggplot() +
  geom_segment(data = ShadeData, aes(x = xmin, y = ymin, xend = xmin, yend = ymax), color = 'black') +
  geom_segment(data = ShadeData, aes(x = xmax, y = ymin, xend = xmax, yend = ymax), color = 'black') +
  geom_segment(data = ShadeData, aes(x = xmin, y = ymin, xend = xmax, yend = ymin), color = 'black') +
  geom_segment(data = ShadeData, aes(x = xmin, y = ymax, xend = xmax, yend = ymax), color = 'black') +
  geom_rect(data = ShadeData,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                fill = Region),
            inherit.aes = FALSE,
            alpha = 0.4) +
  xlab("Longitude") + ylab("Latitude")

plot

newTree = generate_grid_tree(0.25,c(0,0,1,1))
baseTree = generate_grid_tree(0.25,c(0,0,1,1))
sampleNewModels = sample_models_one_pass(newTree, num_models = 3, baseWeight = 0.1)[,c(6,7)]

updatedNewBoundaries = updateCompGridTransitions(8/9, as.matrix(newTree$boundaries),4, as.matrix(sampleNewModels), 3)

knownCornersBase = rnorm(25, mean = 1, sd = 0.1)
knownParXBase = rnorm(25, mean = 0, sd = 0.1)
knownParYBase = rnorm(25, mean = 0, sd = 0.1)
knownParXYBase = rnorm(25, mean = 0, sd = 0.1)

updateCornerQuantities(knownCornersBase, knownParXBase, knownParYBase, knownParXYBase, baseBoundaries = as.matrix(newTree$boundaries), updatedBoundaries = updatedNewBoundaries, model_num = 1)

newTree$models = updatedNewBoundaries[,6:8]
newTree$coefs = list(calculateSurface_KnownCorners(newTree$boundaries, knownCornersBase, knownParXBase, knownParYBase, knownParXYBase),
                     calculateSurface_KnownCorners(newTree$boundaries, knownCornersBase, knownParXBase, knownParYBase, knownParXYBase),
                     calculateSurface_KnownCorners(newTree$boundaries, knownCornersBase, knownParXBase, knownParYBase, knownParXYBase))

updatedCornerQuantities1 = updateCornerQuantities(knownCornersBase, knownParXBase, knownParYBase, knownParXYBase, baseBoundaries = as.matrix(newTree$boundaries), updatedBoundaries = updatedNewBoundaries, model_num = 1)
updatedCornerQuantities2 = updateCornerQuantities(knownCornersBase, knownParXBase, knownParYBase, knownParXYBase, baseBoundaries = as.matrix(newTree$boundaries), updatedBoundaries = updatedNewBoundaries, model_num = 2)
updatedCornerQuantities3 = updateCornerQuantities(knownCornersBase, knownParXBase, knownParYBase, knownParXYBase, baseBoundaries = as.matrix(newTree$boundaries), updatedBoundaries = updatedNewBoundaries, model_num = 3)

updatedCoefs1 = calculateSurface_KnownCorners(data.frame(L1 = updatedNewBoundaries[,1],
                                                         L2 = updatedNewBoundaries[,2],
                                                         U1 = updatedNewBoundaries[,3],
                                                         U2 = updatedNewBoundaries[,4]), GridValues = updatedCornerQuantities1[,1], updatedCornerQuantities1[,2], updatedCornerQuantities1[,3], updatedCornerQuantities1[,4])

updatedCoefs2 = calculateSurface_KnownCorners(data.frame(L1 = updatedNewBoundaries[,1],
                                                         L2 = updatedNewBoundaries[,2],
                                                         U1 = updatedNewBoundaries[,3],
                                                         U2 = updatedNewBoundaries[,4]), GridValues = updatedCornerQuantities2[,1], updatedCornerQuantities2[,2], updatedCornerQuantities2[,3], updatedCornerQuantities2[,4])

updatedCoefs3 = calculateSurface_KnownCorners(data.frame(L1 = updatedNewBoundaries[,1],
                                                         L2 = updatedNewBoundaries[,2],
                                                         U1 = updatedNewBoundaries[,3],
                                                         U2 = updatedNewBoundaries[,4]), GridValues = updatedCornerQuantities3[,1], updatedCornerQuantities3[,2], updatedCornerQuantities3[,3], updatedCornerQuantities3[,4])

updatedTree = list(boundaries = data.frame(L1 = updatedNewBoundaries[,1],
                                           L2 = updatedNewBoundaries[,2],
                                           U1 = updatedNewBoundaries[,3],
                                           U2 = updatedNewBoundaries[,4]),
                   coefs = list(updatedCoefs1, updatedCoefs2, updatedCoefs3),
                   models = updatedNewBoundaries[,6:8], border = c(0,0,1,1))

plotFullTree(newTree)
plotUpdatedTree = plotFullTree(updatedTree)

plotUpdatedTree$ModelSurfaces[[2]]



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

GMM_means = matrix(c(0,0,0,0), nrow = 2)
GMM_cov = array(c(1,0,0,1,1,0,0,1), dim = c(2,2,2))
GMM_weights = c(1,0)

Data_t = rep(0,100)
Data_pos = matrix(rnorm(200), nrow = 100)
Data_vel = matrix(nrow = 100, ncol = 2)

for(i in 1:100){

  Data_vel[i,] = rowSums(baseVectorFields(0, Data_pos[i,])) + rnorm(2,0,0.1)

}

Stan_Data = matrix(c(Data_t, Data_pos, Data_vel), nrow = 100, ncol = 5, byrow = F)

stan_data = list(
  N_data = 100,                 # Number of data points
  comp_res = 1,                 # Computational Grid Resolution
  Data = Stan_Data,                 # Particle Velocities with Positions for now (t, x, y, v_x, v_y)
  GMM_num = 2,              # Number of Gaussian Mixtures in the Transformation
  GMM_means = GMM_means,     # Means of all Gaussian Mixtures;
  GMM_cov = GMM_cov,        # Covariance matrices of all Mixtures
  GMM_weights = GMM_weights #Weights for Mixtures
)

# --- 3. Compile the Stan Model ---
# This will take a minute, but it uses the REAL compiler
mod <- rstan::stan_model("STAN_Files/SurfaceTrajectoryOnlyFunctions.stan")
mod <- rstan::stan_model("STAN_Files/test_functions.stan")
mod <- rstan::stan_model("STAN_Files/SurfaceTrajectory.stan")


# --- 4. Run the Model to Test the Function ---
# This part is very fast
fit <- rstan::sampling(
  object = mod,
  data = stan_data,
  algorithm = "Fixed_param",  # <-- The magic trick
  iter = 1,                   # We only need one "draw"
  warmup = 0,
  chains = 1
)

# --- 5. Extract and Inspect the Result ---
fit_extract <- rstan::extract(fit)

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
