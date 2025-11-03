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

# Now you can call them just like R functions!
result <- my_add_one(10)
print(result)
# [1] 11

my_vec <- c(1, 2, 3)
vec_result <- my_vector_func(my_vec)
print(vec_result)

evaluateCubicPatchParY(coef = tree$coefs[[1]][2,], border = c(0,0,1,1), curPos = c(0.5,0.4))
evaluateCubicPatchParXY()

evaluateSampledTreeValue(tree = tree, coefs = tree$coefs[[1]], curPos = c(0.5,0.4))
evaluateSampledSurfaceValue(as.matrix(tree$boundaries), coefs = tree$coefs[[1]], curPos = c(0.5,0.4))

calculatePatch_KnownDerivates(c(0,0,1,1), c(1,1,1,1), c(1,1,1,1), c(2,2,2,2), c(-1,-1,-1,-1))

knownCorners = rnorm(49, mean = 1, sd = 0.1)
knownParX = rnorm(49, mean = 0, sd = 0.1)
knownParY = rnorm(49, mean = 0, sd = 0.1)
knownParXY = rnorm(49, mean = 0, sd = 0.1)

baseTree = generate_grid_tree(0.5, c(0,0,1,1))
baseModels = sample_models_one_pass(tree = baseTree, num_models = 4, baseWeight = 0.1)[,c(6,7)]

stan_surface = calculateSurface_KnownCorners(updatedBoundaries[,1:4], knownCorners, knownParX, knownParY, knownParXY)
r_surface = calculateSurface_KnownCorners(data.frame(L1 = updatedBoundaries[,1],
                                                     L2 = updatedBoundaries[,2],
                                                     U1 = updatedBoundaries[,3],
                                                     U2 = updatedBoundaries[,4]), knownCorners, knownParX, knownParY, knownParXY)

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
