source("AllRequiredFunctions.R")
library(mclust)

testData = sampleFromFullPrior(vel_error_params = c(1,10), pos_error_params = c(1,100000), n_particles = 100)

testData$Hyperparameters
ggplot(testData$ObsData, aes(x = X1, y = X2, color = Particle)) + geom_point()

testObsData_GMM = testObsData = testData$ObsData[,c(2,3)]
names(testObsData) = c("X","Y")
testObsData$Model = NA

testGMM = Mclust(testObsData_GMM, G = 1000, modelNames = "EEV")
mclustModelNames

testGMM

testGMM$modelName
plot(log(testGMM$parameters$pro))

testGMM_Right = list(Mean = t(testGMM$parameters$mean), Cov = apply(testGMM$parameters$variance$sigma, MARGIN = 3, FUN = identity, simplify = FALSE), Weights = testGMM$parameters$pro)

testCompData = compSpaceData(GMM = testGMM_Right, phySpace_Data = testObsData)

testCompData$Particle = testData$ObsData$Particle

ggplot(testCompData, aes(x = X, y = Y, color = Particle)) + geom_point() + theme(legend.position = "none")

