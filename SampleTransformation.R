### First need to sample Gaussian mixture model

library(MCMCpack)
library(ggplot2)


sampleGMM = function(phy_space_border){

  phy_x_width = phy_space_border[3] - phy_space_border[1]
  phy_y_width = phy_space_border[4] - phy_space_border[2]

  n_comp = rpois(n = 1,lambda = 19) + 1

  mu_x = runif(n_comp, min = phy_space_border[1] + (1/10)*phy_x_width, max = phy_space_border[3] - (1/10)*phy_x_width)
  mu_y = runif(n_comp, min = phy_space_border[2] + (1/10)*phy_y_width, max = phy_space_border[4] - (1/10)*phy_y_width)

  mu = matrix(c(mu_x, mu_y), byrow = F, ncol = 2)

  Sigmas = replicate(n_comp, riwish(4, matrix(c((0.1*phy_x_width)^2, 0,0,(0.1*phy_y_width)^2), nrow = 2)), F)

  wts = as.vector(rdirichlet(n = 1, alpha = rep(1,n_comp)))


  GMM = list(Mean = mu, Cov = Sigmas, Weights = wts)

  GMM

}

sim_one_GMM_point = function(GMM){

  sel_model = sample(1:length(GMM$Weights), size = 1, prob = GMM$Weights)

  mvrnorm(mu = GMM$Mean[sel_model,], Sigma = GMM$Cov[[sel_model]])


}

simulateGMMData = function(GMM, n_points){

  sample = data.frame(t(replicate(n_points, sim_one_GMM_point(GMM))))

  names(sample) = c("X","Y")
  sample
}

GMM = sampleGMM(c(0,0,1,1))
GMM_Data = simulateGMMData(GMM, 10000)

ggplot(GMM_Data, aes(x = X, y = Y)) + geom_point()

