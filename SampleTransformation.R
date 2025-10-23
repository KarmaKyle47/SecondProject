### First need to sample Gaussian mixture model



sampleGMM = function(phy_space_border){

  phy_x_width = phy_space_border[3] - phy_space_border[1]
  phy_y_width = phy_space_border[4] - phy_space_border[2]

  n_comp = rpois(n = 1,lambda = 3) + 1

  mu_x = runif(n_comp, min = (9*phy_space_border[1] + phy_space_border[3])/10, max = (phy_space_border[1] + 9*phy_space_border[3])/10)
  mu_y = runif(n_comp, min = (9*phy_space_border[2] + phy_space_border[4])/10, max = (phy_space_border[2] + 9*phy_space_border[4])/10)

  Sigmas =



}
