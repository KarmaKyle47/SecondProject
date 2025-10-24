### First need to sample Gaussian mixture model

library(MCMCpack)
library(ggplot2)
library(ggpubr)

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

  c(mvrnorm(mu = GMM$Mean[sel_model,], Sigma = GMM$Cov[[sel_model]]), sel_model)


}

simulateGMMData = function(GMM, n_points){

  sample = data.frame(t(replicate(n_points, sim_one_GMM_point(GMM))))

  names(sample) = c("X","Y","Model")
  sample
}

get_compSpace_pos = function(GMM, curPos_phy){

  n_models = length(GMM$Weights)

  x_cdf = sum(apply(matrix(1:n_models),MARGIN = 1, FUN = function(model){
    pnorm(curPos_phy[1], mean = GMM$Mean[model, 1], sd = sqrt(GMM$Cov[[model]][1,1]))
  }) * GMM$Weights)

  cond_log_weights_uw = apply(matrix(1:n_models),MARGIN = 1, FUN = function(model){
    dnorm(curPos_phy[1], mean = GMM$Mean[model, 1], sd = sqrt(GMM$Cov[[model]][1,1]), log = T) + log(GMM$Weights[model])
  })

  max_cond_log_weights_uw = max(cond_log_weights_uw)

  cond_log_weights = cond_log_weights_uw - (max_cond_log_weights_uw + log(sum(exp(cond_log_weights_uw - max_cond_log_weights_uw))))

  cond_weights = exp(cond_log_weights)

  y_given_x_cdf = sum(apply(matrix(1:n_models), MARGIN = 1, FUN = function(model){
    pnorm(curPos_phy[2], mean = GMM$Mean[model, 2] + GMM$Cov[[model]][2,1]*(1/GMM$Cov[[model]][1,1])*(curPos_phy[1] - GMM$Mean[model, 1]), sd = sqrt(GMM$Cov[[model]][2,2] - GMM$Cov[[model]][2,1]*(1/GMM$Cov[[model]][1,1])*GMM$Cov[[model]][1,2]))
  }) * cond_weights)

  c(x_cdf, y_given_x_cdf)

}

get_phySpace_pos = function(GMM, curPos_comp, phySpaceBorder){

  n_models = length(GMM$Weights)

  invCDF_x_finder = function(x){

    sum(apply(matrix(1:n_models),MARGIN = 1, FUN = function(model){
      pnorm(x, mean = GMM$Mean[model, 1], sd = sqrt(GMM$Cov[[model]][1,1]))
    }) * GMM$Weights) -  curPos_comp[1]

  }

  x_Invcdf = try(uniroot(invCDF_x_finder, interval = c(phySpaceBorder[1] - (phySpaceBorder[3] - phySpaceBorder[1]), phySpaceBorder[3] + (phySpaceBorder[3] - phySpaceBorder[1]))), silent = TRUE)

  if(inherits(x_Invcdf, "try-error")) {
    warning("Could not find inverse root for x.")
    return(c(x = NA, y = NA))
  }
  x_Invcdf = x_Invcdf$root

  invCDF_y_given_x_finder = function(y){

    cond_log_weights_uw = apply(matrix(1:n_models),MARGIN = 1, FUN = function(model){
      dnorm(x_Invcdf, mean = GMM$Mean[model, 1], sd = sqrt(GMM$Cov[[model]][1,1]), log = T) + log(GMM$Weights[model])
    })

    max_cond_log_weights_uw = max(cond_log_weights_uw)

    cond_log_weights = cond_log_weights_uw - (max_cond_log_weights_uw + log(sum(exp(cond_log_weights_uw - max_cond_log_weights_uw))))

    cond_weights = exp(cond_log_weights)

    y_given_x_cdf = sum(apply(matrix(1:n_models), MARGIN = 1, FUN = function(model){
      pnorm(y, mean = GMM$Mean[model, 2] + GMM$Cov[[model]][2,1]*(1/GMM$Cov[[model]][1,1])*(x_Invcdf - GMM$Mean[model, 1]), sd = sqrt(GMM$Cov[[model]][2,2] - GMM$Cov[[model]][2,1]*(1/GMM$Cov[[model]][1,1])*GMM$Cov[[model]][1,2]))
    }) * cond_weights) - curPos_comp[2]

  }

  y_given_x_Invcdf = try(uniroot(invCDF_y_given_x_finder, interval = c(phySpaceBorder[2] - (phySpaceBorder[4] - phySpaceBorder[2]), phySpaceBorder[4] + (phySpaceBorder[4] - phySpaceBorder[2]))), silent = TRUE)

  if(inherits(y_given_x_Invcdf, "try-error")) {
    warning("Could not find inverse root for y.")
    return(c(x = NA, y = NA))
  }
  y_given_x_Invcdf = y_given_x_Invcdf$root

  c(x_Invcdf, y_given_x_Invcdf)

}

compSpaceData = function(GMM, phySpace_Data){

  compSpace = data.frame(t(apply(matrix(1:nrow(phySpace_Data)), MARGIN = 1, FUN = function(i){get_compSpace_pos(GMM, c(phySpace_Data$X[i], phySpace_Data$Y[i]))})), phySpace_Data$Model)
  names(compSpace) = c("X","Y","Model")

  compSpace

}

phySpaceData = function(GMM, compSpace_Data, phySpaceBorder){

  phySpace = data.frame(t(apply(matrix(1:nrow(compSpace_Data)), MARGIN = 1, FUN = function(i){
    svMisc::progress(i, nrow(compSpace_Data))
    get_phySpace_pos(GMM, curPos_comp = c(compSpace_Data$X[i], compSpace_Data$Y[i]), phySpaceBorder)})), compSpace_Data$Model)
  names(phySpace) = c("X","Y","Model")

  phySpace

}

visualizeSampledTransformation = function(tree, phySpaceBorder, n_data_points = 10000, boundary_grid_size = 0.01){

  sampled_GMM = sampleGMM(phySpaceBorder)
  phy_Data = simulateGMMData(sampled_GMM, n_data_points)
  comp_Data = compSpaceData(sampled_GMM, phy_Data)

  boundary_data_list = list()

  for(i in 1:nrow(tree$boundaries)){

    L1 = tree$boundaries[i,]$L1
    if(L1 == 0){
      L1 = 0.0001
    }
    L2 = tree$boundaries[i,]$L2
    if(L2 == 0){
      L2 = 0.0001
    }
    U1 = tree$boundaries[i,]$U1
    if(U1 == 1){
      U1 = 0.9999
    }
    U2 = tree$boundaries[i,]$U2
    if(U2 == 1){
      U2 = 0.9999
    }

    X_data = seq(L1, U1, boundary_grid_size)
    if(!(U1 %in% X_data)){
      X_data = c(X_data, U1)
    }

    line_bottom = data.frame(X = X_data, Y = L2, line_id = str_c("box", i, "_b"))
    line_top = data.frame(X = X_data, Y = U2, line_id = str_c("box", i, "_t"))
    line_left = data.frame(X = L1, Y = c(L2,U2), line_id = str_c("box", i, "_l"))
    line_right = data.frame(X = U1, Y = c(L2,U2), line_id = str_c("box", i, "_r"))

    boundary_data_list[[i]] = rbind(line_bottom, line_top, line_left, line_right)

  }

  BoundaryData_Comp = data.frame(do.call(rbind, boundary_data_list))
  names(BoundaryData_Comp) = c("X","Y","Model")

  BoundaryData_Phy = phySpaceData(sampled_GMM, compSpace_Data = BoundaryData_Comp, phySpaceBorder = phySpaceBorder)

  phySpacePlot = ggplot() + geom_point(data = phy_Data, aes(x = X, y = Y, color = as.factor(Model))) +
    geom_line(data = BoundaryData_Phy, aes(x = X, y = Y, group = Model), color = 'black', size = 0.1) +
    ggtitle("Physical Space") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

  compSpacePlot = ggplot() + geom_point(data = comp_Data, aes(x = X, y = Y, color = as.factor(Model))) +
    geom_segment(data = tree$boundaries, aes(x = L1, y = L2, xend = L1, yend = U2), color = 'black') +
    geom_segment(data = tree$boundaries, aes(x = U1, y = L2, xend = U1, yend = U2), color = 'black') +
    geom_segment(data = tree$boundaries, aes(x = L1, y = L2, xend = U1, yend = L2), color = 'black') +
    geom_segment(data = tree$boundaries, aes(x = L1, y = U2, xend = U1, yend = U2), color = 'black') +
    ggtitle("Computational Space") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

  ggarrange(phySpacePlot, compSpacePlot, nrow = 1)

}

tree = generate_grid_tree(0.25, border = c(0,0,1,1))

visualizeSampledTransformation(tree, phySpaceBorder, n_data_points = 1000, boundary_grid_size = 0.001)

phySpaceBorder = c(0,0,10,10)

GMM = sampleGMM(phySpaceBorder)
GMM_Data = simulateGMMData(GMM, 10000)
GMM_Data_Comp = compSpaceData(GMM, GMM_Data)
GMM_Data_Phy = phySpaceData(GMM, compSpace_Data = GMM_Data_Comp, phySpaceBorder)


ggplot(GMM_Data, aes(x = X, y = Y, color = as.factor(Model))) + geom_point()
ggplot(GMM_Data_Comp, aes(x = X, y = Y, color = as.factor(Model))) + geom_point()
ggplot(GMM_Data_Phy, aes(x = X, y = Y, color = as.factor(Model))) + geom_point()

