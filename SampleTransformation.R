### First need to sample Gaussian mixture model

library(MCMCpack)
library(ggplot2)
library(ggpubr)

sampleGMM = function(phy_space_border, lambda = 19, border_buffer_perc = 0.1, dim_sigma = 0.1){

  phy_x_width = phy_space_border[3] - phy_space_border[1]
  phy_y_width = phy_space_border[4] - phy_space_border[2]

  n_comp = rpois(n = 1,lambda = lambda) + 1

  mu_x = runif(n_comp, min = phy_space_border[1] + border_buffer_perc*phy_x_width, max = phy_space_border[3] - border_buffer_perc*phy_x_width)
  mu_y = runif(n_comp, min = phy_space_border[2] + border_buffer_perc*phy_y_width, max = phy_space_border[4] - border_buffer_perc*phy_y_width)

  mu = matrix(c(mu_x, mu_y), byrow = F, ncol = 2)

  Sigmas = replicate(n_comp, riwish(4, matrix(c((dim_sigma*phy_x_width)^2, 0,0,(dim_sigma*phy_y_width)^2), nrow = 2)), F)

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

  x_Invcdf = try(uniroot(invCDF_x_finder, interval = c(phySpaceBorder[1] - 3*(phySpaceBorder[3] - phySpaceBorder[1]), phySpaceBorder[3] + 3*(phySpaceBorder[3] - phySpaceBorder[1]))), silent = TRUE)

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

  y_given_x_Invcdf = try(uniroot(invCDF_y_given_x_finder, interval = c(phySpaceBorder[2] - 3*(phySpaceBorder[4] - phySpaceBorder[2]), phySpaceBorder[4] + 3*(phySpaceBorder[4] - phySpaceBorder[2]))), silent = TRUE)

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

sampled_GMM = testGMM_Right
tree = Traj1
phy_Data = testObsData
comp_Data = testCompData
phySpaceBorder = c(-5,-5,5,5)

visualizeSampledTransformation = function(tree, phySpaceBorder, n_data_points = 10000, boundary_grid_size = 0.01, sampled_GMM, phy_Data){

  if(n_data_points == 0){
    if(missing(sampled_GMM)){sampled_GMM = sampleGMM(phySpaceBorder)}

    x_splits = unique(c(tree$boundaries$L1, tree$boundaries$U1))
    x_splits = x_splits[2:(length(x_splits)-1)]

    y_splits = unique(c(tree$boundaries$L2, tree$boundaries$U2))
    y_splits = y_splits[2:(length(y_splits)-1)]

    x_seq = compSpaceData(sampled_GMM, data.frame(X = unique(c(seq(phySpaceBorder[1], phySpaceBorder[3], boundary_grid_size), phySpaceBorder[3])), Y = mean(phySpaceBorder[c(2,4)]), Model = NA))$X

    Y_BoundaryData_Comp = data.frame(X = rep(x_seq, length(y_splits)), Y = rep(y_splits, each = length(x_seq)), Model = rep(1:length(y_splits), each = length(x_seq)))
    Y_BoundaryData_Phy = phySpaceData(sampled_GMM, compSpace_Data = Y_BoundaryData_Comp, phySpaceBorder = phySpaceBorder)

    min_Y_Boundary = min(Y_BoundaryData_Phy$Y, na.rm = T) - boundary_grid_size
    max_Y_Boundary = max(Y_BoundaryData_Phy$Y, na.rm = T) + boundary_grid_size

    outerBorder = c(phySpaceBorder[1], min(c(phySpaceBorder[2], min_Y_Boundary)), phySpaceBorder[3], max(c(phySpaceBorder[4], max_Y_Boundary)))

    X_BoundaryData_Phy = phySpaceData(sampled_GMM, compSpace_Data = data.frame(X = x_splits, Y = y_splits[1], Model = NA), phySpaceBorder = phySpaceBorder)

    phySpacePlot = ggplot() +
      geom_line(data = Y_BoundaryData_Phy, aes(x = X, y = Y, group = Model), color = 'black', size = 0.1) +
      geom_segment(data = X_BoundaryData_Phy, aes(x = X, y = outerBorder[2], xend = X, yend = outerBorder[4]), color = 'black') +
      geom_segment(aes(x = outerBorder[1], y = outerBorder[2], xend = outerBorder[1], yend = outerBorder[4]), color = 'black') +
      geom_segment(aes(x = outerBorder[3], y = outerBorder[2], xend = outerBorder[3], yend = outerBorder[4]), color = 'black') +
      geom_segment(aes(x = outerBorder[1], y = outerBorder[2], xend = outerBorder[3], yend = outerBorder[2]), color = 'black') +
      geom_segment(aes(x = outerBorder[1], y = outerBorder[4], xend = outerBorder[3], yend = outerBorder[4]), color = 'black') +
      ggtitle("Physical Space") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

    compSpacePlot = ggplot() +
      geom_segment(data = tree$boundaries, aes(x = L1, y = L2, xend = L1, yend = U2), color = 'black') +
      geom_segment(data = tree$boundaries, aes(x = U1, y = L2, xend = U1, yend = U2), color = 'black') +
      geom_segment(data = tree$boundaries, aes(x = L1, y = L2, xend = U1, yend = L2), color = 'black') +
      geom_segment(data = tree$boundaries, aes(x = L1, y = U2, xend = U1, yend = U2), color = 'black') +
      ggtitle("Computational Space") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + xlab("X") + ylab("Y")

    ggarrange(phySpacePlot, compSpacePlot, nrow = 1)
  } else{
    if(missing(sampled_GMM)){sampled_GMM = sampleGMM(phySpaceBorder)}
    if(missing(phy_Data)){phy_Data = simulateGMMData(sampled_GMM, n_data_points)}
    names(phy_Data) = c("X","Y","Model")

    comp_Data = compSpaceData(sampled_GMM, phy_Data)

    x_splits = unique(c(tree$boundaries$L1, tree$boundaries$U1))
    x_splits = x_splits[2:(length(x_splits)-1)]

    y_splits = unique(c(tree$boundaries$L2, tree$boundaries$U2))
    y_splits = y_splits[2:(length(y_splits)-1)]

    x_seq = compSpaceData(sampled_GMM, data.frame(X = unique(c(seq(min(phy_Data$X, phySpaceBorder[1]), max(phy_Data$X, phySpaceBorder[3]), boundary_grid_size), max(phy_Data$X, phySpaceBorder[3]))), Y = mean(phySpaceBorder[c(2,4)]), Model = NA))$X

    Y_BoundaryData_Comp = data.frame(X = rep(x_seq, length(y_splits)), Y = rep(y_splits, each = length(x_seq)), Model = rep(1:length(y_splits), each = length(x_seq)))
    Y_BoundaryData_Phy = phySpaceData(sampled_GMM, compSpace_Data = Y_BoundaryData_Comp, phySpaceBorder = phySpaceBorder)

    min_Y_Boundary = min(Y_BoundaryData_Phy$Y, na.rm = T) - boundary_grid_size
    max_Y_Boundary = max(Y_BoundaryData_Phy$Y, na.rm = T) + boundary_grid_size

    outerBorder = c(min(phy_Data$X, phySpaceBorder[1]), min(phy_Data$Y, phySpaceBorder[2], min_Y_Boundary), max(phy_Data$X, phySpaceBorder[3]), max(phy_Data$Y, phySpaceBorder[4], max_Y_Boundary))

    X_BoundaryData_Phy = phySpaceData(sampled_GMM, compSpace_Data = data.frame(X = x_splits, Y = y_splits[1], Model = NA), phySpaceBorder = phySpaceBorder)

    phySpacePlot = ggplot() + geom_point(data = phy_Data, aes(x = X, y = Y, color = as.factor(Model))) +
      geom_line(data = Y_BoundaryData_Phy, aes(x = X, y = Y, group = Model), color = 'black', size = 0.1) +
      geom_segment(data = X_BoundaryData_Phy, aes(x = X, y = outerBorder[2], xend = X, yend = outerBorder[4]), color = 'black') +
      geom_segment(aes(x = outerBorder[1], y = outerBorder[2], xend = outerBorder[1], yend = outerBorder[4]), color = 'black') +
      geom_segment(aes(x = outerBorder[3], y = outerBorder[2], xend = outerBorder[3], yend = outerBorder[4]), color = 'black') +
      geom_segment(aes(x = outerBorder[1], y = outerBorder[2], xend = outerBorder[3], yend = outerBorder[2]), color = 'black') +
      geom_segment(aes(x = outerBorder[1], y = outerBorder[4], xend = outerBorder[3], yend = outerBorder[4]), color = 'black') +
      ggtitle("Physical Space") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

    compSpacePlot = ggplot() + geom_point(data = comp_Data, aes(x = X, y = Y, color = as.factor(Model))) +
      geom_segment(data = tree$boundaries, aes(x = L1, y = L2, xend = L1, yend = U2), color = 'black') +
      geom_segment(data = tree$boundaries, aes(x = U1, y = L2, xend = U1, yend = U2), color = 'black') +
      geom_segment(data = tree$boundaries, aes(x = L1, y = L2, xend = U1, yend = L2), color = 'black') +
      geom_segment(data = tree$boundaries, aes(x = L1, y = U2, xend = U1, yend = U2), color = 'black') +
      ggtitle("Computational Space") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

    ggarrange(phySpacePlot, compSpacePlot, nrow = 1)
  }



}

plotTransformedRegions = function(GMM, treeBoundaries, types, phySpaceBorder, boundary_grid_size = 0.001, title = "") {

  x_splits = sort(unique(c(treeBoundaries$L1, treeBoundaries$U1)))
  x_splits = x_splits[2:(length(x_splits)-1)]

  y_splits = sort(unique(c(treeBoundaries$L2, treeBoundaries$U2)))
  y_splits = y_splits[2:(length(y_splits)-1)]

  x_seq = compSpaceData(GMM, data.frame(X = unique(c(seq(phySpaceBorder[1], phySpaceBorder[3], boundary_grid_size), phySpaceBorder[3])), Y = mean(phySpaceBorder[c(2,4)]), Model = NA))$X

  Y_BoundaryData_Comp = data.frame(X = rep(x_seq, length(y_splits)), Y = rep(y_splits, each = length(x_seq)), Model = rep(1:length(y_splits), each = length(x_seq)))
  Y_BoundaryData_Phy = phySpaceData(GMM, compSpace_Data = Y_BoundaryData_Comp, phySpaceBorder = phySpaceBorder)

  min_Y_Boundary = min(Y_BoundaryData_Phy$Y, na.rm = T) - boundary_grid_size
  max_Y_Boundary = max(Y_BoundaryData_Phy$Y, na.rm = T) + boundary_grid_size

  outerBorder = c(phySpaceBorder[1], min(c(phySpaceBorder[2], min_Y_Boundary)), phySpaceBorder[3], max(c(phySpaceBorder[4], max_Y_Boundary)))

  X_BoundaryData_Phy = phySpaceData(GMM, compSpace_Data = data.frame(X = x_splits, Y = y_splits[1], Model = NA), phySpaceBorder = phySpaceBorder)

  RegionRect_List = list()

  x_grid = sort(unique(Y_BoundaryData_Phy$X))

  for(i in 0:max(Y_BoundaryData_Comp$Model)){

    for(j in 1:(length(x_grid)-1)){

      if(i == 0){

        cur_CompPos = compSpaceData(GMM, data.frame(X = mean(c(x_grid[j], x_grid[j+1])), Y = mean(c(outerBorder[2], Y_BoundaryData_Phy$Y[Y_BoundaryData_Phy$X == x_grid[j] & Y_BoundaryData_Phy$Model == i+1])), Model = NA))[c(1,2)]
        cur_type = types[treeBoundaries$L1 <= cur_CompPos[1,1] & treeBoundaries$U1 >= cur_CompPos[1,1] & treeBoundaries$L2 <= cur_CompPos[1,2] & treeBoundaries$U2 >= cur_CompPos[1,2]]

        RegionRect_List[[length(RegionRect_List)+1]] = data.frame(xmin = x_grid[j], xmax = x_grid[j+1], ymin = outerBorder[2], ymax = Y_BoundaryData_Phy$Y[Y_BoundaryData_Phy$X == x_grid[j] & Y_BoundaryData_Phy$Model == i+1], type = cur_type)

      } else if(i == max(Y_BoundaryData_Comp$Model)){

        cur_CompPos = compSpaceData(GMM, data.frame(X = mean(c(x_grid[j], x_grid[j+1])), Y = mean(c(outerBorder[4], Y_BoundaryData_Phy$Y[Y_BoundaryData_Phy$X == x_grid[j] & Y_BoundaryData_Phy$Model == i])), Model = NA))[c(1,2)]
        cur_type = types[treeBoundaries$L1 <= cur_CompPos[1,1] & treeBoundaries$U1 >= cur_CompPos[1,1] & treeBoundaries$L2 <= cur_CompPos[1,2] & treeBoundaries$U2 >= cur_CompPos[1,2]]

        RegionRect_List[[length(RegionRect_List)+1]] = data.frame(xmin = x_grid[j], xmax = x_grid[j+1], ymin = Y_BoundaryData_Phy$Y[Y_BoundaryData_Phy$X == x_grid[j] & Y_BoundaryData_Phy$Model == i], ymax = outerBorder[4], type = cur_type)

      } else{

        cur_CompPos = compSpaceData(GMM, data.frame(X = mean(c(x_grid[j], x_grid[j+1])), Y = mean(c(Y_BoundaryData_Phy$Y[Y_BoundaryData_Phy$X == x_grid[j] & Y_BoundaryData_Phy$Model == i], Y_BoundaryData_Phy$Y[Y_BoundaryData_Phy$X == x_grid[j] & Y_BoundaryData_Phy$Model == (i+1)])), Model = NA))[c(1,2)]
        cur_type = types[treeBoundaries$L1 <= cur_CompPos[1,1] & treeBoundaries$U1 >= cur_CompPos[1,1] & treeBoundaries$L2 <= cur_CompPos[1,2] & treeBoundaries$U2 >= cur_CompPos[1,2]]

        RegionRect_List[[length(RegionRect_List)+1]] = data.frame(xmin = x_grid[j], xmax = x_grid[j+1], ymin = Y_BoundaryData_Phy$Y[Y_BoundaryData_Phy$X == x_grid[j] & Y_BoundaryData_Phy$Model == i], ymax = Y_BoundaryData_Phy$Y[Y_BoundaryData_Phy$X == x_grid[j] & Y_BoundaryData_Phy$Model == (i+1)], type = cur_type)

      }

    }


  }

  RegionRect_df = do.call(rbind, RegionRect_List)


  phySpacePlot = ggplot() +
    geom_rect(data = RegionRect_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = type), alpha = 0.4) +
    geom_line(data = Y_BoundaryData_Phy, aes(x = X, y = Y, group = Model), color = 'black', size = 0.1) +
    geom_segment(data = X_BoundaryData_Phy, aes(x = X, y = outerBorder[2], xend = X, yend = outerBorder[4]), color = 'black') +
    geom_segment(aes(x = outerBorder[1], y = outerBorder[2], xend = outerBorder[1], yend = outerBorder[4]), color = 'black') +
    geom_segment(aes(x = outerBorder[3], y = outerBorder[2], xend = outerBorder[3], yend = outerBorder[4]), color = 'black') +
    geom_segment(aes(x = outerBorder[1], y = outerBorder[2], xend = outerBorder[3], yend = outerBorder[2]), color = 'black') +
    geom_segment(aes(x = outerBorder[1], y = outerBorder[4], xend = outerBorder[3], yend = outerBorder[4]), color = 'black') +
    ggtitle(title) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

  phySpacePlot
}

plotTransformedFullTree = function(fullTree, GMM, phySpaceBorder, surface_grid_size = 0.01, region_grid_size = 0.001){

  num_models = length(fullTree$coefs)
  border = fullTree$border

  typePlots = list()
  surfacePlots = list()

  x_phy_seq = unique(c(seq(phySpaceBorder[1], phySpaceBorder[3], surface_grid_size), phySpaceBorder[3]))
  y_phy_seq = unique(c(seq(phySpaceBorder[2], phySpaceBorder[4], surface_grid_size), phySpaceBorder[4]))

  phy_grid = expand.grid(x_phy_seq, y_phy_seq)

  comp_grid = compSpaceData(GMM, data.frame(X = phy_grid[,1], Y = phy_grid[,2], Model = NA))[,c(1,2)]

  for(i in 1:num_models){

    typePlots[[i]] = plotTransformedRegions(GMM, fullTree$boundaries, fullTree$models[,i], phySpaceBorder, region_grid_size, str_c("Model ", i, " Regions"))

    PatchValues = apply(comp_grid, 1, FUN = evaluateSampledTreeValue, tree = fullTree, coefs = fullTree$coefs[[i]])

    plottingGrid = data.frame(phy_grid, PatchValues)
    names(plottingGrid) = c('X','Y','Value')

    Value <- xtabs(Value ~ Y + X, data = plottingGrid)
    X <- as.numeric(rownames(Value))
    Y <- as.numeric(colnames(Value))

    surfacePlot = plot_ly(x = X, y = Y, z = Value, type = "surface")

    surfacePlot = surfacePlot %>% layout(
      title = str_c("Model ", i, " Surface"),
      scene = list(
        xaxis = list(title = 'X'),
        yaxis = list(title = 'Y'),
        zaxis = list(
          title = 'Value'
        )
      )
    )

    surfacePlots[[i]] = surfacePlot

    svMisc::progress(i, num_models)

  }


  names(typePlots) = names(surfacePlots) = str_c("Model",1:num_models)

  list(ModelRegions = typePlots, ModelSurfaces = surfacePlots)

}



tree = generate_grid_tree(0.1, border = c(0,0,1,1))

visualizeSampledTransformation(tree, phySpaceBorder, n_data_points = 1000, boundary_grid_size = 0.01)

phySpaceBorder = c(0,0,1,1)

GMM = sampleGMM(phySpaceBorder)
GMM_Data = simulateGMMData(GMM, 10000)
GMM_Data_Comp = compSpaceData(GMM, GMM_Data)
GMM_Data_Phy = phySpaceData(GMM, compSpace_Data = GMM_Data_Comp, phySpaceBorder)


ggplot(GMM_Data, aes(x = X, y = Y, color = as.factor(Model))) + geom_point()
ggplot(GMM_Data_Comp, aes(x = X, y = Y, color = as.factor(Model))) + geom_point()
ggplot(GMM_Data_Phy, aes(x = X, y = Y, color = as.factor(Model))) + geom_point()

plotCompPatch = plotFullTree(fullTree)
plotTransPatch = plotTransformedFullTree(fullTree = fullTree, GMM = GMM, phySpaceBorder, surface_grid_size = 0.01, region_grid_size = 0.001)

plotCompPatch$ModelRegions[[10]]

plotTransPatch$ModelRegions[[10]]


plotTransformedRegions(GMM, treeBoundaries = fullTree$boundaries, types = fullTree$models[,2], phySpaceBorder, boundary_grid_size = 0.005)



