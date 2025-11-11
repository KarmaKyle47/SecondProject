sampleFullCubicSurface_AllModels = function(border_length = 1, cells_per_dim = 10, num_models = 3, k = 1, base_weight = 0.1, trans_prop = 0.9, prior_mean = 1){

  tree = generate_grid_tree(grid_size = border_length/cells_per_dim , c(0,0,border_length, border_length))

  sampledModels = sample_models_one_pass(tree = tree, num_models = num_models, baseWeight = base_weight)[,6:7]

  sampledCoefs = list()
  regionTypes = c()

  for(i in 1:num_models){

    cur_SampledTree = sampleTransitionSurface_ByGrid(tree, sampledModels, model = i, trans_prop, k, prior_mean)
    sampledCoefs[[i]] = cur_SampledTree$coefs
    regionTypes = cbind(regionTypes, cur_SampledTree$boundaries$Type)

  }

  finalTree = cur_SampledTree

  finalTree$boundaries = finalTree$boundaries[,1:4]
  finalTree$coefs = sampledCoefs
  finalTree$models = data.frame(regionTypes)
  names(finalTree$models) = str_c("Model",1:num_models)

  finalTree

}

plotFullTree = function(fullTree, surface_grid_size = 0.01){

  num_models = length(fullTree$coefs)
  border = fullTree$border

  typePlots = list()
  surfacePlots = list()

  for(i in 1:num_models){

    typePlots[[i]] = plotTransitionRegions(fullTree$boundaries, fullTree$models[,i], str_c("Model ", i, " Regions"))


    surface_grid = expand.grid(seq(border[1],border[3], surface_grid_size), seq(border[2],border[4],surface_grid_size))
    PatchValues = apply(surface_grid, 1, FUN = evaluateSampledTreeValue, tree = fullTree, coefs = fullTree$coefs[[i]])

    plottingGrid = data.frame(surface_grid, PatchValues)
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

plotFullTreeNew = function(fullTree, surface_grid_size = 0.01){

  num_models = length(fullTree$coefs)
  border = fullTree$border

  typePlots = list()
  surfacePlots = list()

  for(i in 1:num_models){

    typePlots[[i]] = visualizeNewModelExistence(Boundaries = fullTree$baseBoundaries, ModelLogits = fullTree$logits, model = i)


    surface_grid = expand.grid(seq(border[1],border[3], surface_grid_size), seq(border[2],border[4],surface_grid_size))
    PatchValues = apply(surface_grid, 1, FUN = evaluateSampledTreeValue, tree = fullTree, coefs = fullTree$coefs[[i]])

    plottingGrid = data.frame(surface_grid, PatchValues)
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

fullTree = sampleFullCubicSurface_AllModels(k = 0.1, num_models = 10, border_length = 1, cells_per_dim = 10, base_weight = 0.1, trans_prop = 8/9)
fullTree_Plots = plotFullTree(fullTree)

fullTree_Plots$ModelRegions[[5]]
fullTree_Plots$ModelSurfaces[[5]]

fullTree_Plots$ModelSurfaces[[5]]
