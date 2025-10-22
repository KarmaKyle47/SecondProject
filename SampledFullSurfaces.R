sampleFullCubicSurface_AllModels = function(border_length = 1, cells_per_dim = 10, num_models = 3, k = 1, base_weight = 0.1, trans_prop = 0.9){

  tree = generate_grid_tree(grid_size = border_length/cells_per_dim , c(0,0,border_length, border_length))

  sampledModels = sample_models_one_pass(tree = tree, num_models = num_models, baseWeight = base_weight)[,6:7]

  rawModelCoefs = replicate(n = num_models, samplePatch_FullTree(tree = tree, k = k)$coefs, simplify = F)

  updatedModelCoefs = list()
  regionTypes = c()

  for(i in 1:num_models){

    cur_sampledTree = tree
    cur_sampledTree$coefs = rawModelCoefs[[i]]
    cur_UpdatedTree = sampleTransitions(sampledTree = cur_sampledTree, sampledModels = sampledModels, model = i, trans_prop = trans_prop)
    updatedModelCoefs[[i]] = cur_UpdatedTree$coefs
    regionTypes = cbind(regionTypes, cur_UpdatedTree$boundaries$Type)

  }

  finalTree = cur_UpdatedTree

  finalTree$boundaries = finalTree$boundaries[,1:4]
  finalTree$coefs = updatedModelCoefs
  finalTree$models = data.frame(regionTypes)
  names(finalTree$models) = str_c("Model",1:num_models)

  finalTree

}

fullTree = finalTree

plotFullTree = function(fullTree, surface_grid_size = 0.01){

  num_models = length(fullTree$coefs)
  border = fullTree$border

  typePlots = list()
  surfacePlots = list()
i=3 # This still needs to be fixed
  for(i in 1:num_models){

    typePlots[[i]] = plotTransitionRegions(fullTree$boundaries, fullTree$models[,i], str_c("Model ", i, " Regions"))


    surface_grid = expand.grid(seq(border[1],border[3], surface_grid_size), seq(border[2],border[4],surface_grid_size))
    PatchValues = apply(surface_grid, 1, FUN = evaluateSampledTreeValue, tree = fullTree, coefs = fullTree$coefs[[i]])

    plottingGrid = data.frame(surface_grid, PatchValues)
    names(plottingGrid) = c('X','Y','Value')

    Value <- xtabs(Value ~ Y + X, data = plottingGrid)
    X <- as.numeric(rownames(Value))
    Y <- as.numeric(colnames(Value))

    surfacePlot = plot_ly(x = ~X, y = ~Y, z = ~Value, type = "surface")

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



  typePlots[[3]]
  surfacePlots[[10]]

  names(typePlots) = names(surfacePlots) = str_c("Model",1:num_models)

  list(ModelRegions = typePlots, ModelSurfaces = surfacePlots)

}

fullTree = sampleFullCubicSurface_AllModels(k = 0.5, num_models = 10, border_length = 1, cells_per_dim = 10, base_weight = 0.1, trans_prop = 0.99)
fullTree_Plots = plotFullTree(fullTree)

fullTree_Plots$ModelSurfaces[[2]]
