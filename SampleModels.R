find_neighbors <- function(all_boundaries, index) {

  # Coordinates of the current boundary
  b <- all_boundaries[index, ]

  # Find boundaries that share an edge
  # Note: The small tolerance `1e-9` handles potential floating-point inaccuracies.
  right_neighbors <- which(abs(all_boundaries$L1 - b$U1) < 1e-9 & all_boundaries$L2 < b$U2 & all_boundaries$U2 > b$L2)
  left_neighbors  <- which(abs(all_boundaries$U1 - b$L1) < 1e-9 & all_boundaries$L2 < b$U2 & all_boundaries$U2 > b$L2)
  up_neighbors    <- which(abs(all_boundaries$L2 - b$U2) < 1e-9 & all_boundaries$L1 < b$U1 & all_boundaries$U1 > b$L1)
  down_neighbors  <- which(abs(all_boundaries$U2 - b$L2) < 1e-9 & all_boundaries$L1 < b$U1 & all_boundaries$U1 > b$L1)

  return(unique(c(right_neighbors, left_neighbors, up_neighbors, down_neighbors)))
}

sample_models_per_region = function(tree, num_models, baseWeight = 1){

  treeBoundaries = data.frame(orderBoundaries_GeminiCleaned(tree)[[1]], model1 = 0, model2 = 0)

  baseTransMat = (baseWeight/num_models)*(matrix(rep(1, num_models^2), nrow = num_models) - diag(num_models))

  for(i in 1:nrow(treeBoundaries)){

    curIndex = which(treeBoundaries$order == i)

    neighbors = treeBoundaries[find_neighbors(treeBoundaries, curIndex),]
    definedNeighbors = neighbors[neighbors$order < i,]

    curTransMat = baseTransMat + table(factor(definedNeighbors$model1, levels = 1:num_models), factor(definedNeighbors$model2, levels = 1:num_models)) +
      table(factor(definedNeighbors$model2, levels = 1:num_models), factor(definedNeighbors$model1, levels = 1:num_models))

    model1 = sample.int(num_models, size = 1, prob = rowSums(curTransMat))
    model2 = sample.int(num_models, size = 1, prob = curTransMat[model1,])

    treeBoundaries$model1[curIndex] = model1
    treeBoundaries$model2[curIndex] = model2


  }

  treeBoundaries

}

visualizeModelExistence = function(sampledTree, model){

  ShadeData = data.frame(xmin = sampledTree$L1,
                         xmax = sampledTree$U1,
                         ymin = sampledTree$L2,
                         ymax = sampledTree$U2,
                         Region = str_c("Region", 1:nrow(sampledTree)))

  ggplot() +
    geom_segment(data = sampledTree, aes(x = L1, y = L2, xend = L1, yend = U2), color = 'black') +
    geom_segment(data = sampledTree, aes(x = U1, y = L2, xend = U1, yend = U2), color = 'black') +
    geom_segment(data = sampledTree, aes(x = L1, y = L2, xend = U1, yend = L2), color = 'black') +
    geom_segment(data = sampledTree, aes(x = L1, y = U2, xend = U1, yend = U2), color = 'black') +
    geom_rect(data = ShadeData,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                  fill = factor(sampledTree$model1 == model | sampledTree$model2 == model, levels = c(T,F))),
              inherit.aes = FALSE,
              alpha = 0.4) +
    scale_fill_manual(
      name = "Model Match", # Optional: adds a title to the legend
      values = c("TRUE" = "#2E7D32", "FALSE" = "#C62828")
    )  +
    xlab("Longitude") + ylab("Latitude") + theme(legend.position = "none")

}

tree = generateTree(199, c(0,0,1,1))

orderBoundaries_GeminiCleaned(tree)

plotTreeGrid(tree)

tree_sampled = sample_models_per_region(tree, 5, 0.1)

visualizeModelExistence(tree_sampled, 5)


