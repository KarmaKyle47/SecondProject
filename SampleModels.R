find_neighbors <- function(all_boundaries, index, sep = F, include_corners = F) {

  # Coordinates of the current boundary
  b <- all_boundaries[index, ]

  # Find boundaries that share an edge
  # Note: The small tolerance `1e-9` handles potential floating-point inaccuracies.
  if(include_corners){

    right_neighbors <- which(abs(all_boundaries$L1 - b$U1) <= 0 & all_boundaries$L2 <= b$U2 & all_boundaries$U2 >= b$L2)
    left_neighbors  <- which(abs(all_boundaries$U1 - b$L1) <= 0 & all_boundaries$L2 <= b$U2 & all_boundaries$U2 >= b$L2)
    up_neighbors    <- which(abs(all_boundaries$L2 - b$U2) <= 0 & all_boundaries$L1 <= b$U1 & all_boundaries$U1 >= b$L1)
    down_neighbors  <- which(abs(all_boundaries$U2 - b$L2) <= 0 & all_boundaries$L1 <= b$U1 & all_boundaries$U1 >= b$L1)

  } else{
    right_neighbors <- which(abs(all_boundaries$L1 - b$U1) < 1e-9 & all_boundaries$L2 < b$U2 & all_boundaries$U2 > b$L2)
    left_neighbors  <- which(abs(all_boundaries$U1 - b$L1) < 1e-9 & all_boundaries$L2 < b$U2 & all_boundaries$U2 > b$L2)
    up_neighbors    <- which(abs(all_boundaries$L2 - b$U2) < 1e-9 & all_boundaries$L1 < b$U1 & all_boundaries$U1 > b$L1)
    down_neighbors  <- which(abs(all_boundaries$U2 - b$L2) < 1e-9 & all_boundaries$L1 < b$U1 & all_boundaries$U1 > b$L1)
  }

  right_edges = apply(cbind(all_boundaries$L2[right_neighbors],
                      all_boundaries$U2[right_neighbors],
                      rep(b$L2, length(right_neighbors)),
                      rep(b$U2, length(right_neighbors))), MARGIN = 1, FUN = function(row){diff(sort(row)[c(2,3)])})
  left_edges = apply(cbind(all_boundaries$L2[left_neighbors],
                           all_boundaries$U2[left_neighbors],
                           rep(b$L2, length(left_neighbors)),
                           rep(b$U2, length(left_neighbors))), MARGIN = 1, FUN = function(row){diff(sort(row)[c(2,3)])})
  up_edges = apply(cbind(all_boundaries$L1[up_neighbors],
                            all_boundaries$U1[up_neighbors],
                            rep(b$L1, length(up_neighbors)),
                            rep(b$U1, length(up_neighbors))), MARGIN = 1, FUN = function(row){diff(sort(row)[c(2,3)])})
  down_edges = apply(cbind(all_boundaries$L1[down_neighbors],
                           all_boundaries$U1[down_neighbors],
                           rep(b$L1, length(down_neighbors)),
                           rep(b$U1, length(down_neighbors))), MARGIN = 1, FUN = function(row){diff(sort(row)[c(2,3)])})

  all_neighbors = c(right_neighbors, left_neighbors, up_neighbors, down_neighbors)
  all_edges = c(right_edges, left_edges, up_edges, down_edges)

  if(sep){

    return(list(Right = cbind(right_neighbors, right_edges),
                Left = cbind(left_neighbors, left_edges),
                Up = cbind(up_neighbors, up_edges),
                Down = cbind(down_neighbors, down_edges)))

  } else{
    return(cbind(all_neighbors, all_edges))
  }


}

sample_models_one_pass = function(tree, num_models, baseWeight = 0.1){

  treeBoundaries = data.frame(orderBoundaries_GeminiCleaned(tree)[[1]], model1 = 0, model2 = 0)

  baseTransMat = (baseWeight/num_models)*(matrix(rep(1, num_models^2), nrow = num_models) - diag(num_models))

  for(i in 1:nrow(treeBoundaries)){

    curTransMat = baseTransMat

    curIndex = which(treeBoundaries$order == i)
    neighbors_index_edge = find_neighbors(treeBoundaries, curIndex)

    definedNeighbors = which(treeBoundaries$order[neighbors_index_edge[,1]] < i)

    sum_of_defined_edges = sum(neighbors_index_edge[definedNeighbors, 2])

    if(sum_of_defined_edges > 0){

      for(j in definedNeighbors){

        curTransMat[treeBoundaries$model1[neighbors_index_edge[j,1]], treeBoundaries$model2[neighbors_index_edge[j,1]]] = curTransMat[treeBoundaries$model2[neighbors_index_edge[j,1]], treeBoundaries$model1[neighbors_index_edge[j,1]]] = curTransMat[treeBoundaries$model1[neighbors_index_edge[j,1]], treeBoundaries$model2[neighbors_index_edge[j,1]]] + neighbors_index_edge[j,2]/sum_of_defined_edges


      }

    }


    model1 = sample.int(num_models, size = 1, prob = rowSums(curTransMat))
    model2 = sample.int(num_models, size = 1, prob = curTransMat[model1,])

    treeBoundaries$model1[curIndex] = model1
    treeBoundaries$model2[curIndex] = model2


  }

  treeBoundaries

}

Energy_Pair = function(node1_models, node2_models){

  # models_similar = sum(node1_models %in% node2_models)
  #
  # if(models_similar == 2){
  #   0
  # } else if(models_similar == 1){
  #   1/2
  # } else if(models_similar == 0){
  #   1
  # }

  1 - (sum(node1_models %in% node2_models)/length(unique(c(node1_models,node2_models))))

}

Energy_Self = function(node_models, penatly){

  penatly*(node_models[1] == node_models[2])

}

sample_models_gibbs = function(tree, num_models, num_passes, self_penalty = 10000, temperature = 1){

  init_models = t(replicate(nrow(tree$split)+1, sample.int(num_models, 2)))

  treeBoundaries = data.frame(treeBorders(tree), model1 = init_models[,1], model2 = init_models[,2])

  tree_energy = c(calculate_tree_energy(treeBoundaries, self_penalty = self_penalty, temperature = temperature))

  one_pass = function(order, treeBoundaries){

    for(i in order){

      cur_models = c(treeBoundaries$model1[i], treeBoundaries$model2[i])
      cand_models = sample.int(num_models, 2)
      neighbors = find_neighbors(treeBoundaries, i)

      delta_self = Energy_Self(cand_models, penatly = self_penalty) - Energy_Self(cur_models, penatly = self_penalty)

      cur_edges_energy = 0
      cand_edges_energy = 0

      for(j in 1:nrow(neighbors)){

        cur_neighbor_models = c(treeBoundaries$model1[neighbors[j,1]], treeBoundaries$model2[neighbors[j,1]])
        cur_edge_length = neighbors[j,2]

        cur_edges_energy = cur_edges_energy + Energy_Pair(cur_models, cur_neighbor_models)*cur_edge_length
        cand_edges_energy = cand_edges_energy + Energy_Pair(cand_models, cur_neighbor_models)*cur_edge_length

      }

      delta_pairs = cand_edges_energy - cur_edges_energy

      delta_E = (delta_self + delta_pairs)/temperature

      target_prob = min(1, exp(-delta_E))
      prob = runif(1)

      if(prob <= target_prob){

        treeBoundaries$model1[i] = cand_models[1]
        treeBoundaries$model2[i] = cand_models[2]

      }

    }

    treeBoundaries

  }

  for(pass in 1:num_passes){

    order = sample(1:nrow(treeBoundaries))

    cur_pass = one_pass(order, treeBoundaries)

    treeBoundaries = cur_pass
    tree_energy = c(tree_energy, calculate_tree_energy(treeBoundaries, self_penalty, 1))

    svMisc::progress(pass, num_passes)

  }

  list(treeBoundaries, tree_energy)

}
sampledTree = baseTreeSampled

calculate_tree_energy = function(sampledTree, self_penalty = 10000, temperature = 1){

  self_energy = 0
  pair_energy = 0

  for(i in 1:nrow(sampledTree)){

    cur_models = c(sampledTree$model1[i], sampledTree$model2[i])
    neighbors = find_neighbors(sampledTree, i)

    self_energy = self_energy + Energy_Self(cur_models, penatly = self_penalty)

    for(j in 1:nrow(neighbors)){

      cur_neighbor_models = c(sampledTree$model1[neighbors[j,1]], sampledTree$model2[neighbors[j,1]])
      cur_edge_length = neighbors[j,2]

      pair_energy = pair_energy + Energy_Pair(cur_models, cur_neighbor_models)

    }

  }

  unname((self_energy + pair_energy/2)/temperature)

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
    xlab("X") + ylab("Y") + theme(legend.position = "none")

}

test_models = sample_models_one_pass(tree, 3)

model1 = visualizeModelExistence(test_models, 1) + ggtitle("Model 1") + theme(plot.title = element_text(hjust = 0.5))
model2 = visualizeModelExistence(test_models, 2) + ggtitle("Model 2") + theme(plot.title = element_text(hjust = 0.5))
model3 = visualizeModelExistence(test_models, 3) + ggtitle("Model 3") + theme(plot.title = element_text(hjust = 0.5))

sampleTransitionSurface_ByGrid(tree, sampledModels = test_models[,c(6,7)], model = 1, trans_prop = 8/9, k = 0.1, prior_mean = 1)

ggarrange(model1, model2, model3, nrow = 1)

tree = generateTree(99, c(0,0,1,1))

treeBoundaries = orderBoundaries_GeminiCleaned(tree)[[1]]

plotTreeGrid(tree)

one_pass_dist = c()

for(i in 1:10000){

  one_pass_dist = c(one_pass_dist, calculate_tree_energy(sample_models_per_region(tree, 5, 0.1), self_penalty = 10000, temperature = 1))

  svMisc::progress(i, 10000)
}

hist(one_pass_dist)


tree_sampled_one_pass = sample_models_per_region(tree, 5, 0.1)

tree_sampled_gibbs = sample_models_gibbs(tree, 5, num_passes = 1000, self_penalty = 1000, temperature = 0.05)


visualizeModelExistence(tree_sampled_gibbs[[1]], 5)
visualizeModelExistence(tree_sampled_one_pass, 1)

calculate_tree_energy(tree_sampled_one_pass, self_penalty = 10000, temperature = 0.04)
calculate_tree_energy(tree_sampled_gibbs[[1]], self_penalty = 10000, temperature = 0.04)


plot(tree_sampled_gibbs[[2]])
hist(tree_sampled_gibbs[[2]][100:1000])

length(tree_sampled_gibbs[[2]])

kindaBoundardies$Model1 = "Off"
kindaBoundardies$Model2 = "Off"
kindaBoundardies$Model3 = "Off"

for(i in 1:nrow(kindaBoundardies)){

  cur_trueIndex = which(tree$boundaries$L1 == kindaBoundardies$L1[i] & tree$boundaries$L2 == kindaBoundardies$L2[i])

  if(tree$models$Model1[cur_trueIndex] == "On"){
    kindaBoundardies$Model1[i] = "On"
  }
  if(tree$models$Model2[cur_trueIndex] == "On"){
    kindaBoundardies$Model2[i] = "On"
  }
  if(tree$models$Model3[cur_trueIndex] == "On"){
    kindaBoundardies$Model3[i] = "On"
  }


}

model1 = plotTransitionRegions(kindaBoundardies, types = kindaBoundardies$Model1, title = "Model 1") +     scale_fill_manual(
  name = "Model Match", # Optional: adds a title to the legend
  values = c("On" = "#2E7D32", "Off" = "#C62828")
) + xlab("X") + ylab("Y")

model2 = plotTransitionRegions(kindaBoundardies, types = kindaBoundardies$Model2, title = "Model 2") +     scale_fill_manual(
  name = "Model Match", # Optional: adds a title to the legend
  values = c("On" = "#2E7D32", "Off" = "#C62828")
) + xlab("X") + ylab("Y")

model3 = plotTransitionRegions(kindaBoundardies, types = kindaBoundardies$Model3, title = "Model 3") +     scale_fill_manual(
  name = "Model Match", # Optional: adds a title to the legend
  values = c("On" = "#2E7D32", "Off" = "#C62828")
) + xlab("X") + ylab("Y")

ggarrange(model1, model2, model3, nrow = 1)

model1_Trans = plotTransitionRegions(tree$boundaries, types = tree$models$Model1, title = "Model 1")
model2_Trans = plotTransitionRegions(tree$boundaries, types = tree$models$Model2, title = "Model 2")
model3_Trans = plotTransitionRegions(tree$boundaries, types = tree$models$Model3, title = "Model 3")

ggarrange(model1_Trans, model2_Trans, model3_Trans, nrow = 1)

plotTree = plotFullTree(tree, 0.01)
plotTree$ModelSurfaces[[2]]

ggplot(testPrior$ObsData, aes(x = X1, y = X2, color = Particle)) + geom_point() +theme(legend.position = "none")



## New sample models
# Instead of sampling which 2 models are on (discrete) each model gets a probability of being on (encoded with a logit)
# Still using the boltzmann-type distribution
#Self energy is the L2 norm of the sum of the probabilities and 2 (encourage 2 models to be "on" in expectation)
#Pairwise energy is the L2 norms of the logits

library(MASS)
num_models = 5
sample_new_models_one_pass = function(tree, num_models, logit_sd = 0.1){

  treeBoundaries = data.frame(orderBoundaries_GeminiCleaned(tree)[[1]])

  RegionLogits = matrix(nrow = nrow(treeBoundaries), ncol = num_models)

  if(num_models == 2){
    baseLogits = c(1000, 1000)
  } else{
    baseLogits = c(log(9),log(9), rep(log((0.2/(num_models - 2))/(1-0.2/(num_models - 2))), num_models - 2))
  }


  for(i in 1:nrow(treeBoundaries)){

    curIndex = which(treeBoundaries$order == i)
    neighbors_index_edge = find_neighbors(treeBoundaries, curIndex)

    definedNeighbors = neighbors_index_edge[which(treeBoundaries$order[neighbors_index_edge[,1]] < i),1]

    if(length(definedNeighbors) > 0){

      defined_logits = matrix(RegionLogits[definedNeighbors,], ncol = num_models)
      avg_neighbor_logits = colMeans(defined_logits)


      RegionLogits[curIndex,] = mvrnorm(1, mu = avg_neighbor_logits, Sigma = logit_sd^2*diag(num_models))

    } else{

      RegionLogits[curIndex,] = mvrnorm(1, mu = sample(baseLogits), Sigma = logit_sd^2*diag(num_models))

    }


  }

  RegionLogits

}

visualizeNewModelExistence = function(Boundaries, ModelLogits, model){

  ShadeData = data.frame(xmin = Boundaries$L1,
                         xmax = Boundaries$U1,
                         ymin = Boundaries$L2,
                         ymax = Boundaries$U2,
                         Region = str_c("Region", 1:nrow(Boundaries)),
                         Prob = 1/(1 + exp(-1*ModelLogits[,model])))

  ggplot() +
    geom_segment(data = Boundaries, aes(x = L1, y = L2, xend = L1, yend = U2), color = 'black') +
    geom_segment(data = Boundaries, aes(x = U1, y = L2, xend = U1, yend = U2), color = 'black') +
    geom_segment(data = Boundaries, aes(x = L1, y = L2, xend = U1, yend = L2), color = 'black') +
    geom_segment(data = Boundaries, aes(x = L1, y = U2, xend = U1, yend = U2), color = 'black') +
    geom_rect(data = ShadeData,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Prob),
              inherit.aes = FALSE,
              alpha = 0.4) +
    scale_fill_gradient(
      low = "#C62828", high = "#2E7D32", limits = c(0,1),
      name = str_c("Probability of Model", model)
    )  +
    xlab("X") + ylab("Y")# + theme(legend.position = "none")

}



tree = generate_grid_tree(0.1, c(0,0,1,1))

ModelLogits = sample_new_models_one_pass(tree, 5, 5)

visualizeNewModelExistence(Boundaries = tree$boundaries, ModelLogits = ModelLogits, model = 4)

Energy_Pair_New = function(cell1_logits, cell2_logits){

  mean((cell1_logits - cell2_logits)^2)

}

Energy_Self_New = function(cell_logits, penatly){

  cell_probs = 1/(1+exp(-1*cell_logits))

  penatly*mean((sum(cell_probs) - 2)^2)

}

calculate_tree_energy_new = function(treeBoundaries, cellLogits, self_penalty = 10, temperature = 1){

  self_energy = 0
  pair_energy = 0

  for(i in 1:nrow(treeBoundaries)){

    neighbors = find_neighbors(treeBoundaries, i)

    self_energy = self_energy + Energy_Self_New(cellLogits[i,], penatly = self_penalty)

    for(j in 1:nrow(neighbors)){

      pair_energy = pair_energy + Energy_Pair_New(cellLogits[i,], cellLogits[j,])

    }

  }

  unname((self_energy + pair_energy/2)/temperature)

}



test_energy = c()

for(i in 1:1000){

  test_energy = c(test_energy, calculate_tree_energy_new(treeBoundaries = tree$boundaries,
                                                         cellLogits = sample_new_models_one_pass(tree, 5, 5), self_penalty = 10, temperature = 1))

  svMisc::progress(i, 1000)
}

plot(test_energy, type = 'l')



