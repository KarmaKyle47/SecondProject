

updatePatchesByModelExistence = function(sampledTree, sampledModelsBoundaries, model){

  curCoefs = sampledTree$coefs

  for(i in 1:nrow(sampledTree$boundaries)){

    if(sampledModelsBoundaries$model1[i] != model & sampledModelsBoundaries$model2[i] != model){

      curCoefs[i,] = rep(0, 16)

    }


  }

  updatedTree = sampledTree
  updatedTree$coefs = curCoefs

  updatedTree
}

sampleTransitions = function(modelUpdatedTree, trans_prop = 0.3){

  #Constructing all possible boundary regions

  if(missing(trans_prop)){
    trans_prop = rbeta(1, 10, 90)
  }

  x_grid_size = modelUpdatedTree$boundaries[1,"U1"] - modelUpdatedTree$boundaries[1,"L1"]
  y_grid_size = modelUpdatedTree$boundaries[1,"U2"] - modelUpdatedTree$boundaries[1,"L2"]

  trans_length = (x_grid_size + y_grid_size - sqrt((x_grid_size + y_grid_size)^2 - 4*trans_prop*x_grid_size*y_grid_size))/4

  newBoundaries = modelUpdatedTree$boundaries
  newBoundaries$Type = "On"
  newBoundaries$Type[rowSums(modelUpdatedTree$coefs == 0) == 16] = "Off"
  newBoundaries$OGIndex = 1:nrow(newBoundaries)

  newSplits = modelUpdatedTree$split

  x_cells_num = (modelUpdatedTree$border[3] - modelUpdatedTree$border[1])/x_grid_size
  y_cells_num = (modelUpdatedTree$border[4] - modelUpdatedTree$border[2])/x_grid_size

  for(i in 1:x_cells_num){

    curSplitLocLow = (i-1)*x_grid_size + trans_length
    curSplitLocHigh = i*x_grid_size - trans_length
    cellsToSplit = which(newBoundaries$L1 < curSplitLocLow & newBoundaries$U1 > curSplitLocHigh)
    cellTypes = newBoundaries[cellsToSplit,]$Type
    cellOGIndex = newBoundaries[cellsToSplit,]$OGIndex
    n_cells = length(cellsToSplit)

    curSplits = data.frame(splits = c(rep(curSplitLocLow, n_cells), rep(curSplitLocHigh, n_cells)), dim = 1, split_part = c(cellsToSplit, 1:n_cells + nrow(newBoundaries)))
    newSplits = rbind(newSplits, curSplits)

    newBoundaries[cellsToSplit,]$U1 = curSplitLocLow
    newBoundaries[cellsToSplit,]$Type = "T"

    curBoundaries = data.frame(L1 = c(rep(curSplitLocLow, n_cells), rep(curSplitLocHigh, n_cells)), L2 = rep(newBoundaries[cellsToSplit,]$L2,2), U1 = c(rep(curSplitLocHigh, n_cells), rep(curSplitLocHigh + trans_length, n_cells)), U2 = rep(newBoundaries[cellsToSplit,]$U2,2), order = NA, Type = c(cellTypes, rep("T", n_cells)), OGIndex = rep(cellOGIndex, 2))

    newBoundaries = rbind(newBoundaries, curBoundaries)

  }


  for(i in 1:y_cells_num){

    curSplitLocLow = (i-1)*y_grid_size + trans_length
    curSplitLocHigh = i*y_grid_size - trans_length
    cellsToSplit = which(newBoundaries$L2 < curSplitLocLow & newBoundaries$U2 > curSplitLocHigh)
    cellTypes = newBoundaries[cellsToSplit,]$Type
    cellOGIndex = newBoundaries[cellsToSplit,]$OGIndex
    n_cells = length(cellsToSplit)

    curSplits = data.frame(splits = c(rep(curSplitLocLow, n_cells), rep(curSplitLocHigh, n_cells)), dim = 2, split_part = c(cellsToSplit, 1:n_cells + nrow(newBoundaries)))
    newSplits = rbind(newSplits, curSplits)

    newBoundaries[cellsToSplit,]$U2 = curSplitLocLow
    newBoundaries[cellsToSplit,]$Type = "T"

    curBoundaries = data.frame(L1 = rep(newBoundaries[cellsToSplit,]$L1,2), L2 = c(rep(curSplitLocLow, n_cells), rep(curSplitLocHigh, n_cells)), U1 = rep(newBoundaries[cellsToSplit,]$U1,2), U2 = c(rep(curSplitLocHigh, n_cells), rep(curSplitLocHigh + trans_length, n_cells)), order = NA, Type = c(cellTypes, rep("T", n_cells)), OGIndex = rep(cellOGIndex, 2))

    newBoundaries = rbind(newBoundaries, curBoundaries)

  }

  #Defining true boundary regions

  transition_index = which(newBoundaries$Type == "T")

  neighbors_FirstPass = apply(matrix(transition_index), MARGIN = 1,
                              FUN = function(i){
                                curNeighbors = find_neighbors(all_boundaries = newBoundaries, index = i, sep = F, include_corners = T)
                                neighborTypes = newBoundaries$Type[curNeighbors]

                                if(sum(neighborTypes == "On") > 0){
                                  "On"
                                } else{
                                  "T"
                                }

                                })

  newBoundaries$Type[transition_index] = neighbors_FirstPass

  transition_index = which(newBoundaries$Type == "T")

  neighbors_SecondPass = apply(matrix(transition_index), MARGIN = 1,
                              FUN = function(i){
                                curNeighbors = find_neighbors(all_boundaries = newBoundaries, index = i, sep = F, include_corners = T)
                                neighborTypes = newBoundaries$Type[curNeighbors]

                                if(sum(neighborTypes == "On") == 0){
                                  "Off"
                                } else{
                                  "T"
                                }

                              })
  newBoundaries$Type[transition_index] = neighbors_SecondPass

  #Defining new coefficients

  newCoefs = matrix(0, nrow = nrow(newBoundaries), ncol = 16)

  on_Index = which(newBoundaries$Type == "On")
  t_Index = which(newBoundaries$Type == "T")


  for(i in on_Index){

    curBorder = unname(unlist(newBoundaries[i,1:4]))
    k=0

    cur_BLcoefs = cur_BRcoefs = cur_TLcoefs = cur_TRcoefs = unname(unlist(modelUpdatedTree$coefs[newBoundaries$OGIndex[i],]))
    cur_BLborder = cur_BRborder = cur_TLborder = curTRborder = unname(unlist(modelUpdatedTree$boundaries[newBoundaries$OGIndex[i],1:4]))

    newCoefs[i,] = samplePatch_Known_Corners(curBorder, k, cur_BLcoefs, cur_BLborder, cur_BRcoefs, cur_BRborder, cur_TLcoefs, cur_TLborder, cur_TRcoefs, curTRborder)

  }

  for(i in t_Index){

    curBorder = unname(unlist(newBoundaries[i,1:4]))
    k=0

    curNeighbors = find_neighbors(newBoundaries, index = i, sep = T, include_corners = T)
    curNieghbors_NoT = lapply(curNeighbors, FUN = function(mat){
      matrix(mat[newBoundaries$Type[mat[,1]] != "T",], ncol = 2)
    })

    curBL_index = c(curNieghbors_NoT$Down[,1], curNieghbors_NoT$Left[,1])[1]
    curTL_index = c(curNieghbors_NoT$Up[,1], curNieghbors_NoT$Left[,1])[1]
    curBR_index = c(curNieghbors_NoT$Down[,1], curNieghbors_NoT$Right[,1])[1]
    curTR_index = c(curNieghbors_NoT$Up[,1], curNieghbors_NoT$Right[,1])[1]

    cur_BLborder = unname(unlist(newBoundaries[curBL_index,1:4]))
    cur_TLborder = unname(unlist(newBoundaries[curTL_index,1:4]))
    cur_BRborder = unname(unlist(newBoundaries[curBR_index,1:4]))
    cur_TRborder = unname(unlist(newBoundaries[curTR_index,1:4]))

    cur_BLcoefs = unname(unlist(newCoefs[curBL_index,]))
    cur_TLcoefs = unname(unlist(newCoefs[curTL_index,]))
    cur_BRcoefs = unname(unlist(newCoefs[curBR_index,]))
    cur_TRcoefs = unname(unlist(newCoefs[curTR_index,]))

    newCoefs[i,] = samplePatch_Known_Corners(curBorder, k, cur_BLcoefs, cur_BLborder, cur_BRcoefs, cur_BRborder, cur_TLcoefs, cur_TLborder, cur_TRcoefs, cur_TRborder)

  }



  TransitionTree = list(split = newSplits,
                        boundaries = newBoundaries,
                        dimension = modelUpdatedTree$dimension,
                        border = modelUpdatedTree$border,
                        coefs = newCoefs)
  TransitionTree
}





tree = generate_grid_tree(0.1,c(0,0,1,1))
plotTreeGrid(tree)

sampledTree = samplePatch_FullTree(tree, 0.001)
plotCubicPatch3D(sampledTree,grid_size = 0.01, z_limit = c(-5,6))

sampledModelsBoundaries = sample_models_one_pass(tree, 5, baseWeight = 0.1)

updatedTree_1 = updatePatchesByModelExistence(sampledTree, sampledModelsBoundaries, model = 1)
updatedTree_2 = updatePatchesByModelExistence(sampledTree, sampledModelsBoundaries, model = 2)
updatedTree_3 = updatePatchesByModelExistence(sampledTree, sampledModelsBoundaries, model = 3)
updatedTree_4 = updatePatchesByModelExistence(sampledTree, sampledModelsBoundaries, model = 4)
updatedTree_5 = updatePatchesByModelExistence(sampledTree, sampledModelsBoundaries, model = 5)
visualizeModelExistence(sampledModelsBoundaries, model = 1)

plotCubicPatch3D(transitionTree,grid_size = 0.001)

modelUpdatedTree = updatedTree_1

transitionTree = sampleTransitions(updatedTree_1, trans_prop = 0.9)

treeBorders(transitionTree)

plotTreeGrid(transitionTree)

Borders = transitionTree$boundaries

ShadeData = data.frame(xmin = Borders$L1,
                       xmax = Borders$U1,
                       ymin = Borders$L2,
                       ymax = Borders$U2,
                       Region = Borders$Type)

plot = ggplot() +
  geom_segment(data = Borders, aes(x = L1, y = L2, xend = L1, yend = U2), color = 'black') +
  geom_segment(data = Borders, aes(x = U1, y = L2, xend = U1, yend = U2), color = 'black') +
  geom_segment(data = Borders, aes(x = L1, y = L2, xend = U1, yend = L2), color = 'black') +
  geom_segment(data = Borders, aes(x = L1, y = U2, xend = U1, yend = U2), color = 'black') +
  geom_rect(data = ShadeData,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                fill = Region),
            inherit.aes = FALSE,
            alpha = 0.4) +
  xlab("Longitude") + ylab("Latitude") + theme(legend.position = "none")

plot

