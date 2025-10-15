

updatePatchesByModelExistence = function(sampledTree, sampledModelsBoundaries, model){

  curCoefs = sampledTree$coefs

  for(i in 1:nrow(sampledTree$boundaries)){

    if(sampledModelsBoundaries$model1[i] == model | sampledModelsBoundaries$model2[i] == model){

      curCoefs[i,] = rep(0, 16)

    }


  }

  updatedTree = sampledTree
  updatedTree$coefs = curCoefs

  updatedTree
}





tree = generate_grid_tree(0.1,c(0,0,1,1))
plotTreeGrid(tree)

sampledTree = samplePatch_FullTree(tree, 0.1)
plotCubicPatch3D(sampledTree,grid_size = 0.01, z_limit = c(-5,6))

sampledModelsBoundaries = sample_models_one_pass(tree, 5, baseWeight = 0.1)

updatedTree_1 = updatePatchesByModelExistence(sampledTree, sampledModelsBoundaries, model = 1)
updatedTree_2 = updatePatchesByModelExistence(sampledTree, sampledModelsBoundaries, model = 2)
updatedTree_3 = updatePatchesByModelExistence(sampledTree, sampledModelsBoundaries, model = 3)
updatedTree_4 = updatePatchesByModelExistence(sampledTree, sampledModelsBoundaries, model = 4)
updatedTree_5 = updatePatchesByModelExistence(sampledTree, sampledModelsBoundaries, model = 5)
plotCubicPatch3D(updatedTree_4,grid_size = 0.01, z_limit = c(-5,6))
