baseVectorFields = function(t, curPos){

  f1 = c(curPos[2],-1*curPos[1])/sqrt(sum(c(curPos[1],curPos[2])^2))
  f2 = c(curPos[1],curPos[2])/sqrt(sum(c(curPos[1],curPos[2])^2))
  matrix(c(f1,f2), nrow = 2, byrow = F)

}

compPatchTree = sampleFullCubicSurface_AllModels(k = 0.5, num_models = 2, border_length = 1, cells_per_dim = 10, base_weight = 0.1, trans_prop = 8/9, prior_mean = 1)

phySpaceBorder = c(-5,-5,5,5)
GMM = sampleGMM(phySpaceBorder)
plotCompPatch = plotFullTree(compPatchTree)
plotTransPatch = plotTransformedFullTree(fullTree = compPatchTree, GMM = GMM, phySpaceBorder, surface_grid_size = 0.1, region_grid_size = 0.01)

plotTransPatch$ModelSurfaces[[2]]

curPos = c(-3.4, 2.3)
t = 0

TrajWeightedBaseVectorFields = function(t, curPos, baseVectorFields, compPatchTree, GMM){

  cur_CompPos = as.numeric(compSpaceData(GMM, data.frame(X = curPos[1],
                                                         Y = curPos[2],
                                                         Model = NA))[c(1,2)])

  cur_traj_value = as.numeric(lapply(X = compPatchTree$coefs, FUN = evaluateSampledTreeValue, tree = compPatchTree, curPos = cur_CompPos))

  cur_ModelVel = baseVectorFields(t, curPos)

  t(cur_ModelVel %*% matrix(cur_traj_value))

}

RungeKutta = function(startTime, startPos, baseVectorFields, compPatchTree, vel_sigma = 0.1, GMM, endTime, t_step = 0.01){

  n_dim = length(startPos)

  t_sim = seq(startTime, endTime, t_step)

  pos_sim = matrix(startPos, nrow = 1, byrow = T)

  for(i in 1:(length(t_sim)-1)){

    cur_t = t_sim[i]
    cur_pos = pos_sim[i,]

    k1 = TrajWeightedBaseVectorFields(cur_t, cur_pos, baseVectorFields, compPatchTree, GMM)
    k2 = TrajWeightedBaseVectorFields(cur_t + 0.5*t_step, cur_pos + t_step*0.5*k1, baseVectorFields, compPatchTree, GMM)
    k3 = TrajWeightedBaseVectorFields(cur_t + 0.5*t_step, cur_pos + t_step*0.5*k2, baseVectorFields, compPatchTree, GMM)
    k4 = TrajWeightedBaseVectorFields(cur_t + t_step, cur_pos + t_step*k3, baseVectorFields, compPatchTree, GMM)

    pos_sim = rbind(pos_sim, cur_pos + t_step*(k1+2*k2+2*k3+k4)/6 + rnorm(n_dim, mean = 0, sd = vel_sigma)*sqrt(t_step))

  }

  full_sim = data.frame(cbind(t_sim, pos_sim))

  names(full_sim) = c('t', stringr::str_c('X', 1:n_dim))

  full_sim
}

EulerMaruyama = function(startTime, startPos, baseVectorFields, compPatchTree, vel_sigma = 0.1, GMM, endTime, t_step = 0.01){

  n_dim = length(startPos)

  t_sim = seq(startTime, endTime, t_step)

  pos_sim = matrix(startPos, nrow = 1, byrow = T)

  for(i in 1:(length(t_sim)-1)){

    cur_t = t_sim[i]
    cur_pos = pos_sim[i,]

    drift = TrajWeightedBaseVectorFields(cur_t, cur_pos, baseVectorFields, compPatchTree, GMM)
    diffusion = rnorm(n_dim, mean = 0, sd = vel_sigma)

    pos_sim = rbind(pos_sim, cur_pos + drift*t_step + diffusion*sqrt(t_step))

  }

  full_sim = data.frame(cbind(t_sim, pos_sim))

  names(full_sim) = c('t', stringr::str_c('X', 1:n_dim))

  full_sim
}

samplePhySpaceParticles = function(n_particles, startTime, endTime, phySpaceBorder, phySpaceBorderBuffer = 0.1, baseVectorFields, compPatchTree, GMM, t_step = 0.01, vel_sigma = 0.1, pos_sigma = 0.1){

  x_phy_width = phySpaceBorder[3] - phySpaceBorder[1]
  y_phy_width = phySpaceBorder[4] - phySpaceBorder[2]

  startPos = data.frame(X = runif(n_particles, min = phySpaceBorder[1] + phySpaceBorderBuffer*x_phy_width, max = phySpaceBorder[3] - phySpaceBorderBuffer*x_phy_width),
                        Y = runif(n_particles, min = phySpaceBorder[2] + phySpaceBorderBuffer*y_phy_width, max = phySpaceBorder[4] - phySpaceBorderBuffer*y_phy_width))

  particleData_List = list()

  for(i in 1:n_particles){

    particleData_List[[i]] = cbind(RungeKutta(startTime = startTime, startPos = c(startPos$X[i], startPos$Y[i]), baseVectorFields = baseVectorFields, compPatchTree = compPatchTree, vel_sigma = vel_sigma, GMM = GMM, endTime = endTime, t_step = t_step), str_c("Particle",i))

    svMisc::progress(i, n_particles)
  }


  particleData_True = data.frame(do.call(rbind, particleData_List))

  names(particleData_True) = c('t', 'X1','X2', 'Particle')

  particleData_PosError = matrix(rnorm(2*nrow(particleData_True), mean = 0, sd = pos_sigma), ncol = 2)

  particleData_Obs = particleData_True
  particleData_Obs[,c(2:3)] = particleData_True[,c(2:3)] + particleData_PosError

  particleData_Obs

}

GMM = sampleGMM(c(-5,-5,5,5))

t1 = Sys.time()
testData_RK = RungeKutta(0, c(1,3), baseVectorFields, compPatchTree = compPatchTree, GMM = GMM, endTime = 10, t_step = 0.01, vel_sigma = 0)
t2 = Sys.time()
testData_EM = EulerMaruyama(0, c(1,3), baseVectorFields, compPatchTree = compPatchTree, GMM = GMM, endTime = 10, t_step = 0.01, vel_sigma = 0)
t3 = Sys.time()

ggplot(testData_RK, aes(x = X1, y = X2)) + geom_point() + xlim(0,15) + ylim(-4,4)
ggplot(testData_EM, aes(x = X1, y = X2)) + geom_point() + xlim(0,15) + ylim(-4,4)


testMultiPart = samplePhySpaceParticles(n_particles = 100, startTime = 0, endTime = 10,
                                        phySpaceBorder = phySpaceBorder, phySpaceBorderBuffer = 0.1,
                                        baseVectorFields = baseVectorFields, compPatchTree = compPatchTree, GMM = GMM, t_step = 0.1, vel_sigma = 0.1, pos_sigma = 0.1)

ggplot(testMultiPart, aes(x = X1, y = X2, color = Particle)) + geom_point()
