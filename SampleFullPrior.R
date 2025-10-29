sampleFromFullPrior = function(phySpaceBorder, baseVectorFields, n_particles, regions_per_dim_CompSpace = 10, traj_k_param = 5, base_weight_params = c(100,900),
                               trans_proportion_params = c(80,10), traj_mean = 1,
                               error_k_param = 10000, n_GMM_mixtures_param = 19, GMM_buffer_perc_params = c(10,40), GMM_cov_params = c(10,40), startTime = 0, endTime = 5,
                               particle_t_step = 0.01, pos_error_param = 10000, vel_error_param = 10000){

  #Generate Computational Space Trajectories

  n_models = ncol(baseVectorFields(t = 0, curPos = c(0.01,0.01))) #Fixed

  #regions_per_dim_CompSpace -> maybe add extra hyperprior

  traj_k = sqrt(rinvgamma(1, shape = traj_k_param)) #std on the value and first derivatives of the trajectory surfaces at the corners of the grid
  error_k = sqrt(rinvgamma(1, shape = error_k_param)) #std on the value and first derivatives of the error surfaces at the corners of the grid (not sure about this one)
  base_weight_model_transition = rbeta(1, base_weight_params[1], base_weight_params[2]) #controls the stickiness of the models in the computational space
  trans_proportion = rbeta(1, trans_proportion_params[1], trans_proportion_params[2]) #the percent of each grid cell to be allocated to be a transition from On to Off

  cat("Sampling Vector Field Trajectory Surfaces...")

  modelTraj = sampleFullCubicSurface_AllModels(border_length = 1, cells_per_dim = regions_per_dim_CompSpace,
                                               num_models = n_models, k = traj_k, base_weight = base_weight_model_transition, #Trajectory Surfaces for the Vector Fields
                                               trans_prop = trans_proportion, prior_mean = traj_mean)

  #Generate Physical Space Transformation

  cat("Sampling Gaussian Mixture Model for Transition...")

  phySpaceBorder_Buffer = rbeta(1, GMM_buffer_perc_params[1], GMM_buffer_perc_params[2])/2 #the percent of the physical space border to treat as a buffer
  GMM_Marginal_Variance_Perc = rbeta(1, GMM_cov_params[1], GMM_cov_params[2])/2 #the std for the marginal variance of each mixture component (proportion of the physical space width)

  phySpaceGMM = sampleGMM(phy_space_border = phySpaceBorder, lambda = n_GMM_mixtures_param,
                          border_buffer_perc = phySpaceBorder_Buffer, dim_sigma = GMM_Marginal_Variance_Perc) #Sampled Gaussian Mixture Model to act as the non-uniform density in the physical space


  #Generate Particle Data

  cat("Sampling Particle Data...")

  PosErrorVar = rinvgamma(1, pos_error_param) #variance for the positional errors in the data (not sure about this either)
  VelErrorVar = rinvgamma(1, vel_error_param) #variance for the velocity stochastic errors in the data (not sure about this either)


  ParticleData = samplePhySpaceParticles(n_particles = n_particles, startTime = startTime, endTime = endTime,
                                      phySpaceBorder = phySpaceBorder, phySpaceBorderBuffer = phySpaceBorder_Buffer, #Samples true particle locations
                                      baseVectorFields = baseVectorFields,
                                      compPatchTree = modelTraj, GMM = phySpaceGMM, t_step = particle_t_step, vel_sigma = sqrt(VelErrorVar), pos_sigma = sqrt(PosErrorVar))
  ParticleData = na.omit(ParticleData)

  hyperparameters = c(traj_k = traj_k,
                      base_weight_model_transition = base_weight_model_transition,
                      transition_proportion = trans_proportion,
                      phy_space_border_buffer = phySpaceBorder_Buffer,
                      gmm_marginal_variance_proportion = GMM_Marginal_Variance_Perc,
                      pos_error_variance = PosErrorVar, pos_error_variance = VelErrorVar)

  list(ObsData = ParticleData, TrajectorySurfaces = modelTraj, TransformationGMM = phySpaceGMM, Hyperparameters = hyperparameters)

}


FullPriorSample = sampleFromFullPrior(phySpaceBorder = c(-10,-10,10,10), baseVectorFields = baseVectorFields, n_particles = 10, regions_per_dim_CompSpace = 10, startTime = 0, endTime = 10, particle_t_step = 0.01)
ggplot(FullPriorSample$ObsData, aes(x = X1, y = X2, color = Particle)) + geom_point() + theme(legend.position = "none")

FullPriorSample$Hyperparameters
