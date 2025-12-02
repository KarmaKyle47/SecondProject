library(plotly)
library(MCMCpack)
library(stringr)

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

compSpaceData = function(GMM, phySpace_Data){

  compSpace = data.frame(t(apply(matrix(1:nrow(phySpace_Data)), MARGIN = 1, FUN = function(i){get_compSpace_pos(GMM, c(phySpace_Data$X[i], phySpace_Data$Y[i]))})), phySpace_Data$Model)
  names(compSpace) = c("X","Y","Model")

  compSpace

}

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

evaluateHSGP = function(z, k, l, M, curPos){

  omega = (0:M)*pi

  spec_den = sqrt(2*pi)*l*exp(-0.5*l^2*omega^2)
  beta = k*diag(sqrt(spec_den)) %*% z %*% diag(sqrt(spec_den))

  phi_x = c(1, sqrt(2)*cos(omega[1:M+1]*curPos[1]))
  phi_y = c(1, sqrt(2)*cos(omega[1:M+1]*curPos[2]))

  as.numeric(phi_x %*% beta %*% phi_y)

}

invlogit <- function(x) {
  1 / (1 + exp(-x))
}

grid_res = 100
phySpaceBorder = c(-5,-5,5,5)
sampledGMM = sampleGMM(phySpaceBorder)
z_log = sampledHSGP$log_z[[1]]
k_log = sampledHSGP$log_k[1]
l_log = sampledHSGP$log_l[1]
z_logit = sampledHSGP$logit_z[[1]]
k_logit = sampledHSGP$logit_k[1]
l_logit = sampledHSGP$logit_l[1]
logit_prior_mean = 10
color_limits = c(0,2)

plotHSGP = function(grid_res, z_log, k_log, l_log, z_logit, k_logit, l_logit, logit_prior_mean, M, color_limits, GMM, phySpaceBorder){

  x_seq_comp = y_seq_comp = seq(0,1,length.out = grid_res+1)
  comp_grid = expand.grid(x_seq_comp, y_seq_comp)

  x_seq_phy = seq(phySpaceBorder[1],phySpaceBorder[3],length.out = grid_res+1)
  y_seq_phy = seq(phySpaceBorder[2],phySpaceBorder[4],length.out = grid_res+1)
  phy_grid = expand.grid(x_seq_phy, y_seq_phy)
  phy_grid_comp = compSpaceData(GMM, data.frame(X = phy_grid[,1], Y = phy_grid[,2], Model = NA))[,c(1,2)]

  values_coef_comp = matrix(exp(apply(comp_grid, MARGIN = 1, FUN = evaluateHSGP, z = z_log, k=k_log, l=l_log, M=M)), byrow = F, nrow = grid_res+1)
  values_weight_comp = matrix(invlogit(apply(comp_grid, MARGIN = 1, FUN = evaluateHSGP, z = z_logit, k=k_logit, l=l_logit, M=M) + logit_prior_mean), byrow = F, nrow = grid_res+1)
  values_traj_comp = values_coef_comp*values_weight_comp

  coef_surface_comp = plot_ly(x = x_seq_comp, y = y_seq_comp, z = values_coef_comp, type = "surface", color_scale = "Viridis")
  weight_surface_comp = plot_ly(x = x_seq_comp, y = y_seq_comp, z = values_weight_comp, type = "surface")
  traj_surface_comp = plot_ly(x = x_seq_comp, y = y_seq_comp, z = values_traj_comp, type = "surface", color_scale = "Viridis", cmin = color_limits[1], cmax = color_limits[2])

  values_coef_phy = matrix(exp(apply(phy_grid_comp, MARGIN = 1, FUN = evaluateHSGP, z = z_log, k=k_log, l=l_log, M=M)), byrow = F, nrow = grid_res+1)
  values_weight_phy = matrix(invlogit(apply(phy_grid_comp, MARGIN = 1, FUN = evaluateHSGP, z = z_logit, k=k_logit, l=l_logit, M=M) + logit_prior_mean), byrow = F, nrow = grid_res+1)
  values_traj_phy = values_coef_phy*values_weight_phy

  coef_surface_phy = plot_ly(x = x_seq_phy, y = y_seq_phy, z = values_coef_phy, type = "surface", color_scale = "Viridis")
  weight_surface_phy = plot_ly(x = x_seq_phy, y = y_seq_phy, z = values_weight_phy, type = "surface")
  traj_surface_phy = plot_ly(x = x_seq_phy, y = y_seq_phy, z = values_traj_phy, type = "surface", color_scale = "Viridis", cmin = color_limits[1], cmax = color_limits[2])



  list("Computational" = list("Coefficient" = coef_surface_comp, "Weight" = weight_surface_comp, "Trajectory" = traj_surface_comp),
       "Physical" = list("Coefficient" = coef_surface_phy, "Weight" = weight_surface_phy, "Trajectory" = traj_surface_phy))

}

N_models = 2
M=10

sampleFullTrajectoriesHSGP = function(N_models, M, log_k_alpha = 9, log_k_beta = 2,
                                      log_l_alpha = 4, log_l_beta = 1,
                                      logit_k_alpha = 2, logit_k_beta = 10,
                                      logit_l_alpha = 4, logit_l_beta = 1){

  log_ks = rinvgamma(n = N_models, shape = log_k_alpha, scale = log_k_beta)
  log_ls = rinvgamma(n = N_models, shape = log_l_alpha, scale = log_l_beta)

  logit_ks = rinvgamma(n = N_models, shape = logit_k_alpha, scale = logit_k_beta)
  logit_ls = rinvgamma(n = N_models, shape = logit_l_alpha, scale = logit_l_beta)

  log_zs = list()
  logit_zs = list()

  for(i in 1:N_models){

    log_zs[[i]] = matrix(rnorm((M+1)^2), nrow = M+1)
    logit_zs[[i]] = matrix(rnorm((M+1)^2), nrow = M+1)

  }

  list(log_z = log_zs, logit_z = logit_zs, log_k = log_ks, log_l = log_ls, logit_k = logit_ks, logit_l = logit_ls)


}

sampledHSGP = sampleFullTrajectoriesHSGP(2, 10)

plotFullTrajectoriesHSGP = function(sampledHSGP, grid_res, color_limits, GMM, phySpaceBorder){

  N_models = length(sampledHSGP$log_z)
  M = nrow(sampledHSGP$log_z[[1]])-1

  if(N_models == 2){
    logit_prior_mean = 10
  } else{
    logit_prior_mean = log(N_models/(N_models - 2))
  }

  plots = list()
  for(i in 1:N_models){

    plots[[i]] = plotHSGP(grid_res = grid_res, z_log = sampledHSGP$log_z[[i]], k_log = sampledHSGP$log_k[i], l_log = sampledHSGP$log_l[i],
                          z_logit = sampledHSGP$logit_z[[i]], k_logit = sampledHSGP$logit_k[i], l_logit = sampledHSGP$logit_l[i],
                          logit_prior_mean = logit_prior_mean, M = M, color_limits = color_limits, GMM = GMM, phySpaceBorder = phySpaceBorder)


  }

  plots

}

sampledHSGP$logit_k

sampledHSGP_plots = plotFullTrajectoriesHSGP(sampledHSGP, 100, color_limits = c(0,2), GMM = GMM, phySpaceBorder = phySpaceBorder)

sampledHSGP_plots[[2]]$Computational$Trajectory
sampledHSGP_plots[[2]]$Physical$Trajectory


M=10
z_log = matrix(rnorm((M+1)*(M+1)), nrow = M+1, ncol = M+1)
k_log = 0.175
l_log = 0.2
z_logit = matrix(rnorm((M+1)*(M+1)), nrow = M+1, ncol = M+1)
k_logit = 5
l_logit = 0.1


test_HSGP = plotHSGP(grid_res = 100, z_log, k_log, l_log, z_logit, k_logit, l_logit, M = 10)

test_HSGP$Coefficient
test_HSGP$Weight
test_HSGP$Trajectory

plotHSGP(100, coefs, prior_mean, k, l = 0.01, M)



TrajWeightedBaseVectorFields_HSGP = function(t, curPos, baseVectorFields,
                                             sampledHSGP,
                                             logit_prior_mean, M,
                                             GMM){
  N_models = length(sampledHSGP$log_z)
  M = nrow(sampledHSGP$log_z[[1]])-1

  if(N_models == 2){
    logit_prior_mean = 10
  } else{
    logit_prior_mean = log(N_models/(N_models - 2))
  }


  cur_CompPos = as.numeric(compSpaceData(GMM, data.frame(X = curPos[1],
                                                         Y = curPos[2],
                                                         Model = NA))[c(1,2)])

  cur_log_value = c()
  cur_logit_value = c()

  for(i in 1:N_models){

    cur_log_value = c(cur_log_value, evaluateHSGP(z = sampledHSGP$log_z[[i]], k = sampledHSGP$log_k[i], l = sampledHSGP$log_l[i], M = M, curPos = cur_CompPos))
    cur_logit_value = c(cur_logit_value, evaluateHSGP(z = sampledHSGP$logit_z[[i]], k = sampledHSGP$logit_k[i], l = sampledHSGP$logit_l[i], M = M, curPos = cur_CompPos) + logit_prior_mean)

  }

  cur_traj_value = exp(cur_log_value)*invlogit(cur_logit_value)

  cur_ModelVel = baseVectorFields(t, curPos)

  t(cur_ModelVel %*% matrix(cur_traj_value))

}

EulerMaruyama = function(startTime, startPos, baseVectorFields,
                         sampledHSGP,
                         logit_prior_mean, M,
                         vel_sigma = 0.1, GMM, n_obs = 100, t_step_mean = 0.1){

  n_dim = length(startPos)

  t_sim = c(startTime, startTime + cumsum(rexp(n_obs-1, rate = 1/t_step_mean)))

  pos_sim = matrix(startPos, nrow = 1, byrow = T)

  for(i in 1:(n_obs-1)){

    cur_t = t_sim[i]
    cur_t_step = t_sim[i+1] - t_sim[i]
    cur_pos = pos_sim[i,]

    drift = TrajWeightedBaseVectorFields_HSGP(cur_t, cur_pos, baseVectorFields,
                                              sampledHSGP,
                                              logit_prior_mean, M, GMM)

    diffusion = rnorm(n_dim, mean = 0, sd = vel_sigma)

    pos_sim = rbind(pos_sim, cur_pos + drift*cur_t_step + diffusion*sqrt(cur_t_step))

  }

  full_sim = data.frame(cbind(t_sim, pos_sim))

  names(full_sim) = c('t', stringr::str_c('X', 1:n_dim))

  full_sim
}

test_particle = EulerMaruyama(startTime = 0, startPos = c(-1,1), baseVectorFields, sampledHSGP, logit_prior_mean = 10, M = 10, vel_sigma = 0.1, GMM = GMM, n_obs = 100, t_step_mean = 0.1)

plot(test_particle[,2:3])


samplePhySpaceParticles = function(n_particles, startTime, n_obs, phySpaceBorder, phySpaceBorderBuffer = 0.1, baseVectorFields,
                                   sampledHSGP,
                                   logit_prior_mean, M, GMM, t_step_mean = 0.01, vel_sigma = 0.1, pos_sigma = 0.01){

  x_phy_width = phySpaceBorder[3] - phySpaceBorder[1]
  y_phy_width = phySpaceBorder[4] - phySpaceBorder[2]

  startPos = data.frame(X = runif(n_particles, min = phySpaceBorder[1] + phySpaceBorderBuffer*x_phy_width, max = phySpaceBorder[3] - phySpaceBorderBuffer*x_phy_width),
                        Y = runif(n_particles, min = phySpaceBorder[2] + phySpaceBorderBuffer*y_phy_width, max = phySpaceBorder[4] - phySpaceBorderBuffer*y_phy_width))

  particleData_List = list()

  for(i in 1:n_particles){

    particleData_List[[i]] = cbind(EulerMaruyama(startTime = startTime, startPos = c(startPos$X[i], startPos$Y[i]), baseVectorFields = baseVectorFields,
                                                 sampledHSGP = sampledHSGP, logit_prior_mean = logit_prior_mean, M = M,
                                                 vel_sigma = vel_sigma, GMM = GMM, n_obs = n_obs, t_step_mean = t_step_mean), str_c("Particle",i))

    svMisc::progress(i, n_particles)
  }


  particleData_True = data.frame(do.call(rbind, particleData_List))

  names(particleData_True) = c('t', 'X1','X2', 'Particle')

  particleData_PosError = matrix(rnorm(2*nrow(particleData_True), mean = 0, sd = pos_sigma), ncol = 2)

  particleData_Obs = particleData_True
  particleData_Obs[,c(2:3)] = particleData_True[,c(2:3)] + particleData_PosError

  particleData_Obs

}




sampledParticles = samplePhySpaceParticles(100, startTime = 0, n_obs = 10, phySpaceBorder = c(-5,-5,5,5), phySpaceBorderBuffer = 0.2, baseVectorFields, sampledHSGP,
                        logit_prior_mean = 10, M = 10, GMM = sampledGMM, t_step_mean = 0.01, vel_sigma = 0.1, pos_sigma = 0.01)


ggplot(sampledParticles, aes(x = X1, y = X2, color = Particle)) + geom_point()









