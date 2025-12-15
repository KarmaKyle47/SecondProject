library(plotly)
library(MCMCpack)
library(stringr)

evaluateHSGP = function(z, k, l, M, border, curPos){

  Lx = border[3] - border[1]
  Ly = border[4] - border[2]

  omega = (0:M)*pi

  spec_den = sqrt(2*pi)*l*exp(-0.5*l^2*omega^2)

  beta = k*diag(sqrt(spec_den)) %*% z %*% diag(sqrt(spec_den))

  phi_x = c(1, sqrt(2)*cos(omega[1:M+1]*(curPos[1] - border[1])/Lx))
  phi_y = c(1, sqrt(2)*cos(omega[1:M+1]*(curPos[2] - border[2])/Ly))

  as.numeric(phi_x %*% beta %*% phi_y)

}

curTime = 1.1

evaluateHSGP_1D = function(z, k, l, M, boundary, curTime, prior_mean){

  Lt = diff(boundary)

  omega = (0:M)*pi

  spec_den = sqrt(2*pi)*l*exp(-0.5*l^2*omega^2)

  beta = k*sqrt(spec_den)*z

  phi = c(1, sqrt(2)*cos(omega[1:M+1]*(curTime - boundary[1])/Lt))

  sum(phi * beta) + prior_mean

}

evaluateHSGP_1D_Derivative = function(z, k, l, M, boundary, curTime){

  Lt = diff(boundary)

  omega = (1:M)*pi

  spec_den = sqrt(2*pi)*l*exp(-0.5*l^2*omega^2)

  beta = k*sqrt(spec_den)*z[1:M+1] * (-1*omega/Lt)

  phi = sqrt(2)*sin(omega[1:M]*(curTime - boundary[1])/Lt)

  sum(phi * beta)

}

FullTrajectoriesHSGP = sampleFullTrajectoriesHSGP(2, 10)
border = c(-5,-5,5,5)

t_res = 100

getMeanCoefsEst = function(FullTrajectoriesHSGP, sampledPath, baseVectorFields, border, boundary, t_res){

  M_Path = length(sampledPath$x_zs) - 1
  M_Traj = nrow(FullTrajectoriesHSGP$log_z[[1]])-1

  Lt = diff(boundary)

  t_grid = seq(boundary[1], boundary[2], length.out = t_res + 1)

  delta_t = t_grid[2] - t_grid[1]

  x_seq = apply(matrix(t_grid), MARGIN = 1, FUN = evaluateHSGP_1D, z = sampledPath$x_zs, k = sampledPath$x_k, l = sampledPath$x_l, M = M_Path, boundary = boundary)
  y_seq = apply(matrix(t_grid), MARGIN = 1, FUN = evaluateHSGP_1D, z = sampledPath$y_zs, k = sampledPath$y_k, l = sampledPath$y_l, M = M_Path, boundary = boundary)

  all_pos_mat = cbind(t_grid, x_seq, y_seq)

  mu_v = apply(all_pos_mat, MARGIN = 1, FUN = function(row){
    TrajWeightedBaseVectorFields_HSGP(t = row[1], curPos = row[2:3], baseVectorFields = baseVectorFields, sampledHSGP = FullTrajectoriesHSGP, M = M_Traj, border = border)
  })

  sin_b = sqrt(2) * sin(outer((t_grid - boundary[1])/Lt, 1:M_Path * pi))

  c = t(mu_v %*% sin_b * delta_t) / Lt

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

plotHSGP = function(grid_res, z_log, k_log, l_log, M, color_limits, border){

  x_seq = seq(border[1],border[3],length.out = grid_res+1)
  y_seq = seq(border[2],border[4],length.out = grid_res+1)
  plot_grid = expand.grid(x_seq, y_seq)

  values_coef = matrix(exp(apply(plot_grid, MARGIN = 1, FUN = evaluateHSGP, z = z_log, k=k_log, l=l_log, M=M, border=border)), byrow = F, nrow = grid_res+1)

  coef_surface = plot_ly(x = x_seq, y = y_seq, z = values_coef, type = "surface", color_scale = "Viridis")


  list("Coefficient" = coef_surface, "CoefValues" = values_coef)

}

N_models = 2
M=10

sampleFullTrajectoriesHSGP = function(N_models, M, log_k_alpha = 9, log_k_beta = 2,
                                      log_l_alpha = 4, log_l_beta = 1,
                                      logit_k_alpha = 4, logit_k_beta = 1,
                                      logit_l_alpha = 10, logit_l_beta = 1){

  # log_ks = rinvgamma(n = N_models, shape = log_k_alpha, scale = log_k_beta)
  # log_ls = rinvgamma(n = N_models, shape = log_l_alpha, scale = log_l_beta)

  log_ks = rep(0.175, N_models)
  log_ls = rep(0.2, N_models)


  log_zs = list()

  for(i in 1:N_models){

    log_zs[[i]] = matrix(rnorm((M+1)^2), nrow = M+1)

  }

  list(log_z = log_zs, log_k = log_ks, log_l = log_ls)


}

samplePathHSGP = function(M, x_k_alpha = 4, x_k_beta = 10,
                               x_l_alpha = 20, x_l_beta = 1,
                               y_k_alpha = 4, y_k_beta = 10,
                               y_l_alpha = 20, y_l_beta = 1){

  x_k = rinvgamma(n=1,shape = x_k_alpha, scale = x_k_beta)
  x_l = rinvgamma(n=1,shape = x_l_alpha, scale = x_l_beta)

  y_k = rinvgamma(n=1,shape = y_k_alpha, scale = y_k_beta)
  y_l = rinvgamma(n=1,shape = y_l_alpha, scale = y_l_beta)


  x_zs = rnorm(M+1)
  y_zs = rnorm(M+1)

  list(x_zs = x_zs, y_zs = y_zs, x_k = x_k, x_l= x_l, y_k = y_k, y_l = y_l)

}

sampledPath = samplePathHSGP(10)

boundary = c(1, 10)
t_res = 1000

plotPathHSGP = function(sampledPath, t_res, boundary){

  M = length(sampledPath$x_zs) - 1

  t_grid = seq(boundary[1], boundary[2], length.out = t_res + 1)

  x_pos = apply(matrix(t_grid), MARGIN = 1, FUN = evaluateHSGP_1D, z = sampledPath$x_zs, k = sampledPath$x_k, l = sampledPath$x_l, M = M, boundary = boundary)
  y_pos = apply(matrix(t_grid), MARGIN = 1, FUN = evaluateHSGP_1D, z = sampledPath$y_zs, k = sampledPath$y_k, l = sampledPath$y_l, M = M, boundary = boundary)

  pos_data = data.frame(t = t_grid, x = x_pos, y = y_pos)

  ggplot(data = pos_data, aes(x = x, y = y)) + geom_point()

}


plotFullTrajectoriesHSGP = function(sampledHSGP, grid_res, color_limits, border){

  N_models = length(sampledHSGP$log_z)
  M = nrow(sampledHSGP$log_z[[1]])-1

  plots = list()
  for(i in 1:N_models){

    plots[[i]] = plotHSGP(grid_res = grid_res, z_log = sampledHSGP$log_z[[i]], k_log = sampledHSGP$log_k[i], l_log = sampledHSGP$log_l[i],
                          M = M, color_limits = color_limits, border = border)


  }

  plots

}

sampledHSGP = sampleFullTrajectoriesHSGP(2, 10)
border = c(-5,-5,5,5)

sampledHSGP_plots = plotFullTrajectoriesHSGP(sampledHSGP, 100, color_limits = c(0,2), border = border)

sampledHSGP_plots[[1]]$Coefficient
sampledHSGP_plots[[5]]$Trajectory

sampledHSGP$log_z[[1]]


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

baseVectorFields = function(t, curPos){

  f1 = c(curPos[2],-1*curPos[1])/sqrt(sum(c(curPos[1],curPos[2])^2))
  f2 = c(curPos[1],curPos[2])/sqrt(sum(c(curPos[1],curPos[2])^2))
  # f3 = c(sqrt(2),sqrt(2))
  # f4 = c(-1*curPos[2],curPos[1])/sqrt(sum(c(curPos[1],curPos[2])^2))
  # f5 = c(-1*curPos[1],-1*curPos[2])/sqrt(sum(c(curPos[1],curPos[2])^2))
  matrix(c(f1,f2), nrow = 2, byrow = F)

}

TrajWeightedBaseVectorFields_HSGP = function(t, curPos, baseVectorFields,
                                             sampledHSGP,
                                             M, border){
  N_models = length(sampledHSGP$log_z)
  M = nrow(sampledHSGP$log_z[[1]])-1

  cur_log_value = c()

  for(i in 1:N_models){

    cur_log_value = c(cur_log_value, evaluateHSGP(z = sampledHSGP$log_z[[i]], k = sampledHSGP$log_k[i], l = sampledHSGP$log_l[i], M = M, curPos = curPos, border = border))

  }

  cur_traj_value = exp(cur_log_value)

  cur_ModelVel = baseVectorFields(t, curPos)

  t(cur_ModelVel %*% matrix(cur_traj_value))

}

EulerMaruyama = function(startTime, startPos, baseVectorFields,
                         sampledHSGP,
                         M, border,
                         vel_sigma = 0.1, n_obs = 100, t_step_mean = 0.1){

  n_dim = length(startPos)

  t_sim = c(startTime, startTime + cumsum(rexp(n_obs-1, rate = 1/t_step_mean)))

  pos_sim = matrix(startPos, nrow = 1, byrow = T)
  vel_sim = c()

  for(i in 1:(n_obs-1)){

    cur_t = t_sim[i]
    cur_t_step = t_sim[i+1] - t_sim[i]
    cur_pos = pos_sim[i,]

    drift = TrajWeightedBaseVectorFields_HSGP(cur_t, cur_pos, baseVectorFields,
                                              sampledHSGP, M, border)

    diffusion = rnorm(n_dim, mean = 0, sd = vel_sigma)

    pos_sim = rbind(pos_sim, cur_pos + drift*cur_t_step + diffusion*sqrt(cur_t_step))
    vel_sim = rbind(vel_sim, drift + diffusion)

  }

  vel_sim = rbind(vel_sim, TrajWeightedBaseVectorFields_HSGP(t_sim[n_obs], pos_sim[n_obs,], baseVectorFields,
                                                             sampledHSGP, M, border) + rnorm(n_dim, mean = 0, sd = vel_sigma))



  full_sim = data.frame(cbind(t_sim, pos_sim, vel_sim))

  names(full_sim) = c('t', stringr::str_c('X', 1:n_dim), stringr::str_c('X', 1:n_dim,'v'))

  full_sim
}

border = c(-5,-5,5,5)

test_particle = EulerMaruyama(startTime = 0, startPos = c(-2,2), baseVectorFields = baseVectorFields, sampledHSGP = sampledHSGP, M = 10, vel_sigma = 0.1, border = border, n_obs = 100, t_step_mean = 0.05)

plot(test_particle[,2:3])


samplePhySpaceParticles = function(n_particles, startTime, n_obs, border, borderBuffer = 0.1, baseVectorFields,
                                   sampledHSGP,
                                   M, t_step_mean = 0.01, vel_sigma = 0.1, pos_sigma = 0.01){

  Lx = border[3] - border[1]
  Ly = border[4] - border[2]

  startPos = data.frame(X = runif(n_particles, min = border[1] + borderBuffer*Lx, max = border[3] - borderBuffer*Lx),
                        Y = runif(n_particles, min = border[2] + borderBuffer*Ly, max = border[4] - borderBuffer*Ly))

  particleData_List = list()

  for(i in 1:n_particles){

    particleData_List[[i]] = cbind(EulerMaruyama(startTime = startTime, startPos = c(startPos$X[i], startPos$Y[i]), baseVectorFields = baseVectorFields,
                                                 sampledHSGP = sampledHSGP, M = M,
                                                 vel_sigma = vel_sigma, border = border, n_obs = n_obs, t_step_mean = t_step_mean), str_c("Particle",i))

    svMisc::progress(i, n_particles)
  }


  particleData_True = data.frame(do.call(rbind, particleData_List))

  names(particleData_True) = c('t', 'X1','X2','X1v','X2v', 'Particle')

  particleData_PosError = matrix(rnorm(2*nrow(particleData_True), mean = 0, sd = pos_sigma), ncol = 2)

  particleData_Obs = particleData_True
  particleData_Obs[,c(2:3)] = particleData_True[,c(2:3)] + particleData_PosError

  particleData_Obs

}




sampledParticles = samplePhySpaceParticles(1, startTime = 0, n_obs = 1000, border = border, borderBuffer = 0.2, baseVectorFields, sampledHSGP,
                        M = 10, t_step_mean = 0.1, vel_sigma = 0, pos_sigma = 0.25)

ggplot(sampledParticles, aes(x = X1, y = X2, color = Particle)) + geom_point()
ggplot(sampledParticles[abs(sampledParticles$X1)<=5 & abs(sampledParticles$X2)<=5, ], aes(x = X1, y = X2, color = Particle)) + geom_point()

sampledParticles_Subset = sampledParticles[abs(sampledParticles$X1)<=5 & abs(sampledParticles$X2)<=5, ]


sampledParticles_Subset = sampledParticles[sampledParticles$X1 >= border[1] & sampledParticles$X1 <= border[3] &
                                             sampledParticles$X2 >= border[2] & sampledParticles$X2 <= border[4],]
ggplot(sampledParticles_Subset, aes(x = X1, y = X2, color = Particle)) + geom_point()

testBorder = c(-5,-5,5,5)
sampledHSGP = sampleFullTrajectoriesHSGP(2, 10)

testData = cbind(0, runif(10000, testBorder[1], testBorder[3]), runif(10000, testBorder[2], testBorder[4]))

testData_V = t(apply(testData[,2:3], MARGIN = 1, FUN = TrajWeightedBaseVectorFields_HSGP, t = 0, baseVectorFields = baseVectorFields, sampledHSGP = sampledHSGP, logit_prior_mean = 10, M = 10, border = testBorder)) + matrix(rnorm(20000, 0, 0.01), nrow = 10000)

plot(testData_V)

testDataFull = cbind(testData, testData_V)
