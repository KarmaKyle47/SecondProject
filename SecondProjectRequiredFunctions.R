library(plotly)
library(MCMCpack)
library(stringr)
library(cmdstanr)
library(MASS)
library(ggplot2)

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

evaluateHSGP_1D = function(z, k, l, M, boundary, curTime, prior_mean){
  
  Lt = diff(boundary)
  
  omega = (0:M)*pi
  
  spec_den = sqrt(2*pi)*l*exp(-0.5*l^2*omega^2)
  
  beta = k*sqrt(spec_den)*z
  
  phi = c(1, sqrt(2)*cos(omega[1:M+1]*(curTime - boundary[1])/Lt))
  
  sum(phi * beta) + prior_mean
  
}

evaluateHSGP_1D_Times = function(z, k, l, M, boundary, times, prior_mean){
  
  Lt = diff(boundary)
  
  omega = (0:M)*pi
  
  spec_den = sqrt(2*pi)*l*exp(-0.5*l^2*omega^2)
  
  beta = k*sqrt(spec_den)*z
  
  phi = cbind(1, sqrt(2) * cos(outer((times - boundary[1])/Lt, omega[1:M+1])))

  as.numeric(phi %*% beta + prior_mean)
  
}


evaluateHSGP_1D_Derivative = function(z, k, l, M, boundary, curTime){
  
  Lt = diff(boundary)
  
  omega = (1:M)*pi
  
  spec_den = sqrt(2*pi)*l*exp(-0.5*l^2*omega^2)
  
  beta = k*sqrt(spec_den)*z[1:M+1] * (-1*omega/Lt)
  
  phi = sqrt(2)*sin(omega[1:M]*(curTime - boundary[1])/Lt)
  
  sum(phi * beta)
  
}

evaluateHSGP_1D_Derivative_Times = function(z, k, l, M, boundary, times){
  
  Lt = diff(boundary)
  
  omega = (1:M)*pi
  
  spec_den = sqrt(2*pi)*l*exp(-0.5*l^2*omega^2)
  
  beta = k*sqrt(spec_den)*z[1:M+1] * (-1*omega/Lt)
  
  phi = sqrt(2) * sin(outer((times - boundary[1])/Lt, omega[1:M]))
  
  as.numeric(phi %*% beta)
  
}

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

plotHSGP = function(grid_res, z_log, k_log, l_log, M, color_limits, border, grid_border){
  
  x_seq = seq(grid_border[1],grid_border[3],length.out = grid_res+1)
  y_seq = seq(grid_border[2],grid_border[4],length.out = grid_res+1)
  plot_grid = expand.grid(x_seq, y_seq)
  
  values_coef = matrix(exp(apply(plot_grid, MARGIN = 1, FUN = evaluateHSGP, z = z_log, k=k_log, l=l_log, M=M, border=border)), byrow = F, nrow = grid_res+1)
  
  coef_surface = plot_ly(x = x_seq, y = y_seq, z = values_coef, type = "surface", color_scale = "Viridis")
  
  
  list("Coefficient" = coef_surface, "CoefValues" = values_coef)
  
}

sampleFullTrajectoriesHSGP = function(N_models, M, log_k = 0.35, log_l = 0.2){
  
  # log_ks = rinvgamma(n = N_models, shape = log_k_alpha, scale = log_k_beta)
  # log_ls = rinvgamma(n = N_models, shape = log_l_alpha, scale = log_l_beta)
  
  log_ks = rep(log_k, N_models)
  log_ls = rep(log_l, N_models)
  
  
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

plotPathHSGP = function(sampledPath, t_res, boundary){
  
  M = length(sampledPath$x_zs) - 1
  
  t_grid = seq(boundary[1], boundary[2], length.out = t_res + 1)
  
  x_pos = apply(matrix(t_grid), MARGIN = 1, FUN = evaluateHSGP_1D, z = sampledPath$x_zs, k = sampledPath$x_k, l = sampledPath$x_l, M = M, boundary = boundary)
  y_pos = apply(matrix(t_grid), MARGIN = 1, FUN = evaluateHSGP_1D, z = sampledPath$y_zs, k = sampledPath$y_k, l = sampledPath$y_l, M = M, boundary = boundary)
  
  pos_data = data.frame(t = t_grid, x = x_pos, y = y_pos)
  
  ggplot(data = pos_data, aes(x = x, y = y)) + geom_point()
  
}

plotFullTrajectoriesHSGP = function(sampledHSGP, grid_res, color_limits, border, grid_border){
  
  N_models = length(sampledHSGP$log_z)
  M = nrow(sampledHSGP$log_z[[1]])-1
  
  plots = list()
  for(i in 1:N_models){
    
    plots[[i]] = plotHSGP(grid_res = grid_res, z_log = sampledHSGP$log_z[[i]], k_log = sampledHSGP$log_k[i], l_log = sampledHSGP$log_l[i],
                          M = M, color_limits = color_limits, border = border, grid_border = grid_border)
    
    
  }
  
  plots
  
}

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

getGridStanHSGP = function(fit, grid_res, M, n_iters, border, grid_border){
  
  Lx = border[3] - border[1]
  Ly = border[4] - border[2]
  
  # --- 1. Setup Constants ---
  # Number of coefficients per model
  N_basis = (M+1)^2
  omega = (0:M) * pi
  
  # Create Grid Vectors (0 to 1)
  x_seq = seq(grid_border[1], grid_border[3], length.out = grid_res + 1)
  y_seq = seq(grid_border[2], grid_border[4], length.out = grid_res + 1)
  
  # --- 2. Pre-calculate Basis Matrices (The Speedup) ---
  # We do this ONCE, not inside the loop.
  # Phi_x: [(grid_res+1) x (M+1)]
  
  # Helper to make cosine matrix
  make_phi <- function(coords, start, L) {
    # n=0 is column of 1s
    # n>0 is sqrt(2)*cos(...)
    # outer() creates the matrix of coords vs frequencies
    waves <- sqrt(2) * cos(outer((coords - start)/L, omega[2:(M+1)]))
    return(cbind(1, waves))
  }
  
  Phi_x <- make_phi(x_seq, border[1], Lx)
  Phi_y <- make_phi(y_seq, border[2], Ly)
  
  # --- 3. Extract Hyperparameters ---
  # We need the specific k and l for every draw to scale the Zs correctly
  # Assuming N_models = 2
  
  # Coefficient (Log) Hypers
  
  log_ks = fit$draws('log_ks', format = 'draws_matrix')
  log_ls = fit$draws('log_ls', format = 'draws_matrix')
  
  log_ks_1 <- as.numeric(log_ks[,1])
  log_ls_1 <- as.numeric(log_ls[,1])
  log_ks_2 <- as.numeric(log_ks[,2])
  log_ls_2 <- as.numeric(log_ls[,2])
  
  
  # --- 4. Extract Raw Z Matrices ---
  # Safer extraction using variable names
  
  log_draws <- fit$draws("log_zs", format = "draws_matrix")
  
  log_draws_1 = array(dim = c(M+1, M+1, n_iters))
  log_draws_2 = array(dim = c(M+1, M+1, n_iters))
  
  for(i in 1:n_iters){
    
    log_draws_1[,,i] = matrix(log_draws[i,1:(N_basis)*2 - 1], nrow = M+1, ncol = M+1, byrow = F)
    log_draws_2[,,i] = matrix(log_draws[i,1:(N_basis)*2], nrow = M+1, ncol = M+1, byrow = F)
  }
  
  # --- 5. Initialize Output Arrays ---
  traj_grid_draws_1 = array(dim = c(grid_res+1, grid_res+1, n_iters))
  traj_grid_draws_2 = array(dim = c(grid_res+1, grid_res+1, n_iters))

  
  # Helper for Spectral Density calculation
  get_beta_matrix <- function(z_flat, k, l, M) {
    # Calculate Spectral Density for this specific length scale
    # Squared Exponential Kernel
    # Note: Using Unit Variance first (k=1.0 inside), multiply by k at end
    spd_unit <- sqrt(2*pi) * l * exp(-0.5 * l^2 * omega^2)
    scale_vec <- sqrt(spd_unit)
    
    # Reshape Z to matrix
    z_mat <- matrix(z_flat, nrow = M+1, ncol = M+1)
    
    # Apply Scaling: D * Z * D
    beta_unit <- diag(scale_vec) %*% z_mat %*% diag(scale_vec)
    
    return(k * beta_unit)
  }
  
  # --- 6. The Main Loop (Now Fast) ---
  # Uses matrix multiplication instead of 'apply'
  
  print("Reconstructing surfaces...")
  pb <- txtProgressBar(min = 0, max = n_iters, style = 3)
  
  for(i in 1:n_iters){
    
    # --- Model 1 ---
    # A. Coefficients
    beta_log_1 <- get_beta_matrix(log_draws_1[,,i], log_ks_1[i], log_ls_1[i], M)
    # Tensor Projection: Phi_x * Beta * Phi_y'
    surf_log_1 <- Phi_x %*% beta_log_1 %*% t(Phi_y)
    traj_grid_draws_1[,,i] <- exp(surf_log_1) # Exponential Lin
    
    # --- Model 2 ---
    # A. Coefficients
    beta_log_2 <- get_beta_matrix(log_draws_2[,,i], log_ks_2[i], log_ls_2[i], M)
    surf_log_2 <- Phi_x %*% beta_log_2 %*% t(Phi_y)
    traj_grid_draws_2[,,i] <- exp(surf_log_2)
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Return list of results
  return(list(
    traj_1 = traj_grid_draws_1,
    traj_2 = traj_grid_draws_2,
    x = x_seq,
    y = y_seq
  ))
}

getGridStanHSGP_wLK = function(fit, grid_res, M, n_iters, border, grid_border, log_ks, log_ls){
  
  Lx = border[3] - border[1]
  Ly = border[4] - border[2]
  
  # --- 1. Setup Constants ---
  # Number of coefficients per model
  N_basis = (M+1)^2
  omega = (0:M) * pi
  
  # Create Grid Vectors (0 to 1)
  x_seq = seq(grid_border[1], grid_border[3], length.out = grid_res + 1)
  y_seq = seq(grid_border[2], grid_border[4], length.out = grid_res + 1)
  
  # --- 2. Pre-calculate Basis Matrices (The Speedup) ---
  # We do this ONCE, not inside the loop.
  # Phi_x: [(grid_res+1) x (M+1)]
  
  # Helper to make cosine matrix
  make_phi <- function(coords, start, L) {
    # n=0 is column of 1s
    # n>0 is sqrt(2)*cos(...)
    # outer() creates the matrix of coords vs frequencies
    waves <- sqrt(2) * cos(outer((coords - start)/L, omega[2:(M+1)]))
    return(cbind(1, waves))
  }
  
  Phi_x <- make_phi(x_seq, border[1], Lx)
  Phi_y <- make_phi(y_seq, border[2], Ly)
  
  # --- 3. Extract Hyperparameters ---
  # We need the specific k and l for every draw to scale the Zs correctly
  # Assuming N_models = 2
  
  log_ks_1 <- rep(log_ks[1], n_iters)
  log_ls_1 <- rep(log_ls[1], n_iters)
  log_ks_2 <- rep(log_ks[2], n_iters)
  log_ls_2 <- rep(log_ls[2], n_iters)
  
  
  # --- 4. Extract Raw Z Matrices ---
  # Safer extraction using variable names
  
  log_draws <- fit$draws("log_zs", format = "draws_matrix")
  
  log_draws_1 = array(dim = c(M+1, M+1, n_iters))
  log_draws_2 = array(dim = c(M+1, M+1, n_iters))
  
  for(i in 1:n_iters){
    
    log_draws_1[,,i] = matrix(log_draws[i,1:(N_basis)*2 - 1], nrow = M+1, ncol = M+1, byrow = F)
    log_draws_2[,,i] = matrix(log_draws[i,1:(N_basis)*2], nrow = M+1, ncol = M+1, byrow = F)
  }
  
  # --- 5. Initialize Output Arrays ---
  traj_grid_draws_1 = array(dim = c(grid_res+1, grid_res+1, n_iters))
  traj_grid_draws_2 = array(dim = c(grid_res+1, grid_res+1, n_iters))
  
  
  # Helper for Spectral Density calculation
  get_beta_matrix <- function(z_flat, k, l, M) {
    # Calculate Spectral Density for this specific length scale
    # Squared Exponential Kernel
    # Note: Using Unit Variance first (k=1.0 inside), multiply by k at end
    spd_unit <- sqrt(2*pi) * l * exp(-0.5 * l^2 * omega^2)
    scale_vec <- sqrt(spd_unit)
    
    # Reshape Z to matrix
    z_mat <- matrix(z_flat, nrow = M+1, ncol = M+1)
    
    # Apply Scaling: D * Z * D
    beta_unit <- diag(scale_vec) %*% z_mat %*% diag(scale_vec)
    
    return(k * beta_unit)
  }
  
  # --- 6. The Main Loop (Now Fast) ---
  # Uses matrix multiplication instead of 'apply'
  
  print("Reconstructing surfaces...")
  pb <- txtProgressBar(min = 0, max = n_iters, style = 3)
  
  for(i in 1:n_iters){
    
    # --- Model 1 ---
    # A. Coefficients
    beta_log_1 <- get_beta_matrix(log_draws_1[,,i], log_ks_1[i], log_ls_1[i], M)
    # Tensor Projection: Phi_x * Beta * Phi_y'
    surf_log_1 <- Phi_x %*% beta_log_1 %*% t(Phi_y)
    traj_grid_draws_1[,,i] <- exp(surf_log_1) # Exponential Lin
    
    # --- Model 2 ---
    # A. Coefficients
    beta_log_2 <- get_beta_matrix(log_draws_2[,,i], log_ks_2[i], log_ls_2[i], M)
    surf_log_2 <- Phi_x %*% beta_log_2 %*% t(Phi_y)
    traj_grid_draws_2[,,i] <- exp(surf_log_2)
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Return list of results
  return(list(
    traj_1 = traj_grid_draws_1,
    traj_2 = traj_grid_draws_2,
    x = x_seq,
    y = y_seq
  ))
}


getPointsStanHSGP = function(fit, points, M, n_iters, border){
  
  # points: N x 2 matrix or dataframe (col 1 = x, col 2 = y)
  # fit: The Stan fit object
  # M: Number of basis functions
  # border: c(xmin, ymin, xmax, ymax) used for normalization
  
  Lx = border[3] - border[1]
  Ly = border[4] - border[2]
  
  # --- 1. Setup Constants ---
  N_basis = (M+1)^2
  omega = (0:M) * pi
  
  # Extract point coordinates
  # Ensure it's a matrix for safe indexing
  points_mat <- as.matrix(points)
  x_coords <- points_mat[,1]
  y_coords <- points_mat[,2]
  N_points <- nrow(points_mat)
  
  # --- 2. Pre-calculate Basis Matrices ---
  # Phi_x and Phi_y are now [N_points x (M+1)]
  
  make_phi <- function(coords, start, L) {
    # n=0 is column of 1s
    # n>0 is sqrt(2)*cos(...)
    waves <- sqrt(2) * cos(outer((coords - start)/L, omega[2:(M+1)]))
    return(cbind(1, waves))
  }
  
  Phi_x <- make_phi(x_coords, border[1], Lx)
  Phi_y <- make_phi(y_coords, border[2], Ly)
  
  # --- 3. Extract Hyperparameters ---
  # (Identical to previous version)
  
  log_ks = fit$draws('log_ks', format = 'draws_matrix')
  log_ls = fit$draws('log_ls', format = 'draws_matrix')
  
  log_ks_1 <- as.numeric(log_ks[,1])
  log_ls_1 <- as.numeric(log_ls[,1])
  log_ks_2 <- as.numeric(log_ks[,2])
  log_ls_2 <- as.numeric(log_ls[,2])
  
  # --- 4. Extract Raw Z Matrices ---
  # (Identical to previous version)
  
  log_draws <- fit$draws("log_zs", format = "draws_matrix")
  
  log_draws_1 = array(dim = c(M+1, M+1, n_iters))
  log_draws_2 = array(dim = c(M+1, M+1, n_iters))
  
  for(i in 1:n_iters){
    log_draws_1[,,i] = matrix(log_draws[i,1:(N_basis)*2 - 1], nrow = M+1, ncol = M+1, byrow = F)
    log_draws_2[,,i] = matrix(log_draws[i,1:(N_basis)*2], nrow = M+1, ncol = M+1, byrow = F)
  }
  
  # --- 5. Initialize Output Arrays ---
  # CHANGE: Shape is now [N_points x n_iters] instead of a grid
  traj_points_draws_1 = matrix(NA, nrow = N_points, ncol = n_iters)
  traj_points_draws_2 = matrix(NA, nrow = N_points, ncol = n_iters)
  
  # Helper for Spectral Density calculation (Identical)
  get_beta_matrix <- function(z_flat, k, l, M) {
    spd_unit <- sqrt(2*pi) * l * exp(-0.5 * l^2 * omega^2)
    scale_vec <- sqrt(spd_unit)
    z_mat <- matrix(z_flat, nrow = M+1, ncol = M+1)
    beta_unit <- diag(scale_vec) %*% z_mat %*% diag(scale_vec)
    return(k * beta_unit)
  }
  
  # --- 6. The Main Loop ---
  
  print(paste("Evaluating trajectory for", N_points, "points..."))
  pb <- txtProgressBar(min = 0, max = n_iters, style = 3)
  
  for(i in 1:n_iters){
    
    # --- Model 1 ---
    beta_log_1 <- get_beta_matrix(log_draws_1[,,i], log_ks_1[i], log_ls_1[i], M)
    
    # MATH CHANGE: Point-wise evaluation
    # 1. Project Phi_x through Beta: (N x M) * (M x M) = (N x M)
    # 2. Element-wise multiply with Phi_y and sum rows
    # This calculates Sum(Phi_x_ik * Beta_kl * Phi_y_il) efficiently
    
    surf_log_1 <- rowSums((Phi_x %*% beta_log_1) * Phi_y)
    traj_points_draws_1[,i] <- exp(surf_log_1) 
    
    # --- Model 2 ---
    beta_log_2 <- get_beta_matrix(log_draws_2[,,i], log_ks_2[i], log_ls_2[i], M)
    
    surf_log_2 <- rowSums((Phi_x %*% beta_log_2) * Phi_y)
    traj_points_draws_2[,i] <- exp(surf_log_2)
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Return list of results
  # Returns matrices of (Points x Iterations)
  return(list(
    traj_1 = traj_points_draws_1,
    traj_2 = traj_points_draws_2,
    points = points_mat
  ))
}

getGridStanHSGP_Beta = function(surface_beta_draws, grid_res, M, n_iters, border, grid_border){
  
  Lx = border[3] - border[1]
  Ly = border[4] - border[2]
  
  # --- 1. Setup Constants ---
  # Number of coefficients per model
  N_basis = (M+1)^2
  omega = (0:M) * pi
  
  # Create Grid Vectors (0 to 1)
  x_seq = seq(grid_border[1], grid_border[3], length.out = grid_res + 1)
  y_seq = seq(grid_border[2], grid_border[4], length.out = grid_res + 1)
  
  # --- 2. Pre-calculate Basis Matrices (The Speedup) ---
  # We do this ONCE, not inside the loop.
  # Phi_x: [(grid_res+1) x (M+1)]
  
  # Helper to make cosine matrix
  make_phi <- function(coords, start, L) {
    # n=0 is column of 1s
    # n>0 is sqrt(2)*cos(...)
    # outer() creates the matrix of coords vs frequencies
    waves <- sqrt(2) * cos(outer((coords - start)/L, omega[2:(M+1)]))
    return(cbind(1, waves))
  }
  
  Phi_x <- make_phi(x_seq, border[1], Lx)
  Phi_y <- make_phi(y_seq, border[2], Ly)
  
  
  # --- 4. Extract Beta Matrices ---
  # Safer extraction using variable names
  
  # beta_draws <- fit$draws("surface_betas", format = "draws_matrix")
  
  beta_draws = surface_beta_draws
  
  beta_draws_1 = array(dim = c(M+1, M+1, n_iters))
  beta_draws_2 = array(dim = c(M+1, M+1, n_iters))
  
  for(i in 1:n_iters){
    
    beta_draws_1[,,i] = matrix(beta_draws[i,1:(N_basis)*2 - 1], nrow = M+1, ncol = M+1, byrow = F)
    beta_draws_2[,,i] = matrix(beta_draws[i,1:(N_basis)*2], nrow = M+1, ncol = M+1, byrow = F)
  }
  
  # --- 5. Initialize Output Arrays ---
  traj_grid_draws_1 = array(dim = c(grid_res+1, grid_res+1, n_iters))
  traj_grid_draws_2 = array(dim = c(grid_res+1, grid_res+1, n_iters))
  
  # --- 6. The Main Loop (Now Fast) ---
  # Uses matrix multiplication instead of 'apply'
  
  print("Reconstructing surfaces...")
  pb <- txtProgressBar(min = 0, max = n_iters, style = 3)
  
  for(i in 1:n_iters){
    
    # Tensor Projection: Phi_x * Beta * Phi_y'
    surf_log_1 <- Phi_x %*% beta_draws_1[,,i] %*% t(Phi_y)
    traj_grid_draws_1[,,i] <- exp(surf_log_1) # Exponential Lin
    
    surf_log_2 <- Phi_x %*% beta_draws_2[,,i] %*% t(Phi_y)
    traj_grid_draws_2[,,i] <- exp(surf_log_2)
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Return list of results
  return(list(
    traj_1 = traj_grid_draws_1,
    traj_2 = traj_grid_draws_2,
    x = x_seq,
    y = y_seq
  ))
}

getPointsStanHSGP_Beta = function(surface_beta_draws, points, M, n_iters, border){
  
  # points: N x 2 matrix or dataframe (col 1 = x, col 2 = y)
  # fit: The Stan fit object
  # M: Number of basis functions
  # border: c(xmin, ymin, xmax, ymax) used for normalization
  
  Lx = border[3] - border[1]
  Ly = border[4] - border[2]
  
  # --- 1. Setup Constants ---
  N_basis = (M+1)^2
  omega = (0:M) * pi
  
  # Extract point coordinates
  # Ensure it's a matrix for safe indexing
  points_mat <- as.matrix(points)
  x_coords <- points_mat[,1]
  y_coords <- points_mat[,2]
  N_points <- nrow(points_mat)
  
  # --- 2. Pre-calculate Basis Matrices ---
  # Phi_x and Phi_y are now [N_points x (M+1)]
  
  make_phi <- function(coords, start, L) {
    # n=0 is column of 1s
    # n>0 is sqrt(2)*cos(...)
    waves <- sqrt(2) * cos(outer((coords - start)/L, omega[2:(M+1)]))
    return(cbind(1, waves))
  }
  
  Phi_x <- make_phi(x_coords, border[1], Lx)
  Phi_y <- make_phi(y_coords, border[2], Ly)
  
  
  beta_draws_1 = array(dim = c(M+1, M+1, n_iters))
  beta_draws_2 = array(dim = c(M+1, M+1, n_iters))
  
  for(i in 1:n_iters){
    
    beta_draws_1[,,i] = matrix(surface_beta_draws[i,1:(N_basis)*2 - 1], nrow = M+1, ncol = M+1, byrow = F)
    beta_draws_2[,,i] = matrix(surface_beta_draws[i,1:(N_basis)*2], nrow = M+1, ncol = M+1, byrow = F)
  }
  
  # --- 5. Initialize Output Arrays ---
  # CHANGE: Shape is now [N_points x n_iters] instead of a grid
  traj_points_draws_1 = matrix(NA, nrow = N_points, ncol = n_iters)
  traj_points_draws_2 = matrix(NA, nrow = N_points, ncol = n_iters)

  # --- 6. The Main Loop ---
  
  print(paste("Evaluating trajectory for", N_points, "points..."))
  pb <- txtProgressBar(min = 0, max = n_iters, style = 3)
  
  for(i in 1:n_iters){
    
    # --- Model 1 ---
    beta_log_1 <- beta_draws_1[,,i]
    
    # MATH CHANGE: Point-wise evaluation
    # 1. Project Phi_x through Beta: (N x M) * (M x M) = (N x M)
    # 2. Element-wise multiply with Phi_y and sum rows
    # This calculates Sum(Phi_x_ik * Beta_kl * Phi_y_il) efficiently
    
    surf_log_1 <- rowSums((Phi_x %*% beta_log_1) * Phi_y)
    traj_points_draws_1[,i] <- exp(surf_log_1) 
    
    # --- Model 2 ---
    beta_log_2 <- beta_draws_2[,,i]
    
    surf_log_2 <- rowSums((Phi_x %*% beta_log_2) * Phi_y)
    traj_points_draws_2[,i] <- exp(surf_log_2)
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Return list of results
  # Returns matrices of (Points x Iterations)
  return(list(
    traj_1 = traj_points_draws_1,
    traj_2 = traj_points_draws_2,
    points = points_mat
  ))
}


fast_mean_3d <- function(arr_3d) {
  d <- dim(arr_3d)
  dim(arr_3d) <- c(d[1] * d[2], d[3])
  means <- rowMeans(arr_3d)
  dim(means) <- c(d[1], d[2])
  return(means)
}

fast_quantile_3d = function(arr_3d, q){
  d <- dim(arr_3d)
  dim(arr_3d) <- c(d[1] * d[2], d[3])
  qs <- apply(arr_3d, 1, FUN = quantile, q)
  dim(qs) <- c(d[1], d[2])
  return(qs)
}

