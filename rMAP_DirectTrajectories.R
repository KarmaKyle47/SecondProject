library(plotly)
library(stringr)
library(MASS)
library(ggplot2)
library(minqa)
library(nloptr)
library(mvtnorm)
library(deSolve)
library(Matrix)
library(pracma)


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

SplitStep_RK4 = function(startTime, startPos, baseVectorFields,
                         sampledHSGP,
                         M, border,
                         vel_sigma = 0.1, n_obs = 100, t_step_mean = 0.1){

  n_dim = length(startPos)

  t_sim = c(startTime, startTime + cumsum(rexp(n_obs-1, rate = 1/t_step_mean)))

  pos_sim = matrix(startPos, nrow = 1, byrow = T)
  vel_sim = c()

  for(i in 1:(n_obs-1)){

    cur_t = t_sim[i]
    h = t_sim[i+1] - t_sim[i] # Step size (delta t)
    cur_pos = pos_sim[i,]

    # --- 1. Deterministic RK4 Drift ---
    # Evaluate at current position
    k1 = TrajWeightedBaseVectorFields_HSGP(cur_t, cur_pos, baseVectorFields,
                                           sampledHSGP, M, border)

    # Evaluate at halfway point using k1
    k2 = TrajWeightedBaseVectorFields_HSGP(cur_t + h/2, cur_pos + k1 * (h/2), baseVectorFields,
                                           sampledHSGP, M, border)

    # Evaluate at halfway point using k2
    k3 = TrajWeightedBaseVectorFields_HSGP(cur_t + h/2, cur_pos + k2 * (h/2), baseVectorFields,
                                           sampledHSGP, M, border)

    # Evaluate at full step using k3
    k4 = TrajWeightedBaseVectorFields_HSGP(cur_t + h, cur_pos + k3 * h, baseVectorFields,
                                           sampledHSGP, M, border)

    # Weighted average of the tangents
    rk4_drift = (k1 + 2*k2 + 2*k3 + k4) / 6

    # --- 2. Additive Stochastic Diffusion ---
    diffusion = rnorm(n_dim, mean = 0, sd = vel_sigma)

    # --- 3. Combine ---
    pos_sim = rbind(pos_sim, cur_pos + rk4_drift*h + diffusion*sqrt(h))
    vel_sim = rbind(vel_sim, rk4_drift + diffusion)

  }

  # Final step velocity calc
  final_drift = TrajWeightedBaseVectorFields_HSGP(t_sim[n_obs], pos_sim[n_obs,], baseVectorFields,
                                                  sampledHSGP, M, border)
  vel_sim = rbind(vel_sim, final_drift + rnorm(n_dim, mean = 0, sd = vel_sigma))

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

    particleData_List[[i]] = cbind(SplitStep_RK4(startTime = startTime, startPos = c(startPos$X[i], startPos$Y[i]), baseVectorFields = baseVectorFields,
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

curPos

beta_mat = cbind(rep(0,4), c(-100,0,0,0))

exp(evaluate2DCosine_fast(beta_mat, pos_mat, border))

evaluate2DCosine_fast = function(beta_mat, pos_mat, border){

  if(is.null(dim(pos_mat))) {
    pos_mat = matrix(pos_mat, nrow = 1, ncol = 2)
  }

  Lx = border[3] - border[1]
  Ly = border[4] - border[2]

  M = sqrt(nrow(beta_mat))-1
  N_points = nrow(pos_mat)

  omega = (0:M)*pi

  # 1. Compute all spatial frequencies simultaneously
  X_scaled = (pos_mat[,1] - border[1]) / Lx
  Y_scaled = (pos_mat[,2] - border[2]) / Ly

  # outer() automatically creates a matrix of all combinations
  cos_X = cos(outer(X_scaled, omega))
  cos_Y = cos(outer(Y_scaled, omega))

  # 2. Apply the sqrt(2) scaling (1 for the first column, sqrt(2) for the rest)
  scale_vec = c(1, rep(sqrt(2), M))
  cos_X = sweep(cos_X, 2, scale_vec, `*`)
  cos_Y = sweep(cos_Y, 2, scale_vec, `*`)

  # 3. Create the Phi combinations instantly using R's fast vector recycling
  idx_i = rep(1:(M+1), times = M+1)
  idx_j = rep(1:(M+1), each = M+1)

  Phi = cos_X[, idx_i] * cos_Y[, idx_j]

  # Matrix multiply
  return(Phi %*% beta_mat)
}

evaluate2DCosine_part_x_fast = function(beta_mat, pos_mat, border){

  if(is.null(dim(pos_mat))) {
    pos_mat = matrix(pos_mat, nrow = 1, ncol = 2)
  }

  Lx = border[3] - border[1]
  Ly = border[4] - border[2]

  M = sqrt(nrow(beta_mat))-1

  omega = (0:M)*pi

  # 1. Compute all spatial frequencies
  X_scaled = (pos_mat[,1] - border[1]) / Lx
  Y_scaled = (pos_mat[,2] - border[2]) / Ly

  # d/dx(cos) -> -sin for X, Y remains cos
  sin_X = -sin(outer(X_scaled, omega))
  cos_Y = cos(outer(Y_scaled, omega))

  # 2. Apply scaling: sqrt(2) logic AND the chain rule for X (omega / Lx)
  c_scale = c(1, rep(sqrt(2), M))
  scale_vec_X = c_scale * (omega / Lx)

  sin_X = sweep(sin_X, 2, scale_vec_X, `*`)
  cos_Y = sweep(cos_Y, 2, c_scale, `*`)

  # 3. Create the Phi combinations
  idx_i = rep(1:(M+1), times = M+1)
  idx_j = rep(1:(M+1), each = M+1)

  Phi_dx = sin_X[, idx_i] * cos_Y[, idx_j]

  # Matrix multiply
  return(Phi_dx %*% beta_mat)
}

evaluate2DCosine_part_y_fast = function(beta_mat, pos_mat, border){

  if(is.null(dim(pos_mat))) {
    pos_mat = matrix(pos_mat, nrow = 1, ncol = 2)
  }

  Lx = border[3] - border[1]
  Ly = border[4] - border[2]

  M = sqrt(nrow(beta_mat))-1

  omega = (0:M)*pi

  # 1. Compute all spatial frequencies
  X_scaled = (pos_mat[,1] - border[1]) / Lx
  Y_scaled = (pos_mat[,2] - border[2]) / Ly

  # X remains cos, d/dy(cos) -> -sin for Y
  cos_X = cos(outer(X_scaled, omega))
  sin_Y = -sin(outer(Y_scaled, omega))

  # 2. Apply scaling: sqrt(2) logic AND the chain rule for Y (omega / Ly)
  c_scale = c(1, rep(sqrt(2), M))
  scale_vec_Y = c_scale * (omega / Ly)

  cos_X = sweep(cos_X, 2, c_scale, `*`)
  sin_Y = sweep(sin_Y, 2, scale_vec_Y, `*`)

  # 3. Create the Phi combinations
  idx_i = rep(1:(M+1), times = M+1)
  idx_j = rep(1:(M+1), each = M+1)

  Phi_dy = cos_X[, idx_i] * sin_Y[, idx_j]

  # Matrix multiply
  return(Phi_dy %*% beta_mat)
}

evaluate2DCosine_part_beta_fast = function(beta_mat, pos_mat, border){

  if(is.null(dim(pos_mat))) {
    pos_mat = matrix(pos_mat, nrow = 1, ncol = 2)
  }

  Lx = border[3] - border[1]
  Ly = border[4] - border[2]

  P = nrow(beta_mat) # Coefficients per surface
  K = ncol(beta_mat) # Number of surfaces

  M_degree = sqrt(P) - 1
  omega = (0:M_degree)*pi

  # 1. Compute spatial frequencies
  X_scaled = (pos_mat[,1] - border[1]) / Lx
  Y_scaled = (pos_mat[,2] - border[2]) / Ly

  cos_X = cos(outer(X_scaled, omega))
  cos_Y = cos(outer(Y_scaled, omega))

  # 2. Apply scaling
  scale_vec = c(1, rep(sqrt(2), M_degree))
  cos_X = sweep(cos_X, 2, scale_vec, `*`)
  cos_Y = sweep(cos_Y, 2, scale_vec, `*`)

  # 3. Create Phi (Dimensions: N x P)
  idx_i = rep(1:(M_degree+1), times = M_degree+1)
  idx_j = rep(1:(M_degree+1), each = M_degree+1)
  Phi = cos_X[, idx_i] * cos_Y[, idx_j]

  return(Phi)
}

TrajWeightedBaseVectorFields_2D_Cosine = function(pos_t_mat, beta_mat, baseVectorFields_Vec, border){

  if(is.null(dim(pos_t_mat))) {
    pos_t_mat = matrix(pos_t_mat, nrow = 1, ncol = 3)
  }

  log_Traj = evaluate2DCosine_fast(beta_mat = beta_mat, pos_mat = pos_t_mat[,c(2,3)], border = border)
  log_Traj_mat = cbind(log_Traj, log_Traj)

  baseVF = baseVectorFields_Vec(pos_t_mat[,])

  Weighted_Vel = exp(log_Traj_mat) * baseVF

  cbind(Weighted_Vel[,1] + Weighted_Vel[,2], Weighted_Vel[,3] + Weighted_Vel[,4])


}

baseVectorFields = function(t, curPos){

  c = sqrt(sum(curPos^2))

  f1 = c(curPos[2],-1*curPos[1]) / c
  f2 = c(curPos[1],curPos[2]) / c

  matrix(c(f1,f2), nrow = 2, byrow = F)

}

baseVectorFields_Vec = function(pos_t_mat){

  if(is.null(dim(pos_t_mat))) {
    pos_t_mat = matrix(pos_t_mat, nrow = 1, ncol = 3)
    pos_mat = matrix(pos_t_mat[2:3], nrow = 1, ncol = 2)
  } else{
    pos_mat = pos_t_mat[,c(2,3)]
  }

  c = sqrt(rowSums(pos_mat^2))

  f1x = pos_t_mat[,3] / c
  f1y = -1*pos_t_mat[,2] / c

  f2x = pos_t_mat[,2] / c
  f2y = pos_t_mat[,3] / c

  cbind(f1x,f2x, f1y, f2y)

}

baseVectorFields_Jacobian_Vec = function(pos_t_mat){

  if(is.null(dim(pos_t_mat))) {
    pos_t_mat = matrix(pos_t_mat, nrow = 1, ncol = 3)
    pos_mat = matrix(pos_t_mat[2:3], nrow = 1, ncol = 2)
  } else{
    pos_mat = pos_t_mat[,c(2,3)]
  }

  x = pos_t_mat[,2]
  y = pos_t_mat[,3]

  c = rowSums(pos_mat^2)^(3/2)

  part_x_f1x = -(x*y) / c
  part_y_f1x = (x^2) / c

  part_x_f1y = -(y^2) / c
  part_y_f1y = (x*y) / c

  part_x_f2x = (y^2) / c
  part_y_f2x = -(x*y) / c

  part_x_f2y = -(x*y) / c
  part_y_f2y = (x^2) / c

  cbind(part_x_f1x, part_y_f1x, part_x_f1y, part_y_f1y, part_x_f2x, part_y_f2x, part_x_f2y, part_y_f2y)

}

calculate_Jacobian_f_wrt_y = function(pos_t_mat, beta_mat, border){

  if(is.null(dim(pos_t_mat))) {
    pos_t_mat = matrix(pos_t_mat, nrow = 1, ncol = 2)
  }

  pos_mat = pos_t_mat[,2:3]

  traj = exp(evaluate2DCosine_fast(beta_mat, pos_mat, border))
  traj_part_x = traj*evaluate2DCosine_part_x_fast(beta_mat, pos_mat, border)
  traj_part_y = traj*evaluate2DCosine_part_y_fast(beta_mat, pos_mat, border)
  VF = baseVectorFields_Vec(pos_t_mat)
  VF_part = baseVectorFields_Jacobian_Vec(pos_t_mat)

  T1 = traj[,1]
  T2 = traj[,2]

  T1_part_x = traj_part_x[,1]
  T1_part_y = traj_part_y[,1]

  T2_part_x = traj_part_x[,2]
  T2_part_y = traj_part_y[,2]

  VF1_x = VF[,1]
  VF2_x = VF[,2]
  VF1_y = VF[,3]
  VF2_y = VF[,4]

  VF1_x_part_x = VF_part[,1]
  VF1_x_part_y = VF_part[,2]
  VF1_y_part_x = VF_part[,3]
  VF1_y_part_y = VF_part[,4]
  VF2_x_part_x = VF_part[,5]
  VF2_x_part_y = VF_part[,6]
  VF2_y_part_x = VF_part[,7]
  VF2_y_part_y = VF_part[,8]

  part_x_dx = T1_part_x * VF1_x + T1 * VF1_x_part_x + T2_part_x * VF2_x + T2 * VF2_x_part_x
  part_y_dx = T1_part_y * VF1_x + T1 * VF1_x_part_y + T2_part_y * VF2_x + T2 * VF2_x_part_y
  part_x_dy = T1_part_x * VF1_y + T1 * VF1_y_part_x + T2_part_x * VF2_y + T2 * VF2_y_part_x
  part_y_dy = T1_part_y * VF1_y + T1 * VF1_y_part_y + T2_part_y * VF2_y + T2 * VF2_y_part_y


  J_arr = array(dim = c(2,2,nrow(pos_t_mat)))

  J_arr[1,1,] = part_x_dx
  J_arr[1,2,] = part_y_dx
  J_arr[2,1,] = part_x_dy
  J_arr[2,2,] = part_y_dy

  J_arr

}

calculate_Gradient_f_wrt_beta = function(pos_t_mat, beta_mat, border){

  if(is.null(dim(pos_t_mat))) {
    pos_t_mat = matrix(pos_t_mat, nrow = 1, ncol = 3)
  }

  pos_mat = pos_t_mat[,2:3]

  M = sqrt(nrow(beta_mat))

  traj = exp(evaluate2DCosine_fast(beta_mat, pos_mat, border))
  log_traj_part_beta = evaluate2DCosine_part_beta_fast(beta_mat, pos_mat, border)
  VF = baseVectorFields_Vec(pos_t_mat)

  T1 = traj[,1]
  T2 = traj[,2]

  VF1_x = VF[,1]
  VF2_x = VF[,2]
  VF1_y = VF[,3]
  VF2_y = VF[,4]

  J_arr = array(dim = c(2*M^2,2, nrow(pos_t_mat)))

  for(i in 1:(M^2)){

    J_arr[i,1,] = VF1_x * T1 * log_traj_part_beta[,i]
    J_arr[i,2,] = VF1_y * T1 * log_traj_part_beta[,i]

    J_arr[i+(M^2),1,] = VF2_x * T2 * log_traj_part_beta[,i]
    J_arr[i+(M^2),2,] = VF2_y * T2 * log_traj_part_beta[,i]

  }

  J_arr

}

beta_mat_true

beta_mat = cbind(c(real_traj_test_1$Beta1Posterior[62,], rep(0,3)), c(real_traj_test_1$Beta2Posterior[62,], rep(0,3)))

calculate_part_path_part_beta = function(beta_mat, border, baseVectorFields_Vec, start_t_pos, end_t, N_prop_steps){

  M = sqrt(nrow(beta_mat))

  y_0 = as.numeric(start_t_pos[2:3])
  t_0 = as.numeric(start_t_pos[1])
  g_0 = matrix(rep(0, 2*(2*M^2)), nrow = 2)

  t_step <- (end_t - t_0) / N_prop_steps
  curPos = as.numeric(start_t_pos)

  pos_mat_k1 = matrix(rep(0, 3*N_prop_steps), nrow = N_prop_steps, ncol = 3)
  pos_mat_k2 = matrix(rep(0, 3*N_prop_steps), nrow = N_prop_steps, ncol = 3)
  pos_mat_k3 = matrix(rep(0, 3*N_prop_steps), nrow = N_prop_steps, ncol = 3)
  pos_mat_k4 = matrix(rep(0, 3*N_prop_steps), nrow = N_prop_steps, ncol = 3)

  t = rep(0,N_prop_steps+1)
  x = rep(0,N_prop_steps+1)
  y = rep(0,N_prop_steps+1)

  t[1] = curPos[1]
  x[1] = curPos[2]
  y[1] = curPos[3]

  for(j in 1:N_prop_steps) {

    pos_mat_k1[j,] = curPos

    k1 <- TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = curPos, beta_mat = beta_mat, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    pos_mat_k2[j,] <- curPos
    pos_mat_k2[j, 1] <- pos_mat_k2[j, 1] + t_step / 2
    pos_mat_k2[j, 2] <- pos_mat_k2[j, 2] + k1[, 1] * (t_step / 2)
    pos_mat_k2[j, 3] <- pos_mat_k2[j, 3] + k1[, 2] * (t_step / 2)

    k2 <- TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = pos_mat_k2[j,], beta_mat = beta_mat, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    pos_mat_k3[j,] <- curPos
    pos_mat_k3[j, 1] <- pos_mat_k3[j, 1] + t_step / 2
    pos_mat_k3[j, 2] <- pos_mat_k3[j, 2] + k2[, 1] * (t_step / 2)
    pos_mat_k3[j, 3] <- pos_mat_k3[j, 3] + k2[, 2] * (t_step / 2)

    k3 <- TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = pos_mat_k3[j,], beta_mat = beta_mat, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    pos_mat_k4[j,] <- curPos
    pos_mat_k4[j, 1] <- pos_mat_k4[j, 1] + t_step
    pos_mat_k4[j, 2] <- pos_mat_k4[j, 2] + k3[, 1] * t_step
    pos_mat_k4[j, 3] <- pos_mat_k4[j, 3] + k3[, 2] * t_step

    k4 <- TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = pos_mat_k4[j,], beta_mat = beta_mat, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    rk4_drift_1 <- (k1[, 1] + 2 * k2[, 1] + 2 * k3[, 1] + k4[, 1]) / 6
    rk4_drift_2 <- (k1[, 2] + 2 * k2[, 2] + 2 * k3[, 2] + k4[, 2]) / 6

    curPos[1] <- curPos[1] + t_step
    curPos[2] <- curPos[2] + rk4_drift_1 * t_step
    curPos[3] <- curPos[3] + rk4_drift_2 * t_step

    t[j+1] = curPos[1]
    x[j+1] = curPos[2]
    y[j+1] = curPos[3]
  }


  J_f_array_k1 = calculate_Jacobian_f_wrt_y(pos_t_mat = pos_mat_k1, beta_mat = beta_mat, border = border)
  J_f_array_k2 = calculate_Jacobian_f_wrt_y(pos_t_mat = pos_mat_k2, beta_mat = beta_mat, border = border)
  J_f_array_k3 = calculate_Jacobian_f_wrt_y(pos_t_mat = pos_mat_k3, beta_mat = beta_mat, border = border)
  J_f_array_k4 = calculate_Jacobian_f_wrt_y(pos_t_mat = pos_mat_k4, beta_mat = beta_mat, border = border)

  grad_f_beta_k1 = calculate_Gradient_f_wrt_beta(pos_t_mat = pos_mat_k1, beta_mat = beta_mat, border = border)
  grad_f_beta_k2 = calculate_Gradient_f_wrt_beta(pos_t_mat = pos_mat_k2, beta_mat = beta_mat, border = border)
  grad_f_beta_k3 = calculate_Gradient_f_wrt_beta(pos_t_mat = pos_mat_k3, beta_mat = beta_mat, border = border)
  grad_f_beta_k4 = calculate_Gradient_f_wrt_beta(pos_t_mat = pos_mat_k4, beta_mat = beta_mat, border = border)

  cur_g = g_0

  for(i in 1:N_prop_steps) {

    # ---------------------------------------------------------
    # RK4 Substep 1 (k1 for g)
    # ---------------------------------------------------------
    A1 <- J_f_array_k1[,,i]               # 2x2 matrix
    B1 <- t(grad_f_beta_k1[,,i])          # Transpose to 2xN matrix

    k1_g <- A1 %*% cur_g + B1             # Derivative at k1

    # ---------------------------------------------------------
    # RK4 Substep 2 (k2 for g)
    # ---------------------------------------------------------
    # Propagate g forward by half a step using k1_g
    g2 <- cur_g + k1_g * (t_step / 2)

    A2 <- J_f_array_k2[,,i]
    B2 <- t(grad_f_beta_k2[,,i])

    k2_g <- A2 %*% g2 + B2                # Derivative at k2

    # ---------------------------------------------------------
    # RK4 Substep 3 (k3 for g)
    # ---------------------------------------------------------
    # Propagate g forward by half a step using k2_g
    g3 <- cur_g + k2_g * (t_step / 2)

    A3 <- J_f_array_k3[,,i]
    B3 <- t(grad_f_beta_k3[,,i])

    k3_g <- A3 %*% g3 + B3                # Derivative at k3

    # ---------------------------------------------------------
    # RK4 Substep 4 (k4 for g)
    # ---------------------------------------------------------
    # Propagate g forward by a FULL step using k3_g
    g4 <- cur_g + k3_g * t_step

    A4 <- J_f_array_k4[,,i]
    B4 <- t(grad_f_beta_k4[,,i])

    k4_g <- A4 %*% g4 + B4                # Derivative at k4

    # ---------------------------------------------------------
    # Final RK4 Update for cur_g
    # ---------------------------------------------------------
    # Now we properly combine them and ADD to the previous cur_g
    cur_g <- cur_g + (t_step / 6) * (k1_g + 2*k2_g + 2*k3_g + k4_g)

  }

  return(list(Del_G = cur_g, PropPos = curPos))

}

start_t_pos_mat

beta_mat = cbind(c(real_traj_test_1$Beta1Posterior[62,], rep(0,3)), c(real_traj_test_1$Beta2Posterior[62,], rep(0,3)))
beta_mat = cbind(c(0,0,0,0), c(0,0,0,0))

beta_mat = beta_mat_true

calculate_importance_weight = function(beta_mat, border, baseVectorFields_Vec, start_t_pos_mat, end_t_pos_mat, N_prop_steps, pos_sd, prior_beta_sigma){

  N_data = nrow(start_t_pos_mat)

  del_G_list = list()
  end_t_pos_prop_mat = matrix(0, nrow = N_data, ncol = 2)

  for(i in 1:N_data){

    cur_del_path_del_beta = calculate_part_path_part_beta(beta_mat, border, baseVectorFields_Vec, start_t_pos = start_t_pos_mat[i,], end_t = end_t_pos_mat[i,1], N_prop_steps)

    del_G_list[[i]] = cur_del_path_del_beta$Del_G
    end_t_pos_prop_mat[i,] = cur_del_path_del_beta$PropPos[2:3]

  }

  del_G = do.call(rbind, del_G_list)
  end_t_pos_true = c(t(end_t_pos_mat[,2:3]))
  end_t_pos_prop = c(t(end_t_pos_prop_mat))
  beta_draw = c(beta_mat)

  var_inv <- 1 / pos_sd^2
  C <- blkdiag(prior_beta_sigma, prior_beta_sigma)
  L_inv = diag(var_inv, 2 * N_data)

  K <- ((end_t_pos_prop - end_t_pos_true) + del_G %*% beta_draw) * var_inv
  H <- solve(L_inv + (del_G %*% C %*% t(del_G)) * (var_inv^2))

  J_approx <- diag(length(beta_mat)) + (C %*% t(del_G) %*% del_G) * var_inv

  # Compute log(|J|) safely
  log_det_J <- as.numeric(determinant(J_approx, logarithm = TRUE)$modulus)

  # Compute log(importance)
  log_importance <- -0.5 * as.numeric(t(K) %*% H %*% K) - 0.5 * log_det_J

  return(log_importance)

}

t1 = Sys.time()

log_importance = calculate_importance_weight(beta_mat_fake, border, baseVectorFields_Vec, start_t_pos_mat, end_t_pos_mat, N_prop_steps, pos_sd, prior_beta_sigma)

t2 = Sys.time()

log_weight <- calculate_log_importance_cpp(
  beta = c(beta_mat_true), # Pass the flattened vector
  M_sq = nrow(beta_mat_true),
  start_t_pos_mat = as.matrix(start_t_pos_mat),
  end_t_pos_true_mat = as.matrix(end_t_pos_mat),
  t_steps = end_t_pos_mat[,1] - start_t_pos_mat[,1], # Supply dt array directly
  N_prop_steps = N_prop_steps,
  border = border,
  pos_sd = pos_sd,
  prior_beta_sigma = prior_precision_mat,
  n_threads = 32 # Set to your 32-core capacity!
)


find_one_rMAP_Trajectory = function(sim_data_list, pos_sd = 0.001, vel_sd = 0.1, M = 2, prior_k = 0.35, prior_l = 1, baseVectorFields_Vec, border, N_prop_steps = 10, traj_eval_grid, print_every = 50){

  N_v = as.numeric(lapply(sim_data_list, nrow))
  N = sum(N_v)
  D = length(N_v)

  # Add positional error
  t_err = lapply(N_v, FUN = rep, x = 0)
  x_pos_err = lapply(X = N_v, FUN = rnorm, mean = 0, sd = pos_sd)
  y_pos_err = lapply(X = N_v, FUN = rnorm, mean = 0, sd = pos_sd)

  err_list = lapply(1:D, FUN = function(i, l1, l2, l3){cbind(l1[[i]],l2[[i]],l3[[i]])}, l1 = t_err, l2 = x_pos_err, l3 = y_pos_err)

  aug_data_list = lapply(1:D, FUN = function(i, l1, l2){l1[[i]] + l2[[i]]}, l1 = sim_data_list, l2 = err_list)

  aug_data_starts = do.call(rbind, lapply(1:D, FUN = function(i, l){l[[i]][1:(N_v[i]-1),]}, l = aug_data_list))
  aug_data_ends = do.call(rbind, lapply(1:D, FUN = function(i, l){l[[i]][-1,]}, l = aug_data_list))

  N_advects = nrow(aug_data_starts)

  rand_vel_1 = matrix(rnorm(N_advects*N_prop_steps, 0, vel_sd), nrow = N_prop_steps, ncol = N_advects)
  rand_vel_2 = matrix(rnorm(N_advects*N_prop_steps, 0, vel_sd), nrow = N_prop_steps, ncol = N_advects)

  omega = (1:M-1)*pi
  spec_den = sqrt(2*pi)*prior_l*exp(-0.5*prior_l^2*omega^2)

  if(M == 1){

    prior_beta_sigma = prior_k^2 * spec_den^2

  } else{

    prior_beta_sigma = diag(c(prior_k^2 * diag((spec_den)) %*% matrix(rep(1,M^2), nrow = M) %*% diag((spec_den))))

  }

  start_beta_1 = mvrnorm(1, rep(0, M^2), Sigma = prior_beta_sigma)
  start_beta_2 = mvrnorm(1, rep(0, M^2), Sigma = prior_beta_sigma)

  start_beta = c(start_beta_1, start_beta_2)

  t_steps_vec <- as.numeric(aug_data_ends$t - aug_data_starts$t) / N_prop_steps

  aug_starts_mat <- as.matrix(aug_data_starts)
  aug_ends_mat <- as.matrix(aug_data_ends)
  rand_vel_1_mat <- as.matrix(rand_vel_1)
  rand_vel_2_mat <- as.matrix(rand_vel_2)
  border_vec <- as.numeric(border)
  start_beta_vec <- as.numeric(start_beta)

  # Good catch on the prior! Pre-calculate the precision matrix (inverse of covariance) here
  prior_precision_mat <- ginv(as.matrix(prior_beta_sigma))

  eval_counter <- 0

  cat(sprintf("\n\n--- Starting Optimization for New Sample ---\n"))

  # Optimize for the best w that matched the mixing fields
  rMAP_loss <- function(beta) {

    eval_counter <<- eval_counter + 1

    # The only thing happening here is the jump to C++
    current_loss <- rMAP_loss_cpp(
      beta = beta,
      M_sq = M^2,
      aug_data_starts = aug_starts_mat,
      aug_data_ends = aug_ends_mat,
      t_steps = t_steps_vec,
      rand_vel_1 = rand_vel_1_mat,
      rand_vel_2 = rand_vel_2_mat,
      border = border_vec,
      pos_sd = pos_sd,
      prior_beta_sigma = prior_precision_mat,
      start_beta = start_beta_vec
    )

    if (eval_counter %% print_every == 0) {
      beta_str <- paste(sprintf("%.3f", beta), collapse = ", ")

      # Notice the \n at the VERY BEGINNING here too, to clear the progress bar
      cat(sprintf("\nIter: %4d | Loss: %10.4f | Beta: [%s]",
                  eval_counter, current_loss, beta_str))
    }

    return(current_loss)

  }

    # Run the optimizer
  opt_result <- nloptr(
    x0 = start_beta_vec,               # Starting values
    eval_f = rMAP_loss,            # Your C++ wrapper function
    opts = list(
      "algorithm" = "NLOPT_LN_NEWUOA",  # LN = Local, No-derivative
      "ftol_rel" = 1e-6,                # Stop when parameters stop changing by this fraction
      "maxeval" = 2000,                 # Maximum number of evaluations
      "print_level" = 0                 # 0 = silent, 1 = show progress, 2 = verbose
    )
  )

  final_beta_str <- paste(sprintf("%.3f", opt_result$solution), collapse = ", ")

  # Using "FINAL" so it stands out from the regular 100-step updates
  cat(sprintf("\nFINAL: Iter: %4d | Loss: %10.4f | Beta: [%s]\n",
              opt_result$iterations, opt_result$objective, final_beta_str))

  post_beta_1 = opt_result$solution[1:(M^2)]
  post_beta_2 = opt_result$solution[1:(M^2)+(M^2)]
  post_beta_mat = cbind(post_beta_1, post_beta_2)

  post_traj_grid = exp(evaluate2DCosine_fast(beta_mat = post_beta_mat, pos_mat = traj_eval_grid, border = border))

  # --- Eval Optimal Traj (Updated to RK4) ---

  M_sq <- M^2
  t_steps <- (aug_data_ends$t - aug_data_starts$t) / N_prop_steps
  sqrt_t_steps <- sqrt(t_steps)
  curPos_mat <- aug_data_starts

  for(j in 1:N_prop_steps) {

    k1 <- TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = curPos_mat, beta_mat = post_beta_mat, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    pos_mat_k2 <- curPos_mat
    pos_mat_k2[, 1] <- pos_mat_k2[, 1] + t_steps / 2
    pos_mat_k2[, 2] <- pos_mat_k2[, 2] + k1[, 1] * (t_steps / 2)
    pos_mat_k2[, 3] <- pos_mat_k2[, 3] + k1[, 2] * (t_steps / 2)

    k2 <- TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = pos_mat_k2, beta_mat = post_beta_mat, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    pos_mat_k3 <- curPos_mat
    pos_mat_k3[, 1] <- pos_mat_k3[, 1] + t_steps / 2
    pos_mat_k3[, 2] <- pos_mat_k3[, 2] + k2[, 1] * (t_steps / 2)
    pos_mat_k3[, 3] <- pos_mat_k3[, 3] + k2[, 2] * (t_steps / 2)

    k3 <- TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = pos_mat_k3, beta_mat = post_beta_mat, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    pos_mat_k4 <- curPos_mat
    pos_mat_k4[, 1] <- pos_mat_k4[, 1] + t_steps
    pos_mat_k4[, 2] <- pos_mat_k4[, 2] + k3[, 1] * t_steps
    pos_mat_k4[, 3] <- pos_mat_k4[, 3] + k3[, 2] * t_steps

    k4 <- TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = pos_mat_k4, beta_mat = post_beta_mat, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    rk4_drift_1 <- (k1[, 1] + 2 * k2[, 1] + 2 * k3[, 1] + k4[, 1]) / 6
    rk4_drift_2 <- (k1[, 2] + 2 * k2[, 2] + 2 * k3[, 2] + k4[, 2]) / 6

    cur_diffusion_vec_1 <- rand_vel_1[j,]
    cur_diffusion_vec_2 <- rand_vel_2[j,]

    curPos_mat[, 1] <- curPos_mat[, 1] + t_steps
    curPos_mat[, 2] <- curPos_mat[, 2] + rk4_drift_1 * t_steps + cur_diffusion_vec_1 * sqrt_t_steps
    curPos_mat[, 3] <- curPos_mat[, 3] + rk4_drift_2 * t_steps + cur_diffusion_vec_2 * sqrt_t_steps
  }

  diff_1 <- curPos_mat[, 2] - aug_data_ends$X1
  diff_2 <- curPos_mat[, 3] - aug_data_ends$X2

  post_likelihood_loss_pos <- sum(diff_1^2 + diff_2^2) / (2 * pos_sd^2)

  d_beta1 <- post_beta_1 - start_beta_1
  d_beta2 <- post_beta_2 - start_beta_2

  post_prior_loss <- crossprod(d_beta1, prior_precision_mat %*% d_beta1) +
    crossprod(d_beta2, prior_precision_mat %*% d_beta2)

  cat("\n")

  list(post_beta_mat, post_likelihood_loss_pos, post_prior_loss, post_traj_grid)

}

find_one_rMAP_Trajectory_Multithread = function(sim_data_list, pos_sd = 0.001, vel_sd = 0.1, M = 2, prior_k = 0.35, prior_l = 1, baseVectorFields_Vec, border, N_prop_steps = 10, traj_eval_grid, print_every = 50, n_threads = 1){

  N_v = as.numeric(lapply(sim_data_list, nrow))
  N = sum(N_v)
  D = length(N_v)

  # Add positional error
  t_err = lapply(N_v, FUN = rep, x = 0)
  x_pos_err = lapply(X = N_v, FUN = rnorm, mean = 0, sd = pos_sd)
  y_pos_err = lapply(X = N_v, FUN = rnorm, mean = 0, sd = pos_sd)

  err_list = lapply(1:D, FUN = function(i, l1, l2, l3){cbind(l1[[i]],l2[[i]],l3[[i]])}, l1 = t_err, l2 = x_pos_err, l3 = y_pos_err)

  aug_data_list = lapply(1:D, FUN = function(i, l1, l2){l1[[i]] + l2[[i]]}, l1 = sim_data_list, l2 = err_list)

  aug_data_starts = do.call(rbind, lapply(1:D, FUN = function(i, l){l[[i]][1:(N_v[i]-1),]}, l = aug_data_list))
  aug_data_ends = do.call(rbind, lapply(1:D, FUN = function(i, l){l[[i]][-1,]}, l = aug_data_list))

  data_starts = do.call(rbind, lapply(1:D, FUN = function(i, l){l[[i]][1:(N_v[i]-1),]}, l = sim_data_list))
  data_ends = do.call(rbind, lapply(1:D, FUN = function(i, l){l[[i]][-1,]}, l = sim_data_list))

  N_advects = nrow(aug_data_starts)

  rand_vel_1 = matrix(rnorm(N_advects*N_prop_steps, 0, vel_sd), nrow = N_prop_steps, ncol = N_advects)
  rand_vel_2 = matrix(rnorm(N_advects*N_prop_steps, 0, vel_sd), nrow = N_prop_steps, ncol = N_advects)

  omega = (1:M-1)*pi
  spec_den = sqrt(2*pi)*prior_l*exp(-0.5*prior_l^2*omega^2)

  if(M == 1){

    prior_beta_sigma = prior_k^2 * spec_den^2

  } else{

    prior_beta_sigma = diag(c(prior_k^2 * diag((spec_den)) %*% matrix(rep(1,M^2), nrow = M) %*% diag((spec_den))))

  }

  start_beta_1 = mvrnorm(1, rep(0, M^2), Sigma = prior_beta_sigma)
  start_beta_2 = mvrnorm(1, rep(0, M^2), Sigma = prior_beta_sigma)

  start_beta = c(start_beta_1, start_beta_2)

  t_steps_vec <- as.numeric(aug_data_ends$t - aug_data_starts$t) / N_prop_steps


  data_starts_mat = as.matrix(data_starts)
  data_ends_mat = as.matrix(data_ends)
  aug_starts_mat <- as.matrix(aug_data_starts)
  aug_ends_mat <- as.matrix(aug_data_ends)
  rand_vel_1_mat <- as.matrix(rand_vel_1)
  rand_vel_2_mat <- as.matrix(rand_vel_2)
  border_vec <- as.numeric(border)
  start_beta_vec <- as.numeric(start_beta)

  # Good catch on the prior! Pre-calculate the precision matrix (inverse of covariance) here
  prior_precision_mat <- ginv(as.matrix(prior_beta_sigma))

  eval_counter <- 0

  cat(sprintf("\n\n--- Starting Optimization for New Sample ---\n"))

  # Optimize for the best w that matched the mixing fields
  rMAP_loss <- function(beta) {

    eval_counter <<- eval_counter + 1

    # The only thing happening here is the jump to C++
    current_loss <- rMAP_loss_cpp_multi(
      beta = beta,
      M_sq = M^2,
      aug_data_starts = aug_starts_mat,
      aug_data_ends = aug_ends_mat,
      t_steps = t_steps_vec,
      rand_vel_1 = rand_vel_1_mat,
      rand_vel_2 = rand_vel_2_mat,
      border = border_vec,
      pos_sd = pos_sd,
      prior_beta_sigma = prior_precision_mat,
      start_beta = start_beta_vec,
      n_threads = n_threads
    )

    if (eval_counter %% print_every == 0) {
      beta_str <- paste(sprintf("%.3f", beta), collapse = ", ")

      # Notice the \n at the VERY BEGINNING here too, to clear the progress bar
      cat(sprintf("\nIter: %4d | Loss: %10.4f | Beta: [%s]",
                  eval_counter, current_loss, beta_str))
    }

    return(current_loss)

  }

  # Run the optimizer
  opt_result <- nloptr(
    x0 = start_beta_vec,               # Starting values
    eval_f = rMAP_loss,            # Your C++ wrapper function
    opts = list(
      "algorithm" = "NLOPT_LN_NEWUOA",  # LN = Local, No-derivative
      "ftol_rel" = 1e-6,                # Stop when parameters stop changing by this fraction
      "maxeval" = 5000,                 # Maximum number of evaluations
      "print_level" = 0                 # 0 = silent, 1 = show progress, 2 = verbose
    )
  )

  final_beta_str <- paste(sprintf("%.3f", opt_result$solution), collapse = ", ")

  post_beta_1 = opt_result$solution[1:(M^2)]
  post_beta_2 = opt_result$solution[1:(M^2)+(M^2)]
  post_beta_mat = cbind(post_beta_1, post_beta_2)

  post_traj_grid = exp(evaluate2DCosine_fast(beta_mat = post_beta_mat, pos_mat = traj_eval_grid, border = border))

  # --- Eval Optimal Traj (Updated to RK4) ---

  post_loss = rMAP_loss_cpp_multi(
    beta = opt_result$solution,
    M_sq = M^2,
    aug_data_starts = aug_starts_mat,
    aug_data_ends = aug_ends_mat,
    t_steps = t_steps_vec,
    rand_vel_1 = rand_vel_1_mat,
    rand_vel_2 = rand_vel_2_mat,
    border = border_vec,
    pos_sd = pos_sd,
    prior_beta_sigma = prior_precision_mat,
    start_beta = start_beta_vec,
    n_threads = n_threads
  )

  log_importance = calculate_log_importance_cpp(
    beta = opt_result$solution,
    M_sq = M^2,
    start_t_pos_mat = aug_starts_mat,
    end_t_pos_true_mat = data_ends_mat,
    t_steps = data_ends_mat[,1] - data_starts_mat[,1],
    N_prop_steps = N_prop_steps,
    border = border,
    pos_sd = pos_sd,
    prior_beta_sigma = prior_precision_mat,
    n_threads = n_threads
  )

  # Using "FINAL" so it stands out from the regular 100-step updates
  cat(sprintf("\nFINAL: Iter: %4d | Loss: %10.4f | Beta: [%s] | Log Importance: %.4f\n",
              opt_result$iterations, opt_result$objective, final_beta_str, log_importance))

  cat("\n")

  list(post_beta_mat, post_loss, post_traj_grid, log_importance)

}

run_rMAP_Trajectory = function(N_samples, sim_data_list, pos_sd, vel_sd, M, prior_k, prior_l, baseVectorFields_Vec, border, N_prop_steps, traj_eval_grid, print_every = 50){

  beta_1_post_mat = matrix(nrow = N_samples, ncol = M^2)
  beta_2_post_mat = matrix(nrow = N_samples, ncol = M^2)

  post_traj_eval_1_list = list()
  post_traj_eval_2_list = list()

  post_like_pos = rep(0, N_samples)
  post_like_prior = rep(0, N_samples)

  for(i in 1:N_samples){

    cat(sprintf("=========== PROGRESS: Sample %d of %d ===========\n", i, N_samples))
    flush.console()

    cur_draw = find_one_rMAP_Trajectory(sim_data_list = sim_data_list, pos_sd = pos_sd, vel_sd = vel_sd, M = M, prior_k = prior_k, prior_l = prior_l, baseVectorFields_Vec = baseVectorFields_Vec, border = border, N_prop_steps = N_prop_steps, traj_eval_grid = traj_eval_grid, print_every = print_every)

    beta_1_post_mat[i,] = cur_draw[[1]][,1]
    beta_2_post_mat[i,] = cur_draw[[1]][,2]

    post_like_pos[i] = cur_draw[[2]]
    post_like_prior[i] = cur_draw[[3]]

    post_traj_eval_1_list[[i]] = cur_draw[[4]][,1]
    post_traj_eval_2_list[[i]] = cur_draw[[4]][,2]


  }

  post_draws_traj_eval_1 = do.call(cbind, post_traj_eval_1_list)
  post_draws_traj_eval_2 = do.call(cbind, post_traj_eval_2_list)

  list(Beta1Posterior = beta_1_post_mat, Beta2Posterior = beta_2_post_mat,
       Traj1Posterior = post_draws_traj_eval_1, Traj2Posterior = post_draws_traj_eval_2,
       PosPosteriorNLL = post_like_pos, PriorPosteriorNLL = post_like_prior)

}

run_rMAP_Trajectory_Multithread = function(N_samples, sim_data_list, pos_sd, vel_sd, M, prior_k, prior_l, baseVectorFields_Vec, border, N_prop_steps, traj_eval_grid, print_every = 50, n_threads = 1){

  beta_1_post_mat = matrix(nrow = N_samples, ncol = M^2)
  beta_2_post_mat = matrix(nrow = N_samples, ncol = M^2)

  post_traj_eval_1_list = list()
  post_traj_eval_2_list = list()

  post_like = rep(0, N_samples)
  log_importance = rep(0, N_samples)

  for(i in 1:N_samples){

    cat(sprintf("=========== PROGRESS: Sample %d of %d ===========\n", i, N_samples))
    flush.console()

    cur_draw = find_one_rMAP_Trajectory_Multithread(sim_data_list = sim_data_list, pos_sd = pos_sd, vel_sd = vel_sd, M = M, prior_k = prior_k, prior_l = prior_l, baseVectorFields_Vec = baseVectorFields_Vec, border = border, N_prop_steps = N_prop_steps, traj_eval_grid = traj_eval_grid, print_every = print_every, n_threads = n_threads)

    beta_1_post_mat[i,] = cur_draw[[1]][,1]
    beta_2_post_mat[i,] = cur_draw[[1]][,2]

    post_like[i] = cur_draw[[2]]
    log_importance[i] = cur_draw[[4]]

    post_traj_eval_1_list[[i]] = cur_draw[[3]][,1]
    post_traj_eval_2_list[[i]] = cur_draw[[3]][,2]

  }

  post_draws_traj_eval_1 = do.call(cbind, post_traj_eval_1_list)
  post_draws_traj_eval_2 = do.call(cbind, post_traj_eval_2_list)

  list(Beta1Posterior = beta_1_post_mat, Beta2Posterior = beta_2_post_mat,
       Traj1Posterior = post_draws_traj_eval_1, Traj2Posterior = post_draws_traj_eval_2,
       PosteriorNLL = post_like,
       LogImportanceWeights = log_importance)

}


get_starting_mode_trajectories = function(traj_samples, PC_importance_threshold, max_clusters){

  N_coefs = ncol(traj_samples$Beta1Posterior)

  traj_coef_pca = prcomp(cbind(traj_samples$Beta1Posterior, traj_samples$Beta2Posterior), center = T, scale. = F)
  pca_summary = summary(traj_coef_pca)

  N_PCs = as.numeric(which(pca_summary$importance[3,] > 0.99)[1])

  PCA_sub = traj_coef_pca$x[,1:(N_PCs)]
  W_pca = traj_coef_pca$rotation[,1:(N_PCs)]
  mu_pca = traj_coef_pca$center

  k_means = run_k_means(PCA_sub = PCA_sub, max_k = max_clusters)

  K_means_centers = k_means$K_means$centers
  K_means_K = k_means$K

  PCA_center_beta = K_means_centers %*% t(W_pca) + matrix(rep(mu_pca, K_means_K), nrow = K_means_K, byrow = T)

  PCA_center_beta1 = PCA_center_beta[,1:N_coefs]
  PCA_center_beta2 = PCA_center_beta[,1:N_coefs + N_coefs]

  list(Center_Beta1 = PCA_center_beta1, Center_Beta2 = PCA_center_beta2)

}

find_one_mode_location_traj = function(center_beta_start_1, center_beta_start_2, sim_data_list, pos_sd = 0.001, vel_sd = 0.1, M = 2, prior_k = 0.35, prior_l = 1, baseVectorFields_Vec, border, N_prop_steps = 10, traj_eval_grid, print_every = 50){

  center_beta_start_full = c(center_beta_start_1, center_beta_start_2)

  N_v = as.numeric(lapply(sim_data_list, nrow))
  N = sum(N_v)
  D = length(N_v)

  data_starts = do.call(rbind, lapply(1:D, FUN = function(i, l){l[[i]][1:(N_v[i]-1),]}, l = sim_data_list))
  data_ends = do.call(rbind, lapply(1:D, FUN = function(i, l){l[[i]][-1,]}, l = sim_data_list))

  N_advects = nrow(data_starts)

  rand_vel_1 = matrix(rep(0, N_advects*N_prop_steps), nrow = N_prop_steps, ncol = N_advects)
  rand_vel_2 = matrix(rep(0, N_advects*N_prop_steps), nrow = N_prop_steps, ncol = N_advects)

  omega = (1:M-1)*pi
  spec_den = sqrt(2*pi)*prior_l*exp(-0.5*prior_l^2*omega^2)

  prior_beta_sigma = diag(c(prior_k^2 * diag(spec_den) %*% matrix(rep(1,4), nrow = 2) %*% diag(spec_den)))

  start_beta = c(rep(0,2*M^2))

  t_steps_vec <- as.numeric(data_ends$t - data_starts$t) / N_prop_steps

  starts_mat <- as.matrix(data_starts)
  ends_mat <- as.matrix(data_ends)
  rand_vel_1_mat <- as.matrix(rand_vel_1)
  rand_vel_2_mat <- as.matrix(rand_vel_2)
  border_vec <- as.numeric(border)
  start_beta_vec <- as.numeric(start_beta)

  # Good catch on the prior! Pre-calculate the precision matrix (inverse of covariance) here
  prior_precision_mat <- ginv(as.matrix(prior_beta_sigma))

  eval_counter <- 0

  cat(sprintf("\n\n--- Starting Optimization for New Sample ---\n"))

  # Optimize for the best w that matched the mixing fields
  rMAP_loss <- function(beta) {

    eval_counter <<- eval_counter + 1

    # The only thing happening here is the jump to C++
    current_loss <- rMAP_loss_cpp(
      beta = beta,
      M_sq = M^2,
      aug_data_starts = starts_mat,
      aug_data_ends = ends_mat,
      t_steps = t_steps_vec,
      rand_vel_1 = rand_vel_1_mat,
      rand_vel_2 = rand_vel_2_mat,
      border = border_vec,
      pos_sd = pos_sd,
      prior_beta_sigma = prior_precision_mat,
      start_beta = start_beta_vec
    )

    if (eval_counter %% print_every == 0) {
      beta_str <- paste(sprintf("%.3f", beta), collapse = ", ")

      # Notice the \n at the VERY BEGINNING here too, to clear the progress bar
      cat(sprintf("\nIter: %4d | Loss: %10.4f | Beta: [%s]",
                  eval_counter, current_loss, beta_str))
    }

    return(current_loss)

  }

  # Run the optimizer
  opt_result <- nloptr(
    x0 = center_beta_start_full,               # Starting values
    eval_f = rMAP_loss,            # Your C++ wrapper function
    opts = list(
      "algorithm" = "NLOPT_LN_NEWUOA",  # LN = Local, No-derivative
      "ftol_rel" = 1e-6,                # Stop when parameters stop changing by this fraction
      "maxeval" = 2000,                 # Maximum number of evaluations
      "print_level" = 0                 # 0 = silent, 1 = show progress, 2 = verbose
    )
  )

  final_beta_str <- paste(sprintf("%.3f", opt_result$solution), collapse = ", ")

  # Using "FINAL" so it stands out from the regular 100-step updates
  cat(sprintf("\nFINAL: Iter: %4d | Loss: %10.4f | Beta: [%s]\n",
              opt_result$iterations, opt_result$objective, final_beta_str))

  post_beta_1 = opt_result$solution[1:(M^2)]
  post_beta_2 = opt_result$solution[1:(M^2)+(M^2)]
  post_beta_mat = cbind(post_beta_1, post_beta_2)

  post_traj_grid = exp(evaluate2DCosine_fast(beta_mat = post_beta_mat, pos_mat = traj_eval_grid, border = border))

  # --- Eval Optimal Traj (Updated to RK4) ---

  M_sq <- M^2
  t_steps <- (data_ends$t - data_starts$t) / N_prop_steps
  sqrt_t_steps <- sqrt(t_steps)
  curPos_mat <- data_starts

  for(j in 1:N_prop_steps) {

    k1 <- TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = curPos_mat, beta_mat = post_beta_mat, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    pos_mat_k2 <- curPos_mat
    pos_mat_k2[, 1] <- pos_mat_k2[, 1] + t_steps / 2
    pos_mat_k2[, 2] <- pos_mat_k2[, 2] + k1[, 1] * (t_steps / 2)
    pos_mat_k2[, 3] <- pos_mat_k2[, 3] + k1[, 2] * (t_steps / 2)

    k2 <- TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = pos_mat_k2, beta_mat = post_beta_mat, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    pos_mat_k3 <- curPos_mat
    pos_mat_k3[, 1] <- pos_mat_k3[, 1] + t_steps / 2
    pos_mat_k3[, 2] <- pos_mat_k3[, 2] + k2[, 1] * (t_steps / 2)
    pos_mat_k3[, 3] <- pos_mat_k3[, 3] + k2[, 2] * (t_steps / 2)

    k3 <- TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = pos_mat_k3, beta_mat = post_beta_mat, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    pos_mat_k4 <- curPos_mat
    pos_mat_k4[, 1] <- pos_mat_k4[, 1] + t_steps
    pos_mat_k4[, 2] <- pos_mat_k4[, 2] + k3[, 1] * t_steps
    pos_mat_k4[, 3] <- pos_mat_k4[, 3] + k3[, 2] * t_steps

    k4 <- TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = pos_mat_k4, beta_mat = post_beta_mat, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    rk4_drift_1 <- (k1[, 1] + 2 * k2[, 1] + 2 * k3[, 1] + k4[, 1]) / 6
    rk4_drift_2 <- (k1[, 2] + 2 * k2[, 2] + 2 * k3[, 2] + k4[, 2]) / 6

    cur_diffusion_vec_1 <- rand_vel_1[j,]
    cur_diffusion_vec_2 <- rand_vel_2[j,]

    curPos_mat[, 1] <- curPos_mat[, 1] + t_steps
    curPos_mat[, 2] <- curPos_mat[, 2] + rk4_drift_1 * t_steps + cur_diffusion_vec_1 * sqrt_t_steps
    curPos_mat[, 3] <- curPos_mat[, 3] + rk4_drift_2 * t_steps + cur_diffusion_vec_2 * sqrt_t_steps
  }

  diff_1 <- curPos_mat[, 2] - data_ends$X1
  diff_2 <- curPos_mat[, 3] - data_ends$X2

  post_likelihood_loss_pos <- sum(diff_1^2 + diff_2^2) / (2 * pos_sd^2)

  d_beta1 <- post_beta_1
  d_beta2 <- post_beta_2

  post_prior_loss <- as.numeric(crossprod(d_beta1, prior_precision_mat %*% d_beta1) +
    crossprod(d_beta2, prior_precision_mat %*% d_beta2))

  cat("\n")

  list(OptimizedBeta = post_beta_mat, PostNLL_Pos = post_likelihood_loss_pos, PostNLL_Prior = post_prior_loss, OptimizedTrajVals = post_traj_grid)

}

get_posterior_traj_modes = function(center_beta_mat_start_1, center_beta_mat_start_2, sim_data_list, pos_sd = 0.001, vel_sd = 0.1, M = 2, prior_k = 0.35, prior_l = 1, baseVectorFields_Vec, border, N_prop_steps = 10, traj_eval_grid, print_every = 50, unique_tol = 0.001){

  # Define matrices for mode coefficients, mode path evaluations, NLL evaluations
  N_traj_modes = nrow(center_beta_mat_start_1)
  N_beta = ncol(center_beta_mat_start_1)
  N_evals = nrow(traj_eval_grid)

  traj_mode_beta1 = matrix(nrow = N_traj_modes, ncol = N_beta)
  traj_mode_beta2 = matrix(nrow = N_traj_modes, ncol = N_beta)

  traj_mode_eval1 = matrix(nrow = N_traj_modes, ncol = N_evals)
  traj_mode_eval2 = matrix(nrow = N_traj_modes, ncol = N_evals)

  traj_mode_NLL_pos = rep(0, N_traj_modes)
  traj_mode_NLL_prior = rep(0, N_traj_modes)

  # Tracker for the number of disjoint modes discovered
  unique_count = 0

  # Loops through all mode starting locations
  for(i in 1:N_traj_modes){

    cur_traj_start_1 = center_beta_mat_start_1[i,]
    cur_traj_start_2 = center_beta_mat_start_2[i,]

    cur_trajMode = find_one_mode_location_traj(
      center_beta_start_1 = cur_traj_start_1,
      center_beta_start_2 = cur_traj_start_2,
      sim_data_list = sim_data_list, pos_sd = pos_sd, vel_sd = vel_sd,
      M = M, prior_k = prior_k, prior_l = prior_l,
      baseVectorFields_Vec = baseVectorFields_Vec, border = border,
      N_prop_steps = N_prop_steps, traj_eval_grid = traj_eval_grid,
      print_every = print_every
    )

    # Extract the optimized joint parameter vector
    opt_beta1 = cur_trajMode$OptimizedBeta[,1]
    opt_beta2 = cur_trajMode$OptimizedBeta[,2]
    opt_joint = c(opt_beta1, opt_beta2)

    # -------------------------------------------------------------
    # Deduplication Filter
    # -------------------------------------------------------------
    is_duplicate = FALSE

    if(unique_count > 0) {
      for(j in 1:unique_count) {
        # Reconstruct the joint vector of previously accepted unique modes
        saved_joint = c(traj_mode_beta1[j,], traj_mode_beta2[j,])

        # Calculate Euclidean distance in the joint parameter space
        dist = sqrt(sum((opt_joint - saved_joint)^2))

        if (dist < unique_tol) {
          is_duplicate = TRUE
          break # Immediately halt check if it fell into a known valley
        }
      }
    }

    # Store only if it represents a novel physical route
    if(!is_duplicate) {
      unique_count = unique_count + 1

      traj_mode_beta1[unique_count,] = opt_beta1
      traj_mode_beta2[unique_count,] = opt_beta2

      traj_mode_eval1[unique_count,] = cur_trajMode$OptimizedTrajVals[,1]
      traj_mode_eval2[unique_count,] = cur_trajMode$OptimizedTrajVals[,2]

      traj_mode_NLL_pos[unique_count] = cur_trajMode$PostNLL_Pos
      traj_mode_NLL_prior[unique_count] = cur_trajMode$PostNLL_Prior
    }
  }

  # Trim pre-allocated structures down to only the unique modes found
  # drop = FALSE ensures matrices don't collapse to vectors if unique_count == 1
  if (unique_count > 0) {
    traj_mode_beta1 = traj_mode_beta1[1:unique_count, , drop = FALSE]
    traj_mode_beta2 = traj_mode_beta2[1:unique_count, , drop = FALSE]
    traj_mode_eval1 = traj_mode_eval1[1:unique_count, , drop = FALSE]
    traj_mode_eval2 = traj_mode_eval2[1:unique_count, , drop = FALSE]
    traj_mode_NLL_pos = traj_mode_NLL_pos[1:unique_count]
    traj_mode_NLL_prior = traj_mode_NLL_prior[1:unique_count]
  }

  return(list(
    Beta1 = traj_mode_beta1,
    Beta2 = traj_mode_beta2,
    TrajVals1 = traj_mode_eval1,
    TrajVals2 = traj_mode_eval2,
    NLL_Pos = traj_mode_NLL_pos,
    NLL_Prior = traj_mode_NLL_prior,
    Total_Unique_Modes = unique_count
  ))
}

test_posterior_traj_modes = get_posterior_traj_modes(center_beta_mat_start_1 = center_beta_mat_start_1, center_beta_mat_start_2 = center_beta_mat_start_2, sim_data_list = sim_data_list, pos_sd = 0.001, vel_sd = 0.1, M = 2, prior_k = 0.35, prior_l = 1, baseVectorFields_Vec = baseVectorFields_Vec, border = c(-2,-2,2,2), N_prop_steps = 10, traj_eval_grid = traj_eval_grid, print_every = 50, unique_tol = 0.001)

posterior_traj_mode_NLL = test_posterior_traj_modes$NLL_Pos + test_posterior_traj_modes$NLL_Prior

posterior_traj_mode_1 = test_posterior_traj_modes$Beta1[1,]
posterior_traj_mode_2 = test_posterior_traj_modes$Beta2[1,]

get_one_traj_hessian = function(posterior_traj_mode_1, posterior_traj_mode_2, sim_data_list, pos_sd = 0.001, vel_sd = 0.1, M = 2, prior_k = 0.35, prior_l = 1, baseVectorFields_Vec, border, N_prop_steps = 10, traj_eval_grid, print_every = 50){

  posterior_traj_mode_full = c(posterior_traj_mode_1, posterior_traj_mode_2)

  N_v = as.numeric(lapply(sim_data_list, nrow))
  N = sum(N_v)
  D = length(N_v)

  data_starts = do.call(rbind, lapply(1:D, FUN = function(i, l){l[[i]][1:(N_v[i]-1),]}, l = sim_data_list))
  data_ends = do.call(rbind, lapply(1:D, FUN = function(i, l){l[[i]][-1,]}, l = sim_data_list))

  N_advects = nrow(data_starts)

  rand_vel_1 = matrix(rep(0, N_advects*N_prop_steps), nrow = N_prop_steps, ncol = N_advects)
  rand_vel_2 = matrix(rep(0, N_advects*N_prop_steps), nrow = N_prop_steps, ncol = N_advects)

  omega = (1:M-1)*pi
  spec_den = sqrt(2*pi)*prior_l*exp(-0.5*prior_l^2*omega^2)

  prior_beta_sigma = diag(c(prior_k^2 * diag(spec_den) %*% matrix(rep(1,4), nrow = 2) %*% diag(spec_den)))

  start_beta = c(rep(0,2*M^2))

  t_steps_vec <- as.numeric(data_ends$t - data_starts$t) / N_prop_steps

  starts_mat <- as.matrix(data_starts)
  ends_mat <- as.matrix(data_ends)
  rand_vel_1_mat <- as.matrix(rand_vel_1)
  rand_vel_2_mat <- as.matrix(rand_vel_2)
  border_vec <- as.numeric(border)
  start_beta_vec <- as.numeric(posterior_traj_mode_full)

  # Good catch on the prior! Pre-calculate the precision matrix (inverse of covariance) here
  prior_precision_mat <- ginv(as.matrix(prior_beta_sigma))

  eval_counter <- 0

  cat(sprintf("\n\n--- Starting Optimization for New Sample ---\n"))

  # Optimize for the best w that matched the mixing fields
  rMAP_loss <- function(beta) {

    eval_counter <<- eval_counter + 1

    # The only thing happening here is the jump to C++
    current_loss <- rMAP_loss_cpp(
      beta = beta,
      M_sq = M^2,
      aug_data_starts = starts_mat,
      aug_data_ends = ends_mat,
      t_steps = t_steps_vec,
      rand_vel_1 = rand_vel_1_mat,
      rand_vel_2 = rand_vel_2_mat,
      border = border_vec,
      pos_sd = pos_sd,
      prior_beta_sigma = prior_precision_mat,
      start_beta = start_beta_vec
    )

    if (eval_counter %% print_every == 0) {
      beta_str <- paste(sprintf("%.3f", beta), collapse = ", ")

      # Notice the \n at the VERY BEGINNING here too, to clear the progress bar
      cat(sprintf("\nIter: %4d | Loss: %10.4f | Beta: [%s]",
                  eval_counter, current_loss, beta_str))
    }

    return(current_loss)

  }

  H_matrix = numDeriv::hessian(func = rMAP_loss, x = start_beta_vec)

  H_matrix

}

posterior_traj_mode_mat_1 = test_posterior_traj_modes$Beta1
posterior_traj_mode_mat_2 = test_posterior_traj_modes$Beta2

get_posterior_traj_hessians = function(posterior_traj_mode_mat_1, posterior_traj_mode_mat_2, sim_data_list, pos_sd = 0.001, vel_sd = 0.1, M = 2, prior_k = 0.35, prior_l = 1, baseVectorFields_Vec, border, N_prop_steps = 10, traj_eval_grid, print_every = 50){

  #Define matrices for mode coefficients, mode path evaluations, NLL evaluations

  N_traj_modes = nrow(posterior_traj_mode_mat_1)
  N_beta = ncol(posterior_traj_mode_mat_1)
  N_evals = nrow(traj_eval_grid)

  traj_hessians = list()

  #Loops through all mode starting locations

  for(i in 1:N_traj_modes){

    cur_traj_mode_1 = posterior_traj_mode_mat_1[i,]
    cur_traj_mode_2 = posterior_traj_mode_mat_2[i,]

    cur_trajMode_H = get_one_traj_hessian(posterior_traj_mode_1 = cur_traj_mode_1, posterior_traj_mode_2 = cur_traj_mode_2, sim_data_list = sim_data_list, pos_sd = pos_sd, vel_sd = vel_sd, M = M, prior_k = prior_k, prior_l = prior_l, baseVectorFields_Vec = baseVectorFields_Vec, border = border, N_prop_steps = N_prop_steps, traj_eval_grid = traj_eval_grid, print_every = print_every)

    traj_hessians[[i]] = cur_trajMode_H

  }

  return(traj_hessians)

}

N_draws = 100
traj_mode = c(posterior_traj_mode_mat_1[1,], posterior_traj_mode_mat_2[1,])
traj_hessian = traj_hessians[[1]]
temp = 1

run_RWMH_one_traj_mode = function(N_draws, traj_mode, traj_hessian, temp,  sim_data_list, pos_sd = 0.001, vel_sd = 0.1, M = 2, prior_k = 0.35, prior_l = 1, baseVectorFields_Vec, border, N_prop_steps = 10){

  N_v = as.numeric(lapply(sim_data_list, nrow))
  N = sum(N_v)
  D = length(N_v)

  data_starts = do.call(rbind, lapply(1:D, FUN = function(i, l){l[[i]][1:(N_v[i]-1),]}, l = sim_data_list))
  data_ends = do.call(rbind, lapply(1:D, FUN = function(i, l){l[[i]][-1,]}, l = sim_data_list))

  N_advects = nrow(data_starts)

  omega = (1:M-1)*pi
  spec_den = sqrt(2*pi)*prior_l*exp(-0.5*prior_l^2*omega^2)

  prior_beta_sigma = diag(c(prior_k^2 * diag(spec_den) %*% matrix(rep(1,M^2), nrow = M) %*% diag(spec_den)))
  prior_precision_mat <- ginv(as.matrix(prior_beta_sigma))

  t_steps_vec <- as.numeric(data_ends$t - data_starts$t) / N_prop_steps

  starts_mat <- as.matrix(data_starts)
  ends_mat <- as.matrix(data_ends)
  border_vec <- as.numeric(border)

  #Find gaussian covariance estimate

  traj_estSigma = ginv(traj_hessian)

  #Calculate covariance of the proposal distributions
  #    The 2.382^2/d should yield optimal acceptance rate
  #    Adding the temperature for added flexibility

  d = length(traj_mode)

  traj_Proposal_Cov = traj_estSigma * ((2.382)^2/d) * temp

  traj_coefs = matrix(nrow = N_draws, ncol = d)
  LL = rep(0, N_draws)
  accept = rep(0, N_draws)

  prev_model_LL = -Inf
  prev_traj = traj_mode

  prior_beta = rep(0, d)

  svMisc::progress(0, N_draws)

  for(i in 1:N_draws){

    cur_traj = mvrnorm(n=1, mu = prev_traj, Sigma = traj_Proposal_Cov)

    cur_vel_err_1 = matrix(rnorm(N_advects*N_prop_steps, 0, vel_sd), nrow = N_prop_steps, ncol = N_advects)
    cur_vel_err_2 = matrix(rnorm(N_advects*N_prop_steps, 0, vel_sd), nrow = N_prop_steps, ncol = N_advects)

    cur_model_LL = -1*rMAP_loss_cpp(
      beta = cur_traj,
      M_sq = M^2,
      aug_data_starts = starts_mat,
      aug_data_ends = ends_mat,
      t_steps = t_steps_vec,
      rand_vel_1 = cur_vel_err_1,
      rand_vel_2 = cur_vel_err_2,
      border = border_vec,
      pos_sd = pos_sd,
      prior_beta_sigma = prior_precision_mat,
      start_beta = prior_beta
    )

    A = (cur_model_LL) - (prev_model_LL)

    log_u = log(runif(1))

    if(log_u <= A){

      prev_traj = cur_traj
      prev_model_LL = cur_model_LL

      traj_coefs[i,] = cur_traj
      LL[i] = cur_model_LL
      accept[i] = 1

    } else{

      traj_coefs[i,] = prev_traj
      LL[i] = prev_model_LL

    }

    svMisc::progress(i, N_draws)

  }

  return(list(Traj_Draws = traj_coefs, LogLikelihoods = LL, Acceptance = accept))

}

RWMH_test = run_RWMH_one_traj_mode(N_draws = 100, traj_mode = traj_mode, traj_hessian = traj_hessian, temp = 1, sim_data_list = sim_data_list, pos_sd = 0.001, vel_sd = 0.1, M = 2, prior_k = 0.35, prior_l = 1, baseVectorFields_Vec = baseVectorFields_Vec, border = border, N_prop_steps = 10)

plot(RWMH_test$LogLikelihoods)
mean(RWMH_test$Acceptance)


RWMH_test$Traj_Draws

### Tests

M=2
true_l = 1
true_k = 0.35

trueHSGP_1 = sampleFullTrajectoriesHSGP(2, M = M-1, log_k = true_k, log_l = true_l)
trueHSGP_2 = sampleFullTrajectoriesHSGP(2, M = M-1, log_k = true_k, log_l = true_l)

trueHSGP$log_z[[2]] = matrix(rep(-100, 4), nrow = 2)
trueHSGP$log_z[[1]] = matrix(rep(0, 4), nrow = 2)

omega = (1:M-1)*pi
spec_den = sqrt(2*pi)*true_l*exp(-0.5*true_l^2*omega^2)

prior_beta_sigma = diag(c(true_k^2 * diag(spec_den) %*% matrix(rep(1,M^2), nrow = M) %*% diag(spec_den)))

true_beta_1 = c(true_k*diag(sqrt(spec_den)) %*% trueHSGP$log_z[[1]] %*% diag(sqrt(spec_den)))
true_beta_2 = c(true_k*diag(sqrt(spec_den)) %*% trueHSGP$log_z[[2]] %*% diag(sqrt(spec_den)))

true_beta_1_1 = c(true_k*diag(sqrt(spec_den)) %*% trueHSGP_1$log_z[[1]] %*% diag(sqrt(spec_den)))
true_beta_1_2 = c(true_k*diag(sqrt(spec_den)) %*% trueHSGP_1$log_z[[2]] %*% diag(sqrt(spec_den)))

true_beta_mat = cbind(true_beta_1, true_beta_2)

true_beta_mat_1 = cbind(true_beta_1_1, true_beta_1_2)

sampledParticles = samplePhySpaceParticles(n_particles = 100, startTime = 0, n_obs = 20*100, border = c(-2,-2,2,2), borderBuffer = 0.2, baseVectorFields, trueHSGP,
                                           M = 2, t_step_mean = 0.1, vel_sigma = 0.001, pos_sigma = 0.001)

sampledParticles_1 = samplePhySpaceParticles(n_particles = 100, startTime = 0, n_obs = 20*10, border = c(-2,-2,2,2), borderBuffer = 0.2, baseVectorFields, trueHSGP_1,
                                             M = 2, t_step_mean = 0.001, vel_sigma = 0.001, pos_sigma = 0.001)

sampledParticles = rbind(sampledParticles_1, sampledParticles_2)

sampledParticles$Particle = rep(str_c('Particle', 1:100), each = 200)

sampledParticles_Sub = sampledParticles_1[1:(100*20) * 10 - (10-1),]
ggplot(sampledParticles_Sub, aes(x = X1, y = X2, color = Particle)) + geom_point() + theme(legend.position = 'None')

sim_data_list = list()

for(d in 1:100){

  curDrifter = sampledParticles_Sub[sampledParticles_Sub$Particle == str_c('Particle',d),]

  sim_data_list[[d]] = curDrifter[,c(1,2,3)]


}

sim_data_list = list(data.frame(t = 0:19*(9*pi), X1 = rep(0,20), X2 = rep(c(1,-1),10)))


traj_eval_grid = expand.grid(seq(-2,2, length.out = 100), seq(-2,2, length.out = 100))

TrueTrajVals = exp(evaluate2DCosine_fast(true_beta_mat_1, pos_mat = traj_eval_grid, border = c(-2,-2,2,2)))

t1 = Sys.time()

real_traj_test_1 = run_rMAP_Trajectory_Multithread(N_samples = 100, sim_data_list = sim_data_list,
                                       pos_sd = 0.001, vel_sd = 0, M = 1, prior_k = 0.35,
                                       prior_l = 1, baseVectorFields_Vec = baseVectorFields_Vec,
                                       border = c(-2,-2,2,2), N_prop_steps = 1000, traj_eval_grid = traj_eval_grid, print_every = 50, n_threads = 32)

t2 = Sys.time()

t2-t1

max_log_w = max(real_traj_test_1$LogImportanceWeights)

log_w_scaled = real_traj_test_1$LogImportanceWeights - max_log_w

importance_weights = exp(log_w_scaled) / sum(exp(log_w_scaled))

which.max(importance_weights)
which.min(real_traj_test_1$PosteriorNLL)

exp(c(real_traj_test_1$Beta1Posterior[20,], real_traj_test_1$Beta2Posterior[20,]))

real_traj_test_1$Traj1Posterior

sample(exp(real_traj_test_1$Beta1Posterior)*9, size = 1000, replace = T, prob = importance_weights)

sort(importance_weights)

PostTraj1_Mean = rowMeans(real_traj_test_1$Traj1Posterior)
PostTraj2_Mean = rowMeans(real_traj_test_1$Traj2Posterior)

hist((PostTraj1_Mean - TrueTrajVals[,1])^2)
hist((PostTraj2_Mean - TrueTrajVals[,2])^2)

CI_traj_1 = t(apply(real_traj_test_1$Traj1Posterior, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975)))
CI_traj_2 = t(apply(real_traj_test_1$Traj2Posterior, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975)))

mean(CI_traj_1[,1] <= TrueTrajVals[,1] & CI_traj_1[,2] >= TrueTrajVals[,1])
mean(CI_traj_2[,1] <= TrueTrajVals[,2] & CI_traj_2[,2] >= TrueTrajVals[,2])

hist(CI_traj_1[,2] - CI_traj_1[,1])
hist(CI_traj_2[,2] - CI_traj_2[,1])


ggplot() + geom_histogram(aes(x = (exp(real_traj_test_1$Beta1Posterior)*9)[exp(real_traj_test_1$Beta1Posterior)*9 < 20]), binwidth = 1/10)


plot(exp(real_traj_test_1$Beta1Posterior)*9, real_traj_test_1$LogImportanceWeights)

plot(real_traj_test_1$LogImportanceWeights)

which(real_traj_test_1$Beta1Posterior > 2)

colMeans(real_traj_test_1$Beta1Posterior)
colMeans(real_traj_test_1$Beta2Posterior)

plot(exp(real_traj_test_1$Beta1Posterior)*9, exp(real_traj_test_1$Beta2Posterior))

post_loss = rMAP_loss_cpp_multi(
  beta = c(real_traj_test_1$Beta1Posterior[62,], real_traj_test_1$Beta2Posterior[62,]),
  M_sq = M^2,
  aug_data_starts = aug_starts_mat,
  aug_data_ends = data_ends_mat,
  t_steps = t_steps_vec,
  rand_vel_1 = rand_vel_1_mat,
  rand_vel_2 = rand_vel_2_mat,
  border = border_vec,
  pos_sd = pos_sd,
  prior_beta_sigma = prior_precision_mat,
  start_beta = start_beta_vec,
  n_threads = n_threads
)

M=1

log_importance = calculate_log_importance_cpp(
  beta = c(real_traj_test_1$Beta1Posterior[62,], real_traj_test_1$Beta2Posterior[62,]),
  M_sq = M^2,
  start_t_pos_mat = aug_starts_mat,
  end_t_pos_true_mat = data_ends_mat,
  t_steps = data_ends_mat[,1] - data_starts_mat[,1],
  N_prop_steps = N_prop_steps,
  border = border,
  pos_sd = pos_sd,
  prior_beta_sigma = prior_precision_mat,
  n_threads = n_threads
)
