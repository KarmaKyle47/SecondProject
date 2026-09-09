library(plotly)
library(stringr)
library(MASS)
library(ggplot2)
library(minqa)
library(nloptr)
library(splines2)
library(parallel)
library(pbapply)


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

evaluate2DCosine_fast = function(beta_mat, pos_mat, border){

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

TrajWeightedBaseVectorFields_2D_Cosine = function(pos_t_mat, beta_mat, baseVectorFields_Vec, border){


  log_Traj = evaluate2DCosine_fast(beta_mat = beta_mat, pos_mat = pos_t_mat[,c(2,3)], border = border)
  log_Traj_mat = cbind(log_Traj, log_Traj)

  baseVF = baseVectorFields_Vec(pos_t_mat)

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

  c = sqrt(rowSums(pos_t_mat[,c(2,3)]^2))

  f1x = pos_t_mat[,3] / c
  f1y = -1*pos_t_mat[,2] / c

  f2x = pos_t_mat[,2] / c
  f2y = pos_t_mat[,3] / c

  cbind(f1x,f2x, f1y, f2y)

}


######## Start Building the Path Mode Finder Function #####################

trajectoryPost = cbind(real_traj_test_1$Beta1Posterior, real_traj_test_1$Beta2Posterior)
data = sim_data_list[[1]]
pos_sd = 0.001
pos_selection_sd = 0.001
vel_sd = 0.001

full_nodes = seq(0, max(data$t), length.out = 252)[-c(1,252)]
N_quad = 1000
N_prop_steps = 100

find_one_rMAP_Path = function(data, pos_sd, vel_sd, pos_selection_sd, trajectoryPost, border, N_quad, baseVectorFields_Vec, full_nodes, N_prop_steps, print_every, plot){

  # Step 1: Sample Trajectory - trajectory Post is a data.frame or matrix with each posterior draw as a row

  sampledTrajectory = matrix(trajectoryPost[sample(1:nrow(trajectoryPost), size = 1),], ncol = 2, byrow = F)

  # Step 2: Augment Data with positional error
  # data is Nx3 with columns t, X1, X2

  N = nrow(data)
  aug_data = data + matrix(c(rep(0, N), rnorm(N, 0, pos_sd), rnorm(N, 0, pos_sd)), byrow = F, ncol = 3)

  aug_data_starts = aug_data[1:(N-1),]
  aug_data_ends = aug_data[2:N,]

  t_steps <- (aug_data_ends$t - aug_data_starts$t) / N_prop_steps
  f_curPos_mat <- aug_data_starts
  b_curPos_mat = aug_data_ends

  f_all_steps = data.frame(t = rep(0,N_prop_steps*(N-1) + 1), X1 = rep(0,N_prop_steps*(N-1) + 1), X2 = rep(0,N_prop_steps*(N-1) + 1))
  b_all_steps = data.frame(t = rep(0,N_prop_steps*(N-1) + 1), X1 = rep(0,N_prop_steps*(N-1) + 1), X2 = rep(0,N_prop_steps*(N-1) + 1))

  f_all_steps[1,] = f_curPos_mat[1,]
  b_all_steps[N_prop_steps*(N-1) + 1,] = b_curPos_mat[N-1,]

  for(j in 1:N_prop_steps) {

    #Forward

    f_k1 <- TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = f_curPos_mat, beta_mat = sampledTrajectory, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    f_pos_mat_k2 <- f_curPos_mat
    f_pos_mat_k2[, 1] <- f_pos_mat_k2[, 1] + t_steps / 2
    f_pos_mat_k2[, 2] <- f_pos_mat_k2[, 2] + f_k1[, 1] * (t_steps / 2)
    f_pos_mat_k2[, 3] <- f_pos_mat_k2[, 3] + f_k1[, 2] * (t_steps / 2)

    f_k2 <- TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = f_pos_mat_k2, beta_mat = sampledTrajectory, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    f_pos_mat_k3 <- f_curPos_mat
    f_pos_mat_k3[, 1] <- f_pos_mat_k3[, 1] + t_steps / 2
    f_pos_mat_k3[, 2] <- f_pos_mat_k3[, 2] + f_k2[, 1] * (t_steps / 2)
    f_pos_mat_k3[, 3] <- f_pos_mat_k3[, 3] + f_k2[, 2] * (t_steps / 2)

    f_k3 <- TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = f_pos_mat_k3, beta_mat = sampledTrajectory, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    f_pos_mat_k4 <- f_curPos_mat
    f_pos_mat_k4[, 1] <- f_pos_mat_k4[, 1] + t_steps
    f_pos_mat_k4[, 2] <- f_pos_mat_k4[, 2] + f_k3[, 1] * t_steps
    f_pos_mat_k4[, 3] <- f_pos_mat_k4[, 3] + f_k3[, 2] * t_steps

    f_k4 <- TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = f_pos_mat_k4, beta_mat = sampledTrajectory, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    f_rk4_drift_1 <- (f_k1[, 1] + 2 * f_k2[, 1] + 2 * f_k3[, 1] + f_k4[, 1]) / 6
    f_rk4_drift_2 <- (f_k1[, 2] + 2 * f_k2[, 2] + 2 * f_k3[, 2] + f_k4[, 2]) / 6

    f_curPos_mat[, 1] <- f_curPos_mat[, 1] + t_steps
    f_curPos_mat[, 2] <- f_curPos_mat[, 2] + f_rk4_drift_1 * t_steps
    f_curPos_mat[, 3] <- f_curPos_mat[, 3] + f_rk4_drift_2 * t_steps

    f_all_steps[0:(N-2) * N_prop_steps + j + 1,1] = f_curPos_mat[,1]
    f_all_steps[0:(N-2) * N_prop_steps + j + 1,2] = f_curPos_mat[,2]
    f_all_steps[0:(N-2) * N_prop_steps + j + 1,3] = f_curPos_mat[,3]

    #Backward

    b_k1 <- -1*TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = b_curPos_mat, beta_mat = sampledTrajectory, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    b_pos_mat_k2 <- b_curPos_mat
    b_pos_mat_k2[, 1] <- b_pos_mat_k2[, 1] - t_steps / 2
    b_pos_mat_k2[, 2] <- b_pos_mat_k2[, 2] + b_k1[, 1] * (t_steps / 2)
    b_pos_mat_k2[, 3] <- b_pos_mat_k2[, 3] + b_k1[, 2] * (t_steps / 2)

    b_k2 <- -1*TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = b_pos_mat_k2, beta_mat = sampledTrajectory, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    b_pos_mat_k3 <- b_curPos_mat
    b_pos_mat_k3[, 1] <- b_pos_mat_k3[, 1] - t_steps / 2
    b_pos_mat_k3[, 2] <- b_pos_mat_k3[, 2] + b_k2[, 1] * (t_steps / 2)
    b_pos_mat_k3[, 3] <- b_pos_mat_k3[, 3] + b_k2[, 2] * (t_steps / 2)

    b_k3 <- -1*TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = b_pos_mat_k3, beta_mat = sampledTrajectory, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    b_pos_mat_k4 <- b_curPos_mat
    b_pos_mat_k4[, 1] <- b_pos_mat_k4[, 1] - t_steps
    b_pos_mat_k4[, 2] <- b_pos_mat_k4[, 2] + b_k3[, 1] * t_steps
    b_pos_mat_k4[, 3] <- b_pos_mat_k4[, 3] + b_k3[, 2] * t_steps

    b_k4 <- -1*TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = b_pos_mat_k4, beta_mat = sampledTrajectory, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    b_rk4_drift_1 <- (b_k1[, 1] + 2 * b_k2[, 1] + 2 * b_k3[, 1] + b_k4[, 1]) / 6
    b_rk4_drift_2 <- (b_k1[, 2] + 2 * b_k2[, 2] + 2 * b_k3[, 2] + b_k4[, 2]) / 6

    b_curPos_mat[, 1] <- b_curPos_mat[, 1] - t_steps
    b_curPos_mat[, 2] <- b_curPos_mat[, 2] + b_rk4_drift_1 * t_steps
    b_curPos_mat[, 3] <- b_curPos_mat[, 3] + b_rk4_drift_2 * t_steps

    b_all_steps[1:(N-1) * N_prop_steps + 1 - j,1] = b_curPos_mat[,1]
    b_all_steps[1:(N-1) * N_prop_steps + 1 - j,2] = b_curPos_mat[,2]
    b_all_steps[1:(N-1) * N_prop_steps + 1 - j,3] = b_curPos_mat[,3]

  }

  Prop_Avg_Pos_mat = (f_all_steps + b_all_steps) / 2

  # Step 3: Get prior path using a bSpline regression with prior_nodes

  Phi_prior = bSpline(Prop_Avg_Pos_mat$t, knots = full_nodes, intercept = T)
  Phi_prior_d = bSpline(Prop_Avg_Pos_mat$t, knots = full_nodes, intercept = T, derivs = 1)

  Phi_prior_H = solve(t(Phi_prior) %*% Phi_prior) %*% t(Phi_prior)
  prior_c_x = Phi_prior_H %*% Prop_Avg_Pos_mat$X1
  prior_c_y = Phi_prior_H %*% Prop_Avg_Pos_mat$X2

  # Step 4: Get full path basis to add onto prior with full nodes

  t_quad = seq(min(aug_data$t), max(aug_data$t), length.out = N_quad)

  Phi_data_full = bSpline(aug_data$t, knots = full_nodes, intercept = T)
  Phi_full = bSpline(t_quad, knots = full_nodes, intercept = T)
  Phi_full_d = bSpline(t_quad, knots = full_nodes, intercept = T, derivs = 1)
  N_param = ncol(Phi_full)

  # Step 5: Sample random initial path
  Lt = max(aug_data$t) - min(aug_data$t)

  helper_M = t(Phi_full) %*% Phi_full
  helper_M_inv = ginv(helper_M)

  c_prior_sigma = (N_quad / Lt * pos_selection_sd^2) * helper_M_inv
  inv_c_prior_sigma = (Lt / (N_quad * pos_selection_sd^2)) * helper_M

  start_c_x = prior_c_x + mvrnorm(mu = rep(0, N_param), Sigma = c_prior_sigma)
  start_c_y = prior_c_y + mvrnorm(mu = rep(0, N_param), Sigma = c_prior_sigma)

  start_c = c(start_c_x, start_c_y)

  rand_vel_x = rnorm(N_quad, 0, sd = vel_sd)
  rand_vel_y = rnorm(N_quad, 0, sd = vel_sd)

  # =================================================================
  # Step 6: Define loss function
  # =================================================================

  aug_X_matrix = as.matrix(aug_data[, c("X1", "X2")])
  eval_counter <- 0

  # Clean visual header for a new optimization run
  cat("\n=======================================================\n")
  cat("       Starting Optimization for New Sample            \n")
  cat("=======================================================\n\n")

  rMAP_loss = function(c) {

    eval_counter <<- eval_counter + 1

    # Call the compiled C++ function
    loss = spline_rMAP_loss_cpp(
      c_params          = c,
      Phi_data          = Phi_data_full,
      Phi_quad          = Phi_full,
      Phi_deriv_quad    = Phi_full_d,
      aug_X             = aug_X_matrix,
      sampledTrajectory = sampledTrajectory,
      border            = border,
      pos_sd            = pos_sd,
      vel_sd            = vel_sd,
      rand_vel_x        = rand_vel_x,
      rand_vel_y        = rand_vel_y,
      prior_c           = start_c,
      inv_c_prior_sigma = inv_c_prior_sigma,
      Lt                = Lt
    )

    if (eval_counter %% print_every == 0) {

      # Truncate C array so it doesn't word-wrap in the console
      n_c = length(c)
      if (n_c > 6) {
        c_str = paste0(sprintf("%.3f, %.3f, %.3f", c[1], c[2], c[3]),
                       ", ... , ",
                       sprintf("%.3f, %.3f, %.3f", c[n_c-2], c[n_c-1], c[n_c]))
      } else {
        c_str = paste(sprintf("%.3f", c), collapse = ", ")
      }

      # Formatted output with trailing newline
      cat(sprintf("  Iter: %4d  |  Loss: %12.4f  |  C: [%s]\n",
                  eval_counter, loss, c_str))
    }

    return(loss)
  }

  # =================================================================
  # Step 7: Run the optimization
  # =================================================================

  opt_result = nloptr(
    x0 = start_c,
    eval_f = rMAP_loss,
    opts = list(
      "algorithm"   = "NLOPT_LN_NEWUOA",
      "ftol_rel"    = 1e-6,
      "maxeval"     = 5000,
      "print_level" = 0
    )
  )

  # Format the Final output identically
  final_c = opt_result$solution
  n_fc = length(final_c)
  if (n_fc > 6) {
    final_c_str = paste0(sprintf("%.3f, %.3f, %.3f", final_c[1], final_c[2], final_c[3]),
                         ", ... , ",
                         sprintf("%.3f, %.3f, %.3f", final_c[n_fc-2], final_c[n_fc-1], final_c[n_fc]))
  } else {
    final_c_str = paste(sprintf("%.3f", final_c), collapse = ", ")
  }

  # Add visual footer to close out the optimization block
  cat("\n-------------------------------------------------------\n")
  cat(sprintf("  FINAL Iter: %4d  |  Loss: %12.4f  |  C: [%s]\n",
              opt_result$iterations, opt_result$objective, final_c_str))
  cat("-------------------------------------------------------\n\n")

  # Step 8: Extract the optimized control points
  optimized_c = opt_result$solution

  opt_c_x = optimized_c[1:N_param]
  opt_c_y = optimized_c[1:N_param + N_param]

  # Eval Optimized Path

  opt_path_x = Phi_full %*% opt_c_x
  opt_path_y = Phi_full %*% opt_c_y

  opt_path_data_x = Phi_data_full %*% opt_c_x
  opt_path_data_y = Phi_data_full %*% opt_c_y

  opt_path_x_d = Phi_full_d %*% opt_c_x
  opt_path_y_d = Phi_full_d %*% opt_c_y

  opt_path_traj_vel = TrajWeightedBaseVectorFields_2D_Cosine(cbind(t_quad, opt_path_x, opt_path_y), sampledTrajectory, baseVectorFields_Vec, border)

  opt_path_traj_vel_x = opt_path_traj_vel[,1]
  opt_path_traj_vel_y = opt_path_traj_vel[,2]


  ## Positional Loss

  data_diff_x = data$X1 - opt_path_data_x
  data_diff_y = data$X2 - opt_path_data_y

  opt_pos_NLL = sum(data_diff_x^2 + data_diff_y^2) / (2 * pos_sd^2)

  ## Velocity Loss

  vel_diff_x = opt_path_x_d - opt_path_traj_vel_x
  vel_diff_y = opt_path_y_d - opt_path_traj_vel_y

  opt_vel_NLL = sum(vel_diff_x^2 + vel_diff_y^2) / (2 * vel_sd^2) * (Lt/N_quad)

  ## Prior Loss

  opt_prior_NLL = as.numeric((t(opt_c_x - start_c_x) %*% inv_c_prior_sigma %*% (opt_c_x - start_c_x) + t(opt_c_y - start_c_y) %*% inv_c_prior_sigma %*% (opt_c_y - start_c_y)) / 2)

  opt_NLL = opt_pos_NLL + opt_vel_NLL + opt_prior_NLL

  cat("\n")

  if(plot){

    ggplot() + geom_path(aes(x = opt_path_x, y = opt_path_y), size = 0.75) + geom_point(data = aug_data, aes(x = X1, y = X2), color = 'red')

  }

  list(Optimized_C = cbind(opt_c_x, opt_c_y), Optimized_Path = cbind(opt_path_x, opt_path_y), Posterior_NLL_Position = opt_pos_NLL, Posterior_NLL_Velocity = opt_vel_NLL, Posterior_NLL_Prior = opt_prior_NLL)

}

find_one_rMAP_Path_Multithread = function(data, pos_sd, vel_sd, pos_selection_sd, trajectoryPost, border, N_quad, baseVectorFields_Vec, full_nodes, N_prop_steps, print_every, plot, n_threads){

  # Step 1: Sample Trajectory - trajectory Post is a data.frame or matrix with each posterior draw as a row

  sampledTrajectory = matrix(trajectoryPost[sample(1:nrow(trajectoryPost), size = 1),], ncol = 2, byrow = F)

  # Step 2: Augment Data with positional error
  # data is Nx3 with columns t, X1, X2

  N = nrow(data)
  aug_data = data + matrix(c(rep(0, N), rnorm(N, 0, pos_sd), rnorm(N, 0, pos_sd)), byrow = F, ncol = 3)

  aug_data_starts = aug_data[1:(N-1),]
  aug_data_ends = aug_data[2:N,]

  t_steps <- (aug_data_ends$t - aug_data_starts$t) / N_prop_steps
  f_curPos_mat <- aug_data_starts
  b_curPos_mat = aug_data_ends

  f_all_steps = data.frame(t = rep(0,N_prop_steps*(N-1) + 1), X1 = rep(0,N_prop_steps*(N-1) + 1), X2 = rep(0,N_prop_steps*(N-1) + 1))
  b_all_steps = data.frame(t = rep(0,N_prop_steps*(N-1) + 1), X1 = rep(0,N_prop_steps*(N-1) + 1), X2 = rep(0,N_prop_steps*(N-1) + 1))

  f_all_steps[1,] = f_curPos_mat[1,]
  b_all_steps[N_prop_steps*(N-1) + 1,] = b_curPos_mat[N-1,]

  for(j in 1:N_prop_steps) {

    #Forward

    f_k1 <- TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = f_curPos_mat, beta_mat = sampledTrajectory, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    f_pos_mat_k2 <- f_curPos_mat
    f_pos_mat_k2[, 1] <- f_pos_mat_k2[, 1] + t_steps / 2
    f_pos_mat_k2[, 2] <- f_pos_mat_k2[, 2] + f_k1[, 1] * (t_steps / 2)
    f_pos_mat_k2[, 3] <- f_pos_mat_k2[, 3] + f_k1[, 2] * (t_steps / 2)

    f_k2 <- TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = f_pos_mat_k2, beta_mat = sampledTrajectory, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    f_pos_mat_k3 <- f_curPos_mat
    f_pos_mat_k3[, 1] <- f_pos_mat_k3[, 1] + t_steps / 2
    f_pos_mat_k3[, 2] <- f_pos_mat_k3[, 2] + f_k2[, 1] * (t_steps / 2)
    f_pos_mat_k3[, 3] <- f_pos_mat_k3[, 3] + f_k2[, 2] * (t_steps / 2)

    f_k3 <- TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = f_pos_mat_k3, beta_mat = sampledTrajectory, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    f_pos_mat_k4 <- f_curPos_mat
    f_pos_mat_k4[, 1] <- f_pos_mat_k4[, 1] + t_steps
    f_pos_mat_k4[, 2] <- f_pos_mat_k4[, 2] + f_k3[, 1] * t_steps
    f_pos_mat_k4[, 3] <- f_pos_mat_k4[, 3] + f_k3[, 2] * t_steps

    f_k4 <- TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = f_pos_mat_k4, beta_mat = sampledTrajectory, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    f_rk4_drift_1 <- (f_k1[, 1] + 2 * f_k2[, 1] + 2 * f_k3[, 1] + f_k4[, 1]) / 6
    f_rk4_drift_2 <- (f_k1[, 2] + 2 * f_k2[, 2] + 2 * f_k3[, 2] + f_k4[, 2]) / 6

    f_curPos_mat[, 1] <- f_curPos_mat[, 1] + t_steps
    f_curPos_mat[, 2] <- f_curPos_mat[, 2] + f_rk4_drift_1 * t_steps
    f_curPos_mat[, 3] <- f_curPos_mat[, 3] + f_rk4_drift_2 * t_steps

    f_all_steps[0:(N-2) * N_prop_steps + j + 1,1] = f_curPos_mat[,1]
    f_all_steps[0:(N-2) * N_prop_steps + j + 1,2] = f_curPos_mat[,2]
    f_all_steps[0:(N-2) * N_prop_steps + j + 1,3] = f_curPos_mat[,3]

    #Backward

    b_k1 <- -1*TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = b_curPos_mat, beta_mat = sampledTrajectory, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    b_pos_mat_k2 <- b_curPos_mat
    b_pos_mat_k2[, 1] <- b_pos_mat_k2[, 1] - t_steps / 2
    b_pos_mat_k2[, 2] <- b_pos_mat_k2[, 2] + b_k1[, 1] * (t_steps / 2)
    b_pos_mat_k2[, 3] <- b_pos_mat_k2[, 3] + b_k1[, 2] * (t_steps / 2)

    b_k2 <- -1*TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = b_pos_mat_k2, beta_mat = sampledTrajectory, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    b_pos_mat_k3 <- b_curPos_mat
    b_pos_mat_k3[, 1] <- b_pos_mat_k3[, 1] - t_steps / 2
    b_pos_mat_k3[, 2] <- b_pos_mat_k3[, 2] + b_k2[, 1] * (t_steps / 2)
    b_pos_mat_k3[, 3] <- b_pos_mat_k3[, 3] + b_k2[, 2] * (t_steps / 2)

    b_k3 <- -1*TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = b_pos_mat_k3, beta_mat = sampledTrajectory, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    b_pos_mat_k4 <- b_curPos_mat
    b_pos_mat_k4[, 1] <- b_pos_mat_k4[, 1] - t_steps
    b_pos_mat_k4[, 2] <- b_pos_mat_k4[, 2] + b_k3[, 1] * t_steps
    b_pos_mat_k4[, 3] <- b_pos_mat_k4[, 3] + b_k3[, 2] * t_steps

    b_k4 <- -1*TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = b_pos_mat_k4, beta_mat = sampledTrajectory, baseVectorFields_Vec = baseVectorFields_Vec, border = border
    )

    b_rk4_drift_1 <- (b_k1[, 1] + 2 * b_k2[, 1] + 2 * b_k3[, 1] + b_k4[, 1]) / 6
    b_rk4_drift_2 <- (b_k1[, 2] + 2 * b_k2[, 2] + 2 * b_k3[, 2] + b_k4[, 2]) / 6

    b_curPos_mat[, 1] <- b_curPos_mat[, 1] - t_steps
    b_curPos_mat[, 2] <- b_curPos_mat[, 2] + b_rk4_drift_1 * t_steps
    b_curPos_mat[, 3] <- b_curPos_mat[, 3] + b_rk4_drift_2 * t_steps

    b_all_steps[1:(N-1) * N_prop_steps + 1 - j,1] = b_curPos_mat[,1]
    b_all_steps[1:(N-1) * N_prop_steps + 1 - j,2] = b_curPos_mat[,2]
    b_all_steps[1:(N-1) * N_prop_steps + 1 - j,3] = b_curPos_mat[,3]

  }

  Prop_Avg_Pos_mat = (f_all_steps + b_all_steps) / 2

  # Step 3: Get prior path using a bSpline regression with prior_nodes

  Phi_prior = bSpline(Prop_Avg_Pos_mat$t, knots = full_nodes, intercept = T)
  Phi_prior_d = bSpline(Prop_Avg_Pos_mat$t, knots = full_nodes, intercept = T, derivs = 1)

  Phi_prior_H = solve(t(Phi_prior) %*% Phi_prior) %*% t(Phi_prior)
  prior_c_x = Phi_prior_H %*% Prop_Avg_Pos_mat$X1
  prior_c_y = Phi_prior_H %*% Prop_Avg_Pos_mat$X2

  # Step 4: Get full path basis to add onto prior with full nodes

  t_quad = seq(min(aug_data$t), max(aug_data$t), length.out = N_quad)

  Phi_data_full = bSpline(aug_data$t, knots = full_nodes, intercept = T)
  Phi_full = bSpline(t_quad, knots = full_nodes, intercept = T)
  Phi_full_d = bSpline(t_quad, knots = full_nodes, intercept = T, derivs = 1)
  N_param = ncol(Phi_full)

  # Step 5: Sample random initial path
  Lt = max(aug_data$t) - min(aug_data$t)

  helper_M = t(Phi_full) %*% Phi_full
  helper_M_inv = ginv(helper_M)

  c_prior_sigma = (N_quad / Lt * pos_selection_sd^2) * helper_M_inv
  inv_c_prior_sigma = (Lt / (N_quad * pos_selection_sd^2)) * helper_M

  start_c_x = prior_c_x + mvrnorm(mu = rep(0, N_param), Sigma = c_prior_sigma)
  start_c_y = prior_c_y + mvrnorm(mu = rep(0, N_param), Sigma = c_prior_sigma)

  start_c = c(start_c_x, start_c_y)

  rand_vel_x = rnorm(N_quad, 0, sd = vel_sd)
  rand_vel_y = rnorm(N_quad, 0, sd = vel_sd)

  # =================================================================
  # Step 6: Define loss function
  # =================================================================

  aug_X_matrix = as.matrix(aug_data[, c("X1", "X2")])
  eval_counter <- 0

  # Clean visual header for a new optimization run
  cat("\n=======================================================\n")
  cat("       Starting Optimization for New Sample            \n")
  cat("=======================================================\n\n")

  rMAP_loss = function(c) {

    eval_counter <<- eval_counter + 1

    # Call the compiled C++ function
    loss = spline_rMAP_loss_cpp_multi(
      c_params          = c,
      Phi_data          = Phi_data_full,
      Phi_quad          = Phi_full,
      Phi_deriv_quad    = Phi_full_d,
      aug_X             = aug_X_matrix,
      sampledTrajectory = sampledTrajectory,
      border            = border,
      pos_sd            = pos_sd,
      vel_sd            = vel_sd,
      rand_vel_x        = rand_vel_x,
      rand_vel_y        = rand_vel_y,
      prior_c           = start_c,
      inv_c_prior_sigma = inv_c_prior_sigma,
      Lt                = Lt,
      n_threads         = n_threads
    )

    if (eval_counter %% print_every == 0) {

      # Truncate C array so it doesn't word-wrap in the console
      n_c = length(c)
      if (n_c > 6) {
        c_str = paste0(sprintf("%.3f, %.3f, %.3f", c[1], c[2], c[3]),
                       ", ... , ",
                       sprintf("%.3f, %.3f, %.3f", c[n_c-2], c[n_c-1], c[n_c]))
      } else {
        c_str = paste(sprintf("%.3f", c), collapse = ", ")
      }

      # Formatted output with trailing newline
      cat(sprintf("  Iter: %4d  |  Loss: %12.4f  |  C: [%s]\n",
                  eval_counter, loss, c_str))
    }

    return(loss)
  }

  # =================================================================
  # Step 7: Run the optimization
  # =================================================================

  opt_result = nloptr(
    x0 = start_c,
    eval_f = rMAP_loss,
    opts = list(
      "algorithm"   = "NLOPT_LN_NEWUOA",
      "ftol_rel"    = 1e-6,
      "maxeval"     = 5000,
      "print_level" = 0
    )
  )

  # Format the Final output identically
  final_c = opt_result$solution
  n_fc = length(final_c)
  if (n_fc > 6) {
    final_c_str = paste0(sprintf("%.3f, %.3f, %.3f", final_c[1], final_c[2], final_c[3]),
                         ", ... , ",
                         sprintf("%.3f, %.3f, %.3f", final_c[n_fc-2], final_c[n_fc-1], final_c[n_fc]))
  } else {
    final_c_str = paste(sprintf("%.3f", final_c), collapse = ", ")
  }

  # Add visual footer to close out the optimization block
  cat("\n-------------------------------------------------------\n")
  cat(sprintf("  FINAL Iter: %4d  |  Loss: %12.4f  |  C: [%s]\n",
              opt_result$iterations, opt_result$objective, final_c_str))
  cat("-------------------------------------------------------\n\n")

  # Step 8: Extract the optimized control points
  optimized_c = opt_result$solution

  opt_c_x = optimized_c[1:N_param]
  opt_c_y = optimized_c[1:N_param + N_param]

  # Eval Optimized Path

  opt_path_x = Phi_full %*% opt_c_x
  opt_path_y = Phi_full %*% opt_c_y

  opt_path_data_x = Phi_data_full %*% opt_c_x
  opt_path_data_y = Phi_data_full %*% opt_c_y

  opt_path_x_d = Phi_full_d %*% opt_c_x
  opt_path_y_d = Phi_full_d %*% opt_c_y

  opt_path_traj_vel = TrajWeightedBaseVectorFields_2D_Cosine(cbind(t_quad, opt_path_x, opt_path_y), sampledTrajectory, baseVectorFields_Vec, border)

  opt_path_traj_vel_x = opt_path_traj_vel[,1]
  opt_path_traj_vel_y = opt_path_traj_vel[,2]


  ## Positional Loss

  data_diff_x = data$X1 - opt_path_data_x
  data_diff_y = data$X2 - opt_path_data_y

  opt_pos_NLL = sum(data_diff_x^2 + data_diff_y^2) / (2 * pos_sd^2)

  ## Velocity Loss

  vel_diff_x = opt_path_x_d - opt_path_traj_vel_x
  vel_diff_y = opt_path_y_d - opt_path_traj_vel_y

  opt_vel_NLL = sum(vel_diff_x^2 + vel_diff_y^2) / (2 * vel_sd^2) * (Lt/N_quad)

  ## Prior Loss

  opt_prior_NLL = as.numeric((t(opt_c_x - start_c_x) %*% inv_c_prior_sigma %*% (opt_c_x - start_c_x) + t(opt_c_y - start_c_y) %*% inv_c_prior_sigma %*% (opt_c_y - start_c_y)) / 2)

  opt_NLL = opt_pos_NLL + opt_vel_NLL + opt_prior_NLL

  cat("\n")

  if(plot){

    ggplot() + geom_path(aes(x = opt_path_x, y = opt_path_y), size = 0.75) + geom_point(data = aug_data, aes(x = X1, y = X2), color = 'red')

  }

  list(Optimized_C = cbind(opt_c_x, opt_c_y), Optimized_Path = cbind(opt_path_x, opt_path_y), Posterior_NLL_Position = opt_pos_NLL, Posterior_NLL_Velocity = opt_vel_NLL, Posterior_NLL_Prior = opt_prior_NLL)

}

run_rMAP_Path = function(N_samples, data, pos_sd, vel_sd, pos_selection_sd, trajectoryPost, border, N_quad, baseVectorFields_Vec, full_nodes, N_prop_steps, print_every = 50, plot = F){

  N_full_nodes = length(full_nodes)

  C_X_Samples = matrix(nrow = N_samples, ncol = N_full_nodes+4)
  C_Y_Samples = matrix(nrow = N_samples, ncol = N_full_nodes+4)

  Path_X_Samples = matrix(nrow = N_samples, ncol = N_quad)
  Path_Y_Samples = matrix(nrow = N_samples, ncol = N_quad)

  NLL_Pos = rep(0, N_samples)
  NLL_Vel = rep(0, N_samples)
  NLL_Prior = rep(0, N_samples)

  for(i in 1:N_samples){

    cat(sprintf("=========== PROGRESS: Sample %d of %d ===========\n", i, N_samples))
    flush.console()

    cur_sample = find_one_rMAP_Path(data = data, pos_sd = pos_sd, vel_sd = vel_sd, pos_selection_sd = pos_selection_sd, trajectoryPost = trajectoryPost, border = c(-2,-2,2,2),
                                      N_quad = N_quad, baseVectorFields_Vec = baseVectorFields_Vec, full_nodes = full_nodes, N_prop_steps = N_prop_steps, print_every = print_every, plot = plot)

    C_X_Samples[i,] = cur_sample$Optimized_C[,1]
    C_Y_Samples[i,] = cur_sample$Optimized_C[,2]

    Path_X_Samples[i,] = cur_sample$Optimized_Path[,1]
    Path_Y_Samples[i,] = cur_sample$Optimized_Path[,2]

    NLL_Pos[i] = cur_sample$Posterior_NLL_Position
    NLL_Vel[i] = cur_sample$Posterior_NLL_Velocity
    NLL_Prior[i] = cur_sample$Posterior_NLL_Prior

  }

  list(C_Draws_X = C_X_Samples, C_Draws_Y = C_Y_Samples, Path_X_Samples = Path_X_Samples, Path_Y_Samples = Path_Y_Samples,
       NLL_Pos = NLL_Pos, NLL_Vel = NLL_Vel, NLL_Prior = NLL_Prior)

}

run_rMAP_Path_Parellel <- function(N_samples, data, pos_sd, vel_sd, pos_selection_sd,
                                  trajectoryPost, border, N_quad, baseVectorFields_Vec,
                                  full_nodes, N_prop_steps, num_cores = 8,
                                  cpp_code_string) { # Pass your cpp_code string here

  # 1. Set up the Windows cluster
  cat(sprintf("Setting up cluster with %d cores...\n", num_cores))
  cl <- makeCluster(num_cores)

  # 2a. Export LOCAL variables (from inside this function's arguments)
  clusterExport(cl, varlist = c(
    "cpp_code_string",
    "data",
    "pos_sd",
    "vel_sd",
    "pos_selection_sd",
    "trajectoryPost",
    "border",
    "N_quad",
    "baseVectorFields_Vec",
    "full_nodes",
    "N_prop_steps"
  ), envir = environment())

  # 2b. Export GLOBAL custom functions (from your main R script)
  clusterExport(cl, varlist = c(
    "find_one_rMAP_Path",
    "TrajWeightedBaseVectorFields_2D_Cosine",
    "evaluate2DCosine_fast"
    # Add ANY other custom functions find_one_rMAP_Path uses here!
  ), envir = .GlobalEnv)

  # 3. Initialize the workers (Load packages and compile C++)
  # This takes a few seconds but only happens ONCE when the cluster starts
  cat("Compiling C++ code on worker nodes...\n")
  clusterEvalQ(cl, {
    library(Rcpp)
    library(RcppArmadillo)
    library(plotly)
    library(stringr)
    library(MASS)
    library(ggplot2)
    library(minqa)
    library(nloptr)
    library(splines2)
    library(parallel)
    # Workers compile the C++ code so they have valid memory pointers
    sourceCpp(code = cpp_code_string)
  })

  # 4. Define the worker function with a debugger
  worker_func <- function(i) {
    tryCatch({

      # Your normal function call
      cur_sample = find_one_rMAP_Path(
        data = data, pos_sd = pos_sd, vel_sd = vel_sd,
        pos_selection_sd = pos_selection_sd, trajectoryPost = trajectoryPost,
        border = c(-2,-2,2,2), N_quad = N_quad,
        baseVectorFields_Vec = baseVectorFields_Vec, full_nodes = full_nodes,
        N_prop_steps = N_prop_steps, print_every = 99999999999999, plot = FALSE
      )

      list(
        C_X       = cur_sample$Optimized_C[,1],
        C_Y       = cur_sample$Optimized_C[,2],
        Path_X    = cur_sample$Optimized_Path[,1],
        Path_Y    = cur_sample$Optimized_Path[,2],
        NLL_Pos   = cur_sample$Posterior_NLL_Position,
        NLL_Vel   = cur_sample$Posterior_NLL_Velocity,
        NLL_Prior = cur_sample$Posterior_NLL_Prior
      )

    }, error = function(e) {
      # If it crashes, return the EXACT error and traceback
      paste("Worker Failed. Error:", e$message)
    })
  }

  # 5. Run the loop in parallel
  cat("Running optimizations...\n")
  pboptions(type = "timer") # Gives a nice estimated time remaining
  results <- pblapply(1:N_samples, worker_func, cl = cl)

  # 6. Shut down the cluster
  stopCluster(cl)

  # ========================================================
  # NEW: Intercept and print errors before unpacking!
  # ========================================================
  is_error <- sapply(results, is.character)
  if (any(is_error)) {
    cat("\n=========================================\n")
    cat("FATAL ERROR ON WORKER NODES DETECTED:\n")
    # Print the exact error message from the first failed node
    print(results[[which(is_error)[1]]])
    cat("=========================================\n")
    stop("Execution halted to prevent unpacking crash.")
  }

  closeAllConnections()

  # 7. Unpack and return (only runs if everything succeeded)
  list(
    C_Draws_X      = do.call(rbind, lapply(results, `[[`, "C_X")),
    C_Draws_Y      = do.call(rbind, lapply(results, `[[`, "C_Y")),
    Path_X_Samples = do.call(rbind, lapply(results, `[[`, "Path_X")),
    Path_Y_Samples = do.call(rbind, lapply(results, `[[`, "Path_Y")),
    NLL_Pos        = sapply(results, `[[`, "NLL_Pos"),
    NLL_Vel        = sapply(results, `[[`, "NLL_Vel"),
    NLL_Prior      = sapply(results, `[[`, "NLL_Prior")
  )
}

run_rMAP_Path_Multithread = function(N_samples, data, pos_sd, vel_sd, pos_selection_sd, trajectoryPost, border, N_quad, baseVectorFields_Vec, full_nodes, N_prop_steps, print_every = 50, plot = F, n_threads){

  N_full_nodes = length(full_nodes)

  C_X_Samples = matrix(nrow = N_samples, ncol = N_full_nodes+4)
  C_Y_Samples = matrix(nrow = N_samples, ncol = N_full_nodes+4)

  Path_X_Samples = matrix(nrow = N_samples, ncol = N_quad)
  Path_Y_Samples = matrix(nrow = N_samples, ncol = N_quad)

  NLL_Pos = rep(0, N_samples)
  NLL_Vel = rep(0, N_samples)
  NLL_Prior = rep(0, N_samples)

  for(i in 1:N_samples){

    cat(sprintf("=========== PROGRESS: Sample %d of %d ===========\n", i, N_samples))
    flush.console()

    cur_sample = find_one_rMAP_Path_Multithread(data = data, pos_sd = pos_sd, vel_sd = vel_sd, pos_selection_sd = pos_selection_sd, trajectoryPost = trajectoryPost, border = c(-2,-2,2,2),
                                    N_quad = N_quad, baseVectorFields_Vec = baseVectorFields_Vec, full_nodes = full_nodes, N_prop_steps = N_prop_steps, print_every = print_every, plot = plot, n_threads = n_threads)

    C_X_Samples[i,] = cur_sample$Optimized_C[,1]
    C_Y_Samples[i,] = cur_sample$Optimized_C[,2]

    Path_X_Samples[i,] = cur_sample$Optimized_Path[,1]
    Path_Y_Samples[i,] = cur_sample$Optimized_Path[,2]

    NLL_Pos[i] = cur_sample$Posterior_NLL_Position
    NLL_Vel[i] = cur_sample$Posterior_NLL_Velocity
    NLL_Prior[i] = cur_sample$Posterior_NLL_Prior

  }

  list(C_Draws_X = C_X_Samples, C_Draws_Y = C_Y_Samples, Path_X_Samples = Path_X_Samples, Path_Y_Samples = Path_Y_Samples,
       NLL_Pos = NLL_Pos, NLL_Vel = NLL_Vel, NLL_Prior = NLL_Prior)

}


## Testing ##

trajectoryPost = cbind(real_traj_test_1$Beta1Posterior, real_traj_test_1$Beta2Posterior)
data = sim_data_list[[1]]
pos_sd = 0.001
pos_selection_sd = 0.01
vel_sd = 0.001
N_prior_nodes = 10
prior_nodes = seq(0, max(data$t), length.out = N_prior_nodes+2)[-c(1,N_prior_nodes+2)]
N_full_nodes = 250
full_nodes = seq(0, max(data$t), length.out = N_full_nodes+2)[-c(1,N_full_nodes+2)]
N_quad = 1000
N_prop_steps = 100

N_samples = 100

t1 = Sys.time()

test_path_samples = run_rMAP_Path_Parellel(N_samples = N_samples, data = data, pos_sd = pos_sd, vel_sd = vel_sd, pos_selection_sd = pos_selection_sd,
                                              trajectoryPost = trajectoryPost, border = c(-2,-2,2,2), N_quad = N_quad,
                                              baseVectorFields_Vec = baseVectorFields_Vec, full_nodes = full_nodes, N_prop_steps = N_prop_steps, num_cores = 32, cpp_code_string = cpp_code)
t2 = Sys.time()

test_path_samples = run_rMAP_Path_Multithread(N_samples = N_samples, data = data, pos_sd = pos_sd, vel_sd = vel_sd, pos_selection_sd = pos_selection_sd,
                                        trajectoryPost = trajectoryPost, border = c(-2,-2,2,2), N_quad = N_quad,
                                        baseVectorFields_Vec = baseVectorFields_Vec, full_nodes = full_nodes, N_prop_steps = N_prop_steps, print_every = 500, plot = T, n_threads = 1)

t3 = Sys.time()

t3-t2

total_NLL = test_path_samples$NLL_Pos + test_path_samples$NLL_Vel + test_path_samples$NLL_Prior

path_samples_plotting = data.frame(t = rep(seq(min(data$t),max(data$t), length.out = N_quad), N_samples), X1 = c(t(test_path_samples$Path_X_Samples)), X2 = c(t(test_path_samples$Path_Y_Samples)), Sample = rep(1:N_samples, each = N_quad))

ggplot(data = path_samples_plotting, aes(x = t, y = X1, group = Sample)) + geom_path(alpha = 0.2)

?prcomp

PCA_Path_Samples = prcomp(cbind(test_path_samples$C_Draws_X, test_path_samples$C_Draws_Y), center = T, scale. = F)
PCA_Path_Summary = summary(PCA_Path_Samples)
N_comp = 2

ggplot(PCA_Path_Samples$x, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = log(test_path_samples$NLL_Pos+test_path_samples$NLL_Vel+test_path_samples$NLL_Prior))) +
  labs(color = "NLL")


## K-Means on the top principle components

library(cluster)

N_important_PC = as.numeric(which(PCA_Path_Summary$importance[3,] > 0.99)[1])

PCA_sub = PCA_Path_Samples$x[,1:(N_important_PC)]
W_pca = PCA_Path_Samples$rotation[,1:(N_important_PC)]
mu_pca = PCA_Path_Samples$center

run_k_means = function(PCA_sub, max_k){

  ## Find optimal k

  gap_stat <- clusGap(PCA_sub, FUNcluster = kmeans, K.max = max_k, B = 1000)
  optimal_K <- maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"], method = "firstSEmax")

  opt_k_means = kmeans(PCA_sub, centers = optimal_K)

  list(K = optimal_K, K_means = opt_k_means)
}

path_samples = test_path_samples
PC_importance_threshold = 0.99
max_clusters = 10

## Get the starting paths for finding the true modes

get_starting_mode_paths = function(path_samples, PC_importance_threshold, max_clusters){

  N_coefs = ncol(path_samples$C_Draws_X)

  path_coef_pca = prcomp(cbind(path_samples$C_Draws_X, path_samples$C_Draws_Y), center = T, scale. = F)
  pca_summary = summary(path_coef_pca)

  N_PCs = as.numeric(which(pca_summary$importance[3,] > 0.99)[1])

  PCA_sub = path_coef_pca$x[,1:(N_PCs)]
  W_pca = path_coef_pca$rotation[,1:(N_PCs)]
  mu_pca = path_coef_pca$center

  k_means = run_k_means(PCA_sub = PCA_sub, max_k = max_clusters)

  K_means_centers = k_means$K_means$centers
  K_means_K = k_means$K

  PCA_center_c = K_means_centers %*% t(W_pca) + matrix(rep(mu_pca, K_means_K), nrow = K_means_K, byrow = T)

  PCA_center_c_x = PCA_center_c[,1:N_coefs]
  PCA_center_c_y = PCA_center_c[,1:N_coefs + N_coefs]

  list(Center_c_x = PCA_center_c_x, Center_c_y = PCA_center_c_y)

}

Mode_Opt_Paths = get_starting_mode_paths(path_samples = test_path_samples, PC_importance_threshold = 0.99, max_clusters = 10)

t_quad = seq(min(data$t), max(data$t), length.out = N_quad)
Phi_full = bSpline(t_quad, knots = full_nodes, intercept = T)

PCA_center_c_x = Mode_Opt_Paths$Center_c_x
PCA_center_c_y = Mode_Opt_Paths$Center_c_y

PCA_center_path_x = Phi_full %*% t(PCA_center_c_x)
PCA_center_path_y = Phi_full %*% t(PCA_center_c_y)

K_means_K = nrow(PCA_center_c_x)

PCA_center_plotting_df = data.frame(t = rep(t_quad, K_means_K), X1 = c(PCA_center_path_x), X2 = c(PCA_center_path_y), Center = rep(1:3, each = N_quad))

ggplot(data = PCA_center_plotting_df, aes(x = X1, y = X2, group = Center)) + geom_path()

## Optimize from the pca found centers

center_c_start = cbind(PCA_center_c_x[1,], PCA_center_c_y[2,])

find_one_mode_location_path = function(center_c_start, data, pos_sd, vel_sd, pos_selection_sd, trajectoryMode, border, N_quad, baseVectorFields_Vec, full_nodes, print_every, plot){

  center_c_start_x = center_c_start[,1]
  center_c_start_y = center_c_start[,2]

  center_c_start_full = c(center_c_start)

  # Step 2: Augment Data with positional error
  # data is Nx3 with columns t, X1, X2

  N = nrow(data)

  # Step 3: Get prior path using a bSpline regression with prior_nodes

  t_quad = seq(min(data$t), max(data$t), length.out = N_quad)

  # Step 4: Get full path basis to add onto prior with full nodes

  Phi_data_full = bSpline(data$t, knots = full_nodes, intercept = T)
  Phi_full = bSpline(t_quad, knots = full_nodes, intercept = T)
  Phi_full_d = bSpline(t_quad, knots = full_nodes, intercept = T, derivs = 1)
  N_param = ncol(Phi_full)

  # Step 5: Sample random initial path
  Lt = max(data$t) - min(data$t)

  helper_M = t(Phi_full) %*% Phi_full
  helper_M_inv = ginv(helper_M)

  c_prior_sigma = (N_quad / Lt * pos_selection_sd^2) * helper_M_inv
  inv_c_prior_sigma = (Lt / (N_quad * pos_selection_sd^2)) * helper_M

  # =================================================================
  # Step 6: Define loss function
  # =================================================================

  aug_X_matrix = as.matrix(data[, c("X1", "X2")])
  eval_counter <- 0

  # Clean visual header for a new optimization run
  cat("\n=======================================================\n")
  cat("       Starting Optimization for New Mode            \n")
  cat("=======================================================\n\n")

  rMAP_loss = function(c) {

    eval_counter <<- eval_counter + 1

    # Call the compiled C++ function
    loss = spline_rMAP_loss_cpp(
      c_params          = c,
      Phi_data          = Phi_data_full,
      Phi_quad          = Phi_full,
      Phi_deriv_quad    = Phi_full_d,
      aug_X             = aug_X_matrix,
      sampledTrajectory = trajectoryMode,
      border            = border,
      pos_sd            = pos_sd,
      vel_sd            = vel_sd,
      rand_vel_x        = rep(0,N_quad),
      rand_vel_y        = rep(0,N_quad),
      prior_c           = center_c_start_full,
      inv_c_prior_sigma = inv_c_prior_sigma,
      Lt                = Lt
    )

    if (eval_counter %% print_every == 0) {

      # Truncate C array so it doesn't word-wrap in the console
      n_c = length(c)
      if (n_c > 6) {
        c_str = paste0(sprintf("%.3f, %.3f, %.3f", c[1], c[2], c[3]),
                       ", ... , ",
                       sprintf("%.3f, %.3f, %.3f", c[n_c-2], c[n_c-1], c[n_c]))
      } else {
        c_str = paste(sprintf("%.3f", c), collapse = ", ")
      }

      # Formatted output with trailing newline
      cat(sprintf("  Iter: %4d  |  Loss: %12.4f  |  C: [%s]\n",
                  eval_counter, loss, c_str))
    }

    return(loss)
  }

  # =================================================================
  # Step 7: Run the optimization
  # =================================================================

  opt_result = nloptr(
    x0 = center_c_start_full,
    eval_f = rMAP_loss,
    opts = list(
      "algorithm"   = "NLOPT_LN_NEWUOA",
      "ftol_rel"    = 1e-6,
      "maxeval"     = 2000,
      "print_level" = 0
    )
  )

  # Format the Final output identically
  final_c = opt_result$solution
  n_fc = length(final_c)
  if (n_fc > 6) {
    final_c_str = paste0(sprintf("%.3f, %.3f, %.3f", final_c[1], final_c[2], final_c[3]),
                         ", ... , ",
                         sprintf("%.3f, %.3f, %.3f", final_c[n_fc-2], final_c[n_fc-1], final_c[n_fc]))
  } else {
    final_c_str = paste(sprintf("%.3f", final_c), collapse = ", ")
  }

  # Add visual footer to close out the optimization block
  cat("\n-------------------------------------------------------\n")
  cat(sprintf("  FINAL Iter: %4d  |  Loss: %12.4f  |  C: [%s]\n",
              opt_result$iterations, opt_result$objective, final_c_str))
  cat("-------------------------------------------------------\n\n")

  # Step 8: Extract the optimized control points
  optimized_c = opt_result$solution

  opt_c_x = optimized_c[1:N_param]
  opt_c_y = optimized_c[1:N_param + N_param]

  # Eval Optimized Path

  opt_path_x = Phi_full %*% opt_c_x
  opt_path_y = Phi_full %*% opt_c_y

  opt_path_data_x = Phi_data_full %*% opt_c_x
  opt_path_data_y = Phi_data_full %*% opt_c_y

  opt_path_x_d = Phi_full_d %*% opt_c_x
  opt_path_y_d = Phi_full_d %*% opt_c_y

  opt_path_traj_vel = TrajWeightedBaseVectorFields_2D_Cosine(cbind(t_quad, opt_path_x, opt_path_y), trajectoryMode, baseVectorFields_Vec, border)

  opt_path_traj_vel_x = opt_path_traj_vel[,1]
  opt_path_traj_vel_y = opt_path_traj_vel[,2]


  ## Positional Loss

  data_diff_x = data$X1 - opt_path_data_x
  data_diff_y = data$X2 - opt_path_data_y

  opt_pos_NLL = sum(data_diff_x^2 + data_diff_y^2) / (2 * pos_sd^2)

  ## Velocity Loss

  vel_diff_x = opt_path_x_d - opt_path_traj_vel_x
  vel_diff_y = opt_path_y_d - opt_path_traj_vel_y

  opt_vel_NLL = sum(vel_diff_x^2 + vel_diff_y^2) / (2 * vel_sd^2) * (Lt/N_quad)

  ## Prior Loss

  opt_prior_NLL = as.numeric((t(opt_c_x - center_c_start_x) %*% inv_c_prior_sigma %*% (opt_c_x - center_c_start_x) + t(opt_c_y - center_c_start_y) %*% inv_c_prior_sigma %*% (opt_c_y - center_c_start_y)) / 2)

  opt_NLL = opt_pos_NLL + opt_vel_NLL + opt_prior_NLL

  cat("\n")

  if(plot){

    mode_plot = ggplot() + geom_path(aes(x = opt_path_x, y = opt_path_y), size = 0.75) + geom_point(data = data, aes(x = X1, y = X2), color = 'red')

  } else{
    mode_plot = NULL
  }

  return(list(Optimized_C = cbind(opt_c_x, opt_c_y), Optimized_Path = cbind(opt_path_x, opt_path_y), Posterior_NLL_Position = opt_pos_NLL, Posterior_NLL_Velocity = opt_vel_NLL, Posterior_NLL_Prior = opt_prior_NLL, Plot = mode_plot))

}

center_c_start_mat_x = Mode_Opt_Paths$Center_c_x
center_c_start_mat_y = Mode_Opt_Paths$Center_c_y

trajectoryModes = list(matrix(trajectoryPost[37,], ncol = 2, byrow = F), matrix(trajectoryPost[56,], ncol = 2, byrow = F), matrix(trajectoryPost[100,], ncol = 2, byrow = F))

get_posterior_path_modes = function(center_c_start_mat_x, center_c_start_mat_y, data, pos_sd, vel_sd, pos_selection_sd, trajectoryModes, border, N_quad, baseVectorFields_Vec, full_nodes, print_every, plot){

  #Define matrices for mode coefficients, mode path evaluations, NLL evaluations

  N_path_modes = nrow(center_c_start_mat_x)
  N_traj_modes = length(trajectoryModes)
  N_combos = N_path_modes * N_traj_modes
  N_c = ncol(center_c_start_mat_x)

  path_mode_c_x = matrix(nrow = N_combos, ncol = N_c)
  path_mode_c_y = matrix(nrow = N_combos, ncol = N_c)

  path_mode_eval_x = matrix(nrow = N_combos, ncol = N_quad)
  path_mode_eval_y = matrix(nrow = N_combos, ncol = N_quad)

  path_mode_NLL_pos = rep(0, N_combos)
  path_mode_NLL_vel = rep(0, N_combos)
  path_mode_NLL_prior = rep(0, N_combos)

  if(plot){
    path_mode_plots = list()
  } else{
    path_mode_plots = NULL
  }

  #Loops through all mode starting locations

  for(i in 1:N_traj_modes){

    cur_trajMode = trajectoryModes[[i]]

    for(j in 1:N_path_modes){

      cur_path_start = cbind(center_c_start_mat_x[j,], center_c_start_mat_y[j,])

      cur_pathMode = find_mode_location(center_c_start = cur_path_start, data = data, pos_sd = pos_sd, vel_sd = vel_sd, pos_selection_sd = pos_selection_sd, trajectoryMode = cur_trajMode, border = border, N_quad = N_quad, baseVectorFields_Vec = baseVectorFields_Vec, full_nodes = full_nodes, print_every = print_every, plot = plot)

      path_mode_c_x[(N_traj_modes)*(i-1) + j,] = cur_pathMode$Optimized_C[,1]
      path_mode_c_y[(N_traj_modes)*(i-1) + j,] = cur_pathMode$Optimized_C[,2]

      path_mode_eval_x[(N_traj_modes)*(i-1) + j,] = cur_pathMode$Optimized_Path[,1]
      path_mode_eval_y[(N_traj_modes)*(i-1) + j,] = cur_pathMode$Optimized_Path[,2]

      path_mode_NLL_pos[(N_traj_modes)*(i-1) + j] = cur_pathMode$Posterior_NLL_Position
      path_mode_NLL_vel[(N_traj_modes)*(i-1) + j] = cur_pathMode$Posterior_NLL_Velocity
      path_mode_NLL_prior[(N_traj_modes)*(i-1) + j] = cur_pathMode$Posterior_NLL_Prior

      if(plot){
        path_mode_plots[[(N_traj_modes)*(i-1) + j]] = cur_pathMode$Plot
      }

    }

  }

  return(list(C_X = path_mode_c_x, C_Y = path_mode_c_y, Path_X = path_mode_eval_x, Path_Y = path_mode_eval_y, NLL_Pos = path_mode_NLL_pos, NLL_Vel = path_mode_NLL_vel, NLL_Prior = path_mode_NLL_prior, Plots = path_mode_plots))

}

test_posterior_path_modes = get_posterior_path_modes(center_c_start_mat_x = center_c_start_mat_x, center_c_start_mat_y = center_c_start_mat_y,
                                                     data = data, pos_sd = pos_sd, vel_sd = vel_sd, pos_selection_sd = pos_selection_sd, trajectoryModes = trajectoryModes, border = border, N_quad = N_quad, baseVectorFields_Vec = baseVectorFields_Vec, full_nodes = full_nodes, print_every = print_every, plot = plot)


