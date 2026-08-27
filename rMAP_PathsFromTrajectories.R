library(plotly)
library(stringr)
library(MASS)
library(ggplot2)
library(minqa)
library(nloptr)
library(splines2)


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
data = sim_data_list[[18]]
pos_sd = 0.001
pos_selection_sd = 0.001
vel_sd = 0.001

prior_nodes = seq(0, max(data$t), length.out = 6)[-c(1,12)]
full_nodes = seq(0, max(data$t), length.out = 22)[-c(1,52)]
N_quad = 1000

find_rMAP_Path_Modes = function(data, pos_sd, vel_sd, pos_selection_sd, trajectoryPost, border, N_quad, baseVectorFields_Vec, prior_nodes, full_nodes, print_every, plot){

  # Step 1: Sample Trajectory - trajectory Post is a data.frame or matrix with each posterior draw as a row

  sampledTrajectory = matrix(trajectoryPost[sample(1:nrow(trajectoryPost), size = 1),], ncol = 2, byrow = F)

  # Step 2: Augment Data with positional error
  # data is Nx3 with columns t, X1, X2

  N = nrow(data)
  aug_data = data + matrix(c(rep(0, N), rnorm(N, 0, pos_sd), rnorm(N, 0, pos_sd)), byrow = F, ncol = 3)

  # Step 3: Get prior path using a bSpline regression with prior_nodes

  t_quad = seq(min(aug_data$t), max(aug_data$t), length.out = N_quad)

  Phi_data_prior = bSpline(aug_data$t, knots = prior_nodes, intercept = T)
  Phi_prior = bSpline(t_quad, knots = prior_nodes, intercept = T)
  Phi_prior_d = bSpline(t_quad, knots = prior_nodes, intercept = T, derivs = 1)

  Phi_data_prior_H = solve(t(Phi_data_prior) %*% Phi_data_prior) %*% t(Phi_data_prior)
  prior_c_x = Phi_data_prior_H %*% aug_data$X1
  prior_c_y = Phi_data_prior_H %*% aug_data$X2

  # Step 4: Get full path basis to add onto prior with full nodes

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

  # Step 5.5: Project prior_c to the full dimension

  full_Proj_M = helper_M_inv %*% t(Phi_full) %*% Phi_prior

  prior_full_c_x = full_Proj_M %*% prior_c_x
  prior_full_c_y = full_Proj_M %*% prior_c_y

  start_c_x = prior_full_c_x + mvrnorm(mu = rep(0, N_param), Sigma = c_prior_sigma)
  start_c_y = prior_full_c_y + mvrnorm(mu = rep(0, N_param), Sigma = c_prior_sigma)

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

  opt_path_traj_vel = TrajWeightedBaseVectorFields_2D_Cosine(cbind(t_quad, opt_path_x, opt_path_y), sampledTrajectory, baseVectorFields_Vec, border) + cbind(rand_vel_x, rand_vel_y)

  opt_path_traj_vel_x = opt_path_traj_vel[,1]
  opt_path_traj_vel_y = opt_path_traj_vel[,2]


  ## Positional Loss

  data_diff_x = aug_data$X1 - opt_path_data_x
  data_diff_y = aug_data$X2 - opt_path_data_y

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

run_rMAP_Trajectory = function(N_samples, data, pos_sd, vel_sd, pos_selection_sd, trajectoryPost, border, N_quad, baseVectorFields_Vec, prior_nodes, full_nodes, print_every = 50, plot = F){

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

    cur_sample = find_rMAP_Path_Modes(data = data, pos_sd = pos_sd, vel_sd = vel_sd, pos_selection_sd = pos_selection_sd, trajectoryPost = trajectoryPost, border = c(-2,-2,2,2),
                                      N_quad = N_quad, baseVectorFields_Vec = baseVectorFields_Vec, prior_nodes = prior_nodes, full_nodes = full_nodes, print_every = print_every, plot = T)

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
data = sim_data_list[[18]]
pos_sd = 0.001
pos_selection_sd = 0.001
vel_sd = 0.001
N_prior_nodes = 4
prior_nodes = seq(0, max(data$t), length.out = N_prior_nodes+2)[-c(1,N_prior_nodes+2)]
N_full_nodes = 20
full_nodes = seq(0, max(data$t), length.out = N_full_nodes+2)[-c(1,N_full_nodes+2)]
N_quad = 1000

N_samples = 100

test_path_samples = run_rMAP_Trajectory(N_samples = 100, data = data, pos_sd = pos_sd, vel_sd = vel_sd, pos_selection_sd = pos_selection_sd,
                                        trajectoryPost = trajectoryPost, border = c(-2,-2,2,2), N_quad = N_quad,
                                        baseVectorFields_Vec = baseVectorFields_Vec, prior_nodes = prior_nodes, full_nodes = full_nodes, print_every = 500, plot = T)



PCA_Path_Samples = prcomp(cbind(test_path_samples$C_Draws_X, test_path_samples$C_Draws_Y), center = T, scale. = T)

N_comp = 2

ggplot(PCA_Path_Samples$x[,1:N_comp], aes(x = PC1, y = PC2)) + geom_point(aes(color = log(test_path_samples$NLL_Pos+test_path_samples$NLL_Vel+test_path_samples$NLL_Prior)))



