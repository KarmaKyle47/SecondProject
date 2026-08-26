library(plotly)
library(stringr)
library(MASS)
library(ggplot2)
library(minqa)
library(nloptr)


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

sample_rMAP_trajectories = function(sim_data_list, pos_sd = 0.001, vel_sd = 0.1, M = 2, prior_k = 0.35, prior_l = 1, baseVectorFields_Vec, border, N_prop_steps = 10, traj_eval_grid, print_every = 50){

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

  prior_beta_sigma = diag(c(prior_k^2 * diag(spec_den) %*% matrix(rep(1,4), nrow = 2) %*% diag(spec_den)))

  upper_beta = rep(sqrt(diag(prior_beta_sigma))*3,2)
  lower_beta = -1*upper_beta

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
    x0 = start_beta,               # Starting values
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

  post_likelihood_loss_pos <- sum(diff_1^2 + diff_2^2) / (2 * (N_advects) * pos_sd^2)

  d_beta1 <- post_beta_1 - start_beta_1
  d_beta2 <- post_beta_2 - start_beta_2

  post_prior_loss <- crossprod(d_beta1, prior_precision_mat %*% d_beta1) +
    crossprod(d_beta2, prior_precision_mat %*% d_beta2)

  cat("\n")

  list(post_beta_mat, post_likelihood_loss_pos, post_prior_loss, post_traj_grid)

}

run_rMAP_Trajectory = function(N_samples, sim_data_list, pos_sd, vel_sd, M, prior_k, prior_l, baseVectorFields_Vec, border, N_prop_steps, traj_eval_grid, print_every = 50){

  beta_1_post_mat = matrix(nrow = N_samples, ncol = M^2)
  beta_2_post_mat = matrix(nrow = N_samples, ncol = M^2)

  post_traj_eval_1_list = list()
  post_traj_eval_2_list = list()

  post_like_pos = rep(0, N_samples)
  post_like_prior = rep(0, N_samples)
  accept = rep(0, N_samples)

  for(i in 1:N_samples){

    cat(sprintf("=========== PROGRESS: Sample %d of %d ===========\n", i, N_samples))
    flush.console()

    cur_draw = sample_rMAP_trajectories(sim_data_list = sim_data_list, pos_sd = pos_sd, vel_sd = vel_sd, M = M, prior_k = prior_k, prior_l = prior_l, baseVectorFields_Vec = baseVectorFields_Vec, border = border, N_prop_steps = N_prop_steps, traj_eval_grid = traj_eval_grid, print_every = print_every)

    # Metropolis
    if(i > 1){

      prev_NLL = post_like_pos[i-1] + post_like_prior[i-1]
      cur_NLL = cur_draw[[2]] + cur_draw[[3]]

      log_acc_prob = -1* (cur_NLL - prev_NLL)

      u = runif(1)

      if(log(u) < log_acc_prob){

        beta_1_post_mat[i,] = cur_draw[[1]][,1]
        beta_2_post_mat[i,] = cur_draw[[1]][,2]

        post_like_pos[i] = cur_draw[[2]]
        post_like_prior[i] = cur_draw[[3]]

        post_traj_eval_1_list[[i]] = cur_draw[[4]][,1]
        post_traj_eval_2_list[[i]] = cur_draw[[4]][,2]

        accept[i] = 1

        cat('\n*** Accepted ***\n\n')

      } else{

        beta_1_post_mat[i,] = beta_1_post_mat[i-1,]
        beta_2_post_mat[i,] = beta_2_post_mat[i-1,]

        post_like_pos[i] = post_like_pos[i-1]
        post_like_prior[i] = post_like_prior[i-1]

        post_traj_eval_1_list[[i]] = post_traj_eval_1_list[[i-1]]
        post_traj_eval_2_list[[i]] = post_traj_eval_2_list[[i-1]]

        cat('\n*** Rejected ***\n\n')

      }

    } else{

      beta_1_post_mat[i,] = cur_draw[[1]][,1]
      beta_2_post_mat[i,] = cur_draw[[1]][,2]

      post_like_pos[i] = cur_draw[[2]]
      post_like_prior[i] = cur_draw[[3]]

      post_traj_eval_1_list[[i]] = cur_draw[[4]][,1]
      post_traj_eval_2_list[[i]] = cur_draw[[4]][,2]

      accept[i] = 1

      cat('\n*** Initial Accepted ***\n\n')

    }

  }

  post_draws_traj_eval_1 = do.call(cbind, post_traj_eval_1_list)
  post_draws_traj_eval_2 = do.call(cbind, post_traj_eval_2_list)

  list(Beta1Posterior = beta_1_post_mat, Beta2Posterior = beta_2_post_mat,
       Traj1Posterior = post_draws_traj_eval_1, Traj2Posterior = post_draws_traj_eval_2,
       PosPosteriorNLL = post_like_pos, PriorPosteriorNLL = post_like_prior, MH_Acceptance = accept)

}

### Tests

M=2
true_l = 1
true_k = 0.35

trueHSGP = sampleFullTrajectoriesHSGP(2, M = M-1, log_k = true_k, log_l = true_l)

omega = (1:M-1)*pi
spec_den = sqrt(2*pi)*true_l*exp(-0.5*true_l^2*omega^2)

prior_beta_sigma = diag(c(true_k^2 * diag(spec_den) %*% matrix(rep(1,M^2), nrow = M) %*% diag(spec_den)))

true_beta_1 = c(true_k*diag(sqrt(spec_den)) %*% trueHSGP$log_z[[1]] %*% diag(sqrt(spec_den)))
true_beta_2 = c(true_k*diag(sqrt(spec_den)) %*% trueHSGP$log_z[[2]] %*% diag(sqrt(spec_den)))

true_beta_mat = cbind(true_beta_1, true_beta_2)


sampledParticles = samplePhySpaceParticles(n_particles = 100, startTime = 0, n_obs = 20*10, border = c(-2,-2,2,2), borderBuffer = 0.2, baseVectorFields, trueHSGP,
                                           M = 2, t_step_mean = 0.001, vel_sigma = 0, pos_sigma = 0)

sampledParticles_Sub = sampledParticles[1:(100*20) * 10 - (10-1),]
ggplot(sampledParticles_Sub, aes(x = X1, y = X2, color = Particle)) + geom_point() + theme(legend.position = 'None')

sim_data_list = list()

for(d in 1:100){

  curDrifter = sampledParticles_Sub[sampledParticles_Sub$Particle == str_c('Particle',d),]

  sim_data_list[[d]] = curDrifter[,c(1,2,3)]


}


traj_eval_grid = expand.grid(seq(-2,2, length.out = 100), seq(-2,2, length.out = 100))

TrueTrajVals = exp(evaluate2DCosine_fast(true_beta_mat, pos_mat = traj_eval_grid, border = c(-2,-2,2,2)))


real_traj_test_1 = run_rMAP_Trajectory(N_samples = 100, sim_data_list = sim_data_list,
                                       pos_sd = 0.001, vel_sd = 0, M = 2, prior_k = 0.35,
                                       prior_l = 1, baseVectorFields_Vec = baseVectorFields_Vec,
                                       border = c(-2,-2,2,2), N_prop_steps = 5, traj_eval_grid = traj_eval_grid, print_every = 50)
??nloptr

mean(real_traj_test_1$MH_Acceptance)

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


hist(real_traj_test_1$Beta1Posterior[,1])

colMeans(real_traj_test_1$Beta1Posterior)
colMeans(real_traj_test_1$Beta2Posterior)


