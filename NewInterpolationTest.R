source('SecondProjectRequiredFunctions.R')
options(mc.cores = parallel::detectCores())
max_cores = parallel::detectCores()


baseVectorFields = function(t, curPos){

  window = cos(pi*curPos[1]/20) * cos(pi*curPos[2]/20)

  f1 = c(curPos[2],-1*curPos[1]) * window
  f2 = c(curPos[1],curPos[2])*sin(t/8) * window
  matrix(c(f1,f2), nrow = 2, byrow = F)

}

true_log_k = 0.35
true_log_l = 0.05
M = 10

sampledHSGP = sampleFullTrajectoriesHSGP(2, M, true_log_k, true_log_l)

rand_res = 100

sampledParticles_V2 = samplePhySpaceParticles(100, startTime = 0, n_obs = 100*rand_res, border = c(-10,-10,10,10), borderBuffer = 0.2, baseVectorFields, sampledHSGP,
                                           M = M, t_step_mean = 0.01, vel_sigma = 0, pos_sigma = 0)

sampledParticles_Sub = sampledParticles[1:10000 * rand_res - (rand_res-1),]

sampledParticles_Sub = sampledParticles_Sub[sampledParticles_Sub$X1 >= -10 & sampledParticles_Sub$X1 <= 10 &
                                              sampledParticles_Sub$X2 >= -10 & sampledParticles_Sub$X2 <= 10,]

sampledParticles_Sub = sim_data_sub

ggplot(sampledParticles_Sub, aes(x = X1, y = X2, color = Particle)) + geom_point()


drifter_names = unique(sampledParticles_Sub$Particle)
drifter_num_points = c()
i=1
for(i in 1:length(drifter_names)){

  cur_drifter = sampledParticles_Sub[sampledParticles_Sub$Particle == drifter_names[i],]

  drifter_num_points = c(drifter_num_points, nrow(cur_drifter))

}

N_quad = 1000

D = length(drifter_names)
N = drifter_num_points
M_drifter = 20
K = N+M_drifter

c_part_x_list = list()
c_part_y_list = list()
Z_flat_full_list = list()
Phi_List = list()
Phi_d_List = list()
t_grid_list = list()
Lts = rep(0,D)
drifter_boundaries = matrix(nrow = D, ncol = 2)
d=1
for(d in 1:D){

  curDrifter = sampledParticles_Sub[sampledParticles_Sub$Particle == drifter_names[d],]
  cur_boundary = c(min(curDrifter$t),max(curDrifter$t))
  cur_quad_grid = seq(cur_boundary[1], cur_boundary[2], length.out = N_quad)
  cur_K = K[d]

  cur_P = cos(outer((curDrifter$t - cur_boundary[1])/diff(cur_boundary), pi * (1:cur_K - 1)))

  cur_Z = Null(t(cur_P))[,1:M_drifter]

  cur_Phi = cos(outer((cur_quad_grid - cur_boundary[1])/diff(cur_boundary), pi * (1:cur_K - 1)))
  cur_Phi_d = sin(outer((cur_quad_grid - cur_boundary[1])/diff(cur_boundary), pi * (1:cur_K - 1))) * outer(rep(1, N_quad), -1*(pi*(1:cur_K - 1)/diff(cur_boundary)))

  LI_x = spline(x = curDrifter$t, y = curDrifter$X1, xout = cur_quad_grid)$y
  LI_y = spline(x = curDrifter$t, y = curDrifter$X2, xout = cur_quad_grid)$y

  cur_c_part_x = solve(t(cur_Phi) %*% cur_Phi) %*% t(cur_Phi) %*% LI_x
  cur_c_part_y = solve(t(cur_Phi) %*% cur_Phi) %*% t(cur_Phi) %*% LI_y

  c_part_x_list[[d]] = cur_c_part_x
  c_part_y_list[[d]] = cur_c_part_y

  Z_flat_full_list[[d]] = c(cur_Z)

  Phi_List[[d]] = cur_Phi
  Phi_d_List[[d]] = cur_Phi_d

  t_grid_list[[d]] = cur_quad_grid
  Lts[d] = diff(cur_boundary)

  drifter_boundaries[d,] = cur_boundary

  svMisc::progress(d, D)

}
d=1
plotPriorDrifter = function(d){

  curDrifter = sampledParticles_Sub[sampledParticles_Sub$Particle == drifter_names[d],]

  prior_pos_x = as.numeric(Phi_List[[d]] %*% c_part_x_list[[d]])
  prior_pos_y = as.numeric(Phi_List[[d]] %*% c_part_y_list[[d]])

  ggplot() + geom_line(aes(x = t_grid_list[[d]], y = prior_pos_x), color = 'blue') +
             geom_point(data = curDrifter, aes(x = t, y = X1), color = 'black', size = 1)

  ggplot() + geom_line(aes(x = t_grid_list[[d]], y = prior_pos_y), color = 'blue') +
             geom_point(data = curDrifter, aes(x = t, y = X2), color = 'black', size = 0.1)

  ggplot() + geom_path(aes(x = prior_pos_x, y = prior_pos_y), color = 'blue') +
    geom_point(data = curDrifter, aes(x = X1, y = X2), color = 'black', size = 1)
}

plotPriorDrifter(1)

unlist(lapply(Z_flat_full_list, length)) == K*M_drifter
length(Z_flat_full_list[[1]])
do.call(, Z_flat_full_list)

unlist(c_part_x_list)

full_data = list(
  D = D,
  N_quad = N_quad,

  M = rep(M_drifter, D),
  total_M = M_drifter*D,

  K = K,
  total_K = sum(K),
  total_Z = sum(K*M_drifter),

  c_part_x_flat = unlist(c_part_x_list),
  c_part_y_flat = unlist(c_part_y_list),

  drifter_boundaries = drifter_boundaries,
  drifter_times = sampledParticles_Sub$t,

  N_models = 2,
  M_Surface = 5,

  border = c(-12,-12,12,12),
  fixed_ks = c(0.35,0.35)#,
  #fixed_ls = c(0.05, 0.05)#,

  # size_hidden = 10,
  # W1 = list(GeoSubset_emulator$weights[[1]], EkmanSubset_emulator$weights[[1]]),
  # W2 = list(GeoSubset_emulator$weights[[3]], EkmanSubset_emulator$weights[[3]]),
  # #W3 = list(GeoSubset_emulator$weights[[5]], EkmanSubset_emulator$weights[[5]]),
  # #W4 = list(GeoSubset_emulator$weights[[7]], EkmanSubset_emulator$weights[[7]]),
  # #W5 = list(GeoSubset_emulator$weights[[9]], EkmanSubset_emulator$weights[[9]]),
  # W_out = list(GeoSubset_emulator$weights[[5]], EkmanSubset_emulator$weights[[5]]),
  #
  # B1 = list(GeoSubset_emulator$weights[[2]], EkmanSubset_emulator$weights[[2]]),
  # B2 = list(GeoSubset_emulator$weights[[4]], EkmanSubset_emulator$weights[[4]]),
  # #B3 = list(GeoSubset_emulator$weights[[6]], EkmanSubset_emulator$weights[[6]]),
  # #B4 = list(GeoSubset_emulator$weights[[8]], EkmanSubset_emulator$weights[[8]]),
  # #B5 = list(GeoSubset_emulator$weights[[10]], EkmanSubset_emulator$weights[[10]]),
  # B_out = list(GeoSubset_emulator$weights[[6]], EkmanSubset_emulator$weights[[6]]),
  #
  # VF_means = prior_means,
  # VF_sds = prior_sds

)

new_interpolator_model = cmdstan_model('NewInterpolationModel.stan', cpp_options = list(stan_threads = T))

full_fit = new_interpolator_model$variational(
  data = full_data,
  threads = max_cores/4,
  iter = 10000,
  draws = 1000, show_messages = T, init = 0#, algorithm = "fullrank"
)

PF_sum = pure_fit$summary()

PF_sum[PF_sum$variable == 'fixed_ls',]

surface_beta_draws = full_fit$draws('surface_betas')
w_x_draws = full_fit$draws('w_x')
w_y_draws = full_fit$draws('w_y')
sigma_vel_draws = full_fit$draws('sigma_vel')

c(full_fit$draws('post_vel_mse'))
c(full_fit$draws('prior_vel_mse'))

save(sampledHSGP, sampledParticles_Sub,
     surface_beta_draws, w_x_draws, w_y_draws, sigma_vel_draws,
     file = '/home/kdp2abu/NewInterpolatorVI_1.RData')


hist(sigma_vel_draws)

truth = plotFullTrajectoriesHSGP(sampledHSGP, grid_res = 100, color_limits = c(0,4), border = c(-12,-12,12,12), grid_border = c(-10,-10,10,10))

truth[[1]]$Coefficient
truth[[2]]$Coefficient

FullTrajPost = getGridStanHSGP_Beta(surface_beta_draws = surface_beta_draws, grid_res = 100, M = 5, n_iters = 1000, border = c(-12,-12,12,12), grid_border = c(-10,-10,10,10))
FullTraj1_PostMean = fast_mean_3d(FullTrajPost$traj_1)
FullTraj1_Lower = fast_quantile_3d(FullTrajPost$traj_1, 0.025)
FullTraj1_Upper = fast_quantile_3d(FullTrajPost$traj_1, 0.975)

FullTraj2_PostMean = fast_mean_3d(FullTrajPost$traj_2)
FullTraj2_Lower = fast_quantile_3d(FullTrajPost$traj_2, 0.025)
FullTraj2_Upper = fast_quantile_3d(FullTrajPost$traj_2, 0.975)

FullTraj1_MSE = median((truth[[1]]$CoefValues - FullTraj1_PostMean)^2)
FullTraj2_MSE = median((truth[[2]]$CoefValues - FullTraj2_PostMean)^2)
FullTraj1_CR = mean(truth[[1]]$CoefValues > FullTraj1_Lower & truth[[1]]$CoefValues < FullTraj1_Upper)
FullTraj2_CR = mean(truth[[2]]$CoefValues > FullTraj2_Lower & truth[[2]]$CoefValues < FullTraj2_Upper)

FullTraj1_Avg_CI_Width = median(FullTraj1_Upper - FullTraj1_Lower)
FullTraj2_Avg_CI_Width = median(FullTraj2_Upper - FullTraj2_Lower)


x_seq = FullTrajPost$x
y_seq = FullTrajPost$y

c_min = 0
c_max = 4

plot_ly() %>%

  # # Layer 1: The Lower Bound (Floor)
  add_surface(x = x_seq, y = y_seq, z = FullTraj1_Lower,
              opacity = 0.3,           # <--- Make it ghostly
              colorscale = "Viridis",
              cmin = c_min, cmax = c_max,
              showscale = FALSE) %>%   # Hide legend for this layer
  #
  # # Layer 2: The Upper Bound (Ceiling)
  add_surface(x = x_seq, y = y_seq, z = FullTraj1_Upper,
              opacity = 0.3,           # <--- Make it ghostly
              colorscale = "Viridis",
              cmin = c_min, cmax = c_max,
              showscale = FALSE) %>%

  # # Layer 3: The Posterior Mean (The Core)
  # add_surface(x = x_seq, y = y_seq, z = FullTraj1_PostMean,
  #             opacity = 0.8,           # <--- Solid
  #             colorscale = "Viridis",
  #             cmin = c_min, cmax = c_max,
  #             colorbar = list(title = "Trajectory")) %>%

  # Layer 4: The Truth
  add_surface(x = x_seq, y = y_seq, z = truth[[1]]$CoefValues,
              opacity = 1.0,           # <--- Solid
              colorscale = "Viridis",
              cmin = c_min, cmax = c_max,
              colorbar = list(title = "Trajectory")) %>%

  layout(title = "Truth with 95% Credible Envelope",
         scene = list(
           zaxis = list(title = "Trajectory Value")
         ))


plot_ly() %>%

  # Layer 1: The Lower Bound (Floor)
  add_surface(x = x_seq, y = y_seq, z = FullTraj2_Lower,
              opacity = 0.3,           # <--- Make it ghostly
              colorscale = "Viridis",
              cmin = c_min, cmax = c_max,
              showscale = FALSE) %>%   # Hide legend for this layer

  # Layer 2: The Upper Bound (Ceiling)
  add_surface(x = x_seq, y = y_seq, z = FullTraj2_Upper,
              opacity = 0.3,           # <--- Make it ghostly
              colorscale = "Viridis",
              cmin = c_min, cmax = c_max,
              showscale = FALSE) %>%

  # # Layer 3: The Posterior Mean (The Core)
  # add_surface(x = x_seq, y = y_seq, z = FullTraj2_PostMean,
  #             opacity = 0.8,           # <--- Solid
  #             colorscale = "Viridis",
  #             cmin = c_min, cmax = c_max,
  #             colorbar = list(title = "Trajectory")) %>%

  # Layer 4: The Truth
  add_surface(x = x_seq, y = y_seq, z = truth[[2]]$CoefValues,
              opacity = 1.0,           # <--- Solid
              colorscale = "Viridis",
              cmin = c_min, cmax = c_max,
              colorbar = list(title = "Trajectory")) %>%

  layout(title = "Truth with 95% Credible Envelope",
         scene = list(
           zaxis = list(title = "Trajectory Value")
         ))

exporter_model = cmdstan_model('NewInterpolationModel_Exporter.stan', cpp_options = list(stan_threads = T))

stan_export <- exporter_model$sample(
  data = full_data,    # Your exact same data list
  fixed_param = TRUE,  # Turns off MCMC/Variational Inference
  iter_sampling = 1,   # Only do this once!
  iter_warmup = 0,
  threads_per_chain = max_cores/4,
  chains = 1
)

Z_elements = stan_export$draws(variables = 'Z_elements', format = 'draws_matrix')[1,]
Z_elements_v = as.vector(Z_elements)


Z_per_drifter = unlist(lapply(Z_flat_full_list, length))
Z_flat_full_list_STAN = list()
z_start = 1
d=1
for(d in 1:D){

  cur_z = Z_per_drifter[d]
  z_end = z_start + cur_z - 1

  Z_flat_full_list_STAN[[d]] = Z_elements_v[z_start:z_end]

  z_start = z_end + 1

}

plotPostDrifter = function(d){

  curDrifter = sampledParticles_Sub[sampledParticles_Sub$Particle == drifter_names[d],]

  post_draws_coef_x = matrix(rep(as.numeric(c_part_x_list[[d]]), 1000), nrow = 1000, ncol = K[d], byrow = T) + w_x_draws[,1:M_drifter + (d-1)*M_drifter] %*% t(matrix(Z_flat_full_list_STAN[[d]], ncol = M_drifter, byrow = F))
  post_draws_coef_y = matrix(rep(as.numeric(c_part_y_list[[d]]), 1000), nrow = 1000, ncol = K[d], byrow = T) + w_y_draws[,1:M_drifter + (d-1)*M_drifter] %*% t(matrix(Z_flat_full_list_STAN[[d]], ncol = M_drifter, byrow = F))

  post_draws_pos_x = post_draws_coef_x %*% t(Phi_List[[d]])
  post_draws_pos_y = post_draws_coef_y %*% t(Phi_List[[d]])

  post_mean_path_x = as.numeric(colMeans(post_draws_pos_x))
  post_mean_path_y = as.numeric(colMeans(post_draws_pos_y))

  post_lower_path_x = as.numeric(apply(post_draws_pos_x, MARGIN = 2, FUN = quantile, p = 0.025))
  post_lower_path_y = as.numeric(apply(post_draws_pos_y, MARGIN = 2, FUN = quantile, p = 0.025))

  post_upper_path_x = as.numeric(apply(post_draws_pos_x, MARGIN = 2, FUN = quantile, p = 0.975))
  post_upper_path_y = as.numeric(apply(post_draws_pos_y, MARGIN = 2, FUN = quantile, p = 0.975))

  ggplot() + geom_line(aes(x = t_grid_list[[d]], y = post_mean_path_x), color = 'blue') +
             geom_point(data = curDrifter, aes(x = t, y = X1), color = 'black', size = 0.1)

  ggplot() + geom_line(aes(x = t_grid_list[[d]], y = post_mean_path_y), color = 'blue') +
             geom_point(data = curDrifter, aes(x = t, y = X2), color = 'black', size = 0.1)

  ggplot() + geom_path(aes(x = post_mean_path_x, y = post_mean_path_y), color = 'blue') +
    geom_point(data = curDrifter, aes(x = X1, y = X2), color = 'black', size = 1)

}

plotPriorDrifter(4)
plotPostDrifter(4)
d=1

save(list = ls(all.names = TRUE), file = "April30NewInterpolationTesting.RData")
d=1
checkDrifterVel = function(d){

  curDrifter = sampledParticles_Sub[sampledParticles_Sub$Particle == drifter_names[d],]

  post_draws_coef_x = matrix(rep(as.numeric(c_part_x_list[[d]]), 1000), nrow = 1000, ncol = K[d], byrow = T) + w_x_draws[,1:M_drifter + (d-1)*M_drifter] %*% t(matrix(Z_flat_full_list_STAN[[d]], ncol = M_drifter, byrow = F))
  post_draws_coef_y = matrix(rep(as.numeric(c_part_y_list[[d]]), 1000), nrow = 1000, ncol = K[d], byrow = T) + w_y_draws[,1:M_drifter + (d-1)*M_drifter] %*% t(matrix(Z_flat_full_list_STAN[[d]], ncol = M_drifter, byrow = F))

  post_draws_pos_x = post_draws_coef_x %*% t(Phi_List[[d]])
  post_draws_pos_y = post_draws_coef_y %*% t(Phi_List[[d]])

  post_mean_pos_x = as.numeric(colMeans(post_draws_pos_x))
  post_mean_pos_y = as.numeric(colMeans(post_draws_pos_y))

  post_draws_vel_x = post_draws_coef_x %*% t(Phi_d_List[[d]])
  post_draws_vel_y = post_draws_coef_y %*% t(Phi_d_List[[d]])

  post_mean_vel_x = as.numeric(colMeans(post_draws_vel_x))
  post_mean_vel_y = as.numeric(colMeans(post_draws_vel_y))

  post_mean_pos_mat = matrix(c(post_mean_pos_x, post_mean_pos_y), nrow = 1000, byrow = F)

  All_Traj_Vals = getPointsStanHSGP_Beta(surface_beta_draws, post_mean_pos_mat, M = 5, n_iters = 1000, border = c(-12,-12,12,12))

  Traj_1_mean = rowMeans(All_Traj_Vals$traj_1)
  Traj_2_mean = rowMeans(All_Traj_Vals$traj_2)

  quad_grid = seq(min(curDrifter$t), max(curDrifter$t), length.out = N_quad)

  eval_points = matrix(c(quad_grid, post_mean_pos_x, post_mean_pos_y), nrow = N_quad, byrow = F)

  ## Geo NN Eval

  baseVFs = t(apply(eval_points, MARGIN = 1, FUN = function(row){baseVectorFields(row[1], row[-1])}))

  ## Est VF

  Est_VF_long = baseVFs[,1] * Traj_1_mean + baseVFs[,3] * Traj_2_mean
  Est_VF_lat = baseVFs[,2] * Traj_1_mean + baseVFs[,4] * Traj_2_mean

  cur_post_vel_MSE = mean(c((post_mean_vel_x - Est_VF_long)^2, (post_mean_vel_y - Est_VF_lat)^2))


  #### Base Path

  prior_pos_x = as.numeric(Phi_List[[d]] %*% c_part_x_list[[d]])
  prior_pos_y = as.numeric(Phi_List[[d]] %*% c_part_y_list[[d]])

  prior_vel_x = as.numeric(Phi_d_List[[d]] %*% c_part_x_list[[d]])
  prior_vel_y = as.numeric(Phi_d_List[[d]] %*% c_part_y_list[[d]])

  prior_eval_points = matrix(c(quad_grid, prior_pos_x, prior_pos_y), nrow = N_quad, byrow = F)

  prior_baseVFs = t(apply(prior_eval_points, MARGIN = 1, FUN = function(row){baseVectorFields(row[1], row[-1])}))

  ## Est VF

  prior_Est_VF_long =  prior_baseVFs[,1] + prior_baseVFs[,3]
  prior_Est_VF_lat = prior_baseVFs[,2] + prior_baseVFs[,4]

  cur_prior_vel_MSE = mean(c((prior_vel_x - prior_Est_VF_long)^2, (prior_vel_y - prior_Est_VF_lat)^2))

  c(cur_post_vel_MSE, cur_prior_vel_MSE)

}

checkDrifterVel(1)

checkPriorPost = matrix(ncol = 2, nrow = 16)

for(d in 1:D){

  checkPriorPost[d,] = checkDrifterVel(d)

}

priorPerDraw = rowMeans(checkPriorPost[,1:100 * 2 - 1])

median((checkPriorPost[,2] - checkPriorPost[,1])/checkPriorPost[,2] * 100)


plot(checkPriorPost)
abline(a = 0, b = 1)

colMeans(checkPriorPost)


prior_x_pos_STAN = stan_export$draws('p_x_base', format = 'draws_matrix')
prior_y_pos_STAN = stan_export$draws('p_y_base', format = 'draws_matrix')

colnames(prior_x_pos_STAN)

prior_x_pos_STAN_m = matrix(prior_x_pos_STAN, nrow = 100, byrow = F)
prior_y_pos_STAN_m = matrix(prior_y_pos_STAN, nrow = 100, byrow = F)

prior_x_pos_m = matrix(ncol = 5000, nrow = 100)
prior_y_pos_m = matrix(ncol = 5000, nrow = 100)

for(d in 1:D){
  prior_x_pos_m[d,] = as.numeric(Phi_List[[d]] %*% c_part_x_list[[d]])
  prior_y_pos_m[d,] = as.numeric(Phi_List[[d]] %*% c_part_y_list[[d]])
}

plot(prior_x_pos_STAN_m)
as.vector(prior_x_pos_STAN)


mean(abs(prior_x_pos_m - prior_x_pos_STAN_m))
prior_x_pos_STAN_m[1,]


prior_x_vel_STAN = stan_export$draws('v_x_base', format = 'draws_matrix')
prior_y_vel_STAN = stan_export$draws('v_y_base', format = 'draws_matrix')

colnames(prior_x_pos_STAN)

prior_x_vel_STAN_m = matrix(prior_x_vel_STAN, nrow = 100, byrow = F)
prior_y_vel_STAN_m = matrix(prior_y_vel_STAN, nrow = 100, byrow = F)

prior_x_vel_m = matrix(ncol = 5000, nrow = 100)
prior_y_vel_m = matrix(ncol = 5000, nrow = 100)

for(d in 1:D){
  prior_x_vel_m[d,] = as.numeric(Phi_d_List[[d]] %*% c_part_x_list[[d]])
  prior_y_vel_m[d,] = as.numeric(Phi_d_List[[d]] %*% c_part_y_list[[d]])
}


mean(abs(prior_x_vel_m - prior_x_vel_STAN_m))


Phi_Z_flat_STAN = stan_export$draws('Phi_Z_flat', format = 'draws_matrix')

colnames(Phi_Z_flat_STAN)[5000:6000]

Phi_Z_flat_STAN_m = matrix(Phi_Z_flat_STAN, nrow = 5000, byrow = F)

Phi_Z_list_STAN = list()
d=1
for(d in 1:D){

  Phi_Z_list_STAN[[d]] = Phi_Z_flat_STAN_m[,1:20 + (d-1)*20]


}

post_mean_w_x = colMeans(w_x_draws)
post_mean_w_y = colMeans(w_y_draws)

prior_x_pos_m[1,] + Phi_Z_list_STAN[[1]] %*% post_mean_w_x[1:20]

post_x_pos_m = matrix(nrow = 100, ncol = 5000)
post_y_pos_m = matrix(nrow = 100, ncol = 5000)

for(d in 1:D){

  curDrifter = sampledParticles_Sub[sampledParticles_Sub$Particle == drifter_names[d],]

  post_draws_coef_x = matrix(rep(as.numeric(c_part_x_list[[d]]), 1000), nrow = 1000, ncol = K[d], byrow = T) + w_x_draws[,1:M_drifter + (d-1)*M_drifter] %*% t(matrix(Z_flat_full_list_STAN[[d]], ncol = M_drifter, byrow = F))
  post_draws_coef_y = matrix(rep(as.numeric(c_part_y_list[[d]]), 1000), nrow = 1000, ncol = K[d], byrow = T) + w_y_draws[,1:M_drifter + (d-1)*M_drifter] %*% t(matrix(Z_flat_full_list_STAN[[d]], ncol = M_drifter, byrow = F))

  post_draws_pos_x = post_draws_coef_x %*% t(Phi_List[[d]])
  post_draws_pos_y = post_draws_coef_y %*% t(Phi_List[[d]])

  post_x_pos_m[d,] = as.numeric(colMeans(post_draws_pos_x))
  post_y_pos_m[d,] = as.numeric(colMeans(post_draws_pos_y))
}

post_x_pos_STAN_m = matrix(nrow = 100, ncol = 5000)
post_y_pos_STAN_m = matrix(nrow = 100, ncol = 5000)

for(d in 1:D){

  post_x_pos_STAN_m[d,] = prior_x_pos_STAN_m[d,] + Phi_Z_list_STAN[[d]] %*% post_mean_w_x[1:20 + (d-1)*20]
  post_y_pos_STAN_m[d,] = prior_y_pos_STAN_m[d,] + Phi_Z_list_STAN[[d]] %*% post_mean_w_y[1:20 + (d-1)*20]

}


mean(abs(post_x_pos_m - post_x_pos_STAN_m))
mean(abs(post_y_pos_m - post_y_pos_STAN_m))



Phi_d_Z_flat_STAN = stan_export$draws('Phi_d_Z_flat', format = 'draws_matrix')

Phi_d_Z_flat_STAN_m = matrix(Phi_d_Z_flat_STAN, nrow = 5000, byrow = F)

Phi_d_Z_list_STAN = list()
d=1
for(d in 1:D){

  Phi_d_Z_list_STAN[[d]] = Phi_d_Z_flat_STAN_m[,1:20 + (d-1)*20]


}


post_x_vel_m = matrix(nrow = 100, ncol = 5000)
post_y_vel_m = matrix(nrow = 100, ncol = 5000)

for(d in 1:D){

  curDrifter = sampledParticles_Sub[sampledParticles_Sub$Particle == drifter_names[d],]

  post_draws_coef_x = matrix(rep(as.numeric(c_part_x_list[[d]]), 1000), nrow = 1000, ncol = K[d], byrow = T) + w_x_draws[,1:M_drifter + (d-1)*M_drifter] %*% t(matrix(Z_flat_full_list_STAN[[d]], ncol = M_drifter, byrow = F))
  post_draws_coef_y = matrix(rep(as.numeric(c_part_y_list[[d]]), 1000), nrow = 1000, ncol = K[d], byrow = T) + w_y_draws[,1:M_drifter + (d-1)*M_drifter] %*% t(matrix(Z_flat_full_list_STAN[[d]], ncol = M_drifter, byrow = F))

  post_draws_vel_x = post_draws_coef_x %*% t(Phi_d_List[[d]])
  post_draws_vel_y = post_draws_coef_y %*% t(Phi_d_List[[d]])

  post_x_vel_m[d,] = as.numeric(colMeans(post_draws_vel_x))
  post_y_vel_m[d,] = as.numeric(colMeans(post_draws_vel_y))
}

post_x_vel_STAN_m = matrix(nrow = 100, ncol = 5000)
post_y_vel_STAN_m = matrix(nrow = 100, ncol = 5000)

for(d in 1:D){

  post_x_vel_STAN_m[d,] = prior_x_vel_STAN_m[d,] + Phi_d_Z_list_STAN[[d]] %*% post_mean_w_x[1:20 + (d-1)*20]
  post_y_vel_STAN_m[d,] = prior_y_vel_STAN_m[d,] + Phi_d_Z_list_STAN[[d]] %*% post_mean_w_y[1:20 + (d-1)*20]

}


mean(abs(post_x_vel_m - post_x_vel_STAN_m))
mean(abs(post_y_vel_m - post_y_vel_STAN_m))

### Testing Surface Eval and VFs

otherQuantities_model = cmdstan_model('/home/kdp2abu/NewInterpolationModel_OtherQuantities.stan', cpp_options = list(stan_threads = T))

d = 30

plot(post_x_pos_m[d,], post_y_pos_m[d,])

matrix(surface_beta_draws[500,1:121 * 2 - 1], nrow = 11, byrow = F)

otherQuantities_data = list(

  N_quad = 5000,
  N_models = 2,

  t = t_grid_list[[d]],
  x = post_x_pos_m[d,],
  y = post_y_pos_m[d,],

  M_Surface = 10,
  surface_betas = matrix(surface_beta_draws[500,1:121 * 2 - 1], nrow = 11, byrow = F),
  omega_surface = (0:10)*pi,

  Lx = 25,
  Ly = 25,
  border = c(-12.5,-12.5,12.5,12.5)


)

OQ_export <- otherQuantities_model$sample(
  data = otherQuantities_data,    # Your exact same data list
  fixed_param = TRUE,  # Turns off MCMC/Variational Inference
  iter_sampling = 1,   # Only do this once!
  iter_warmup = 0,
  threads_per_chain = max_cores/4,
  chains = 1
)

STAN_Coef1 = OQ_export$draws('Coef1', format = 'draws_matrix')

R_Coef1 = getPointsStanHSGP_Beta(surface_beta_draws = surface_beta_draws[500,], points = matrix(c(post_x_pos_m[d,], post_y_pos_m[d,]), byrow = F, ncol = 2),
                       M = 10, n_iters = 1, border = c(-12.5, -12.5, 12.5, 12.5))

mean(abs(R_Coef1$traj_1 - c(OQ_export$draws('Coef1', format = 'draws_matrix'))))

STAN_VFs = OQ_export$draws('VFs', format = 'draws_matrix')

R_VFs = apply(matrix(c(t_grid_list[[d]], post_x_pos_m[d,], post_y_pos_m[d,]), ncol = 3, byrow = F), MARGIN = 1, FUN = function(row){baseVectorFields(row[1], row[-1])})

STAN_VFs_M1_x = STAN_VFs[1:5000 * 2 - 1]
R_VFs_M1_x = R_VFs[1,]

(STAN_VFs_M1_x - R_VFs_M1_x)
