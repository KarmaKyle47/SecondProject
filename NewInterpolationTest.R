source('SecondProjectRequiredFunctions.R')
options(mc.cores = parallel::detectCores())
max_cores = parallel::detectCores()
library(splines2)


baseVectorFields = function(t, curPos){

  window = cos(pi*curPos[1]/20) * cos(pi*curPos[2]/20)

  f1 = c(curPos[2],-1*curPos[1]) * window
  f2 = c(curPos[1],curPos[2])*sin(t/8) * window
  matrix(c(f1,f2), nrow = 2, byrow = F)

}

baseVectorFields = function(t, curPos){

  f1 = c(1,0)
  f2 = c(0,1)
  matrix(c(f1,f2), nrow = 2, byrow = F)

}

baseVectorFields = function(t, curPos) {
  x = curPos[1]
  y = curPos[2]

  # 1. Original window function
  window = cos(pi * x / 20) * cos(pi * y / 20)

  # 2. Original base fields (multiplied by window)
  f1 = c(y, -1 * x) * window
  f2 = c(x, y) * sin(t / 8) * window

  # 3. Boundary Repulsion Setup
  boundary_limit = 10  # Matches the zeroes of your cosine function
  k = 0.1             # Repulsion strength (tune this!)
  eps = 0.01           # Small buffer to prevent division by zero

  # Inverse square repulsion:
  # As 'x' approaches 'boundary_limit', the second term explodes negatively (pushes left).
  # As 'x' approaches '-boundary_limit', the first term explodes positively (pushes right).
  repel_x = k * (1 / (x + boundary_limit + eps)^2 - 1 / (boundary_limit - x + eps)^2)
  repel_y = k * (1 / (y + boundary_limit + eps)^2 - 1 / (boundary_limit - y + eps)^2)

  repulsion = c(repel_x, repel_y)

  # 4. Add the repulsion to your fields (outside the window multiplier)
  f1 = f1 + repulsion
  f2 = f2 + repulsion

  matrix(c(f1, f2), nrow = 2, byrow = FALSE)
}

true_log_k = 0.35
true_log_l = 1
M = 4

sampledHSGP = sampleFullTrajectoriesHSGP(2, M, true_log_k, true_log_l)

rand_res = 100

sampledParticles_V2 = samplePhySpaceParticles(10, startTime = 0, n_obs = 100*rand_res, border = c(-10,-10,10,10), borderBuffer = 0.2, baseVectorFields, sampledHSGP,
                                              M = M, t_step_mean = 0.01, vel_sigma = 0, pos_sigma = 0)

sampledParticles_Sub = sampledParticles_V2[1:1000 * rand_res - (rand_res-1),]

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
K_drifter = 100

p_x_base = list()
p_y_base = list()
v_x_base = list()
v_y_base = list()

Phi_List = list()
Phi_d_List = list()
Phi_data_List = list()

t_grid_list = list()
Lts = rep(0,D)
drifter_boundaries = matrix(nrow = D, ncol = 2)
d=1
for(d in 1:D){

  curDrifter = sampledParticles_Sub[sampledParticles_Sub$Particle == drifter_names[d],]
  cur_boundary = c(min(curDrifter$t),max(curDrifter$t))
  cur_quad_grid = seq(cur_boundary[1], cur_boundary[2], length.out = N_quad)

  cur_x_spline = splinefun(x = curDrifter$t, y = curDrifter$X1, method = 'fmm')
  p_x_base[[d]] = cur_x_spline(cur_quad_grid)
  v_x_base[[d]] = cur_x_spline(cur_quad_grid, deriv = 1)

  cur_y_spline = splinefun(x = curDrifter$t, y = curDrifter$X2, method = 'fmm')
  p_y_base[[d]] = cur_y_spline(cur_quad_grid)
  v_y_base[[d]] = cur_y_spline(cur_quad_grid, deriv = 1)

  Phi_List[[d]] = bSpline(x = cur_quad_grid, df = K_drifter, intercept = F)
  Phi_d_List[[d]] = bSpline(x = cur_quad_grid, df = K_drifter, intercept = F, derivs = 1)

  Phi_data_List[[d]] = bSpline(x = curDrifter$t, df = K_drifter, intercept = F)

  t_grid_list[[d]] = cur_quad_grid
  Lts[d] = diff(cur_boundary)

  drifter_boundaries[d,] = cur_boundary

  svMisc::progress(d, D)

}
d=1
plotPriorDrifter = function(d){

  curDrifter = sampledParticles_Sub[sampledParticles_Sub$Particle == drifter_names[d],]

  prior_pos_x = p_x_base[[d]]
  prior_pos_y = p_y_base[[d]]

  ggplot() + geom_line(aes(x = t_grid_list[[d]], y = prior_pos_x), color = 'blue') +
             geom_point(data = curDrifter, aes(x = t, y = X1), color = 'black', size = 1)

  # ggplot() + geom_line(aes(x = t_grid_list[[d]], y = prior_pos_y), color = 'blue') +
  #            geom_point(data = curDrifter, aes(x = t, y = X2), color = 'black', size = 0.1)
  #
  ggplot() + geom_path(aes(x = prior_pos_x, y = prior_pos_y), color = 'blue') +
    geom_point(data = curDrifter, aes(x = X1, y = X2), color = 'black', size = 1)
}

plotPriorDrifter(2)

unlist(lapply(Z_flat_full_list, length)) == K*M_drifter
length(Z_flat_full_list[[1]])
do.call(, Z_flat_full_list)

unlist(c_part_x_list)

full_data = list(
  D = D,
  N_quad = N_quad,

  K = rep(K_drifter, D),
  total_K = K_drifter*D,

  N_drifter = drifter_num_points,
  total_data = sum(drifter_num_points),

  total_Z = sum(K_drifter*drifter_num_points),

  p_x_base = p_x_base,
  p_y_base = p_y_base,
  v_x_base = v_x_base,
  v_y_base = v_y_base,

  Phi_flat = do.call(cbind, Phi_List),
  Phi_d_flat = do.call(cbind, Phi_d_List),

  Phi_data_elements = unlist(Phi_data_List),

  drifter_boundaries = drifter_boundaries,

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

spline_interpolator_model = cmdstan_model('BaseSplineInterpolationModel.stan', cpp_options = list(stan_threads = T), force_recompile = T)

spline_fit = spline_interpolator_model$variational(
  data = full_data,
  threads = max_cores/4,
  iter = 10000,
  draws = 1000, show_messages = T, init = 0#, algorithm = "fullrank"
)

PF_sum = spline_fit$summary()

PF_sum[PF_sum$variable == 'fixed_ls',]

surface_beta_draws = spline_fit$draws('surface_betas')
w_x_draws = spline_fit$draws('w_x')
w_y_draws = spline_fit$draws('w_y')
sigma_mag_draws = spline_fit$draws('sigma_mag')
sigma_angle_draws = spline_fit$draws('sigma_angle')
sigma_pos_draws = spline_fit$draws('sigma_pos')
rho_mag_draws = spline_fit$draws('rho_mag')
rho_angle_draws = spline_fit$draws('rho_angle')

hist(rho_angle_draws)

c(full_fit$draws('post_vel_mse'))
c(full_fit$draws('prior_vel_mse'))

save(sampledHSGP, sampledParticles_Sub,
     surface_beta_draws, w_x_draws, w_y_draws, sigma_vel_draws,
     file = 'UpdatedModelwL5_28.RData')

load('UpdatedModelwL5_28.RData')


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

  post_draws_pos_x = w_x_draws[,1:K_drifter + (d-1)*K_drifter] %*% t(Phi_List[[d]])
  post_draws_pos_y = w_y_draws[,1:K_drifter + (d-1)*K_drifter] %*% t(Phi_List[[d]])

  post_mean_path_x = p_x_base[[d]] + as.numeric(colMeans(post_draws_pos_x))
  post_mean_path_y = p_y_base[[d]] + as.numeric(colMeans(post_draws_pos_y))

  post_lower_path_x = p_x_base[[d]] + as.numeric(apply(post_draws_pos_x, MARGIN = 2, FUN = quantile, p = 0.025))
  post_lower_path_y = p_y_base[[d]] + as.numeric(apply(post_draws_pos_y, MARGIN = 2, FUN = quantile, p = 0.025))

  post_upper_path_x = p_x_base[[d]] + as.numeric(apply(post_draws_pos_x, MARGIN = 2, FUN = quantile, p = 0.975))
  post_upper_path_y = p_y_base[[d]] + as.numeric(apply(post_draws_pos_y, MARGIN = 2, FUN = quantile, p = 0.975))

  ggplot() + geom_line(aes(x = t_grid_list[[d]], y = post_mean_path_x), color = 'blue') +
    geom_line(aes(x = t_grid_list[[d]], y = post_lower_path_x), color = 'red', alpha = 0.2) +
    geom_line(aes(x = t_grid_list[[d]], y = post_upper_path_x), color = 'red', alpha = 0.2) +
             geom_point(data = curDrifter, aes(x = t, y = X1), color = 'black', size = 1)


  # ggplot() + geom_line(aes(x = t_grid_list[[d]], y = post_mean_path_y), color = 'blue') +
  #            geom_point(data = curDrifter, aes(x = t, y = X2), color = 'black', size = 0.1)
  #
  ggplot() + geom_path(aes(x = as.numeric(post_mean_path_x), y = as.numeric(post_mean_path_y)), color = 'blue') +
    geom_point(data = curDrifter, aes(x = X1, y = X2), color = 'black', size = 1)

}

plotPriorDrifter(1)
plotPostDrifter(1)
d=1

save(list = ls(all.names = TRUE), file = "June4SplineInterpolationTesting_FirstMagAngle_wAR1.RData")
d=1
load("May31SplineInterpolationTesting.RData")

checkDrifterVel = function(d){

  curDrifter = sampledParticles_Sub[sampledParticles_Sub$Particle == drifter_names[d],]

  post_draws_pos_x = w_x_draws[,1:K_drifter + (d-1)*K_drifter] %*% t(Phi_List[[d]])
  post_draws_pos_y = w_y_draws[,1:K_drifter + (d-1)*K_drifter] %*% t(Phi_List[[d]])

  post_draws_vel_x = w_x_draws[,1:K_drifter + (d-1)*K_drifter] %*% t(Phi_d_List[[d]])
  post_draws_vel_y = w_y_draws[,1:K_drifter + (d-1)*K_drifter] %*% t(Phi_d_List[[d]])

  post_mean_pos_x = p_x_base[[d]] + as.numeric(colMeans(post_draws_pos_x))
  post_mean_pos_y = p_y_base[[d]] + as.numeric(colMeans(post_draws_pos_y))

  post_mean_vel_x = v_x_base[[d]] + as.numeric(colMeans(post_draws_vel_x))
  post_mean_vel_y = v_y_base[[d]] + as.numeric(colMeans(post_draws_vel_y))

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

  prior_pos_x = p_x_base[[d]]
  prior_pos_y = p_y_base[[d]]

  prior_vel_x = v_x_base[[d]]
  prior_vel_y = v_y_base[[d]]

  prior_eval_points = matrix(c(quad_grid, prior_pos_x, prior_pos_y), nrow = N_quad, byrow = F)

  prior_baseVFs = t(apply(prior_eval_points, MARGIN = 1, FUN = function(row){baseVectorFields(row[1], row[-1])}))

  ## Est VF

  prior_Est_VF_long =  prior_baseVFs[,1] + prior_baseVFs[,3]
  prior_Est_VF_lat = prior_baseVFs[,2] + prior_baseVFs[,4]

  cur_prior_vel_MSE = mean(c((prior_vel_x - prior_Est_VF_long)^2, (prior_vel_y - prior_Est_VF_lat)^2))

  c(cur_post_vel_MSE, cur_prior_vel_MSE)

}

checkDrifterVel(1)

checkPriorPost = matrix(ncol = 2, nrow = D)

for(d in 1:D){

  checkPriorPost[d,] = checkDrifterVel(d)

}

priorPerDraw = rowMeans(checkPriorPost[,1:100 * 2 - 1])

median((checkPriorPost[,2] - checkPriorPost[,1])/checkPriorPost[,2] * 100)


plot(checkPriorPost)
abline(a = 0, b = 1)

mean(checkPriorPost[,2] > checkPriorPost[,1])

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




#Super Simple 1 to 1 model

log_coef1 = rnorm(1,0,0.35)
log_coef2 = rnorm(1,0,0.35)

# 2. Create X and Y vectors
# seq() generates a sequence of numbers from -10 to 10
x <- seq(-10, 10, length.out = 50)
y <- seq(-10, 10, length.out = 50)

# 3. Create a Z matrix filled with the constant value
# We create a 50x50 matrix where every single cell is 5
z <- matrix(exp(log_coef1), nrow = length(y), ncol = length(x))

# 4. Create the surface plot
fig <- plot_ly(
  x = ~x,
  y = ~y,
  z = ~z,
  type = "surface",
  colorscale = "Viridis"
)

plotHSGP

# 5. Update layout for better viewing
fig <- fig %>% layout(
  title = paste("Constant Plane Surface (Z =", constant_z, ")"),
  scene = list(
    xaxis = list(title = "X Axis"),
    yaxis = list(title = "Y Axis"),
    zaxis = list(
      title = "Z Axis",
      # Set the z-axis range so the plane doesn't look like a flat box
      range = c(0, 10)
    )
  )
)

# 6. Show the plot
fig


N_drifter = 100
K_drifter = 5

start_pos_x = runif(N_drifter, min = -10, max = 10)
start_pos_y = runif(N_drifter, min = -10, max = 10)

sigma_vel = 0.5

w_x_matrix = matrix(rnorm(n = N_drifter*K_drifter,0,sigma_vel), nrow = K_drifter)
w_y_matrix = matrix(rnorm(n = N_drifter*K_drifter,0,sigma_vel), nrow = K_drifter)

t_start = 0
t_end = 10
t_res = 100
t_grid_data = seq(t_start,t_end,length.out = t_res)

p_x_base = exp(log_coef1) * t_grid_data
p_y_base = exp(log_coef2) * t_grid_data

p_x_base_mat = matrix(nrow = t_res, ncol = N_drifter)
p_y_base_mat = matrix(nrow = t_res, ncol = N_drifter)

for(d in 1:N_drifter){

  p_x_base_mat[,d] = p_x_base + start_pos_x[d]
  p_y_base_mat[,d] = p_y_base + start_pos_y[d]


}

Phi_data = bSpline(x = t_grid_data, df = K_drifter, intercept = F)
#
p_x_mat = p_x_base_mat + Phi_data %*% w_x_matrix
p_y_mat = p_y_base_mat + Phi_data %*% w_y_matrix

sigma_pos = 0.1

p_x_mat_data = p_x_mat + matrix(rnorm(N_drifter*t_res, 0, sigma_pos), nrow = t_res)
p_y_mat_data = p_y_mat + matrix(rnorm(N_drifter*t_res, 0, sigma_pos), nrow = t_res)


Data = data.frame(t = rep(t_grid_data, N_drifter), x = c(p_x_mat_data), y = c(p_y_mat_data),
                  Particle = rep(str_c("Particle",1:N_drifter), each = t_res))

ggplot(Data, aes(x = x, y = y, color = Particle)) + geom_point() + theme(legend.position = "none")


data_stan = list(
  D = N_drifter,
  K = K_drifter,
  N_drifter = t_res,
  Phi = Phi_data,

  t_grid_data = t_grid_data,
  x_pos = p_x_mat_data,
  y_pos = p_y_mat_data
)

simple_interpolator_model = cmdstan_model('SuperSimpleBaseSplineInterpolator.stan')
line_simple_model = cmdstan_model('StraightLineSimple.stan')


simple_fit = simple_interpolator_model$variational(
  data = data_stan,
  iter = 100000,
  draws = 1000, show_messages = T, init = 0#, algorithm = "fullrank"
)

straight_fit_VI = line_simple_model$variational(
  data = data_stan,
  iter = 100000,
  draws = 1000, show_messages = T, init = 0#, algorithm = "fullrank"
)

straight_fit_HMC = straight_fit

straight_fit = line_simple_model$sample(
  data = data_stan,
  chains = 4,
  parallel_chains = 4,
  refresh = 10
)

simple_fit = simple_interpolator_model$sample(
  data = data_stan,
  chains = 4,
  parallel_chains = 4,
  refresh = 10
)

simple_fit = straight_fit_VI

SF_sum = simple_fit$summary()

SF_sum$rhat

sigma_pos_draws = simple_fit$draws('sigma_pos')
sigma_pos_mean = mean(sigma_pos_draws)
sigma_pos_lower = as.numeric(quantile(sigma_pos_draws,0.025))
sigma_pos_upper = as.numeric(quantile(sigma_pos_draws,0.975))

sigma_pos_mean - sigma_pos
sigma_pos_lower <= sigma_pos & sigma_pos_upper >= sigma_pos

sigma_vel_draws = simple_fit$draws('sigma_vel')
sigma_vel_mean = mean(sigma_vel_draws)
sigma_vel_lower = as.numeric(quantile(sigma_vel_draws,0.025))
sigma_vel_upper = as.numeric(quantile(sigma_vel_draws,0.975))

sigma_vel_mean - sigma_vel
sigma_vel_lower <= sigma_vel & sigma_vel_upper >= sigma_vel

log_coef1_draws = simple_fit$draws('logcoef1')
log_coef1_mean = mean(log_coef1_draws)
log_coef1_lower = as.numeric(quantile(log_coef1_draws,0.025))
log_coef1_upper = as.numeric(quantile(log_coef1_draws,0.975))

log_coef1_mean - log_coef1
log_coef1_lower <= log_coef1 & log_coef1_upper >= log_coef1

log_coef2_draws = simple_fit$draws('logcoef2')
log_coef2_mean = mean(log_coef2_draws)
log_coef2_lower = as.numeric(quantile(log_coef2_draws,0.025))
log_coef2_upper = as.numeric(quantile(log_coef2_draws,0.975))

log_coef2_mean - log_coef2
log_coef2_lower <= log_coef2 & log_coef2_upper >= log_coef2

start_pos_x_draws = simple_fit$draws('start_pos_x')
start_pos_x_mean = colMeans(start_pos_x_draws)
start_pos_x_lower = as.numeric(apply(start_pos_x_draws, MARGIN = 2, FUN = quantile, p = 0.025))
start_pos_x_upper = as.numeric(apply(start_pos_x_draws, MARGIN = 2, FUN = quantile, p = 0.975))

start_pos_x_mean - start_pos_x
start_pos_x_lower <= start_pos_x & start_pos_x_upper >= start_pos_x
mean(start_pos_x_lower <= start_pos_x & start_pos_x_upper >= start_pos_x)

start_pos_y_draws = simple_fit$draws('start_pos_y')
start_pos_y_mean = colMeans(start_pos_y_draws)
start_pos_y_lower = as.numeric(apply(start_pos_y_draws, MARGIN = 2, FUN = quantile, p = 0.025))
start_pos_y_upper = as.numeric(apply(start_pos_y_draws, MARGIN = 2, FUN = quantile, p = 0.025))

start_pos_y_mean - start_pos_y
start_pos_y_lower <= start_pos_y & start_pos_y_upper >= start_pos_y
mean(start_pos_y_lower <= start_pos_y & start_pos_y_upper >= start_pos_y)


w_x_draws = simple_fit$draws('w_x')
w_x_mean = matrix(colMeans(w_x_draws), nrow = K_drifter, byrow = F)
w_x_lower = matrix(apply(w_x_draws, MARGIN = 2, FUN = quantile, p = 0.025), nrow = K_drifter, byrow = F)
w_x_upper = matrix(apply(w_x_draws, MARGIN = 2, FUN = quantile, p = 0.975), nrow = K_drifter, byrow = F)

w_x_mean - w_x_matrix
w_x_lower <= w_x_matrix & w_x_upper >= w_x_matrix

mean(w_x_lower <= w_x_matrix & w_x_upper >= w_x_matrix)

w_y_draws = simple_fit$draws('w_y')
w_y_mean = matrix(colMeans(w_y_draws), nrow = K_drifter, byrow = F)
w_y_lower = matrix(apply(w_y_draws, MARGIN = 2, FUN = quantile, p = 0.025), nrow = K_drifter, byrow = F)
w_y_upper = matrix(apply(w_y_draws, MARGIN = 2, FUN = quantile, p = 0.975), nrow = K_drifter, byrow = F)

w_y_mean - w_y_matrix
w_y_lower <= w_y_matrix & w_y_upper >= w_y_matrix

mean(w_y_lower <= w_y_matrix & w_y_upper >= w_y_matrix)


p_x_base_post = exp(log_coef1_mean) * t_grid_data
p_y_base_post = exp(log_coef2_mean) * t_grid_data

p_x_base_mat_post = matrix(nrow = t_res, ncol = N_drifter)
p_y_base_mat_post = matrix(nrow = t_res, ncol = N_drifter)

for(d in 1:N_drifter){

  p_x_base_mat_post[,d] = p_x_base_post + start_pos_x_mean[d]
  p_y_base_mat_post[,d] = p_y_base_post + start_pos_y_mean[d]


}

p_x_mat_post = p_x_base_mat_post + Phi_data %*% w_x_mean
p_y_mat_post = p_y_base_mat_post + Phi_data %*% w_y_mean


Data_post = data.frame(t = rep(t_grid_data, N_drifter), x = c(p_x_mat_post), y = c(p_y_mat_post),
                  Particle = rep(str_c("Particle",1:N_drifter), each = t_res))

ggplot(Data_post, aes(x = x, y = y, color = Particle)) + geom_point() + theme(legend.position = "none")



Phi_data %*% w_x_mean
Phi_data %*% w_y_mean


SF_sum = straight_fit$summary()

straight_fit = straight_fit_VI

max(SF_sum$rhat)

a_x_draws = straight_fit$draws('a_x')
b_x_draws = straight_fit$draws('b_x', format = 'draws_matrix')

hist(b_x_draws)

b_x_lower = apply(b_x_draws, MARGIN = 2, FUN = quantile, p = 0.025)
b_x_upper = apply(b_x_draws, MARGIN = 2, FUN = quantile, p = 0.975)

(mean(b_x_draws) - exp(log_coef1))^2


b_x_upper >= exp(log_coef1) & b_x_lower <= exp(log_coef1)

mean(b_x_upper >= exp(log_coef1))
mean(b_x_lower <= exp(log_coef1))
mean(b_x_upper >= exp(log_coef1) & b_x_lower <= exp(log_coef1))

exp(log_coef1)

a_y_draws = straight_fit$draws('a_y')
b_y_draws = straight_fit$draws('b_y')

b_y_draws = straight_fit$draws('b_y', format = 'draws_matrix')
b_y_lower = apply(b_y_draws, MARGIN = 2, FUN = quantile, p = 0.025)
b_y_upper = apply(b_y_draws, MARGIN = 2, FUN = quantile, p = 0.975)

hist(b_x_draws[,7])
hist(straight_fit_VI$draws('b_x')[,7])


mean(b_y_upper >= exp(log_coef2) & b_y_lower <= exp(log_coef2))

exp(log_coef2)

hist(b_y_draws)

log_coef1_draws = straight_fit$draws('logcoef1', format = 'draws_matrix')
log_coef1_lower = as.numeric(quantile(exp(log_coef1_draws), p = 0.025))
log_coef1_upper = as.numeric(quantile(exp(log_coef1_draws), p = 0.975))


(mean(exp(log_coef1_draws)) - exp(log_coef1))^2


hist(exp(log_coef1_draws))

log_coef2_draws = straight_fit$draws('logcoef2')
hist(exp(log_coef2_draws))

straight_fit$draws('sigma_vel')


run_model_VI = function(){

  log_coef1 = rnorm(1,0,0.35)
  log_coef2 = rnorm(1,0,0.35)

  N_drifter = 100
  K_drifter = 5

  start_pos_x = runif(N_drifter, min = -10, max = 10)
  start_pos_y = runif(N_drifter, min = -10, max = 10)

  sigma_vel = 0.5

  w_x_matrix = matrix(rnorm(n = N_drifter*K_drifter,0,sigma_vel), nrow = K_drifter)
  w_y_matrix = matrix(rnorm(n = N_drifter*K_drifter,0,sigma_vel), nrow = K_drifter)

  t_start = 0
  t_end = 10
  t_res = 100
  t_grid_data = seq(t_start,t_end,length.out = t_res)

  p_x_base = exp(log_coef1) * t_grid_data
  p_y_base = exp(log_coef2) * t_grid_data

  p_x_base_mat = matrix(nrow = t_res, ncol = N_drifter)
  p_y_base_mat = matrix(nrow = t_res, ncol = N_drifter)

  for(d in 1:N_drifter){

    p_x_base_mat[,d] = p_x_base + start_pos_x[d]
    p_y_base_mat[,d] = p_y_base + start_pos_y[d]


  }

  Phi_data = bSpline(x = t_grid_data, df = K_drifter, intercept = F)
  #
  p_x_mat = p_x_base_mat + Phi_data %*% w_x_matrix
  p_y_mat = p_y_base_mat + Phi_data %*% w_y_matrix

  sigma_pos = 0.1

  p_x_mat_data = p_x_mat + matrix(rnorm(N_drifter*t_res, 0, sigma_pos), nrow = t_res)
  p_y_mat_data = p_y_mat + matrix(rnorm(N_drifter*t_res, 0, sigma_pos), nrow = t_res)

  data_stan = list(
    D = N_drifter,
    K = K_drifter,
    N_drifter = t_res,
    Phi = Phi_data,

    t_grid_data = t_grid_data,
    x_pos = p_x_mat_data,
    y_pos = p_y_mat_data
  )

  straight_fit = simple_interpolator_model$variational(
    data = data_stan,
    iter = 1000000,
    draws = 10000, show_messages = T, init = 0, tol_rel_obj = 0.001#, algorithm = "fullrank"
  )


  raw_draws = straight_fit$draws(c('lp__','lp_approx__', 'logcoef1'),format = 'draws_matrix')

  hist(exp(raw_draws[,3]))
  exp(log_coef1)

  raw_draws[,1]

  log_coef1_draws = straight_fit$draws(format = 'draws_matrix')

  # Extract the critical log-density variables evaluated on the unconstrained space
  # lp__         = log p(theta, y) (Target density)
  # lp_approx__  = log q(theta)    (Proposal density)
  log_p <- raw_draws[,1]
  log_q <- raw_draws[,2]

  # 5. Calculate the raw log importance weights
  # log(r_s) = log(p) - log(q)
  raw_log_weights <- log_p - log_q

  # 6. Apply Pareto Smoothed Importance Sampling (PSIS)
  # The loo::psis() function automatically determines the truncation threshold M,
  # fits the GPD to the upper tail, and computes the pareto-k diagnostic.
  psis_result <- loo::psis(log_ratios = raw_log_weights, r_eff = 1)

  # Extract the k-hat diagnostic value
  k_hat <- psis_result$diagnostics$pareto_k
  cat(sprintf("\nPSIS Pareto k-hat diagnostic: %.3f\n", k_hat))

  # Evaluate the reliability of the Variational Approximation
  if (k_hat <= 0.5) {
    cat("Result: k <= 0.5. The variational approximation is reliable and variance is finite.\n")
  } else if (k_hat <= 0.7) {
    cat("Result: 0.5 < k <= 0.7. The approximation has degraded, but the PSIS correction is still capable of unbiased estimation.\n")
  } else {
    warning("Result: k > 0.7. The importance weights exhibit infinite variance. The VI estimate and the PSIS correction are unreliable. Exact HMC sampling is required.")
  }

  # 7. Extract the smoothed, unnormalized log weights from the PSIS object
  # Normalization is deferred to the resampling stage for numerical stability
  smoothed_log_weights <- weights(psis_result, log = TRUE, normalize = FALSE)

  # 8. Correct the Variational Posterior Space
  # Bind the smoothed log weights directly to the raw draws object
  weighted_draws <- posterior::weight_draws(
    x = raw_draws,
    weights = smoothed_log_weights,
    log = TRUE
  )

  # Execute Stratified Resampling
  # This function generates an unweighted draws object that approximates the true
  # posterior much more accurately than the raw ADVI draws, systematically correcting
  # for the variance underestimation caused by the KL-divergence penalty.
  corrected_draws <- posterior::resample_draws(
    x = weighted_draws,
    method = "stratified"
  )

  hist(corrected_draws[,3])

  log_coef1_lower = as.numeric(quantile(exp(log_coef1_draws), p = 0.025))
  log_coef1_upper = as.numeric(quantile(exp(log_coef1_draws), p = 0.975))

  squared_error = (mean(exp(log_coef1_draws)) - exp(log_coef1))^2

  in_CI = log_coef1_upper >= exp(log_coef1) & log_coef1_lower <= exp(log_coef1)

  c(squared_error, in_CI, straight_fit$time()$total)

}

SF_sum = simple_fit$summary()

straight_fit = simple_fit

run_model_HMC = function(){

  log_coef1 = rnorm(1,0,0.35)
  log_coef2 = rnorm(1,0,0.35)

  N_drifter = 100

  start_pos_x = runif(N_drifter, min = -10, max = 10)
  start_pos_y = runif(N_drifter, min = -10, max = 10)

  sigma_vel = 0.1

  b_x_errors = rnorm(N_drifter, 0, sigma_vel)
  b_y_errors = rnorm(N_drifter, 0, sigma_vel)

  t_start = 0
  t_end = 10
  t_res = 100
  t_grid_data = seq(t_start,t_end,length.out = t_res)

  p_x_base = exp(log_coef1) * t_grid_data
  p_y_base = exp(log_coef2) * t_grid_data

  p_x_base_mat = matrix(nrow = t_res, ncol = N_drifter)
  p_y_base_mat = matrix(nrow = t_res, ncol = N_drifter)

  for(d in 1:N_drifter){

    p_x_base_mat[,d] = p_x_base + start_pos_x[d] + b_x_errors[d]*t_grid_data
    p_y_base_mat[,d] = p_y_base + start_pos_y[d] + b_y_errors[d]*t_grid_data


  }

  # Phi_data = bSpline(x = t_grid_data, df = K_drifter, intercept = F)
  #
  p_x_mat = p_x_base_mat# + Phi_data %*% w_x_matrix
  p_y_mat = p_y_base_mat# + Phi_data %*% w_y_matrix

  sigma_pos = 0.1

  p_x_mat_data = p_x_mat + matrix(rnorm(N_drifter*t_res, 0, sigma_pos), nrow = t_res)
  p_y_mat_data = p_y_mat + matrix(rnorm(N_drifter*t_res, 0, sigma_pos), nrow = t_res)

  data_stan = list(
    D = N_drifter,
    K = K_drifter,
    N_drifter = t_res,

    t_grid_data = t_grid_data,
    x_pos = p_x_mat_data,
    y_pos = p_y_mat_data
  )

  straight_fit = line_simple_model$sample(
    data = data_stan,
    chains = 4,
    parallel_chains = 4,
    show_messages = F
  )

  log_coef1_draws = straight_fit$draws('logcoef1', format = 'draws_matrix')
  log_coef1_lower = as.numeric(quantile(exp(log_coef1_draws), p = 0.025))
  log_coef1_upper = as.numeric(quantile(exp(log_coef1_draws), p = 0.975))

  squared_error = (mean(exp(log_coef1_draws)) - exp(log_coef1))^2

  in_CI = log_coef1_upper >= exp(log_coef1) & log_coef1_lower <= exp(log_coef1)

  c(squared_error, in_CI, straight_fit$time()$total)

}

max(SF_sum$rhat)

Spline_V1_Test_VI = data.frame(Sq_Error = rep(0,100), In_CI = rep(0,100), Time = rep(0,100))
i=1
for(i in 2:100){

  Spline_V1_Test_VI[i,] = run_model_VI()

  print(str_c('Done with run ',i, '.'))


}

Spline_V1_Test_VI[1,]
warnings()


mean(Spline_V1_Test_VI$Sq_Error)
mean(Spline_V1_Test_VI$In_CI)
mean(Spline_V1_Test_VI$Time)



