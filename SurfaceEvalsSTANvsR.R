library(cmdstanr)

# 1. Create dummy test data (e.g., 5 simple points)
N_test <- 5
x_test <- c(0, 5, 10, 15, 20)
y_test <- c(0, 5, 10, 15, 20)
t_test <- c(1, 2, 3, 4, 5)

M_Surface <- 10 # or whatever your M is
border <- c(-50, -50, 50, 50) # Your true borders
Lx <- border[3] - border[1]
Ly <- border[4] - border[2]

# Create a fixed matrix of beta coefficients (e.g., all 0.1)
beta_test <- matrix(0.1, nrow = M_Surface + 1, ncol = M_Surface + 1)

test_data <- list(
  N_test = N_test, x_test = x_test, y_test = y_test, t_test = t_test,
  M_Surface = M_Surface, Lx = Lx, Ly = Ly, border = border,
  beta_test = beta_test
)

# 2. Run the Stan Unit Test (1 iteration, no MCMC)
test_model <- cmdstan_model("test_surface.stan")
test_fit <- test_model$sample(data = test_data, fixed_param = TRUE, iter_sampling = 1)

# 3. Extract Stan's Math
hsgp_stan <- as.numeric(test_fit$draws("hsgp_output", format="matrix")[1,])
vf1_stan <- matrix(test_fit$draws("vf_model1_output", format="matrix")[1,], ncol=2)
vf2_stan <- matrix(test_fit$draws("vf_model2_output", format="matrix")[1,], ncol=2)

# =====================================================================
# 4. RUN YOUR R FUNCTIONS HERE
# =====================================================================
hsgp_R <- getPointsStanHSGP_Beta(matrix(rep(0.1,242), nrow = 1), matrix(c(x_test,y_test), ncol = 2, byrow = F), M = 10, n_iters = 1, border = border)
vf1_R <- t(apply(matrix(c(t_test,x_test,y_test), ncol = 3, byrow = F), MARGIN = 1, FUN = function(row){baseVectorFields(row[1], row[-1])}))
c(hsgp_R$traj_1)
# 5. THE MOMENT OF TRUTH
print("Max absolute difference in HSGP:")
print(max(abs(hsgp_stan - hsgp_R$traj_1)))

print("Max absolute difference in Vector Field 1:")
print(max(abs(vf1_stan - vf1_R)))
vf2_stan


##### Pure Surface Test #############

baseVectorFields = function(t, curPos){

  f1 = c(curPos[2],-1*curPos[1])/sqrt(sum(c(curPos[1],curPos[2])^2))
  f2 = c(curPos[1],curPos[2])/sqrt(sum(c(curPos[1],curPos[2])^2))
  matrix(c(f1,f2), nrow = 2, byrow = F)

}


library(cmdstanr)

# =========================================================================
# 1. MAP YOUR SIMULATED DATA HERE
# (Just replace 'sim_data$x' with whatever your true SDE dataframe is called)
# =========================================================================
true_log_k = 0.35
true_log_l = 1
M = 4

sampledHSGP = sampleFullTrajectoriesHSGP(2, M, true_log_k, true_log_l)

sampledHSGP_Plot = plotFullTrajectoriesHSGP(sampledHSGP, grid_res = 100, border = c(-10,-10,10,10), grid_border = c(-10,-10,10,10))

sampledHSGP_Plot[[2]]$Coefficient

sim_points = expand.grid(1, seq(-10,10, length.out = 20), seq(-10,10, length.out = 20))

sim_points = data.frame(t = 1, X1 = sim_data_sub$X1, X2 = sim_data_sub$X2)

sim_vels = apply(sim_points, MARGIN = 1, FUN = function(row){TrajWeightedBaseVectorFields_HSGP(row[1], row[-1], baseVectorFields, sampledHSGP, M = 4, border = c(-10,-10,10,10))})

hist(sim_vels)

rand_res = 100

sim_data = samplePhySpaceParticles(16, startTime = 0, n_obs = 100*rand_res, border = c(-10,-10,10,10), borderBuffer = 0.2, baseVectorFields, sampledHSGP,
                                              M = M, t_step_mean = 0.01, vel_sigma = 0, pos_sigma = 0)

sim_data_sub = sim_data[1:3200 * 50 - (50-1),]

sim_vels = apply(sim_data_sub[,c(1,2,3)], MARGIN = 1, FUN = function(row){TrajWeightedBaseVectorFields_HSGP(row[1], row[-1], baseVectorFields, sampledHSGP, M = 10, border = c(-10,-10,10,10))})

sampledHSGP

plot(sim_points$X1, sim_points$X2)

plot(sim_data_sub$X1, sim_data_sub$X2)

sampledParticles_Sub_V2 = sim_data_sub[sim_data_sub$X1 >= -10 & sim_data_sub$X1 <= 10 &
                                         sim_data_sub$X2 >= -10 & sim_data_sub$X2 <= 10,]

plot(sampledParticles_Sub_V2$X1, sampledParticles_Sub_V2$X2)

true_x <- sim_points[,2]
true_y <- sim_points[,3]
true_vx <- sim_vels[1,] + rnorm(400, 0, 0.2)
true_vy <- sim_vels[2,] + rnorm(400, 0, 0.2)
true_t <- sim_points[,1]

# =========================================================================
# 2. RUN THE TEST MODEL (NO EDITS NEEDED BELOW HERE)
# =========================================================================

# Build the data list (using the exact environment variables you already have)
pure_test_data <- list(
  N_points = length(true_x),
  x = true_x,
  y = true_y,
  v_x = true_vx,
  v_y = true_vy,
  t = true_t,
  N_models = 2,
  M_Surface = 15,  # Uses your existing M_Surface variable
  border = c(-12.5,-12.5,12.5,12.5),        # Uses your existing border variable
  fixed_ks = rep(0.35,2)#,    # Uses your existing fixed_ks variable
  #fixed_ls = rep(0.01,2)     # Uses your existing fixed_ls variable
)

6.4*0.1

pure_test_data <- list(
  N_points = length(true_x),
  x = true_x,
  y = true_y,
  v_x = true_vx,
  v_y = true_vy,
  t = true_t,
  N_models = 2,
  M_Surface = 5,  # Uses your existing M_Surface variable
  border = c(-12,-12,12,12),        # Uses your existing border variable
  fixed_ks = rep(0.35,2)#,    # Uses your existing fixed_ks variable
  #fixed_ls = rep(0.5,2)     # Uses your existing fixed_ls variable
)

pure_test_data_datal <- list(
  N_points = length(true_x),
  x = true_x,
  y = true_y,
  v_x = true_vx,
  v_y = true_vy,
  t = true_t,
  N_models = 2,
  M_Surface = 6,  # Uses your existing M_Surface variable
  border = c(-12,-12,12,12),        # Uses your existing border variable
  fixed_ks = rep(0.35,2),    # Uses your existing fixed_ks variable
  fixed_ls = 0.724      # Uses your existing fixed_ls variable
)

# Compile and run
print("Compiling pure surface test model...")
pure_model_opt <- cmdstan_model("pure_surface_test.stan")

pure_model_full <- cmdstan_model("pure_surface_test_parellel.stan", cpp_options = list(stan_threads = T))
pure_model_full_Datal <- cmdstan_model("pure_surface_test_parellel_Datal.stan", cpp_options = list(stan_threads = T))


pure_fit_optimize = pure_model$optimize(data = pure_test_data, init = 0.01)

pure_fit_optimize$summary()

0.724/sqrt(log(1600))

print("Running ADVI on pure true paths...")
pure_fit <- pure_model_full$variational(
  data = pure_test_data,
  threads = 4,
  iter = 10000,
  draws = 1000, show_messages = T, init = 0
) # Note: Used pathfinder for lightning-fast, highly robust VI

pure_fit_datal <- pure_model_full_Datal$variational(
  data = pure_test_data_datal,
  threads = 4,
  iter = 10000,
  draws = 1000, show_messages = T, init = 0
)

pure_test_data_datal

PF_sum = pure_fit$summary()

PF_sum[PF_sum$variable == 'fixed_ls',]

pure_fit_PF <- pure_model$pathfinder(
  data = pure_test_data,          # Make sure this matches your actual data list name
  num_paths = 4,             # Runs 4 independent optimization paths in parallel
  single_path_draws = 250,   # Draws to take from each path (1000 total)
  max_lbfgs_iters = 2000,    # Gives it plenty of room to find the peak
  history_size = 50          # How much memory the L-BFGS algorithm keeps (helps with complex geometry)
)

pure_fit_MCMC = pure_model$sample(
  data = pure_test_data,          # Make sure this is your filtered, 2-drifter data list
  chains = 4,                # 4 independent Markov chains
  parallel_chains = 4,       # Run them all at the same time on your CPU cores
  iter_warmup = 500,         # Short warmup just to find the typical set
  iter_sampling = 1000       # Short sampling just to get a rough variance estimate

)

# =========================================================================
# 3. EXTRACT AND PLOT
# =========================================================================

# Extract the surface betas
surface_beta_draws_test <- pure_fit$draws(variables = "surface_betas", format = "draws_matrix")
z_surface_beta_draws_test <- pure_fit$draws(variables = "z_surface_betas", format = "draws_matrix")

post_mean_zs = matrix(colMeans(z_surface_beta_draws_test)[1:256 * 2 -1], nrow = 16, byrow = F)
sampledHSGP$log_z[[1]]

sd(pure_fit$draws(variables = "sigma_vel", format = "draws_matrix"))

hist(pure_fit$draws(variables = "sigma_vel", format = "draws_matrix"))

# Use your newly fixed function to generate the surfaces!
# Make sure you have a 'plot_points' grid ready to evaluate the surface on

plot_points = expand.grid(seq(-10,10, length.out = 101), seq(-10,10, length.out = 101))
truth_1 = exp(apply(plot_points, MARGIN = 1, FUN = function(row){evaluateHSGP(z = sampledHSGP$log_z[[1]], k = 0.35, l = 1, M = 4, border = c(-10,-10,10,10), curPos = row)}))
truth_2 = exp(apply(plot_points, MARGIN = 1, FUN = function(row){evaluateHSGP(z = sampledHSGP$log_z[[2]], k = 0.35, l = 1, M = 4, border = c(-10,-10,10,10), curPos = row)}))

truth_1 - truth[[1]]$CoefValues

truth[[1]]$Coefficient

post_mean = exp(apply(plot_points, MARGIN = 1, FUN = function(row){evaluateHSGP(post_mean_zs, k = 0.35, l = 10, M = 15, border = c(-12.5,-12.5,12.5,12.5), curPos = row)}))

median(abs(truth - post_mean))

test_results <- getPointsStanHSGP_Beta(
  surface_beta_draws = surface_beta_draws_test,
  points = plot_points,  # Your spatial grid for plotting
  M = 5,
  n_iters = nrow(surface_beta_draws_test),
  border = c(-12,-12,12,12)
)

trueVF = apply(plot_points, MARGIN = 1, FUN = function(row){TrajWeightedBaseVectorFields_HSGP(1, row, baseVectorFields, sampledHSGP, M = 4, border = c(-10,-10,10,10))})

baseVF = t(apply(plot_points, 1, FUN = function(row){c(baseVectorFields(1, row))}))

post_xv = baseVF[,1] * test_results$traj_1 + baseVF[,3] * test_results$traj_2
post_yv = baseVF[,2] * test_results$traj_1 + baseVF[,4] * test_results$traj_2

mean_xv = rowMeans(post_xv)
lower_xv = apply(post_xv, MARGIN = 1, FUN = quantile, probs = 0.025)
upper_xv = apply(post_xv, MARGIN = 1, FUN = quantile, probs = 0.975)

mean(trueVF[1,] >= lower_xv & trueVF[1,] <= upper_xv)

hist(upper_xv - lower_xv)


hist((mean_xv - trueVF[1,]))

mean_yv = rowMeans(post_yv)
lower_yv = apply(post_yv, MARGIN = 1, FUN = quantile, probs = 0.025)
upper_yv = apply(post_yv, MARGIN = 1, FUN = quantile, probs = 0.975)

mean(trueVF[2,] >= lower_yv & trueVF[2,] <= upper_yv)

hist(upper_yv - lower_yv)
hist((mean_yv - trueVF[2,]))

hist(post_xv)
hist(sim_vels[1,])

c(baseVectorFields(1, test_results$points[10,]))

hist(post_xv - sim_vels[1,])

warnings()
# Calculate the Posterior Mean Surface for Model 1
mean_surface_1 <- rowMeans(test_results$traj_1)

lower_surface_1 = apply(test_results$traj_1, MARGIN = 1, FUN = quantile, probs = 0.025)
upper_surface_1 = apply(test_results$traj_1, MARGIN = 1, FUN = quantile, probs = 0.975)

hist(upper_surface_1 - lower_surface_1)

median(abs(mean_surface_1 - truth))

hist(mean_surface_1 - truth_1)

mean(truth_1 >= lower_surface_1 & truth_1 <= upper_surface_1)

mean(truth_1 < lower_surface_1)

mean(truth_1 > upper_surface_1)



mean_surface_2 <- rowMeans(test_results$traj_2)

lower_surface_2 = apply(test_results$traj_2, MARGIN = 1, FUN = quantile, probs = 0.025)
upper_surface_2 = apply(test_results$traj_2, MARGIN = 1, FUN = quantile, probs = 0.975)

hist(upper_surface_2 - lower_surface_2)

median(abs(mean_surface_2 - truth_2))

hist(mean_surface_2 - truth_2)

mean(truth_2 >= lower_surface_2 & truth_2 <= upper_surface_2)

mean(truth_2 < lower_surface_2)

mean(truth_2 > upper_surface_2)


# Quick and dirty base plot to immediately see if it worked
plot(plot_points[,1], plot_points[,2],
     col = hcl.colors(100, "Viridis")[cut(mean_xv - trueVF[1,], 100)],
     pch = 15, cex = 1.5,
     main = "Pure Surface Test: Model 1 HSGP Field")



hist(mean_surface_2 - truth_2)

truth = plotFullTrajectoriesHSGP(sampledHSGP, grid_res = 100, color_limits = c(0,4), border = c(-10,-10,10,10), grid_border = c(-10,-10,10,10))

FullTrajPost = getGridStanHSGP_Beta(surface_beta_draws = surface_beta_draws_test, grid_res = 100, M = 5, n_iters = 1000, border = c(-12,-12,12,12), grid_border = c(-10,-10,10,10))
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

  # # Layer 1: The Lower Bound (Floor)
  add_surface(x = x_seq, y = y_seq, z = FullTraj2_Lower,
              opacity = 0.3,           # <--- Make it ghostly
              colorscale = "Viridis",
              cmin = c_min, cmax = c_max,
              showscale = FALSE) %>%   # Hide legend for this layer
  #
  # # Layer 2: The Upper Bound (Ceiling)
  add_surface(x = x_seq, y = y_seq, z = FullTraj2_Upper,
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
  add_surface(x = x_seq, y = y_seq, z = truth[[2]]$CoefValues,
              opacity = 1.0,           # <--- Solid
              colorscale = "Viridis",
              cmin = c_min, cmax = c_max,
              colorbar = list(title = "Trajectory")) %>%

  layout(title = "Truth with 95% Credible Envelope",
         scene = list(
           zaxis = list(title = "Trajectory Value")
         ))



### Putting back in the path

sampledParticles = samplePhySpaceParticles(100, startTime = 0, n_obs = 100*rand_res, border = c(-10,-10,10,10), borderBuffer = 0.2, baseVectorFields, sampledHSGP,
                                           M = M, t_step_mean = 0.01, vel_sigma = 0, pos_sigma = 0)

sampledParticles_Sub = sampledParticles[1:10000 * rand_res - (rand_res-1),]

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

  LI_x = approx(x = curDrifter$t, y = curDrifter$X1, xout = cur_quad_grid)$y
  LI_y = approx(x = curDrifter$t, y = curDrifter$X2, xout = cur_quad_grid)$y

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
  M_Surface = 15,

  border = c(-12.5,-12.5,12.5,12.5),
  fixed_ks = c(0.35,0.35),
  fixed_ls = c(0.01, 0.01)

)

new_interpolator_model = cmdstan_model('/home/kdp2abu/NewInterpolationModel.stan', cpp_options = list(stan_threads = T))

full_fit = new_interpolator_model$variational(
  data = full_data,
  threads = max_cores/4,
  iter = 100000,
  draws = 1000, show_messages = T, init = 0#, algorithm = "fullrank"
)

surface_beta_draws_test <- full_fit$draws(variables = "surface_betas", format = "draws_matrix")
z_surface_beta_draws_test <- full_fit$draws(variables = "z_surface_betas", format = "draws_matrix")

post_mean_zs = matrix(colMeans(z_surface_beta_draws_test)[1:256 * 2 -1], nrow = 16, byrow = F)
sampledHSGP$log_z[[1]]

full_fit$draws(variables = "sigma_vel", format = "draws_matrix")

# Use your newly fixed function to generate the surfaces!
# Make sure you have a 'plot_points' grid ready to evaluate the surface on

plot_points = expand.grid(seq(-10,10, length.out = 100), seq(-10,10, length.out = 100))
truth = exp(apply(plot_points, MARGIN = 1, FUN = function(row){evaluateHSGP(sampledHSGP$log_z[[1]], k = 0.35, l = 0.05, M = 10, border = c(-10,-10,10,10), curPos = row)}))
post_mean = exp(apply(plot_points, MARGIN = 1, FUN = function(row){evaluateHSGP(post_mean_zs, k = 0.35, l = 0.01, M = 15, border = c(-12.5,-12.5,12.5,12.5), curPos = row)}))

mean(abs(truth - post_mean))

test_results <- getPointsStanHSGP_Beta(
  surface_beta_draws = surface_beta_draws_test,
  points = plot_points,  # Your spatial grid for plotting
  M = 15,
  n_iters = nrow(surface_beta_draws_test),
  border = c(-12.5,-12.5,12.5,12.5)
)
warnings()
# Calculate the Posterior Mean Surface for Model 1
mean_surface_1 <- rowMeans(test_results$traj_1)
median(abs(mean_surface_1 - truth))

hist(test_results$traj_1[45,])

# Quick and dirty base plot to immediately see if it worked
plot(plot_points[,1], plot_points[,2],
     col = hcl.colors(100, "Viridis")[cut(abs(mean_surface_1 - truth), 100)],
     pch = 15, cex = 1.5,
     main = "Pure Surface Test: Model 1 HSGP Field")


# 1. Throw it into a quick dataframe
plot_df <- data.frame(
  x = plot_points[, 1],
  y = plot_points[, 2],
  error = log((mean_surface_1 - truth)^2)
)

# 2. Plot with the automatic continuous scale
ggplot(plot_df, aes(x = x, y = y, color = error)) +
  geom_point(shape = 15, size = 3) +
  scale_color_viridis_c(option = "viridis") +
  coord_fixed() + # Keeps your spatial map from stretching
  theme_minimal() +
  labs(title = "Model 1 HSGP Field: Absolute Deviation from Truth",
       color = "Absolute\nError")
