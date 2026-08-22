source('SecondProjectRequiredFunctions.R')
options(mc.cores = parallel::detectCores())
max_cores = parallel::detectCores()
library(splines2)
library(VFmixing)


#Gameplan

# 1. Simulate Data with positional error (easy)
# 2. From simulated data, sample errors to add to the OG data
# 3. Optimize for best BASE path (overspecified basis and desiring the lowest velocity path through all the data points)
# 4. Sample prior w's from prior distribution (Just N(0,1) for now)
# 5. With augmented data and sampled priors, find optimal path abiding by the current trajectory distribution (start with 1,1 constant)
# 6. Calculate Metropolis ratio... not sure how exactly ... and accept/reject proposed sample
# 7. Repeat 2-6 for as many samples as possible


dipoleVectorField = function(t, curPos) {
  x = curPos[1]
  y = curPos[2]

  speed = 1
  epsilon = 0.000001

  # Calculate the perfectly circular directional vectors
  dx = -2 * x * y
  dy = x^2 - y^2

  # Normalize to a unit vector, then multiply by desired speed
  magnitude_sq = x^2 + y^2 + epsilon

  vx_dipole = (dx / magnitude_sq) * speed
  vy_dipole = (dy / magnitude_sq) * speed

  vx_jet = 0
  vy_jet = -speed

  blend_weight = exp(-magnitude_sq / 0.1)

  f1 = c(vx_dipole, vy_dipole) * (1 - blend_weight)
  f2 = c(vx_jet, vy_jet) * blend_weight

  matrix(c(f1,f2), nrow = 2, byrow = F)
}

OneTraj = list(Traj = list(matrix(c(0,0,0,1,0,0,0,1), nrow = 2, byrow = T)))


# Step 1

OneTraj = list(log_z = list(matrix(rep(0,4), nrow=2), matrix(rep(0,4), nrow=2)), log_k = c(0,0), log_l = c(100,100))

plotHSGP(grid_res = 100, z_log = OneTraj$log_z[[1]], k_log = 0, l_log = 100, M = 1, color_limits = c(0,2), border = c(0,0,1,1), grid_border = c(0,0,1,1))

rand_res = 100

sampledParticles = samplePhySpaceParticles(n_particles = 1, startTime = 0, n_obs = 100*rand_res, border = c(-1,-1,1,1), borderBuffer = 0.4, dipoleVectorField, OneTraj,
                                              M = 2, t_step_mean = 0.01, vel_sigma = 0, pos_sigma = 0)

test_particle = sampledParticles[1:100 * rand_res - (rand_res-1),1:3]

test_particle = RungeKutta(startTime = 0, startPos = c(-0.1, -0.1), baseVectorFields = dipoleVectorField, TrajList = OneTraj$Traj, TrajTimeSplits = c(0,10), endTime = 10, t_step = 0.01)


ggplot(sampledParticles_Sub, aes(x = X1, y = X2)) + geom_point()

# Step 2 (basically)

true_pos_sd = 0.05

sim_particle = test_particle + matrix(c(rep(0, nrow(test_particle)), rnorm(2*nrow(test_particle), 0, true_pos_sd)), byrow = F, ncol = 3)

N = 20

sim_particle_sub = sim_particle[seq(1, nrow(sim_particle), length.out = N),]
ggplot(sim_particle, aes(x = t, y = X1)) + geom_point()

sim_particle_sub = sim_particle

# Step 3

K = 10
N = nrow(sim_particle_sub)
drifter_boundary = c(min(sim_particle_sub$t), max(sim_particle_sub$t))
Lt = diff(drifter_boundary)
c = 1.25
N_quad = 1000
t_quad = seq(drifter_boundary[1], drifter_boundary[2], length.out = N_quad)

drifter_boundary_extended = c(drifter_boundary[1] - (c-1)/2 * Lt, drifter_boundary[2] + (c-1)/2 * Lt)


Phi_data = matrix(nrow = N, ncol = N+K)
Phi = matrix(nrow = N_quad, ncol = N+K)
Phi_d = matrix(nrow = N_quad, ncol = N+K)
Phi_a = matrix(nrow = N_quad, ncol = N+K)

for(i in 1:(N+K)){

  Phi_data[,i] = cos((i-1)*pi*(sim_particle_sub$t - drifter_boundary_extended[1]) / (Lt*c))
  Phi[,i] = cos((i-1)*pi*(t_quad - drifter_boundary_extended[1]) / (Lt*c))
  Phi_d[,i] = -1*((i-1)*pi / (Lt*c))*sin((i-1)*pi*(t_quad - drifter_boundary_extended[1]) / (Lt*c))
  Phi_a[,i] = -1*((i-1)*pi / (Lt*c))^2*cos((i-1)*pi*(t_quad - drifter_boundary_extended[1]) / (Lt*c))

}

c_base_x = ginv(Phi_data) %*% sim_particle_sub$X1
c_base_y = ginv(Phi_data) %*% sim_particle_sub$X2

Z = Null(t(Phi_data))

c_base_x_wig = c_base_x + Z %*% matrix(rnorm(K, 0, 1))
c_base_y_wig = c_base_y + Z %*% matrix(rnorm(K, 0, 1))

Phi %*% c_base_x

ggplot() + geom_point(data = sim_particle_sub, aes(x = t, y = X2)) + geom_line(aes(x = t_quad, y = Phi %*% (c_base_y + Z %*% matrix(rnorm(K, 0, 1)))))

opt_trans_mat_pos = -1*solve(t(Z) %*% t(Phi) %*% Phi %*% Z) %*% t(Z) %*% t(Phi) %*% Phi
opt_trans_mat_vel = -1*solve(t(Z) %*% t(Phi_d) %*% Phi_d %*% Z) %*% t(Z) %*% t(Phi_d) %*% Phi_d
opt_trans_mat_acc = -1*solve(t(Z) %*% t(Phi_a) %*% Phi_a %*% Z) %*% t(Z) %*% t(Phi_a) %*% Phi_a

w_opt_x_pos = opt_trans_mat_pos %*% c_base_x
w_opt_y_pos = opt_trans_mat_pos %*% c_base_y

w_opt_x_vel = opt_trans_mat_vel %*% c_base_x
w_opt_y_vel = opt_trans_mat_vel %*% c_base_y

w_opt_x_acc = opt_trans_mat_acc %*% c_base_x
w_opt_y_acc = opt_trans_mat_acc %*% c_base_y

c_opt_x_pos = c_base_x + Z %*% w_opt_x_pos
c_opt_y_pos = c_base_y + Z %*% w_opt_y_pos

c_opt_x_vel = c_base_x + Z %*% w_opt_x_vel
c_opt_y_vel = c_base_y + Z %*% w_opt_y_vel

c_opt_x_acc = c_base_x + Z %*% w_opt_x_acc
c_opt_y_acc = c_base_y + Z %*% w_opt_y_acc

ggplot() + geom_point(data = sim_particle_sub, aes(x = t, y = X1)) + geom_line(aes(x = t_quad, y = Phi %*% c_base_x))
ggplot() + geom_point(data = sim_particle_sub, aes(x = t, y = X1)) + geom_line(aes(x = t_quad, y = Phi %*% c_opt_x_pos))
ggplot() + geom_point(data = sim_particle_sub, aes(x = t, y = X1)) + geom_line(aes(x = t_quad, y = Phi %*% c_opt_x_vel))
ggplot() + geom_point(data = sim_particle_sub, aes(x = t, y = X1)) + geom_line(aes(x = t_quad, y = Phi %*% c_opt_x_acc))

ggplot() + geom_point(data = sim_particle_sub, aes(x = t, y = X2)) + geom_line(aes(x = t_quad, y = Phi %*% c_base_y))
ggplot() + geom_point(data = sim_particle_sub, aes(x = t, y = X2)) + geom_line(aes(x = t_quad, y = Phi %*% c_opt_y_pos))
ggplot() + geom_point(data = sim_particle_sub, aes(x = t, y = X2)) + geom_line(aes(x = t_quad, y = Phi %*% c_opt_y_vel))
ggplot() + geom_point(data = sim_particle_sub, aes(x = t, y = X2)) + geom_line(aes(x = t_quad, y = Phi %*% c_opt_y_acc))

ggplot() + geom_point(data = sim_particle_sub, aes(x = X1, y = X2)) + geom_path(aes(x = Phi %*% c_base_x, y = Phi %*% c_base_y))
ggplot() + geom_point(data = sim_particle_sub, aes(x = X1, y = X2)) + geom_path(aes(x = Phi %*% c_opt_x_pos, y = Phi %*% c_opt_y_pos))
ggplot() + geom_point(data = sim_particle_sub, aes(x = X1, y = X2)) + geom_path(aes(x = Phi %*% c_opt_x_vel, y = Phi %*% c_opt_y_vel))
ggplot() + geom_point(data = sim_particle_sub, aes(x = X1, y = X2)) + geom_path(aes(x = Phi %*% c_opt_x_acc, y = Phi %*% c_opt_y_acc))



t(Phi_a %*% c_base_x) %*% (Phi_a %*% c_base_x)
t(Phi_a %*% c_opt_x_pos) %*% (Phi_a %*% c_opt_x_pos)
t(Phi_a %*% c_opt_x_vel) %*% (Phi_a %*% c_opt_x_vel)
t(Phi_a %*% c_opt_x_acc) %*% (Phi_a %*% c_opt_x_acc)


# Step 4

overall_variance = 10
Sigma_Prior = solve(t(Z) %*% t(Phi_a) %*% Phi_a %*% Z)

w_samp_x = mvrnorm(mu = rep(0, K), Sigma = N_quad * (overall_variance / K) * Sigma_Prior)
w_samp_y = mvrnorm(mu = rep(0, K), Sigma = N_quad * (overall_variance / K) * Sigma_Prior)


ggplot() + geom_point(data = sim_particle_sub, aes(x = t, y = X1)) + geom_line(aes(x = t_quad, y = Phi %*% c_opt_x_acc))
ggplot() + geom_point(data = sim_particle_sub, aes(x = t, y = X1)) + geom_line(aes(x = t_quad, y = Phi %*% (c_opt_x_acc + Z %*% w_samp_x)))

ggplot() + geom_point(data = sim_particle_sub, aes(x = t, y = X2)) + geom_line(aes(x = t_quad, y = Phi %*% c_opt_y_acc))
ggplot() + geom_point(data = sim_particle_sub, aes(x = t, y = X2)) + geom_line(aes(x = t_quad, y = Phi %*% (c_opt_y_acc + Z %*% w_samp_y)))

ggplot() + geom_point(data = sim_particle_sub, aes(x = X1, y = X2)) + geom_path(aes(x = Phi %*% c_opt_x_acc, y = Phi %*% c_opt_y_acc))
ggplot() + geom_point(data = sim_particle_sub, aes(x = X1, y = X2)) + geom_path(aes(x = Phi %*% (c_opt_x_acc + Z %*% w_samp_x), y = Phi %*% (c_opt_y_acc + Z %*% w_samp_y)))



overall_variance = 100
added_acc = rep(0,10000)

for(i in 1:10000){

  w_samp = mvrnorm(mu = rep(0, K), Sigma = N_quad * (overall_variance / K) * Sigma_Prior)

  added_acc[i] = t(Phi_a %*% Z %*% w_samp_x) %*% (Phi_a %*% Z %*% w_samp_x) / N_quad


}

hist(added_acc)
mean(added_acc)
sd(added_acc)


# Step 5

sigma_vel = 0.1

dipoleVectorField_vectorized_add = function(input_mat) {
  # input_mat is an n x 3 matrix. We extract all x and y values at once.
  # t is input_mat[, 1], but it is unused in the math.
  x = input_mat[, 2]
  y = input_mat[, 3]

  speed = 1
  epsilon = 0.000001

  # Element-wise math on the entire vectors of x and y
  dx = -2 * x * y
  dy = x^2 - y^2

  magnitude_sq = x^2 + y^2 + epsilon

  vx_dipole = (dx / magnitude_sq) * speed
  vy_dipole = (dy / magnitude_sq) * speed

  blend_weight = exp(-magnitude_sq / 0.1)

  # Calculate f1 for all n points
  f1_x = vx_dipole * (1 - blend_weight)
  f1_y = vy_dipole * (1 - blend_weight)

  # Calculate f2 for all n points
  # vx_jet is 0, so f2_x is always 0
  f2_x = rep(0, length(x))
  f2_y = -speed * blend_weight

  # Combine into an n x 2 x 2 array.
  # For any point i, the result is accessed via result[i, , ]

  matrix(c(f1_x + f2_x, f1_y + f2_y), byrow = F, ncol = 1)
}


plot_VF_grid = expand.grid(0, seq(-0.5,0.5,length.out = 100), seq(-0.5,0.5,length.out = 100))

VF_values = dipoleVectorField_vectorized_add(plot_VF_grid)

plot_VF_values = matrix(VF_values, ncol = 2, byrow = F)

ggplot() + geom_point(aes(x = plot_VF_grid[,2], y = plot_VF_grid[,3], color = sqrt(plot_VF_values[,1]^2 + plot_VF_values[,2]^2)))



c_rand_start_x = c_opt_x_acc #+ Z %*% w_samp_x
c_rand_start_y = c_opt_y_acc #+ Z %*% w_samp_y

path_vel = rbind(Phi_d %*% c_rand_start_x, Phi_d %*% c_rand_start_y)

traj_vel = dipoleVectorField_vectorized_add(cbind(t_quad, Phi %*% c_rand_start_x, Phi %*% c_rand_start_y))


likelihood_loss = (1 / sigma_vel^2) * (t(path_vel - traj_vel) %*% (path_vel - traj_vel))

inv_Sigma_Prior = (K / (N_quad * overall_variance)) * (t(Z) %*% t(Phi_a) %*% Phi_a %*% Z)

prior_loss = t(w_samp_x) %*% inv_Sigma_Prior %*% w_samp_x + t(w_samp_y) %*% inv_Sigma_Prior %*% w_samp_y



sigma_vel = 0.1
inv_Sigma_Prior = (K / (N_quad * overall_variance)) * (t(Z) %*% t(Phi_a) %*% Phi_a %*% Z)

rand_w_x = mvrnorm(mu = rep(0, K), Sigma = N_quad * (overall_variance / K) * Sigma_Prior)
rand_w_y = mvrnorm(mu = rep(0, K), Sigma = N_quad * (overall_variance / K) * Sigma_Prior)

c_rand_start_x = c_opt_x_acc + Z %*% rand_w_x
c_rand_start_y = c_opt_y_acc + Z %*% rand_w_y

rMAP_loss = function(w){

  w_x = w[1:K]
  w_y = w[1:K + K]

  c_x = c_rand_start_x + Z %*% w_x
  c_y = c_rand_start_y + Z %*% w_y

  path_vel = rbind(Phi_d %*% c_x, Phi_d %*% c_y)
  traj_vel = dipoleVectorField_vectorized_add(cbind(t_quad, Phi %*% c_x, Phi %*% c_y))

  likelihood_loss = (1 / sigma_vel^2) * (t(path_vel - traj_vel) %*% (path_vel - traj_vel))
  prior_loss = t(w_x - rand_w_x) %*% inv_Sigma_Prior %*% (w_x - rand_w_x) + t(w_y - rand_w_y) %*% inv_Sigma_Prior %*% (w_y - rand_w_y)

  as.numeric(likelihood_loss + prior_loss)

}


w_init <- c(rand_w_x, rand_w_y)

# Run the optimizer
opt_result <- optim(
  par = w_init,             # Starting values
  fn = rMAP_loss,           # Your loss function
  method = "BFGS",          # Unconstrained quasi-Newton
  control = list(
    maxit = 1000,           # Increase max iterations (default is 100)
    trace = 0,              # Prints progress to the console
    REPORT = 10             # How often to print progress
  )
)


w_opt_full = opt_result$value

w_opt_x = w_opt_full[1:K]
w_opt_y = w_opt_full[1:K + K]

c_opt_x = c_rand_start_x + Z %*% w_opt_x
c_opt_y = c_rand_start_y + Z %*% w_opt_y


opt_pos_x = Phi %*% c_opt_x
opt_pos_y = Phi %*% c_opt_y

ggplot() + geom_point(data = sim_particle_sub, aes(x = X1, y = X2)) + geom_path(aes(x = opt_pos_x, y = opt_pos_y))


path_vel = rbind(Phi_d %*% c_opt_x, Phi_d %*% c_opt_y)
traj_vel = dipoleVectorField_vectorized_add(cbind(t_quad, Phi %*% c_opt_x, Phi %*% c_opt_y))

(t(path_vel - traj_vel) %*% (path_vel - traj_vel)) / (2 * N_quad)



# Create Function to do all steps

sim_data = sim_particle_sub

sample_rMAP = function(sim_data, pos_sd = 0.1, K = 10, N_quad = 1000, add_acc_var = 10, vel_sd = 0.1){

  N = nrow(sim_data)

  # Add positional error

  x_pos_err = rnorm(N, 0, pos_sd)
  y_pos_err = rnorm(N, 0, pos_sd)

  aug_data = sim_data + matrix(c(rep(0, N), x_pos_err, y_pos_err), byrow = F, ncol = 3)

  # Get optimal base path

  drifter_boundary = c(min(aug_data$t), max(aug_data$t))
  Lt = diff(drifter_boundary)
  c = 1.25

  t_quad = seq(drifter_boundary[1], drifter_boundary[2], length.out = N_quad)
  drifter_boundary_extended = c(drifter_boundary[1] - (c-1)/2 * Lt, drifter_boundary[2] + (c-1)/2 * Lt)

  # Propagate the basis matrices

  Phi_data = matrix(nrow = N, ncol = N+K)
  Phi = matrix(nrow = N_quad, ncol = N+K)
  Phi_d = matrix(nrow = N_quad, ncol = N+K)
  Phi_a = matrix(nrow = N_quad, ncol = N+K)

  for(i in 1:(N+K)){

    Phi_data[,i] = cos((i-1)*pi*(aug_data$t - drifter_boundary_extended[1]) / (Lt*c))
    Phi[,i] = cos((i-1)*pi*(t_quad - drifter_boundary_extended[1]) / (Lt*c))
    Phi_d[,i] = -1*((i-1)*pi / (Lt*c))*sin((i-1)*pi*(t_quad - drifter_boundary_extended[1]) / (Lt*c))
    Phi_a[,i] = -1*((i-1)*pi / (Lt*c))^2*cos((i-1)*pi*(t_quad - drifter_boundary_extended[1]) / (Lt*c))

  }

  # Get generic interpolator from ginv

  c_base_x = ginv(Phi_data) %*% aug_data$X1
  c_base_y = ginv(Phi_data) %*% aug_data$X2

  # Compute the null space

  Z = Null(t(Phi_data))

  # Find matrix to transform the generic interpolator to the one with the lowest acceleration

  helper_M = t(Z) %*% t(Phi_a) %*% Phi_a %*% Z

  helper_M_inv = ginv(helper_M)

  opt_trans_mat_acc = -1*helper_M_inv %*% t(Z) %*% t(Phi_a) %*% Phi_a

  # Find the optimal w

  w_opt_x_acc = opt_trans_mat_acc %*% c_base_x
  w_opt_y_acc = opt_trans_mat_acc %*% c_base_y

  #Get optimal c

  c_opt_x_acc = c_base_x + Z %*% w_opt_x_acc
  c_opt_y_acc = c_base_y + Z %*% w_opt_y_acc

  ggplot() + geom_point(data = aug_data, aes(x = t, y = X1)) + geom_path(aes(x = t_quad, y = Phi %*% c_opt_x_acc))


  # Sample random w as a particular starting point from the prior
  # Sigma comes from the added acceleration being distribution Chi-Squared

  w_prior_sigma = (N_quad * (add_acc_var / K)) * helper_M_inv
  inv_w_prior_sigma = (K / (N_quad * add_acc_var)) * helper_M

  rand_w_x = mvrnorm(mu = rep(0, K), Sigma = w_prior_sigma)
  rand_w_y = mvrnorm(mu = rep(0, K), Sigma = w_prior_sigma)

  c_rand_start_x = c_opt_x_acc + Z %*% rand_w_x
  c_rand_start_y = c_opt_y_acc + Z %*% rand_w_y

  # Optimize for the best w that matched the mixing fields

  rMAP_loss = function(w){

    w_x = w[1:K]
    w_y = w[1:K + K]

    c_x = c_rand_start_x + Z %*% w_x
    c_y = c_rand_start_y + Z %*% w_y

    path_vel = rbind(Phi_d %*% c_x, Phi_d %*% c_y)
    traj_vel = dipoleVectorField_vectorized_add(cbind(t_quad, Phi %*% c_x, Phi %*% c_y))

    likelihood_loss = (1 / vel_sd^2) * (t(path_vel - traj_vel) %*% (path_vel - traj_vel))
    prior_loss = t(w_x - rand_w_x) %*% inv_w_prior_sigma %*% (w_x - rand_w_x) + t(w_y - rand_w_y) %*% inv_w_prior_sigma %*% (w_y - rand_w_y)

    as.numeric(likelihood_loss + prior_loss)

  }

  w_init <- c(rand_w_x, rand_w_y)

  # Run the optimizer
  opt_result <- optim(
    par = w_init,             # Starting values
    fn = rMAP_loss,           # Your loss function
    method = "BFGS",          # Unconstrained quasi-Newton
    control = list(
      maxit = 1000,           # Increase max iterations (default is 100)
      trace = 1,              # Prints progress to the console
      REPORT = 10             # How often to print progress
    )
  )

  post_w_x = opt_result$par[1:K]
  post_w_y = opt_result$par[1:K + K]

  pos_x_post = Phi %*% (c_rand_start_x + Z %*% post_w_x)
  pos_y_post = Phi %*% (c_rand_start_y + Z %*% post_w_y)

  list(pos_x_post, pos_y_post, opt_result$value)

}


pos_x_post = matrix(nrow = N_quad, ncol = 100)
pos_y_post = matrix(nrow = N_quad, ncol = 100)
post_like = rep(0,100)

svMisc::progress(0,100)

for(i in 1:100){

  cur_draw = sample_rMAP(sim_particle_sub)

  pos_x_post[,i] = cur_draw[[1]]
  pos_y_post[,i] = cur_draw[[2]]

  post_like[i] = cur_draw[[3]]

  svMisc::progress(i,100)

}


plotting_Post_df = data.frame(t = rep(t_quad, 100), X1 = c(pos_x_post), X2 = c(pos_y_post), run = rep(1:100, each = N_quad))

ggplot() + geom_point(data = test_particle, aes(x = X1, y = X2), color = 'red') + geom_path(data = plotting_Post_df, aes(x = X1, y = X2, group = run), alpha = 0.2)


ggplot() + geom_path(data = test_particle, aes(x = X1, y = X2), color = 'red') + geom_path(aes(x = rowMeans(pos_x_post), y = rowMeans(pos_y_post)))

ggplot() + geom_line(data = test_particle, aes(x = t, y = X1), color = 'red') + #geom_line(aes(x = t_quad, y = rowMeans(pos_x_post)), color = 'black') +
  geom_line(aes(x = t_quad, y = apply(pos_x_post, 1, quantile, p = 0.025)), color = 'blue') + geom_line(aes(x = t_quad, y = apply(pos_x_post, 1, quantile, p = 0.975)), color = 'blue')

ggplot() + geom_line(data = test_particle, aes(x = t, y = X2), color = 'red') + #geom_line(aes(x = t_quad, y = rowMeans(pos_y_post)), color = 'black') +
  geom_line(aes(x = t_quad, y = apply(pos_y_post, 1, quantile, p = 0.025)), color = 'blue') + geom_line(aes(x = t_quad, y = apply(pos_y_post, 1, quantile, p = 0.975)), color = 'blue')



ggplot() + geom_path(data = test_particle, aes(x = X1, y = X2), color = 'red') + #geom_line(aes(x = t_quad, y = rowMeans(pos_x_post)), color = 'black') +
  geom_path(aes(x = apply(pos_x_post, 1, quantile, p = 0.025), y = apply(pos_y_post, 1, quantile, p = 0.025)), color = 'blue') +
  geom_path(aes(x = apply(pos_x_post, 1, quantile, p = 0.975), y = apply(pos_y_post, 1, quantile, p = 0.975)), color = 'blue')








## Test with all origin points

N_quad = 1000
end_t = 20

all_origin = data.frame(t = seq(0,end_t, length.out = 20), X1 = rep(0,20), X2 = rep(0,20))


pos_x_post_origin = matrix(nrow = N_quad, ncol = 100)
pos_y_post_origin = matrix(nrow = N_quad, ncol = 100)
post_like_origin = rep(0,100)

svMisc::progress(0,100)

for(i in 1:100){

  cur_draw = sample_rMAP_nointerp(all_origin, pos_sd = 0.001, K_base = 1, K_full = 30, N_quad = N_quad, vel_sd = 0.1, pos_selection_sd = 0.1)

  pos_x_post_origin[,i] = cur_draw[[1]]
  pos_y_post_origin[,i] = cur_draw[[2]]

  post_like_origin[i] = cur_draw[[3]] + cur_draw[[4]] + cur_draw[[5]]

  svMisc::progress(i,100)

}

t_quad_origin = seq(0,end_t,length.out = N_quad)

plotting_Post_Origin_df = data.frame(t = rep(t_quad_origin, 100), X1 = c(pos_x_post_origin), X2 = c(pos_y_post_origin), run = rep(1:100, each = N_quad))

ggplot() + geom_path(data = plotting_Post_Origin_df, aes(x = X1, y = X2, group = run), alpha = 0.2) + geom_point(data = all_origin, aes(x = X1, y = X2), color = 'red')# + xlim(-0.5,0.5) + ylim(-0.5,0.5)
ggplot() + geom_path(data = plotting_Post_Origin_df, aes(x = t, y = X1, group = run), alpha = 0.2) + geom_point(data = all_origin, aes(x = t, y = X1), color = 'red')
ggplot() + geom_path(data = plotting_Post_Origin_df, aes(x = t, y = X2, group = run), alpha = 0.2) + geom_point(data = all_origin, aes(x = t, y = X2), color = 'red')



ggplot() + geom_path(data = test_particle, aes(x = X1, y = X2), color = 'red') + geom_path(aes(x = rowMeans(pos_x_post), y = rowMeans(pos_y_post)))
ggplot() + geom_path(data = test_particle, aes(x = t, y = X1), color = 'red') + geom_path(aes(x = t_quad, y = rowMeans(pos_x_post)))
ggplot() + geom_path(data = test_particle, aes(x = t, y = X2), color = 'red') + geom_path(aes(x = t_quad, y = rowMeans(pos_y_post)))

ggplot() + geom_line(data = test_particle, aes(x = t, y = X1), color = 'red') + #geom_line(aes(x = t_quad, y = rowMeans(pos_x_post)), color = 'black') +
  geom_line(aes(x = t_quad, y = apply(pos_x_post, 1, quantile, p = 0.025)), color = 'blue') + geom_line(aes(x = t_quad, y = apply(pos_x_post, 1, quantile, p = 0.975)), color = 'blue')

ggplot() + geom_line(data = test_particle, aes(x = t, y = X2), color = 'red') + #geom_line(aes(x = t_quad, y = rowMeans(pos_y_post)), color = 'black') +
  geom_line(aes(x = t_quad, y = apply(pos_y_post, 1, quantile, p = 0.025)), color = 'blue') + geom_line(aes(x = t_quad, y = apply(pos_y_post, 1, quantile, p = 0.975)), color = 'blue')



ggplot() + geom_path(data = test_particle, aes(x = X1, y = X2), color = 'red') + #geom_line(aes(x = t_quad, y = rowMeans(pos_x_post)), color = 'black') +
  geom_path(aes(x = apply(pos_x_post, 1, quantile, p = 0.025), y = apply(pos_y_post, 1, quantile, p = 0.025)), color = 'blue') +
  geom_path(aes(x = apply(pos_x_post, 1, quantile, p = 0.975), y = apply(pos_y_post, 1, quantile, p = 0.975)), color = 'blue')




K_base = 50
K_full = 50
# Real Tests

sample_rMAP_nointerp = function(sim_data, pos_sd = 0.1, K_base = 10, K_full = 20, N_quad = 1000, vel_sd = 0.1, pos_selection_sd = 0.1){

  N = nrow(sim_data)

  # Add positional error

  x_pos_err = rnorm(N, 0, pos_sd)
  y_pos_err = rnorm(N, 0, pos_sd)

  aug_data = sim_data + matrix(c(rep(0, N), x_pos_err, y_pos_err), byrow = F, ncol = 3)

  # Get optimal base path

  drifter_boundary = c(min(aug_data$t), max(aug_data$t))
  Lt = diff(drifter_boundary)
  c = 1.25

  t_quad = seq(drifter_boundary[1], drifter_boundary[2], length.out = N_quad)
  drifter_boundary_extended = c(drifter_boundary[1] - (c-1)/2 * Lt, drifter_boundary[2] + (c-1)/2 * Lt)

  # Propagate the basis matrices

  Phi_data = matrix(nrow = N, ncol = K_full)
  Phi = matrix(nrow = N_quad, ncol = K_full)
  Phi_d = matrix(nrow = N_quad, ncol = K_full)

  for(i in 1:K_full){

    Phi_data[,i] = cos((i-1)*pi*(aug_data$t - drifter_boundary_extended[1]) / (Lt*c))
    Phi[,i] = cos((i-1)*pi*(t_quad - drifter_boundary_extended[1]) / (Lt*c))
    Phi_d[,i] = -1*((i-1)*pi / (Lt*c))*sin((i-1)*pi*(t_quad - drifter_boundary_extended[1]) / (Lt*c))

  }

  # Get least-squares estimate from first K_base parameters

  Phi_data_base = Phi_data[,1:K_base]

  c_base_x = c(ginv(t(Phi_data_base) %*% Phi_data_base) %*% t(Phi_data_base) %*% aug_data$X1, rep(0, K_full-K_base))
  c_base_y = c(ginv(t(Phi_data_base) %*% Phi_data_base) %*% t(Phi_data_base) %*% aug_data$X2, rep(0, K_full-K_base))

  # ggplot() + geom_point(data = aug_data, aes(x = t, y = X1)) + geom_path(aes(x = t_quad, y = Phi %*% c_base_x))


  # Sample random w as a particular starting point from the prior
  # Sigma comes from the added acceleration being distribution Chi-Squared

  helper_M = t(Phi) %*% Phi
  helper_M_inv = ginv(helper_M)

  c_prior_sigma = (N_quad * pos_selection_sd^2 / Lt) * helper_M_inv
  inv_c_prior_sigma = (Lt / (N_quad * pos_selection_sd^2)) * helper_M

  rand_c_x = c_base_x + mvrnorm(mu = rep(0, K_full), Sigma = c_prior_sigma)
  rand_c_y = c_base_y + mvrnorm(mu = rep(0, K_full), Sigma = c_prior_sigma)

  # Optimize for the best w that matched the mixing fields

  rMAP_loss = function(c){

    c_x = c[1:K_full]
    c_y = c[1:K_full + K_full]

    full_c_x = rand_c_x + c_x
    full_c_y = rand_c_y + c_y

    path_vel = rbind(Phi_d %*% full_c_x, Phi_d %*% full_c_y)
    traj_vel = dipoleVectorField_vectorized_add(cbind(t_quad, Phi %*% full_c_x, Phi %*% full_c_y))

    path_pos_data_x = Phi_data %*% full_c_x
    path_pos_data_y = Phi_data %*% full_c_y


    likelihood_loss_pos = (1 / pos_sd^2) * (t(path_pos_data_x - aug_data$X1) %*% (path_pos_data_x - aug_data$X1) + t(path_pos_data_y - aug_data$X2) %*% (path_pos_data_y - aug_data$X2)) / (2*N)
    likelihood_loss_vel = (1 / vel_sd^2) * (t(path_vel - traj_vel) %*% (path_vel - traj_vel)) / (2*N_quad)
    prior_loss = t(c_x) %*% inv_c_prior_sigma %*% (c_x) + t(c_y) %*% inv_c_prior_sigma %*% (c_y)

    as.numeric(likelihood_loss_pos + likelihood_loss_vel + prior_loss)

  }

  w_init <- c(rand_c_x, rand_c_y)

  # Run the optimizer
  opt_result <- optim(
    par = w_init,             # Starting values
    fn = rMAP_loss,           # Your loss function
    method = "BFGS",          # Unconstrained quasi-Newton
    control = list(
      maxit = 1000,           # Increase max iterations (default is 100)
      trace = 0,              # Prints progress to the console
      REPORT = 10             # How often to print progress
    )
  )

  post_c_x = opt_result$par[1:K_full]
  post_c_y = opt_result$par[1:K_full + K_full]

  post_c_x_full = rand_c_x + post_c_x
  post_c_y_full = rand_c_y + post_c_y

  post_path_vel = rbind(Phi_d %*% post_c_x_full, Phi_d %*% post_c_y_full)
  post_traj_vel = dipoleVectorField_vectorized_add(cbind(t_quad, Phi %*% post_c_x_full, Phi %*% post_c_y_full))

  post_path_pos_data_x = Phi_data %*% post_c_x_full
  post_path_pos_data_y = Phi_data %*% post_c_y_full


  post_likelihood_loss_pos = (1 / pos_sd^2) * (t(post_path_pos_data_x - aug_data$X1) %*% (post_path_pos_data_x - aug_data$X1) + t(post_path_pos_data_y - aug_data$X2) %*% (post_path_pos_data_y - aug_data$X2)) / (2*N)
  post_likelihood_loss_vel = (1 / vel_sd^2) * (t(post_path_vel - post_traj_vel) %*% (post_path_vel - post_traj_vel)) / (2*N_quad)
  post_prior_loss = t(post_c_x) %*% inv_c_prior_sigma %*% (post_c_x) + t(post_c_y) %*% inv_c_prior_sigma %*% (post_c_y)

  pos_x_post = Phi %*% post_c_x_full
  pos_y_post = Phi %*% post_c_y_full

  list(pos_x_post, pos_y_post, post_likelihood_loss_pos, post_likelihood_loss_vel, post_prior_loss)

}

sample_rMAP_bspline = function(sim_data, pos_sd = 0.1, K_base = 10, K_full = 20, N_quad = 1000, vel_sd = 0.1, pos_selection_sd = 0.1){

  N = nrow(sim_data)

  # Add positional error

  x_pos_err = rnorm(N, 0, pos_sd)
  y_pos_err = rnorm(N, 0, pos_sd)

  aug_data = sim_data + matrix(c(rep(0, N), x_pos_err, y_pos_err), byrow = F, ncol = 3)

  # Get optimal base path

  drifter_boundary = c(min(aug_data$t), max(aug_data$t))
  Lt = diff(drifter_boundary)
  t_quad = seq(drifter_boundary[1], drifter_boundary[2], length.out = N_quad)

  # Propagate the basis matrices

  Phi_data_base = bSpline(aug_data$t, df = K_base, intercept = T)
  Phi_base = bSpline(t_quad, df = K_base, intercept = T)
  Phi_d_base = bSpline(t_quad, df = K_base, intercept = T, derivs = 1)

  Phi_data = bSpline(aug_data$t, df = K_full, intercept = T)
  Phi = bSpline(t_quad, df = K_full, intercept = T)
  Phi_d = bSpline(t_quad, df = K_full, intercept = T, derivs = 1)

  # Get least-squares estimate from first K_base parameters

  c_base_x = c(ginv(t(Phi_data_base) %*% Phi_data_base) %*% t(Phi_data_base) %*% aug_data$X1)
  c_base_y = c(ginv(t(Phi_data_base) %*% Phi_data_base) %*% t(Phi_data_base) %*% aug_data$X2)

  # ggplot() + geom_point(data = aug_data, aes(x = t, y = X1)) + geom_path(aes(x = t_quad, y = Phi_base %*% c_base_x))


  # Sample random w as a particular starting point from the prior
  # Sigma comes from the added acceleration being distribution Chi-Squared

  helper_M = t(Phi) %*% Phi
  helper_M_inv = ginv(helper_M)

  c_prior_sigma = (N_quad * pos_selection_sd^2 / Lt) * helper_M_inv
  inv_c_prior_sigma = (Lt / (N_quad * pos_selection_sd^2)) * helper_M

  rand_c_x = mvrnorm(mu = rep(0, K_full), Sigma = c_prior_sigma)
  rand_c_y = mvrnorm(mu = rep(0, K_full), Sigma = c_prior_sigma)

  rand_vel = rnorm(2*N_quad, 0, sd = vel_sd)

  # Optimize for the best w that matched the mixing fields

  rMAP_loss = function(c){

    c_x = c[1:K_full]
    c_y = c[1:K_full + K_full]

    path_vel = rbind(Phi_d_base %*% c_base_x + Phi_d %*% c_x, Phi_d_base %*% c_base_y + Phi_d %*% c_y)
    traj_vel = dipoleVectorField_vectorized_add(cbind(t_quad, Phi_base %*% c_base_x + Phi %*% c_x, Phi_base %*% c_base_y + Phi %*% c_y))

    path_pos_data_x = Phi_data_base %*% c_base_x + Phi_data %*% c_x
    path_pos_data_y = Phi_data_base %*% c_base_y + Phi_data %*% c_y


    likelihood_loss_pos = (1 / pos_sd^2) * (t(path_pos_data_x - aug_data$X1) %*% (path_pos_data_x - aug_data$X1) + t(path_pos_data_y - aug_data$X2) %*% (path_pos_data_y - aug_data$X2)) / (2*N)
    likelihood_loss_vel = (1 / vel_sd^2) * (t(path_vel - traj_vel - rand_vel) %*% (path_vel - traj_vel - rand_vel)) / (2*N_quad)
    prior_loss = t(c_x - rand_c_x) %*% inv_c_prior_sigma %*% (c_x - rand_c_x) + t(c_y - rand_c_y) %*% inv_c_prior_sigma %*% (c_y - rand_c_y)

    as.numeric(likelihood_loss_pos + likelihood_loss_vel + prior_loss)

  }

  w_init <- c(rand_c_x, rand_c_y)

  # Run the optimizer
  opt_result <- optim(
    par = w_init,             # Starting values
    fn = rMAP_loss,           # Your loss function
    method = "BFGS",          # Unconstrained quasi-Newton
    control = list(
      maxit = 1000,           # Increase max iterations (default is 100)
      trace = 1,              # Prints progress to the console
      REPORT = 10             # How often to print progress
    )
  )

  post_c_x = opt_result$par[1:K_full]
  post_c_y = opt_result$par[1:K_full + K_full]

  post_path_vel = rbind(Phi_d_base %*% c_base_x + Phi_d %*% post_c_x, Phi_d_base %*% c_base_y + Phi_d %*% post_c_y)
  post_traj_vel = dipoleVectorField_vectorized_add(cbind(t_quad, Phi_base %*% c_base_x + Phi %*% post_c_x, Phi_base %*% c_base_y + Phi %*% post_c_y))

  post_path_pos_data_x = Phi_data_base %*% c_base_x + Phi_data %*% post_c_x
  post_path_pos_data_y = Phi_data_base %*% c_base_y + Phi_data %*% post_c_y


  post_likelihood_loss_pos = (1 / pos_sd^2) * (t(post_path_pos_data_x - aug_data$X1) %*% (post_path_pos_data_x - aug_data$X1) + t(post_path_pos_data_y - aug_data$X2) %*% (post_path_pos_data_y - aug_data$X2)) / (2*N)
  post_likelihood_loss_vel = (1 / vel_sd^2) * (t(post_path_vel - post_traj_vel - rand_vel) %*% (post_path_vel - post_traj_vel - rand_vel)) / (2*N_quad)
  post_prior_loss = t(post_c_x - rand_c_x) %*% inv_c_prior_sigma %*% (post_c_x - rand_c_x) + t(post_c_y - rand_c_y) %*% inv_c_prior_sigma %*% (post_c_y - rand_c_y)

  pos_x_post = Phi_base %*% c_base_x + Phi %*% post_c_x
  pos_y_post = Phi_base %*% c_base_y + Phi %*% post_c_y

  plot(pos_x_post, pos_y_post)

  list(pos_x_post, pos_y_post, post_likelihood_loss_pos, post_likelihood_loss_vel, post_prior_loss)

}


#Maybe some aliasing

OneTraj = list(log_z = list(matrix(rep(0,4), nrow=2), matrix(rep(0,4), nrow=2)), log_k = c(0,0), log_l = c(100,100))

plotHSGP(grid_res = 100, z_log = OneTraj$log_z[[1]], k_log = 0, l_log = 100, M = 1, color_limits = c(0,2), border = c(0,0,1,1), grid_border = c(0,0,1,1))

rand_res = 100

sampledParticles = samplePhySpaceParticles(n_particles = 1, startTime = 0, n_obs = 100*rand_res, border = c(-0.5,-0.5,0.5, 0.5), borderBuffer = 0.4, dipoleVectorField, OneTraj,
                                           M = 2, t_step_mean = 0.01, vel_sigma = 0.1, pos_sigma = 0.1)

test_particle = sampledParticles[1:100 * rand_res - (rand_res-1),1:3]

test_particle = RungeKutta(startTime = 0, startPos = runif(2,0,0.1), baseVectorFields = dipoleVectorField, TrajList = OneTraj$Traj, TrajTimeSplits = c(0,50), endTime = 10, t_step = 0.01)

true_pos_sd = 0.01

sim_particle = test_particle + matrix(c(rep(0, nrow(test_particle)), rnorm(2*nrow(test_particle), 0, true_pos_sd)), byrow = F, ncol = 3)

N = 50

sim_particle_sub = sim_particle[seq(1, nrow(sim_particle), length.out = N),]
ggplot(test_particle, aes(x = t, y = X1)) + geom_point()

N_quad = 1000

pos_x_post = matrix(nrow = N_quad, ncol = 100)
pos_y_post = matrix(nrow = N_quad, ncol = 100)
post_like_pos = rep(0,100)
post_like_vel = rep(0,100)
post_like_prior = rep(0,100)

svMisc::progress(0,100)

sim_particle_sub = test_particle

for(i in 1:100){

  cur_draw = sample_rMAP_bspline(sim_data = sim_particle_sub, pos_sd = 0.1, K_base = 75, K_full = 75, N_quad = N_quad, vel_sd = 0.1, pos_selection_sd = 2)

  #Metropolis

  if(i > 1){

    prev_like = post_like_pos[i-1] + post_like_vel[i-1] + post_like_prior[i-1]
    cur_like = cur_draw[[3]] + cur_draw[[4]] + cur_draw[[5]]

    a = cur_like / prev_like

    u = runif(1)

    if(u < a){

      pos_x_post[,i] = cur_draw[[1]]
      pos_y_post[,i] = cur_draw[[2]]

      post_like_pos[i] = cur_draw[[3]]
      post_like_vel[i] = cur_draw[[4]]
      post_like_prior[i] = cur_draw[[5]]

      print('Accepted')

    } else{

      pos_x_post[,i] = pos_x_post[,i-1]
      pos_y_post[,i] = pos_y_post[,i-1]

      post_like_pos[i] = post_like_pos[i-1]
      post_like_vel[i] = post_like_pos[i-1]
      post_like_prior[i] = post_like_pos[i-1]

      print('Rejected')

    }

  } else{

    pos_x_post[,i] = cur_draw[[1]]
    pos_y_post[,i] = cur_draw[[2]]

    post_like_pos[i] = cur_draw[[3]]
    post_like_vel[i] = cur_draw[[4]]
    post_like_prior[i] = cur_draw[[5]]

  }

  svMisc::progress(i,100)

}

plot(post_like_pos+post_like_vel+post_like_prior)

median(post_like_vel)

t_quad = seq(min(sim_particle_sub$t),max(sim_particle_sub$t), length.out = N_quad)


plotting_Post_df = data.frame(t = rep(t_quad, 100), X1 = c(pos_x_post), X2 = c(pos_y_post), run = rep(1:100, each = N_quad))

test_particle = sampledParticles

ggplot() + geom_path(data = test_particle, aes(x = X1, y = X2), color = 'red', size = 1) + geom_path(data = plotting_Post_df, aes(x = X1, y = X2, group = run), alpha = 0.2) + geom_point(data = sim_particle_sub, aes(x = X1, y = X2), color = 'red')
ggplot() + geom_path(data = test_particle, aes(x = t, y = X1), color = 'red') + geom_path(data = plotting_Post_df, aes(x = t, y = X1, group = run), alpha = 0.2) + geom_point(data = sim_particle_sub, aes(x = t, y = X1), color = 'red')
ggplot() + geom_path(data = test_particle, aes(x = t, y = X2), color = 'red') + geom_path(data = plotting_Post_df, aes(x = t, y = X2, group = run), alpha = 0.2) + geom_point(data = sim_particle_sub, aes(x = t, y = X2), color = 'red')



ggplot() + geom_path(data = test_particle, aes(x = X1, y = X2), color = 'red') + geom_path(aes(x = rowMeans(pos_x_post), y = rowMeans(pos_y_post)))
ggplot() + geom_path(data = test_particle, aes(x = t, y = X1), color = 'red') + geom_path(aes(x = t_quad, y = rowMeans(pos_x_post)))
ggplot() + geom_path(data = test_particle, aes(x = t, y = X2), color = 'red') + geom_path(aes(x = t_quad, y = rowMeans(pos_y_post)))

ggplot() + geom_line(data = test_particle, aes(x = t, y = X1), color = 'red') + #geom_line(aes(x = t_quad, y = rowMeans(pos_x_post)), color = 'black') +
  geom_line(aes(x = t_quad, y = apply(pos_x_post, 1, quantile, p = 0.025)), color = 'blue') + geom_line(aes(x = t_quad, y = apply(pos_x_post, 1, quantile, p = 0.975)), color = 'blue')

ggplot() + geom_line(data = test_particle, aes(x = t, y = X2), color = 'red') + #geom_line(aes(x = t_quad, y = rowMeans(pos_y_post)), color = 'black') +
  geom_line(aes(x = t_quad, y = apply(pos_y_post, 1, quantile, p = 0.025)), color = 'blue') + geom_line(aes(x = t_quad, y = apply(pos_y_post, 1, quantile, p = 0.975)), color = 'blue')



ggplot() + geom_path(data = test_particle, aes(x = X1, y = X2), color = 'red') + #geom_line(aes(x = t_quad, y = rowMeans(pos_x_post)), color = 'black') +
  geom_path(aes(x = apply(pos_x_post, 1, quantile, p = 0.025), y = apply(pos_y_post, 1, quantile, p = 0.025)), color = 'blue') +
  geom_path(aes(x = apply(pos_x_post, 1, quantile, p = 0.975), y = apply(pos_y_post, 1, quantile, p = 0.975)), color = 'blue')



######### Simple VFs

evaluate2DCosine = function(beta_mat, pos_mat, border){

  Lx = border[3] - border[1]
  Ly = border[4] - border[2]

  M = sqrt(nrow(beta_mat))-1
  N_points = nrow(pos_mat)

  omega = (0:M)*pi

  Phi = matrix(nrow = N_points, ncol = (M+1)^2)

  for(i in 0:M){
    for(j in 0:M){

      Phi[,i*(M+1) + j + 1] = sqrt(2)^(i!=0) * sqrt(2)^(j!=0) * cos(omega[i+1]*(pos_mat[,1] - border[1])/Lx) * cos(omega[j+1]*(pos_mat[,2] - border[2])/Ly)

    }
  }


  Phi %*% beta_mat

}

pos_mat = expand.grid(seq(0.01,1, length.out = 10), seq(0.01,1, length.out = 10))
pos_t_mat = cbind(rep(0,100), pos_mat)
evaluate2DCosine(beta_mat = matrix(c(1,0.3,0.6,-2,
                          2,0.8,-3,4), nrow = 4, byrow = F), pos_mat, border = c(0,0,1,1))


TrajWeightedBaseVectorFields_2D_Cosine = function(pos_t_mat, beta_mat, baseVectorFields_Vec, border){


  log_Traj = evaluate2DCosine(beta_mat = beta_mat, pos_mat = pos_t_mat[,c(2,3)], border = border)
  log_Traj_mat = cbind(log_Traj, log_Traj)

  baseVF = baseVectorFields_Vec(pos_t_mat)

  Weighted_Vel = exp(log_Traj_mat) * baseVF

  cbind(Weighted_Vel[,1] + Weighted_Vel[,2], Weighted_Vel[,3] + Weighted_Vel[,4])


}



baseVectorFields_Vec = function(pos_t_mat){

  c = sqrt(rowSums(pos_t_mat[,c(2,3)]^2))

  f1x = pos_t_mat[,3] / c
  f1y = -1*pos_t_mat[,2] / c

  f2x = pos_t_mat[,2] / c
  f2y = pos_t_mat[,3] / c

  cbind(f1x,f2x, f1y, f2y)

}


sample_rMAP_bspline = function(sim_data, pos_sd = 0.1, K_base = 10, K_full = 20, N_quad = 1000, vel_sd = 0.1, pos_selection_sd = 0.1, baseVectorFields_Vec, beta_mat, border){

  N = nrow(sim_data)

  # Add positional error

  x_pos_err = rnorm(N, 0, pos_sd)
  y_pos_err = rnorm(N, 0, pos_sd)

  aug_data = sim_data + matrix(c(rep(0, N), x_pos_err, y_pos_err), byrow = F, ncol = 3)

  # Get optimal base path

  drifter_boundary = c(min(aug_data$t), max(aug_data$t))
  Lt = diff(drifter_boundary)
  t_quad = seq(drifter_boundary[1], drifter_boundary[2], length.out = N_quad)

  # Propagate the basis matrices

  Phi_data_base = bSpline(aug_data$t, df = K_base, intercept = T)
  Phi_base = bSpline(t_quad, df = K_base, intercept = T)
  Phi_d_base = bSpline(t_quad, df = K_base, intercept = T, derivs = 1)

  Phi_data = bSpline(aug_data$t, df = K_full, intercept = T)
  Phi = bSpline(t_quad, df = K_full, intercept = T)
  Phi_d = bSpline(t_quad, df = K_full, intercept = T, derivs = 1)

  # Get least-squares estimate from first K_base parameters

  c_base_x = c(ginv(t(Phi_data_base) %*% Phi_data_base) %*% t(Phi_data_base) %*% aug_data$X1)
  c_base_y = c(ginv(t(Phi_data_base) %*% Phi_data_base) %*% t(Phi_data_base) %*% aug_data$X2)

  ggplot() + geom_point(data = aug_data, aes(x = t, y = X1)) + geom_path(aes(x = t_quad, y = Phi_base %*% c_base_x))


  # Sample random w as a particular starting point from the prior
  # Sigma comes from the added acceleration being distribution Chi-Squared

  helper_M = t(Phi) %*% Phi
  helper_M_inv = ginv(helper_M)

  c_prior_sigma = (N_quad * pos_selection_sd^2 / Lt) * helper_M_inv
  inv_c_prior_sigma = (Lt / (N_quad * pos_selection_sd^2)) * helper_M

  rand_c_x = mvrnorm(mu = rep(0, K_full), Sigma = c_prior_sigma)
  rand_c_y = mvrnorm(mu = rep(0, K_full), Sigma = c_prior_sigma)

  rand_vel = rnorm(2*N_quad, 0, sd = vel_sd)

  # Optimize for the best w that matched the mixing fields
  rMAP_loss = function(c){

    c_x = c[1:K_full]
    c_y = c[1:K_full + K_full]

    path_pos_data_x = Phi_data_base %*% c_base_x + Phi_data %*% c_x
    path_pos_data_y = Phi_data_base %*% c_base_y + Phi_data %*% c_y

    path_vel = rbind(Phi_d_base %*% c_base_x + Phi_d %*% c_x, Phi_d_base %*% c_base_y + Phi_d %*% c_y)

    traj_vel = c(TrajWeightedBaseVectorFields_2D_Cosine(pos_t_mat = cbind(t_quad, Phi_base %*% c_base_x + Phi %*% c_x, Phi_base %*% c_base_y + Phi %*% c_y),
                                                      baseVectorFields_Vec = baseVectorFields_Vec, beta_mat = beta_mat, border = border))

    likelihood_loss_pos = (1 / pos_sd^2) * (t(path_pos_data_x - aug_data$X1) %*% (path_pos_data_x - aug_data$X1) + t(path_pos_data_y - aug_data$X2) %*% (path_pos_data_y - aug_data$X2)) / (2*N)
    likelihood_loss_vel = (1 / vel_sd^2) * (t(path_vel - traj_vel - rand_vel) %*% (path_vel - traj_vel - rand_vel)) / (2*N_quad)
    prior_loss = t(c_x - rand_c_x) %*% inv_c_prior_sigma %*% (c_x - rand_c_x) + t(c_y - rand_c_y) %*% inv_c_prior_sigma %*% (c_y - rand_c_y)

    as.numeric(likelihood_loss_pos + likelihood_loss_vel + prior_loss)

  }

  w_init <- c(rand_c_x, rand_c_y)

  # Run the optimizer
  opt_result <- optim(
    par = w_init,             # Starting values
    fn = rMAP_loss,           # Your loss function
    method = "BFGS",          # Unconstrained quasi-Newton
    control = list(
      maxit = 1000,           # Increase max iterations (default is 100)
      trace = 1,              # Prints progress to the console
      REPORT = 10             # How often to print progress
    )
  )

  post_c_x = opt_result$par[1:K_full]
  post_c_y = opt_result$par[1:K_full + K_full]

  post_path_vel = rbind(Phi_d_base %*% c_base_x + Phi_d %*% post_c_x, Phi_d_base %*% c_base_y + Phi_d %*% post_c_y)
  post_traj_vel = c(TrajWeightedBaseVectorFields_2D_Cosine(pos_t_mat = cbind(t_quad, Phi_base %*% c_base_x + Phi %*% post_c_x, Phi_base %*% c_base_y + Phi %*% post_c_y),
                                                           baseVectorFields_Vec = baseVectorFields_Vec, beta_mat = beta_mat, border = border))

  post_path_pos_data_x = Phi_data_base %*% c_base_x + Phi_data %*% post_c_x
  post_path_pos_data_y = Phi_data_base %*% c_base_y + Phi_data %*% post_c_y


  post_likelihood_loss_pos = (1 / pos_sd^2) * (t(post_path_pos_data_x - aug_data$X1) %*% (post_path_pos_data_x - aug_data$X1) + t(post_path_pos_data_y - aug_data$X2) %*% (post_path_pos_data_y - aug_data$X2)) / (2*N)
  post_likelihood_loss_vel = (1 / vel_sd^2) * (t(post_path_vel - post_traj_vel - rand_vel) %*% (post_path_vel - post_traj_vel - rand_vel)) / (2*N_quad)
  post_prior_loss = t(post_c_x - rand_c_x) %*% inv_c_prior_sigma %*% (post_c_x - rand_c_x) + t(post_c_y - rand_c_y) %*% inv_c_prior_sigma %*% (post_c_y - rand_c_y)

  pos_x_post = Phi_base %*% c_base_x + Phi %*% post_c_x
  pos_y_post = Phi_base %*% c_base_y + Phi %*% post_c_y

  plot(pos_x_post, pos_y_post)

  list(pos_x_post, pos_y_post, post_likelihood_loss_pos, post_likelihood_loss_vel, post_prior_loss)

}





baseVectorFields = function(t, curPos){

  c = sqrt(sum(curPos^2))

  f1 = c(curPos[2],-1*curPos[1]) / c
  f2 = c(curPos[1],curPos[2]) / c

  matrix(c(f1,f2), nrow = 2, byrow = F)

}

rand_res = 100
OneTraj = list(log_z = list(matrix(rep(0,4), nrow=2), matrix(rep(0,4), nrow=2)), log_k = c(0,0), log_l = c(100,100))

sqrt(0.01)

sampledParticles = samplePhySpaceParticles(n_particles = 1, startTime = 0, n_obs = 100*10, border = c(-5,-5,5,5), borderBuffer = 0.4, baseVectorFields, OneTraj,
                                           M = 2, t_step_mean = 0.001, vel_sigma = 0.1, pos_sigma = 0)

ggplot(sampledParticles, aes(x = X1, y = X2, color = Particle)) + geom_point()

sampledParticles_Sub = sampledParticles[1:100 * 10 - (10-1),]

ggplot(sampledParticles_Sub, aes(x = X1, y = X2)) + geom_point()


N_quad = 1000

pos_x_post = matrix(nrow = N_quad, ncol = 100)
pos_y_post = matrix(nrow = N_quad, ncol = 100)
post_like_pos = rep(0,100)
post_like_vel = rep(0,100)
post_like_prior = rep(0,100)

svMisc::progress(0,100)

beta_list = OneTraj$log_z

for(i in 1:100){

  cur_draw = sample_rMAP_bspline(sim_data = sampledParticles_Sub[,1:3], pos_sd = 0.01, K_base = 10, K_full = 10, N_quad = N_quad, vel_sd = 0.01, pos_selection_sd = 0.1, beta_mat = matrix(rep(0,8), nrow = 4), baseVectorFields_Vec = baseVectorFields_Vec, border = c(-5,-5,5,5))

  #Metropolis

  pos_x_post[,i] = cur_draw[[1]]
  pos_y_post[,i] = cur_draw[[2]]

  post_like_pos[i] = cur_draw[[3]]
  post_like_vel[i] = cur_draw[[4]]
  post_like_prior[i] = cur_draw[[5]]


  svMisc::progress(i,100)

}


plot(post_like_pos+post_like_vel+post_like_prior)

median(post_like_pos)

t_quad = seq(min(sampledParticles_Sub$t),max(sampledParticles_Sub$t), length.out = N_quad)


plotting_Post_df = data.frame(t = rep(t_quad, 100), X1 = c(pos_x_post), X2 = c(pos_y_post), run = rep(1:100, each = N_quad))

test_particle = sampledParticles
sim_particle_sub = sampledParticles_Sub

ggplot() + geom_path(data = test_particle, aes(x = X1, y = X2), color = 'red', size = 1) + geom_path(data = plotting_Post_df, aes(x = X1, y = X2, group = run), alpha = 0.2) + geom_point(data = sim_particle_sub, aes(x = X1, y = X2), color = 'red')
ggplot() + geom_path(data = test_particle, aes(x = t, y = X1), color = 'red') + geom_path(data = plotting_Post_df, aes(x = t, y = X1, group = run), alpha = 0.2) + geom_point(data = sim_particle_sub, aes(x = t, y = X1), color = 'red')
ggplot() + geom_path(data = test_particle, aes(x = t, y = X2), color = 'red') + geom_path(data = plotting_Post_df, aes(x = t, y = X2, group = run), alpha = 0.2) + geom_point(data = sim_particle_sub, aes(x = t, y = X2), color = 'red')



ggplot() + geom_path(data = test_particle, aes(x = X1, y = X2), color = 'red') + geom_path(aes(x = rowMeans(pos_x_post), y = rowMeans(pos_y_post)))
ggplot() + geom_path(data = test_particle, aes(x = t, y = X1), color = 'red') + geom_path(aes(x = t_quad, y = rowMeans(pos_x_post)))
ggplot() + geom_path(data = test_particle, aes(x = t, y = X2), color = 'red') + geom_path(aes(x = t_quad, y = rowMeans(pos_y_post)))

ggplot() + geom_line(data = test_particle, aes(x = t, y = X1), color = 'red') + #geom_line(aes(x = t_quad, y = rowMeans(pos_x_post)), color = 'black') +
  geom_line(aes(x = t_quad, y = apply(pos_x_post, 1, quantile, p = 0.025)), color = 'blue') + geom_line(aes(x = t_quad, y = apply(pos_x_post, 1, quantile, p = 0.975)), color = 'blue')

ggplot() + geom_line(data = test_particle, aes(x = t, y = X2), color = 'red') + #geom_line(aes(x = t_quad, y = rowMeans(pos_y_post)), color = 'black') +
  geom_line(aes(x = t_quad, y = apply(pos_y_post, 1, quantile, p = 0.025)), color = 'blue') + geom_line(aes(x = t_quad, y = apply(pos_y_post, 1, quantile, p = 0.975)), color = 'blue')



ggplot() + geom_path(data = test_particle, aes(x = X1, y = X2), color = 'red') + #geom_line(aes(x = t_quad, y = rowMeans(pos_x_post)), color = 'black') +
  geom_path(aes(x = apply(pos_x_post, 1, quantile, p = 0.025), y = apply(pos_y_post, 1, quantile, p = 0.025)), color = 'blue') +
  geom_path(aes(x = apply(pos_x_post, 1, quantile, p = 0.975), y = apply(pos_y_post, 1, quantile, p = 0.975)), color = 'blue')





#rMAP but we go for the trajectories themselves

#1. Add positional error to the data
#2. Get the basis to evaluate the trajectory surface at a given position
#3. Add velocity error to the mixing fields by randomly sampling the diffusion terms, but keeping them fixed for optimization purposes
#4. The optimization function is going to be only positional error differences
#     > How well can the mixing vector field advect the particle to the next point in space
#5. The starting trajectory comes from an HSGP, but then we immediately extract the beta
#   and work in the beta space
#6. Optimize for the betas with the rMAP approach
#7. After getting multiple trajectory samples, we can get optimal paths which interpolate
#   (maybe not interpolate) that align with trajectory functions from the posterior
#8. This is a trajectory then path approach, so no need to iterate
#9. Lord help me and let this work PLEASE


sampledParticles = samplePhySpaceParticles(n_particles = 100, startTime = 0, n_obs = 20*10, border = c(-5,-5,5,5), borderBuffer = 0.4, baseVectorFields, OneTraj,
                                           M = 2, t_step_mean = 0.001, vel_sigma = 0.1, pos_sigma = 0)

sampledParticles_Sub = sampledParticles[1:(100*20) * 10 - (10-1),]
ggplot(sampledParticles_Sub, aes(x = X1, y = X2, color = Particle)) + geom_point()

sim_data_list = list()

for(d in 1:100){

  curDrifter = sampledParticles_Sub[sampledParticles_Sub$Particle == str_c('Particle',d),]

  sim_data_list[[d]] = curDrifter[,c(1,2,3)]


}

sim_data = sampledParticles_Sub[,1:3]
pos_sd = 0.01
vel_sd = 0.1
N_prop_steps = 10
border = c(-5,-5,5,5)

traj_eval_grid = expand.grid(seq(-5,5, length.out = 100), seq(-5,5, length.out = 100))

#Should work and is pretty quick...

#Need to sample lots of drifters around the border to get a better reading on the trajectory


sample_rMAP_trajectories = function(sim_data_list, pos_sd = 0.001, vel_sd = 0.1, M = 2, prior_k = 0.35, prior_l = 1, baseVectorFields_Vec, border, N_prop_steps = 10, traj_eval_grid){

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

  start_beta_1 = mvrnorm(1, rep(0, M^2), Sigma = prior_beta_sigma)
  start_beta_2 = mvrnorm(1, rep(0, M^2), Sigma = prior_beta_sigma)

  start_beta = c(start_beta_1, start_beta_2)

  # Optimize for the best w that matched the mixing fields
  rMAP_loss = function(beta){

    M_sq <- M^2
    beta_1 <- beta[1:M_sq]
    beta_2 <- beta[(M_sq + 1):(2 * M_sq)]
    beta_mat <- cbind(beta_1, beta_2)

    # Pre-calculate time step vectors for all i
    t_steps <- (aug_data_ends$t - aug_data_starts$t) / N_prop_steps
    sqrt_t_steps <- sqrt(t_steps)

    # Initialize the starting positions for all (N-1) segments at once
    # curPos_mat is an (N-1) x 3 matrix
    curPos_mat <- aug_data_starts

    # Inner loop: Propagate all segments simultaneously
    for(j in 1:N_prop_steps) {

      # Evaluate drift for all (N-1) positions simultaneously
      # Ensure this function is updated to return an (N-1) x 2 matrix
      cur_drift_mat <- TrajWeightedBaseVectorFields_2D_Cosine(
        pos_t_mat = curPos_mat,
        beta_mat = beta_mat,
        baseVectorFields_Vec = baseVectorFields_Vec,
        border = border
      )

      # Extract diffusion for all N-1 segments for step j
      cur_diffusion_vec_1 <- rand_vel_1[j,]
      cur_diffusion_vec_2 <- rand_vel_2[j,]

      # Update Time (Column 1)
      curPos_mat[, 1] <- curPos_mat[, 1] + t_steps

      # Update X and Y positions (Columns 2 and 3)
      curPos_mat[, 2] <- curPos_mat[, 2] + cur_drift_mat[, 1] * t_steps + cur_diffusion_vec_1 * sqrt_t_steps
      curPos_mat[, 3] <- curPos_mat[, 3] + cur_drift_mat[, 2] * t_steps + cur_diffusion_vec_2 * sqrt_t_steps
    }

    # Fast Likelihood calculation using sum of squares instead of matrix multiplication
    diff_1 <- curPos_mat[, 2] - aug_data_ends$X1
    diff_2 <- curPos_mat[, 3] - aug_data_ends$X2

    likelihood_loss_pos <- sum(diff_1^2 + diff_2^2) / (2 * (N_advects) * pos_sd^2)

    # Fast Prior calculation using crossprod
    d_beta1 <- beta_1 - start_beta_1
    d_beta2 <- beta_2 - start_beta_2

    prior_loss <- crossprod(d_beta1, prior_beta_sigma %*% d_beta1) +
      crossprod(d_beta2, prior_beta_sigma %*% d_beta2)

    return(as.numeric(likelihood_loss_pos + prior_loss))

  }

  # Run the optimizer
  opt_result <- optim(
    par = start_beta,             # Starting values
    fn = rMAP_loss,           # Your loss function
    method = "BFGS",          # Unconstrained quasi-Newton
    control = list(
      maxit = 1000,           # Increase max iterations (default is 100)
      trace = 1,              # Prints progress to the console
      REPORT = 10             # How often to print progress
    )
  )

  post_beta_1 = opt_result$par[1:(M^2)]
  post_beta_2 = opt_result$par[1:(M^2)+(M^2)]
  post_beta_mat = cbind(post_beta_1, post_beta_2)

  post_traj_grid = exp(evaluate2DCosine(beta_mat = post_beta_mat, pos_mat = traj_eval_grid, border = border))

  #Eval Optimal Traj

  M_sq <- M^2

  # Pre-calculate time step vectors for all i
  t_steps <- (aug_data_ends$t - aug_data_starts$t) / N_prop_steps
  sqrt_t_steps <- sqrt(t_steps)

  # Initialize the starting positions for all (N-1) segments at once
  # curPos_mat is an (N-1) x 3 matrix
  curPos_mat <- aug_data_starts

  # Inner loop: Propagate all segments simultaneously
  for(j in 1:N_prop_steps) {

    # Evaluate drift for all (N-1) positions simultaneously
    # Ensure this function is updated to return an (N-1) x 2 matrix
    cur_drift_mat <- TrajWeightedBaseVectorFields_2D_Cosine(
      pos_t_mat = curPos_mat,
      beta_mat = post_beta_mat,
      baseVectorFields_Vec = baseVectorFields_Vec,
      border = border
    )

    # Extract diffusion for all N-1 segments for step j
    cur_diffusion_vec_1 <- rand_vel_1[j,]
    cur_diffusion_vec_2 <- rand_vel_2[j,]

    # Update Time (Column 1)
    curPos_mat[, 1] <- curPos_mat[, 1] + t_steps

    # Update X and Y positions (Columns 2 and 3)
    curPos_mat[, 2] <- curPos_mat[, 2] + cur_drift_mat[, 1] * t_steps + cur_diffusion_vec_1 * sqrt_t_steps
    curPos_mat[, 3] <- curPos_mat[, 3] + cur_drift_mat[, 2] * t_steps + cur_diffusion_vec_2 * sqrt_t_steps
  }

  # Fast Likelihood calculation using sum of squares instead of matrix multiplication
  diff_1 <- curPos_mat[, 2] - aug_data_ends$X1
  diff_2 <- curPos_mat[, 3] - aug_data_ends$X2

  post_likelihood_loss_pos <- sum(diff_1^2 + diff_2^2) / (2 * (N_advects) * pos_sd^2)

  # Fast Prior calculation using crossprod
  d_beta1 <- post_beta_1 - start_beta_1
  d_beta2 <- post_beta_2 - start_beta_2

  post_prior_loss <- crossprod(d_beta1, prior_beta_sigma %*% d_beta1) +
    crossprod(d_beta2, prior_beta_sigma %*% d_beta2)

  list(post_beta_mat, post_likelihood_loss_pos, post_prior_loss, post_traj_grid)

}

traj_eval_grid = expand.grid(seq(-2,2, length.out = 100), seq(-2,2, length.out = 100))


test_sample = sample_rMAP_trajectories(sim_data = sim_data_list, pos_sd = 0.001, vel_sd = 0.1, M = 2, prior_k = 0.35, prior_l = 1, baseVectorFields_Vec = baseVectorFields_Vec, border = c(-3,-3,3,3), N_prop_steps = 10, traj_eval_grid = traj_eval_grid)

post_traj_grid_plot = cbind(traj_eval_grid, test_sample[[4]])

test_sample[[4]]

ggplot() + geom_point(data = post_traj_grid_plot, aes(x=  Var1, y = Var2, color = abs(log(post_beta_1)) < 0.01), size = 5) + geom_point(data = sim_data, aes(x = X1, y = X2))
ggplot() + geom_point(data = post_traj_grid_plot, aes(x=  Var1, y = Var2, color = abs(log(post_beta_2)) < 0.01), size = 5) + geom_point(data = sim_data, aes(x = X1, y = X2))

#running test

M=2

beta_1_post_mat = matrix(nrow = 100, ncol = M^2)
beta_2_post_mat = matrix(nrow = 100, ncol = M^2)

post_traj_eval_1_list = list()
post_traj_eval_2_list = list()

post_like_pos = rep(0,100)
post_like_prior = rep(0,100)

svMisc::progress(0,100)

for(i in 1:100){

  cur_draw = sample_rMAP_trajectories(sim_data = sim_data_list, pos_sd = sqrt(sqrt(var(post_like_full)/2.83))*0.001, vel_sd = 0.1, M = 2, prior_k = 0.35, prior_l = 1, baseVectorFields_Vec = baseVectorFields_Vec, border = c(-2,-2,2,2), N_prop_steps = 5, traj_eval_grid = traj_eval_grid)

  #Metropolis

  if(i > 1){

    prev_LL =

  }

  beta_1_post_mat_V2[i,] = cur_draw[[1]][,1]
  beta_2_post_mat_V2[i,] = cur_draw[[1]][,2]

  post_like_pos_V2[i] = cur_draw[[2]]
  post_like_prior_V2[i] = cur_draw[[3]]

  post_traj_eval_1_list_V2[[i]] = cur_draw[[4]][,1]
  post_traj_eval_2_list_V2[[i]] = cur_draw[[4]][,2]


  svMisc::progress(i,100)

}

plot(post_like_pos_V2 + post_like_prior_V2)

post_like_full = post_like_pos_V2 + post_like_prior_V2

acc_prob = rep(1,100)

for(i in 2:100){

  cur_like = post_like_full[i]
  prev_like = post_like_full[i-1]

  acc_prob[i] = min(1, exp(cur_like - prev_like))

}

var(post_like_full)

mean(acc_prob)

hist(beta_1_post_mat_V2[,4])
quantile(beta_1_post_mat_V2[,4], probs = c(0.025, 0.975))


test_sample






exp(0.05)

sqrt(sqrt(var(post_like_full)/2.83))



post_draws_traj_eval_1_V2 = do.call(cbind, post_traj_eval_1_list_V2)
post_draws_traj_eval_2_V2 = do.call(cbind, post_traj_eval_2_list_V2)

CI_traj_1_V2 = t(apply(post_draws_traj_eval_1_V2, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975)))
CI_traj_2_V2 = t(apply(post_draws_traj_eval_2_V2, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975)))

mean(CI_traj_1_V2[,1] <= 1 & CI_traj_1_V2[,2] >= 1)

mean(CI_traj_2_V2[,1] <= 1 & CI_traj_2_V2[,2] >= 1)





