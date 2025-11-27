library(plotly)
M = 20
k=0.5
l = 0.05

coefs = matrix(rnorm((M+1)*(M+1)), nrow = M+1, ncol = M+1)

curPos = c(.6,-0.4)

evaluateHSGP = function(z, k, l, M, curPos){

  omega = (0:M)*pi

  spec_den = sqrt(2*pi)*l*exp(-0.5*l^2*omega^2)
  beta = k*diag(sqrt(spec_den)) %*% z %*% diag(sqrt(spec_den))

  phi_x = c(1, sqrt(2)*cos(omega[1:M+1]*curPos[1]))
  phi_y = c(1, sqrt(2)*cos(omega[1:M+1]*curPos[2]))

  as.numeric(phi_x %*% beta %*% phi_y)

}

invlogit <- function(x) {
  1 / (1 + exp(-x))
}
grid_res = 100
plotHSGP = function(grid_res, z_log, k_log, l_log, z_logit, k_logit, l_logit, logit_prior_mean, M, color_limits){

  plot_grid = expand.grid(seq(0,1,length.out = grid_res+1), seq(0,1,length.out = grid_res+1))

  x_seq = y_seq = seq(0,1,length.out = grid_res+1)

  values_coef = matrix(exp(apply(plot_grid, MARGIN = 1, FUN = evaluateHSGP, z = z_log, k=k_log, l=l_log, M=M)), byrow = F, nrow = grid_res+1)
  values_weight = matrix(invlogit(apply(plot_grid, MARGIN = 1, FUN = evaluateHSGP, z = z_logit, k=k_logit, l=l_logit, M=M) + 10), byrow = F, nrow = grid_res+1)

  values_traj = values_coef*values_weight

  coef_surface = plot_ly(x = x_seq, y = y_seq, z = values_coef, type = "surface", color_scale = "Viridis")
  weight_surface = plot_ly(x = x_seq, y = y_seq, z = values_weight, type = "surface")

  traj_surface = plot_ly(x = x_seq, y = y_seq, z = values_traj, type = "surface", color_scale = "Viridis", cmin = color_limits[1], cmax = color_limits[2])

  list("Coefficient" = coef_surface, "Weight" = weight_surface, "Trajectory" = traj_surface)

}

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
