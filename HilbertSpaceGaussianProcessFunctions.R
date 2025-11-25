library(plotly)
M = 20
k=0.5
l = 0.05

coefs = matrix(rnorm((M+1)*(M+1)), nrow = M+1, ncol = M+1)

curPos = c(.6,-0.4)

evaluateHSGP = function(coefs, prior_mean, k, l, M, curPos){

  omega = (0:M)*pi

  spec_den = (k^2)*sqrt(2*pi)*l*exp(-0.5*l^2*omega^2)

  beta = diag(sqrt(spec_den)) %*% coefs %*% diag(sqrt(spec_den))

  phi_x = c(1, sqrt(2)*cos(omega[1:M+1]*curPos[1]))
  phi_y = c(1, sqrt(2)*cos(omega[1:M+1]*curPos[2]))

  as.numeric(phi_x %*% beta %*% phi_y) + prior_mean

}

plotHSGP = function(grid_res, coefs, prior_mean, k, l, M){

  plot_grid = expand.grid(seq(0,1,length.out = grid_res+1), seq(0,1,length.out = grid_res+1))

  x_seq = y_seq = seq(0,1,length.out = grid_res+1)

  values = matrix(apply(plot_grid, MARGIN = 1, FUN = evaluateHSGP, coefs = coefs, prior_mean = prior_mean, k=k, l=l, M=M), byrow = F, nrow = grid_res+1)

  plot_ly(x = x_seq, y = y_seq, z = values, type = "surface")

}

plotHSGP(100, coefs, prior_mean, k, l = 0.01, M)
