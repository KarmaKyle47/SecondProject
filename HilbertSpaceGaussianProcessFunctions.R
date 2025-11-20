
M = 20
k=0.5
l = 0.3

coefs = matrix(c(rnorm(1,1,1), rnorm(21*21-1)), nrow = 21, ncol = 21)

curPos = c(.6,-0.4)

evaluateHSGP = function(coefs, k, l, M, curPos){

  omega = (0:M)*pi

  spec_den = (k^2)*sqrt(2*pi)*l*exp(-0.5*l^2*freqs^2)

  beta = diag(sqrt(spec_den)) %*% coefs %*% diag(sqrt(spec_den))

  phi_x = c(1, sqrt(2)*cos(omega[1:M+1]*curPos[1]))
  phi_y = c(1, sqrt(2)*cos(omega[1:M+1]*curPos[2]))

  phi_x %*% beta %*% phi_y

  sum(outer(phi_x, phi_y) * beta)


}
