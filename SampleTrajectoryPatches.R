library(MASS)
library(plotly)

evaluateCubicPatchValue = function(coef, border, curPos){

  x_0 = border[1]
  y_0 = border[2]
  x = curPos[1]
  y = curPos[2]

  polyTerms = c(1, (y-y_0), (x-x_0), (x-x_0)*(y-y_0),
                (x-x_0)^2, (x-x_0)^2*(y-y_0), (x-x_0)^3, (x-x_0)^3*(y-y_0),
                (y-y_0)^2, (y-y_0)^3, (x-x_0)*(y-y_0)^2, (x-x_0)*(y-y_0)^3,
                (x-x_0)^2*(y-y_0)^2, (x-x_0)^2*(y-y_0)^3, (x-x_0)^3*(y-y_0)^2, (x-x_0)^3*(y-y_0)^3)

  sum(coef * polyTerms)

}
evaluateCubicPatchParX = function(coef, border, curPos){

  x_0 = border[1]
  y_0 = border[2]
  x = curPos[1]
  y = curPos[2]

  polyTerms = c(0, 0, 1, (y-y_0),
                2*(x-x_0), 2*(x-x_0)*(y-y_0), 3*(x-x_0)^2, 3*(x-x_0)^2*(y-y_0),
                0, 0, (y-y_0)^2, (y-y_0)^3,
                2*(x-x_0)*(y-y_0)^2, 2*(x-x_0)*(y-y_0)^3, 3*(x-x_0)^2*(y-y_0)^2, 3*(x-x_0)^2*(y-y_0)^3)

  sum(coef * polyTerms)
}
evaluateCubicPatchParY = function(coef, border, curPos){

  x_0 = border[1]
  y_0 = border[2]
  x = curPos[1]
  y = curPos[2]

  polyTerms = c(0, 1, 0, (x-x_0),
                0, (x-x_0)^2, 0, (x-x_0)^3,
                2*(y-y_0), 3*(y-y_0)^2, 2*(x-x_0)*(y-y_0), 3*(x-x_0)*(y-y_0)^2,
                2*(x-x_0)^2*(y-y_0), 3*(x-x_0)^2*(y-y_0)^2, 2*(x-x_0)^3*(y-y_0), 3*(x-x_0)^3*(y-y_0)^2)

  sum(coef * polyTerms)


}
evaluateCubicPatchParXY = function(coef, border, curPos){

  x_0 = border[1]
  y_0 = border[2]
  x = curPos[1]
  y = curPos[2]

  polyTerms = c(0, 0, 0, 1,
                0, 2*(x-x_0), 0, 3*(x-x_0)^2,
                0, 0, 2*(y-y_0), 3*(y-y_0)^2,
                4*(x-x_0)*(y-y_0), 6*(x-x_0)*(y-y_0)^2, 6*(x-x_0)^2*(y-y_0), 9*(x-x_0)^2*(y-y_0)^2)

  sum(coef * polyTerms)
}

samplePatch_Unrestricted = function(border, k, prior_mean = 1){

  dx = border[3] - border[1]
  dy = border[4] - border[2]

  A = matrix(c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
               1,0,dx,0,dx^2,0,dx^3,0,0,0,0,0,0,0,0,0,
               1,dy,0,0,0,0,0,0,dy^2,dy^3,0,0,0,0,0,0,
               1,dy,dx,dx*dy,dx^2,dx^2*dy,dx^3,dx^3*dy,dy^2,dy^3,dx*dy^2,dx*dy^3,dx^2*dy^2,dx^2*dy^3,dx^3*dy^2,dx^3*dy^3,
               0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
               0,0,1,0,2*dx,0,3*dx^2,0,0,0,0,0,0,0,0,0,
               0,0,1,dy,0,0,0,0,0,0,dy^2,dy^3,0,0,0,0,
               0,0,1,dy,2*dx,2*dx*dy,3*dx^2,3*dx^2*dy,0,0,dy^2,dy^3,2*dx*dy^2,2*dx*dy^3,3*dx^2*dy^2,3*dx^2*dy^3,
               0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
               0,1,0,dx,0,dx^2,0,dx^3,0,0,0,0,0,0,0,0,
               0,1,0,0,0,0,0,0,2*dy,3*dy^2,0,0,0,0,0,0,
               0,1,0,dx,0,dx^2,0,dx^3,2*dy,3*dy^2,2*dx*dy,3*dx*dy^2,2*dx^2*dy,3*dx^2*dy^2,2*dx^3*dy,3*dx^3*dy^2,
               0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
               0,0,0,1,0,2*dx,0,3*dx^2,0,0,0,0,0,0,0,0,
               0,0,0,1,0,0,0,0,0,0,2*dy,3*dy^2,0,0,0,0,
               0,0,0,1,0,2*dx,0,3*dx^2,0,0,2*dy,3*dy^2,4*dx*dy,6*dx*dy^2,6*dx^2*dy,9*dx^2*dy^2), nrow = 16, ncol = 16, byrow = T)

  A_inv = solve(A)

  mu = A_inv %*% c(rep(prior_mean,4), rep(0, 12))

  sigma = k^2*(A_inv %*% t(A_inv))

  coef = mvrnorm(mu = mu, Sigma = sigma)

  coef

}

samplePatch_Known_Bottom = function(border, k, BL_coefs, BL_border, BR_coefs, BR_border, prior_mean = 1){

  dx = border[3] - border[1]
  dy = border[4] - border[2]

  A = matrix(c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
               1,0,dx,0,dx^2,0,dx^3,0,0,0,0,0,0,0,0,0,
               1,dy,0,0,0,0,0,0,dy^2,dy^3,0,0,0,0,0,0,
               1,dy,dx,dx*dy,dx^2,dx^2*dy,dx^3,dx^3*dy,dy^2,dy^3,dx*dy^2,dx*dy^3,dx^2*dy^2,dx^2*dy^3,dx^3*dy^2,dx^3*dy^3,
               0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
               0,0,1,0,2*dx,0,3*dx^2,0,0,0,0,0,0,0,0,0,
               0,0,1,dy,0,0,0,0,0,0,dy^2,dy^3,0,0,0,0,
               0,0,1,dy,2*dx,2*dx*dy,3*dx^2,3*dx^2*dy,0,0,dy^2,dy^3,2*dx*dy^2,2*dx*dy^3,3*dx^2*dy^2,3*dx^2*dy^3,
               0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
               0,1,0,dx,0,dx^2,0,dx^3,0,0,0,0,0,0,0,0,
               0,1,0,0,0,0,0,0,2*dy,3*dy^2,0,0,0,0,0,0,
               0,1,0,dx,0,dx^2,0,dx^3,2*dy,3*dy^2,2*dx*dy,3*dx*dy^2,2*dx^2*dy,3*dx^2*dy^2,2*dx^3*dy,3*dx^3*dy^2,
               0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
               0,0,0,1,0,2*dx,0,3*dx^2,0,0,0,0,0,0,0,0,
               0,0,0,1,0,0,0,0,0,0,2*dy,3*dy^2,0,0,0,0,
               0,0,0,1,0,2*dx,0,3*dx^2,0,0,2*dy,3*dy^2,4*dx*dy,6*dx*dy^2,6*dx^2*dy,9*dx^2*dy^2), nrow = 16, ncol = 16, byrow = T)

  c00 = evaluateCubicPatchValue(coef = BL_coefs, border = BL_border, curPos = c(border[1], border[2]))
  c01 = evaluateCubicPatchParY(coef = BL_coefs, border = BL_border, curPos = c(border[1], border[2]))
  c10 = evaluateCubicPatchParX(coef = BL_coefs, border = BL_border, curPos = c(border[1], border[2]))
  c11 = evaluateCubicPatchParXY(coef = BL_coefs, border = BL_border, curPos = c(border[1], border[2]))

  value_BR = evaluateCubicPatchValue(coef = BR_coefs, border = BR_border, curPos = c(border[3], border[2]))
  par_X_BR = evaluateCubicPatchParX(coef = BR_coefs, border = BR_border, curPos = c(border[3], border[2]))
  par_Y_BR = evaluateCubicPatchParY(coef = BR_coefs, border = BR_border, curPos = c(border[3], border[2]))
  par_XY_BR = evaluateCubicPatchParXY(coef = BR_coefs, border = BR_border, curPos = c(border[3], border[2]))

  c20 = (-1*par_X_BR*dx + 3*value_BR - 3*c00 - 2*c10*dx)/(dx^2)
  c21 = (-1*par_XY_BR*dx + 3*par_Y_BR - 3*c01 - 2*c11*dx)/(dx^2)
  c30 = (par_X_BR*dx - 2*value_BR + 2*c00 + c10*dx)/(dx^3)
  c31 = (par_XY_BR*dx - 2*par_Y_BR + 2*c01 + c11*dx)/(dx^3)


  known_coefs = c(c00, c01, c10, c11, c20, c21, c30, c31)

  mu_trans = c(rep(prior_mean,2), rep(0,6)) - (A[c(3,4,7,8,11,12,15,16), 1:8] %*% known_coefs)

  A_inv = solve(A[c(3,4,7,8,11,12,15,16), 9:16])

  mu = A_inv %*% mu_trans
  sigma = k^2*(A_inv %*% t(A_inv))

  coef = mvrnorm(mu = mu, Sigma = sigma)

  c(known_coefs, coef)

}

samplePatch_Known_Left = function(border, k, BL_coefs, BL_border, TL_coefs, TL_border, prior_mean = 1){

  dx = border[3] - border[1]
  dy = border[4] - border[2]

  A = matrix(c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
               1,0,dx,0,dx^2,0,dx^3,0,0,0,0,0,0,0,0,0,
               1,dy,0,0,0,0,0,0,dy^2,dy^3,0,0,0,0,0,0,
               1,dy,dx,dx*dy,dx^2,dx^2*dy,dx^3,dx^3*dy,dy^2,dy^3,dx*dy^2,dx*dy^3,dx^2*dy^2,dx^2*dy^3,dx^3*dy^2,dx^3*dy^3,
               0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
               0,0,1,0,2*dx,0,3*dx^2,0,0,0,0,0,0,0,0,0,
               0,0,1,dy,0,0,0,0,0,0,dy^2,dy^3,0,0,0,0,
               0,0,1,dy,2*dx,2*dx*dy,3*dx^2,3*dx^2*dy,0,0,dy^2,dy^3,2*dx*dy^2,2*dx*dy^3,3*dx^2*dy^2,3*dx^2*dy^3,
               0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
               0,1,0,dx,0,dx^2,0,dx^3,0,0,0,0,0,0,0,0,
               0,1,0,0,0,0,0,0,2*dy,3*dy^2,0,0,0,0,0,0,
               0,1,0,dx,0,dx^2,0,dx^3,2*dy,3*dy^2,2*dx*dy,3*dx*dy^2,2*dx^2*dy,3*dx^2*dy^2,2*dx^3*dy,3*dx^3*dy^2,
               0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
               0,0,0,1,0,2*dx,0,3*dx^2,0,0,0,0,0,0,0,0,
               0,0,0,1,0,0,0,0,0,0,2*dy,3*dy^2,0,0,0,0,
               0,0,0,1,0,2*dx,0,3*dx^2,0,0,2*dy,3*dy^2,4*dx*dy,6*dx*dy^2,6*dx^2*dy,9*dx^2*dy^2), nrow = 16, ncol = 16, byrow = T)

  c00 = evaluateCubicPatchValue(coef = BL_coefs, border = BL_border, curPos = c(border[1], border[2]))
  c01 = evaluateCubicPatchParY(coef = BL_coefs, border = BL_border, curPos = c(border[1], border[2]))
  c10 = evaluateCubicPatchParX(coef = BL_coefs, border = BL_border, curPos = c(border[1], border[2]))
  c11 = evaluateCubicPatchParXY(coef = BL_coefs, border = BL_border, curPos = c(border[1], border[2]))

  value_TL = evaluateCubicPatchValue(coef = TL_coefs, border = TL_border, curPos = c(border[1], border[4]))
  par_X_TL = evaluateCubicPatchParX(coef = TL_coefs, border = TL_border, curPos = c(border[1], border[4]))
  par_Y_TL = evaluateCubicPatchParY(coef = TL_coefs, border = TL_border, curPos = c(border[1], border[4]))
  par_XY_TL = evaluateCubicPatchParXY(coef = TL_coefs, border = TL_border, curPos = c(border[1], border[4]))

  c02 = (-1*par_Y_TL*dy + 3*value_TL - 3*c00 - 2*c01*dy)/(dy^2)
  c12 = (-1*par_XY_TL*dy + 3*par_X_TL - 3*c10 - 2*c11*dy)/(dy^2)
  c03 = (par_Y_TL*dy - 2*value_TL + 2*c00 + c01*dy)/(dy^3)
  c13 = (par_XY_TL*dy - 2*par_X_TL + 2*c10 + c11*dy)/(dy^3)


  known_coefs = c(c00, c01, c10, c11, c02, c03, c12, c13)

  mu_trans = c(rep(prior_mean,2), rep(0,6)) - (A[c(2,4,6,8,10,12,14,16), c(1:4,9:12)] %*% known_coefs)

  A_inv = solve(A[c(2,4,6,8,10,12,14,16), c(5:8,13:16)])

  mu = A_inv %*% mu_trans
  sigma = k^2*(A_inv %*% t(A_inv))

  coef = mvrnorm(mu = mu, Sigma = sigma)

  c(known_coefs[1:4], coef[1:4], known_coefs[5:8], coef[5:8])


}

samplePatch_Known_Bottom_and_Left = function(border, k, BL_coefs, BL_border, BR_coefs, BR_border, TL_coefs, TL_border, prior_mean = 1){

  dx = border[3] - border[1]
  dy = border[4] - border[2]

  A = matrix(c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
               1,0,dx,0,dx^2,0,dx^3,0,0,0,0,0,0,0,0,0,
               1,dy,0,0,0,0,0,0,dy^2,dy^3,0,0,0,0,0,0,
               1,dy,dx,dx*dy,dx^2,dx^2*dy,dx^3,dx^3*dy,dy^2,dy^3,dx*dy^2,dx*dy^3,dx^2*dy^2,dx^2*dy^3,dx^3*dy^2,dx^3*dy^3,
               0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
               0,0,1,0,2*dx,0,3*dx^2,0,0,0,0,0,0,0,0,0,
               0,0,1,dy,0,0,0,0,0,0,dy^2,dy^3,0,0,0,0,
               0,0,1,dy,2*dx,2*dx*dy,3*dx^2,3*dx^2*dy,0,0,dy^2,dy^3,2*dx*dy^2,2*dx*dy^3,3*dx^2*dy^2,3*dx^2*dy^3,
               0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
               0,1,0,dx,0,dx^2,0,dx^3,0,0,0,0,0,0,0,0,
               0,1,0,0,0,0,0,0,2*dy,3*dy^2,0,0,0,0,0,0,
               0,1,0,dx,0,dx^2,0,dx^3,2*dy,3*dy^2,2*dx*dy,3*dx*dy^2,2*dx^2*dy,3*dx^2*dy^2,2*dx^3*dy,3*dx^3*dy^2,
               0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
               0,0,0,1,0,2*dx,0,3*dx^2,0,0,0,0,0,0,0,0,
               0,0,0,1,0,0,0,0,0,0,2*dy,3*dy^2,0,0,0,0,
               0,0,0,1,0,2*dx,0,3*dx^2,0,0,2*dy,3*dy^2,4*dx*dy,6*dx*dy^2,6*dx^2*dy,9*dx^2*dy^2), nrow = 16, ncol = 16, byrow = T)

  c00 = evaluateCubicPatchValue(coef = BL_coefs, border = BL_border, curPos = c(border[1], border[2]))
  c01 = evaluateCubicPatchParY(coef = BL_coefs, border = BL_border, curPos = c(border[1], border[2]))
  c10 = evaluateCubicPatchParX(coef = BL_coefs, border = BL_border, curPos = c(border[1], border[2]))
  c11 = evaluateCubicPatchParXY(coef = BL_coefs, border = BL_border, curPos = c(border[1], border[2]))

  value_BR = evaluateCubicPatchValue(coef = BR_coefs, border = BR_border, curPos = c(border[3], border[2]))
  par_X_BR = evaluateCubicPatchParX(coef = BR_coefs, border = BR_border, curPos = c(border[3], border[2]))
  par_Y_BR = evaluateCubicPatchParY(coef = BR_coefs, border = BR_border, curPos = c(border[3], border[2]))
  par_XY_BR = evaluateCubicPatchParXY(coef = BR_coefs, border = BR_border, curPos = c(border[3], border[2]))

  c20 = (-1*par_X_BR*dx + 3*value_BR - 3*c00 - 2*c10*dx)/(dx^2)
  c21 = (-1*par_XY_BR*dx + 3*par_Y_BR - 3*c01 - 2*c11*dx)/(dx^2)
  c30 = (par_X_BR*dx - 2*value_BR + 2*c00 + c10*dx)/(dx^3)
  c31 = (par_XY_BR*dx - 2*par_Y_BR + 2*c01 + c11*dx)/(dx^3)

  value_TL = evaluateCubicPatchValue(coef = TL_coefs, border = TL_border, curPos = c(border[1], border[4]))
  par_X_TL = evaluateCubicPatchParX(coef = TL_coefs, border = TL_border, curPos = c(border[1], border[4]))
  par_Y_TL = evaluateCubicPatchParY(coef = TL_coefs, border = TL_border, curPos = c(border[1], border[4]))
  par_XY_TL = evaluateCubicPatchParXY(coef = TL_coefs, border = TL_border, curPos = c(border[1], border[4]))

  c02 = (-1*par_Y_TL*dy + 3*value_TL - 3*c00 - 2*c01*dy)/(dy^2)
  c12 = (-1*par_XY_TL*dy + 3*par_X_TL - 3*c10 - 2*c11*dy)/(dy^2)
  c03 = (par_Y_TL*dy - 2*value_TL + 2*c00 + c01*dy)/(dy^3)
  c13 = (par_XY_TL*dy - 2*par_X_TL + 2*c10 + c11*dy)/(dy^3)

  known_coefs = c(c00, c01, c10, c11, c20, c21, c30, c31, c02, c03, c12, c13)

  mu_trans = c(prior_mean,0,0,0) - (A[c(4,8,12,16), 1:12] %*% known_coefs)

  A_inv = solve(A[c(4,8,12,16), 13:16])

  mu = A_inv %*% mu_trans
  sigma = k^2*(A_inv %*% t(A_inv))

  coef = mvrnorm(mu = mu, Sigma = sigma)

  c(known_coefs, coef)


}

samplePatch_Known_Corners = function(border, k, BL_coefs, BL_border, BR_coefs, BR_border, TL_coefs, TL_border, TR_coefs, TR_border){

  dx = border[3] - border[1]
  dy = border[4] - border[2]

  c00 = evaluateCubicPatchValue(coef = BL_coefs, border = BL_border, curPos = c(border[1], border[2]))
  c01 = evaluateCubicPatchParY(coef = BL_coefs, border = BL_border, curPos = c(border[1], border[2]))
  c10 = evaluateCubicPatchParX(coef = BL_coefs, border = BL_border, curPos = c(border[1], border[2]))
  c11 = evaluateCubicPatchParXY(coef = BL_coefs, border = BL_border, curPos = c(border[1], border[2]))

  value_BR = evaluateCubicPatchValue(coef = BR_coefs, border = BR_border, curPos = c(border[3], border[2]))
  par_X_BR = evaluateCubicPatchParX(coef = BR_coefs, border = BR_border, curPos = c(border[3], border[2]))
  par_Y_BR = evaluateCubicPatchParY(coef = BR_coefs, border = BR_border, curPos = c(border[3], border[2]))
  par_XY_BR = evaluateCubicPatchParXY(coef = BR_coefs, border = BR_border, curPos = c(border[3], border[2]))

  c20 = (-1*par_X_BR*dx + 3*value_BR - 3*c00 - 2*c10*dx)/(dx^2)
  c21 = (-1*par_XY_BR*dx + 3*par_Y_BR - 3*c01 - 2*c11*dx)/(dx^2)
  c30 = (par_X_BR*dx - 2*value_BR + 2*c00 + c10*dx)/(dx^3)
  c31 = (par_XY_BR*dx - 2*par_Y_BR + 2*c01 + c11*dx)/(dx^3)

  value_TL = evaluateCubicPatchValue(coef = TL_coefs, border = TL_border, curPos = c(border[1], border[4]))
  par_X_TL = evaluateCubicPatchParX(coef = TL_coefs, border = TL_border, curPos = c(border[1], border[4]))
  par_Y_TL = evaluateCubicPatchParY(coef = TL_coefs, border = TL_border, curPos = c(border[1], border[4]))
  par_XY_TL = evaluateCubicPatchParXY(coef = TL_coefs, border = TL_border, curPos = c(border[1], border[4]))

  c02 = (-1*par_Y_TL*dy + 3*value_TL - 3*c00 - 2*c01*dy)/(dy^2)
  c12 = (-1*par_XY_TL*dy + 3*par_X_TL - 3*c10 - 2*c11*dy)/(dy^2)
  c03 = (par_Y_TL*dy - 2*value_TL + 2*c00 + c01*dy)/(dy^3)
  c13 = (par_XY_TL*dy - 2*par_X_TL + 2*c10 + c11*dy)/(dy^3)

  value_TR = evaluateCubicPatchValue(coef = TR_coefs, border = TR_border, curPos = c(border[3], border[4]))
  par_X_TR = evaluateCubicPatchParX(coef = TR_coefs, border = TR_border, curPos = c(border[3], border[4]))
  par_Y_TR = evaluateCubicPatchParY(coef = TR_coefs, border = TR_border, curPos = c(border[3], border[4]))
  par_XY_TR = evaluateCubicPatchParXY(coef = TR_coefs, border = TR_border, curPos = c(border[3], border[4]))

  temp_coefs = c(c00, c01, c10, c11, c20, c21, c30, c31, c02, c03, c12, c13, 0, 0, 0, 0)

  Cz = evaluateCubicPatchValue(coef = temp_coefs, border = border, curPos = c(border[3], border[4]))
  Cx = evaluateCubicPatchParX(coef = temp_coefs, border = border, curPos = c(border[3], border[4]))
  Cy = evaluateCubicPatchParY(coef = temp_coefs, border = border, curPos = c(border[3], border[4]))
  Cxy = evaluateCubicPatchParXY(coef = temp_coefs, border = border, curPos = c(border[3], border[4]))

  c22 = ((par_XY_TR - Cxy)*dx*dy - 3*(par_Y_TR - Cy)*dy - 3*(par_X_TR - Cx)*dx + 9*(value_TR - Cz))/(dx^2*dy^2)
  c23 = (-1*(par_XY_TR - Cxy)*dx*dy + 3*(par_Y_TR - Cy)*dy + 2*(par_X_TR - Cx)*dx - 6*(value_TR - Cz))/(dx^2*dy^3)
  c32 = (-1*(par_XY_TR - Cxy)*dx*dy + 2*(par_Y_TR - Cy)*dy + 3*(par_X_TR - Cx)*dx - 6*(value_TR - Cz))/(dx^3*dy^2)
  c33 = ((par_XY_TR - Cxy)*dx*dy - 2*(par_Y_TR - Cy)*dy - 2*(par_X_TR - Cx)*dx + 4*(value_TR - Cz))/(dx^3*dy^3)


  known_coefs = c(c00, c01, c10, c11, c20, c21, c30, c31, c02, c03, c12, c13, c22, c23, c32, c33)

  known_coefs

}

#CornerVectors -> BL, BR, TL, TR

calculatePatch_KnownDerivates = function(border, CornerValues, CornerParXs, CornerParYs, CornerParXYs){

  dx = border[3] - border[1]
  dy = border[4] - border[2]

  c00 = CornerValues[1]
  c01 = CornerParYs[1]
  c10 = CornerParXs[1]
  c11 = CornerParXYs[1]

  c20 = (-1*CornerParXs[2]*dx + 3*CornerValues[2] - 3*c00 - 2*c10*dx)/(dx^2)
  c21 = (-1*CornerParXYs[2]*dx + 3*CornerParYs[2] - 3*c01 - 2*c11*dx)/(dx^2)
  c30 = (CornerParXs[2]*dx - 2*CornerValues[2] + 2*c00 + c10*dx)/(dx^3)
  c31 = (CornerParXYs[2]*dx - 2*CornerParYs[2] + 2*c01 + c11*dx)/(dx^3)

  c02 = (-1*CornerParYs[3]*dy + 3*CornerValues[3] - 3*c00 - 2*c01*dy)/(dy^2)
  c12 = (-1*CornerParXYs[3]*dy + 3*CornerParXs[3] - 3*c10 - 2*c11*dy)/(dy^2)
  c03 = (CornerParYs[3]*dy - 2*CornerValues[3] + 2*c00 + c01*dy)/(dy^3)
  c13 = (CornerParXYs[3]*dy - 2*CornerParXs[3] + 2*c10 + c11*dy)/(dy^3)

  temp_coefs = c(c00, c01, c10, c11, c20, c21, c30, c31, c02, c03, c12, c13, 0, 0, 0, 0)

  Cz = evaluateCubicPatchValue(coef = temp_coefs, border = border, curPos = c(border[3], border[4]))
  Cx = evaluateCubicPatchParX(coef = temp_coefs, border = border, curPos = c(border[3], border[4]))
  Cy = evaluateCubicPatchParY(coef = temp_coefs, border = border, curPos = c(border[3], border[4]))
  Cxy = evaluateCubicPatchParXY(coef = temp_coefs, border = border, curPos = c(border[3], border[4]))

  c22 = ((CornerParXYs[4] - Cxy)*dx*dy - 3*(CornerParYs[4] - Cy)*dy - 3*(CornerParXs[4] - Cx)*dx + 9*(CornerValues[4] - Cz))/(dx^2*dy^2)
  c23 = (-1*(CornerParXYs[4] - Cxy)*dx*dy + 3*(CornerParYs[4] - Cy)*dy + 2*(CornerParXs[4] - Cx)*dx - 6*(CornerValues[4] - Cz))/(dx^2*dy^3)
  c32 = (-1*(CornerParXYs[4] - Cxy)*dx*dy + 2*(CornerParYs[4] - Cy)*dy + 3*(CornerParXs[4] - Cx)*dx - 6*(CornerValues[4] - Cz))/(dx^3*dy^2)
  c33 = ((CornerParXYs[4] - Cxy)*dx*dy - 2*(CornerParYs[4] - Cy)*dy - 2*(CornerParXs[4] - Cx)*dx + 4*(CornerValues[4] - Cz))/(dx^3*dy^3)


  known_coefs = c(c00, c01, c10, c11, c20, c21, c30, c31, c02, c03, c12, c13, c22, c23, c32, c33)

  known_coefs

}

#GridVectors -> BL corner of entire grid moving L to R then bottom to top

calculateSurface_KnownCorners = function(boundaries, GridValues, GridParXs, GridParYs, GridParXYs){

  x_grid = sort(unique(c(boundaries$L1, boundaries$U1)))
  y_grid = sort(unique(c(boundaries$L2, boundaries$U2)))

  x_grid_len = length(x_grid)
  y_grid_len = length(y_grid)

  Coefs = matrix(nrow = nrow(boundaries), ncol = 16)

  for(i in 1:nrow(boundaries)){

    cur_border = as.numeric(boundaries[i,1:4])

    x_index = which(x_grid == cur_border[1])
    y_index = which(y_grid == cur_border[2])

    cur_corner_indices = c(0:1 + (y_index-1)*x_grid_len + x_index, 0:1 + (y_index)*x_grid_len + x_index)

    cur_CornerValues = GridValues[cur_corner_indices]
    cur_CornerParXs = GridParXs[cur_corner_indices]
    cur_CornerParYs = GridParYs[cur_corner_indices]
    cur_CornerParXYs = GridParXYs[cur_corner_indices]

    Coefs[i,] = calculatePatch_KnownDerivates(cur_border, cur_CornerValues, cur_CornerParXs, cur_CornerParYs, cur_CornerParXYs)

  }

  Coefs

}

samplePatch_FullTree = function(tree, k, prior_mean = 1){

  treeBoundaries = orderBoundaries_GeminiCleaned(tree, F)
  label_order = treeBoundaries$order

  sampledCoefs = matrix(nrow = nrow(treeBoundaries), ncol = 16)

  for(i in 1:nrow(sampledCoefs)){

    cur_index = which(label_order == i)
    cur_border = as.numeric(treeBoundaries[cur_index, 1:4])

    cur_neighbors = find_neighbors(treeBoundaries, cur_index, T)

    pot_BL = unname(label_order[c(cur_neighbors$Down[,1], cur_neighbors$Left[,1])])
    pot_BL = pot_BL[pot_BL < i]

    pot_BR = unname(label_order[c(cur_neighbors$Down[,1], cur_neighbors$Right[,1])])
    pot_BR = pot_BR[pot_BR < i]

    pot_TL = unname(label_order[c(cur_neighbors$Left[,1], cur_neighbors$Up[,1])])
    pot_TL = pot_TL[pot_TL < i]

    pot_TR = unname(label_order[c(cur_neighbors$Right[,1], cur_neighbors$Up[,1])])
    pot_TR = pot_TR[pot_TR < i]

    num_corners_defined = (length(pot_BL) > 0) + (length(pot_BR) > 0) + (length(pot_TL) > 0) + (length(pot_TR) > 0)

    if(num_corners_defined == 0){

      curCoefs = samplePatch_Unrestricted(border = cur_border, k = k, prior_mean = prior_mean)
      sampledCoefs[cur_index,] = curCoefs

    } else if(num_corners_defined == 2){

      if(length(pot_BR) > 0){

        BL_order = min(c(suppressWarnings(min(label_order[cur_neighbors$Down[,1]])),suppressWarnings(min(label_order[cur_neighbors$Left[,1]]))))
        BL_index = which(label_order == BL_order)
        BL_coefs = sampledCoefs[BL_index,]
        BL_border = as.numeric(treeBoundaries[BL_index, 1:4])

        BR_order = min(c(suppressWarnings(max(label_order[cur_neighbors$Down[,1]])),suppressWarnings(min(label_order[cur_neighbors$Right[,1]]))))
        BR_index = which(label_order == BR_order)
        BR_coefs = sampledCoefs[BR_index,]
        BR_border = as.numeric(treeBoundaries[BR_index, 1:4])

        curCoefs = samplePatch_Known_Bottom(cur_border, k = k, BL_coefs, BL_border, BR_coefs, BR_border, prior_mean = prior_mean)
        sampledCoefs[cur_index,] = curCoefs

      } else{

        BL_order = min(c(suppressWarnings(min(label_order[cur_neighbors$Down[,1]])),suppressWarnings(min(label_order[cur_neighbors$Left[,1]]))))
        BL_index = which(label_order == BL_order)
        BL_coefs = sampledCoefs[BL_index,]
        BL_border = as.numeric(treeBoundaries[BL_index, 1:4])

        TL_order = min(c(suppressWarnings(max(label_order[cur_neighbors$Left[,1]])),suppressWarnings(min(label_order[cur_neighbors$Up[,1]]))))
        TL_index = which(label_order == TL_order)
        TL_coefs = sampledCoefs[TL_index,]
        TL_border = as.numeric(treeBoundaries[TL_index, 1:4])

        curCoefs = samplePatch_Known_Left(cur_border, k = k, BL_coefs, BL_border, TL_coefs, TL_border, prior_mean = prior_mean)
        sampledCoefs[cur_index,] = curCoefs

      }

    } else if(num_corners_defined == 3){

      BL_order = min(c(suppressWarnings(min(label_order[cur_neighbors$Down[,1]])),suppressWarnings(min(label_order[cur_neighbors$Left[,1]]))))
      BL_index = which(label_order == BL_order)
      BL_coefs = sampledCoefs[BL_index,]
      BL_border = as.numeric(treeBoundaries[BL_index, 1:4])

      BR_order = min(c(suppressWarnings(max(label_order[cur_neighbors$Down[,1]])),suppressWarnings(min(label_order[cur_neighbors$Right[,1]]))))
      BR_index = which(label_order == BR_order)
      BR_coefs = sampledCoefs[BR_index,]
      BR_border = as.numeric(treeBoundaries[BR_index, 1:4])

      TL_order = min(c(suppressWarnings(max(label_order[cur_neighbors$Left[,1]])),suppressWarnings(min(label_order[cur_neighbors$Up[,1]]))))
      TL_index = which(label_order == TL_order)
      TL_coefs = sampledCoefs[TL_index,]
      TL_border = as.numeric(treeBoundaries[TL_index, 1:4])

      curCoefs = samplePatch_Known_Bottom_and_Left(cur_border, k = k, BL_coefs, BL_border, BR_coefs, BR_border, TL_coefs, TL_border, prior_mean = prior_mean)
      sampledCoefs[cur_index,] = curCoefs


    } else{

      BL_order = min(c(suppressWarnings(min(label_order[cur_neighbors$Down[,1]])),suppressWarnings(min(label_order[cur_neighbors$Left[,1]]))))
      BL_index = which(label_order == BL_order)
      BL_coefs = sampledCoefs[BL_index,]
      BL_border = as.numeric(treeBoundaries[BL_index, 1:4])

      BR_order = min(c(suppressWarnings(max(label_order[cur_neighbors$Down[,1]])),suppressWarnings(min(label_order[cur_neighbors$Right[,1]]))))
      BR_index = which(label_order == BR_order)
      BR_coefs = sampledCoefs[BR_index,]
      BR_border = as.numeric(treeBoundaries[BR_index, 1:4])

      TL_order = min(c(suppressWarnings(max(label_order[cur_neighbors$Left[,1]])),suppressWarnings(min(label_order[cur_neighbors$Up[,1]]))))
      TL_index = which(label_order == TL_order)
      TL_coefs = sampledCoefs[TL_index,]
      TL_border = as.numeric(treeBoundaries[TL_index, 1:4])

      TR_order = min(c(suppressWarnings(max(label_order[cur_neighbors$Right[,1]])),suppressWarnings(max(label_order[cur_neighbors$Up[,1]]))))
      TR_index = which(label_order == TR_order)
      TR_coefs = sampledCoefs[TR_index,]
      TR_border = as.numeric(treeBoundaries[TR_index, 1:4])

      curCoefs = samplePatch_Known_Corners(cur_border, k = k, BL_coefs, BL_border, BR_coefs, BR_border, TL_coefs, TL_border, TR_coefs, TR_border)
      sampledCoefs[cur_index,] = curCoefs

    }


  }

  tree$boundaries = treeBoundaries
  tree$coefs = sampledCoefs

  tree
}

evaluateSampledTreeValue = function(tree, coefs, curPos){

  cur_index = which(curPos[1] >= tree$boundaries$L1 & curPos[1] <= tree$boundaries$U1 & curPos[2] >= tree$boundaries$L2 & curPos[2] <= tree$boundaries$U2)[1]

  cur_border = as.numeric(tree$boundaries[cur_index, 1:4])
  cur_coefs = as.numeric(coefs[cur_index,])

  evaluateCubicPatchValue(cur_coefs, cur_border, curPos)

}

plotCubicPatch3D = function(sampledTree, grid_size = 0.1, z_limit = c(-5,6)){

  border = sampledTree$border

  test_grid = expand.grid(seq(border[1],border[3], grid_size), seq(border[2],border[4],grid_size))
  PatchValues = apply(test_grid, 1, FUN = evaluateSampledTreeValue, tree = sampledTree, coefs = sampledTree$coefs)

  plottingGrid = data.frame(test_grid, PatchValues)
  names(plottingGrid) = c('X','Y','Value')

  Value <- xtabs(Value ~ Y + X, data = plottingGrid)
  X <- as.numeric(rownames(Value))
  Y <- as.numeric(colnames(Value))

  plot = plot_ly(x = X, y = Y, z = Value, type = "surface")

  plot %>% layout(
    title = "3D Surface with Z-Axis Limit",
    scene = list(
      xaxis = list(title = 'X'),
      yaxis = list(title = 'Y'),
      zaxis = list(
        title = 'Value',
        range = z_limit # <-- This sets the Z-axis limit
      )
    )
  )

}

#####

tree = generate_grid_tree(0.1,c(0,0,1,1))
plotTreeGrid(tree)

sampledTree = samplePatch_FullTree(tree, 0.1)

plotCubicPatch3D(sampledTree,grid_size = 0.01, z_limit = c(-5,6))

plotTreeGrid(sampledTree)

init_coef = samplePatch_Unrestricted(border_init, 0.5)

test_grid = expand.grid(seq(0,2, 0.1), seq(0,3,0.1))

PatchValues = apply(test_grid, 1, FUN = evaluateCubicPatchValue, coef = init_coef, border = border_init)

plottingGrid = data.frame(test_grid, PatchValues)
names(plottingGrid) = c('X','Y','Value')

ggplot(plottingGrid, aes(x = X, y = Y, color = Value)) + geom_point(shape = 15, size = 20)

plot_ly(plottingGrid, x = ~X, y = ~Y, z = ~Value, type = "scatter3d", mode = "markers")

grid_size = 0.1


border_init = c(0,0,2,3)

border_connect = c(2,0.5,3,2.5)

init_coef = samplePatch_Unrestricted(border_init, 0.5)
connect_coef = samplePatch_Known_Left(border = border_connect, k = 0.5, left_coefs = init_coef, left_border = border_init)




plotCubicPatch3D(init_coef, border_init, 0.1)



test_grid_init = expand.grid(seq(0,2, 0.1), seq(0,3,0.1))
PatchValues = apply(test_grid_init, 1, FUN = evaluateCubicPatchValue, coef = init_coef, border = border_init)

plottingGrid = data.frame(test_grid, PatchValues)
names(plottingGrid) = c('X','Y','Value')

test_grid_connect = expand.grid(seq(2,3, 0.1), seq(0.5,2.5,0.1))
PatchValues_Connect = apply(test_grid_connect, 1, FUN = evaluateCubicPatchValue, coef = connect_coef, border = border_connect)

PlottingGrid_Connect = data.frame(test_grid_connect, PatchValues_Connect)
names(PlottingGrid_Connect) = c('X','Y','Value')

PlottingGrid_Combined = rbind(plottingGrid, PlottingGrid_Connect)

Value <- xtabs(Value ~ X + Y, data = PlottingGrid_Combined)
X <- as.numeric(rownames(Value))
Y <- as.numeric(colnames(Value))

plot = plot_ly(x = ~X, y = ~Y, z = ~Value, type = "surface")

plot %>% layout(
  title = "3D Surface with Z-Axis Limit",
  scene = list(
    xaxis = list(title = 'X'),
    yaxis = list(title = 'Y'),
    zaxis = list(
      title = 'Value',
      range = c(-5, 6) # <-- This sets the Z-axis limit
    )
  )
)

plot_ly(PlottingGrid_Combined, x = ~X, y = ~Y, z = ~Value, type = "scatter3d", mode = "markers")



## Test Patches

border1 = c(0,0,5,1)
border2 = c(0,1,2,5)
border3 = c(2,1,5,3)
border3.5 = c(2,1,4,3)
border4 = c(4,1,5,3)
border5 = c(2,3,5,5)



coef1 = samplePatch_Unrestricted(border1, 0.2)
coef2 = samplePatch_Known_Bottom(border2, 0.2, coef1, border1, BR_coefs = coef1, BR_border = border1)
coef3 = samplePatch_Known_Bottom_and_Left(border3, 0.2, coef1, border1, coef1, border1, coef2, border2)
coef5 = samplePatch_Known_Bottom_and_Left(border5, 0.2, coef3, border3, coef3, border3, coef2, border2)

coef3.5 = samplePatch_Known_Corners(border3.5, k = 0.2, coef1, border1, coef2, border2, coef5, border5)
coef4 = samplePatch_Known_Corners(border4, k = 0.2, coef1, border1, coef3.5, border3.5, coef5, border5)

grid_size = 0.1

grid1 = expand.grid(seq(border1[1],border1[3], grid_size), seq(border1[2],border1[4], grid_size))
grid2 = expand.grid(seq(border2[1],border2[3], grid_size), seq(border2[2],border2[4], grid_size))
grid3 = expand.grid(seq(border3[1],border3[3], grid_size), seq(border3[2],border3[4], grid_size))
grid5 = expand.grid(seq(border5[1],border5[3], grid_size), seq(border5[2],border5[4], grid_size))

grid3.5 = expand.grid(seq(border3.5[1],border3.5[3], grid_size), seq(border3.5[2],border3.5[4], grid_size))
grid4 = expand.grid(seq(border4[1],border4[3], grid_size), seq(border4[2],border4[4], grid_size))

PatchValues1 = apply(grid1, 1, FUN = evaluateCubicPatchValue, coef = coef1, border = border1)
PatchValues2 = apply(grid2, 1, FUN = evaluateCubicPatchValue, coef = coef2, border = border2)
PatchValues3 = apply(grid3, 1, FUN = evaluateCubicPatchValue, coef = coef3, border = border3)
PatchValues5 = apply(grid5, 1, FUN = evaluateCubicPatchValue, coef = coef5, border = border5)

PatchValues3.5 = apply(grid3.5, 1, FUN = evaluateCubicPatchValue, coef = coef3.5, border = border3.5)
PatchValues4 = apply(grid4, 1, FUN = evaluateCubicPatchValue, coef = coef4, border = border4)


plottingGrid = data.frame(round(rbind(grid1, grid2,grid3,grid5), 1), c(PatchValues1, PatchValues2, PatchValues3, PatchValues5))
plottingGrid_2 = data.frame(round(rbind(grid1, grid2, grid3.5, grid4,grid5),1), c(PatchValues1, PatchValues2, PatchValues3.5, PatchValues4, PatchValues5))

names(plottingGrid) = c('X','Y','Value')
names(plottingGrid_2) = c('X','Y','Value')

duplicate_rows <- duplicated(plottingGrid[,c(1,2)])
duplicate_rows2 <- duplicated(plottingGrid_2[,c(1,2)])

plottingGrid = plottingGrid[!duplicate_rows,]
plottingGrid_2 = plottingGrid_2[!duplicate_rows2,]

plot_ly(plottingGrid, x = ~X, y = ~Y, z = ~Value, type = "scatter3d", mode = "markers")
plot_ly(plottingGrid_2, x = ~X, y = ~Y, z = ~Value, type = "scatter3d", mode = "markers")
plot_ly(plottingGrid_Compare, x = ~X, y = ~Y, z = ~Value, type = "scatter3d", mode = "markers")


mean(abs(plottingGrid[order(plottingGrid$X, plottingGrid$Y),]$Value - plottingGrid_2[order(plottingGrid_2$X, plottingGrid_2$Y),]$Value) < 0.1)

plottingGrid_Sorted = plottingGrid[order(plottingGrid$X, plottingGrid$Y),]
plottingGrid_Sorted2 = plottingGrid_2[order(plottingGrid_2$X, plottingGrid_2$Y),]

plottingGrid_Compare = plottingGrid_Sorted
plottingGrid_Compare$Value = plottingGrid_Sorted$Value - plottingGrid_Sorted2$Value

Value <- xtabs(Value ~ X+Y, data = plottingGrid)
X <- as.numeric(rownames(Value))
Y <- as.numeric(colnames(Value))

Value2 <- xtabs(Value ~ X+Y, data = plottingGrid_2)
X2 <- as.numeric(rownames(Value2))
Y2 <- as.numeric(colnames(Value2))

plot_ly(x = ~X, y = ~Y, z = ~Value, type = "surface")
plot_ly(x = ~X2, y = ~Y2, z = ~Value2, type = "surface")

mean(plottingGrid == plottingGrid_2)

### Sampling trajectories with the new approach
# In Every region in the base grid, each model has a probability of being "on"
# This corresponds to a mixture model where with the given probability the corner quantities are drawn from on "on" distribution and the other an "off" distribution

#The 3rd dimension corresponds to a cell, in the grid.
#The first dimension corresponds to a corner in the cell. (BL, BR, TL, TR)
#The second dimension corresponds to a quantity. (Value, ParX, ParY, ParXY)
sample_NewPatches = function(ModelLogits, k = 0.1){

  num_models = ncol(ModelLogits)
  num_cells = nrow(ModelLogits)

  basePatches = list()

  for(m in 1:num_models){

    cur_modelPatches = array(dim = c(4,4,num_cells))

    for(c in 1:num_cells){

      cur_prob = 1/(1 + exp(-1*ModelLogits[c,m]))
      is_on = sample(c(T,F), 1, prob = c(cur_prob, 1 - cur_prob))

      if(is_on){
        cur_modelPatches[,,c] = mvrnorm(4, mu = c(1,0,0,0), Sigma = k^2*diag(4))
      } else{
        cur_modelPatches[,,c] = mvrnorm(4, mu = c(0,0,0,0), Sigma = (0.01)^2*diag(4))
      }


    }

    basePatches[[m]] = cur_modelPatches


  }

  basePatches


}

