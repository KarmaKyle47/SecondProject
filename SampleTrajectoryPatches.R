library(MASS)
library(plotly)

samplePatch_Unrestricted = function(border, k){

  dx = border[3] - border[1]
  dy = border[4] - border[2]

  A = matrix(c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
               0,0,dx,0,dx^2,0,dx^3,0,0,0,0,0,0,0,0,0,
               0,dy,0,0,0,0,0,0,dy^2,dy^3,0,0,0,0,0,0,
               0,0,0,dx*dy,0,dx^2*dy,0,dx^3*dy,0,0,dx*dy^2,dx*dy^3,dx^2*dy^2,dx^2*dy^3,dx^3*dy^2,dx^3*dy^3,
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

  mu = A_inv %*% c(1, rep(0, 15))

  sigma = k^2*(A_inv %*% t(A_inv))

  coef = mvrnorm(mu = mu, Sigma = sigma)

  coef

}

samplePatch_Known_Bottom = function(border, k, bottom_coefs, bottom_border){

  dx = border[3] - border[1]
  dy = border[4] - border[2]

  A = matrix(c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
               0,0,dx,0,dx^2,0,dx^3,0,0,0,0,0,0,0,0,0,
               0,dy,0,0,0,0,0,0,dy^2,dy^3,0,0,0,0,0,0,
               0,0,0,dx*dy,0,dx^2*dy,0,dx^3*dy,0,0,dx*dy^2,dx*dy^3,dx^2*dy^2,dx^2*dy^3,dx^3*dy^2,dx^3*dy^3,
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

  c00 = evaluateCubicPatchValue(coef = bottom_coefs, border = bottom_border, curPos = c(border[1], border[2]))
  c01 = evaluateCubicPatchParY(coef = bottom_coefs, border = bottom_border, curPos = c(border[1], border[2]))
  c10 = evaluateCubicPatchParX(coef = bottom_coefs, border = bottom_border, curPos = c(border[1], border[2]))
  c11 = evaluateCubicPatchParXY(coef = bottom_coefs, border = bottom_border, curPos = c(border[1], border[2]))

  value_BR = evaluateCubicPatchValue(coef = bottom_coefs, border = bottom_border, curPos = c(border[3], border[2]))
  par_X_BR = evaluateCubicPatchParX(coef = bottom_coefs, border = bottom_border, curPos = c(border[3], border[2]))
  par_Y_BR = evaluateCubicPatchParY(coef = bottom_coefs, border = bottom_border, curPos = c(border[3], border[2]))
  par_XY_BR = evaluateCubicPatchParXY(coef = bottom_coefs, border = bottom_border, curPos = c(border[3], border[2]))

  c20 = (-1*par_X_BR*dx + 3*value_BR - 3*c00 - 2*c10*dx)/(dx^2)
  c21 = (-1*par_XY_BR*dx + 3*par_Y_BR - 3*c01 - 2*c11*dx)/(dx^2)
  c30 = (par_X_BR*dx - 2*value_BR + 2*c00 + c10*dx)/(dx^3)
  c31 = (par_XY_BR*dx - 2*par_Y_BR + 2*c01 + c11*dx)/(dx^3)


  known_coefs = c(c00, c01, c10, c11, c20, c21, c30, c31)

  mu_trans = -1*(A[c(3,4,7,8,11,12,15,16), 1:8] %*% known_coefs)

  A_inv = solve(A[c(3,4,7,8,11,12,15,16), 9:16])

  mu = A_inv %*% mu_trans
  sigma = k^2*(A_inv %*% t(A_inv))

  coef = mvrnorm(mu = mu, Sigma = sigma)

  c(known_coefs, coef)

}

samplePatch_Known_Left = function(border, k, left_coefs, left_border){

  dx = border[3] - border[1]
  dy = border[4] - border[2]

  A = matrix(c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
               0,0,dx,0,dx^2,0,dx^3,0,0,0,0,0,0,0,0,0,
               0,dy,0,0,0,0,0,0,dy^2,dy^3,0,0,0,0,0,0,
               0,0,0,dx*dy,0,dx^2*dy,0,dx^3*dy,0,0,dx*dy^2,dx*dy^3,dx^2*dy^2,dx^2*dy^3,dx^3*dy^2,dx^3*dy^3,
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

  c00 = evaluateCubicPatchValue(coef = left_coefs, border = left_border, curPos = c(border[1], border[2]))
  c01 = evaluateCubicPatchParY(coef = left_coefs, border = left_border, curPos = c(border[1], border[2]))
  c10 = evaluateCubicPatchParX(coef = left_coefs, border = left_border, curPos = c(border[1], border[2]))
  c11 = evaluateCubicPatchParXY(coef = left_coefs, border = left_border, curPos = c(border[1], border[2]))

  value_BR = evaluateCubicPatchValue(coef = left_coefs, border = left_border, curPos = c(border[1], border[4]))
  par_X_BR = evaluateCubicPatchParX(coef = left_coefs, border = left_border, curPos = c(border[1], border[4]))
  par_Y_BR = evaluateCubicPatchParY(coef = left_coefs, border = left_border, curPos = c(border[1], border[4]))
  par_XY_BR = evaluateCubicPatchParXY(coef = left_coefs, border = left_border, curPos = c(border[1], border[4]))

  c02 = (-1*par_Y_BR*dy + 3*value_BR - 3*c00 - 2*c01*dy)/(dy^2)
  c12 = (-1*par_XY_BR*dy + 3*par_X_BR - 3*c10 - 2*c11*dy)/(dy^2)
  c03 = (par_Y_BR*dy - 2*value_BR + 2*c00 + c01*dy)/(dy^3)
  c13 = (par_XY_BR*dy - 2*par_X_BR + 2*c10 + c11*dy)/(dy^3)


  known_coefs = c(c00, c01, c10, c11, c02, c03, c12, c13)

  mu_trans = -1*(A[c(2,4,6,8,10,12,14,16), c(1:4,9:12)] %*% known_coefs)

  A_inv = solve(A[c(2,4,6,8,10,12,14,16), c(5:8,13:16)])

  mu = A_inv %*% mu_trans
  sigma = k^2*(A_inv %*% t(A_inv))

  coef = mvrnorm(mu = mu, Sigma = sigma)

  c(known_coefs[1:4], coef[1:4], known_coefs[5:8], coef[5:8])


}

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

border_init = c(0,0,2,3)

init_coef = samplePatch_Unrestricted(border_init, 0.5)

test_grid = expand.grid(seq(0,2, 0.1), seq(0,3,0.1))

PatchValues = apply(test_grid, 1, FUN = evaluateCubicPatchValue, coef = init_coef, border = border_init)

plottingGrid = data.frame(test_grid, PatchValues)
names(plottingGrid) = c('X','Y','Value')

ggplot(plottingGrid, aes(x = X, y = Y, color = Value)) + geom_point(shape = 15, size = 20)

plot_ly(plottingGrid, x = ~X, y = ~Y, z = ~Value, type = "scatter3d", mode = "markers")

grid_size = 0.1

plotCubicPatch3D = function(coef, border, grid_size){

  test_grid = expand.grid(seq(border[1],border[3], grid_size), seq(border[2],border[4],grid_size))
  PatchValues = apply(test_grid, 1, FUN = evaluateCubicPatchValue, coef = coef, border = border)

  plottingGrid = data.frame(test_grid, PatchValues)
  names(plottingGrid) = c('X','Y','Value')

  Value <- xtabs(Value ~ X + Y, data = plottingGrid)
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

}

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






