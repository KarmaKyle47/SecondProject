
library(MASS)
library(plotly)

border = c(0,0,2,3)

k=0.5

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

curPos = c(1,1)


evaluateCubicPatchValue = function(coef, border, curPos){

  x_0 = border[1]
  y_0 = border[2]
  x = curPos[1]
  y = curPos[2]

  polyTerms = c(1, (y-y_0), (x-x_0), (x-x_0)^2,
                (x-x_0)^2, (x-x_0)^2*(y-y_0), (x-x_0)^3, (x-x_0)^3*(y-y_0),
                (y-y_0)^2, (y-y_0)^3, (x-x_0)*(y-y_0)^2, (x-x_0)*(y-y_0)^3,
                (x-x_0)^2*(y-y_0)^2, (x-x_0)^2*(y-y_0)^3, (x-x_0)^3*(y-y_0)^2, (x-x_0)^3*(y-y_0)^3)

  sum(coef * polyTerms)

}


test_grid = expand.grid(seq(0,2, 0.1), seq(0,3,0.1))

PatchValues = apply(test_grid, 1, FUN = evaluateCubicPatchValue, coef = coef, border = border)

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
        range = c(-200, 200) # <-- This sets the Z-axis limit
      )
    )
  )

}

rand_coef = samplePatch_Unrestricted(border, 10)
plotCubicPatch3D(rand_coef, border, 0.1)

