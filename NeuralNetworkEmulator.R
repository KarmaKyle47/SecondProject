# Load ggplot2 for the visualization step at the end
library(ggplot2)

#' Generate a 3D Gridded Synthetic Vector Field
#'
#' @param nx Number of spatial grid points in x
#' @param ny Number of spatial grid points in y
#' @param nt Number of time steps
#' @param x_range Spatial boundaries for x axis
#' @param y_range Spatial boundaries for y axis
#' @param t_range Time boundaries
#' @return A data.frame with columns: t, x, y, u, v
generate_synthetic_vf <- function(nx = 40, ny = 40, nt = 10,
                                  x_range = c(-10, 10),
                                  y_range = c(-10, 10),
                                  t_range = c(0, 5)) {

  # 1. Create the 3D Grid
  x_seq <- seq(x_range[1], x_range[2], length.out = nx)
  y_seq <- seq(y_range[1], y_range[2], length.out = ny)
  t_seq <- seq(t_range[1], t_range[2], length.out = nt)

  df <- expand.grid(x = x_seq, y = y_seq, t = t_seq)

  # 2. Background Flow: A meandering current that pulses over time
  # u_bg oscillates based on y and t
  # v_bg oscillates based on x
  U_bg <- 1.0 + 0.5 * cos(pi * df$y / 5) * cos(2 * pi * df$t / max(t_seq))
  V_bg <- 0.5 * sin(pi * df$x / 5)

  # 3. Eddy 1: Cyclonic (Counter-Clockwise), moving Eastward
  # Its center (x0, y0) shifts based on time t
  eddy1_x <- -5 + 2 * df$t
  eddy1_y <- 2
  r2_1 <- (df$x - eddy1_x)^2 + (df$y - eddy1_y)^2
  sigma2 <- 4.0   # Controls the size/spread of the eddies
  A1 <- 3.0       # Amplitude of the rotation

  # Cartesian curl equations for a Gaussian vortex
  U_e1 <- -A1 * (df$y - eddy1_y) * exp(-r2_1 / (2 * sigma2))
  V_e1 <-  A1 * (df$x - eddy1_x) * exp(-r2_1 / (2 * sigma2))

  # 4. Eddy 2: Anticyclonic (Clockwise), moving Westward
  eddy2_x <- 5 - 1.5 * df$t
  eddy2_y <- -3
  r2_2 <- (df$x - eddy2_x)^2 + (df$y - eddy2_y)^2
  A2 <- -2.5      # Negative amplitude reverses rotation

  U_e2 <- -A2 * (df$y - eddy2_y) * exp(-r2_2 / (2 * sigma2))
  V_e2 <-  A2 * (df$x - eddy2_x) * exp(-r2_2 / (2 * sigma2))

  # 5. Combine all components
  df$u <- U_bg + U_e1 + U_e2
  df$v <- V_bg + V_e1 + V_e2

  return(df)
}

# --- Generate the Data ---
# Creates a dataframe with 16,000 rows (40 x 40 x 10)
grid_data <- generate_synthetic_vf()


# Filter for just the first time step
plot_data <- subset(grid_data, t == 3 + 1/3)

# Calculate speed for color mapping
plot_data$speed <- sqrt(plot_data$u^2 + plot_data$v^2)

# Plot the vector field using segments
ggplot(plot_data, aes(x = x, y = y)) +
  # Add arrows. We scale down u and v slightly so the arrows don't overlap too much
  geom_segment(aes(xend = x + u * 0.3, yend = y + v * 0.3, color = speed),
               arrow = arrow(length = unit(0.1, "cm")), size = 0.5) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(title = "Synthetic Vector Field (t = 0)",
       subtitle = "Notice the two distinct eddies and the background flow",
       x = "X Position", y = "Y Position", color = "Speed") +
  coord_fixed()


# Install gganimate and transformr if you haven't already
# install.packages(c("gganimate", "transformr"))

library(ggplot2)
library(gganimate)
install.packages(c("gifski", "av"))

# 1. Calculate the speed for the color mapping (if not done already)
grid_data$speed <- sqrt(grid_data$u^2 + grid_data$v^2)

# 2. Build the base plot (notice we don't filter for t == 0 anymore)
p <- ggplot(grid_data, aes(x = x, y = y)) +
  geom_segment(aes(xend = x + u * 0.3, yend = y + v * 0.3, color = speed),
               arrow = arrow(length = unit(0.1, "cm")), size = 0.5) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  coord_fixed() +
  labs(title = "Synthetic Ocean Currents",
       # The {frame_time} variable dynamically updates the title during animation!
       subtitle = "Time: {round(frame_time, 2)}",
       x = "X Position", y = "Y Position", color = "Speed")

# 3. Add the animation layer
# transition_time tells gganimate which variable controls the frames
anim <- p +
  transition_time(t) +
  ease_aes('linear')

# 4. Render the animation
# nframes controls smoothness; fps controls speed
animated_gif <- animate(anim, nframes = 100, fps = 20, width = 600, height = 600)

# Display it in the viewer
animated_gif

# Optional: Save it to your working directory
# anim_save("synthetic_eddies.gif", animation = animated_gif)



### Emulator



install.packages("keras3")
library(keras3)

install_keras(backend = "tensorflow", version = "2.10")

install_keras()

# 1. Standardize your data (same as before)
t_mean <- mean(grid_data$t); t_sd <- sd(grid_data$t)
x_mean <- mean(grid_data$x); x_sd <- sd(grid_data$x)
y_mean <- mean(grid_data$y); y_sd <- sd(grid_data$y)

t_stand = (grid_data$t - t_mean) / t_sd
x_stand = (grid_data$x - x_mean) / x_sd
y_stand = (grid_data$y - y_mean) / y_sd

X_train = matrix(data = c(t_stand, x_stand, y_stand), ncol = 3, byrow = F)
Y_train <- as.matrix(grid_data[,c('u','v')])

# 2. Build the Model
model <- keras_model_sequential() %>%
  layer_dense(units = 32, activation = "tanh", input_shape = c(3)) %>%
  layer_dense(units = 32, activation = "tanh") %>%
  layer_dense(units = 2) # Velocity u, v (Linear output)

# 3. Compile
model %>% compile(
  optimizer = optimizer_adam(learning_rate = 0.01),
  loss = "mse"
)

# 4. Train
history <- model %>% fit(
  X_train, Y_train,
  epochs = 500,
  batch_size = 128,
  view_metrics = TRUE # Opens a live plot in RStudio
)





