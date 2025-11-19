library(plotly)

# --- 1. Generator Function (Random State) ---
generate_hsgp_surface <- function(n_grid = 60, L_factor = 2, M = 20,
                                  length_scale = 0.2, magnitude = 1.0) {

  # Define Grid (0 to 1)
  x <- seq(0, 1, length.out = n_grid)
  y <- seq(0, 1, length.out = n_grid)

  # HSGP Hyperparameters
  L <- L_factor
  omega <- ((1:M) * pi) / (2 * L)

  # Spectral Density (Power Spectrum)
  spd <- magnitude^2 * sqrt(2 * pi) * length_scale * exp(-0.5 * (length_scale * omega)^2)

  # Random Coefficients (Standard Normals)
  # This is where the variation comes from
  z_matrix <- matrix(rnorm(M * M), nrow = M, ncol = M)

  # Scale coefficients by spectral density
  beta <- diag(sqrt(spd)) %*% z_matrix %*% diag(sqrt(spd))

  # Basis Functions (Sine Waves)
  scale_factor <- 1 / sqrt(L)
  phi_x <- outer(x - 0.5, 1:M, function(x, m) sin(m * pi * (x + L) / (2 * L))) * scale_factor
  phi_y <- outer(y - 0.5, 1:M, function(y, m) sin(m * pi * (y + L) / (2 * L))) * scale_factor

  # Tensor Product
  surface <- phi_x %*% beta %*% t(phi_y)

  return(list(x = x, y = y, z = surface))
}

# --- 2. The Plotting Function ---
get_hsgp_plots <- function(n_grid = 100) {

  # --- Generate Random Data ---

  # 1. Physics (Coefficient) Surface
  # High length_scale (0.3) = Smooth features
  sim_phys <- generate_hsgp_surface(n_grid, length_scale = 0.10, magnitude = 0.3, M = 50)
  mat_phys <- sim_phys$z + 1.0 # Centered at 1.0

  # Grab the coordinate vectors (0 to 1) for plotting
  x_vec <- sim_phys$x
  y_vec <- sim_phys$y

  # 2. Logit Surface (The Switch)
  # Low length_scale (0.15) = Sharp features / Fronts
  sim_logit <- generate_hsgp_surface(n_grid, length_scale = 0.10, magnitude = 3.0, M = 50)
  mat_logit <- sim_logit$z - 2.0 # Centered at -2.0 (Default Off)

  # 3. Calculate Weight and Final
  mat_weight <- 1 / (1 + exp(-mat_logit)) # Sigmoid
  mat_final <- mat_phys * mat_weight      # Element-wise multiplication

  # --- Plotting Setup ---
  ax <- list(title = "", showgrid = FALSE, showbackground = FALSE, showticklabels = TRUE)
  layout_settings <- list(scene = list(xaxis = ax, yaxis = ax, zaxis = ax))

  # --- Create Individual Plots ---

  # Plot 1: Physics (Viridis)
  p1 <- plot_ly(x = ~x_vec, y = ~y_vec, z = ~mat_phys, type = "surface", colorscale = "Viridis") %>%
    layout(title = "1. Coefficient Surface (Physics)", scene = layout_settings$scene)

  # Plot 2: Logit (RdBu)
  p2 <- plot_ly(x = ~x_vec, y = ~y_vec, z = ~mat_logit, type = "surface", colorscale = "RdBu") %>%
    layout(title = "2. Logit Surface (The Switch)", scene = layout_settings$scene)

  # Plot 3: Weight (Blues)
  p3 <- plot_ly(x = ~x_vec, y = ~y_vec, z = ~mat_weight, type = "surface", colorscale = "Blues") %>%
    layout(title = "3. Weight Surface (0 to 1)", scene = layout_settings$scene)

  # Plot 4: Final Result (Plasma)
  p4 <- plot_ly(x = ~x_vec, y = ~y_vec, z = ~mat_final, type = "surface", colorscale = "Plasma") %>%
    layout(title = "4. Final Surface (Result)", scene = layout_settings$scene)

  return(list(physics = p1, logit = p2, weight = p3, final = p4))
}

# --- Usage ---
# Run this line repeatedly to see new random surfaces
my_plots <- get_hsgp_plots()

# View them
my_plots$physics
my_plots$logit
my_plots$weight
my_plots$final
