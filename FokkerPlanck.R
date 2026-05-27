library(ggplot2)
library(gganimate)
library(dplyr)
library(tidyr)
library(gifski)

# ==========================================
# 1. Configuration
# ==========================================

nx <- 1000; ny <- 1000
x_lim <- c(-5, 5); y_lim <- c(-5, 5)
T_max <- 2.0
dt <- 0.001
D <- 0.001

n_particles <- 100000

dx <- (x_lim[2] - x_lim[1]) / (nx - 1)
dy <- (y_lim[2] - y_lim[1]) / (ny - 1)
x_vals <- seq(x_lim[1], x_lim[2], length.out = nx)
y_vals <- seq(y_lim[1], y_lim[2], length.out = ny)

# ==========================================
# 2. Vector Field & Helper Functions
# ==========================================

# get_u <- function(x, y) { -pi * sin(pi * x) * cos(pi * y) }
# get_v <- function(x, y) {  pi * cos(pi * x) * sin(pi * y) }

get_u <- function(x, y) { x/sqrt(x^2 + y^2) }
get_v <- function(x, y) { y/sqrt(x^2 + y^2) }

# --- Chang-Cooper Helper: Bernoulli Function ---
# B(x) = x / (exp(x) - 1)
# Used to weight the flux based on local Peclet number
Bernoulli <- function(x) {
  # Avoid division by zero for very small x using Taylor expansion
  # (Standard numerical trick for x/(e^x - 1))
  res <- x / (expm1(x))
  res[abs(x) < 1e-5] <- 1 - x[abs(x) < 1e-5]/2
  return(res)
}

# ==========================================
# 3. Initialization
# ==========================================

# Grid P
P <- matrix(0, nrow = nx, ncol = ny)
grid_x <- outer(x_vals, rep(1, ny))
grid_y <- outer(rep(1, nx), y_vals)
mu_x <- 1; mu_y <- 1; sigma <- 0.08
P <- exp(-((grid_x - mu_x)^2 + (grid_y - mu_y)^2) / (2 * sigma^2))
P <- P / (sum(P) * dx * dy)

# Particles
set.seed(42)
p_x <- rnorm(n_particles, mean = mu_x, sd = sigma)
p_y <- rnorm(n_particles, mean = mu_y, sd = sigma)

# ==========================================
# 4. Pre-Calculate Interface Velocities
# ==========================================
# Chang-Cooper requires velocity at the "walls" between cells
# X-interfaces: at x + dx/2 (Size: nx-1 x ny)
# Y-interfaces: at y + dy/2 (Size: nx x ny-1)

x_int <- x_vals[-nx] + dx/2
y_int <- y_vals[-ny] + dy/2

# Velocity at X-interfaces (between columns)
U_int <- outer(x_int, y_vals, get_u)

# Velocity at Y-interfaces (between rows)
V_int <- outer(x_vals, y_int, get_v)

# ==========================================
# 5. Simulation Loop
# ==========================================

save_interval <- 15
n_steps <- floor(T_max / dt)
combined_results <- list()

cat("Running Chang-Cooper FPE vs Monte Carlo...\n")

for (t in 1:n_steps) {

  # --- STEP 1: Chang-Cooper Solver ---

  # A. Flux in X Direction
  # J_{i+1/2} = (D/dx) * [ B(w)*P_i - B(-w)*P_{i+1} ]
  # w = u_{i+1/2} * dx / D

  w_x <- U_int * dx / D
  B_w <- Bernoulli(-w_x)
  B_mw <- Bernoulli(w_x) # B(-w) = B(w) + w

  # P_left is P[1:(n-1), ], P_right is P[2:n, ]
  P_left <- P[-nx, ]
  P_right <- P[-1, ]

  # Flux between cells (Inner fluxes)
  Jx_inner <- (D / dx) * (B_w * P_left - B_mw * P_right)

  # Pad with 0 for Reflective Boundaries (No flux at walls)
  Jx_full <- rbind(0, Jx_inner, 0)

  # B. Flux in Y Direction
  w_y <- V_int * dy / D
  B_wy <- Bernoulli(-w_y)
  B_mwy <- Bernoulli(w_y)

  P_down <- P[, -ny]
  P_up <- P[, -1]

  Jy_inner <- (D / dy) * (B_wy * P_down - B_mwy * P_up)

  # Pad with 0 for Reflective Boundaries
  Jy_full <- cbind(0, Jy_inner, 0)

  # C. Divergence Update
  # dP/dt = - (dJx/dx + dJy/dy)
  # Note: Chang-Cooper computes Flux J directly, so we just take Div(J).
  # We do NOT add D*Laplacian separately; diffusion is built into J.

  div_J <- (diff(Jx_full) / dx) + (t(diff(t(Jy_full))) / dy)

  P <- P - div_J * dt

  # (No need for "P[P<0] <- 0" anymore! Chang-Cooper handles it.)

  # --- STEP 2: Monte Carlo (Particles) ---

  u_p <- get_u(p_x, p_y)
  v_p <- get_v(p_x, p_y)

  noise <- sqrt(2 * D * dt)
  p_x <- p_x + u_p * dt + rnorm(n_particles) * noise
  p_y <- p_y + v_p * dt + rnorm(n_particles) * noise

  # Reflective Bounds
  p_x <- ifelse(p_x < x_lim[1], 2*x_lim[1] - p_x, p_x)
  p_x <- ifelse(p_x > x_lim[2], 2*x_lim[2] - p_x, p_x)
  p_y <- ifelse(p_y < y_lim[1], 2*y_lim[1] - p_y, p_y)
  p_y <- ifelse(p_y > y_lim[2], 2*y_lim[2] - p_y, p_y)

  # --- STEP 3: Save & Binning ---
  if (t %% save_interval == 0 || t == 1) {
    curr_time <- t * dt

    # Grid Data
    df_fpe <- expand.grid(x = x_vals, y = y_vals)
    df_fpe$density <- as.vector(P)
    df_fpe$method <- "1. Chang-Cooper FPE"
    df_fpe$time <- curr_time

    # Particle Binning
    x_bins <- cut(p_x, breaks = seq(x_lim[1], x_lim[2], length.out = nx + 1), labels = FALSE)
    y_bins <- cut(p_y, breaks = seq(y_lim[1], y_lim[2], length.out = ny + 1), labels = FALSE)

    valid <- !is.na(x_bins) & !is.na(y_bins)
    mat_mc <- matrix(0, nx, ny)
    if(any(valid)) {
      tbl <- table(factor(x_bins[valid], levels=1:nx),
                   factor(y_bins[valid], levels=1:ny))
      mat_mc <- as.matrix(tbl)
    }

    df_mc <- expand.grid(x = x_vals, y = y_vals)
    df_mc$density <- as.vector(mat_mc) / n_particles / (dx * dy)
    df_mc$method <- "2. Monte Carlo Empirical"
    df_mc$time <- curr_time

    combined_results[[length(combined_results) + 1]] <- rbind(df_fpe, df_mc)
  }

  svMisc::progress(t, n_steps)
}

sim_data <- do.call(rbind, combined_results)
sim_data <- sim_data[sim_data$density > 1e-5, ]


# ==========================================
# 6. Visualization
# ==========================================

p <- ggplot(sim_data, aes(x = x, y = y, fill = density)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", name = "Density") +

  # geom_contour(aes(z = density, group = interaction(time, method)),
  #              color = "white", alpha = 0.3, bins = 5) +

  facet_wrap(~method) +

  labs(title = 'Method Comparison: t = {frame_time}',
       subtitle = 'Chang-Cooper guarantees positivity without clamping',
       x = "Longitude", y = "Latitude") +

  theme_minimal() +
  theme(panel.background = element_rect(fill = "gray10"),
        plot.background = element_rect(fill = "gray10"),
        text = element_text(color = "white"),
        axis.text = element_text(color = "gray"),
        strip.text = element_text(size = 12, face = "bold", color = "white")) +

  transition_time(time)

animate(p, nframes = 100, fps = 15, width = 800, height = 400, renderer = gifski_renderer())





