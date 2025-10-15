# Title: Uniform Grid Transformation with Leader-Line Numbering
# Author: Gemini
# Date: 2025-10-11
# Description: This definitive script uses a robust leader-line labeling
# system to prevent overlapping numbers in crowded regions of the physical
# domain, ensuring a perfectly clear visualization of cell correspondence.

# --- 1. Install and Load Required Libraries ---
# install.packages(c("MASS", "akima"))
library(MASS)
library(akima)

# --- 2. Generate Irregularly Spaced Data ---
set.seed(42)
n_points_per_cluster <- 250
x1 <- rnorm(n_points_per_cluster, mean = 0.2, sd = 0.08)
y1 <- rnorm(n_points_per_cluster, mean = 0.2, sd = 0.08)
x2 <- rnorm(n_points_per_cluster, mean = 0.8, sd = 0.08)
y2 <- rnorm(n_points_per_cluster, mean = 0.8, sd = 0.08)
data <- data.frame(x = c(x1, x2), y = c(y1, y2))
data <- subset(data, x >= 0 & x <= 1 & y >= 0 & y <= 1)

# --- 3. Estimate Data Density and Calculate CDFs ---
grid_size_cdf <- 256
kde <- kde2d(data$x, data$y, n = grid_size_cdf, lims = c(0, 1, 0, 1))
rho <- kde$z
x_coords_cdf <- kde$x
y_coords_cdf <- kde$y

rho <- rho / sum(rho)
rho_x <- colSums(rho)
cdf_x_vec <- cumsum(rho_x)
eta <- matrix(rep(cdf_x_vec, each = grid_size_cdf), nrow = grid_size_cdf)
rho_y_given_x <- sweep(rho, 2, rho_x, FUN = "/")
rho_y_given_x[is.na(rho_y_given_x)] <- 0
xi <- apply(rho_y_given_x, 2, cumsum)

# --- 4. Visualization with Leader-Line Numbering ---

n_vis_cells <- 5
levels_to_draw <- seq(0, 1, length.out = n_vis_cells + 1)

# Create a matrix that maps each physical point to a cell ID
cell_id_matrix <- matrix(NA, nrow = grid_size_cdf, ncol = grid_size_cdf)
cell_number <- 1
for (j in 1:n_vis_cells) {
  for (i in 1:n_vis_cells) {
    eta_min <- levels_to_draw[j]; eta_max <- levels_to_draw[j + 1]
    xi_min <- levels_to_draw[i]; xi_max <- levels_to_draw[i + 1]
    indices <- which(eta >= eta_min & eta <= eta_max & xi >= xi_min & xi <= xi_max, arr.ind = TRUE)
    cell_id_matrix[indices] <- cell_number
    cell_number <- cell_number + 1
  }
}

par(mfrow = c(1, 2), pty = "s", mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 2, 0))

# Plot 1: Original Domain
plot(NULL, asp = 1, xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1), # Expanded limits for leader lines
     xlab = "x (Physical)", ylab = "y (Physical)", main = "Original Domain & Adaptive Grid")
contour(x_coords_cdf, y_coords_cdf, cell_id_matrix,
        levels = 1:(n_vis_cells^2), add = TRUE, col = "gray70", drawlabels = FALSE)
points(data$x, data$y, pch = 20, col = "gray30", cex = 0.5)

# Plot 2: Transformed Domain
plot(NULL, asp = 1, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "η (Computational)", ylab = "ξ (Computational)", main = "Transformed Domain & Uniform Grid")
abline(h = levels_to_draw, col = "gray70")
abline(v = levels_to_draw, col = "gray70")
transformed_x <- interp2(x_coords_cdf, y_coords_cdf, eta, data$x, data$y)
transformed_y <- interp2(x_coords_cdf, y_coords_cdf, xi, data$x, data$y)
points(transformed_x, transformed_y, pch = 20, col = "gray30", cex = 0.5)

# --- Smart Labeling Logic ---
# Define a threshold for what is considered a "small" cell
# (based on the number of high-res grid points it contains)
area_threshold <- 100
leader_line_length <- 0.15

cell_number <- 1
for (j in 1:n_vis_cells) {
  for (i in 1:n_vis_cells) {
    # Labeling for Transformed (Right) Plot is always simple
    eta_center <- (levels_to_draw[j] + levels_to_draw[j + 1]) / 2
    xi_center <- (levels_to_draw[i] + levels_to_draw[i + 1]) / 2
    text(eta_center, xi_center, cell_number, col = "blue", font = 2)

    # Smart Labeling for Original (Left) Plot
    indices <- which(cell_id_matrix == cell_number, arr.ind = TRUE)
    if (nrow(indices) > 0) {
      cell_area <- nrow(indices)
      x_center_phys <- mean(x_coords_cdf[indices[, 2]])
      y_center_phys <- mean(y_coords_cdf[indices[, 1]])

      if (cell_area < area_threshold) {
        # --- For SMALL cells, use leader lines ---
        points(x_center_phys, y_center_phys, pch = 21, bg = "white", col="blue", cex=0.7) # Add a dot

        # Calculate angle to push label away from center of plot
        angle <- atan2(y_center_phys - 0.5, x_center_phys - 0.5)

        # Calculate leader line end point
        x_end <- x_center_phys + leader_line_length * cos(angle)
        y_end <- y_center_phys + leader_line_length * sin(angle)

        segments(x_center_phys, y_center_phys, x_end, y_end, col = "blue", lwd = 0.8)
        text(x_end, y_end, cell_number, col = "blue", font = 2, pos = ifelse(abs(cos(angle)) > 0.5,
                                                                             ifelse(cos(angle) > 0, 4, 2),
                                                                             ifelse(sin(angle) > 0, 3, 1)))
      } else {
        # --- For LARGE cells, label the center directly ---
        text(x_center_phys, y_center_phys, cell_number, col = "blue", font = 2)
      }
    }
    cell_number <- cell_number + 1
  }
}
mtext("Numbered Cell Correspondence with Leader Lines", outer = TRUE, cex = 1.5)

