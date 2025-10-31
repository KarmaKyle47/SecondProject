install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
install.packages(c("bayesplot", "loo", "posterior", "cmdstanr"))

# Automatically cache compiled Stan models to avoid recompilation
rstan_options(auto_write = TRUE)

# Use multiple cores for parallel chain execution
options(mc.cores = parallel::detectCores())


## Linear Regression Test

# Set seed for reproducibility
set.seed(42)

# Simulate data
N <- 100
alpha_true <- 1.5
beta_true <- 2.0
sigma_true <- 3.0
x <- rnorm(N, 5, 2)
y <- alpha_true + beta_true * x + rnorm(N, 0, sigma_true)

# Prepare data for Stan as a named list
stan_data <- list(N = N, x = x, y = y)

fit_linear <- stan(
  file = "STAN_Files/linear_regression.stan", # Path to the Stan program
  data = stan_data,               # Named list of data
  chains = 32,                     # Number of Markov chains
  iter = 2000,                    # Total iterations per chain
  warmup = 1000,                  # Burn-in iterations
  cores = 32                       # Number of cores for parallel execution
)

print(fit_linear, pars = c("alpha", "beta", "sigma"), probs = c(0.025, 0.5, 0.975))


library(cmdstanr)

# --- 1. Create Mock Data for Testing ---
# Let's test with a simple 2x2 grid (comp_res = 2)
# and 3 total models (numModels = 3)

comp_res_test <- 10
num_models_test <- 5
N_test <- comp_res_test^2 # N = 4

# Create the 4 base cells:
# (0,0) to (1,1)
# (1,0) to (2,1)
# (0,1) to (1,2)
# (1,1) to (2,2)

baseTree = generate_grid_tree(0.1, c(0,0,1,1))
baseTree$boundaries = baseTree$boundaries[order(baseTree$boundaries[,2], baseTree$boundaries[,1]),]

base_boundaries = baseTree$boundaries

base_models = sample_models_one_pass(baseTree, num_models = 5, baseWeight = 0.1)[,c(6,7)]

# --- 2. Assemble the Stan Data List ---
# (This section is identical)
stan_data <- list(
  comp_res = comp_res_test,
  numModels = num_models_test,
  trans_prop = 8/9,
  N = N_test,
  baseCompGridBoundaries = base_boundaries,
  models = base_models,
  y_dummy = 1.0
)

# --- 3. Compile the Stan Model ---
# rstan separates compiling and running. This just compiles the C++.
# (This only needs to be done once)
file <- file.path(getwd(), "STAN_Files/test_functions.stan")
mod <- rstan::stan_model(file)

# --- 4. Run the Model to Test the Function ---
# We use algorithm = "Fixed_param"
# This is the rstan equivalent of cmdstanr's fixed_param = TRUE
fit <- rstan::sampling(
  object = mod,         # The compiled model object
  data = stan_data,
  algorithm = "Fixed_param", # The key change!
  iter = 1,             # We only need one "draw"
  warmup = 0,
  chains = 1
)

# --- 5. Extract and Inspect the Result ---
# rstan uses the extract() function
fit_extract <- rstan::extract(fit)

# 'final_grid' will be in the extracted list as a 3D array: [iteration, row, col]
# We just want the first (and only) iteration
final_matrix <- fit_extract$final_grid[1, , ]

# Print the final matrix
# Dimensions should be:
# (N*9) = 36 rows
# (5 + numModels) = 8 columns
print(dim(final_matrix))
print(final_matrix)

# Print the model columns (6, 7, 8)
print("Model Columns:")
print(final_matrix[, 6:8])

Borders = data.frame(L1 = final_matrix[,1],
                     L2 = final_matrix[,2],
                     U1 = final_matrix[,3],
                     U2 = final_matrix[,4])

ShadeData = data.frame(xmin = final_matrix[,1],
                       xmax = final_matrix[,3],
                       ymin = final_matrix[,2],
                       ymax = final_matrix[,4],
                       Region = final_matrix[,7])

plot = ggplot() +
  geom_segment(data = Borders, aes(x = L1, y = L2, xend = L1, yend = U2), color = 'black') +
  geom_segment(data = Borders, aes(x = U1, y = L2, xend = U1, yend = U2), color = 'black') +
  geom_segment(data = Borders, aes(x = L1, y = L2, xend = U1, yend = L2), color = 'black') +
  geom_segment(data = Borders, aes(x = L1, y = U2, xend = U1, yend = U2), color = 'black') +
  geom_rect(data = ShadeData,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                fill = factor(Region)),
            inherit.aes = FALSE,
            alpha = 0.4) +
  xlab("Longitude") + ylab("Latitude")

plot


plotTransitionRegions(treeBoundaries = newBoundaries[,1:4], types = newBoundaries$Type)

















