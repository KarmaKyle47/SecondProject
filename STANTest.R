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
