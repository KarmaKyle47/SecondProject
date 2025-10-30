// File: linear_regression.stan
data {
  int<lower=0> N;      // Number of observations
  vector[N] x;         // Predictor variable
  vector[N] y;         // Outcome variable
}

parameters {
  real alpha;          // Intercept
  real beta;           // Slope
  real<lower=0> sigma; // Residual standard deviation
}

model {
  // Priors
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
  sigma ~ cauchy(0, 5); // A weakly informative prior for scale parameters

  // Likelihood
  y ~ normal(alpha + beta * x, sigma);
}

generated quantities {
  vector[N] y_rep; // For posterior predictive checks
  for (n in 1:N) {
    y_rep[n] = normal_rng(alpha + beta * x[n], sigma);
  }
}