library(Rcpp)
library(RcppArmadillo)

cpp_code <- '
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// -------------------------------------------------------
// Helper to calculate the combined weighted vector field
// -------------------------------------------------------
arma::mat calc_drift(const arma::mat& pos_t, const arma::mat& beta_mat, const arma::vec& border) {
    int N = pos_t.n_rows;
    arma::mat drift(N, 2);
    double Lx = border[2] - border[0];
    double Ly = border[3] - border[1];

    int num_bases = beta_mat.n_rows;
    int M = std::round(std::sqrt(num_bases)) - 1;
    arma::vec omega = arma::regspace(0, M) * arma::datum::pi;

    // 1. Calculate Phi matrix
    arma::mat Phi(N, num_bases, arma::fill::ones);

    for(int i = 0; i <= M; ++i) {
        for(int j = 0; j <= M; ++j) {
            double scale = 1.0;
            if(i != 0) scale *= std::sqrt(2.0);
            if(j != 0) scale *= std::sqrt(2.0);

            int col_idx = j * (M + 1) + i;
            for(int r = 0; r < N; ++r) {
                double cx = std::cos(omega[i] * (pos_t(r, 1) - border[0]) / Lx);
                double cy = std::cos(omega[j] * (pos_t(r, 2) - border[1]) / Ly);
                Phi(r, col_idx) = scale * cx * cy;
            }
        }
    }

    arma::mat log_Traj = Phi * beta_mat;

    // 2. Base Vector Field + Weighting
    for(int r = 0; r < N; ++r) {
        double x = pos_t(r, 1);
        double y = pos_t(r, 2);
        // Added 1e-6 epsilon to prevent division by zero near origin
        double c = std::sqrt(x*x + y*y) + 1e-6;

        double f1x = y / c;
        double f1y = -x / c;
        double f2x = x / c;
        double f2y = y / c;

        double exp_L1 = std::exp(log_Traj(r, 0));
        double exp_L2 = std::exp(log_Traj(r, 1));

        drift(r, 0) = exp_L1 * f1x + exp_L2 * f2x;
        drift(r, 1) = exp_L1 * f1y + exp_L2 * f2y;
    }

    return drift;
}

// -------------------------------------------------------
// Main rMAP Loss Function
// -------------------------------------------------------
// [[Rcpp::export]]
double spline_rMAP_loss_cpp(
    arma::vec c_params,           // Optimizable control points (length 2 * N_param)
    arma::mat Phi_data,           // Basis evaluated at obs times (N_obs x N_param)
    arma::mat Phi_quad,           // Basis at dense grid (N_quad x N_param)
    arma::mat Phi_deriv_quad,     // Deriv basis at dense grid (N_quad x N_param)
    arma::mat aug_X,              // Augmented positional data (N_obs x 2)
    arma::mat sampledTrajectory,  // Eulerian weight surface (M_sq x 2 matrix)
    arma::vec border,             // Border for HSGP scaling (length 4)
    double pos_sd,                // Positional error SD
    double vel_sd,                // Velocity error SD
    arma::vec rand_vel_x,         // Frozen stochastic velocity errors (X)
    arma::vec rand_vel_y,         // Frozen stochastic velocity errors (Y)
    arma::vec prior_c,            // The fixed projected prior control points
    arma::mat inv_c_prior_sigma,   // The inverse covariance for the control points
    double Lt
) {

    int N_param = Phi_data.n_cols;
    int N_obs = Phi_data.n_rows;
    int N_quad = Phi_quad.n_rows;

    // 1. Unpack Control Points
    arma::vec c_x = c_params.subvec(0, N_param - 1);
    arma::vec c_y = c_params.subvec(N_param, 2 * N_param - 1);
    arma::mat C = arma::join_rows(c_x, c_y);

    // =======================================================
    // PART 1: POSITIONAL LOSS (At data observation times)
    // =======================================================

    arma::mat X_pred_data = Phi_data * C;

    double nll_pos = 0.0;
    for(int i = 0; i < N_obs; i++) {
        double err_x = X_pred_data(i, 0) - aug_X(i, 0);
        double err_y = X_pred_data(i, 1) - aug_X(i, 1);
        nll_pos += (err_x * err_x + err_y * err_y) / (2.0 * pos_sd * pos_sd);
    }

    // =======================================================
    // PART 2: VELOCITY/PHYSICS LOSS (At dense quadrature times)
    // =======================================================

    arma::mat X_quad_pred = Phi_quad * C;
    arma::mat V_spline_pred = Phi_deriv_quad * C;

    // Format pos_t for calc_drift (N_quad x 3 matrix where col 1 is X, col 2 is Y)
    // (Col 0 is initialized to 0.0 as time is unused by the spatial field)
    arma::mat pos_t(N_quad, 3, arma::fill::zeros);
    pos_t.col(1) = X_quad_pred.col(0);
    pos_t.col(2) = X_quad_pred.col(1);

    // Compute the physical vector field directly using your helper function
    // Calculates velocities for all N_quad points simultaneously
    arma::mat V_RK4 = calc_drift(pos_t, sampledTrajectory, border);

    double nll_vel = 0.0;
    for(int i = 0; i < N_quad; i++) {
        // Add the frozen stochastic error for this specific optimizer run
        double target_vx = V_RK4(i, 0) + rand_vel_x(i);
        double target_vy = V_RK4(i, 1) + rand_vel_y(i);

        // Penalize the spline derivative for fighting the targeted physical flow
        double diff_vx = V_spline_pred(i, 0) - target_vx;
        double diff_vy = V_spline_pred(i, 1) - target_vy;

        nll_vel += (diff_vx * diff_vx + diff_vy * diff_vy) / (2.0 * vel_sd * vel_sd);
    }

    nll_vel = nll_vel * (Lt / N_quad);

    // =======================================================
    // PART 3: PRIOR PENALTY LOSS
    // =======================================================

    arma::vec d_cx = c_x - prior_c.subvec(0, N_param - 1);
    arma::vec d_cy = c_y - prior_c.subvec(N_param, 2 * N_param - 1);

    // Quadratic form for the prior penalty
    double prior_loss = arma::dot(d_cx, inv_c_prior_sigma * d_cx) +
                        arma::dot(d_cy, inv_c_prior_sigma * d_cy);

    prior_loss = 0.5 * prior_loss;

    return nll_pos + nll_vel + prior_loss;
}
'

# Compile the C++ code into R
Rcpp::sourceCpp(code = cpp_code)
