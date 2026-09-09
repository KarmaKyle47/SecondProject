library(Rcpp)
library(RcppArmadillo)

cpp_code <- '
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// -------------------------------------------------------
// Optimized Helper to calculate drift for a SINGLE point
// -------------------------------------------------------
inline void single_drift(double x, double y, const arma::mat& beta_mat,
                         const arma::vec& border, const arma::vec& omega,
                         int M, double Lx, double Ly, double& dx, double& dy) {

    double log_T1 = 0.0;
    double log_T2 = 0.0;
    double cx[128]; // Assuming M < 128. Increase size if you use massive basis sets.
    double cy[128];

    // O(M) Precalculation
    for(int i = 0; i <= M; ++i) {
        double scale = (i == 0) ? 1.0 : std::sqrt(2.0);
        cx[i] = scale * std::cos(omega[i] * (x - border[0]) / Lx);
        cy[i] = scale * std::cos(omega[i] * (y - border[1]) / Ly);
    }

    // O(M^2) dot product without matrix allocation
    for(int j = 0; j <= M; ++j) {
        for(int i = 0; i <= M; ++i) {
            double phi_val = cx[i] * cy[j];
            int col_idx = j * (M + 1) + i;
            log_T1 += phi_val * beta_mat(col_idx, 0);
            log_T2 += phi_val * beta_mat(col_idx, 1);
        }
    }

    double c = std::sqrt(x*x + y*y) + 1e-6;
    double f1x = y / c;
    double f1y = -x / c;
    double f2x = x / c;
    double f2y = y / c;

    double exp_L1 = std::exp(log_T1);
    double exp_L2 = std::exp(log_T2);

    dx = exp_L1 * f1x + exp_L2 * f2x;
    dy = exp_L1 * f1y + exp_L2 * f2y;
}

// -------------------------------------------------------
// Main rMAP Loss Function
// -------------------------------------------------------
// [[Rcpp::export]]
double spline_rMAP_loss_cpp(
    const arma::vec& c_params,            // Passed by const reference
    const arma::mat& Phi_data,
    const arma::mat& Phi_quad,
    const arma::mat& Phi_deriv_quad,
    const arma::mat& aug_X,
    const arma::mat& sampledTrajectory,
    const arma::vec& border,
    double pos_sd,
    double vel_sd,
    const arma::vec& rand_vel_x,
    const arma::vec& rand_vel_y,
    const arma::vec& prior_c,
    const arma::mat& inv_c_prior_sigma,
    double Lt
) {

    int N_param = Phi_data.n_cols;
    int N_obs = Phi_data.n_rows;
    int N_quad = Phi_quad.n_rows;

    arma::vec c_x = c_params.subvec(0, N_param - 1);
    arma::vec c_y = c_params.subvec(N_param, 2 * N_param - 1);
    arma::mat C = arma::join_rows(c_x, c_y);

    // =======================================================
    // PART 1: POSITIONAL LOSS
    // =======================================================
    arma::mat X_pred_data = Phi_data * C;
    double nll_pos = 0.0;

    // Quick scalar accumulation
    for(int i = 0; i < N_obs; i++) {
        double err_x = X_pred_data(i, 0) - aug_X(i, 0);
        double err_y = X_pred_data(i, 1) - aug_X(i, 1);
        nll_pos += (err_x * err_x + err_y * err_y);
    }
    nll_pos /= (2.0 * pos_sd * pos_sd);

    // =======================================================
    // PART 2: VELOCITY/PHYSICS LOSS
    // =======================================================
    arma::mat X_quad_pred = Phi_quad * C;
    arma::mat V_spline_pred = Phi_deriv_quad * C;

    double Lx = border[2] - border[0];
    double Ly = border[3] - border[1];
    int num_bases = sampledTrajectory.n_rows;
    int M = std::round(std::sqrt(num_bases)) - 1;
    arma::vec omega = arma::regspace(0, M) * arma::datum::pi;

    double nll_vel = 0.0;

    // Standard single-thread loop
    for(int i = 0; i < N_quad; i++) {
        double dx, dy;

        // Calculate physics flow purely locally (no matrices created)
        single_drift(X_quad_pred(i, 0), X_quad_pred(i, 1), sampledTrajectory,
                     border, omega, M, Lx, Ly, dx, dy);

        double diff_vx = V_spline_pred(i, 0) - (dx + rand_vel_x(i));
        double diff_vy = V_spline_pred(i, 1) - (dy + rand_vel_y(i));

        nll_vel += (diff_vx * diff_vx + diff_vy * diff_vy);
    }

    nll_vel /= (2.0 * vel_sd * vel_sd);
    nll_vel *= (Lt / N_quad);

    // =======================================================
    // PART 3: PRIOR PENALTY LOSS
    // =======================================================
    arma::vec d_cx = c_x - prior_c.subvec(0, N_param - 1);
    arma::vec d_cy = c_y - prior_c.subvec(N_param, 2 * N_param - 1);

    double prior_loss = arma::dot(d_cx, inv_c_prior_sigma * d_cx) +
                        arma::dot(d_cy, inv_c_prior_sigma * d_cy);

    return nll_pos + nll_vel + (0.5 * prior_loss);
}
'

# Compile the C++ code into R
Rcpp::sourceCpp(code = cpp_code)
